/*
    Copyright (C) 2002 MIT Center For Space Research

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include "config.h"

#include <stdio.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <string.h>

#include "jdfits.h"

static void usage (void)
{
   fprintf (stderr, "fitsdump v%s\n",
	    JDFITS_VERSION);

   fprintf (stderr, "\
Usage: fitsdump [Options] [-e <extension name>] <fitsfile>\n\
 Options:\n\
  -t                    test fits file integrity\n\
  -r                    dump raw data, do not scale\n\
  -s                    apply full scaling\n\
  -h                    dump headers\n\
  -H                    dump only headers\n\
  -e <extension name/number>\n\
                        dump out a particular extension.  The default is all.\n\
  -c <column name>      dump out a particular column with specified name\n"
	    );
   exit (1);
}

static char *Extension_Name;
static unsigned int Option_Flags;

#define RAW_FLAG 0x1
#define FULL_SCALE 0x2
#define DUMP_HEADERS 0x4
#define HEADERS_ONLY 0x8
#define TEST_INTEGRITY 0x10

static void parse_options (char *str)
{
   char ch;
   while ((ch = *str++) != 0)
     {
	switch (ch)
	  {
	   case 'r':
	     Option_Flags |= RAW_FLAG;
	     break;

	   case 's':
	     Option_Flags |= FULL_SCALE;
	     break;

	   case 'H':
	     Option_Flags |= (HEADERS_ONLY | DUMP_HEADERS);
	     break;

	   case 'h':
	     Option_Flags |= DUMP_HEADERS;
	     break;

	   case 't':
	     Option_Flags |= (TEST_INTEGRITY | DUMP_HEADERS);
	     break;

	   default:
	     fprintf (stderr, "fitsdump: Undefined option: %c\n", ch);
	     usage ();
	  }
     }
}

static int fitsdump_file (char *file, int ext_to_dump, char *colname, int prefix_with_filename)
{
   JDFits_Type *ft;
   JDFits_Keyword_Type *kw;
   JDFits_Read_Keyword_Type *r;
   int ok_to_dump, extnum;
   char *header_fmt = NULL;

   ft = jdfits_open_file (file, JDFITS_READ_MODE);

   if (ft == NULL)
     {
	fprintf (stderr, "Unable to open %s as fits file.\n", file);
	return -1;
     }

   extnum = 0;
   if ((Option_Flags & TEST_INTEGRITY) == 0)
     if (Extension_Name == NULL)
       header_fmt = "# EXTNAME[%d]: %s\n";

   while (1)
     {
	char *extname = NULL;
	if (Extension_Name == NULL) ok_to_dump = 1;
	else ok_to_dump = (ext_to_dump == extnum);

	if (ft->header->type == JDFITS_BINTABLE_HEADER)
	  {

	     if (-1 == jdfits_bintable_parse_headers (ft)) break;

	     extname = ft->header->ext.bintable->extname;

	     if (Extension_Name != NULL)
	       {
		  if (ext_to_dump == -1)
		    {
		       if ((extname != NULL) && !strcmp (Extension_Name, extname))
			 ok_to_dump = 1;
		    }
	       }
	  }

	if (Option_Flags & TEST_INTEGRITY) ok_to_dump = 1;

	if (extname == NULL)
	  extname = "(none)";

	if (header_fmt != NULL)
	  fprintf (stdout, header_fmt, extnum, extname);

	if (ok_to_dump
	    && (Option_Flags & DUMP_HEADERS)
	    && (NULL != (r = jdfits_open_keywords (ft))))
	  {
	     while ((kw = jdfits_read_keyword (r)) != NULL)
	       {
		  int ret = 0;

		  if (Option_Flags & TEST_INTEGRITY) continue;

		  if (prefix_with_filename
		      && (fprintf (stdout, "%s:", file) < 0))
		    exit (1);

		  if (fprintf (stdout, "%-8s", kw->name) < 0)
		    exit (1);

		  switch (kw->type)
		    {
		     case JDFITS_FLOAT64_TYPE:
		       ret = fprintf (stdout, "=%20.16g", kw->v.dval);
		       break;

		     case JDFITS_FLOAT32_TYPE:
		       ret = fprintf (stdout, "=%20.16g", kw->v.fval);
		       break;

		     case JDFITS_INT32_TYPE:
		       ret = fprintf (stdout, "=%20ld", (long) kw->v.lval);
		       break;

		     case JDFITS_INT_TYPE:
		       ret = fprintf (stdout, "=%20d", kw->v.ival);
		       break;
		     case JDFITS_BOOL_TYPE:
		       ret = fprintf (stdout, "=%20c", kw->v.ival ? 'T' : 'F');
		       break;

		     case JDFITS_INT16_TYPE:
		       ret = fprintf (stdout, "=%20hd", kw->v.hval);
		       break;

		     case JDFITS_STRING_TYPE:
		       ret = fprintf (stdout, "= '%-8s'", kw->v.sval);
		       break;

		     case JDFITS_COMMENT_TYPE:
		       ret = fprintf (stdout, "%s", kw->v.sval);
		       break;

		     default:
		       break;
		    }

		  if (ret < 0)
		    exit (1);

		  if ((kw->comment != NULL) && (kw->type != JDFITS_COMMENT_TYPE))
		    {
		       putc (' ', stdout);
		       fwrite ((char *) kw->comment, 1, kw->comment_len, stdout);
		    }
		  if (EOF == putc ('\n', stdout))
		    exit (1);
	       }

	     jdfits_close_keywords (r);
	  }
	if (Option_Flags & TEST_INTEGRITY) ok_to_dump = 0;

	if (ok_to_dump && (ft->header->type == JDFITS_BINTABLE_HEADER)
	    && ((Option_Flags & HEADERS_ONLY) == 0))
	  {
	     int scale = 0;
	     if (Option_Flags & RAW_FLAG) scale = -1;
	     else if (Option_Flags & FULL_SCALE) scale = 1;

	     if (-1 == jdfits_bintable_dump_data (ft, scale, stdout, colname))
	       {
		  jdfits_close_file (ft);
		  exit (1);
	       }
	  }
	else if (jdfits_skip_to_next_header (ft) == -1)
	  break;

	if (jdfits_read_header (ft) == -1)
	  break;
	/* if (jdfits_parse_header (ft) == -1) break; */

	extnum++;
     }
   jdfits_close_file (ft);
   return 0;
}

int main (int argc, char **argv)
{
   int i;
   char *colname = NULL;
   int ext_to_dump = -1;

   /* Parse command line arguments */
   for (i = 1; i < argc; i++)
     {
	if (!strcmp ("-e", argv[i]))
	  {
	     i++;
	     if (i == argc) usage ();
	     Extension_Name = argv[i];
	     if (1 != sscanf (Extension_Name, "%d", &ext_to_dump))
	       ext_to_dump = -1;
	  }

	else if (!strcmp ("-c", argv[i]))
	  {
	     i++;
	     if (i == argc) usage ();
	     colname = argv[i];
	  }

	else if (*(argv[i]) == '-') parse_options (*(argv + i) + 1);
	else break;
     }

   if (i == argc)
     usage ();

   argc -= i;
   argv += i;

#ifdef IDL_HACK
   return jdfits_to_idl (argv[0], NULL);
#endif

   for (i = 0; i < argc; i++)
     {
	char *file = argv[i];

	if (argc > 1)
	  fprintf (stdout, "==>>> %s:\n", file);

	if (-1 == fitsdump_file (file, ext_to_dump, colname,
				 (argc>1) && (Option_Flags & HEADERS_ONLY)))
	  fprintf (stderr, "Error dumping %s\n", file);
     }

   return 0;
}

