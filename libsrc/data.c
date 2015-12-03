/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2002-2018 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research
    under contract SV1-61010 from the Smithsonian Institution.

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
#include "marx-feat.h"

/*{{{ System headers */

#include <stdio.h>
#include <string.h>
#include <limits.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <stdarg.h>
#include <ctype.h>

/*}}}*/

#include "marx.h"
#include "_marx.h"

static char *Data_Directory;

/* This routine returns a malloced quantity.  It is up to the calling routine
 * to free it.
 */
char *marx_make_data_file_name (char *f) /*{{{*/
{
   char *file;
   char *dir;

   dir = Data_Directory;
   if (dir == NULL)
     {
	dir = getenv ("MARX_DATA_DIR");
#ifdef MARX_DATA_DIR
	if (dir == NULL)
	  dir = MARX_DATA_DIR;
#endif
     }

   if (dir == NULL) dir = ".";

   if ((NULL == (file = marx_find_file_in_path (dir, f, ':')))
       && (NULL == (file = marx_find_file_in_path (dir, f, ' '))))
     {
	marx_error ("Unable to locate data file %s\n", f);
	marx_error ("Check MARX_DATA_DIR environment variable");
     }
   return file;
}

/*}}}*/

int marx_set_data_directory (char *dir) /*{{{*/
{
   if (dir == NULL)
     {
#ifdef MARX_DATA_DIR
	dir = MARX_DATA_DIR;
	if (dir == NULL)
#endif
	  return -1;
     }

   if (*dir == '$')
     {
	char *env = getenv (dir + 1);
	if (env == NULL)
	  {
#ifdef MARX_DATA_DIR
	     env = MARX_DATA_DIR;
#endif
	     if (env == NULL)
	       {
		  marx_error ("Unable to set Data Search Path.\nEnvironment variable %s does not exist.",
			      dir + 1);
		  return -1;
	       }
	  }
	dir = env;
     }

   if (*dir == 0)
     {
	marx_error ("Data Path is empty!");
	return -1;
     }

   if (Data_Directory != NULL)
     free (Data_Directory);

   Data_Directory = (char *)marx_malloc (strlen (dir) + 1);
   if (Data_Directory == NULL)
     return -1;

   strcpy (Data_Directory, dir);

   return 0;
}

/*}}}*/

int marx_f_read_bdat (char *file, unsigned int *nrowsptr,
		      unsigned int ncols, float **col0, ...) /*{{{*/
{
   JDMBData_File_Type *bf;
   unsigned int num_cols;
   unsigned int num_rows;
   float **columns;
   float32 *row_data;
   va_list ap;
   unsigned int i, j;
   float **col;
   FILE *fp;

   bf = JDMbdata_open_file (file);
   if (bf == NULL)
     {
	marx_error ("Error reading binary file %s", file);
	return -1;
     }

   if (bf->data_type != 'E')
     {
	marx_error ("%s has unsupported data type", file);
	JDMbdata_close_file (bf);
	return -1;
     }

   num_rows = bf->nrows;
   num_cols = bf->ncols;

   if (ncols > num_cols)
     {
	marx_error ("%s does not have enough columns", file);
	JDMbdata_close_file (bf);
	return -1;
     }

   columns = (float **) marx_malloc (ncols * sizeof (float *));
   if (columns == NULL)
     {
	JDMbdata_close_file (bf);
	return -1;
     }
   memset ((char *) columns, 0, ncols * sizeof(float *));

   if (NULL == (row_data = (float32 *) marx_malloc (num_cols * sizeof(float32))))
     goto return_error;

   if (col0 != NULL)
     {
	if (NULL == (columns[0] = JDMfloat_vector (num_rows)))
	  goto return_error;
     }

   va_start (ap, col0);
   for (i = 1; i < ncols; i++)
     {
	col = va_arg(ap, float **);
	if (col == NULL)
	  continue;

	if (NULL == (columns[i] = JDMfloat_vector (num_rows)))
	  {
	     va_end (ap);
	     goto return_error;
	  }
     }
   va_end (ap);

   fp = bf->fp;

   for (i = 0; i < num_rows; i++)
     {
	if (num_cols != JDMread_float32 (row_data, num_cols, fp))
	  {
	     marx_error ("Read error while reading %s", file);
	     goto return_error;
	  }

	for (j = 0; j < ncols; j++)
	  {
	     if (columns[j] == NULL)
	       continue;

	     columns[j][i] = (float) row_data[j];
	  }
     }

   JDMbdata_close_file (bf);

   if (col0 != NULL) *col0 = columns[0];

   va_start (ap, col0);
   for (i = 1; i < ncols; i++)
     {
	col = va_arg(ap, float **);
	if (col == NULL)
	  continue;

	*col = columns[i];
     }
   va_end (ap);
   marx_free ((char *) columns);
   marx_free ((char *) row_data);

   *nrowsptr = num_rows;
   return 0;

   return_error:

   if (row_data != NULL) marx_free ((char *) row_data);
   if (bf != NULL) JDMbdata_close_file (bf);
   for (i = 0; i < ncols; i++)
     {
	if (columns[i] != NULL)
	  JDMfree_float_vector (columns[i]);
     }
   marx_free ((char *) columns);

   return -1;
}

/*}}}*/

/*{{{ Simple ASCII data file interface */

static _Marx_Simple_Data_Type *
find_data_item (_Marx_Simple_Data_Type *table, char *name)
{
   char ch;

   ch = *name;
   while (table->name != NULL)
     {
	if ((*name == table->name[0])
	    && (0 == strcmp (name, table->name)))
	  return table;

	table++;
     }
   return NULL;
}

static char *skip_whitespace (char *p)
{
   char ch;

   while (((ch = *p) != 0)
	  && isspace (ch))
     p++;

   return p;
}

static char *skip_non_whitespace (char *p)
{
   char ch;

   while (((ch = *p) != 0)
	  && (0 == isspace (ch)))
     p++;

   return p;
}

static int process_data_item (_Marx_Simple_Data_Type *t, char *p0)
{
   char *p1;
   unsigned int n, i;
   double *data;
   double scale;

   n = t->num_elements;
   data = t->data;
   scale = t->scale;

   for (i = 0; i < n; i++)
     {
	double d;

	p0 = skip_whitespace (p0);
	p1 = skip_non_whitespace (p0);

	if (*p1 != 0)
	  *p1++ = 0;

	if (1 != sscanf (p0, "%lf", &d))
	  {
	     marx_error ("Expecting %u numeric fields", i+1);
	     return -1;
	  }

	data[i] = d * scale;
	p0 = p1;
     }

   t->processed = 1;
   return 0;
}

int _marx_read_simple_data_file (char *file, _Marx_Simple_Data_Type *table)
{
   FILE *fp;
   char line[0x3FFF];
   unsigned int linenum;
   _Marx_Simple_Data_Type *t;
   int ret;

   t = table;
   while (t->name != NULL)
     {
	t->processed = 0;
	t++;
     }

   if (NULL == (fp = fopen (file, "r")))
     {
	marx_error ("Unable to open %s", file);
	return -1;
     }

   linenum = 0;
   while (NULL != (fgets (line, sizeof (line), fp)))
     {
	char *p0, *p1;

	linenum++;
	p0 = skip_whitespace (line);
	if (*p0 == '#')
	  continue;

	p1 = skip_non_whitespace (p0);
	if (*p1 != 0) *p1++ = 0;

	if (NULL == (t = find_data_item (table, p0)))
	  continue;

	if (-1 == process_data_item (t, p1))
	  {
	     marx_error ("Error processing line %u of %s", linenum, file);
	     fclose (fp);
	     return -1;
	  }
     }

   fclose (fp);

   ret = 0;
   t = table;
   while (t->name != NULL)
     {
	if (t->processed == 0)
	  {
	     marx_error ("***ERROR: %s not in %s", t->name, file);
	     ret = -1;
	  }
	t++;
     }
   return ret;
}

/*}}}*/
