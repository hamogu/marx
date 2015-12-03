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

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#include <string.h>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <limits.h>

#ifndef PATH_MAX
# define PATH_MAX 1024
#endif

#define MAX_PATH_LEN PATH_MAX

#include "marx.h"
#include "_marx.h"

double _Marx_TStart_Yrs;
double _Marx_TStart_MJDsecs;

static void do_malloc_error (void) /*{{{*/
{
   marx_error ("Memory allocation failure.");
}

/*}}}*/

void *marx_malloc (unsigned int nbytes) /*{{{*/
{
   void *s;

   s = malloc (nbytes);
   if (s == NULL)
     do_malloc_error ();
   return s;
}

/*}}}*/

void *marx_calloc (unsigned int n, unsigned int size) /*{{{*/
{
   void *s;

   s = calloc (n, size);
   if (s == NULL)
     do_malloc_error ();
   return s;
}

/*}}}*/

void *marx_realloc (void *ptr, unsigned int size)
{
   if (ptr == NULL)
     return marx_malloc (size);

   ptr = realloc (ptr, size);
   if (ptr == NULL)
     do_malloc_error ();
   return ptr;
}

int marx_free (void *s) /*{{{*/
{
   if (s == NULL) return -1;
   free (s);
   return 0;
}

/*}}}*/

/* Concatenate dir and file to produce dir/file.  This routine works with or
 * without the trailing / on dir.  If file is "", then dir is simply
 * returned.
 */
char *marx_dircat (char *dir, char *file) /*{{{*/
{
   unsigned int dirlen;
   unsigned int filelen;
   char *filename;

   if ((dir == NULL) && (file == NULL))
     return NULL;

   dirlen = filelen = 0;

   if (dir != NULL) dirlen = strlen (dir);
   if (file != NULL) filelen = strlen (file);

   if (NULL == (filename = (char *)malloc (dirlen + filelen + 2)))
     {
	do_malloc_error ();
	return NULL;
     }

   if (dirlen)
     {
	strcpy (filename, dir);
	/* Add final / if it is not already there. */
	if (filename[dirlen - 1] != '/')
	  {
	     filename[dirlen] = '/';
	     dirlen++;
	     filename[dirlen] = 0;
	  }
     }
   if (filelen) strcpy (filename + dirlen, file);
   return filename;
}

/*}}}*/

int marx_file_exists (char *file)
{
   struct stat st;
   int m;

#ifdef _S_IFDIR
# ifndef S_IFDIR
#  define S_IFDIR _S_IFDIR
# endif
#endif

#ifndef S_ISDIR
# ifdef S_IFDIR
#  define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
# else
#  define S_ISDIR(m) 0
# endif
#endif

   if (stat(file, &st) < 0) return 0;
   m = st.st_mode;

   if (S_ISDIR(m)) return 2;
   return 1;
}

char *marx_find_file_in_path (char *dirs, char *file, char colon)
{
   char dir_buf [MAX_PATH_LEN + 1];
   char *dir0, *dir1;
   char *filename;

   /* If file is absolute, try to find it. */
   if (file == NULL) return NULL;

   if (*file == '/')
     {
	if (1 == marx_file_exists (file))
	  return marx_dircat (NULL, file);
     }

   if (dirs == NULL)
     return NULL;

   dir0 = dirs;

   while (*dir0)
     {
	unsigned int len;

	dir1 = strchr (dir0, colon);
	if (dir1 == NULL)
	  {
	     dir1 = dir0 + strlen (dir0);
	  }
	if (dir0 == dir1)
	  {
	     dir0++;
	     continue;
	  }
	if ((len = (dir1 - dir0)) <= MAX_PATH_LEN)
	  {
	     strncpy (dir_buf, dir0, len);
	     dir_buf[len] = 0;

	     filename = marx_dircat (dir_buf, file);
	     if (filename == NULL)
	       return NULL;

	     if (1 == marx_file_exists (filename))
	       return filename;

	     free (filename);
	  }

	dir0 = dir1;
	if (*dir0) dir0++;
     }

   return NULL;
}

int marx_set_time (double tyrs, double mjdsecs)
{
   _Marx_TStart_Yrs = tyrs;
   _Marx_TStart_MJDsecs = mjdsecs;
   return 0;
}

Param_File_Type *marx_pf_parse_cmd_line
  (char *file, char *mode, int argc, char **argv)
{
   Param_File_Type *pf;

   pf = pf_parse_cmd_line_no_exit (file, mode, argc, argv);
   if (pf != NULL)
     return pf;
#ifdef MARX_PFILE_DIR
   if (PF_Errno == PF_FILE_NOT_FOUND)
     {
	fprintf (stderr, "*****\n\
Unable to locate a %s parameter file.  Marx parameter files may be\n\
found in:\n\
 %s/\n\
*****\n",
		 file, MARX_PFILE_DIR);
     }
#endif
   exit (1);
}

int _marx_check_monotonicity_f (float *p, unsigned int n)
{
   float x, *pmax;

   pmax = p + n;
   x = *p++;
   while (p < pmax)
     {
	if (*p < x)
	  return -1;
	x = *p++;
     }
   return 0;
}

int _marx_check_monotonicity_d (double *p, unsigned int n)
{
   double x, *pmax;

   pmax = p + n;
   x = *p++;
   while (p < pmax)
     {
	if (*p < x)
	  return -1;
	x = *p++;
     }
   return 0;
}


/* Range syntax:
 *
 *   ranges:
 *      range
 *      range "," ranges
 *
 *   range:
 *      number_optional
 *      lowerbndchar_optional number_optional ":" number_optional upperbndchar_optional
 *
 *   lowerbndchar:
 *     "(" | "["
 *   upperbndchar:
 *     "]" | ")"
 *
 *  examples:
 *       2:5   ==> 2 <= x < 5
 *       2:5]  ==> 2 <= x <= 5
 *      (2:    ==> 2 < x
 *       :5]   ==> x <= 5
 *       :     null range
 *      2:)    ==> 2 <= x
 *      (2)    ==> 2 < x < 2 (not 2)
 *
 *  Note: Infinity is NOT handled.  
 */

struct _Marx_Range_Type
{
   double a, b;
   int type;
#define RANGE_TYPE_EMPTY	0
#define RANGE_TYPE_LT_LT	1      /* "(a,b)"   a < x < y  */
#define RANGE_TYPE_LT_LE	2      /* "(a,b]"   a < x <= y */
#define RANGE_TYPE_LE_LT	3      /* "a:b"     a <= x < y */
#define RANGE_TYPE_LE_LE	4      /* "a:b]"    a <= x <= y */
#define RANGE_TYPE_LT		5      /* "(a:"     a < x */
#define RANGE_TYPE_LE		6      /* "[a:"     a <= x */
#define RANGE_TYPE_GT		7      /* ":a"      a > x */
#define RANGE_TYPE_GE		8      /* ":a]"     a >= x */
   struct _Marx_Range_Type *next;
};

int _marx_is_in_range (Marx_Range_Type *r, double x)
{
   while (r != NULL)
     {
	switch (r->type)
	  {
	   case RANGE_TYPE_EMPTY: return 1;
	   case RANGE_TYPE_LT_LT: if ((r->a < x) && (x < r->b)) return 1; break;
	   case RANGE_TYPE_LT_LE: if ((r->a < x) && (x <= r->b)) return 1; break;
	   case RANGE_TYPE_LE_LT: if ((r->a <= x) && (x < r->b)) return 1; break;
	   case RANGE_TYPE_LE_LE: if ((r->a <= x) && (x <= r->b)) return 1; break;
	   case RANGE_TYPE_LT: if (r->a < x) return 1; break;
	   case RANGE_TYPE_LE: if (r->a <= x) return 1; break;
	   case RANGE_TYPE_GT: if (r->a > x) return 1; break;
	   case RANGE_TYPE_GE: if (r->a <= x) return 1; break;
	  }
	r = r->next;
     }
   return 0;
}

double _marx_compute_range_length (Marx_Range_Type *r)
{
   double sum = 0.0;

#define UNBOUNDED_RANGE -1.0

   /* Note: this does not handle overlapping sectors! */
   while (r != NULL)
     {
	switch (r->type)
	  {
	   case RANGE_TYPE_LT_LT:
	   case RANGE_TYPE_LT_LE:
	   case RANGE_TYPE_LE_LT:
	   case RANGE_TYPE_LE_LE:
	     if (r->b > r->a)
	       sum += r->b - r->a;
	     break;
	   case RANGE_TYPE_LT:
	   case RANGE_TYPE_LE:
	   case RANGE_TYPE_GT:
	   case RANGE_TYPE_GE:
	   case RANGE_TYPE_EMPTY:
	     return UNBOUNDED_RANGE;
	  }
	r = r->next;
     }
   return sum;
}


void _marx_free_range_type (Marx_Range_Type *r)
{
   while (r != NULL)
     {
	Marx_Range_Type *next = r->next;
	marx_free ((char *)r);
	r = next;
     }
}

static int parse_number (char **strp, double *xp)
{
   char *s0, *s1;
   double x;

   s0 = *strp;
   while (isspace (*s0)) s0++;

   x = strtod (s0, &s1);
   *strp = s1;
   *xp = x;

   if (s0 == s1)
     return 0;

   return 1;
}

Marx_Range_Type *_marx_parse_range_string (char *str)
{
   Marx_Range_Type *root = NULL, *tail;

   if (str != NULL) while (1)
     {
	char ch;
	int type;
#define RANGEOPEN 1
#define RANGECLOSED 2
	int a_type, b_type, has_a, has_b;
	double a = 0, b = 0;
	Marx_Range_Type *r;

	ch = *str;
	if (ch == 0)
	  break;

	if (isspace (ch) || (ch == ','))
	  {
	     str++;
	     continue;
	  }

	/* at the start of a new range */
	if (ch == '(')
	  {
	     a_type = RANGEOPEN;
	     str++;
	  }
	else if (ch == '[')
	  {
	     a_type = RANGECLOSED;
	     str++;
	  }
	else a_type = RANGECLOSED;

	has_a = parse_number (&str, &a);

	while (isspace (*str))
	  str++;

	if (*str == ':')
	  {
	     str++;
	     has_b = parse_number (&str, &b);
	     while (isspace (*str))
	       str++;
	  }
	else has_b = 0;

	ch = *str;
	if (ch == ')')
	  {
	     b_type = RANGEOPEN;
	     str++;
	     while (isspace (*str)) str++;
	  }
	else if (ch == ']')
	  {
	     b_type = RANGECLOSED;
	     str++;
	     while (isspace (*str)) str++;
	  }
	else b_type = RANGEOPEN;

	if ((*str != 0) && (*str != ','))
	  {
	     marx_error ("While parsing the range specification, expected a comma or end of string, found: '%c'", *str);
	     goto return_error;
	  }

	if (has_a)
	  {
	     if (has_b)
	       {
		  if (a_type == RANGEOPEN)
		    {
		       if (b_type == RANGEOPEN)
			 type = RANGE_TYPE_LT_LT;
		       else
			 type = RANGE_TYPE_LT_LE;
		    }
		  else if (b_type == RANGEOPEN)
		    type = RANGE_TYPE_LE_LT;
		  else
		    type = RANGE_TYPE_LE_LE;
	       }
	     else if (a_type == RANGEOPEN)
	       type = RANGE_TYPE_LT;
	     else
	       type = RANGE_TYPE_LE;
	  }
	else if (has_b)
	  {
	     if (b_type == RANGEOPEN)   /* ":a)" x < a ==> a>x */
	       type = RANGE_TYPE_GT;
	     else 		       /* ":a]" ==> x <= a */
	       type = RANGE_TYPE_GE;
	     a = b;
	  }
	else  /* empty range */
	  type = RANGE_TYPE_EMPTY;

	r = (Marx_Range_Type *)marx_malloc(sizeof (Marx_Range_Type));
	if (r == NULL)
	  goto return_error;
	r->next = NULL;
	r->type = type;
	r->a = a;
	r->b = b;

	if (root == NULL)
	  root = r;
	else
	  tail->next = r;

	tail = r;
     }

   if (root == NULL)
     {
	root = (Marx_Range_Type *)marx_malloc(sizeof (Marx_Range_Type));
	if (root == NULL)
	  return NULL;
	root->type = RANGE_TYPE_EMPTY;
	root->next = NULL;
     }
   return root;

return_error:
   _marx_free_range_type (root);
   return NULL;
}
