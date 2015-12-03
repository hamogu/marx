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
/* Read an RDB file. */
/* -*- mode: C; mode: fold; -*- */
#include "config.h"
#include "marx-feat.h"

#include <stdio.h>
#include <string.h>

#include "marx.h"
#include "_marx.h"

/* I have no idea how applicable these routines are with respect to
 * generic RDB files.  The main purpose of the routines in this file is to
 * read the mirror definition files.
 */

typedef struct
{
   char **fields;
   unsigned int num_fields;
   char *field_data;		       /* fields contain pointers to this
					* area.
					*/
}
RDB_Row_Type;

struct _Marx_RDB_File_Type
{
   char *filename;
   FILE *fp;

   RDB_Row_Type column_names_row;
   RDB_Row_Type units_row;

   unsigned int num_rows;
   RDB_Row_Type *rows;
};

/* Returns 0 if at end of file, 1 if line read, -1 upon error */
static int read_line (FILE *fp, char **line_ptr)
{
   char *line;
   unsigned int total_len, space;
   char buf[4096];
   unsigned int padding;

   *line_ptr = NULL;

   while (1)
     {
	unsigned int len;
	int ch;

	if (NULL == fgets (buf, sizeof (buf), fp))
	  return 0;

	if ((*buf == '\n') || (*buf == 0))
	  continue;

	if (*buf != '#')
	  break;

	/* Comment.  Skip past it. */

	len = strlen (buf);
	if (buf[len - 1] != '\n')
	  {
	     while (((ch = getc (fp)) != EOF)
		    && (ch != '\n'))
	       ;
	  }
     }

   space = total_len = 0;
   line = NULL;

   /* If sizeof(buf) is too small, an increase in padding might be a good
    * idea for efficiency.  In any case, 0 is ok.
    */
   padding = 0;

   do
     {
	unsigned int len, needed_len;

	len = strlen (buf);

	needed_len = total_len + len + 1;

	if (needed_len > space)
	  {
	     char *new_line;

	     space = needed_len + padding;

	     if (NULL == (new_line = (char *)marx_realloc (line, space)))
	       {
		  if (line != NULL) marx_free (line);
		  return -1;
	       }
	     line = new_line;
	  }

	strcpy (line + total_len, buf);
	total_len += len;

	if (total_len && (line [total_len - 1] == '\n'))
	  {
	     total_len -= 1;
	     line [total_len] = 0;
	     break;
	  }
     }
   while (NULL != fgets (buf, sizeof (buf), fp));

   *line_ptr = line;
   return 1;
}

static void free_rdb_row_type (RDB_Row_Type *row)
{
   if (row->field_data != NULL)
     marx_free (row->field_data);
   if (row->fields != NULL)
     marx_free ((char *) row->fields);
}

static int read_rdb_row (FILE *fp, RDB_Row_Type *row)
{
   int status;
   char *line, **fields;
   unsigned int num_fields;
   unsigned char ch;

   memset ((char *) row, 0, sizeof (RDB_Row_Type));

   status = read_line (fp, &line);
   if (status != 1)
     return status;

   row->field_data = line;

   num_fields = 1;
   while ((ch = *line) != 0)
     {
	if (ch == '\t')
	  num_fields++;
	line++;
     }

   if (NULL == (fields = (char **)marx_malloc (num_fields * sizeof (char *))))
     {
	free_rdb_row_type (row);
	return -1;
     }
   row->fields = fields;
   row->num_fields = num_fields;

   num_fields = 0;
   line = row->field_data;

   fields[num_fields] = line;
   num_fields++;

   while ((ch = *line) != 0)
     {
	if (ch == '\t')
	  {
	     *line = 0;
	     fields [num_fields] = line + 1;
	     num_fields++;
	  }
	line++;
     }

   return 1;
}

void marx_close_rdb_file (Marx_RDB_File_Type *rdb)
{
   unsigned int r, num_rows;

   if (rdb == NULL)
     return;

   if (rdb->filename != NULL) marx_free (rdb->filename);
   free_rdb_row_type (&rdb->column_names_row);
   free_rdb_row_type (&rdb->units_row);

   if (rdb->rows != NULL)
     {
	num_rows = rdb->num_rows;
	for (r = 0; r < num_rows; r++)
	  free_rdb_row_type (rdb->rows + r);
	marx_free ((char *) rdb->rows);
     }

   if (rdb->fp != NULL) fclose (rdb->fp);
   marx_free ((char *) rdb);
}

Marx_RDB_File_Type *marx_open_rdb_file (char *file)
{
   FILE *fp;
   Marx_RDB_File_Type *rdb;
   unsigned int space;

   if (file == NULL) return NULL;

   fp = fopen (file, "r");
   if (fp == NULL)
     {
	marx_error ("Error opening RDB file %s", file);
	return NULL;
     }

   if (NULL == (rdb = (Marx_RDB_File_Type *) marx_malloc (sizeof (Marx_RDB_File_Type))))
     {
	fclose (fp);
	return NULL;
     }

   memset ((char *)rdb, 0, sizeof (Marx_RDB_File_Type));

   rdb->fp = fp;

   if (NULL == (rdb->filename = (char *)marx_malloc (1 + strlen (file))))
     {
	marx_close_rdb_file (rdb);
	return NULL;
     }
   strcpy (rdb->filename, file);

   if (1 != read_rdb_row (fp, &rdb->column_names_row))
     {
	marx_error ("RDB file %s does not have column names", file);
	marx_close_rdb_file (rdb);
	return NULL;
     }

   if (1 != read_rdb_row (fp, &rdb->units_row))
     {
	marx_error ("RDB file %s does not have a unit column", file);
	marx_close_rdb_file (rdb);
	return NULL;
     }

   space = 0;

   while (1)
     {
	int status;

	if (rdb->num_rows == space)
	  {
	     RDB_Row_Type *new_rows;

	     space += 1024;
	     new_rows = (RDB_Row_Type *) marx_realloc ((char *)rdb->rows,
							space * sizeof (RDB_Row_Type));
	     if (new_rows == NULL)
	       {
		  marx_close_rdb_file (rdb);
		  return NULL;
	       }

	     rdb->rows = new_rows;
	  }

	status = read_rdb_row (fp, rdb->rows + rdb->num_rows);
	if (status == -1)
	  {
	     marx_close_rdb_file (rdb);
	     return NULL;
	  }

	if (status == 0)
	  break;

	rdb->num_rows += 1;
     }

   return rdb;
}

int marx_rdb_get_col (Marx_RDB_File_Type *rdb, char *name)
{
   unsigned int c;
   unsigned int num;
   char ch;
   char **column_names;

   column_names = rdb->column_names_row.fields;
   num = rdb->column_names_row.num_fields;

   ch = *name;
   c = 0;

   while (c < num)
     {
	if ((ch == column_names[c][0])
	    && (0 == strcmp (name, column_names[c])))
	  return (int) c;

	c++;
     }
   marx_error ("File %s has no column with name %s",
	       rdb->filename, name);
   return -1;
}

int marx_rdb_get_row (Marx_RDB_File_Type *rdb, char *name, char *value)
{
   unsigned int c;
   RDB_Row_Type *row, *row_max;

   if (rdb == NULL)
     return -1;

   c = marx_rdb_get_col (rdb, name);
   if (-1 == (int) c)
     return -1;

   /* Now find the row with the specified value */
   row = rdb->rows;
   row_max = row + rdb->num_rows;

   while (row < row_max)
     {
	if ((c < row->num_fields)
	    && (0 == strcmp (value, row->fields[c])))
	  return (int) (row - rdb->rows);

	row++;
     }

   marx_error ("RDB file %s has no row with value %s",
	       rdb->filename, value);
   return -1;
}

char *marx_rdb_get_value (Marx_RDB_File_Type *rdb,
			  unsigned int row, unsigned int col)
{
   RDB_Row_Type *r;

   if (rdb == NULL)
     return NULL;

   if (rdb->num_rows <= row)
     {
	marx_error ("RDB file %s does not have %u rows",
		    rdb->filename, row);
	return NULL;
     }

   r = rdb->rows + row;
   if (r->num_fields <= col)
     {
	marx_error ("RDB file %s row %u does not have %u columns",
		    rdb->filename, row, col);
	return NULL;
     }

   return r->fields[col];
}

#if 0				       /* TESTING */
int main (int argc, char **argv)
{
   Marx_RDB_File_Type *rdb;
   char *value;
   int row, col;

   if (argc == 1)
     return 1;

   rdb = marx_open_rdb_file (argv[1]);
   if (rdb == NULL)
     return 1;

   row = marx_rdb_get_row (rdb, "mirror", "p3");
   col = marx_rdb_get_col (rdb, "az_mis");

   if ((row == -1) || (col == -1))
     {
	marx_close_rdb_file (rdb);
	return 1;
     }

   value = marx_rdb_get_value (rdb, row, col);
   if (value == NULL)
     {
	marx_close_rdb_file (rdb);
	return 1;
     }

   fprintf (stdout, "Value[%d][%d] = %s\n", row, col, value);

   marx_close_rdb_file (rdb);
   return 0;
}
#endif
