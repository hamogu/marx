/*
    This file is part of MARX

    Copyright (C) 1999 Massachusetts Institute of Technology

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
#include <stdio.h>
#include <jdmath.h>

#include <string.h>

static int dump_float32 (JDMBData_File_Type *bf);

int main (int argc, char **argv)
{
   char *file;
   JDMBData_File_Type *bf;
   char *pgm;
   
   pgm = argv[0];

   if (argc != 2)
     {
	fprintf (stderr, "Usage: %s FILENAME\n", pgm);
	return 1;
     }
   
   file = argv[1];
   
   bf = JDMbdata_open_file (file);
   if (bf == NULL)
     {
	fprintf (stderr, "Unable to open %s\n", file);
	return 1;
     }
   
   fprintf (stdout, "#Data Type: '%c'\n", bf->data_type);
   fprintf (stdout, "#Num Rows: %u\n", bf->nrows);
   fprintf (stdout, "#Num Cols: %u\n", bf->ncols);
   fprintf (stdout, "#Comment: %s\n", bf->comment);
   
   switch (bf->data_type)
     {
      case 'E':
	if (-1 == dump_float32 (bf))
	  return 1;
	break;
	
      default:
	fprintf (stderr, "*** %s: Support for this data type has not been implemented.\n",
		 pgm);
	return 1;
     }

   return 0;
}

static int dump_float32 (JDMBData_File_Type *bf)
{
   unsigned int nrows, ncols;
   unsigned int i, j;
   FILE *fp;
   
   fp = bf->fp;
   nrows = bf->nrows;
   ncols = bf->ncols;
   
   for (i = 0; i < nrows; i++)
     {
	for (j = 0; j < ncols; j++)
	  {
	     float32 x;
	     
	     if (1 != JDMread_float32 (&x, 1, fp))
	       {
		  fprintf (stderr, "*** Read error (row %u, col %u)\n", i, j);
		  return -1;
	       }
	     
	     fprintf (stdout, "%e\t", (double) x);
	  }
	fputc ('\n', stdout);
     }
   
   return 0;
}

	     
