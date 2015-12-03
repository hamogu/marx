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
#include <stdlib.h>

#include <jdmath.h>
#include <string.h>

#define MAX_LINE_LEN 64 * 1024


static char *add_to_comment (char *comment, char *addition)
{
   unsigned int len0, len1;
   
   len1 = strlen (addition);
   if (comment == NULL) 
     {
	len0 = 0;
	comment = malloc (len1 + 1);
     }
   else 
     {
	len0 = strlen (comment);
	comment = realloc (comment, len0 + len1 + 1);
     }
   
   if (comment == NULL)
     {
	fprintf (stderr, "Out of memory.\n");
	return NULL;
     }
   
   strcpy (comment + len0, addition);
   return comment;
}

static char *add_pgm_comments (char *comment, int argc, char **argv)
{
   if (NULL == (comment = add_to_comment (comment, 
					  "#This file was generated using the command line:\n#")))
     return NULL;
   
   while (argc)
     {
	argc--;
	if ((NULL == (comment = add_to_comment (comment, *argv)))
	    || (NULL == (comment = add_to_comment (comment, " "))))
	  return NULL;
	
	argv++;
     }
   return add_to_comment (comment, "\n");
}

	  
static void usage (char *pgm)
{   
   fprintf (stderr, "Usage: %s INPUT_FILE OUTPUT_FILE\n", pgm);
}

char Input_Line [MAX_LINE_LEN];

int main (int argc, char **argv)
{
   JDMBData_File_Type *bf;
   unsigned int ncols, nrows, num_cols, linenum;
   int data_type;
   FILE *fpin, *fpout;
   char *pgm;
   char *infile, *outfile;
   char *comment = NULL;

   pgm = argv[0];
   
   if (argc != 3)
     {
	usage (pgm);
	return 1;
     }
   
   infile = argv[1];
   outfile = argv[2];
   
   fpin = fopen (infile, "r");
   if (fpin == NULL)
     {
	fprintf (stderr, "Unable to open %s\n", infile);
	return 1;
     }
   
   
   data_type = 'E';

   comment = NULL;
   while (NULL != fgets (Input_Line, sizeof(Input_Line), fpin))
     {
	if (*Input_Line != '#')
	  {
	     if ((*Input_Line == '!') || (*Input_Line == '%')
		 || (*Input_Line == ';'))
	       *Input_Line = '#';
	     else if (*Input_Line == '\n')
	       continue;
	     else
	       break;
	  }
	
	if (NULL == (comment = add_to_comment (comment, Input_Line)))
	  {
	     fclose (fpin);
	     return 1;
	  }
     }

   if (NULL == (comment = add_pgm_comments (comment, argc, argv)))
     {
	fclose (fpin);
	return 1;
     }

   bf = JDMbdata_create_file (outfile, data_type, 0, 0, comment);

   if (bf == NULL)
     {
	JDMmsg_error ("Error creating output file");
	free (comment);
	fclose (fpin);
	return 1;
     }
   fpout = bf->fp;
   
   linenum = 0;
   nrows = num_cols = 0;

   do
     {
	char *nptr, *endptr;
	float32 val;

	linenum++;

	ncols = 0;
	endptr = Input_Line;
	while (1)
	  {
	     nptr = endptr;
	     val = (float32) strtod (nptr, &endptr);
	     if (nptr == endptr)
	       break;
	     
	     ncols++;

	     if ((ncols > num_cols) 
		 && (num_cols != 0))
	       {
		  fprintf (stderr, "Too many columns on line %u\n", linenum);
		  fclose (fpin);
		  JDMbdata_close_file (bf);
		  return 1;
	       }
	     
	     if (1 != JDMwrite_float32 (&val, 1, fpout))
	       {
		  fprintf (stderr, "Write error\n");
		  fclose (fpin);
		  JDMbdata_close_file (bf);
		  return 1;
	       }
	  }
	
	if (ncols == 0)
	  continue;

	nrows++;

	if (ncols == num_cols)
	  continue;
	
	if (ncols < num_cols)
	  {
	     fprintf (stderr, "Expecting %u columns on line %u\n",
		      num_cols, linenum);
	     fclose (fpin);
	     JDMbdata_close_file (bf);
	     return 1;
	  }
	
	num_cols = ncols;
     }
   while (NULL != fgets (Input_Line, sizeof (Input_Line), fpin));

   fclose (fpin);

   bf->nrows = nrows;
   bf->ncols = num_cols;

   
   if (-1 == JDMbdata_close_file (bf))
     {
	fprintf (stderr, "Error closing file.\n");
	return 1;
     }
   
   if (comment != NULL)
     free (comment);

   return 0;
}

   
   
