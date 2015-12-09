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

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <jdmath.h>

   /* Now read in each encircled energy array one by one.  Use the array to
    * interpolate a new set of thetas for a common ee grid.  That is, the
    * file contains a mapping of encircled energies as a function of theta
    * where there is a common theta axis.  Simply transform to the reverse---
    * theta as a function of ee.
    */

static void read_error (void) /*{{{*/
{
   fprintf (stderr, "read error.\n");
   exit (-1);
}

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   FILE *fp;
   char *file = "mirr-ee.dat";
   int num_ee;
   int num_energies;
   float value;
   float *energies;
   float *ee;

   int i, n;

   if (NULL == (fp = fopen (file, "rb"))) return -1;

   if ((1 != JDMread_int32 (&num_energies, 1, fp))
       || (1 != JDMread_int32 (&num_ee, 1, fp)))
     read_error ();

   if ((NULL == (ee = JDMfloat_vector (num_ee)))
       || (NULL == (energies = JDMfloat_vector (num_energies))))
     JDMmsg_error (NULL);

   if ((num_energies != JDMread_float32 (energies, num_energies, fp))
       || (num_ee != JDMread_float32 (ee, num_ee, fp)))
     read_error ();

   /* Now interpolate all thetas to this grid */
   for (n = 0; n < 4; n++)
     {
	for (i = 0; i < num_energies; i++)
	  {
	     int j;
	     fprintf (stdout, "#energy: %e\n", energies[i]);
	     for (j = 0; j < num_ee; j++)
	       {
		  if (1 != JDMread_float32 (&value, 1, fp))
		    read_error ();

		  fprintf (stdout, "%e %e\n", ee[j], value);
	       }
	     fprintf (stdout, "@eod\n");
	  }
	fprintf (stdout, "@pause\n@clear\n");
     }

   fclose (fp);
   return 0;
}

/*}}}*/
