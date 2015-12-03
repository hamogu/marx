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

int main (int argc, char **argv) /*{{{*/
{
   FILE *fp, *fpout;
   char *file = "Mirror_EE.bin";
   float32 theta[564];
   float32 new_theta[564];
   float32 energy[25];
   float32 ee_n[4][25][564];
   int32 num_ee = 564;
   int32 num_energies = 25;
   float32 new_ee[564];
   unsigned int i, n;

   (void) argc; (void) argv;

   if (NULL == (fp = fopen (file, "rb"))) return -1;

   if (25 != JDMread_float32 (energy, 25, fp))
     return -1;

   if (564 != JDMread_float32 (theta, 564, fp))
     return -1;

   if (4 * 564 * 25 != JDMread_float32 ((float32 *)ee_n, 4 * 564 * 25, fp))
     return -1;

   fclose (fp);

   /* Fix up data.  The EEs should be an increasing function ending at 1.0!!
    */
   for (n = 0; n < 4; n++)
     {
	for (i = 0; i < 25; i++)
	  {
	     float min_ee = 0.0;
	     unsigned int j;
	     for (j = 0; j < 563; j++)
	       {
		  float this_ee = ee_n[n][i][j];
		  if (this_ee < min_ee) ee_n[n][i][j] = min_ee;
		  else min_ee = this_ee;
	       }
	     ee_n[n][i][563] = 1.0;
	  }
     }

   if (NULL == (fpout = fopen ("mirr-ee.bin", "wb")))
     {
	fprintf (stderr, "Unable to open output file.\n");
	exit (-1);
     }

   JDMwrite_int32 (&num_energies, 1, fpout);
   JDMwrite_int32 (&num_ee, 1, fpout);

   if (num_energies != (int)JDMwrite_float32 (energy, (unsigned int) num_energies, fpout))
     {
	fprintf (stderr, "Write error.\n");
	return -1;
     }

   /* Create a new ee grid.  Since the function goes up very fast and levels
    * out at 1.0, use a log grid to give more samples near 1.0.
    */
   for (i = 1; i <= (unsigned int) num_ee; i++)
     {
	new_ee[i - 1] = (float) (log ((double)i) / log((double) num_ee));
     }

   if (num_ee != (int)JDMwrite_float32 (new_ee, (unsigned int) num_ee, fpout))
     {
	fprintf (stderr, "Write error.\n");
	return -1;
     }

   /* Now interpolate all thetas to this grid */
   for (n = 0; n < 4; n++)
     {
	for (i = 0; i < 25; i++)
	  {
	     unsigned int j;
	     JDMinterpolate_fvector (new_ee, new_theta, (unsigned int) num_ee,
				     (float *)ee_n[n][i], theta, (unsigned int) num_ee);

	     /* Convert the theta array to radians. */
	     for (j = 0; j < (unsigned int) num_ee; j++)
	       {
		  new_theta[j] = (new_theta[j] / 3600.0) * (PI / 180.0);
		  if (new_theta[j] < 0.0) new_theta[j] = 0.0;
	       }

	     if (num_ee != (int)JDMwrite_float32 (new_theta, (unsigned int) num_ee, fpout))
	       {
		  fprintf (stderr, "Write error.\n");
		  return -1;
	       }
	  }
     }
#if 0
   for (i = 0; i < num_ee; i++)
     {
	int j;
	for (j = 0; j < num_energies; j++)
	  {
	     fprintf (stdout, "%e %e %e %e %e %e\n",

		      new_theta[i],
		      new_ee[0][3][i], new_ee[1][3][i],
		      new_ee[2][3][i], new_ee[3][3][i]);
	  }
	fputs ("\n\n", stdout);
     }

#endif

   fclose (fpout);
   return 0;
}

/*}}}*/
