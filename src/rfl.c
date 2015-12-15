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

#include <marx.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

static void usage (char *pgm) /*{{{*/
{
   fprintf (stderr, "%s: usage:\n", pgm);
   fprintf (stderr, "\t%s opt-const-file.dat energy(KeV)\n", pgm);
   fprintf (stderr, "or\t%s opt-const-file.dat -max-graxing-angle(degrees)\n", pgm);
   exit (1);
}

/*}}}*/

/* Simple command line tool to dump reflectivities from a binary file to stdout.
 */
int main (int argc, char **argv) /*{{{*/
{
   char *file;
   double energy, theta;
   float *energies, *betas, *deltas;
   unsigned int nread;
   int angle_mode = 0;

   if (argc != 3) usage (argv[0]);

   file = argv[1];
   if (*argv[2] == '-')
     {
	if (1 != sscanf (argv[2], "%lf", &theta))
	  usage (argv[0]);
	theta = -theta;
	angle_mode = 1;
     }
   else if (1 != sscanf (argv[2], "%lf", &energy))
     usage (argv[0]);

   /* The optical constant file consists of:
    *   energy (KeV), beta, delta
    */
   if (-1 == marx_f_read_bdat (file, &nread, 3, &energies, &betas, &deltas))
     {
	fprintf (stderr, "Error encountered trying to read %s\n", file);
	return 1;
     }

   if (angle_mode)
     {
	unsigned int i;
	double t, dt;

	dt = theta / 10.0;
	for (i = 0; i < nread; i++)
	  {
	     fprintf (stdout, "%f", energies[i]);
	     for (t = dt; t <= theta; t += dt)
	       {
		  double cos_theta = cos (PI/2 - PI*t/180.0);

		  fprintf (stdout, "\t%e",
			   marx_reflectivity (cos_theta, betas[i], deltas[i]));
	       }
	     putc ('\n', stdout);
	  }
	return 0;
     }

   fprintf (stdout, "# Energy = %f KeV\n#Arc-Min Probability\n", energy);

   for (theta = 0.0; theta < 600.0; theta += 0.1)
     {
	double t = PI/2.0 - theta * (PI/ 180.0/60.0);
	fprintf (stdout, "%f\t%e\n", theta,
		 marx_interp_reflectivity (energy, cos (t), energies, betas, deltas,
					   nread));
     }

   return 0;
}

/*}}}*/
