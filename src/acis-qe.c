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

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>
#include <jdmath.h>

#include <marx.h>
#include <pfile.h>

#include "acissim.c"

static Param_Table_Type Parm_Table [] =
{
     {NULL, 0, NULL}
};

static int acis_qe_initialize (int argc, char **argv)
{
   Param_File_Type *p;

   p = marx_pf_parse_cmd_line ("pileup.par", NULL, argc, argv);

   if (p == NULL)
     {
	fprintf (stderr, "Error opening parameter file.\n");
	return -1;
     }

   if (-1 == acissim_init (p))
     {
	pf_close_parameter_file (p);
	return -1;
     }

   if (-1 == pf_get_parameters (p, Parm_Table))
     {
	pf_error ("Error getting parameters.");
	pf_close_parameter_file (p);
	return -1;
     }

   pf_close_parameter_file (p);

   return 0;
}

int main (int argc, char **argv)
{
   AcisSim_Pixel_Event_Type event;
   AcisSim_Ray_Type ray;
   double energy, emin, emax;
   unsigned int i, npts;
   double ratio, factor, abs_coeff;
   unsigned int num_iterations, num_detected, num_g0_detected;

   if (-1 == acis_qe_initialize (argc, argv))
     return 1;

   memset ((char *) &ray, 0, sizeof (AcisSim_Ray_Type));
   ray.xpixel = 512.5;
   ray.ypixel = 512.5;

   emin = 0.04;
   emax = 12.0;
   ratio = emax / emin;
   npts = 1024;
   factor = 1.0 / (npts - 1);

   num_iterations = 1000;
   /* The 10000 has been chosen because we desire to compute the mean
    * to within 1 percent.  The argument is based on the binomial
    * distribution.  For each energy, n rays are run through the
    * simulator; each one is either detected, or not.  Let p be the
    * probability of detection and let q = 1 - p be the probability
    * that a ray will go undetected. Then after n trials, the mean
    * number detected will be np and the standard deviation will be
    * s = sqrt(npq).  The quantum efficiency is the mean value divided
    * by n and the average sigma per event is s/n.  We want s/n to be
    * less than some value a.  Obviously,
    * @   a < sqrt(npq)/n ==> n > pq/a^2
    * For pq ~ 1 and a ~ 0.01 we find n > 10000.  Note also that max(pq)
    * is <= 1.
    */

   if (Filter_Energies != NULL)
     npts = Num_Filter_Energies;

   for (i = 0; i < npts; i++)
     {
	unsigned int count;
	double t_coeff;

	if (Filter_Energies != NULL)
	  {
	     energy = Filter_Energies[i];
	     t_coeff = Filter_Trans_Coeffs[i];
	  }
	else
	  {
	     energy = emin * pow (ratio, i * factor);
	     t_coeff = 1.0;
	  }

	ray.energy = energy;
	abs_coeff = absorbtion_coeff (energy);

	num_g0_detected = num_detected = 0;
	count = num_iterations;
	while (count)
	  {
	     count--;

	     if ((t_coeff != 1.0)
		 && (JDMrandom () >= t_coeff))
	       continue;

	     if (-1 != process_ray_1 (&ray, &event, abs_coeff))
	       {
		  num_detected++;
		  if ((event.flags & REGULAR_EVENT_OK)
		      && (*(event.regular_island.phas
			    + (event.regular_island.island_size * event.regular_island.island_size)/2)
			  > 0.8 * energy))
		    num_g0_detected++;
	       }
	  }

	fprintf (stdout, "%e\t%e\t%e\n", energy,
		 (double)num_detected/(double)num_iterations,
		 (double)num_g0_detected/(double)num_iterations);
     }
   return 0;
}

