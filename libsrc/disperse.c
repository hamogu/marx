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
#include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include "marx.h"

/* This routine computes the position of the photons in the focal plane
 * due to dispersion by the grating.  The geometric quantities specified
 * in this file are determined by the structure of the AXAF.
 */

double HRMA_To_Detector_Distance = 10069.0;    /* mm */
double Grating_To_Detector_Distance = 8635.0;   /* mm */
double Focus_Detector_Offset = 0.2;             /* mm */
/* Note: the Focus_Detector_Offset should really be a time-dependent
 * random number that varies with the position on the focal plane.  Here,
 * it is simply set to the maximum expected distance (worse case).
 */

int disperse_photons (Photon_Type *pt) /*{{{*/
{
   Photon_Attr_Type *photon_attributes, *at;
   unsigned int n, i, *sorted_index;
   Marx_Grating_Assembly_Type *g, grating_assembly[NUM_MIRRORS];
   double sines[NUM_MIRRORS];
   double cosines[NUM_MIRRORS];
   double detector_blur[NUM_MIRRORS];
   double mdist = HRMA_To_Detector_Distance + Focus_Detector_Offset;
   double gdist = Grating_To_Detector_Distance + Focus_Detector_Offset;
   double factor = HBAR_C * 2.0 * PI;
   double blur_ratio;

   if (Grating_To_Detector_Distance <= 0.0)
     {
	blur_ratio = 1.0e10;
     }
   else blur_ratio = Focus_Detector_Offset / Grating_To_Detector_Distance;

   /* precompute sin/cos for efficiency */
   for (i = 0; i < NUM_MIRRORS; i++)
     {
	diffract_update_geometry (grating_assembly + i, i);

	sines[i] = sin(grating_assembly[i].dispersion_angle);
	cosines[i] = cos(grating_assembly[i].dispersion_angle);

	detector_blur[i] = blur_ratio * grating_assembly[i].grating_radius;
     }

   prune_photons (pt);
   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   for (i = 0; i < n; i++)
     {
	unsigned int mirror_shell;
	double c, s;
	float dy, dz;
	int order;

	at = photon_attributes + sorted_index[i];
	mirror_shell = at->mirror_shell;
	g = &grating_assembly[mirror_shell];

	order = at->order;

	if (order != 0)
	  {
	     dy = order * ((factor / at->energy) / g->grating_period) * gdist;
	     dz = 0.0;
	  }
	else
	  {
	     dy = dz = 0.0;
	  }

	/* Grating blur will have to be added here: */

	/* Here is the blur due to the detector not exactly at the focal
	 * plane.
	 */
	dy += detector_blur[mirror_shell];

	/* Now rotate */
	c = cosines[mirror_shell];
	s = sines[mirror_shell];
	at->y =  dy * c + dz * s;
	at->z = -dy * s + dz * c;

	/* Add in mirror blur terms */
	at->y += at->blur_theta * mdist * cos((double)at->azimuth);
	at->z += at->blur_theta * mdist * sin((double)at->azimuth);
     }
   return 0;
}

/*}}}*/

int dump_disperse_variables (FILE *fp) /*{{{*/
{
   fprintf (fp, "--Geometric Variables--\n");
   fprintf (fp, " HRMA_To_Detector_Distance: %f\n", HRMA_To_Detector_Distance);
   fprintf (fp, " Grating_To_Detector_Distance: %f\n", Grating_To_Detector_Distance);
   fprintf (fp, " Focus_Detector_Offset: %f\n", Focus_Detector_Offset);
   return 0;
}

/*}}}*/

