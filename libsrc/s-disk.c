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

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"
#include "source.def"

/* Width of source (arc-seconds) */
static double Disk_Theta;
static double Disk_Min_Theta;

static int disk_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int disk_close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int disk_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   int (*efun) (Marx_Spectrum_Type *, double *);
   JDMVector_Type normal, *p;
   double disk_theta, disk_theta_min;
   Marx_Photon_Attr_Type *at;
   double x0, x1;

   /* Convert to radians from arc-secs. */
   disk_theta = Disk_Theta * (1.0 / 3600.0 *  PI / 180.0);    /* 1 arc sec */
   disk_theta_min = Disk_Min_Theta * (1.0 / 3600.0 *  PI / 180.0);    /* 1 arc sec */

   at = pt->attributes;
   efun = st->spectrum.energy_function;

   /* The disk source is a source whose flux is constant. */

   /* Step 1 */
   normal = st->p_normal;
   p = &st->p;

   x0 = disk_theta_min / disk_theta; x0 = x0 * x0;
   x1 = 1.0 - x0;

   for (i = 0; i < num; i++)
     {
	double rnd;

	if (-1 == (*efun) (&st->spectrum, &at->energy))
	  return -1;

	/* Step 2.
	 * Rotate about p.
	 */
	normal = JDMv_rotate_unit_vector (normal, *p, 2.0 * PI * JDMrandom ());

	/* Step 3. */
	rnd = disk_theta * sqrt (x0 + x1 * JDMrandom ());

	at->p = JDMv_rotate_unit_vector (*p, normal, rnd);
	at++;
     }

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_disk_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   (void) source_id;
   st->open_source = disk_open_source;
   st->create_photons = disk_create_photons;
   st->close_source = disk_close_source;

   if ((-1 == pf_get_double (p, "S-DiskTheta1", &Disk_Theta))
       || (-1 == pf_get_double (p, "S-DiskTheta0", &Disk_Min_Theta)))
     return -1;

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

