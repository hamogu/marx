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

/* The gaussian distribution use here is defined to be one in which the
 * BRIGHTNESS is distributed with a probability distribution of the form
 * g(r) = exp (-r^2/s^s) where s (sigma) controls the width.
 * Note that g(r)dr is NOT the number of rays with r between r and r + dr.
 * The number of rays with r between r and r + dr is r g(r) dr and the total
 * number of rays within radius r is
 * @
 * @  N(r) = \int_0^r dr r g(r) = 1 - exp(-r^2/s^2)
 * @
 * Inverting this yields:
 * @
 * @  r(N) = s sqrt(log (1/(1 - N(r))))
 * @
 * Note that N(r) is normalized to 1; that is, it is the fractional number.
 */

/* Half width of source (arc-seconds) */
static double Sigma_Theta = 0.0;

static int gauss_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int gauss_close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int gauss_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				 unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   int (*efun) (Marx_Spectrum_Type *, double *);
   JDMVector_Type normal, *p;
   double sigma_theta;
   Marx_Photon_Attr_Type *at;

   /* Convert to radians from arc-secs. */
   sigma_theta = Sigma_Theta * (1.0 / 3600.0 *  PI / 180.0);    /* 1 arc sec */

   at = pt->attributes;
   efun = st->spectrum.energy_function;

   /* The gaussian point source consists of rays whose directions are
    * normally distributed about the direction to the origin.  To compute
    * the direction, we need proceed in 3 steps:
    *
    *    1.  Construct a unit vector normal to the unit vector that
    *        points to the origin (st->p).  This step should have already
    *        been performed by calling routine (st->p_normal).
    *
    *    2.  Rotate this vector by a random angle about st->p to produce
    *        a new normal.
    *    3.  Rotate st->p about the normal of step 2 by a normally
    *        distributed random angle.
    *
    * Only the first step is independent of the ray.  The other steps
    * must be performed on every ray.
    */

   /* Step 1 */
   normal = st->p_normal;
   p = &st->p;

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

	do
	  {
	     rnd = JDMrandom ();
	  }
	while (rnd == 0.0);

	rnd = sqrt (-log (rnd));

	at->p = JDMv_rotate_unit_vector (*p, normal,
					 sigma_theta * rnd);
	at++;
     }

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_gauss_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   (void) source_id;
   st->open_source = gauss_open_source;
   st->create_photons = gauss_create_photons;
   st->close_source = gauss_close_source;

   if (-1 == pf_get_double (p, "S-GaussSigma", &Sigma_Theta))
     return -1;

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

