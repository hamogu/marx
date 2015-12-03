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

static int beta_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
#if 0
   unsigned int i;
   for (i = 0; i < st->num_radial_points; i++)
     {
	fprintf (stdout, "%e %e\n", st->cum_radial_rvalues[i],
		 st->cum_radial_dist[i]);
     }
#endif
   return 0;
}

/*}}}*/

static int beta_close_source (Marx_Source_Type *st) /*{{{*/
{
   if (st->cum_radial_dist != NULL)
     {
	JDMfree_double_vector (st->cum_radial_dist);
	st->cum_radial_dist = NULL;
     }
   if (st->cum_radial_rvalues != NULL)
     {
	JDMfree_double_vector (st->cum_radial_rvalues);
	st->cum_radial_rvalues = NULL;
     }

   return 0;
}

/*}}}*/

static double Core_Radius;	       /* radians */
static double Beta;
static double Alpha;

static int beta_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				 unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   double rnd;
   JDMVector_Type normal, *p;
   Marx_Photon_Attr_Type *at;
   int (*efun) (Marx_Spectrum_Type *, double *);
   double xpon;

   efun = st->spectrum.energy_function;
   at = pt->attributes;

   /* See comment in s-gauss.c regarding the procedure there.  Only the
    * details in step 3 differ.
    */
   normal = st->p_normal;
   p = &st->p;

   xpon = 1.0 / (1.0 - Alpha);

   for (i = 0; i < num; i++)
     {
	if (-1 == (*efun) (&st->spectrum, &at->energy))
	  return -1;

	/* Step 2.
	 * Rotate normal about p.
	 */
	normal = JDMv_rotate_unit_vector (normal, *p, 2.0 * PI * JDMrandom ());

	/* step 3. Rotate p about this normal by an angle according to
	 * the beta distribution.  This is accomplished by getting a random
	 * number between 0 and 1 and looking up the corresponding x value
	 * in the cum_radial distribution.
	 */
	do
	  {
	     rnd = JDMrandom ();
	  }
	while (rnd == 0.0);

	rnd = Core_Radius * sqrt (pow (rnd, xpon) - 1.0);

	/* Step 3. */
	at->p = JDMv_rotate_unit_vector (*p, normal, rnd);
	at++;
     }

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_beta_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			     char *name, unsigned int source_id)
{
   (void) source_id;
   st->open_source = beta_open_source;
   st->create_photons = beta_create_photons;
   st->close_source = beta_close_source;

   if (-1 == pf_get_double (p, "S-BetaCoreRadius", &Core_Radius))
     return -1;

   if (-1 == pf_get_double (p, "S-BetaBeta", &Beta))
     return -1;

   Alpha = 3 * Beta - 0.5;

   if ((Alpha <= 1.0) || (Core_Radius <= 0.0))
     {
	marx_error ("Beta parameters out of range.  Must be greater than 1/2.");
	return -1;
     }

   /* Convert core radius from seconds to radians */
   Core_Radius = Core_Radius * (1.0 / 3600.0) * (PI / 180.0);

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

