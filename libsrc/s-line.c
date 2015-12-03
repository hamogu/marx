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

static double Line_Phi = 0.0;	       /* degrees */
static double Line_Theta = 1.0;	       /* arc seconds */

static int line_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int line_close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int line_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				 unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   int (*efun) (Marx_Spectrum_Type *, double *);
   double theta;
   double cos_phi, sin_phi;
   double alpha;
   Marx_Photon_Attr_Type *at;
   JDMVector_Type normal;

   at = pt->attributes;
   efun = st->spectrum.energy_function;

   cos_phi = cos (Line_Phi);
   sin_phi = sin (Line_Phi);

   alpha = JDMv_find_rotation_axis (JDMv_vector(-1,0,0), st->p, &normal);

   for (i = 0; i < num; i++)
     {
	double sin_theta;
	JDMVector_Type p;

	if (-1 == (*efun) (&st->spectrum, &at->energy))
	  return -1;

	theta = Line_Theta * (-1.0 + 2.0 * JDMrandom ());

	sin_theta = -sin (theta);
	p.x = -cos (theta);
	p.y = sin_theta * cos_phi;
	p.z = sin_theta * sin_phi;

	at->p = JDMv_rotate_unit_vector (p, normal, alpha);
	at++;
     }

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_line_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   (void) source_id;
   st->open_source = line_open_source;
   st->create_photons = line_create_photons;
   st->close_source = line_close_source;

   if ((-1 == pf_get_double (p, "S-LinePhi", &Line_Phi))
       || (-1 == pf_get_double (p, "S-LineTheta", &Line_Theta)))
     return -1;

   /* Convert to radians */
   Line_Theta = Line_Theta * (1.0 / 3600.0 *  PI / 180.0);    /* 1 arc sec */
   Line_Phi = Line_Phi * (PI / 180.0);    /* degrees --> radians */

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

