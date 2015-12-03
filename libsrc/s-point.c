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

static int point_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int point_close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int point_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				 unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   Marx_Photon_Attr_Type *at;
   int (*efun) (Marx_Spectrum_Type *, double *);

   at = pt->attributes;
   efun = st->spectrum.energy_function;

   for (i = 0; i < num; i++)
     {
	if (-1 == (*efun) (&st->spectrum, &at->energy))
	  return -1;

	/* For a point source, the calculation of the vector to the origin
	 * is trivial. */
	at->p = st->p;

	at++;
     }

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_point_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   (void) source_id;

   st->open_source = point_open_source;
   st->create_photons = point_create_photons;
   st->close_source = point_close_source;

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

