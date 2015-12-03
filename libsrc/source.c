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
/* Source controller module. */

#include "config.h"
#include "marx-feat.h"

#include <stdio.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <string.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"

typedef struct /*{{{*/
{
   char *name;
   int (*select_source) (Marx_Source_Type *, Param_File_Type *, char *, unsigned int);
   int dither_flags;
   int source_id;
}

/*}}}*/
Source_Object_Type;

#define MARX_SOURCE_C_FILE 1
#include "source.def"

/* If the source is at infinity, Source_Distance is <= 0.
 */
static double Source_Distance = 0.0;
static double Source_Ra, Source_Dec;
static double Source_Azimuth, Source_Elevation;

static char *Source_Name;

static Param_Table_Type Source_Parm_Table [] = /*{{{*/
{
   {"SourceDistance",		PF_REAL_TYPE,	&Source_Distance},
   {"SourceRA",			PF_REAL_TYPE,	&Source_Ra},
   {"SourceDec",		PF_REAL_TYPE,	&Source_Dec},
#if 0
   {"SourceOffsetZ",		PF_REAL_TYPE,	&Source_Elevation},
   {"SourceOffsetY",		PF_REAL_TYPE,	&Source_Azimuth},
#endif
   {"SourceType",		PF_STRING_TYPE,	&Source_Name},
   {NULL, 			0,		NULL}
};

/*}}}*/

static Source_Object_Type *find_source (char *name)
{
   Source_Object_Type *s = Sources;
   int source_id = 0;

   while (s->name != NULL)
     {
	if (!strcmp (name, s->name))
	  {
	     s->source_id = source_id;
	     return s;
	  }
	source_id++;
	s++;
     }
   marx_error ("Unable to locate source \"%s\".", name);
   return NULL;
}

/* This function computes the azimuth and elevation of the source with
 * respect to the spacecraft pointing axis, which is generally different
 * from the optical axis.
 */
extern int marx_compute_elaz (double Ra, double Dec, double *azp, double *elp)
{
   JDMVector_Type src, pnt;
   double ra_pnt, dec_pnt, roll_pnt;
   double az, el;

   src = JDMv_spherical_to_vector (1.0, PI/2.0 - Dec, Ra);
   if (-1 == marx_get_pointing (&ra_pnt, &dec_pnt, &roll_pnt))
     return -1;
   pnt = JDMv_spherical_to_vector (1.0, PI/2.0 - dec_pnt, ra_pnt);
   src = JDMv_rotate_unit_vector (src, pnt, -roll_pnt);
   JDMv_unit_vector_to_spherical (src, &el, &az);
   el = PI/2.0 - el;

   marx_compute_ra_dec_offsets (ra_pnt, dec_pnt, az, el, &az, &el);

   *elp = el;
   *azp = az;
   return 0;
}

static int select_source (Marx_Source_Type *st, Param_File_Type *pf, char *name) /*{{{*/
{
   Source_Object_Type *s;
   double yoff, zoff;
   double ra_nom, dec_nom, roll_nom;
   JDMVector_Type p;

   s = find_source (name);
   if (s == NULL)
     return -1;

   if (-1 == _marx_init_dither (pf, s->dither_flags, &yoff, &zoff))
     return -1;

   if (-1 == marx_get_nominal_pointing (&ra_nom, &dec_nom, &roll_nom))
     return -1;

   if (-1 == marx_compute_elaz (Source_Ra, Source_Dec, &Source_Azimuth, &Source_Elevation))
     return -1;

   p = JDMv_spherical_to_vector (1.0, 0.5*PI-zoff, yoff);

   // SIMPUT source gives absolute RA, DEC
   // MARX sources give relative to source position
   if (strcmp(name, "SIMPUT") != 0){
     /* Now add offsets via the proper rotations */
     p = JDMv_rotate_unit_vector (p, JDMv_vector (0, -1, 0), Source_Elevation);
     p = JDMv_rotate_unit_vector (p, JDMv_vector (0, 0, 1), Source_Azimuth);

     /* Finally roll it so that this point will be invariant under roll.  That is,
      * the dither transformation will (on the average) undo this rotation.
      * See the apply_dither function.
      */
     p = JDMv_rotate_unit_vector (p, JDMv_vector (1, 0, 0), roll_nom);

     /* This vector must point FROM source TO origin. */
     st->p.x = -p.x;
     st->p.y = -p.y;
     st->p.z = -p.z;

     /* Create a vector orthogonal to above.  This will save the source
      * routine the effort required to do this.
      *
      * Since st->p is more or less oriented along the negative x direction,
      * st->p.x will be non-zero.  In fact, this will be required.
      * With this requirement, a normal in the x-y plane is trival to construct.
      */

     if (p.x <= 0.0)
       {
	 marx_error ("Source rays will not hit the telescope.");
	 return -1;
       }

     st->p_normal.z = 0.0;
     st->p_normal.y = 1.0;
     st->p_normal.x = -p.y / p.x;
     JDMv_normalize (&st->p_normal);
   }
   else{
     st->p = p;
   }

   st->distance = Source_Distance;

   marx_message ("Initializing source type %s...\n", name);

   if (-1 == (*s->select_source)(st, pf, name, s->source_id))
     return -1;

   return 0;
}

/*}}}*/

int marx_close_source (Marx_Source_Type *st) /*{{{*/
{
   if (st == NULL)
     return -1;

#if MARX_HAS_DITHER
   if (_Marx_Dither_Mode != _MARX_DITHER_MODE_NONE)
     _marx_close_dither ();
#endif

   if (st->spectrum.close_spectrum != NULL)
     st->spectrum.close_spectrum (&st->spectrum);

   if (st->close_source != NULL)
     {
	(void) (*st->close_source)(st);
     }

   memset ((char *)st, 0, sizeof (Marx_Source_Type));
   marx_free ((char *)st);
   return 0;
}

/*}}}*/

int marx_open_source (Marx_Source_Type *st) /*{{{*/
{
   if (st == NULL)
     return -1;

   if (st->open_source == NULL)
     return -1;

   if (-1 == (*st->open_source) (st))
     return -1;

   return 0;
}

/*}}}*/

static double compute_mean_time (double flux) /*{{{*/
{
   if (flux <= 0.0) return 0.0;
   return 1.0 / flux / Marx_Mirror_Geometric_Area;
}

/*}}}*/

int marx_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
			 unsigned int num, unsigned int *num_collected,
			 double *exposure_time)
{
   unsigned int n, i;
   double *sorted_energies;
   unsigned int *sorted_index;
   Marx_Photon_Attr_Type *at;

   *num_collected = 0;

   if ((st == NULL) || (pt == NULL) || (st->create_photons == NULL))
     return -1;

   pt->source_distance = st->distance;
   pt->history = 0;

   pt->start_time += pt->total_time;

   memset ((char *) pt->attributes, 0, num * sizeof (Marx_Photon_Attr_Type));

   if (-1 == (*st->create_photons) (st, pt, num, &n))
     return -1;

   sorted_energies = pt->sorted_energies;
   at = pt->attributes;

   if (pt->history & MARX_TIME_OK)
     {
	if (exposure_time != NULL)
	  {
	     double tmax = _Marx_TStart_MJDsecs + (*exposure_time);
	     for (i = 0; i < n; i++)
	       {
		  sorted_energies[i] = at[i].energy;
		  if (at[i].arrival_time >= tmax)
		    {
		       n = i + 1;
		       break;
		    }
	       }
	  }
	else for (i = 0; i < n; i++)
	  sorted_energies[i] = at[i].energy;
     }
   else
     {
	double t, mt;

	mt = compute_mean_time (st->spectrum.total_flux);
	t = 0.0;

	for (i = 0; i < n; i++)
	  {
	     sorted_energies[i] = at->energy;
	     /* For independent poisson events, the distribution of times
	      * between events is exponentially distributed.
	      */
	     t += mt * JDMexpn_random ();

	     at->arrival_time = t;
	     at->flags = 0;
	     if ((exposure_time != NULL) && (*exposure_time <= t))
	       {
		  n = i + 1;
		  break;
	       }
	     at++;
	  }

	pt->history |= (MARX_ENERGY_OK
			| MARX_TIME_OK
			| MARX_X_VECTOR_OK
			| MARX_P_VECTOR_OK);
     }

   if (n)
     {
	unsigned int tag_start = pt->tag_start;
	at = pt->attributes;
	pt->num_sorted = pt->n_photons = n;
	for (i = 0; i < n; i++)
	  {
	     at[i].tag = tag_start + i;
	  }
	pt->history |= MARX_TAG_OK;
	pt->tag_start = tag_start + n;
	pt->total_time = at[n - 1].arrival_time;
#if MARX_HAS_DITHER
	if (-1 == _marx_dither_photons (pt, &n))
	  return -1;
#endif
     }

   /* Now sort the energies for later reference */
   if (pt->sorted_index != NULL)
     JDMfree_integer_vector ((int *) pt->sorted_index);

   sorted_index = pt->sorted_index = JDMsort_doubles (sorted_energies, n);
   if (sorted_index == NULL)
     return -1;

   at = pt->attributes;
   /* re-arrange according to index */
   for (i = 0; i < n; i++)
     {
	sorted_energies[i] = at[sorted_index[i]].energy;
     }

   pt->num_sorted = pt->n_photons = n;

   if (n == 0) pt->total_time = 0.0;
   else pt->total_time = at[n - 1].arrival_time;
   *num_collected = n;

   return 0;
}

/*}}}*/

Marx_Source_Type *marx_create_source (Param_File_Type *p) /*{{{*/
{
   Marx_Source_Type *st;

   if (-1 == pf_get_parameters (p, Source_Parm_Table))
     return NULL;

   st = (Marx_Source_Type *) marx_malloc (sizeof (Marx_Source_Type));
   if (st == NULL) return st;

   memset ((char *) st, 0, sizeof (Marx_Source_Type));

   /* Scale to mm */
   Source_Distance = Source_Distance * 1000.0;

   /* Convert from minutes to radians */
   Source_Azimuth = Source_Azimuth * (PI / 180.0 / 60.0);
   Source_Elevation = Source_Elevation * (PI / 180.0 / 60.0);

   /* Convert from degrees to radians */
   Source_Ra *= (PI / 180.0);
   Source_Dec *= (PI / 180.0);

   if (-1 == select_source (st, p, Source_Name))
     {
	marx_free ((char *) st);
	return NULL;
     }
   return st;
}

/*}}}*/

