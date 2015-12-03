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
#include <string.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <jdmath.h>

#include "chandra.h"
#include "marx.h"
#include "_marx.h"

double Marx_Focal_Length = 10065.5;    /* mm */

static Marx_FP_Coord_Type Focal_Plane_Coord_Info [] =
{
     {
	"AXAF-FP-1.1",		       /* fp_system_name */
	0.492 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
	4096.5,			       /* fp x pixel of on axis ray */
	4096.5			       /* fp y pixel of on axis ray */
     },
     {
	"AXAF-FP-2.1",		       /* fp_system_name */
	0.13175 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
	16384.5,		       /* fp x pixel of on axis ray */
	16384.5			       /* fp y pixel of on axis ray */
     },
     {
	"AXAF-FP-2.3",		       /* fp_system_name */
	0.1318 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
	32768.5,		       /* fp x pixel of on axis ray */
	32768.5			       /* fp y pixel of on axis ray */
     },
     {
	NULL, 0, 0, 0
     }
};

Marx_FP_Coord_Type *marx_get_fp_system_info (char *name)
{
   Marx_FP_Coord_Type *fp;

   fp = Focal_Plane_Coord_Info;
   while (fp->fp_system_name != NULL)
     {
	if (0 == strcmp (name, fp->fp_system_name))
	  return fp;
	fp++;
     }

   marx_error ("Focal plane system %s not implemented", name);
   return NULL;
}

typedef struct
{
   char *name;
   int id;
}
Det_Name_Type;

static Det_Name_Type Detectors [] =
{
   {"HRC-S", MARX_DETECTOR_HRC_S},
   {"HRC-I", MARX_DETECTOR_HRC_I},
   {"ACIS-S", MARX_DETECTOR_ACIS_S},
   {"ACIS-I", MARX_DETECTOR_ACIS_I},
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_DETECTOR_IXO_CATGS_CCD
   {"IXOCCD", MARX_DETECTOR_IXO_CATGS_CCD},
# endif
# ifdef MARX_DETECTOR_IXO_XMS
   {"IXOXMS", MARX_DETECTOR_IXO_XMS},
# endif
#endif
   {NULL, 0}
};

Marx_Detector_Type *
marx_get_detector_info (char *detname) /*{{{*/
{
   Marx_Detector_Type *d;
   Det_Name_Type *dn;
   int det;

   if (detname == NULL)
     return NULL;

   dn = Detectors;
   det = -1;
   while (dn->name != NULL)
     {
	if (0 == strcmp (detname, dn->name))
	  {
	     det = dn->id;
	     break;
	  }
	dn++;
     }

   switch (det)
     {
      default:
	marx_error ("Detector %s is unknown", detname);
	return NULL;

      case MARX_DETECTOR_ACIS_S:
	d = _marx_get_acis_s_detector ();
	break;

      case MARX_DETECTOR_ACIS_I:
	d = _marx_get_acis_i_detector ();
	break;

      case MARX_DETECTOR_HRC_S:
	d = _marx_get_hrc_s_detector ();
	break;

      case MARX_DETECTOR_HRC_I:
	d = _marx_get_hrc_i_detector ();
	break;
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_DETECTOR_IXO_CATGS_CCD
      case MARX_DETECTOR_IXO_CATGS_CCD:
	d = _marx_get_ixo_ccd_detector ();
	break;
# endif
# ifdef MARX_DETECTOR_IXO_XMS
      case MARX_DETECTOR_IXO_XMS:
	d = _marx_get_ixo_xms_detector ();
	break;
# endif
#endif
     }

   return d;
}

/*}}}*/

int
marx_compute_tiled_pixel (Marx_Detector_Type *d, /*{{{*/
			      int chip, unsigned int x, unsigned int y,
			      unsigned int *xp, unsigned int *yp)
{
   Marx_Detector_Geometry_Type *g;

   if (d == NULL)
     return -1;

   if (d->tiled_pixel_map_fun == NULL)
     {
	marx_error ("marx_compute_tiled_pixel: detector not supported");
	return -1;
     }

   g = d->facet_list;

   while (g != NULL)
     {
	if (g->id == chip)
	  return (*d->tiled_pixel_map_fun) (d, g, chip, x, y, xp, yp);

	g = g->next;
     }

   marx_error ("chip = %d is not appropriate for this detector", chip);
   return -1;
}

/*}}}*/

int marx_mnc_to_fpc (Marx_FP_Coord_Type *fc,
		     JDMVector_Type *mnc,
		     double *xp, double *yp)
{
   double mnc_x;
   double factor;

   if (fc == NULL)
     return -1;

   mnc_x = mnc->x;
   if (mnc_x == 0.0)
     {
	marx_error ("marx_mnc_to_fpc: mnc.x is 0");
	return -1;
     }

   factor = 1.0 / (fc->fp_delta_s0 * mnc_x);

   /* Consider a ray with a negative MNC->x coordinate, but positive y and
    * z coordinates.  We want this to map to a positive FP_X and a negative
    * FP_Y value.
    */
   *xp = fc->fp_x0 - factor * mnc->y;
   *yp = fc->fp_y0 + factor * mnc->z;

   return 0;
}

