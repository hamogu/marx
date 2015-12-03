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
#include <jdmath.h>

#include "acis.h"
#include "marx.h"
#include "_marx.h"

static Marx_Detector_Type ACIS_I_Detector;
static Marx_Detector_Type ACIS_S_Detector;

static Marx_Detector_Geometry_Type ACIS_S_Geom [_MARX_NUM_ACIS_S_CHIPS];
static double ACIS_S_TDet_Xoffs[_MARX_NUM_ACIS_S_CHIPS] =
  {791, 1833, 2875, 3917, 4959, 6001};
static double ACIS_S_TDet_Yoffs[_MARX_NUM_ACIS_S_CHIPS] =
  {1702, 1702, 1702, 1702, 1702, 1702};
static double ACIS_S_Pixel_Size = 0.024;   /* mm */

static Marx_Detector_Geometry_Type ACIS_I_Geom [_MARX_NUM_ACIS_I_CHIPS];
static double ACIS_I_TDet_Xoffs[_MARX_NUM_ACIS_I_CHIPS] =
  {3061, 5131, 3061, 5131};
static double ACIS_I_TDet_Yoffs[_MARX_NUM_ACIS_I_CHIPS] =
  {5131, 4107, 4085, 3061};
static double ACIS_I_Pixel_Size = 0.024;   /* mm */

/* This is necessary because the corners defined in the pixlib parameter
 * file do not produce a detector whose size is 1024 * 0.024.
 */
static void tweak_acis_pixel_size (Marx_Detector_Type *det)
{
   Marx_Detector_Geometry_Type *d;

   d = det->facet_list;

   while (d != NULL)
     {
	d->x_pixel_size = d->xlen / 1024.0;
	d->y_pixel_size = d->ylen / 1024.0;

	d = d->next;
     }
}

static int print_info (Marx_Detector_Type *det, FILE *fp)
{
   (void) fprintf (fp, "STT-LSI offset: (% 10.4e, % 10.4e, % 10.4e)\n",
		   det->stt_lsi_offset.x,
		   det->stt_lsi_offset.y,
		   det->stt_lsi_offset.z);
   (void) fprintf (fp, "STF-STT offset: (% 10.4e, % 10.4e, % 10.4e)\n",
		   det->stf_stt_offset.x,
		   det->stf_stt_offset.y,
		   det->stf_stt_offset.z);
   return 0;
}

static int post_init_acis_detector (Marx_Detector_Type *d)
{
   if (-1 == _marx_caldb_patch_acis_geom (d))
     return -1;

   if (-1 == _marx_caldb_patch_aimpoint (d))
     return -1;

   if (-1 == _marx_compute_detector_basis (d))
     return -1;

   tweak_acis_pixel_size (d);

   if (NULL == (d->fp_coord_info = marx_get_fp_system_info (d->fp_system_name)))
     return -1;

   d->print_info = print_info;

   d->is_initialized = 1;
   return 0;
}

/* 01
 * 23
 * Note that the (0,0) pixel is not at the lower left of all these chips.
 * For chips 0 and 2, it is at the upper left with the xpixel running down
 * and the y pixel running to the right.
 */

static int
acis_i_to_tiled (Marx_Detector_Type *det,
		 Marx_Detector_Geometry_Type *g,
		 int chip, unsigned int x, unsigned int y,
		 unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) det;
   switch (chip)
     {
      case 0:
      case 2:
	xf = g->tdet_xoff + y;
	yf = g->tdet_yoff - x;
	break;

      case 1:
      case 3:
      default:
	xf = g->tdet_xoff - y;
	yf = g->tdet_yoff + x;
	break;
     }

   if (xf < 0.0) xf = 0.0;
   if (yf < 0.0) yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;

   return 0;
}

Marx_Detector_Type *_marx_get_acis_i_detector (void)
{
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *g;
   int ccd_id;

   d = &ACIS_I_Detector;
   if (d->is_initialized)
     return d;

   d->detector_type = MARX_DETECTOR_ACIS_I;
   d->tiled_pixel_map_fun = &acis_i_to_tiled;
   d->facet_list = _marx_link_detector_facet_list (ACIS_I_Geom, _MARX_NUM_ACIS_I_CHIPS, sizeof(Marx_Detector_Geometry_Type));
   d->fp_system_name = "AXAF-FP-1.1";
   d->first_facet_id = 0;
   d->last_facet_id = 3;

   g = d->facet_list;
   for (ccd_id = 0; ccd_id <= 3; ccd_id++)
     {
	g->id = ccd_id;
	g->tdet_xoff = ACIS_I_TDet_Xoffs[ccd_id];
	g->tdet_yoff = ACIS_I_TDet_Yoffs[ccd_id];
	g->x_pixel_size = ACIS_I_Pixel_Size;
	g->y_pixel_size = ACIS_I_Pixel_Size;
	g->num_x_pixels = 1024;
	g->num_y_pixels = 1024;
	g++;
     }

   if (-1 == post_init_acis_detector (d))
     return NULL;

   return d;
}

static int acis_s_to_tiled (Marx_Detector_Type *det,
			    Marx_Detector_Geometry_Type *g,
			    int chip, unsigned int x, unsigned int y,
			    unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) det;
   (void) chip;

   xf = x + g->tdet_xoff;
   if (xf < 0.0)
     xf = 0.0;

   yf = y + g->tdet_yoff;
   if (yf < 0.0)
     yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;
   return 0;
}

Marx_Detector_Type *_marx_get_acis_s_detector (void)
{
   int ccd_id;
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *g;

   d = &ACIS_S_Detector;
   if (d->is_initialized)
     return d;

   d->detector_type = MARX_DETECTOR_ACIS_S;
   d->tiled_pixel_map_fun = &acis_s_to_tiled;
   d->facet_list = _marx_link_detector_facet_list (ACIS_S_Geom, _MARX_NUM_ACIS_S_CHIPS, sizeof(Marx_Detector_Geometry_Type));
   d->fp_system_name = "AXAF-FP-1.1";
   d->first_facet_id = 4;
   d->last_facet_id = 9;

   g = d->facet_list;
   for (ccd_id = 4; ccd_id <= 9; ccd_id++)
     {
	g->id = ccd_id;
	g->tdet_xoff = ACIS_S_TDet_Xoffs[ccd_id-4];
	g->tdet_yoff = ACIS_S_TDet_Yoffs[ccd_id-4];
	g->x_pixel_size = ACIS_S_Pixel_Size;
	g->y_pixel_size = ACIS_S_Pixel_Size;
	g->num_x_pixels = 1024;
	g->num_y_pixels = 1024;
	g++;
     }

   if (-1 == post_init_acis_detector (d))
     return NULL;

   return d;
}

