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

#include "hrc.h"
#include "marx.h"
#include "_marx.h"

static double X_Pixel_Size;
static double Y_Pixel_Size;
static double LL_XYZ[3];
static double Center_XYZ[3];
static double LL_to_LR_Vector[3];
static double Center_CXCY[2];
static double Aimpoint_CXCY[2];

/* These quantities are computed */
static JDMVector_Type LL_Vector;
static JDMVector_Type LR_Vector;
static JDMVector_Type UL_Vector;
static JDMVector_Type UR_Vector;
static double LL_CXCY[2];

static int setup_coordinate_xforms (void)
{
   double len;
   JDMVector_Type center, ll, lr, a, b, dx;

   a = JDMv_vector (LL_to_LR_Vector[0],
		    LL_to_LR_Vector[1],
		    LL_to_LR_Vector[2]);

   len = JDMv_length (a);
   if (len == 0.0)
     {
	marx_error ("Length of vector from LL to LR must be non-zero");
	return -1;
     }
   a = JDMv_smul (1.0/len, a);

   center = JDMv_vector (Center_XYZ[0], Center_XYZ[1], Center_XYZ[2]);

   LL_Vector = JDMv_vector (LL_XYZ[0], LL_XYZ[1], LL_XYZ[2]);
   ll = JDMv_diff (LL_Vector, center);

   if (0.0 == JDMv_length (ll))
     {
	marx_error ("The LL corner must not be the center");
	return -1;
     }

   lr = JDMv_ax1_bx2 (1.0, ll,
		      -2.0 * JDMv_dot_prod (a, ll), a);

   LR_Vector = JDMv_sum (center, lr);
   UR_Vector = JDMv_sum (center, JDMv_smul (-1.0, ll));
   UL_Vector = JDMv_sum (center, JDMv_smul (-1.0, lr));

   /* Now determine the pixel coord of the LL corner */

   if ((X_Pixel_Size <= 0.0) || (Y_Pixel_Size <= 0.0))
     {
	marx_error ("X and Y Pixel sizes must be larger than 0");
	return -1;
     }

   b = JDMv_unit_vector (JDMv_diff (UL_Vector, LL_Vector));

   LL_CXCY[0] = Center_CXCY[0] + JDMv_dot_prod (ll, a)/X_Pixel_Size;
   LL_CXCY[1] = Center_CXCY[1] + JDMv_dot_prod (ll, b)/Y_Pixel_Size;

   /* Finally we need to translate the system so that the aimpoint is
    * at the origin.
    */

   dx = JDMv_diff (JDMv_ax1_bx2 ((LL_CXCY[0] - Aimpoint_CXCY[0])*X_Pixel_Size, a,
				 (LL_CXCY[1] - Aimpoint_CXCY[1])*Y_Pixel_Size, b),
		   LL_Vector);

   LL_Vector = JDMv_sum (LL_Vector, dx);
   LR_Vector = JDMv_sum (LR_Vector, dx);
   UR_Vector = JDMv_sum (UR_Vector, dx);
   UL_Vector = JDMv_sum (UL_Vector, dx);

   return 0;
}

static int patch_hrc_i_geom (Marx_Detector_Type *d)
{
   Marx_Detector_Geometry_Type *g;

   if (-1 == _marx_hrc_i_geom_init (NULL))
     return -1;

   /* d->y_pixel_size = Y_Pixel_Size; */
   /* d->x_pixel_size = X_Pixel_Size; */

   g = d->facet_list;

   g->id = 0;

   g->x_ll = LL_Vector;
   g->x_ul = UL_Vector;
   g->x_ur = UR_Vector;
   g->x_lr = LR_Vector;

   /* For AXAF-HRC-2.4I */

   g->tdet_xoff = -1.0;
   g->tdet_yoff = 0.0;

   g->xpixel_offset = LL_CXCY[0];
   g->ypixel_offset = LL_CXCY[1];

   return 0;
}

/* dx, dy are measured from the LL corner */
int _marx_hrc_i_compute_pixel (double dx, double dy,
			       double *xpixel, double *ypixel)
{
   dx /= X_Pixel_Size;
   dy /= Y_Pixel_Size;

   *xpixel = LL_CXCY[0] + dx;
   *ypixel = LL_CXCY[1] + dy;

   return 0;
}

int _marx_hrc_i_get_pixel_size (double *dx, double *dy)
{
   *dx = X_Pixel_Size;
   *dy = Y_Pixel_Size;
   return 0;
}

static _Marx_Simple_Data_Type Array_Data_Table [] =
{
   {"HRC_I_X_Pixel_Size",	1, &X_Pixel_Size,	1.0, 0},
   {"HRC_I_Y_Pixel_Size",	1, &Y_Pixel_Size,	1.0, 0},

   {"HRC_I_LL_XYZ",		3, LL_XYZ,		1.0, 0},
   {"HRC_I_Center_XYZ",		3, Center_XYZ,		1.0, 0},
   {"HRC_I_LL_to_LR_Vector",	3, LL_to_LR_Vector,	1.0, 0},

   {"HRC_I_Center_CXCY",	2, Center_CXCY,		1.0, 0},
   {"HRC_I_Aimpoint_CXCY",	2, Aimpoint_CXCY,	1.0, 0},

   {NULL, 0, NULL, 0.0, 0}
};

/* Note: This may be called with pf == NULL */
int _marx_hrc_i_geom_init (Param_File_Type *pf)
{
   char *file;
   static int initialized = 0;

   (void) pf;

   if (initialized)
     return 0;

   file = "hrc/hrc_i_geom.txt";

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;

   marx_message ("\t%s\n", file);

   if (-1 == _marx_read_simple_data_file (file, Array_Data_Table))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   if (-1 == setup_coordinate_xforms ())
     return -1;

   initialized = 1;
   return 0;
}

static int
hrc_i_to_tiled (Marx_Detector_Type *det,
		Marx_Detector_Geometry_Type *g,
		int chip, unsigned int x, unsigned int y,
		unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) chip;
   (void) det;
   (void) g;

   xf = (x + g->tdet_xoff);
   yf = (y + g->tdet_yoff);

   if (xf < 0.0) xf = 0.0;
   if (yf < 0.0) yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;
   return 0;
}

static Marx_Detector_Type HRC_I_Detector;
static Marx_Detector_Geometry_Type HRC_I_Geom[_MARX_NUM_HRC_I_CHIPS];

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

Marx_Detector_Type *_marx_get_hrc_i_detector (void)
{
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *g;

   d = &HRC_I_Detector;
   if (d->is_initialized)
     return d;

   d->detector_type = MARX_DETECTOR_HRC_I;
   d->tiled_pixel_map_fun = &hrc_i_to_tiled;
   d->facet_list = _marx_link_detector_facet_list (HRC_I_Geom, _MARX_NUM_HRC_I_CHIPS, sizeof(Marx_Detector_Geometry_Type));
   d->fp_system_name = "AXAF-FP-2.1";
   d->first_facet_id = 0;
   d->last_facet_id = 0;
   d->print_info = print_info;

   g = d->facet_list;
   g->id = 0;
   g->tdet_xoff = 0;
   g->tdet_yoff = 0;
   g->num_x_pixels = 16384;
   g->num_y_pixels = 16384;
   g->x_pixel_size = 6.429e-3;
   g->y_pixel_size = 6.429e-3;

   if (-1 == patch_hrc_i_geom (d))
     return NULL;

   if (-1 == _marx_caldb_patch_aimpoint (d))
     return NULL;

   if (-1 == _marx_compute_detector_basis (d))
     return NULL;

   if (NULL == (d->fp_coord_info = marx_get_fp_system_info (d->fp_system_name)))
     return NULL;

   d->is_initialized = 1;

   return d;
}
