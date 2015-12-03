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

#define USE_CALDB_FILES	0	       /* doesn't work yet */
/* All these parameters get set from a file.  See below */

static double V_Active [6] =
{
   1665.0, 16404.0, 17060.0, 32165.0, 32975.0, 47605.0
};
static double V_Physical [6] =
{
   1196, 16444, 16900, 32250, 32930, 48117
};

static double U_Active [2] =
{
   560, 3536
};

static double U_Physical [2] =
{
   0, 4096
};

static double U_Pixel_Size = 0.006250;
static double V_Pixel_Size = 0.0064294;

/* The center of the detector is offset from the aimpoint by 4mm */
static double Aimpoint_Offset [3] =
{
   0, 4, 0
};

static double Left_Theta = 1.4321 * PI/180.0;
static double Right_Theta = 1.2202 * PI/180.0;

/* Fiducial Points on each MCP to define (u,v) <---> (x,y,z) mapping */
static double Left_XYZ [3];
static double Left_UV[2];
static double Middle_XYZ[3];
static double Middle_UV[2];
static double Right_XYZ[3];
static double Right_UV[2];
static double Left_CXCY[2];
static double Right_CXCY[2];
static double Middle_CXCY[2];

static double Detector_Theta;

/* For future changes, allow the most general linear transformation between
 * (x,y,z) and (u, v) coordinates.  Basically, this involves:
 *
 *     X_i = X0_i + (u - u0) alpha_i + (v - v0) beta_i
 *
 * Here alpha_i and beta_i are scale factors.  Furthermore, we assume that
 * that the (u,v) system is orthogonal.  Then consider a small change in the
 * u coordinate holding v fixed:
 *
 *    d(u)X_i = du alpha_i
 *
 * and
 *
 *    d(v)X_i = dv beta_i
 *
 * By orthogonal, we mean that the resulting vectors d(u)X_i and d(v)X_i are
 * orthogonal, i.e.,
 *
 *    0 = d(u)X_i d(v)X_i ==>  alpha_i beta_i = 0
 *
 * where repeated indices are summed over.  Thus, alpha and beta are orthogonal.
 * From this result, it trivally follows that:
 *
 *   u = u_0 + alpha_inv_i (X_i - X0_i)
 *   v = v_0 + beta_inv_i (X_i - X0_i)
 *
 * where alpha_inv_i = alpha_i / |alpha|^2
 */

static double Left_Alpha[3], Left_Alpha_Inv[3];
static double Left_Beta[3], Left_Beta_Inv[3];

static double Middle_Alpha[3], Middle_Alpha_Inv[3];
static double Middle_Beta[3], Middle_Beta_Inv[3];

static double Right_Alpha[3], Right_Alpha_Inv[3];
static double Right_Beta[3], Right_Beta_Inv[3];

static void map_uv_to_xyz (double *uv, double *uv0,
			   double *alpha, double *beta, double *xyz0,
			   double *xyz)
{
   unsigned int i;
   double du, dv;

   du = uv[0] - uv0[0];
   dv = uv[1] - uv0[1];

   for (i = 0; i < 3; i++)
     xyz[i] = xyz0[i] + du * alpha[i] + dv * beta[i];
}
#if 0
static void  map_xyz_to_uv (double *xyz, double *xyz0,
			    double *alpha_inv, double *beta_inv,
			    double *uv0, double *uv)
{
   double u, v;
   unsigned int i;

   u = uv0[0];
   v = uv0[1];
   for (i = 0; i < 3; i++)
     {
	double dx = xyz[i] - xyz0[i];
	u += alpha_inv[i] * dx;
	v += beta_inv[i] * dx;
     }
   uv[0] = u;
   uv[1] = v;
}
#endif
static void compute_inverse_vector (double *a, double *a_inv)
{
   double x;

   x = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];

   a_inv[0] = a[0]/x;
   a_inv[1] = a[1]/x;
   a_inv[2] = a[2]/x;
}

static void setup_alphas_and_betas (void)
{
   /* On the left side, increasing v implies decreases in x and y */
   /* HRC-S1 */
   Left_Alpha[0] = 0.0;
   Left_Alpha[1] = 0.0;
   Left_Alpha[2] = U_Pixel_Size;
   compute_inverse_vector (Left_Alpha, Left_Alpha_Inv);
   Left_Beta[0] = -1.0 * sin(Left_Theta) * V_Pixel_Size;
   Left_Beta[1] = -1.0 * cos(Left_Theta) * V_Pixel_Size;
   Left_Beta[2] = 0.0;
   compute_inverse_vector (Left_Beta, Left_Beta_Inv);

   /* On the middle, increasing v implies decrease in y */
   /* HRC-S2 */
   Middle_Alpha[0] = 0.0;
   Middle_Alpha[1] = 0.0;
   Middle_Alpha[2] = U_Pixel_Size;
   compute_inverse_vector (Middle_Alpha, Middle_Alpha_Inv);
   Middle_Beta[0] = 0.0;
   Middle_Beta[1] = -1.0 * V_Pixel_Size;
   Middle_Beta[2] = 0.0;
   compute_inverse_vector (Middle_Beta, Middle_Beta_Inv);

   /* On the right side, increasing v implies decrease in y and increase in x */
   /* HRC-S3 */
   Right_Alpha[0] = 0.0;
   Right_Alpha[1] = 0.0;
   Right_Alpha[2] = U_Pixel_Size;
   compute_inverse_vector (Right_Alpha, Right_Alpha_Inv);
   Right_Beta[0] = 1.0 * sin(Right_Theta) * V_Pixel_Size;
   Right_Beta[1] = -1.0 * cos(Right_Theta) * V_Pixel_Size;
   Right_Beta[2] = 0.0;
   compute_inverse_vector (Right_Beta, Right_Beta_Inv);
}

static JDMVector_Type uv_to_xyz_vector (double u, double v, double *uv0,
					double *alpha, double *beta, double *xyz0)
{
   double xyz[3];
   double uv[2];
   JDMVector_Type xyz_vector;

   uv[0] = u; uv[1] = v;

   map_uv_to_xyz (uv, uv0, alpha, beta, xyz0, xyz);

   xyz_vector.x = xyz[0];
   xyz_vector.y = xyz[1];
   xyz_vector.z = xyz[2];

   return xyz_vector;
}

static JDMVector_Type left_uv_to_xyz_vector (double u, double v)
{
   return uv_to_xyz_vector (u, v, Left_UV, Left_Alpha, Left_Beta, Left_XYZ);
}

static JDMVector_Type right_uv_to_xyz_vector (double u, double v)
{
   return uv_to_xyz_vector (u, v, Right_UV, Right_Alpha, Right_Beta, Right_XYZ);
}

static JDMVector_Type middle_uv_to_xyz_vector (double u, double v)
{
   return uv_to_xyz_vector (u, v, Middle_UV, Middle_Alpha, Middle_Beta, Middle_XYZ);
}

static void rotate_mcp (Marx_Detector_Geometry_Type *g, double theta)
{
   JDMVector_Type axis;
   double sin_theta, cos_theta;

   axis = JDMv_vector (-1, 0, 0);      /* choose axis from origin along optical axis towards the HRMA */
   sin_theta = sin(theta);
   cos_theta = cos(theta);
   g->x_ll = JDMv_rotate_vector1 (g->x_ll, axis, cos_theta, sin_theta);
   g->x_lr = JDMv_rotate_vector1 (g->x_lr, axis, cos_theta, sin_theta);
   g->x_ul = JDMv_rotate_vector1 (g->x_ul, axis, cos_theta, sin_theta);
   g->x_ur = JDMv_rotate_vector1 (g->x_ur, axis, cos_theta, sin_theta);
}

static int patch_hrc_s_geom (Marx_Detector_Type *d)
{
   Marx_Detector_Geometry_Type *g;
   unsigned int i;

   if (-1 == _marx_hrc_s_geom_init (NULL))
     return -1;

   /* d->y_pixel_size = V_Pixel_Size; */
   /* d->x_pixel_size = U_Pixel_Size; */

   /* Ok to index by arrays in this file */
   g = d->facet_list;

   /* The HRC-S corners are designated as follows:
    *
    *     UR   LR
    *     UL   LL
    *
    * whereas ACIS-S is more sensible:
    *
    *     UL  UR
    *     LL  LR
    *
    * Sigh.
    */

#if USE_CALDB_FILES
   g[0].id = 1;
   g[1].id = 2;
   g[2].id = 3;
   if (-1 == _marx_caldb_patch_hrc_s_geom (d))
     return -1;
#else
   /* Handle LEFT Segment */
   /* HRC-S1 */
   g[0].id = 1;
   g[0].x_ul = left_uv_to_xyz_vector (U_Active[0], V_Active[1]);
   g[0].x_ur = left_uv_to_xyz_vector (U_Active[1], V_Active[1]);
   g[0].x_ll = left_uv_to_xyz_vector (U_Active[0], V_Active[0]);
   g[0].x_lr = left_uv_to_xyz_vector (U_Active[1], V_Active[0]);

   /* MIDDLE */
   g[1].id = 2;
   g[1].x_ul = middle_uv_to_xyz_vector (U_Active[0], V_Active[3]);
   g[1].x_ur = middle_uv_to_xyz_vector (U_Active[1], V_Active[3]);
   g[1].x_ll = middle_uv_to_xyz_vector (U_Active[0], V_Active[2]);
   g[1].x_lr = middle_uv_to_xyz_vector (U_Active[1], V_Active[2]);

   /* RIGHT */
   g[2].id = 3;
   g[2].x_ul = right_uv_to_xyz_vector (U_Active[0], V_Active[5]);
   g[2].x_ur = right_uv_to_xyz_vector (U_Active[1], V_Active[5]);
   g[2].x_ll = right_uv_to_xyz_vector (U_Active[0], V_Active[4]);
   g[2].x_lr = right_uv_to_xyz_vector (U_Active[1], V_Active[4]);
#if 0
   /* These are for AXAF-HRC-2.4S */
   g[0].tdet_xoff = 16384.0;
   g[0].tdet_yoff = 1024.0;
   g[1].tdet_xoff = 16128.0;
   g[1].tdet_yoff = 6144.0;
   g[2].tdet_xoff = 16384.0;
   g[2].tdet_yoff = 11264.0;
#endif
#if 0
   /* These are for AXAF-HRC-2.6S */
   g[0].tdet_xoff = 49368.0;
   g[0].tdet_yoff = 0.0;
   g[1].tdet_xoff = 32912.0;
   g[1].tdet_yoff = 0.0;
   g[2].tdet_xoff = 16456.0;
   g[2].tdet_yoff = 0.0;
#endif
#endif				       /* USE_CALDB_FILES */
#if 1
   /* These are for AXAF-HRC-2.7S */
   g[0].tdet_xoff = 0.0;
   g[0].tdet_yoff = -12.0;
   g[1].tdet_xoff = 0.0;
   g[1].tdet_yoff = 16444.0;
   g[2].tdet_xoff = 0.0;
   g[2].tdet_yoff = 32900.0;
#endif
   for (i = 0; i < 3; i++)
     {
	double u, v;
	(void) _marx_hrc_s_compute_pixel (g[i].id, 0, 0,
					  &g[i].xpixel_offset,
					  &g[i].ypixel_offset,
					  &u, &v);
	rotate_mcp (g+i, Detector_Theta);
     }
   return 0;
}

/* dx, dy are measured from the LL corner, which by construction is the
 * active pixel for that corner.
 */
int _marx_hrc_s_compute_pixel (int id, double dx, double dy,
			       double *xpixel, double *ypixel,
			       double *upixel, double *vpixel)
{
   double u, v;
   double u_0, v_0, x_0, y_0;

   switch (id)
     {
      case 3:
	u = U_Active[0];
	v = V_Active[4];
	u_0 = Right_UV[0];
	v_0 = Right_UV[1];
	x_0 = Right_CXCY[0];
	y_0 = Right_CXCY[1];
	break;

      case 2:
	u = U_Active[0];
	v = V_Active[2];
	u_0 = Middle_UV[0];
	v_0 = Middle_UV[1];
	x_0 = Middle_CXCY[0];
	y_0 = Middle_CXCY[1];
	break;

      case 1:
	u = U_Active[0];
	v = V_Active[0];
	u_0 = Left_UV[0];
	v_0 = Left_UV[1];
	x_0 = Left_CXCY[0];
	y_0 = Left_CXCY[1];
	break;

      default:
	return -1;
     }

   u += dx / U_Pixel_Size;
   v += dy / V_Pixel_Size;

   *upixel = u;
   *vpixel = v;

   *xpixel = x_0 + (u - u_0);
   *ypixel = y_0 + (v - v_0);

   return 0;
}

int _marx_hrc_s_get_pixel_size (double *dx, double *dy)
{
   *dx = U_Pixel_Size;
   *dy = V_Pixel_Size;
   return 0;
}
static _Marx_Simple_Data_Type Array_Data_Table [] =
{
   {"V_Pixel_Size",	1, &V_Pixel_Size,	1.0, 0},
   {"V_Physical",	6, V_Physical,		1.0, 0},
   {"V_Active",		6, V_Active,		1.0, 0},

   {"U_Pixel_Size",	1, &U_Pixel_Size,	1.0, 0},
   {"U_Physical",	2, U_Physical,		1.0, 0},
   {"U_Active",		2, U_Active,		1.0, 0},

   {"Aimpoint_Offset",	3, Aimpoint_Offset,	1.0, 0},

   {"Left_UV",		2, Left_UV,		1.0, 0},
   {"Right_UV",		2, Right_UV,		1.0, 0},
   {"Middle_UV",	2, Middle_UV,		1.0, 0},

   {"Left_XYZ",		3, Left_XYZ,		1.0, 0},
   {"Right_XYZ",	3, Right_XYZ,		1.0, 0},
   {"Middle_XYZ",	3, Middle_XYZ,		1.0, 0},

   {"Left_CXCY",	2, Left_CXCY,		1.0, 0},
   {"Right_CXCY",	2, Right_CXCY,		1.0, 0},
   {"Middle_CXCY",	2, Middle_CXCY,		1.0, 0},

   {"Left_Angle",	1, &Left_Theta, 	PI/180.0, 0},
   {"Right_Angle",	1, &Right_Theta,	PI/180.0, 0},
   {"Detector_Angle",	1, &Detector_Theta,	PI/180.0, 0},
   {NULL, 0, NULL, 0.0, 0}
};

/* Note: This may be called with pf == NULL */
int _marx_hrc_s_geom_init (Param_File_Type *pf)
{
   char *file;
   unsigned int i;
   static int initialized = 0;

   (void) pf;

   if (initialized)
     return 0;

   file = "hrc/hrc_s_geom.txt";

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;
#if 0
   marx_message ("\t%s\n", file);
#endif
   if (-1 == _marx_read_simple_data_file (file, Array_Data_Table))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   /* The XYZ coordinates that have been read in are with respect to the
    * center position.
    */

   for (i = 0; i < 3; i++)
     {
	Left_XYZ[i] += Aimpoint_Offset[i];
	Right_XYZ[i] += Aimpoint_Offset[i];
	Middle_XYZ[i] += Aimpoint_Offset[i];
     }

   setup_alphas_and_betas ();
   /* (void) _marx_set_detector_angle (Detector_Theta); */

   initialized = 1;
   return 0;
}

static int
hrc_s_to_tiled (Marx_Detector_Type *det,
		Marx_Detector_Geometry_Type *g,
		int chip, unsigned int x, unsigned int y,
		unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) chip;
   (void) det;

#if 1
   /* For AXAF-HRC-2.7S */
   /* Note: The coordinate document that I was given appears to be messed up
    * for this coord system.  At least, I think it is, so until I am proven
    * wrong, use this:
    */
# if 1
   yf = y;
   yf = yf + g->tdet_yoff;
   xf = x + g->tdet_xoff;
# else
   /* instead of this: */
   xf = x + g->tdet_xoff;
   yf = y + g->tdet_yoff;
# endif
#else
   /* For AXAF-HRC-2.6S */
   yf = y;
   xf = -yf + g->tdet_xoff;
   yf = x + g->tdet_yoff;
#endif

   if (xf < 0.0) xf = 0.0;
   if (yf < 0.0) yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;
   return 0;
}

static Marx_Detector_Type HRC_S_Detector;
static Marx_Detector_Geometry_Type HRC_S_Geom[_MARX_NUM_HRC_S_CHIPS];

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

Marx_Detector_Type *_marx_get_hrc_s_detector (void)
{
   int id;
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *g;

   d = &HRC_S_Detector;
   if (d->is_initialized)
     return d;

   d->detector_type = MARX_DETECTOR_HRC_S;
   d->tiled_pixel_map_fun = &hrc_s_to_tiled;
   d->facet_list = _marx_link_detector_facet_list (HRC_S_Geom, _MARX_NUM_HRC_S_CHIPS, sizeof(Marx_Detector_Geometry_Type));
   d->fp_system_name = "AXAF-FP-2.3";
   d->first_facet_id = 1;
   d->last_facet_id = 3;
   d->print_info = print_info;

   g = d->facet_list;
   for (id = 0; id < _MARX_NUM_HRC_S_CHIPS; id++)
     {
	g->id = id+1;
	g->tdet_xoff = 0;
	g->tdet_yoff = 0;
	g->num_x_pixels = 4096;
	g->num_y_pixels = 16384;
	g->x_pixel_size = 6.429e-3;
	g->y_pixel_size = 6.429e-3;
	g++;
     }

   if (-1 == patch_hrc_s_geom (d))
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
