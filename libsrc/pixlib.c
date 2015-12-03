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

#include "marx.h"

/* The chip to sky transformation goes from the chip pixel coordinates, on
 * a specific detector plane, to mirror nodal coordinates, and then expressed
 * as an elevation and an azimuth.  The transformations implied by this
 * are:
 *
 *     CHIP -> LSI -> STT -> STF -> FC -> MNC
 *
 * Marx itself handles some of this by expressing the detector geometry
 * in the STF system for the nominal focus.  To go from STF->FC requires
 * knowledge of part of the aspect solution (DY, dZ, and Theta).  Also, if
 * the detector is not at the nominal focal position, that offset will have
 * to be added in here.
 *
 * The FC->MNC is simply a translation by the focal length along the optical
 * axis.
 *
 * The structure defined below handles the transformation from CHIP to MNC
 * for a given aspect and detector offset.
 */

#include <jdmath.h>

struct _Marx_Chip_To_MNC_Type
{
   JDMVector_Type mnc_to_chip_offset;
   /* offset from mnc to chip system.  Include focal length, detector offsets,
    * and aspect corrections.
    */

   /* These are the chip unit vectors expressed in the MNC system.  Actually,
    * they are scaled by the pixel sizes.
    */
   JDMVector_Type chip_e1;
   JDMVector_Type chip_e2;

   /* These have the advantage that chip_e1_inv.chip_e1 = 1. */
   JDMVector_Type chip_e1_inv;     /* scaled by 1/pixel size */
   JDMVector_Type chip_e2_inv;     /* scaled by 1/pixel size */

   /* Orthogonal to above two-- unit vector */
   JDMVector_Type chip_e3_hat;

   double offset_dot_e3;

   double min_x_pixel;
   double max_x_pixel;
   double min_y_pixel;
   double max_y_pixel;

   double xpixel_offset;
   double ypixel_offset;

   Marx_Detector_Type *det;
   Marx_Detector_Geometry_Type *geom;
};

void
marx_free_chip_to_mnc (Marx_Chip_To_MNC_Type *m)
{
   marx_free ((char *)m);
}

Marx_Chip_To_MNC_Type *
marx_allocate_chip_to_mnc (char *name, int id)
{
   Marx_Detector_Type *det;
   Marx_Chip_To_MNC_Type *chip2mnc;
   Marx_Detector_Geometry_Type *g;

   if (NULL == (det = marx_get_detector_info (name)))
     return NULL;

   g = det->facet_list;
   while (g != NULL)
     {
	if (g->id == id)
	  break;
	g = g->next;
     }
   if (g == NULL)
     {
	marx_error ("marx_allocate_chip_to_mnc: chip_id out of range");
	return NULL;
     }

   chip2mnc = (Marx_Chip_To_MNC_Type *) marx_malloc (sizeof (Marx_Chip_To_MNC_Type));
   if (chip2mnc == NULL)
     return NULL;

   memset ((char *) chip2mnc, 0, sizeof (Marx_Chip_To_MNC_Type));

   chip2mnc->det = det;
   chip2mnc->geom = g;

   chip2mnc->min_x_pixel = 0;
   chip2mnc->max_x_pixel = g->num_x_pixels;
   chip2mnc->min_y_pixel = 0;
   chip2mnc->max_y_pixel = g->num_y_pixels;

   return chip2mnc;
}

int
marx_init_chip_to_mnc (Marx_Chip_To_MNC_Type *chip2mnc,
		       double focal_length,
		       double xoff, double yoff, double zoff, double theta)
{
   double cos_theta, sin_theta;
   JDMVector_Type *e1, *e2, *xhat, *yhat;
   Marx_Detector_Type *det;
   Marx_Detector_Geometry_Type *g;
   double x_pixel_size, y_pixel_size;

   det = chip2mnc->det;
   g = chip2mnc->geom;

   /* The g->xhat, g->yhat, and g->x_ll vectors are already expressed in the
    * STF coordinate system.  In particular, g->x_ll specifies the origin of
    * the chip system in STF coordinates.  We need to rotate this by theta
    * about the optical axis (X) into the MNC frame.  The resulting vector
    * is translated via the offset parameters to the MNC origin.
    */

   cos_theta = cos (theta);
   sin_theta = sin (theta);

   xhat = &g->x_ll;

   chip2mnc->mnc_to_chip_offset.x = xoff - focal_length;
   chip2mnc->mnc_to_chip_offset.y = yoff + cos_theta * xhat->y - sin_theta * xhat->z;
   chip2mnc->mnc_to_chip_offset.z = zoff + sin_theta * xhat->y + cos_theta * xhat->z;

   /* These are the pixel values at the LL corner */
   chip2mnc->xpixel_offset = g->xpixel_offset;
   chip2mnc->ypixel_offset = g->ypixel_offset;

   /* Now do the unit vectors. */
   e1 = &chip2mnc->chip_e1;
   e2 = &chip2mnc->chip_e2;

   xhat = &g->xhat;
   yhat = &g->yhat;

   e1->x = xhat->x;
   e1->y = cos_theta * xhat->y - sin_theta * xhat->z;
   e1->z = sin_theta * xhat->y + cos_theta * xhat->z;

   e2->x = yhat->x;
   e2->y = cos_theta * yhat->y - sin_theta * yhat->z;
   e2->z = sin_theta * yhat->y + cos_theta * yhat->z;

   chip2mnc->chip_e3_hat = JDMv_pcross_prod (e1, e2);

   /* Now scale by pixel sizes */
   x_pixel_size = g->x_pixel_size;
   y_pixel_size = g->y_pixel_size;

   chip2mnc->chip_e1_inv = JDMv_smul (1.0 / x_pixel_size, *e1);
   chip2mnc->chip_e2_inv = JDMv_smul (1.0 / y_pixel_size, *e2);
   *e1 = JDMv_smul (x_pixel_size, *e1);
   *e2 = JDMv_smul (y_pixel_size, *e2);

   chip2mnc->offset_dot_e3 = JDMv_pdot_prod (&chip2mnc->mnc_to_chip_offset,
					     &chip2mnc->chip_e3_hat);

   return 0;
}

int
marx_chip_to_mnc (Marx_Chip_To_MNC_Type *chip2mnc,
		  double xpixel, double ypixel,
		  JDMVector_Type *mnc)
{
   JDMVector_Type *ofs, *e1, *e2;

   /* Here, we assume that pixels are measure from the LL corner */
   xpixel -= chip2mnc->xpixel_offset;
   ypixel -= chip2mnc->ypixel_offset;

   ofs = &chip2mnc->mnc_to_chip_offset;
   e1 = &chip2mnc->chip_e1;
   e2 = &chip2mnc->chip_e2;

   mnc->x = ofs->x + e1->x * xpixel + e2->x * ypixel;
   mnc->y = ofs->y + e1->y * xpixel + e2->y * ypixel;
   mnc->z = ofs->z + e1->z * xpixel + e2->z * ypixel;

   JDMv_normalize (mnc);

   return 0;
}

static int generic_to_chip (Marx_Chip_To_MNC_Type *chip2mnc,
			    JDMVector_Type *mnc,
			    JDMVector_Type *ofs,
			    double ofs_dot_e3,
			    double *xpixel, double *ypixel)
{
   JDMVector_Type diff;
   double x, y;
   double scale;

   /* Note that the sign of the mnc vector does not matter. Thus this
    * routine is insensitive to whether mnc points to the sky or from
    * the sky.
    */
   scale = ofs_dot_e3 / JDMv_pdot_prod (mnc, &chip2mnc->chip_e3_hat);
   diff.x = scale * mnc->x - ofs->x;
   diff.y = scale * mnc->y - ofs->y;
   diff.z = scale * mnc->z - ofs->z;

   x = JDMv_pdot_prod (&diff, &chip2mnc->chip_e1_inv);
   y = JDMv_pdot_prod (&diff, &chip2mnc->chip_e2_inv);

   *xpixel = x + chip2mnc->xpixel_offset;
   *ypixel = y + chip2mnc->ypixel_offset;

   if (((x < chip2mnc->min_x_pixel) || (x >= chip2mnc->max_x_pixel))
       || ((y < chip2mnc->min_x_pixel) || (y >= chip2mnc->max_x_pixel)))
     return -1;

   return 0;
}

int
marx_mnc_to_chip (Marx_Chip_To_MNC_Type *chip2mnc,
		  JDMVector_Type *mnc,
		  double *xpixel, double *ypixel)
{
   return generic_to_chip (chip2mnc, mnc, &chip2mnc->mnc_to_chip_offset,
			   chip2mnc->offset_dot_e3,
			   xpixel, ypixel);
}

typedef struct
{
   double period;		       /* microns */
   double dispersion_angle;	       /* radians */
}
Grating_Info_Type;

static Grating_Info_Type HEG_Info =
{
   0.200081,
   -5.18 * (PI / 180.0)
};

static Grating_Info_Type MEG_Info =
{
   0.400141,
   4.75 * (PI / 180.0)
};

static Grating_Info_Type LEG_Info =
{
   0.991249,
   0.0 * (PI / 180.0)
};

static Grating_Info_Type *get_grating_info (int type)
{
   switch (type)
     {
      case MARX_GRATING_HEG:
	return &HEG_Info;

      case MARX_GRATING_MEG:
	return &MEG_Info;

      case MARX_GRATING_LEG:
	return &LEG_Info;
     }

   marx_error ("get_grating_info: Grating type %d is unknown", type);
   return NULL;
}

struct _Marx_Grating_Xform_Type
{
   JDMVector_Type grating_node_offset; /* from MNC origin */
   JDMVector_Type l_hat, d_hat, n_hat; /* unit vectors of grating facet */

   JDMVector_Type mnc_to_chip_offset;
   /* Location of chip origin wrt the MNC */

   /* The rest of the quantities depend upon the ray to be diffracted. */
   JDMVector_Type grating_to_chip_offset;
   /* Offset from diffracted ray at grating to the to chip origin.  Includes
    * focal length, detector offsets, and aspect correction.
    */
   JDMVector_Type p0;		       /* ray incident upon grating */
   double p_dot_d;
   double p_dot_l;
   double period;
};

void marx_free_grating_xform (Marx_Grating_Xform_Type *g)
{
   if (g == NULL) return;
   marx_free ((char *)g);
}

Marx_Grating_Xform_Type *
marx_allocate_grating_xform (Marx_Chip_To_MNC_Type *c2mnc, int type)
{
   Marx_Grating_Xform_Type *g;
   Grating_Info_Type *ginfo;

   if (NULL == (ginfo = get_grating_info (type)))
     return NULL;

   if (c2mnc == NULL)
     {
	marx_error ("marx_allocate_grating_xform: NULL passed");
	return NULL;
     }

   g = (Marx_Grating_Xform_Type *) marx_malloc (sizeof (Marx_Grating_Xform_Type));
   if (g == NULL)
     return NULL;

   memset ((char *) g, 0, sizeof (Marx_Grating_Xform_Type));

   g->l_hat = JDMv_vector (0.0, 0.0, 1.0);
   g->d_hat = JDMv_vector (0.0, 1.0, 0.0);
   g->n_hat = JDMv_vector (-1.0, 0.0, 0.0);

   g->l_hat = JDMv_rotate_unit_vector (g->l_hat, g->n_hat, -1 * ginfo->dispersion_angle);
   g->d_hat = JDMv_rotate_unit_vector (g->d_hat, g->n_hat, -1 * ginfo->dispersion_angle);

   g->grating_node_offset = JDMv_vector (-1431.81,0,0);
   /* From equation 88 of rev 4.2 coordinate memo */

   g->mnc_to_chip_offset = c2mnc->mnc_to_chip_offset;

   g->period = ginfo->period;
   return g;
}

/* p specifies a ray in the mirror nodal system.  That is, it represents
 * a ray at the mirror node traveling to the grating.  Here we compute
 * the intersection on the grating and set up values that will be used
 * later for the diffraction.
 */
int
marx_init_grating_xform (Marx_Grating_Xform_Type *g,
			 JDMVector_Type *p)
{
   double t;

   if (g == NULL)
     return -1;

   t = JDMv_dot_prod (g->n_hat, *p);
   if (t == 0.0)
     return -1;

   t = JDMv_dot_prod (g->grating_node_offset, g->n_hat) / t;

   g->grating_to_chip_offset = JDMv_diff (g->mnc_to_chip_offset,
					  JDMv_smul (t, *p));
   g->p0 = *p;

   g->p_dot_d = JDMv_dot_prod (*p, g->d_hat);
   g->p_dot_l = JDMv_dot_prod (*p, g->l_hat);

   return 0;
}

static int diffract_ray (Marx_Grating_Xform_Type *g,
			 int n, double lambda,
			 JDMVector_Type *pp)
{
   double p_d, p_l, p_n;

   p_d = n * lambda / g->period + g->p_dot_d;
   p_l = g->p_dot_l;

   p_n = 1.0 - p_l * p_l - p_d * p_d;
   if (p_n < 0.0)
     return -1;
   p_n = sqrt (p_n);

   *pp = JDMv_ax1_bx2 (p_d, g->d_hat,
		       1.0, JDMv_ax1_bx2 (p_l, g->l_hat,
					  p_n, g->n_hat));
   JDMv_normalize (pp);
   return 0;
}

int marx_grating_to_chip (Marx_Grating_Xform_Type *g,
			  Marx_Chip_To_MNC_Type *c,
			  double lambda, int n,
			  double *xpixel, double *ypixel)
{
   JDMVector_Type p;

   if (-1 == diffract_ray (g, n, lambda, &p))
     return -1;

   return generic_to_chip (c, &p, &g->grating_to_chip_offset,
			   JDMv_dot_prod (g->grating_to_chip_offset,
					  c->chip_e3_hat),
			   xpixel, ypixel);
}

/* In the following pair of functions, p0 is a unit vector that specifies
 * the location of the tangent plane.  The (x,y,z) components of the vector
 * are such that (0,0,1) is at the celestial pole, (1,0,0) is at ra=dec=0.
 */

/* This routine computes the tangent plane coordinate of a vector p
 * where the unit vector p0 specifies the tangent plane.
 */
int
marx_vector_to_tan_plane (JDMVector_Type *p,
			  JDMVector_Type *p0,
			  double *tx, double *ty)
{
   double p0_perp, p0p_perp;
   double p0_x, p0_y, p0_z;
   double p_x, p_y, p_z;
   double den;

   p0_x = p0->x; p0_y = p0->y; p0_z = p0->z;
   p_x = p->x; p_y = p->y; p_z = p->z;

   p0_perp = p0_x * p0_x + p0_y * p0_y;
   p0p_perp = p0_x * p_x + p0_y * p_y;

   /* If den == 0, then point is at inifinity. */
   den = sqrt(p0_perp) * (p0p_perp + p_z * p0_z);

   *tx = (p0_x * p_y - p0_y * p_x) / den;
   *ty = (p_z * p0_perp - p0_z * p0p_perp)/den;

   return 0;
}

int
marx_tan_plane_to_vector (double tx, double ty, JDMVector_Type *p0, JDMVector_Type *p)
{
   double p0_perp;
   double p0_x, p0_y, p0_z;
   double a, b;

   p0_x = p0->x;
   p0_y = p0->y;
   p0_z = p0->z;

   p0_perp = sqrt (p0_x*p0_x + p0_y*p0_y);

   /* FIXME!!! p0_perp could be 0 */
   a = p0_x/p0_perp;
   b = p0_y/p0_perp;

   p->z = p0_z + ty * p0_perp;
   ty = ty * p0_z;

   p->x = p0_x - b * tx - a * ty;
   p->y = p0_y + a * tx - b * ty;

   JDMv_normalize (p);
   return 0;
}

void
marx_mnc_to_ra_dec (JDMVector_Type *mnc, double *ra, double *dec)
{
   double perp;
   double x, y;

   /* We need to flip the sign of the mirror nodal coords to point to
    * the sky.  This is take care of in the signs below, and in the
    * comparisons.
    */
   x = mnc->x;
   y = mnc->y;

   perp = sqrt (x * x + y * y);
   if (perp > 1.0) perp = 1.0;	       /* roundoff */

   if (mnc->z <= 0) *dec = acos (perp);
   else *dec = -acos (perp);

   /* x is negative.  flip sign */
   perp = -x / perp;

   if (y <= 0.0) *ra = acos (perp);
   else *ra = -acos (perp);
}

/* Suppose that we are given two points on the sphere with coords
 * (ra_0, dec_0) and (ra, dec).  Denote delta_ra as the RA of the second
 * point as seen from the first point at (ra_0, dec_0).  Similarly, delta_dec
 * is the Dec of the second point with respect to the first point.
 * This function computes the so-called aspect offsets, (delta_ra,delta_dec).
 */
void
marx_compute_ra_dec_offsets (double ra_0, double dec_0,
			     double ra, double dec,
			     double *delta_ra, double *delta_dec)
{
   double d_ra, d_dec;
   double c;
   double factor;
   double sin_delta_dec;
   double num, den;

   d_ra = ra - ra_0;
   d_dec = dec - dec_0;

   c = cos (dec);
   factor = c * (1 - cos (d_ra));

   sin_delta_dec = sin (d_dec) + factor * sin (dec_0);
   if (fabs (sin_delta_dec) > 1.0)
     {
	/* Roundoff error */
	if (sin_delta_dec < 0) sin_delta_dec = -1;
	else sin_delta_dec = 1;
     }

   /* asin returns result in (-PI,PI), which is what I want */
   *delta_dec = asin (sin_delta_dec);

   num = c * sin (d_ra);
   den = cos (d_dec) - factor * cos (dec_0);

   /* We want to compute atan (num/den) with the result in (-PI,PI).  I could
    * use atan2, but some systems do not have it.  atan always returns the
    * result in (-PI/2,PI/2).
    */
   if (den >= 0) *delta_ra = atan (num/den);
   else
     {
	if (num >= 0)
	  *delta_ra = atan (num/den) + PI;
	else
	  *delta_ra = atan (num/den) - PI;
     }
}
