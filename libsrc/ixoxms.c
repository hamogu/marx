/*
    This file is part of MARX

    Copyright (C) 2011-2013 Massachusetts Institute of Technology

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

#define MARX_DET_FACET_PRIVATE_DATA \
   double energy_gain; \
   double sigma; \
   struct _Marx_QE_Type *qeinfo;

#include "marx.h"
#include "_marx.h"

#define NUM_XMS_FACETS	5

static double Inner_Pixel_Size = 300*1e-3;   /* mm */
static double Outer_Pixel_Size = 600*1e-3;   /* mm */
static double Inner_Outer_Gap = 0.0;

#define FWHM_TO_SIGMA(X) ((X)/2.3548200450309493)

static unsigned int Inner_Num_Pixels = 40;
static unsigned int Outer_Num_X_Pixels = 36;
static unsigned int Outer_Num_Y_Pixels = 16;

static double Inner_Sigma = FWHM_TO_SIGMA(2.5)*0.001; /* keV */
static double Outer_Sigma = FWHM_TO_SIGMA(10)*0.001; /* keV */
static double Inner_Gain = 1.25*0.001; /* keV/channel */
static double Outer_Gain = 1.25*0.001; /* keV/channel */

static char *Inner_QE_File;
static char *Outer_QE_File;
static Param_Table_Type IXO_CCD_Parm_Table [] =
{
   {"IXO_XMS_Inner_QE_File",	PF_FILE_TYPE,	&Inner_QE_File},
   {"IXO_XMS_Outer_QE_File",	PF_FILE_TYPE,	&Outer_QE_File},
   {"IXO_XMS_Inner_Sigma",	PF_REAL_TYPE,	&Inner_Sigma},
   {"IXO_XMS_Outer_Sigma",	PF_REAL_TYPE,	&Outer_Sigma},
   {NULL, 0, NULL}
};

static int compute_pha (Marx_Detector_Geometry_Type *g, double x, double y,
			double energy, float *pi)
{
   double gain = g->energy_gain;
   int pha;

   (void) x;
   (void) y;

   energy = energy + g->sigma * JDMgaussian_random ();
   if (energy < 0.0) energy = 0.0;

   pha = (1 + energy/gain);
   if (pha <= 0)
     {
	*pi = 0.0;
	return 0;
     }
   energy = (pha - JDMrandom())*gain;
   *pi = (float) energy;
   return pha;
}

static int compute_qe (Marx_Detector_Geometry_Type *g,
			double xpixel, double ypixel, double energy,
			double *qep)
{
   (void) xpixel; (void) ypixel;

   *qep = _marx_qe_compute (g->qeinfo, energy);
   return 0;
}

static int apply_qe_and_pha (Marx_Detector_Geometry_Type *g, Marx_Photon_Attr_Type *at)
{
   if (_Marx_Det_Ideal_Flag == 0)
     {
	double qe;

	if (-1 == (*g->qe_fun)(g, at->y_pixel, at->z_pixel, at->energy, &qe))
	  return -1;

	if (JDMrandom () > qe)
	  {
	     at->flags |= PHOTON_UNDETECTED;
	     return -1;
	  }
     }

   if (-1 == (at->pulse_height = (*g->pha_fun) (g, at->y_pixel, at->z_pixel, at->energy, &at->pi)))
     {
	at->flags |= PHOTON_UNDETECTED;
	return -1;
     }

   return 0;
}

static int read_ixo_ccd_parms (Param_File_Type *p)
{
   if (-1 == pf_get_parameters (p, IXO_CCD_Parm_Table))
     return -1;

   Inner_Sigma = FWHM_TO_SIGMA(Inner_Sigma);
   Outer_Sigma = FWHM_TO_SIGMA(Outer_Sigma);

   return 0;
}

/* The XMS looks like:
 *
 *    ╔══════════════════╦═══════╗ z3
 *    ║               <═*║       ║
 *    ║    3             ║       ║
 *    ╠═══════╦══════════╣       ║ z2
 *    ║*      ║          ║       ║
 *    ║║      ║    0     ║      ^║
 *    ║v      ║          ║   2  ║║
 *    ║   4   ║*═>       ║      *║
 *    ║       ╠══════════╩═══════╣ z1
 *    ║       ║              1   ║
 *    ║       ║*═>               ║
 *    ╚═══════╩══════════════════╝ z0
 *    y0     y1         y2       y3
 *
 * I have arbitrarily labeled the LL locations with a "*" and the -> label the
 * local +X pixel axis.
 */

static Marx_Detector_Geometry_Type *setup_geom (void)
{
   Marx_Detector_Geometry_Type *geom, *d;
   JDMVector_Type x_ll, x_lr, x_ul, x_ur;
   double dx, dy;
   double y_0, y_1, y_2, y_3, z_0, z_1, z_2, z_3;
   double gap = 0.0;

   geom = (Marx_Detector_Geometry_Type *)marx_calloc (NUM_XMS_FACETS, sizeof(Marx_Detector_Geometry_Type));
   if (geom == NULL)
     return NULL;

   /* XMS facet 0 */
   dx = Inner_Num_Pixels * Inner_Pixel_Size;
   dy = Inner_Num_Pixels * Inner_Pixel_Size;

   y_1 = -0.5*dx; y_2 = 0.5*dx;
   z_1 = -0.5*dy; z_2 = 0.5*dy;

   x_ll.x = 0; x_ll.y = y_1; x_ll.z = z_1;
   x_ul.x = 0; x_ul.y = y_1; x_ul.z = z_2;
   x_ur.x = 0; x_ur.y = y_2; x_ur.z = z_2;
   x_lr.x = 0; x_lr.y = y_2; x_lr.z = z_1;
   d = geom+0; d->x_ll = x_ll; d->x_lr = x_lr; d->x_ur = x_ur; d->x_ul = x_ul;
   d->id = 0;
   d->x_pixel_size = Inner_Pixel_Size;
   d->y_pixel_size = Inner_Pixel_Size;
   d->num_x_pixels = Inner_Num_Pixels;
   d->num_y_pixels = Inner_Num_Pixels;

   dx = Outer_Num_X_Pixels * Outer_Pixel_Size;
   dy = Outer_Num_Y_Pixels * Outer_Pixel_Size;
   gap = Inner_Outer_Gap;

   y_0 = y_1 - (dy + gap); y_3 = y_1 + dx;
   z_0 = z_1 - (dy + gap); z_3 = z_1 + dx;
   /* Note: It is assumed that gap is such that y_3=y_2+gap+dy*/

   /* 1 */
   x_ll.x = 0; x_ll.y = y_1; x_ll.z = z_0;
   x_ul.x = 0; x_ul.y = y_1; x_ul.z = z_0+dy;
   x_ur.x = 0; x_ur.y = y_1+dx; x_ur.z = z_0+dy;
   x_lr.x = 0; x_lr.y = y_1+dx; x_lr.z = z_0;
   d = geom+1; d->x_ll = x_ll; d->x_lr = x_lr; d->x_ur = x_ur; d->x_ul = x_ul;
   d->id = 1;
   d->x_pixel_size = Outer_Pixel_Size;
   d->y_pixel_size = Outer_Pixel_Size;
   d->num_x_pixels = Outer_Num_X_Pixels;
   d->num_y_pixels = Outer_Num_Y_Pixels;

   /* 2 */
   x_ll.x = 0; x_ll.y = y_3; x_ll.z = z_1;
   x_ul.x = 0; x_ul.y = y_3-dy; x_ul.z = z_1;
   x_ur.x = 0; x_ur.y = y_3-dy; x_ur.z = z_1+dx;
   x_lr.x = 0; x_lr.y = y_3; x_lr.z = z_1+dx;
   d = geom+2; d->x_ll = x_ll; d->x_lr = x_lr; d->x_ur = x_ur; d->x_ul = x_ul;
   d->id = 2;
   d->x_pixel_size = Outer_Pixel_Size;
   d->y_pixel_size = Outer_Pixel_Size;
   d->num_x_pixels = Outer_Num_X_Pixels;
   d->num_y_pixels = Outer_Num_Y_Pixels;

   /* 3 */
   x_ll.x = 0; x_ll.y = y_2; x_ll.z = z_3;
   x_ul.x = 0; x_ul.y = y_2; x_ul.z = z_3-dy;
   x_ur.x = 0; x_ur.y = y_2-dx; x_ur.z = z_3-dy;
   x_lr.x = 0; x_lr.y = y_2-dx; x_lr.z = z_3;
   d = geom+3; d->x_ll = x_ll; d->x_lr = x_lr; d->x_ur = x_ur; d->x_ul = x_ul;
   d->id = 3;
   d->x_pixel_size = Outer_Pixel_Size;
   d->y_pixel_size = Outer_Pixel_Size;
   d->num_x_pixels = Outer_Num_X_Pixels;
   d->num_y_pixels = Outer_Num_Y_Pixels;

   /* 4 */
   x_ll.x = 0; x_ll.y = y_0; x_ll.z = z_2;
   x_ul.x = 0; x_ul.y = y_0+dy; x_ul.z = z_2;
   x_ur.x = 0; x_ur.y = y_0+dy; x_ur.z = z_2-dx;
   x_lr.x = 0; x_lr.y = y_0; x_lr.z = z_2-dx;
   d = geom+4; d->x_ll = x_ll; d->x_lr = x_lr; d->x_ur = x_ur; d->x_ul = x_ul;
   d->id = 4;
   d->x_pixel_size = Outer_Pixel_Size;
   d->y_pixel_size = Outer_Pixel_Size;
   d->num_x_pixels = Outer_Num_X_Pixels;
   d->num_y_pixels = Outer_Num_Y_Pixels;

   return geom;
}

static Marx_Detector_Type IXO_XMS_Detector;

Marx_Detector_Type *_marx_get_ixo_xms_detector (void)
{
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *g;

   d = &IXO_XMS_Detector;
   if (d->is_initialized)
     return d;

   g = setup_geom ();
   if (g == NULL)
     return NULL;

   d->detector_type = MARX_DETECTOR_IXO_XMS;
   d->tiled_pixel_map_fun = NULL;
   d->facet_list = _marx_link_detector_facet_list (g, NUM_XMS_FACETS, sizeof(Marx_Detector_Geometry_Type));
   d->fp_system_name = NULL;
   d->first_facet_id = 0;
   d->last_facet_id = 4;

   if (-1 == _marx_compute_detector_basis (d))
     {
	marx_free ((char *)d);
	return NULL;
     }
   return d;
}

int _marx_ixoxms_init (Param_File_Type *pf)
{
   Marx_QE_Type *qeinfo;
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *xms;

   if (-1 == read_ixo_ccd_parms (pf))
     return -1;

   if (NULL == (d = _marx_get_ixo_xms_detector ()))
     return -1;

   qeinfo = NULL;
   if (_Marx_Det_Ideal_Flag == 0)
     {
	char *file;

	if (_Marx_Det_Ideal_Flag == 0)
	  marx_message ("Reading IXO XMS QE/Filter Files\n");

	if (NULL == (file = marx_make_data_file_name (Inner_QE_File)))
	  return -1;

	marx_message ("\t%s\n", file);
	qeinfo = _marx_qe_read_file (file, "IXO_XMS_QE", "ENERGY", "QE", NULL);
	marx_free (file);

	if (qeinfo == NULL)
	  return -1;

	xms = d->facet_list;
	xms->qeinfo = qeinfo;

	xms = xms->next;

	if (NULL == (file = marx_make_data_file_name (Outer_QE_File)))
	  return -1;
	marx_message ("\t%s\n", file);
	qeinfo = _marx_qe_read_file (file, "IXO_XMS_QE_OUTER", "ENERGY", "QE", NULL);
	marx_free (file);
	if (qeinfo == NULL)
	  return -1;

	while (xms != NULL)
	  {
	     _marx_qe_inc_ref (qeinfo);
	     xms->qeinfo = qeinfo;
	     xms = xms->next;
	  }
	_marx_qe_free (qeinfo);
     }

   /* inner */
   xms = d->facet_list;
   xms->energy_gain = Inner_Gain;
   xms->sigma = Inner_Sigma;
   xms->pha_fun = compute_pha;
   xms->qe_fun = compute_qe;

   xms = xms->next;
   while (xms != NULL)
     {
	xms->energy_gain = Outer_Gain;
	xms->sigma = Outer_Sigma;
	xms->pha_fun = compute_pha;
	xms->qe_fun = compute_qe;
	xms = xms->next;
     }

   return 0;
}

int _marx_ixoxms_detect (Marx_Photon_Type *pt)
{
   Marx_Photon_Attr_Type *at, *attrs;
   unsigned int n_photons, i;
   unsigned int *sorted_index;
   Marx_Detector_Geometry_Type *facet_list;

   if (pt->history & MARX_DET_NUM_OK)
     return 0;

   pt->history |= (MARX_DET_PIXEL_OK | MARX_DET_NUM_OK
		   | MARX_PULSEHEIGHT_OK | MARX_PI_OK);

   marx_prune_photons (pt);

   attrs = pt->attributes;
   n_photons = pt->num_sorted;
   sorted_index = pt->sorted_index;

   facet_list = IXO_XMS_Detector.facet_list;

   for (i = 0; i < n_photons; i++)
     {
	Marx_Detector_Geometry_Type *d;
	double dx, dy;

	at = attrs + sorted_index[i];

#if MARX_HAS_DITHER
	_marx_dither_detector (&at->dither_state);
#endif
	/* Transform ray into local system */
	_marx_transform_ray (&at->x, &at->p,
			     &_Marx_Det_XForm_Matrix);

	/* See if the photon will hit the detector and if so, which one. */
	d = _marx_intersect_with_detector (at->x, at->p,
					   facet_list,
					   &at->x, &dx, &dy,
					   _Marx_Det_Extend_Flag);
	if (d == NULL)
	  {
	     at->flags |= PHOTON_MISSED_DETECTOR;
	     at->ccd_num = -1;
	  }
	else
	  {
	     at->ccd_num = d->id;
	     at->y_pixel = dx / d->x_pixel_size;
	     at->z_pixel = dy / d->y_pixel_size;

	     (void) apply_qe_and_pha (d, at);
	  }
	/* Transform detected photon back to original system */
	_marx_transform_ray_reverse (&at->x, &at->p,
				     &_Marx_Det_XForm_Matrix);
#if MARX_HAS_DITHER
	_marx_undither_detector (&at->dither_state);
#endif
     }

   return 0;
}
