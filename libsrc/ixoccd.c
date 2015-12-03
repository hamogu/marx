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

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>

#include <jdmath.h>
#include <pfile.h>
#include <jdfits.h>

#include "ixoccd.h"
#include "marx.h"
#include "_marx.h"

#define MARX_HAS_IXO_CCD_STREAK 0
#define MAX_NUM_IXO_CCDS	32
static unsigned int Num_IXO_CCDs = 16;
static int Num_X_Pixels = 1024;
static int Num_Y_Pixels = 1024;
static double X_Pixel_Size = 0.024;    /* 24 um */
static double Y_Pixel_Size = 0.024;    /* 24 um */
static double CCD_Gap = 0.5;	       /* 500 um */
static double Center_CCD_Angle = 3.0 * PI/180.0;

static double CCD_Gain = 3.65*0.001;   /* keV/electron */
static double Fano_Factor = 0.12;
static double Read_Noise = 4;			       /* electrons */

struct _IXO_CCD_QE_Type
{
   unsigned int num_refs;
   float *energies;
   unsigned int num_energies;
   float *qes;
   float *filter;
};

static Marx_Detector_Geometry_Type IXO_CCD_List[MAX_NUM_IXO_CCDS];

static void free_ccd_qe_type (IXO_CCD_QE_Type *qeinfo)
{
   if (qeinfo == NULL)
     return;

   if (qeinfo->num_refs > 1)
     {
	qeinfo->num_refs--;
	return;
     }
   /* NULLs ok */
   marx_free ((char *)qeinfo->energies);
   marx_free ((char *)qeinfo->qes);
   marx_free ((char *)qeinfo->filter);
   marx_free ((char *)qeinfo);
}

static IXO_CCD_QE_Type *alloc_ccd_qe_type (unsigned int num_energies)
{
   IXO_CCD_QE_Type *qeinfo;

   if (NULL == (qeinfo = (IXO_CCD_QE_Type *)marx_malloc (sizeof(IXO_CCD_QE_Type))))
     return NULL;
   memset ((char *)qeinfo, 0, sizeof(IXO_CCD_QE_Type));

   qeinfo->num_refs = 1;
   qeinfo->num_energies = num_energies;
   if ((NULL == (qeinfo->energies = (float *)marx_malloc(num_energies*sizeof(float))))
       || (NULL == (qeinfo->qes = (float *)marx_malloc(num_energies*sizeof(float))))
       || (NULL == (qeinfo->filter = (float *)marx_malloc(num_energies*sizeof(float)))))
     {
	free_ccd_qe_type (qeinfo);
	return NULL;
     }
   return qeinfo;
}

static IXO_CCD_QE_Type *read_ccd_qe_file (char *file)
{
   IXO_CCD_QE_Type *qeinfo = NULL;
   JDFits_Type *ft;
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   unsigned int i, num_rows;

   ft = jdfits_open_binary_table (file, "IXO_CCD_QE");
   if (ft == NULL)
     {
	marx_error ("Unable to open IXC CCD QE file %s", file);
	return NULL;
     }

   r = jdfits_bintable_open_rows (ft, 3, "f:ENERGY", "f:QE", "f:filter");
   if (r == NULL)
     goto return_error;

   num_rows = r->num_rows;
   if (num_rows < 2)
     {
	marx_error ("Expecting more than 2 rows in %s", file);
	goto return_error;
     }

   if (NULL == (qeinfo = alloc_ccd_qe_type (num_rows)))
     goto return_error;

   c = r->col_data;
   for (i = 0; i < num_rows; i++)
     {
	if (1 != jdfits_read_next_row (ft, r))
	  {
	     marx_error ("File %s appears to be corrupt", file);
	     goto return_error;
	  }
	qeinfo->energies[i] = c[0].data.f[0];
	qeinfo->qes[i] = c[1].data.f[0];
	qeinfo->filter[i] = c[2].data.f[0];
     }
   jdfits_bintable_close_rows (r);
   (void) jdfits_close_file (ft);
   return qeinfo;

return_error:

   marx_error ("Error processing %s", file);
   /* NULLs ok */
   free_ccd_qe_type (qeinfo);
   jdfits_bintable_close_rows (r);
   (void) jdfits_close_file (ft);
   return NULL;
}

static int apply_qe_and_pha (Marx_Detector_Geometry_Type *ccd, Marx_Photon_Attr_Type *at)
{
   double en = at->energy;

   if (_Marx_Det_Ideal_Flag == 0)
     {
	IXO_CCD_QE_Type *qeinfo;
	double r, qe, filter;

	qeinfo = ccd->qeinfo;

	r = JDMrandom ();

	qe = JDMinterpolate_f (en, qeinfo->energies, qeinfo->qes, qeinfo->num_energies);
	/* FIXME: Add support for separate filter energy grid */
	filter = JDMinterpolate_f (en, qeinfo->energies, qeinfo->filter, qeinfo->num_energies);
	if (r >= qe*filter)
	  {
	     at->flags |= PHOTON_UNDETECTED;
	     return -1;
	  }
     }

   if (-1 == (at->pulse_height = (*ccd->pha_fun) (ccd, at->y_pixel, at->z_pixel, en, &at->pi)))
     {
	at->flags |= PHOTON_UNDETECTED;
	return -1;
     }

   return 0;
}

int _marx_ixoccd_detect (Marx_Photon_Type *pt)
{
   Marx_Photon_Attr_Type *at, *attrs;
   unsigned int n_photons, i;
   unsigned int *sorted_index;
#if MARX_HAS_IXO_CCD_STREAK
   double tstart;
#endif

   if (pt->history & MARX_DET_NUM_OK)
     return 0;

   pt->history |= (MARX_DET_PIXEL_OK | MARX_DET_NUM_OK
		   | MARX_PULSEHEIGHT_OK | MARX_PI_OK);

   marx_prune_photons (pt);

   attrs = pt->attributes;
   n_photons = pt->num_sorted;
   sorted_index = pt->sorted_index;

#if MARX_HAS_IXO_CCD_STREAK
   tstart = pt->start_time;
#endif

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

	/* See if the photon will hit the CCD and if so, which one. */
	d = _marx_intersect_with_detector (at->x, at->p,
					   IXO_CCD_List,
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

	     if (0 == apply_qe_and_pha (d, at))
	       {
#if MARX_HAS_IXO_CCD_STREAK
		  (void) _marx_ixo_apply_streak (tstart, at, d);
#endif
	       }
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

static int compute_pha (Marx_Detector_Geometry_Type *ccd, double x, double y,
			  double energy, float *pi)
{
   double da;
   double gain = ccd->energy_gain;
   double noise = ccd->read_noise;
   int pha;

   (void) x;
   (void) y;

   if (gain == 0.0)
     {
	*pi = 0;
	return 0;
     }

   /* Eq 2.1 of ACIS-PSU-SOP-01 suggests the following: */
   da = gain * sqrt (noise * noise + ccd->fano_factor * energy / gain);

   energy = energy + da * JDMgaussian_random ();
   if (energy < 0.0) energy = 0.0;

   pha = (short) (1 + energy/gain);
   if (pha <= 0)
     {
	*pi = 0.0;
	return 0;
     }
   energy = (pha - JDMrandom())*gain;
   *pi = (float) energy;
   return pha;
}

static int setup_geom (Marx_Detector_Geometry_Type *geom)
{
   Marx_Detector_Geometry_Type *d;
   JDMVector_Type x_ll, x_lr, x_ul, x_ur, axis, cntr;
   double r, delta;
   double delta_phi, phi0;
   double dx, dy;
   unsigned int i;

   /* The method here is to position a detector plane at the origin, and then
    * rotate it into place.  The rotation axis is located at the center of the
    * torus and is directed along the y axis.
    *
    * Looking from the mirror to the detector plane, a CCD at the origin looks
    * like:
    *
    *     xur+-----dy----+ xlr
    *        |           |
    *        |           |
    *      dx|     O     |
    *        |           |
    *        |           |
    *    xul +-----------+ xll
    *
    */

   dx = Num_X_Pixels * X_Pixel_Size;
   dy = Num_Y_Pixels * Y_Pixel_Size;

   /* Center the template CCD on the origin */
   x_ll.x = 0; x_ll.y =  0.5*dy; x_ll.z = -0.5*dx;
   x_ul.x = 0; x_ul.y = -0.5*dy; x_ul.z = -0.5*dx;
   x_ur.x = 0; x_ur.y = -0.5*dy; x_ur.z =  0.5*dx;
   x_lr.x = 0; x_lr.y =  0.5*dy; x_lr.z =  0.5*dx;

   delta = dx*dx/16.0/Rowland_R;       /* gap between torus and center of CCD */
   r = Rowland_R - delta;	       /* R' in memo, and dx is L */
   delta_phi = 2.0*atan(dx/2.0/r)
     + 2.0*asin(CCD_Gap/(2.0*r*sqrt(1.0+0.25*(dx/r)*(dx/r))));

   /* The center angle is the angle of the center of the detector array with
    * respect to a grating on the rowland surface at the optical axis.  Let
    * phi0 denote the angle with respect to the center of the rowland circle.
    * Then:
    *   phi0 = Rowland_Theta + 2*t
    *   Center_CCD_Angle = Rowland_Theta + t
    * ==>
    *   phi0 = Rowland_Theta + 2*(Center_CCD_Angle-Rowland_Theta)
    *        = 2*Center_CCD_Angle - Rowland_Theta;
    */
   phi0 = 2.0*Center_CCD_Angle - Rowland_Theta;
   /* phi0 is now the desired center of the detector array.  Adjust it to
    * get the desired center of the first CCD.
    *
    *     O     | 1 | 2 | 3 |     (odd number)
    *     O   | 1 | 2 | 3 | 4 |   (even number)
    *
    * For an even number
    */
   phi0 -= ((Num_IXO_CCDs/2) - 0.5*((Num_IXO_CCDs+1) % 2))*delta_phi;

#if 0
   /* Setup the imaging CCD */
   d = geom + 0;
   d->x_ll = x_ll;
   d->x_lr = x_lr;
   d->x_ur = x_ur;
   d->x_ul = x_ul;
   d->x_pixel_size = X_Pixel_Size;
   d->y_pixel_size = Y_Pixel_Size;
   d->num_x_pixels = NUM_X_PIXELS;
   d->num_y_pixels = NUM_X_PIXELS;
   d->id = 0;
#endif

   /* Rowland center */
   cntr.x = Rowland_R * cos(Rowland_Theta);
   cntr.y = 0.0;
   cntr.z = Rowland_R * sin(Rowland_Theta);

   /* Now do the spectroscopic CCDs */
   axis.x = 0.0; axis.y = 1.0; axis.z = 0.0;
   for (i = 0; i < Num_IXO_CCDs; i++)
     {
	double phi, sin_phi, cos_phi;
	JDMVector_Type v;

	d = geom + i;

	phi = phi0 + i*delta_phi;
	sin_phi = sin(phi);
	cos_phi = cos(phi);

	d->x_ll = JDMv_rotate_vector1 (x_ll, axis, cos_phi, sin_phi);
	d->x_lr = JDMv_rotate_vector1 (x_lr, axis, cos_phi, sin_phi);
	d->x_ur = JDMv_rotate_vector1 (x_ur, axis, cos_phi, sin_phi);
	d->x_ul = JDMv_rotate_vector1 (x_ul, axis, cos_phi, sin_phi);

	/* Now push it to the center of the torus, and then to the surface */
	v.x = cntr.x-r*cos_phi; v.y = cntr.y+0.0; v.z = cntr.z+r*sin_phi;

	d->x_ll = JDMv_sum (d->x_ll, v);
	d->x_lr = JDMv_sum (d->x_lr, v);
	d->x_ul = JDMv_sum (d->x_ul, v);
	d->x_ur = JDMv_sum (d->x_ur, v);

	d->x_pixel_size = X_Pixel_Size;
	d->y_pixel_size = Y_Pixel_Size;
	d->num_x_pixels = Num_X_Pixels;
	d->num_y_pixels = Num_Y_Pixels;
	d->id = i;
     }

   return 0;
}

static char *QE_CCD_File;
static Param_Table_Type IXO_CCD_Parm_Table [] =
{
   {"IXO_Center_CCD_Angle",	PF_REAL_TYPE,	&Center_CCD_Angle},
   {"IXO_CCD_QE_File",		PF_FILE_TYPE,	&QE_CCD_File},
   {"IXO_CCD_Gain",		PF_REAL_TYPE,	&CCD_Gain},
   {"IXO_CCD_Fano",		PF_REAL_TYPE,	&Fano_Factor},
   {"IXO_CCD_Noise",		PF_REAL_TYPE,	&Read_Noise},
   {"IXO_Num_CCDs",		PF_UINT_TYPE,	&Num_IXO_CCDs},
   {"IXO_CCD_XPixel_Size",	PF_REAL_TYPE,	&X_Pixel_Size},
   {"IXO_CCD_YPixel_Size",	PF_REAL_TYPE,	&Y_Pixel_Size},
   {"IXO_CCD_Num_XPixels",	PF_UINT_TYPE,	&Num_X_Pixels},
   {"IXO_CCD_Num_YPixels",	PF_UINT_TYPE,	&Num_Y_Pixels},
   {NULL, 0, NULL}
};

/* p could be NULL */
static int read_ixo_ccd_parms (Param_File_Type *p)
{
   if (-1 == pf_get_parameters (p, IXO_CCD_Parm_Table))
     return -1;

   if ((Num_X_Pixels == 0) || (Num_Y_Pixels == 0)
       || (X_Pixel_Size <= 0) || (Y_Pixel_Size <= 0))
     {
	marx_error ("%s", "IXO CCD pixel parameters are invalid");
	return -1;
     }
   if (Num_IXO_CCDs > MAX_NUM_IXO_CCDS)
     {
	marx_error ("The number of IXO CCDs must be at most %u.", MAX_NUM_IXO_CCDS);
	return -1;
     }

   Center_CCD_Angle *= PI/180.0;
   CCD_Gain *= 0.001;
   return 0;
}

static Marx_Detector_Type IXO_CCD_Detector;
Marx_Detector_Type *_marx_get_ixo_ccd_detector (void)
{
   Marx_Detector_Type *d;

   d = &IXO_CCD_Detector;
   if (d->is_initialized)
     return d;

   _marx_catgs_init_variables ();

   d->detector_type = MARX_DETECTOR_IXO_CATGS_CCD;
   d->facet_list = _marx_link_detector_facet_list (IXO_CCD_List, Num_IXO_CCDs, sizeof(Marx_Detector_Geometry_Type));
   if (-1 == setup_geom (d->facet_list))
     return NULL;
   (void) _marx_compute_detector_basis (d);
   d->is_initialized = 1;
   return d;
}

int _marx_ixoccd_init (Param_File_Type *p)
{
   IXO_CCD_QE_Type *qeinfo;
   Marx_Detector_Type *d;
   Marx_Detector_Geometry_Type *g;

   if ((0 == CatGS_Init_Called)
       && (-1 == _marx_catgs_init (p)))/* need this initialized */
     return -1;

   if (-1 == read_ixo_ccd_parms (p))
     return -1;

   if (NULL == (d = _marx_get_ixo_ccd_detector ()))
     return -1;

   qeinfo = NULL;
   if (_Marx_Det_Ideal_Flag == 0)
     {
	char *file;

	if (_Marx_Det_Ideal_Flag == 0)
	  marx_message ("Reading IXO CCD QE/Filter Files\n");

	if (NULL == (file = marx_make_data_file_name (QE_CCD_File)))
	  return -1;

	marx_message ("\t%s\n", file);
	qeinfo = read_ccd_qe_file (file);
	marx_free (file);

	if (qeinfo == NULL)
	  return -1;
     }

   g = d->facet_list;
   while (g != NULL)
     {
	g->pha_fun = compute_pha;
	g->read_noise = Read_Noise;
	g->energy_gain = CCD_Gain;
	g->fano_factor = Fano_Factor;

	if (_Marx_Det_Ideal_Flag == 0)
	  {
	     g->qeinfo = qeinfo;
	     if (g != d->facet_list)
	       qeinfo->num_refs++;/* ugly */
	  }
	g = g->next;
     }

   return 0;
}

