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
#include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"

/* This structure is used for the QE curves */
typedef struct
{
   unsigned int n_energies;
   float *energies;
   float *eff;
}
Detector_Type;

static Marx_Detector_Geometry_Type *ACIS_I_Chips;
static _Marx_Acis_Chip_Type Acis_CCDS [_MARX_NUM_ACIS_I_CHIPS];

#if MARX_HAS_ACIS_STREAK
double Frame_Time;
double Exposure_Time;
double Frame_Transfer_Time;

void _marx_acis_apply_streak (double tstart,
			      Marx_Photon_Attr_Type *at,
			      Marx_Detector_Geometry_Type *d)
{
   double t;
   JDMVector_Type dx;
   double xpixel, ypixel;

   if (Frame_Transfer_Time <= 0.0)
     return;

   t = tstart + at->arrival_time;
   t = fmod (t, Frame_Time);
   if (t <= Exposure_Time)
     return;

   /* Border pixels are excluded because of ACIS event detection */
   at->z_pixel = 1.0 + 1022.0 * JDMrandom ();

   /* Now I need to update the X and P vectors.
    * Define the P vector with respect to the mirror nodal position.
    */

   xpixel = (at->y_pixel - d->xpixel_offset) * d->x_pixel_size;
   ypixel = (at->z_pixel - d->ypixel_offset) * d->y_pixel_size;

   dx = JDMv_ax1_bx2 (xpixel, d->xhat, ypixel, d->yhat);
   at->x = JDMv_sum (d->x_ll, dx);

   at->p = JDMv_diff (at->x, JDMv_vector (Marx_Focal_Length, 0, 0));
   JDMv_normalize (&at->p);
   at->flags |= PHOTON_ACIS_STREAKED;
}

#endif

int _marx_acis_i_detect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *at, *attrs;
   unsigned int n_photons, i;
   unsigned int *sorted_index;
#if MARX_HAS_ACIS_STREAK
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

#if MARX_HAS_ACIS_STREAK
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
	_marx_transform_ray (&at->x, &at->p, &_Marx_Det_XForm_Matrix);

	d = _marx_intersect_with_detector (at->x, at->p,
					   ACIS_I_Chips,
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
  	     _marx_acis_enyz_to_grade(at);
	     pt->history |= (MARX_DET_GRADE_OK);

	     if (0 == _marx_acis_apply_qe_and_pha (&Acis_CCDS[(unsigned int) (d->id)], at))
	       {
#if MARX_HAS_ACIS_STREAK
		  (void) _marx_acis_apply_streak (tstart, at, d);
#endif
	       }

	  }

	/* Transform detected photon back to original system */
	_marx_transform_ray_reverse (&at->x, &at->p, &_Marx_Det_XForm_Matrix);

#if MARX_HAS_DITHER
	_marx_undither_detector (&at->dither_state);
#endif
     }

   return 0;
}

/*}}}*/

static Param_Table_Type ACIS_Parm_Table [] =
{
#if MARX_HAS_ACIS_STREAK
     {"ACIS_Exposure_Time",	PF_REAL_TYPE,	&Exposure_Time},
     {"ACIS_Frame_Transfer_Time",PF_REAL_TYPE,	&Frame_Transfer_Time},
#endif
     {NULL, 0, NULL}
};

int _marx_acis_get_generic_parms (Param_File_Type *pf)
{
   if (-1 == pf_get_parameters (pf, ACIS_Parm_Table))
     return -1;

#if MARX_HAS_ACIS_STREAK
   if (Exposure_Time <= 0)
     {
	marx_error ("ACIS_Exposure_Time must be greater than 0");
	return -1;
     }
   if (Frame_Transfer_Time < 0.0)
     Frame_Transfer_Time = 0.0;

   Frame_Time = Frame_Transfer_Time + Exposure_Time;
#endif

   return 0;
}

static Param_Table_Type Parm_Table [] =
{
     {NULL, 0, NULL}
};

#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
static int get_double_param (Param_File_Type *pf, char *fmt, int i, double *v)
{
   char parm[128];

   sprintf (parm, fmt, i);
   if (-1 == pf_get_double (pf, parm, v))
     {
	marx_error ("Unable to get paramter %s", parm);
	return -1;
     }

   return 0;
}
#endif
static int get_file_param (Param_File_Type *pf, char *fmt, int i, char **v)
{
   char parm[128];
   char file [PF_MAX_LINE_LEN];

   sprintf (parm, fmt, i);
   if (-1 == pf_get_file (pf, parm, file, sizeof (file)))
     {
	marx_error ("Unable to get paramter %s", parm);
	return -1;
     }

   if (NULL == (*v = (char *)marx_malloc (strlen (file) + 1)))
     return -1;

   strcpy (*v, file);
   return 0;
}

static int get_acis_parms (Param_File_Type *p)
{
   int i;

   if (-1 == _marx_acis_get_generic_parms (p))
     return -1;

   if (-1 == pf_get_parameters (p, Parm_Table))
     return -1;

   for (i = 0; i < _MARX_NUM_ACIS_I_CHIPS; i++)
     {
	_Marx_Acis_Chip_Type *ccd;

	ccd = &Acis_CCDS[i];

	ccd->ccd_id = i;
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
	if (-1 == get_double_param (p, "ACIS_CCD%d_Gain", i, &ccd->ccd_gain))
	  return -1;
	if (ccd->ccd_gain <= 0.0)
	  ccd->ccd_gain = 4;
	ccd->ccd_gain /= 1000.0;	       /* Convert to KeV */

	if (-1 == get_double_param (p, "ACIS_CCD%d_Offset", i, &ccd->ccd_offset))
	  return -1;
	ccd->ccd_offset /= 1000.0;	       /* Convert to KeV */

	if (-1 == get_double_param (p, "ACIS-I%d-FanoFactor", i, &ccd->fano_factor))
	  return -1;

	if (-1 == get_double_param (p, "ACIS-I%d-ReadNoise", i, &ccd->read_noise))
	  return -1;

	if (-1 == get_double_param (p, "ACIS-I%d-EnergyGain", i, &ccd->energy_gain))
	  return -1;
#endif
	if (-1 == get_file_param (p, "ACIS-I%d-QEFile", i, &ccd->qe_file))
	  return -1;

	if (-1 == get_file_param (p, "ACIS-I%d-FilterFile", i, &ccd->filter_file))
	  return -1;

     }

   return 0;
}

#if MARX_HAS_ACIS_FEF
static short apply_fef (_Marx_Acis_Chip_Type *c, float x, float y, double en, float *pi)
{
   short pha;

   if (-1 == _marx_apply_acis_rmf (c, x, y, en, pi, &pha))
     {
	pha = -1;
	*pi = 0;
     }
   return pha;
}
#endif

int _marx_acis_i_init (Param_File_Type *p) /*{{{*/
{
   Marx_Detector_Type *acis_i;
   unsigned int i;
   int igrfits[2];

#if MARX_HAS_ACIS_FEF
   if (-1 == marx_init_acis_i_rmf (p))
     return -1;
#endif

   for (i = 0; i < _MARX_NUM_ACIS_I_CHIPS; i++)
     _marx_free_acis_chip_type (&Acis_CCDS[i]);

   if (-1 == get_acis_parms (p))
     return -1;

   if (NULL == (acis_i = marx_get_detector_info ("ACIS-I")))
     return -1;

   ACIS_I_Chips = acis_i->facet_list;

#if MARX_HAS_ACIS_GAIN_MAP
   if (-1 == _marx_init_acis_i_gain_map (p))
     return -1;
#endif

   if (_Marx_Det_Ideal_Flag == 0)
     marx_message ("Reading ACIS-I QE/Filter Files\n");

   for (i = 0; i < _MARX_NUM_ACIS_I_CHIPS; i++)
     {
	_Marx_Acis_Chip_Type *ccd = &Acis_CCDS[i];
#if MARX_HAS_ACIS_GAIN_MAP
	ccd->pha_fun = _marx_apply_acis_gain_map;
#else
# if MARX_HAS_ACIS_FEF
	ccd->pha_fun = apply_fef;
# else
	ccd->pha_fun = _marx_acis_compute_fs_pha;
# endif
#endif
	if (_Marx_Det_Ideal_Flag == 0)
	  {
	     if (-1 == _marx_acis_read_chip_efficiencies (ccd))
	       return -1;
	     if (-1 == _marx_acis_contam_init (p, ccd))
	       return -1;
	  }
     }

   /* AML October 11, 2018
    *
    * if grades FITS files have not been read in, read them now.
    *
    *
    */
   set_print_grfits_diag();
   // unset_print_grfits_diag();
   read_acis_grades_files(igrfits);

   return 0;
}

