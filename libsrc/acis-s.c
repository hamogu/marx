/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2002-2015 Massachusetts Institute of Technology

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

/*{{{ Include Files */

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

/*}}}*/

#define USE_CALDB_FILES

static _Marx_Acis_Chip_Type Acis_CCDS [_MARX_NUM_ACIS_S_CHIPS];
static Marx_Detector_Geometry_Type *ACIS_S_Chips;
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
static char *FsBs_Configuration;
#endif
static Param_Table_Type ACIS_Parm_Table [] = /*{{{*/
{
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
     {"ACIS-S-FsBsConf",	PF_STRING_TYPE, &FsBs_Configuration},
#endif
   {NULL, 0, NULL}
};

/*}}}*/

#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
short _marx_acis_compute_fs_pha (_Marx_Acis_Chip_Type *ccd, float x, float y,
				 double energy, float *pi) /*{{{*/
{
   double da;
   double gain = ccd->energy_gain;
   double noise = ccd->read_noise;
   short pha;

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

   pha = (short) ((energy - ccd->ccd_offset)/ ccd->ccd_gain);
   if (pha < 0)
     pha = 0;

   *pi = (float) energy;
   return pha;
}

/*}}}*/

short _marx_acis_compute_bs_pha (_Marx_Acis_Chip_Type *ccd, float x, float y,
				 double energy, float *pi) /*{{{*/
{
   double da;
   double gain = ccd->energy_gain;
   double noise = ccd->read_noise;
   short pha;

   (void) x;
   (void) y;

   if (gain == 0.0)
     {
	*pi = 0.0;
	return 0.0;
     }

   /* For backside chips, there is no simple value for the energy->pha
    * mapping.  Gregory P. suggests the following hack:
    */
   da = energy;
   if (da < 1.0) da = 1.0;

   /* Eq 2.1 of ACIS-PSU-SOP-01 suggests the following: */
   da = gain * sqrt (noise * noise + ccd->fano_factor * da / gain);

   energy = energy + da * JDMgaussian_random ();
   if (energy < 0.0) energy = 0.0;

   pha = (short) ((energy - ccd->ccd_offset)/ ccd->ccd_gain);
   if (pha < 0)
     pha = 0;

   *pi = (float) energy;
   return pha;
}

/*}}}*/
#endif

int _marx_acis_apply_qe_and_pha (_Marx_Acis_Chip_Type *ccd, Marx_Photon_Attr_Type *at)
{
   double qe, qe_filter, qe_contam;
   double r, en;

   en = at->energy;

   if (_Marx_Det_Ideal_Flag == 0)
     {
	r = JDMrandom ();

	if (ccd->qe_num_energies != 0)
	  qe = JDMinterpolate_f (en, ccd->qe_energies, ccd->qe, ccd->qe_num_energies);
	else
	  qe = 1.0;

	if (ccd->filter_num_energies != 0)
	  qe_filter = JDMinterpolate_f (en, ccd->filter_energies, ccd->filter_qe, ccd->filter_num_energies);
	else
	  qe_filter = 1.0;

	qe_contam = (*ccd->contam_fun)(ccd, en, at->y_pixel, at->z_pixel);

	if (r >= qe * qe_filter * qe_contam)
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

int _marx_acis_s_detect (Marx_Photon_Type *pt) /*{{{*/
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
	_marx_transform_ray (&at->x, &at->p,
			     &_Marx_Det_XForm_Matrix);

	/* See if the photon will hit the CCD and if so, which one. */
	d = _marx_intersect_with_detector (at->x, at->p,
					   ACIS_S_Chips,
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

	     if (0 == _marx_acis_apply_qe_and_pha (&Acis_CCDS[(unsigned int) (d->id - 4)], at))
	       {
#if MARX_HAS_ACIS_STREAK
		  (void) _marx_acis_apply_streak (tstart, at, d);
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

/*}}}*/

void _marx_free_acis_chip_type (_Marx_Acis_Chip_Type *c)
{
   if (c == NULL)
     return;

   JDMfree_float_vector (c->qe_energies);
   JDMfree_float_vector (c->qe);
   JDMfree_float_vector (c->filter_qe);
   JDMfree_float_vector (c->filter_energies);

   marx_free (c->qe_file);
   marx_free (c->filter_file);

   memset ((char *) c, 0, sizeof (_Marx_Acis_Chip_Type));
}

#ifndef USE_CALDB_FILES
static int read_efficiency_file (char *file,
				 float **en, float **eff, unsigned int *num)
{
   if ((file == NULL) || (*file == 0))
     {
	*num = 0;
	*eff = NULL;
	*en = NULL;
	return 0;
     }

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;

   marx_message ("\t%s\n", file);

   if (-1 == marx_f_read_bdat (file, num, 2, en, eff))
     {
	*num = 0;
	*eff = NULL;
	*en = NULL;
	marx_free (file);
	return -1;
     }

   marx_free (file);
   return 0;
}
#endif

int _marx_acis_read_chip_efficiencies (_Marx_Acis_Chip_Type *chip)
{
   if (chip == NULL)
     return 0;
#ifdef USE_CALDB_FILES
   return _marx_read_acis_qe (chip->ccd_id, &chip->qe_energies, &chip->qe, &chip->qe_num_energies);
#else
   if (-1 == read_efficiency_file (chip->qe_file,
				   &chip->qe_energies,
				   &chip->qe,
				   &chip->qe_num_energies))
     return -1;

   if (-1 == read_efficiency_file (chip->filter_file,
				   &chip->filter_energies,
				   &chip->filter_qe,
				   &chip->filter_num_energies))
     return -1;
   return 0;
#endif
}

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

   if (-1 == pf_get_parameters (p, ACIS_Parm_Table))
     return -1;

   for (i = 0; i < _MARX_NUM_ACIS_S_CHIPS; i++)
     {
	_Marx_Acis_Chip_Type *ccd;

	ccd = &Acis_CCDS[i];

	ccd->ccd_id = i+4;
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
	if (-1 == get_double_param (p, "ACIS_CCD%d_Gain", i+4, &ccd->ccd_gain))
	  return -1;
	if (ccd->ccd_gain <= 0.0)
	  ccd->ccd_gain = 4;
	ccd->ccd_gain /= 1000.0;	       /* Convert to KeV */

	if (-1 == get_double_param (p, "ACIS_CCD%d_Offset", i+4, &ccd->ccd_offset))
	  return -1;
	ccd->ccd_offset /= 1000.0;	       /* Convert to KeV */

	if (-1 == get_double_param (p, "ACIS-S%d-FanoFactor", i, &ccd->fano_factor))
	  return -1;

	if (-1 == get_double_param (p, "ACIS-S%d-ReadNoise", i, &ccd->read_noise))
	  return -1;

	if (-1 == get_double_param (p, "ACIS-S%d-EnergyGain", i, &ccd->energy_gain))
	  return -1;
#endif
	if (-1 == get_file_param (p, "ACIS-S%d-QEFile", i, &ccd->qe_file))
	  return -1;

	if (-1 == get_file_param (p, "ACIS-S%d-FilterFile", i, &ccd->filter_file))
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

int _marx_acis_s_init (Param_File_Type *p) /*{{{*/
{
   unsigned int i;
   Marx_Detector_Type *acis_s;

#if MARX_HAS_ACIS_FEF
   if (-1 == marx_init_acis_s_rmf (p))
     return -1;
#endif

   for (i = 0; i < _MARX_NUM_ACIS_S_CHIPS; i++)
     _marx_free_acis_chip_type (&Acis_CCDS[i]);

   if (-1 == get_acis_parms (p))
     return -1;

   if (NULL == (acis_s = marx_get_detector_info ("ACIS-S")))
     return -1;

   ACIS_S_Chips = acis_s->facet_list;
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
   /* Check FsBs configuration */
   if (strlen (FsBs_Configuration) != _MARX_NUM_ACIS_S_CHIPS)
     {
	marx_error ("The FsBs Configuration must be 6 characters for 6 chips");
	return -1;
     }
#endif

#if MARX_HAS_ACIS_GAIN_MAP
   if (-1 == _marx_init_acis_s_gain_map (p))
     return -1;
#endif
   if (_Marx_Det_Ideal_Flag == 0)
     marx_message ("Reading ACIS-S QE/Filter Files\n");

   for (i = 0; i < _MARX_NUM_ACIS_S_CHIPS; i++)
     {
	_Marx_Acis_Chip_Type *ccd = &Acis_CCDS[i];

#if MARX_HAS_ACIS_FEF
	ccd->pha_fun = apply_fef;
#else
# if MARX_HAS_ACIS_GAIN_MAP
	ccd->pha_fun = _marx_apply_acis_gain_map;
# else
	switch (FsBs_Configuration[i])
	  {
	   case 'F':
	   case 'f':
	     ccd->pha_fun = _marx_acis_compute_fs_pha;
	     break;

	   case 'b':
	   case 'B':
	     ccd->pha_fun = _marx_acis_compute_bs_pha;
	     break;

	   default:
	     marx_error ("FsBs Configuration character must be 'b' or 'f'");
	     return -1;
	  }
# endif				       /* MARX_HAS_ACIS_GAIN_MAP */
#endif				       /* MARX_HAS_ACIS_FEF */
	if (_Marx_Det_Ideal_Flag == 0)
	  {
	     if (-1 == _marx_acis_contam_init (p, ccd))
	       return -1;

	     if (-1 == _marx_acis_read_chip_efficiencies (ccd))
	       return -1;
	  }
     }

   return 0;
}

/*}}}*/

