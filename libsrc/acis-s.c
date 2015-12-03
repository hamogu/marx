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

/* grades5.h
 * Alan M. Levine
 * October 10, 2018
 * Heritage: grades4.h
 *
 * Header file for grades5.c
 */


extern int init_gr_p_tab_ptr(int isidecode);
extern void set_image_array_indices(int jgr, int jen, int jx, int jy);
extern void get_grprdim(int grlims[4]);
extern int check_gr_pr_indices();
extern int get_gr_prob(double *val);
extern void print_gr_prob(int nc, FILE *fpo);
extern void check_gr_probs(int nc, double maxdel, int ipr, FILE *fpo);

extern int get_nebins();
extern int get_ebin(double ephot);
extern void print_en_bins(FILE *fpo);
extern void check_gr_en_range();

extern int get_ngrades();
extern int get_grade(int jgr);
extern void print_grade_list(FILE *fpo);
extern int read_grades_fits_file(FILE *gradesfp);
extern double test_random();
extern int get_random_grade(int jen, int jx, int jy);
extern int random_grade(double ephot, double xpixsz, double xpos, double ypixsz,
			double ypos);
extern int acis_grades_init(FILE *fpin);
extern int fsbscode[10];
extern int j_read_fsbs_fits[2];

// fs = 1, bs = 2
int fsbscode[10] = { 1, 1, 1, 1, 1, 2, 1, 2, 1, 1 };
int j_read_fsbs_fits[2] = { 0, 0 };
static int num_grade;
char *grades_data_directory;
char grades_input_file[512];

int print_grfits_diag, extra_stderr_diag;
FILE *fpgrdiag;

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

/****************************************************************************/
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
/****************************************************************************/

/*}}}*/

/****************************************************************************/
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
/****************************************************************************/
/*}}}*/
#endif

/****************************************************************************/
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

/****************************************************************************/
// The following function, _marx_acis_yzen_to_grade (), started by AML April 12, 2018
int _marx_acis_enyz_to_grade (Marx_Photon_Attr_Type *at)
{
  double y, z, en;
  int ccd_num, ipval, ifsbs;
  short gr;

  y = at->y_pixel;
  z = at->z_pixel;
  en = at->energy;
  ccd_num = at->ccd_num;
  ifsbs = fsbscode[ccd_num];

  ipval = init_gr_p_tab_ptr(ifsbs);
  gr = random_grade(en,1.0,y,1.0,z);
  at->grade = gr;
  if(gr >= 0) {
    num_grade++;
  }
  if( (extra_stderr_diag == 1) && (num_grade%100 == 0) ) {
    fprintf(stderr,"grade gr = %d\n",gr);
    return(1);
  } 
  if(gr < 0) {
    fprintf(stderr,"******grade gr = %d\n",gr);
    return(0);
  }
}

/****************************************************************************/
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
  	     _marx_acis_enyz_to_grade(at);
	     pt->history |= (MARX_DET_GRADE_OK);

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

/****************************************************************************/
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
/****************************************************************************/
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

/****************************************************************************/
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
/****************************************************************************/

#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
/****************************************************************************/
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
/****************************************************************************/
#endif

/****************************************************************************/
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

/****************************************************************************/
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

/****************************************************************************/
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
/****************************************************************************/

int _marx_acis_s_init (Param_File_Type *p) /*{{{*/
{
   unsigned int i;
   int igrfits[2];
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

   /* AML October 11, 2018
    *
    * if grades FITS files have not been read in, read them now.
    *
    *
    */
   // set_print_grfits_diag();
   unset_print_grfits_diag();
   unset_extra_stderr_diag();
   read_acis_grades_files(igrfits);

   return 0;
}
/****************************************************************************/

/*}}}*/

/****************************************************************************/
/* grades5.c
 * Alan M. Levine
 * October 10, 2018
 * Heritage: grades4.c
 *
 * Cheap FITS file reader.
 * This is not a general code in any sense.
 * Read a Marx event grades definition FITS file.
 */

/* Requirements on FITS input file:
 *
 * 1) Three records in the following order:
 *      (i) primary image with 4-d grade probability table data (doubles)
 *      (ii) extension with two-column binary table containing the energy bin limits (doubles)
 *      iii) extension (IMAGE) with 1-d grade reindex table (16, 32, or 64 bit integers).
 */

/*
 * October 10, 2018 - Enable two or more grade probability tables.  For two tables,
 *                    one would apply to front-side illuminated chips and the other to 
 *                    back-side chips.
 */

/*
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "grades4.h"
*/

#define GR_ELOW_MAX 0.2   // keV
#define GR_EHI_MIN 10.0   // keV

#define GR5_MAX_N_HDR 6
#define GR5_MAX_N_AX 10
#define GR5_BLKSIZ 2880   // no. of bytes in a FITS file block
#define GR5_MAXCOL 10     // maximum no.of columns in a binary table

// extern double JDMrandom (void);

char hdrstr[80], ttype[GR5_MAXCOL][64], tform[GR5_MAXCOL][64];
/*
   naxes - value of NAXIS
   nax[] - values of NAXISnum
   bitp - value of BITPIX
   iext - =1 if 'EXTEND' value is 'T', =0 otherwise

   ihdr = -1 - not set
        =  0 - primary
        =  1 - image
        =  2 - ASCII table
        =  3 - binary table

   nflds - value of 'TFIELDS' keyword
   numbyt[] - number of bytes in a data item 
   intflt[] - =1 for integer types, =2 for floats or doubles
   ncol -  count of columns in a binary table; after a header is read it
           should be equal to nflds
   irec = record number (starts with 1)
*/
int naxes, nax[GR5_MAX_N_AX], bitp, iext, ihdr, nflds, ncol, irec, iecol[2];
int numbyt[GR5_MAXCOL], intflt[GR5_MAXCOL];
short pshrt;
int pint;
long long psll;
float pflt;
double pdbl;

#define GR5_MAX_GRADES 256
#define GR5_MAX_EBINS 101
#define GR5_MAX_SUBPIX 10   // applies to both x and y

// grprdim[]: max index in each dimension of grade probability table = gr_prob
// order = grade, energy bin, x, y
int bpgrt;
int igr, ien, ix, iy;

// 0 => front side
// 1 => back side
struct grtableset {
  int grprdim[4], grprmod[4], nebins, ngrades;
  int gr_table[GR5_MAX_GRADES];
  int ilim[4][2];
  double *gr_prob, ebin[GR5_MAX_EBINS][2];
} gr_p_tables[2];
struct grtableset *grt_ptr;

/****************************************************************************/

void set_print_grfits_diag()
{
  print_grfits_diag = 1;
  fpgrdiag = fopen("grfits_diag.txt","w");
  if(fpgrdiag == 0) {
    fprintf(stderr,"Fatal error opening grfits_diag.txt for writing.  Exiting.\n");
    exit(-10);
  }
}

/****************************************************************************/

void unset_print_grfits_diag()
{
  print_grfits_diag = 0;
  if(fpgrdiag != 0) {
    fclose(fpgrdiag);
  }
}

/****************************************************************************/

void set_extra_stderr_diag()
{
  extra_stderr_diag = 1;
}

/****************************************************************************/

void unset_extra_stderr_diag()
{
  extra_stderr_diag = 0;
}

/****************************************************************************/
// Run this function before the other functions below.
// It needs to be rerun every time the type of chip, i.e., front or back, changes.
// return value = isidecode if valid, otherwise -1
int init_gr_p_tab_ptr(int isidecode)
{
  int iret;

  if( (isidecode == 1) || (isidecode == 2) ) {
    grt_ptr = gr_p_tables + (isidecode - 1);
    iret = isidecode;
  }
  else {
    iret = -1;
  }
  return(iret);
}

/****************************************************************************/

void set_image_array_indices(int jgr, int jen, int jx, int jy)
{
  igr = jgr;
  ien = jen;
  ix = jx;
  iy = jy;
}
/****************************************************************************/

void get_grprdim(int grlims[4])
{
  int i;

  for(i=0;i<4;++i) {
    grlims[i] = grt_ptr->grprdim[i];
  }
}
/****************************************************************************/

int check_grprdim()
{
  int grdim[4], i;

  for(i=0;i<4;++i)
    grdim[i] = grt_ptr->grprdim[i];

  if( (grdim[0] < 1) || (grdim[0] > GR5_MAX_GRADES) )
    return(-1);
  if( (grdim[1] < 1) || (grdim[1] > GR5_MAX_EBINS) )
    return(-2);
  if( (grdim[2] < 1) || (grdim[2] > GR5_MAX_SUBPIX) )
    return(-3);
  if( (grdim[3] < 1) || (grdim[3] > GR5_MAX_SUBPIX) )
    return(-4);
  return(1);
}
/****************************************************************************/

int get_image_array_bin(int irec, int ibin)
{
  int ib;

  if(irec==0) {    // grade probability table
    igr = ibin % grt_ptr->grprdim[0];
    ien = (ibin/grt_ptr->grprmod[1]) % grt_ptr->grprdim[1];
    ix = (ibin/grt_ptr->grprmod[2]) % grt_ptr->grprdim[2];
    iy = (ibin/grt_ptr->grprmod[3]) % grt_ptr->grprdim[3];
    ib = igr + (ien + (ix + (iy*GR5_MAX_SUBPIX))*GR5_MAX_EBINS)*GR5_MAX_GRADES;
    // fprintf(stderr,"giab(): ibin,igr,ien,ix,iy = %d %d %d %d %d\n",ibin,igr,ien,ix,iy);
    return(ib);
  }
  else if(irec==2) {
  }
  return(-1);
}
/****************************************************************************/

int check_gr_pr_indices()
{

  if( (igr < 0) || (igr >= grt_ptr->grprdim[0]) )
    return(-1);
  if( (ien < 0) || (ien >= grt_ptr->grprdim[1]) )
    return(-2);
  if( (ix < 0) || (ix >= grt_ptr->grprdim[2]) )
    return(-3);
  if( (iy < 0) || (iy >= grt_ptr->grprdim[3]) )
    return(-4);
  return(1);
}

/****************************************************************************/
/* gr_prob is to be used as if it were an array
 *     array[MAX_SUBPIX][MAX_SUBPIX]{MAX_EBINS][MAX_GRADES]
 * where the last index varies most rapidly, and the first index indicates 
 * the 'y' subpixel.
 */

int put_gr_prob(double val)
{
  int ib;

  if(check_gr_pr_indices() < 0) {
    fprintf(stderr,"ERROR: A grade prob. table index is out of range.\n");
    fprintf(stderr,"  The indices are igr=%d, ien=%d, ix=%d, iy=%d\n",igr,ien,ix,iy);
    fprintf(stderr,"This error is fatal.  Exiting.\n");
    exit(-11);
  }
    
  // ib = (igr + ien*GR5_MAX_GRADES) + ix*GR5_MAX_GRADES*GR5_MAX_EBINS + 
  //       iy*GR5_MAX_GRADES*GR5_MAX_EBINS*GR5_MAX_SUBPIX;
  ib = igr + (ien + (ix + (iy*GR5_MAX_SUBPIX))*GR5_MAX_EBINS)*GR5_MAX_GRADES;
  *(grt_ptr->gr_prob + ib) = val;
  return(1);
}

/****************************************************************************/

int get_gr_prob(double *val)
{
  int ib;

  if(check_gr_pr_indices()<0) {
    fprintf(stderr,"ERROR: A grade prob. table index is out of range.\n");
    fprintf(stderr,"  The indices are igr=%d, ien=%d, ix=%d, iy=%d\n",igr,ien,ix,iy);
    fprintf(stderr,"This error is fatal.  Exiting.\n");
    exit(-12);
  }
    
  ib = igr + (ien + (ix + (iy*GR5_MAX_SUBPIX))*GR5_MAX_EBINS)*GR5_MAX_GRADES;
  *val = *(grt_ptr->gr_prob + ib);
  return(1);
}

/****************************************************************************/

void print_gr_prob(int nc, FILE *fpo)
{
  double val, gt;
  int j;

  for(iy=grt_ptr->ilim[3][0];iy<grt_ptr->ilim[3][1];++iy) {
    for(ix=grt_ptr->ilim[2][0];ix<grt_ptr->ilim[2][1];++ix) {
      for(ien=grt_ptr->ilim[1][0];ien<grt_ptr->ilim[1][1];++ien) {
	fprintf(fpo,"iy,ix,ien = %d %d %d\n",iy,ix,ien);
	gt = 0.0;
	for(igr=grt_ptr->ilim[0][0];igr<grt_ptr->ilim[0][1];++igr) {
	  if(igr%nc == 0) {
	    if(igr>0)
	      fprintf(fpo,"\n");
	    fprintf(fpo,"%d  ",igr);
	  }
	  j = get_gr_prob(&val);
	  fprintf(fpo," %8.5f",val);
	  gt += val;
	}
	fprintf(fpo,"\n");
	fprintf(fpo,"subtotal = %f\n",gt);
      }
    }
  }

}

/****************************************************************************/

void check_gr_probs(int nc, double maxdel, int ipr, FILE *fpo)
{
  double val, gt, newval;
  int j, k;

  for(iy=grt_ptr->ilim[3][0];iy<grt_ptr->ilim[3][1];++iy) {
    for(ix=grt_ptr->ilim[2][0];ix<grt_ptr->ilim[2][1];++ix) {
      for(ien=grt_ptr->ilim[1][0];ien<grt_ptr->ilim[1][1];++ien) {
	if(ipr==1) fprintf(fpo,"iy,ix,ien = %d %d %d\n",iy,ix,ien);
	gt = 0.0;
	for(igr=grt_ptr->ilim[0][0];igr<grt_ptr->ilim[0][1];++igr) {
	  if(igr%nc == 0) {
	    if(igr>0)
	      if(ipr==1) fprintf(fpo,"\n");
	    if(ipr==1) fprintf(fpo,"%d  ",igr);
	  }
	  j = get_gr_prob(&val);
	  if(ipr==1) fprintf(fpo," %8.5f",val);
	  gt += val;
	}
	if(ipr==1) fprintf(fpo,"\n");
	if(ipr==1) fprintf(fpo,"subtotal = %f\n",gt);

	if(fabs(gt - 1.0) > maxdel) {
	  fprintf(stderr,"Total grade probability is outside allowed range.  Will exit.\n");
	  fprintf(stderr,"Maximum deviation of the subtotal from 1.0 is %e.\n",maxdel);
	  fprintf(stderr,"For ien, ix, iy = %d %d %d, subtotal prob = %e\n",
		  ien,ix,iy,gt);
	  exit(-13);
	}
	else if(gt != 1.0) {
	  // renormalize
	  for(igr=grt_ptr->ilim[0][0];igr<grt_ptr->ilim[0][1];++igr) {
	    j = get_gr_prob(&val);
	    newval = val/gt;
	    k = put_gr_prob(newval);
	  }
	}

      }  // for(ien=...
    }   // for(ix=...
  } // for(iy=...

}

/****************************************************************************/

void alloc_gradepr()
{
  int nb;

  nb = GR5_MAX_GRADES*GR5_MAX_EBINS*GR5_MAX_SUBPIX*GR5_MAX_SUBPIX;
  grt_ptr->gr_prob = calloc(nb,sizeof(double));
}

/****************************************************************************/

void set_e_columns()
{
  if(nflds != 2) {
    fprintf(stderr,"set_e_columns(): ERROR - nflds = %d is not 2, exiting...\n",nflds);
    exit(-1);
  }

  if( (!strncmp("ENERG_LO",&ttype[0][1],8)) && (!strncmp("ENERG_HI",&ttype[1][1],8)) ) {
    iecol[0] = 0;
    iecol[1] = 1;
  }
  else if( (!strncmp("ENERG_HI",&ttype[0][1],8)) && (!strncmp("ENERG_LO",&ttype[1][1],8)) ) {
    iecol[0] = 1;
    iecol[1] = 0;
  }
  else {
    fprintf(stderr,"set_e_columns(): ERROR - ENERG_LO and ENERG_HI misplaced; exiting...\n");
  }
}

/****************************************************************************/

void print_en_bins(FILE *fpo)
{
  int ie;

  for(ie=0;ie<grt_ptr->nebins;++ie) {
    fprintf(fpo,"%d %f %f\n",ie,grt_ptr->ebin[ie][0],grt_ptr->ebin[ie][1]);
  }
}

/****************************************************************************/

int get_nebins()
{
  return(grt_ptr->nebins);
}
/****************************************************************************/
// Given a photon energy, return the energy bin number.
// Retun iebin = -1 if the photon energy is outside the valid range.
// #define MAX_EBINS 101
// double *gr_prob, ebin[MAX_EBINS][2];

int get_ebin(double ephot)
{
  int i, iebin;

  iebin = -1;
  for(i=0;i<grt_ptr->nebins;++i) {
    if( (ephot >= grt_ptr->ebin[i][0]) && (ephot < grt_ptr->ebin[i][1]) ) {
      iebin = i;
      break;
    }
  }
  if(iebin < 0) {
    fprintf(stderr,"get_ebin() error - energy %f is out of range\n",ephot);
    fprintf(stderr,"  returning -1\n");
  }
  return(iebin);
}
/****************************************************************************/

void check_gr_en_range()
{
  double elo, ehi;

  elo = grt_ptr->ebin[0][0];
  ehi = grt_ptr->ebin[grt_ptr->nebins - 1][1];

  if( (elo > GR_ELOW_MAX) || (ehi < GR_EHI_MIN) ) {
    fprintf(stderr,"The grade probability table energy range is too small.\n");
    fprintf(stderr,"The table's lowest energy (%f keV) may not exceed %f keV\n",
	    elo,GR_ELOW_MAX);
    fprintf(stderr,"The table's highest energy (%f keV) may not be lower than %f keV\n",
	    ehi,GR_EHI_MIN);
    fprintf(stderr,"EXITING.\n");
    exit(-14);
  }
}

/****************************************************************************/

int get_ngrades()
{
  return(grt_ptr->ngrades);
}
/****************************************************************************/
// get grade from grade table
// int bpgrt, gr_table[MAX_GRADES];

int get_grade(int jgr)
{

  if( (jgr >= 0) && (jgr < grt_ptr->ngrades) ) {
    return(grt_ptr->gr_table[jgr]);
  }
  else {
    fprintf(stderr,"get_grade() error - index %d is out of range\n",jgr);
    fprintf(stderr,"  returning -1\n");
    return(-1);
  }
}

/****************************************************************************/

void print_grade_list(FILE *fpo)
{
  int ig;

  for(ig=0;ig<grt_ptr->ngrades;++ig) {
    fprintf(fpo,"%d %d\n",ig,grt_ptr->gr_table[ig]);
  }
}
/****************************************************************************/

double test_random()
{
  return(JDMrandom());  // obtain a random number uniformly distr. on 0 to 1
}
/****************************************************************************/
// get a random grade for specified ien, ix, and iy.
// int bpgrt, gr_table[MAX_GRADES];

int get_random_grade(int jen, int jx, int jy)
{
  double drand, cump, val;
  int jgr, jfl, igrade;

  drand = JDMrandom();  // obtain a random number uniformly distr. on 0 to 1

  cump = 0.0;
  for(jgr=0;jgr<grt_ptr->grprdim[0];++jgr) {
    set_image_array_indices(jgr,jen,jx,jy);
    if( (jfl = get_gr_prob(&val)) == 1) {
      if( (drand >= cump) && (drand <= (cump + val)) ) {
	igrade = get_grade(jgr);
	return(igrade);
      }
      cump += val;
    }
    else {
      fprintf(stderr,"get_random_grade() error return from get_gr_prob()\n");
      return(-1);
    }
  }
  return(-2);
}

/****************************************************************************/

int random_grade(double ephot, double xpixsz, double xpos, double ypixsz,
		 double ypos)
{
  double xp, yp;
  int jx, jy, jen;

  jen = get_ebin(ephot);    // add error check !!!

  // Check that xpixsz and ypixsz are > 0 ?

  xp = (xpos/xpixsz) - floor(xpos/xpixsz);
  jx = xp*grt_ptr->grprdim[2];
  if(jx >= grt_ptr->grprdim[2])
    jx = grt_ptr->grprdim[2] - 1;

  yp = (ypos/ypixsz) - floor(ypos/ypixsz);
  jy = yp*grt_ptr->grprdim[3];
  if(jy >= grt_ptr->grprdim[3])
    jy = grt_ptr->grprdim[3] - 1;

  return(get_random_grade(jen,jx,jy));
}
/****************************************************************************/

void mkdbl(unsigned char a[8], double *d)
{
  int i;
  union myuld {
    unsigned long long ull;
    double ud;
  } myu;

  myu.ull = 0;
  for(i=0;i<8;++i) {
    myu.ull = (myu.ull << 8) + a[i];
  }
  *d = myu.ud;
  return;
}

/****************************************************************************/

void mkflt(unsigned char a[4], float *d)
{
  int i;
  union myuld {
    unsigned long ull;
    float ud;
  } myu;

  myu.ull = 0;
  for(i=0;i<4;++i) {
    myu.ull = (myu.ull << 8) + a[i];
  }
  *d = myu.ud;
  return;
}

/****************************************************************************/

void mksll(unsigned char a[8], signed long long *d)
{
  int i;
  union myuld {
    unsigned long long ull;
    double ud;
  } myu;

  myu.ull = 0;
  for(i=0;i<8;++i) {
    myu.ull = (myu.ull << 8) + a[i];
  }
  *d = myu.ud;
  return;
}

/****************************************************************************/

void mkint(unsigned char a[4], int *d)
{
  int i;
  union myuld {
    unsigned long ull;
    int ud;
  } myu;

  myu.ull = 0;
  for(i=0;i<4;++i) {
    myu.ull = (myu.ull << 8) + a[i];
  }
  *d = myu.ud;
  return;
}

/****************************************************************************/

void mkshort(unsigned char a[2], short *d)
{
  int i;
  union myuld {
    unsigned long ull;
    int ud;
  } myu;

  myu.ull = 0;
  for(i=0;i<2;++i) {
    myu.ull = (myu.ull << 8) + a[i];
  }
  *d = myu.ud;
  return;
}

/****************************************************************************/
/* nblks = number of 2880-byte (=36 x 80 bytes) blocks to read 
 *         (this better be two or more in order for the code to 
 *         determine a sensible value for naxis2)
 */

int read_header(FILE *fpin)
{
  char hdr[81], ach[2];
  int i, j, k, km1, iend, ng, nbyt;

  naxes = -1;
  for(i=0;i<GR5_MAX_N_AX;++i){
    nax[i] = -1;
  }
  iend = 0;
  iext = 0;
  ihdr = -1;

  nbyt = 0;
  for(i=0;i<GR5_MAX_N_HDR;++i){
    if(extra_stderr_diag == 1)
	fprintf(stderr,"starting header block %d\n",i);fflush(stderr);
    if(print_grfits_diag==1) {
      fprintf(fpgrdiag,"starting header block %d\n",i); }
    for(j=0;j<36;++j) {
      /* if(print_grfits_diag==1) {
	 fprintf(fpgrdiag,"j= %d\n",j); } */
      ng = fread(hdr,1,80,fpin);
      if(ng < 80) {
	iend = 1;
	break;
      }
      nbyt += ng;
      hdr[80] = 0;      /* null termination of string */
      if(print_grfits_diag==1) {
	fprintf(fpgrdiag,"%d %d %s\n",i,j,hdr); }
      if (!strncmp("SIMPLE",hdr,6)) {
	// This is the primary header.
	ihdr = 0;
	sscanf(hdr+6," %*s %s",hdrstr);
	if(print_grfits_diag==1) {
	  fprintf(fpgrdiag,"SIMPLE = %s\n",hdrstr); }
      }
      else if(!strncmp("XTENSION",hdr,8)) {
	// This is an extension header.
	sscanf(hdr+8," %*s %s",hdrstr);
        if(!strncmp("IMAGE",hdrstr+1,5)) {
	  // This is an IMAGE extension header.
	  ihdr = 1;
	}
        else if(!strncmp("TABLE",hdrstr+1,5)) {
	  // This is a TABLE (ASCII table) extension header.
	  ihdr = 2;
	}
        else if(!strncmp("BINTABLE",hdrstr+1,8)) {
	  // This is an BINTABLE extension header.
	  ihdr = 3;
	  ncol = 0;  // initialization
	}
	if(print_grfits_diag==1) {
	  fprintf(fpgrdiag,"XTENSION = %s %d\n",hdrstr,ihdr); }
      }
      if (!strncmp("BITPIX",hdr,6)) {
	sscanf(hdr+6," %*s %d",&bitp);
	if(print_grfits_diag==1) {
	  fprintf(fpgrdiag,"BITPIX = %d\n",bitp); }
	if(bitp > 0) {
	  intflt[0] = 1;  // integer type
	  numbyt[0] = bitp/8;
	}
	else {
	  intflt[0] = 2;  // floating type
	  numbyt[0] = -bitp/8;
	}
      }
      if (naxes >= 0) {
	if (!strncmp("NAXIS",hdr,5)) {
	  sscanf(&hdr[5],"%d",&k);
	  sscanf(hdr+6," %*s %d",&nax[k-1]);
	  if(print_grfits_diag==1) {
	    fprintf(fpgrdiag,"NAXIS%d = %d\n",k,nax[k-1]); }
	}
      }
      if (naxes == -1) {
	if (!strncmp("NAXIS ",hdr,6)) {
	  sscanf(hdr,"%*s %*s %d",&naxes);
	  if(print_grfits_diag==1) {
	    fprintf(fpgrdiag,"NAXIS = %d\n",naxes); }
	}
      }
      if (ihdr == 3) {
	if (!strncmp("TFIELDS",hdr,7)) {
	  sscanf(hdr+7," %*s %d",&nflds);
	  if(print_grfits_diag==1) {
	    fprintf(fpgrdiag,"TFIELDS = %d\n",nflds); }
	}
	else if (!strncmp("TTYPE",hdr,5)) {
	  ++ncol;
	  sscanf(hdr+5,"%d",&k);
	  /* if(print_grfits_diag==1) {
	     fprintf(fpgrdiag,"k = %d\n",k); } */
	  sscanf(hdr+7," %*s %s",&ttype[k-1][0]);
	  if(print_grfits_diag==1) {
	    fprintf(fpgrdiag,"TTYPE%d = %s\n",k,&ttype[k-1][0]); }
	}
	else if (!strncmp("TFORM",hdr,5)) {
	  sscanf(hdr+5,"%d",&k);
	  /* if(print_grfits_diag==1) {
	     fprintf(fpgrdiag,"k = %d\n",k); } */
	  km1 = k - 1;
	  sscanf(hdr+7," %*s %s",&tform[km1][0]);
	  if(print_grfits_diag==1) {
	    fprintf(fpgrdiag,"TFORM%d = %s\n",k,&tform[km1][0]); }
	  if(tform[km1][1] == 'D') {
	    intflt[km1] = 2;
	    numbyt[km1] = 8;
	  }
	  else if(tform[km1][1] == 'E') {
	    intflt[km1] = 2;
	    numbyt[km1] = 4;
	  }
	  else if(tform[km1][1] == 'I') {
	    intflt[km1] = 1;
	    numbyt[km1] = 2;
	  }
	  else if(tform[km1][1] == 'J') {
	    intflt[km1] = 1;
	    numbyt[km1] = 4;
	  }
	  else if(tform[km1][1] == 'K') {
	    intflt[km1] = 1;
	    numbyt[km1] = 8;
	  }
	  if(print_grfits_diag==1) {
	    fprintf(fpgrdiag,"k,km1,intflt[km1],numbyt[km1] = %d %d %d %d\n",
		    k,km1,intflt[km1],numbyt[km1]); }
	}
      }
      if (!strncmp("EXTEND",hdr,6)) {
	sscanf(hdr+6," %*s %1s",ach);
	ach[1] = 0;
	if(ach[0] == 'T')
	  iext = 1;
	if(print_grfits_diag==1) {
	  fprintf(fpgrdiag,"EXTEND = %1s; iext = %d\n",ach,iext); }
      }
      if (!strncmp("END",hdr,3)) {
	iend = 1;
	if(extra_stderr_diag == 1)
	  fprintf(stderr,"END found.\n");
	if(print_grfits_diag==1) {
	  fprintf(fpgrdiag,"END found.\n"); }
      }
    }  // end of j loop
    if(iend == 1) break;
  }
  return(nbyt);
}

/****************************************************************************/
/* return value = 0 => error
 *              = 1 => no error
 */

int convert_data(int ifl, int numb, unsigned char inbuf[8], int ipr, FILE *fppr)
{

  if(ifl == 1 ) { // integer types
    if(numb == 2) {
	  mkshort(inbuf,&pshrt);
	  if(ipr)
	    fprintf(fppr,"%d\n",pshrt);
    }
    if(numb == 4) {
	  mkint(inbuf,&pint);
	  if(ipr)
	    fprintf(fppr,"%d\n",pint);
    }
    else if(numb == 8) {
	  mksll(inbuf,&psll);
	  if(ipr)
	    fprintf(fppr,"%lld\n",psll);
    }
    else {
      return(0);
    }
  }
  else if(ifl == 2) {  // floating types
    if(numb == 4) {
	  mkflt(inbuf,&pflt);
	  if(ipr)
	    fprintf(fppr,"%e\n",pflt);
    }
    else if(numb == 8) {
	  mkdbl(inbuf,&pdbl);
	  if(ipr)
	    fprintf(fppr,"%e\n",pdbl);
    }
    else {
      return(0);
    }
  }
  else {
    return(0);
  }
  return(1);
}

/****************************************************************************/

int read_pad(int blkbyts, int bytsrd, FILE *infile)
{
  unsigned char inbuf[8];
  int i, ng, npart, nleft, nrd;

  npart = bytsrd % blkbyts;
  nrd = 0;
  if(npart > 0) {
    nleft = (blkbyts - npart);
    if(print_grfits_diag==1) {
      fprintf(fpgrdiag,"blkbyts,bytsrd,npart,nleft = %d %d %d %d\n",
	      blkbyts,bytsrd,npart,nleft); }
    for(i=0;i<nleft;++i) {
      ng = fread(inbuf,1,1,infile);
      nrd += ng;
    }
  }
  return(nrd);
}

/****************************************************************************/
// number of image entries = nim
// number of bytes per entry = numbyt[0]

int read_image_data(FILE *infile)
{
  unsigned char inbuf[8];
  int i, ng, nim, nbytim, nrd, npad, ib, ic;

  nrd = 0;
  npad = 0;

  if (naxes > 0) {
    nim = 1;
    for(i=0;i<naxes;++i) {
      if(extra_stderr_diag == 1)
	fprintf(stderr,"i,nax = %d %d\n",i,nax[i]);
      nim *= nax[i];
    }
    nbytim = nim*numbyt[0];
    if(print_grfits_diag==1) {
      fprintf(fpgrdiag,"nim,numbyt[0],nbytim = %d %d %d\n\n",nim,numbyt[0],nbytim); }
    if(nbytim > 0) {
      // Change into a 'for' statement for each axis?  not now.
      for(i=0;i<nim;++i) {
	ng = fread(inbuf,1,numbyt[0],infile);
	if (ng != numbyt[0]) {
	  fprintf(stderr,"ERROR - ng = %d for i = %d; numbyt[0] = %d\n",ng,i,numbyt[0]);
	}
	nrd += ng;
	// How does the destination variable which could have any of a set of types get specified?
	// one global variable of each type?
	if( (i<5) && (print_grfits_diag==1) )  {
	    convert_data(intflt[0],numbyt[0],inbuf,1,fpgrdiag);
	}
	else
	  convert_data(intflt[0],numbyt[0],inbuf,0,stderr);
	if(irec==0) {
	  ib = get_image_array_bin(irec,i);
	  ic = put_gr_prob(pdbl);
	}
	if(irec==2) {
	  if(bpgrt==16)
	    grt_ptr->gr_table[i] = pshrt;
	  else if(bpgrt==32)
	    grt_ptr->gr_table[i] = pint;
	  else if(bpgrt==64)
	    grt_ptr->gr_table[i] = psll;
	}
      }
      // read block pad
      npad = 0;
      npad = read_pad(GR5_BLKSIZ,nbytim,infile);
    }
  }
  return(nrd + npad);
}

/****************************************************************************/

int read_bintable_data(FILE *infile)
{
  unsigned char inbuf[8];
  int i, j, sng, nim, nbytim, nrd, npad, ng;

  if(extra_stderr_diag == 1)
    fprintf(stderr,"ncol = %d\n",ncol);
  nrd = 0;
  if (naxes > 0) {  // error if naxes != 2
    for(i=0;i<nax[1];++i) {   // nax[1] = number of rows in the binary table
      if(extra_stderr_diag == 1)
	fprintf(stderr,"i,nax[1] = %d %d\n",i,nax[1]);
      for(j=0;j<ncol;++j) {   // ncol = number of columns in a row
	if( (print_grfits_diag==1) && (i==0) ) {
	  fprintf(fpgrdiag,"j,numbyt[j] = %d %d\n",j,numbyt[j]); }
	ng = fread(inbuf,1,numbyt[j],infile);
	nrd += ng;
	if (ng != numbyt[j]) {
	  fprintf(stderr,"ERROR - ng = %d for i = %d\n",ng,i);
	}
	// How does the destination variable which could have any of a set of types get specified?
	// one global variable of each type?
	if( (i<3) && (print_grfits_diag==1) ) {
	    convert_data(intflt[j],numbyt[j],inbuf,1,fpgrdiag);
	}
	else
	  convert_data(intflt[j],numbyt[j],inbuf,0,stderr);
	if(irec==1) {
	  grt_ptr->ebin[i][iecol[j]] = pdbl;
	}
      }
    }
    // read block pad
    npad = 0;
    npad = read_pad(GR5_BLKSIZ,nrd,infile);
  }
  return(nrd + npad);
}

/****************************************************************************/
/* Grade probabilities should be in the primary HDU (an IMAGE),
 * energy bin definitions in the first extension (a BINTABLE),
 * and grades associated with each grade bin inn the second extension (an IMAGE).
 */

int read_grades_fits_file(FILE *gradesfp)
{
  double pix;
  int i, nb, nba, nrd, npart, nleft, ng, ntot, igrpr, icgr;
  signed long long sll;

  alloc_gradepr();

  irec = 0;
  ntot = 0;
  if(extra_stderr_diag == 1) fprintf(stderr,"\nBeginning read_grades_fits_file()\n");
  if(print_grfits_diag==1) fprintf(fpgrdiag,"\nBeginning read_grades_fits_file()\n");
  while(1) {
    nrd = read_header(gradesfp);
    if(extra_stderr_diag == 1) fprintf(stderr,"nrd = %d\n",nrd);
    if( (nrd <= 0) || ((nrd % 2880) != 0) ) {
      return(-1);
    }
    if(extra_stderr_diag == 1) fprintf(stderr,"\nnaxes = %d\n",naxes);
    ntot += nrd;
    if(extra_stderr_diag == 1) fprintf(stderr,"ntot = %d\n",ntot);
    if(irec==0) {
      for(i=0;i<4;++i) {
	grt_ptr->grprdim[i] = nax[i];
	if(i==0)
	  grt_ptr->grprmod[i] = 1;
	else
	  grt_ptr->grprmod[i] = nax[i-1]*grt_ptr->grprmod[i-1];
	if(extra_stderr_diag == 1) 
	  fprintf(stderr,"i,grprdim[i],grprmod[i] = %d %d %d\n",i,grt_ptr->grprdim[i],grt_ptr->grprmod[i]);
      }
      if((icgr = check_grprdim()) < 0) {
	fprintf(stderr,"ERROR - grade table too big - icgr = %d; exiting ...\n",icgr);
	for(i=0;i<4;++i) {
	  if(extra_stderr_diag == 1) fprintf(stderr,"i,dim[i] = %d  %d\n",i,grt_ptr->grprdim[i]);
	}
	exit(-1);
      }
    }
    else if(irec==1) {
      grt_ptr->nebins = nax[1];
    }
    else if(irec==2) {
      grt_ptr->ngrades = nax[0];
      bpgrt = bitp;
    }

    // Read IMAGE data
    if( (ihdr == 0) || (ihdr == 1) ) {
      nrd = read_image_data(gradesfp);
      if(extra_stderr_diag == 1) fprintf(stderr,"irec, nrd = %d %d\n",irec,nrd);
      ntot += nrd;
      if(extra_stderr_diag == 1) fprintf(stderr,"ntot = %d\n",ntot);
    }

    // Read BINTABLE data
    if(ihdr == 3) {
      if(irec==1) {
	set_e_columns();
      }
      nrd = read_bintable_data(gradesfp);
      if(extra_stderr_diag == 1) fprintf(stderr,"irec, nrd = %d %d\n",irec,nrd);
      ntot += nrd;
      if(extra_stderr_diag == 1) fprintf(stderr,"ntot = %d\n",ntot);
    }
    ++irec;
  } // while loop 

  return(1);
}

/****************************************************************************/

int acis_grades_init(FILE *fpin)
{
  int i, ird;

  ird = read_grades_fits_file(fpin);
  if(print_grfits_diag==1) {
    fprintf(fpgrdiag,"grt_ptr->nebins, grt_ptr->ngrades = %d %d\n",
	    grt_ptr->nebins,grt_ptr->ngrades);
    print_en_bins(fpgrdiag);
    print_grade_list(fpgrdiag);
  }

  for(i=0;i<4;++i) {
    grt_ptr->ilim[i][0] = 0;
    grt_ptr->ilim[i][1] = grt_ptr->grprdim[i];
  }
  return(ird);
}

/****************************************************************************/

void grades_fits_file_path(char *filename)
{
  char *dir, *env;

  if (dir == NULL) {
    dir = getenv("MARX_DATA_DIR");
  }

  if (*dir == '$')
    {
      env = getenv (dir + 1);
      if (env == NULL)
	  {
	    marx_error ("Unable to set Data Search Path.\nEnvironment variable %s does not exist.",
			dir + 1);
	  }
      dir = env;
     }

   if (*dir == 0)
     {
	marx_error ("Data Path is empty!");
     }

  sprintf(grades_input_file,"%s/%s",dir,filename);
}

/****************************************************************************/

void read_acis_grades_files(int ireturn[2])
{
  char filein[512];
  int irdfs, irdbs, ipval;
  FILE *grfits_fpin;

   if(j_read_fsbs_fits[0] == 0) {
     marx_message("Opening file grade_probs_fs.fits for reading\n");
     fprintf(stderr,"Opening file grade_probs_fs.fits for reading\n");
     // read fs grades FITS file
     // open file (file name?)
     grades_fits_file_path("grade_probs_fs.fits");
     fprintf(stdout,"grades_input_file = %s\n",grades_input_file);
     grfits_fpin = fopen(grades_input_file,"r");
     // fprintf(stdout,"fpin = %d\n",(long) grfits_fpin);
     ipval = init_gr_p_tab_ptr(1);
     // CHECK FOR ERROR
     irdfs =  acis_grades_init(grfits_fpin);
     fclose(grfits_fpin);
     // CHECK FOR ERROR
     j_read_fsbs_fits[0] = 1;
     check_gr_probs(16,GR5_MAX_DEL,0,stderr);
     check_gr_en_range();
   }
   if(j_read_fsbs_fits[1] == 0) {
     fprintf(stdout,"Opening file grade_probs_bs.fits for reading\n");
     fprintf(stderr,"Opening file grade_probs_bs.fits for reading\n");
     // read bs grades FITS file
     // open file (file name?)
     grades_fits_file_path("grade_probs_bs.fits");
     fprintf(stdout,"grades_input_file = %s\n",grades_input_file);
     grfits_fpin = fopen(grades_input_file,"r");
     ipval = init_gr_p_tab_ptr(2);
     // CHECK FOR ERROR
     irdbs =  acis_grades_init(grfits_fpin);
     fclose(grfits_fpin);
     // CHECK FOR ERROR
     j_read_fsbs_fits[1] = 1;
     check_gr_probs(16,GR5_MAX_DEL,0,stderr);
     check_gr_en_range();
   }

   ireturn[0] = irdfs;
   ireturn[1] = irdbs;
}

/****************************************************************************/
