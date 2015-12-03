/* This file is part of MARX

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
# include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"

/*
 * The general form of the absorbtion is:
 *  acisabs(E,t,x,y) = \sum_j W_j acisabs_j(E,t,x,y)
 * Here j denotes a component and the above sum is over components.
 * For a 2 component model, there are two terms in the above sum.
 *
 * acisabs_j(E,t,x,y) represents the absorbtion for the jth
 * component.  It is represented by:
 *
 *  acisabs_j(E,t,x,y) = exp(-kappa_j*\sum_i \mu_ij(E)z_ij(x,y,t))
 *
 * Here, kappa_j, \mu_ij(E) and z_ij(x,y,t) depend upon the
 * component, which is why they each have a j index.  The sum over i is
 * over layers.
 *
 * It is also assumed that z_ij(x,y,t) is parametrized as
 *   z_ij(x,y,t) = \tau_ij0(t) + \tau_ij1(t) f_ij(x,y)
 */

/* For the AV/HLM hybrid model defined in 2009, only a single multilayer
 * componenent is needed.  Also, the spatial dependence f_ij(x,y) is the
 * same for all layers.  As a result, we have:
 *
 *    acisabs(E,t,x,y) = exp(-\sum_i\mu_i(E)(\tau_i0(t)+\tau_i1(t)f(x,y)))
 *
 * Here, the spatial dependence f(x,y) will be computed on the fly, and the
 * tables for (E,\mu_i) and (t,\tau_i0,1) will be read in.
 *
 * The time-dependent functions are slowly varying and can be assumed to be
 * a constant during an observation.  So a set of fixed constants {\tau}
 * may be assumed.  The only tables that will need interpolated on an
 * event-by-event basis are for the \mu value.  Hence, for fixed \tau:
 *
 *   -log(acisabs(E,x,y))
 *      = \sum_i \tau_i0\mu_i(E)   +   f(x,y)\sum_i\tau_i1\mu_i(E)
 */

#define MAX_LAYERS 5
typedef struct Single_Component_Contam_Type
{
   unsigned int num_layers;
   double tau_0s[MAX_LAYERS];
   double tau_1s[MAX_LAYERS];
   float *mus[MAX_LAYERS];
   float *energies[MAX_LAYERS];
   unsigned int num_mus[MAX_LAYERS];

   /* if fxy == NULL, use fxy_vals as a lookup table. */
   double (*fxy)(double, double);
   float *fxy_vals[MAX_LAYERS];
   unsigned int blocking_factor;
}
Single_Component_Contam_Type;

static double compute_contamination (Single_Component_Contam_Type *c,
				     double en, double cx, double cy)
{
   double fxy;
   unsigned int i;
   double v;

   if (c->fxy != NULL)
     {
        /* This lines assumes the fxy is the same for every i */
	fxy = (*c->fxy)(cx, cy);

	v = 0.0;
	for (i = 0; i < c->num_layers; i++)
	  {
	     double mu;

	     mu = JDMinterpolate_f (en, c->energies[i], c->mus[i], c->num_mus[i]);
	     v += mu*(c->tau_0s[i] + c->tau_1s[i]*fxy);
	  }
     }
   else
     {
	unsigned int ofs;

	if ((cx < 0) || (cx >= 1024) || (cy < 0) || (cy >= 1024))
	  return 0.0;
	cx /= c->blocking_factor;
	cy /= c->blocking_factor;
	ofs = (1024/c->blocking_factor) * (unsigned int) cy
	  + (unsigned int) cx;
	v = 0.0;

	for (i = 0; i < c->num_layers; i++)
	  {
	     double mu;

	     mu = JDMinterpolate_f (en, c->energies[i], c->mus[i], c->num_mus[i]);
	     fxy = c->fxy_vals[i][ofs];
	     v += mu*(c->tau_0s[i] + c->tau_1s[i]*fxy);
	  }
     }
   return exp (-v);
}

static Single_Component_Contam_Type Single_Component_Contam_Table[10];

static double ArcMin_Per_Pixel = 8.0/1024.0;
static double fxy_acis_i (double x, double y, double x_0, double y_0)
{
   double dx = x - x_0;
   double dy = y - y_0;
   double r = ArcMin_Per_Pixel * sqrt(dx*dx + dy*dy);
#if 0
   /* This was for the N0005 contam model */
   /* for gamma=2 */
   double a = 8.07*8.07, b=3.24*3.24;
   return (r*r - b)/(a-b);
#else
   r /= 8.07;
   return 1.29*r*r;      /* N0006 model */
#endif
}

static double fxy_acis0 (double x, double y)
{
   return fxy_acis_i (x, y, 1024.5, 1024.5);
}

static double fxy_acis1 (double x, double y)
{
   return fxy_acis_i (x, y, 0.5, 1024.5);
}
static double fxy_acis2 (double x, double y)
{
   return fxy_acis_i (x, y, 0.5, 1024.5);
}
static double fxy_acis3 (double x, double y)
{
   return fxy_acis_i (x, y, 1024.5, 1024.5);
}

static double fxy_acis456789 (double x, double y)
{
   double y_0 = 512.0;
   double alpha_1 = 5.5;
   double alpha_2 = 4.5;

   (void) x;

   if (y <= 512.0)
     return pow(fabs((y-y_0)/(64.0-y_0)), alpha_1);
   else
     return pow(fabs((y-y_0)/(964.0-y_0)), alpha_2);
}

#define CONTAM_FUN(name, i) \
   static double name (_Marx_Acis_Chip_Type *ccd, \
		       double en, double cx, double cy) \
   { \
      (void) ccd; \
      return compute_contamination (Single_Component_Contam_Table+(i), \
				    en, cx, cy); \
   }

CONTAM_FUN(contam_0, 0)
CONTAM_FUN(contam_1, 1)
CONTAM_FUN(contam_2, 2)
CONTAM_FUN(contam_3, 3)
CONTAM_FUN(contam_4, 4)
CONTAM_FUN(contam_5, 5)
CONTAM_FUN(contam_6, 6)
CONTAM_FUN(contam_7, 7)
CONTAM_FUN(contam_8, 8)
CONTAM_FUN(contam_9, 9)

#define NUM_CONTAM_COLUMNS	8
#define NUM_CONTAM_OPT_COLUMNS	1
static char *Contam_File_Columns[NUM_CONTAM_COLUMNS+NUM_CONTAM_OPT_COLUMNS] =
{
#define COMPONENT_COLUMN	0
   "i:component",
#define N_ENERGY_COLUMN		1
   "i:n_energy",
#define ENERGY_COLUMN		2
   "f:energy",
#define MU_COLUMN		3
   "f:mu",
#define N_TIME_COLUMN		4
   "i:n_time",
#define TIME_COLUMN		5
   "d:time",
#define TAU0_COLUMN		6
   "f:tau0",
#define TAU1_COLUMN		7
   "f:tau1",

   /* optional columns go here */
#define FXY_COLUMN		8
   "f:fxy",
};

static void compute_tau_values (Single_Component_Contam_Type *contam,
				unsigned int layer, unsigned int ntimes,
				double *times, float *tau0s, float *tau1s)
{
   unsigned int n0, n1;
   double t = _Marx_TStart_MJDsecs;
   double t0, t1;

   n0 = JDMbinary_search_d (t, times, ntimes);
   n1 = n0+1;

   if (n1 == ntimes)
     {
	if (n0 == 0)
	  {
	     contam->tau_0s[layer] = tau0s[0];
	     contam->tau_1s[layer] = tau1s[0];
	     return;
	  }
	n1 = n0-1;
     }

   t0 = times[n0];
   t1 = times[n1];
   t = (t-t0)/(t1-t0);

   contam->tau_0s[layer] = tau0s[n0] + (tau0s[n1] - tau0s[n0]) * t;
   contam->tau_1s[layer] = tau1s[n0] + (tau1s[n1] - tau1s[n0]) * t;
}

static void free_contam_table_entry (int ccd)
{
   unsigned int layer;
   Single_Component_Contam_Type *c = Single_Component_Contam_Table+ccd;

   for (layer = 0; layer < c->num_layers; layer++)
     {
	if (c->mus[layer] != NULL)
	  marx_free ((char *) c->mus[layer]);
	if (c->energies[layer] != NULL)
	  marx_free ((char *) c->energies[layer]);
	if (c->fxy_vals[layer] != NULL)
	  marx_free ((char *) c->fxy_vals[layer]);
     }
   memset ((char *)c, 0, sizeof(Single_Component_Contam_Type));
}

static int find_binary_table_callback (void *cd, JDFits_Type *ft)
{
   int ccd;
   char extname[32];
   char *val;
   int ival;

   ccd = *(int *) cd;

   if (-1 == jdfits_read_keyword_string (ft, "EXTNAME", &val))
     return 0;
   /* Files distribued with marx use ACIS$ccd_CONTAM as extnam.
    * The Chandra CALDB uses AXAF_CONTAM, with CCD_ID keyword in header.
    */
   sprintf (extname, "ACIS%d_CONTAM", ccd);
   if (0 == jdfits_strcasecmp (val, extname))
     {
	jdfits_free (val);
	return 1;
     }

   if (0 != jdfits_strcasecmp (val, "AXAF_CONTAM"))
     {
	jdfits_free (val);
	return 0;
     }
   jdfits_free (val);

   if (-1 == jdfits_read_keyword_int (ft, "CCD_ID", &ival))
     return 0;

   if (ival == ccd)
     return 1;

   return 0;
}


static int read_contam_file_for_ccd (char *file, int ccd)
{
   JDFits_Type *f;
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   unsigned int num_layers, layer, num_columns;
   Single_Component_Contam_Type *contam;
   int has_fxy;

   free_contam_table_entry (ccd);

   contam = Single_Component_Contam_Table+ccd;

   marx_message ("Reading ACIS%d Contamination File\n", ccd);
   marx_message ("\t%s\n", file);

   if (NULL == (f = jdfits_find_binary_table (file, find_binary_table_callback, (void *)&ccd)))
     {
	marx_error ("Unable to open acis contam file %s for CCD %d", file, ccd);
	return -1;
     }

   has_fxy = (1 == jdfits_bintable_column_exists (f, "fxy"));
   num_columns = NUM_CONTAM_COLUMNS;
   if (has_fxy)
     {
	num_columns++;
	contam->fxy = NULL;
     }
   else switch (ccd)
     {
      case 0: contam->fxy = fxy_acis0; break;
      case 1: contam->fxy = fxy_acis1; break;
      case 2: contam->fxy = fxy_acis2; break;
      case 3: contam->fxy = fxy_acis3; break;

      default:
	contam->fxy = fxy_acis456789;
     }

   r = jdfits_bintable_aopen_rows (f, num_columns, Contam_File_Columns);
   if (r == NULL)
     {
	marx_error ("Error processing ACIS contam table %s for CCD=%d", file, ccd);
	jdfits_close_file (f);
	return -1;
     }
   num_layers = r->num_rows;
   if (num_layers > MAX_LAYERS)
     {
	marx_error ("The number of layers in %s for CCD=%d is larger than supported",
		    file, ccd);
	jdfits_close_file (f);
	return -1;
     }
   contam->num_layers = num_layers;

   c = r->col_data;

   for (layer = 0; layer < num_layers; layer++)
     {
	unsigned int n_time, n_energy;
	double *times;
	float *tau0s, *tau1s, *energies, *mus;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Unexpected end of binary table %s", file);
	     goto return_error_bad_row;
	  }
	c = r->col_data;

	/* This implementation assumes a single component */
	if (0 != c[COMPONENT_COLUMN].data.i[0])
	  {
	     marx_error ("Expecting COMPONENT=0 in %s for CCD=%d", file, ccd);
	     goto return_error_bad_row;
	  }
	n_time = (unsigned int) c[N_TIME_COLUMN].data.i[0];
	n_energy = (unsigned int) c[N_ENERGY_COLUMN].data.i[0];

	if ((n_time <= 1) || (n_energy <= 1))
	  {
	     marx_error ("Expecting n_energy or n_time to be greater than 1");
	     goto return_error_bad_row;
	  }
	if ((n_time > c[TIME_COLUMN].repeat)
	    || (n_time > c[TAU0_COLUMN].repeat)
	    || (n_time > c[TAU1_COLUMN].repeat))
	  {
	     marx_error ("n_time value is larger than data arrays");
	     goto return_error_bad_row;
	  }

	if ((n_energy > c[ENERGY_COLUMN].repeat)
	    || (n_energy > c[MU_COLUMN].repeat))
	  {
	     marx_error ("n_energy value is larger than data arrays");
	     goto return_error_bad_row;
	  }

	times = c[TIME_COLUMN].data.d;
	tau0s = c[TAU0_COLUMN].data.f;
	tau1s = c[TAU1_COLUMN].data.f;
	energies = c[ENERGY_COLUMN].data.f;
	mus = c[MU_COLUMN].data.f;

	if (-1 == _marx_check_monotonicity_d (times, n_time))
	  {
	     marx_error ("Expecting the time grid to be monotonically increasing");
	     goto return_error_bad_row;
	  }

	if (-1 == _marx_check_monotonicity_f (energies, n_energy))
	  {
	     marx_error ("Expecting the ENERGY grid to be monotonically increasing");
	     goto return_error_bad_row;
	  }

	compute_tau_values (contam, layer, n_time, times, tau0s, tau1s);

	if ((NULL == (contam->energies[layer] = (float *)marx_malloc (sizeof(float)*n_energy)))
	    || (NULL == (contam->mus[layer] = (float *)marx_malloc (sizeof(float)*n_energy))))
	  goto return_error_bad_row;

	memcpy ((char *)contam->energies[layer], energies, n_energy*sizeof(float));
	memcpy ((char *)contam->mus[layer], mus, n_energy*sizeof(float));
	contam->num_mus[layer] = n_energy;

	if (has_fxy)
	  {
	     JDFits_Keyword_Type *k;
	     unsigned int repeat = 1024*1024;
	     contam->blocking_factor = 1;

	     if (NULL != (k = jdfits_find_keyword (f, "fxyblk")))
	       {
		  int b;
		  if (-1 == jdfits_extract_integer (k, &b))
		    goto return_error_bad_row;
		  if ((b <= 0) || (b * (1024/b) != 1024))
		    {
		       marx_error ("Invalid fxyblk value");
		       goto return_error_bad_row;
		    }
		  contam->blocking_factor = b;
		  repeat = (1024/b) * (1024/b);
	       }

	     if (repeat != c[FXY_COLUMN].repeat)
	       {
		  marx_error ("The FXY column is expected to have %u * %u values",
			      1024/contam->blocking_factor, 1024/contam->blocking_factor);
		  goto return_error_bad_row;
	       }
	     if (NULL == (contam->fxy_vals[layer] = (float *)marx_malloc(sizeof(float)*repeat)))
	       goto return_error_bad_row;
	     memcpy ((char *)contam->fxy_vals[layer],
		     (char *)c[FXY_COLUMN].data.f, repeat*sizeof(float));
	  }
     }

   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   return 0;

return_error_bad_row:

   marx_error ("Error encountered reading row %d of %s for CCD=%d",
	       layer+1, file, ccd);
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   return -1;
}

int _marx_acis_contam_init (Param_File_Type *p, _Marx_Acis_Chip_Type *ccd)
{
   char *file;

   (void) p;

   if (NULL == (file = _marx_caldb_get_file ("ACISCONTAM")))
     return -1;

   if (-1 == read_contam_file_for_ccd (file, ccd->ccd_id))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   switch (ccd->ccd_id)
     {
      case 0:
	ccd->contam_fun = contam_0;
	break;
      case 1:
	ccd->contam_fun = contam_1;
	break;
      case 2:
	ccd->contam_fun = contam_2;
	break;
      case 3:
	ccd->contam_fun = contam_3;
	break;
      case 4:
	ccd->contam_fun = contam_4;
	break;
      case 5:
	ccd->contam_fun = contam_5;
	break;
      case 6:
	ccd->contam_fun = contam_6;
	break;
      case 7:
	ccd->contam_fun = contam_7;
	break;
      case 8:
	ccd->contam_fun = contam_8;
	break;
      case 9:
	ccd->contam_fun = contam_9;
	break;

      default:
	marx_error ("_marx_acis_contam_init: unsupported ccdid: %d", ccd->ccd_id);
	return -1;
     }
   return 0;
}

