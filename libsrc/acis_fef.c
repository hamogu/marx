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

#if MARX_HAS_ACIS_FEF

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include <pfile.h>
#include <jdfits.h>

#include "marx.h"
#include "_marx.h"

/*
 * A single FEF extension contains a function written here as
 * R(h,E,i,j,ccdid) where h is the pha, E is the energy, (i,j)
 * specifies a pixel on the CCD, and ccdid specifies the CCD.  This
 * function is encoded as sum of gaussians:
 *
 *    G(h,a(E,i,j,ccdid))
 *
 * where a represents the gaussian parameters (width, amplitude, and
 * center), and depends upon the region (ccdid,i,j) and energy E. The
 * FEF itself contains the functional form (sum of gaussians) in the
 * header and encodes the parameters a(E,i,j,ccdid) in the binary
 * table.  The most obvious way of storing this information is via a
 * multi-dimensional array (possibly spread over multiple extensions
 * for each CCDID) --- after all, that is what 'a' is.  However the
 * FEF stores it in a very poorly designed manner such that each row
 * pertains to a specific (E,i,j,ccdid).  The most obvious failing of
 * this format is that there is no way of knowing how many energies
 * are used for the energy grid, nor are there any constraints on the
 * ordering of the rows.  For example, are they sorted on energy?  There
 * have been cases of programs that broke when the order of the rows was
 * changed.
 *
 * Hence, this implementation makes 2 passes through the table, with the
 * first determining the number of energies per region and to make sure
 * that there are no duplicate energies.
 *
 * Implementation notes:
 *
 * The regions defined by the CHIPX/Y_MIN/MAX values are not
 * uniform in size.  However, the values are _supposed_ to be multiples
 * of 32 and this is what will be assumed in this code.  Hence, at most
 * there will be 32x32=1024 regions per CCD.
 *
 * Within a region, the energy column is assumed to be in ascending order.
 */
#define MIN_REGION_SIZE	32
#define NUM_REGIONS (1024/MIN_REGION_SIZE)
#define FEF_EXTNAME "FUNCTION"

#define FIRST_CHAN 1

typedef struct
{
   float amp;
   float center;
   float sigma;
   float cum_area;
   int use_tail_dist;
}
Gauss_Parm_Type;

typedef struct
{
   unsigned int num_gaussians;
   unsigned int num_energies;
   float *energies;		       /* [num_energies] */
   float *channels;		       /* [num_energies] */
   Gauss_Parm_Type *gaussians;	       /* [num_energies * num_gaussians] */
   unsigned int num_refs;	       /* reference count */
}
Fef_Type;

#define NUM_ACIS_CCDS	10
typedef struct
{
   Fef_Type *fef_map [NUM_REGIONS][NUM_REGIONS];
}
Fef_Map_Type;

static Fef_Map_Type *Fef_Maps[NUM_ACIS_CCDS];

static void free_fef_type (Fef_Type *f)
{
   if (f == NULL)
     return;

   if (f->num_refs > 1)
     {
	f->num_refs--;
	return;
     }
   if (f->energies != NULL)
     marx_free ((char *)f->energies);
   if (f->channels != NULL)
     marx_free ((char *)f->channels);
   if (f->gaussians != NULL)
     marx_free ((char *)f->gaussians);

   marx_free ((char *) f);
}

static void free_fef_map (Fef_Map_Type *map)
{
   unsigned int i, j;

   if (map == NULL)
     return;

   for (i = 0; i < NUM_REGIONS; i++)
     {
	for (j = 0; j < NUM_REGIONS; j++)
	  {
	     free_fef_type (map->fef_map[i][j]);
	     map->fef_map[i][j] = NULL;
	  }
     }
}

static void free_fef_maps (void)
{
   unsigned int i;

   for (i = 0; i < NUM_ACIS_CCDS; i++)
     {
	free_fef_map (Fef_Maps[i]);
	Fef_Maps[i] = NULL;
     }
}

static int allocate_fef_maps (unsigned int min_ccdid, unsigned max_ccdid)
{
   unsigned int i;

   if (min_ccdid >= NUM_ACIS_CCDS)
     {
	marx_error ("Internal Error: min/max_ccdid > NUM_ACIS_CCDS");
	return -1;
     }

   free_fef_maps ();

   for (i = min_ccdid; i <= max_ccdid; i++)
     {
	Fef_Map_Type *f;

	f = (Fef_Map_Type *) marx_malloc (sizeof (Fef_Map_Type));
	if (f == NULL)
	  {
	     free_fef_maps ();
	     return -1;
	  }
	memset ((char *) f, 0, sizeof (Fef_Map_Type));
	Fef_Maps [i] = f;
     }

   return 0;
}

static int check_region_values (int ccdid, int xmin, int ymin, int xmax, int ymax)
{
   if ((ccdid < 0) || (ccdid >= NUM_ACIS_CCDS))
     {
	marx_error ("Invalid CCD_ID value in gain file: %d", ccdid);
	return -1;
     }

   if ((xmin < 1) || (ymin < 1) || (xmax > 1024) || (ymax > 1024))
     {
	marx_error ("CHIPX/Y_LO/HI value out of range");
	return -1;
     }

   xmin--;			       /* 1->0 based */
   ymin--;

   if ((xmax % MIN_REGION_SIZE)
       || (ymax % MIN_REGION_SIZE)
       || (xmin % MIN_REGION_SIZE)
       || (ymin % MIN_REGION_SIZE))
     {
	marx_error ("Expecting ACIS gain region size to be a multiple of %d",
		    MIN_REGION_SIZE);
	return -1;
     }

   return 0;
}

static int map_fef_to_region (int ccdid, int xmin, int ymin, int xmax, int ymax,
			      Fef_Type *f)
{
   Fef_Map_Type *m;
   int i;

   if (-1 == check_region_values (ccdid, xmin, ymin, xmax, ymax))
     return -1;

   xmin--;			       /* 1->0 based */
   ymin--;

   m = Fef_Maps[ccdid];
   if (m == NULL)
     {
	marx_error ("Internal Error: Fef_Maps[%d] is NULL", ccdid);
	return -1;
     }

   xmin /= MIN_REGION_SIZE;
   xmax /= MIN_REGION_SIZE;
   ymin /= MIN_REGION_SIZE;
   ymax /= MIN_REGION_SIZE;

   for (i = xmin; i < xmax; i++)
     {
	int j;
	for (j = ymin; j < ymax; j++)
	  {
	     if (m->fef_map[i][j] != NULL)
	       {
		  marx_message ("***Warning: Fef Map defines an overlapping region for CCD_ID=%d at (%d, %d)\n",
				ccdid, 1+i*MIN_REGION_SIZE, 1+j*MIN_REGION_SIZE);
	       }
	     m->fef_map[i][j] = f;
	     f->num_refs += 1;
	  }
     }

   return 0;
}

static int check_fef_map_for_holes (int min_ccdid, int max_ccdid)
{
   int ccdid;

   for (ccdid = min_ccdid; ccdid <= max_ccdid; ccdid++)
     {
	Fef_Map_Type *f = Fef_Maps [ccdid];
	unsigned int i, j;

	if (f == NULL)
	  {
	     marx_error ("Internal Error: Fef_Maps[%d] is NULL", ccdid);
	     return -1;
	  }

	for (i = 0; i < NUM_REGIONS; i++)
	  {
	     for (j = 0; j < NUM_REGIONS; j++)
	       {
		  if (f->fef_map[i][j] != NULL)
		    continue;

		  marx_message ("***Warning: FEF for CCDID=%d contains holes\n",
				ccdid);
		  i = NUM_REGIONS;
		  break;
	       }
	  }
     }

   return 0;
}

static int check_fef_validity (Fef_Type *f)
{
   if (-1 == _marx_check_monotonicity_f (f->energies, f->num_energies))
     {
	marx_error ("Fef Map does not have monotonically increasing energies");
	return -1;
     }
   if (-1 == _marx_check_monotonicity_f (f->channels, f->num_energies))
     {
	marx_error ("Fef Map does not have monotonically increasing PHAs");
	return -1;
     }
   return 0;
}

static int check_repeat (JDFits_Col_Data_Type *c, unsigned int n)
{
   unsigned int i;
   for (i = 0; i < n; i++)
     {
	if (c[i].repeat != 1)
	  {
	     marx_error ("Expecting a repeat count of 1 for column %d", i+1);
	     return -1;
	  }
     }
   return 0;
}

#define MAX_FEF_COLUMNS 64
#define CCDID_COLUMN	0
#define CHIPX_LO_COLUMN	1
#define CHIPX_HI_COLUMN	2
#define CHIPY_LO_COLUMN	3
#define CHIPY_HI_COLUMN	4
#define REGNUM_COLUMN	5
#define ENERGY_COLUMN	6
#define CHANNEL_COLUMN	7
#define FWHM1_COLUMN	8
#define MAX_GAUSSIANS	((MAX_FEF_COLUMNS - FWHM1_COLUMN)/3)
static char *Fef_Columns [MAX_FEF_COLUMNS] =
{
   "i:CCD_ID", "i:CHIPX_LO", "i:CHIPX_HI", "i:CHIPY_LO", "i:CHIPY_HI",
     "i:REGNUM", "f:ENERGY", "f:CHANNEL"	       /* rest are gauss parms */
};

static JDFits_Type *open_fef_file (char *file)
{
   JDFits_Type *f;

   if (NULL == (f = jdfits_open_binary_table (file, FEF_EXTNAME)))
     marx_error ("Unable to open ACIS FEF file %s", file);

   return f;
}

static JDFits_Row_Type *open_rows (JDFits_Type *f, unsigned int num_columns)
{
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;

   if (NULL == (r = jdfits_bintable_aopen_rows (f, num_columns, Fef_Columns)))
     return NULL;

   c = r->col_data;
   if (-1 == check_repeat (c, num_columns))
     {
	jdfits_bintable_close_rows (r);
	return NULL;
     }

   return r;
}

#define SQRT_2		1.4142135623730951
#define SQRT_2PI	2.5066282746310002

/* It appears that the ciao processing of the FEF does not include the
 * sqrt(2*PI*sigma^2) in the denominator of the gaussian.  In other words, the
 * amplitude to does include this correction.
 */
static double gaussian_integral (double xmin, double xmax, Gauss_Parm_Type *g)
{
   double sigma, x0;

   if (0.0 == (sigma = g->sigma * SQRT_2))
     return 0.0;

   x0 = g->center;
   /* return 0.5 * g->amp * (erf ((xmax-x0)/sigma) - erf((xmin-x0)/sigma)); */

   return 0.5 * g->amp * (erf ((xmax-x0)/sigma) - erf((xmin-x0)/sigma))
     * (SQRT_2PI*g->sigma);
}

static double evaluate_gaussian (double x, Gauss_Parm_Type *g)
{
   double sigma;

   sigma = g->sigma;
   if (sigma == 0.0)
     return 0.0;

   x = (x - g->center)/sigma;
   /* return g->amp * exp (-0.5*x*x) / (SQRT_2PI * sigma); */
   return g->amp * exp (-0.5*x*x);
}

/* All gaussians have a positive weight. */
static int compute_pha_with_pos_amps (Gauss_Parm_Type *gaussians, unsigned int num_gaussians,
				      double *phap)
{
   while (1)
     {
	double r = JDMrandom ();
	Gauss_Parm_Type *g = gaussians;
	Gauss_Parm_Type *gmax = g + num_gaussians;

	while (g < gmax)
	  {
	     double pha;
	     unsigned int count;

	     if (g->cum_area <= r)
	       {
		  g++;
		  continue;
	       }

	     if (g->use_tail_dist == 0)
	       {
		  count = 0;
		  do
		    {
		       pha = g->center + g->sigma * JDMgaussian_random ();
		       count++;
		    }
		  while ((pha < 0) && (count < 100));
	       }
	     else
	       {
		  /* The algorithm here was adapted from the GNU Scientific Library */
		  double u, v, x;
		  double s = (0 - g->center)/g->sigma;

		  do
		    {
		       u = JDMrandom ();
		       do
			 {
			    v = JDMrandom();
			 }
		       while (v == 0.0);
		       x = sqrt (s * s - 2 * log (v));
		    }
		  while (x * u > s);
		  pha = g->center + x * g->sigma;
	       }

	     if (pha < 0)
	       {
		  fprintf (stderr, "Failed to find a pha value\n");
		  break;
	       }

	     *phap = pha;
	     return 0;
	  }
     }
}

/* This function is a variation of the rejection technique. */
static int compute_pha_with_neg_amps (Gauss_Parm_Type *gaussians, unsigned int num_gaussians,
				      double *phap)
{
   int count = 0;
   while (count < 100)
     {
	double pha;
	Gauss_Parm_Type *g, *gmax;
	double pos_sum, sum;

	if (-1 == compute_pha_with_pos_amps (gaussians, num_gaussians, &pha))
	  return -1;

	pos_sum = sum = 0.0;
	g = gaussians;
	gmax = g + num_gaussians;
	while (g < gmax)
	  {
	     double dsum = evaluate_gaussian (pha, g);
	     sum += dsum;
	     if (dsum > 0)
	       pos_sum += dsum;
	     g++;
	  }
	if (JDMrandom () * pos_sum < sum)
	  {
	     *phap = pha;
	     return 0;
	  }
	count++;
     }
   return -1;
}

#define HAS_NEG_AMP_GAUSSIANS	1
#define HAS_NEG_POS_GAUSSIANS	2
#define HAS_TOTAL_NEG_AREA	4
#define MY_INFINITY 1e37
static int normalize_gaussians (Gauss_Parm_Type *gaussians, unsigned int num,
				int *flagsp, float *first_momentp)
{
   Gauss_Parm_Type *g, *gmax;
   double total_pos_area = 0.0;
   double total_neg_area = 0.0;
   double total_area;
   int flags = 0;
   double ave;

   g = gaussians;
   gmax = g + num;
   ave = 0;
   while (g < gmax)
     {
	double area1, area2;

	area1 = gaussian_integral (0, MY_INFINITY, g);
	area2 = gaussian_integral (-MY_INFINITY, 0, g);

	g->use_tail_dist = 0;

	if (area2 > area1)
	  {
	     double ratio = area1/area2;
	     if (ratio < 0.1)
	       g->use_tail_dist = 1;
	  }
	if (area1 >= 0)
	  total_pos_area += area1;
	else
	  total_neg_area -= area1;

	ave += area1 * g->center;
	g->cum_area = total_pos_area;
	g++;
     }

   if (total_pos_area <= total_neg_area)
     {
#if 0
	marx_message ("***WARNING: Gaussian parameters appear invalid: a region with response <= 0 or absurd parameters has been detected.\n");
#endif
	flags |= HAS_TOTAL_NEG_AREA;
	/* return 0; */
     }

   if (total_neg_area != 0.0)
     flags |= HAS_NEG_AMP_GAUSSIANS;

   g = gaussians;
   if (total_pos_area > 0) while (g < gmax)
     {
	/* Since we interpolate over pairs of gaussians, do not normalize
	 * each member of the pair separately.
	 */
	/* g->amp /= total_pos_area; */
	g->cum_area /= total_pos_area;
	g++;
     }

   total_area = total_pos_area - total_neg_area;
   if (total_area != 0)
     *first_momentp = ave/total_area;
   else
     *first_momentp = 0;

   *flagsp = flags;
   return 0;
}

static int analyse_fef (Fef_Type *f, int *flagsp)
{
   Gauss_Parm_Type *g;
   unsigned int i;
   unsigned int num_energies;
   unsigned int num_gaussians;
   int flags;

   g = f->gaussians;
   num_gaussians = f->num_gaussians;
   num_energies = f->num_energies;

   *flagsp = 0;
   for (i = 0; i < num_energies; i++)
     {
	float mean;

	if (-1 == normalize_gaussians (g, num_gaussians, &flags, &mean))
	  return -1;
	if (flags & HAS_TOTAL_NEG_AREA)
	  *flagsp |= HAS_TOTAL_NEG_AREA;
	/* f->channels[i] = mean; */
	g += num_gaussians;
     }

   return 0;
}

typedef struct _Row_Type
{
   float energy;
   float channel;
   Gauss_Parm_Type *gaussians;
   struct _Row_Type *next;
}
Row_Type;

#define FWHM_TO_SIGMA 0.4246609001440095
static int process_region_rows (Row_Type *r, unsigned int num_energies,
				unsigned int num_gaussians,
				int ccdid, int xmin, int ymin,
				int xmax, int ymax, int region_num)
{
   Fef_Type *f;
   unsigned int i;
   int flags;

   (void) region_num;

   f = (Fef_Type *) marx_malloc (sizeof (Fef_Type));
   if (f == NULL)
     return -1;

   memset ((char *) f, 0, sizeof (Fef_Type));

   if ((NULL == (f->energies = (float *) marx_malloc (num_energies * sizeof (float))))
       || (NULL == (f->channels = (float *) marx_malloc (num_energies * sizeof (float))))
       || (NULL == (f->gaussians = (Gauss_Parm_Type *) marx_malloc (num_gaussians * num_energies * sizeof (Gauss_Parm_Type)))))
     {
	free_fef_type (f);
	return -1;
     }
   f->num_energies = num_energies;
   f->num_gaussians = num_gaussians;

   for (i = 0; i < num_energies; i++)
     {
	Gauss_Parm_Type *g, *g1, *gmax;

	f->energies[i] = r->energy;
	f->channels[i] = r->channel;
	g = f->gaussians + i * num_gaussians;
	gmax = g + num_gaussians;
	g1 = r->gaussians;
	while (g < gmax)
	  {
	     g->amp = g1->amp;
	     g->center = g1->center;
	     g->sigma = g1->sigma;
	     g->use_tail_dist = 0;
	     g++;
	     g1++;
	  }

	r = r->next;
     }

   if (-1 == check_fef_validity (f))
     {
	free_fef_type (f);
	return -1;
     }

   if (-1 == map_fef_to_region (ccdid, xmin, ymin, xmax, ymax, f))
     {
	free_fef_type (f);
	return -1;
     }

   if (-1 == analyse_fef (f, &flags))
     {
	free_fef_type (f);
	return -1;
     }
#if 0
   if (flags & HAS_TOTAL_NEG_AREA)
     marx_message ("*** WARNING: One or more energies in region %d has an invalid response\n", region_num);
#endif
   return 0;
}

static void free_row_chain (Row_Type *r)
{
   while (r != NULL)
     {
	Row_Type *n;

	n = r->next;
	if (r->gaussians != NULL)
	  marx_free ((char *) r->gaussians);
	marx_free ((char *) r);
	r = n;
     }
}

static int process_rows (JDFits_Type *f, unsigned int num_columns,
			 unsigned int num_gaussians,
			 int min_ccdid, int max_ccdid)
{
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   unsigned int i;
   Row_Type *root, *last;
   int last_region_num;
   int last_xmin, last_xmax, last_ymin, last_ymax, last_ccdid;
   unsigned int last_i;
   unsigned int num_energies, num_rows;

   if (NULL == (r = open_rows (f, num_columns)))
     return -1;

   root = last = NULL;
   last_xmin = last_xmax = last_ymin = last_ymax = -1;
   last_region_num = -1;
   num_energies = 0;
   last_i = 0;
   last_ccdid = -1;
   c = r->col_data;
   num_rows = r->num_rows;
   for (i = 0; i < num_rows; i++)
     {
	int regnum;
	Row_Type *rt;
	unsigned int j, k;
	int ccdid;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Unexpected end of file.");
	     goto return_error_bad_row;
	  }

	ccdid = c[CCDID_COLUMN].data.i[0];
	if ((ccdid < min_ccdid) || (ccdid > max_ccdid))
	  continue;

	rt = (Row_Type *) marx_malloc (sizeof (Row_Type));
	if (rt == NULL)
	  goto return_error_bad_row;
	if (NULL == (rt->gaussians = (Gauss_Parm_Type *) marx_malloc (num_gaussians * sizeof(Gauss_Parm_Type))))
	  {
	     marx_free ((char *) rt);
	     goto return_error_bad_row;
	  }

	rt->energy = c[ENERGY_COLUMN].data.f[0];
	rt->channel = c[CHANNEL_COLUMN].data.f[0];
	rt->next = NULL;

	k = FWHM1_COLUMN;
	for (j = 0; j < num_gaussians; j++)
	  {
	     Gauss_Parm_Type *g = rt->gaussians + j;
	     g->sigma = FWHM_TO_SIGMA * c[k++].data.f[0];
	     g->center = c[k++].data.f[0];
	     g->amp = c[k++].data.f[0];
	  }
	regnum = c[REGNUM_COLUMN].data.i[0];
	if ((regnum == last_region_num)
	    && (root != NULL))
	  {
	     last->next = rt;
	     last = rt;
	     num_energies++;
	     continue;
	  }

	if ((root != NULL)
	    && (-1 == process_region_rows (root, num_energies, num_gaussians,
					   last_ccdid, last_xmin, last_ymin,
					   last_xmax, last_ymax, last_region_num)))
	  {
	     free_row_chain (rt);
	     goto return_error_bad_row;
	  }

	free_row_chain (root);
	root = last = rt;
	num_energies = 1;
	last_region_num = regnum;
	last_xmin = c[CHIPX_LO_COLUMN].data.i[0];
	last_ymin = c[CHIPY_LO_COLUMN].data.i[0];
	last_xmax = c[CHIPX_HI_COLUMN].data.i[0];
	last_ymax = c[CHIPY_HI_COLUMN].data.i[0];
	last_ccdid = ccdid;
	last_i = i;
     }

   if (root != NULL)
     {
	int status =  process_region_rows (root, num_energies, num_gaussians,
					   last_ccdid, last_xmin, last_ymin,
					   last_xmax, last_ymax, last_region_num);
	free_row_chain (root);
	if (status == -1)
	  goto return_error_bad_row;
     }

   jdfits_bintable_close_rows (r);

   return check_fef_map_for_holes (min_ccdid, max_ccdid);

   return_error_bad_row:
   marx_error ("Error detected processing rows %u-%u\n", last_i+1, i+1);
   jdfits_bintable_close_rows (r);
   return -1;
}

static int read_fef_file (char *file, int min_ccdid, int max_ccdid)
{
   JDFits_Type *f;
   unsigned int i;
   char gauss_columns [3*MAX_GAUSSIANS][32];
   unsigned int num_gaussians, num_columns;

   free_fef_maps ();
   if (-1 == allocate_fef_maps (min_ccdid, max_ccdid))
     return -1;

   if (NULL == (f = open_fef_file (file)))
     return -1;

   /* Since we do not know how many gaussians there are in the file, it is
    * necessary to probe.
    */
   for (i = 0; i < MAX_FEF_COLUMNS - FWHM1_COLUMN; i++)
     {
	unsigned int j, ii, i3, i1;
	static char *fmts[3] =
	  {
	     "f:G%d_FWHM", "f:G%d_POS", "f:G%d_AMPL"
	  };

	i3 = 3 * i;
	i1 = i + 1;
	ii = i3 + FWHM1_COLUMN;
	for (j = 0; j < 3; j++)
	  {
	     char *colname = gauss_columns [i3 + j];
	     sprintf (colname, fmts[j], i1);

	     if (1 != jdfits_bintable_column_exists (f, colname+2))
	       break;

	     Fef_Columns[ii + j] = colname;
	  }
	if (j != 3)
	  break;
     }
   num_gaussians = i;
   num_columns = FWHM1_COLUMN + 3*num_gaussians;

   if (-1 == process_rows (f, num_columns, num_gaussians,
			   min_ccdid, max_ccdid))
     {
	marx_error ("Error processing ACIS fef file %s", file);
	free_fef_maps ();
	jdfits_close_file (f);
	return -1;
     }

   jdfits_close_file (f);
   return 0;
}

static int read_acis_fef (Param_File_Type *p, int min_ccdid, int max_ccdid)
{
   char *file;
   int status;

   (void) p;

   file = _marx_caldb_get_file ("ACISFEF");
   if (file == NULL)
     return -1;

   marx_message ("Reading ACIS-I/S FEF File\n");
   marx_message ("\t%s\n", file);
   status = read_fef_file (file, min_ccdid, max_ccdid);
   marx_free (file);
   return status;
}

int marx_init_acis_i_rmf (Param_File_Type *p)
{
   return read_acis_fef (p, 0, 3);
}

int marx_init_acis_s_rmf (Param_File_Type *p)
{
   return read_acis_fef (p, 4, 9);
}

#if MARX_HAS_IXO_SUPPORT
int marx_init_ixo_ccd_rmf (Param_File_Type *p)
{
   return read_acis_fef (p, 0, 0);
}
#endif

static unsigned int Num_Gaussians;
static Gauss_Parm_Type *Gaussians;

static Fef_Type *find_fef (int ccd_id, float x, float y)
{
   unsigned int i, j;
   Fef_Map_Type *m;
   Fef_Type *f;

   /* x=308; y=494; ccd_id=7; */

   if (x < 0)
     {
	if (_Marx_Det_Extend_Flag == 0)
	  return NULL;
	x = 0;
     }
   else if (x >= MIN_REGION_SIZE*NUM_REGIONS)
     {
	if (_Marx_Det_Extend_Flag == 0)
	  return NULL;
	x = (MIN_REGION_SIZE*NUM_REGIONS)-1;
     }

   if (y < 0)
     {
	if (_Marx_Det_Extend_Flag == 0)
	  return NULL;
	y = 0;
     }
   else if (y >= MIN_REGION_SIZE*NUM_REGIONS)
     {
	if (_Marx_Det_Extend_Flag == 0)
	  return NULL;
	y = (MIN_REGION_SIZE*NUM_REGIONS)-1;
     }

   i = (unsigned int) (x / MIN_REGION_SIZE);
   j = (unsigned int) (y / MIN_REGION_SIZE);

   if ((i >= NUM_REGIONS) || (j >= NUM_REGIONS))
     return NULL;

   if ((ccd_id < 0) || (ccd_id >= NUM_ACIS_CCDS))
     {
	marx_error ("find_fef: ccdid is out of range\n");
	return NULL;
     }

   m = Fef_Maps [ccd_id];
   if ((m == NULL)
       || (NULL == (f = m->fef_map[i][j])))
     {
	marx_error ("find_fef: No fef for region\n");
	return NULL;
     }

   return f;
}

int marx_apply_acis_rmf (int ccd_id, float x, float y,
			 double energy, float *pip, short *phap)
{
   unsigned int i, j;
   Fef_Type *f;
   double pha;
   unsigned int num_energies, num_gaussians;
   Gauss_Parm_Type *g0, *g1;
   double t;
   int flags, status;
   float mean;			       /* unused */

   if (NULL == (f = find_fef (ccd_id, x, y)))
     return -1;

   num_gaussians = f->num_gaussians;
   if (Num_Gaussians < num_gaussians)
     {
	marx_free ((char *) Gaussians);
	Num_Gaussians = 0;
	Gaussians = (Gauss_Parm_Type *) marx_malloc (num_gaussians * sizeof (Gauss_Parm_Type));
	if (Gaussians == NULL)
	  return -1;
	Num_Gaussians = num_gaussians;
     }

   num_energies = f->num_energies;
   i = JDMbinary_search_f (energy, f->energies, num_energies);
   j = i + 1;
   if (j == num_energies)
     {
	/* extrapolate */
	i--;
	j--;
     }

   t = (energy - f->energies[i])/(f->energies[j] - f->energies[i]);
   g0 = f->gaussians + i * num_gaussians;
   g1 = g0 + num_gaussians;

   for (i = 0; i < num_gaussians; i++)
     {
	double v;
	Gaussians[i].center = g0->center + t * (g1->center - g0->center);
	v = g0->sigma + t * (g1->sigma - g0->sigma);
	if (v <= 0.0)
	  {
	     Gaussians[i].sigma = 0.0;
	     Gaussians[i].amp = 0.0;
	  }
	else
	  {
	     Gaussians[i].sigma = v;
	     v = g0->amp + t * (g1->amp - g0->amp);
	     if ((v < 0.0) && ((g1->amp > 0.0) || (g0->amp > 0.0)))
	       v = 0.0;
	     Gaussians[i].amp = v;
	  }

	g0++;
	g1++;
     }

   if (-1 == normalize_gaussians (Gaussians, num_gaussians, &flags, &mean))
     return -1;

   if (flags & HAS_TOTAL_NEG_AREA)
     {
	marx_message ("Warning occurred for ccdid=%d,chipx=%f,chipy=%f,energy=%f\n",
		      ccd_id, x, y, energy);
	return -1;
	/* flags &= ~HAS_TOTAL_NEG_AREA; */
     }

   if (flags == 0)
     status = compute_pha_with_pos_amps (Gaussians, num_gaussians, &pha);
   else
     status = compute_pha_with_neg_amps (Gaussians, num_gaussians, &pha);

   if (status == -1)
     return -1;

   /* At this point, pha is real-valued.  Now convert it to an integer.
    * The first pha bin has a PHA of 1.  The convention (which I do not favor)
    * is to regard the center of the bin as having a value of 1.0, with the
    * left edge at 0.5.  Hence, if the pha is written as an integer (ipha)
    * plus a fraction f, then pha=ipha+f.  If 0 < 0.5 < f, then the desired
    * integer is ipha.  However, if (0.5<=f<1), then the desired value is
    * ipha+1.
    *
    * So, it would seem that we want *phap = (short)(pha+0.5).  However this
    * leads to a bin-shift when comparing to mkrmf.  Better agreement is
    * obtained using a bin defined by ipha <= pha < ipha+1.
    */
   *phap = (short) pha;

   /* The convention is that the first PHA bin is at ipha=1.  We want to
    * randomize within the bin to compute energy value.  I feel that this
    * calculation should go into marx2fits, but we do it here so that the
    * marx output files can be used.
    */
   pha = *phap - JDMrandom();

   *pip = JDMinterpolate_f (pha, f->channels, f->energies, f->num_energies);
   if (*pip < 0)
     return -1;

   return 0;
}

int _marx_apply_acis_rmf (_Marx_Acis_Chip_Type *c, float x, float y,
			  double energy, float *pip, short *phap)
{
   return marx_apply_acis_rmf (c->ccd_id, x, y, energy, pip, phap);
}

int marx_map_energy_to_acis_pha (int ccd_id, int x, int y, double energy, short *phap)
{
   Fef_Type *f;

   if (NULL == (f = find_fef (ccd_id, x, y)))
     return -1;

   *phap = JDMinterpolate_f (energy, f->energies, f->channels, f->num_energies);
   return 0;
}

#endif 				       /* MARX_HAS_ACIS_FEF */
