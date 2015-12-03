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
#include <jdfits.h>

#include "marx.h"
#include "_marx.h"

#if MARX_HAS_ACIS_GAIN_MAP

/* The ACIS gain table is a single fits extension that consists of columns:
 * CHIPX_MIN/MAX, CHIPY_MIN/MAX, PHA[], ENERGY[], SIGMA[], NPOINTS
 * where NPOINTS specifies the length of the PHA, ENERGY, and SIGMA
 * vectors.
 *
 * Unfortunately, the regions defined by the CHIPX/Y_MIN/MAX values are not
 * uniform in size.  However, the values are _supposed_ to be multiples
 * of 32 and this is what will be assumed in this code.  Hence, at most
 * there will be 32x32=1024 regions per CCD.
 *
 * This file is ugly because with FITS, you never know what you are going to
 * get.  So you have got to be ready for anything.
 */
#define MIN_REGION_SIZE	32
#define NUM_REGIONS (1024/MIN_REGION_SIZE)
#define GAIN_MAP_EXTNAME "AXAF_DETGAIN"

#define NUM_GAIN_COLUMNS	9
static char *Gain_Columns [NUM_GAIN_COLUMNS] =
{
   "i:CCD_ID", "i:NPOINTS", "f:ENERGY", "f:PHA", "f:SIGMA",
   "i:CHIPX_MIN", "i:CHIPX_MAX", "i:CHIPY_MIN", "i:CHIPY_MAX"
};
#define CCDID_COLUMN	0
#define NPOINTS_COLUMN	1
#define ENERGY_COLUMN	2
#define PHA_COLUMN	3
#define SIGMA_COLUMN	4
#define CHIPXMIN_COLUMN	5
#define CHIPXMAX_COLUMN	6
#define CHIPYMIN_COLUMN	7
#define CHIPYMAX_COLUMN	8

typedef struct
{
   unsigned int num_points;
   float *energies;		       /* malloced */
   float *phas;
   float *sigmas;
}
Gain_Type;

#define NUM_ACIS_CCDS	10
static Gain_Type *Gain_Maps;
static unsigned int Num_Gain_Maps;

typedef struct
{
   Gain_Type *map[NUM_REGIONS][NUM_REGIONS];
}
Gain_Map_Type;

static Gain_Map_Type *Gain_Map_Tables[NUM_ACIS_CCDS];

static void free_gain_lookup_tables (unsigned int min_ccdid, unsigned int max_ccdid)
{
   unsigned int i;

   for (i = min_ccdid; i <= max_ccdid; i++)
     {
	if (Gain_Map_Tables[i] == NULL)
	  continue;

	marx_free ((char *) Gain_Map_Tables[i]);
	Gain_Map_Tables[i] = NULL;
     }
}

static int allocate_gain_lookup_tables (unsigned int min_ccdid, unsigned max_ccdid)
{
   unsigned int i;

   if (min_ccdid >= NUM_ACIS_CCDS)
     {
	marx_error ("Internal Error: min/max_ccdid > NUM_ACIS_CCDS");
	return -1;
     }

   free_gain_lookup_tables (min_ccdid, max_ccdid);

   for (i = min_ccdid; i <= max_ccdid; i++)
     {
	Gain_Map_Type *g;

	g = (Gain_Map_Type *) marx_malloc (sizeof (Gain_Map_Type));
	if (g == NULL)
	  {
	     free_gain_lookup_tables (min_ccdid, max_ccdid);
	     return -1;
	  }
	memset ((char *) g, 0, sizeof (Gain_Map_Type));
	Gain_Map_Tables[i] = g;
     }

   return 0;
}

static int allocate_gain_type (Gain_Type *g, unsigned int npoints)
{
   if (NULL == (g->energies = (float *)marx_malloc (3 * npoints * sizeof (float))))
     return -1;

   g->phas = g->energies + npoints;
   g->sigmas = g->phas + npoints;
   g->num_points = npoints;

   return 0;
}

static void free_gain_maps (void)
{
   Gain_Type *g, *gmax;

   if (Gain_Maps == NULL)
     return;
   g = Gain_Maps;
   gmax = g + Num_Gain_Maps;

   while (g < gmax)
     {
	if (g->energies != NULL)
	  marx_free ((char *) g->energies);
	g++;
     }
   marx_free ((char *) Gain_Maps);
   Gain_Maps = NULL;
   Num_Gain_Maps = 0;
}

static int map_gain_to_region (int ccdid, int xmin, int ymin, int xmax, int ymax,
			       Gain_Type *g)
{
   Gain_Map_Type *m;
   int i;

   if ((ccdid < 0) || (ccdid >= NUM_ACIS_CCDS))
     {
	marx_error ("Invalid CCD_ID value in gain file: %d", ccdid);
	return -1;
     }

   if (xmin < 1) xmin = 1;
   if (ymin < 1) ymin = 1;
   if (xmax > 1024) xmax = 1024;
   if (ymax > 1024) ymax = 1024;

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

   m = Gain_Map_Tables[ccdid];
   if (m == NULL)
     {
	marx_error ("Internal Error: Gain_Map_Table[%d] is NULL", ccdid);
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
	     if (m->map[i][j] != NULL)
	       {
		  marx_message ("***Warning: Gain Map defines an overlapping region for CCD_ID=%d at (%d, %d)\n",
				ccdid, 1+i*MIN_REGION_SIZE, 1+j*MIN_REGION_SIZE);
	       }
	     m->map[i][j] = g;
	  }
     }

   return 0;
}

static int check_gain_map_for_holes (int min_ccdid, int max_ccdid)
{
   int ccdid;

   for (ccdid = min_ccdid; ccdid <= max_ccdid; ccdid++)
     {
	Gain_Map_Type *m = Gain_Map_Tables[ccdid];
	unsigned int i, j;

	if (m == NULL)
	  {
	     marx_error ("Internal Error: Gain_Map_Tables[%d] is NULL", ccdid);
	     return -1;
	  }

	for (i = 0; i < NUM_REGIONS; i++)
	  {
	     for (j = 0; j < NUM_REGIONS; j++)
	       {
		  if (m->map[i][j] != NULL)
		    continue;

		  marx_message ("***Warning: Gain Map for CCDID=%d contains holes\n",
				ccdid);
		  i = NUM_REGIONS;
		  break;
	       }
	  }
     }

   return 0;
}

static int check_monotonicity (float *p, unsigned int n)
{
   float x, *pmax;

   pmax = p + n;
   x = *p++;
   while (p < pmax)
     {
	if (*p < x)
	  return -1;
	x = *p++;
     }
   return 0;
}

static int check_gain_validity (Gain_Type *g)
{
   float *p, *pmax;

   if (-1 == check_monotonicity (g->energies, g->num_points))
     {
	marx_error ("Gain Map does not have monotonically increasing energies");
	return -1;
     }
   if (-1 == check_monotonicity (g->phas, g->num_points))
     {
	marx_error ("Gain Map does not have monotonically increasing PHAs");
	return -1;
     }

   p = g->energies;
   pmax = p + g->num_points;
   while (p < pmax)
     {
	*p *= 1e-3;		       /* convert to keV */
	p++;
     }

   return 0;
}

static int check_repeat (JDFits_Col_Data_Type *c, int n)
{
   c += n;
   if (c->repeat != 1)
     {
	marx_error ("Expecting a repeat count of 1 for column %s", Gain_Columns[n]+2);
	return -1;
     }
   return 0;
}

static int read_gain_map (char *file, int min_ccdid, int max_ccdid)
{
   JDFits_Type *f;
   JDFits_Row_Type *r;
   unsigned int i, num_rows;
   JDFits_Col_Data_Type *c;
   unsigned int max_repeat;

   free_gain_maps ();
   if (-1 == allocate_gain_lookup_tables (min_ccdid, max_ccdid))
     return -1;

   if (NULL == (f = jdfits_open_binary_table (file, GAIN_MAP_EXTNAME)))
     {
	marx_error ("Unable to open ACIS gain table %s", file);
	return -1;
     }
   r = jdfits_bintable_aopen_rows (f, NUM_GAIN_COLUMNS, Gain_Columns);
   if (r == NULL)
     {
	marx_error ("Error processing ACIS gain table %s", file);
	jdfits_close_file (f);
	return -1;
     }

   if (0 == (num_rows = r->num_rows))
     {
	marx_error ("ACIS gain table %s contains NO data.", file);
	goto return_error;
     }
   c = r->col_data;

   if ((-1 == check_repeat (c, CCDID_COLUMN))
       || (-1 == check_repeat (c, NPOINTS_COLUMN))
       || (-1 == check_repeat (c, CHIPXMIN_COLUMN))
       || (-1 == check_repeat (c, CHIPXMAX_COLUMN))
       || (-1 == check_repeat (c, CHIPYMIN_COLUMN))
       || (-1 == check_repeat (c, CHIPYMAX_COLUMN)))
     goto return_error;

   max_repeat = c[ENERGY_COLUMN].repeat;
   if (max_repeat < c[PHA_COLUMN].repeat)
     max_repeat = c[PHA_COLUMN].repeat;
   if (max_repeat < c[SIGMA_COLUMN].repeat)
     max_repeat = c[SIGMA_COLUMN].repeat;

   if (NULL == (Gain_Maps = (Gain_Type *)marx_malloc (num_rows * sizeof (Gain_Type))))
     goto return_error;
   memset ((char *) Gain_Maps, 0, num_rows * sizeof (Gain_Type));
   Num_Gain_Maps = num_rows;

   for (i = 0; i < num_rows; i++)
     {
	int ccdid, npoints;
	int xmin, ymin, xmax, ymax;
	Gain_Type *g;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Unexpected end of ACIS gain table %s", file);
	     goto return_error_bad_row;
	  }
	c = r->col_data;

	ccdid = c[CCDID_COLUMN].data.i[0];
	if ((ccdid < min_ccdid)
	    || (ccdid > max_ccdid))
	  continue;

	npoints = c[NPOINTS_COLUMN].data.i[0];
	if (npoints <= 0)
	  continue;

	if (npoints == 1)
	  {
	     marx_error ("Expecting Number of points (NPOINTS) to be greater than 1");
	     goto return_error_bad_row;
	  }

	if (npoints > (int) max_repeat)
	  {
	     marx_error ("NPOINTS is greater than repeat count");
	     goto return_error_bad_row;
	  }

	g = Gain_Maps + i;

	if (-1 == allocate_gain_type (g, npoints))
	  goto return_error_bad_row;

	xmin = c[CHIPXMIN_COLUMN].data.i[0];
	ymin = c[CHIPYMIN_COLUMN].data.i[0];
	xmax = c[CHIPXMAX_COLUMN].data.i[0];
	ymax = c[CHIPYMAX_COLUMN].data.i[0];
	memcpy (g->energies, (char *)c[ENERGY_COLUMN].data.f, npoints * sizeof (float));
	memcpy (g->phas, (char *)c[PHA_COLUMN].data.f, npoints * sizeof (float));
	memcpy (g->sigmas, (char *)c[SIGMA_COLUMN].data.f, npoints * sizeof (float));

	if (-1 == check_gain_validity (g))
	  goto return_error_bad_row;

	if (-1 == map_gain_to_region (ccdid, xmin, ymin, xmax, ymax, g))
	  goto return_error_bad_row;
     }

   if (-1 == check_gain_map_for_holes (min_ccdid, min_ccdid))
     goto return_error;

   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   return 0;

   return_error_bad_row:
   marx_error ("Error encountered while processing row %d of %s\n", i+1, file);
   /* drop */

   return_error:
   free_gain_maps ();
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   return -1;
}

static int read_acis_gain_map (Param_File_Type *p, int min_ccdid, int max_ccdid)
{
   char buf [PF_MAX_LINE_LEN];
   char *file;
   int status;
   char *parm = "ACIS_Gain_Map_File";

   if (-1 == pf_get_file (p, parm, buf, sizeof (buf)))
     {
	marx_error ("Unable to get paramter %s", parm);
	return -1;
     }
   file = marx_make_data_file_name (buf);
   if (file == NULL)
     return -1;

   marx_message ("Reading ACIS-I/S Gain File\n");
   marx_message ("\t%s\n", file);
   status = read_gain_map (file, min_ccdid, max_ccdid);
   marx_free (file);
   return status;
}

int _marx_init_acis_i_gain_map (Param_File_Type *p)
{
   return read_acis_gain_map (p, 0, 3);
}

int _marx_init_acis_s_gain_map (Param_File_Type *p)
{
   return read_acis_gain_map (p, 4, 9);
}

short _marx_apply_acis_gain_map (_Marx_Acis_Chip_Type *c, float x, float y, double energy, float *pi)
{
   unsigned int i, j;
   Gain_Map_Type *m;
   Gain_Type *g;
   float pha, sigma;

   i = (unsigned int) (x / MIN_REGION_SIZE);
   j = (unsigned int) (y / MIN_REGION_SIZE);

   if ((i >= NUM_REGIONS) || (j >= NUM_REGIONS))
     return -1;

   if ((c->ccd_id < 0) || (c->ccd_id >= NUM_ACIS_CCDS))
     {
	marx_message ("_marx_apply_acis_gain_map: ccdid is out of range\n");
	return -1;
     }

   m = Gain_Map_Tables[c->ccd_id];
   if (m == NULL)
     return -1;

   g = m->map[i][j];
   if (g == NULL)
     return -1;

   pha = JDMinterpolate_f (energy, g->energies, g->phas, g->num_points);
   sigma = JDMinterpolate_f (energy, g->energies, g->sigmas, g->num_points);

   pha += sigma * JDMgaussian_random ();

   if (pha < 0.0)
     return -1;

   energy = JDMinterpolate_f (pha, g->phas, g->energies, g->num_points);
   if (energy < 0.0)
     return -1;

   *pi = energy;
   return (short) pha;
}

#endif				       /* MARX_HAS_ACIS_GAIN_MAP */
