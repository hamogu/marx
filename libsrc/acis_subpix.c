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
#include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include <pfile.h>
#include <jdfits.h>

#include "marx.h"
#include "_marx.h"

typedef struct
{
   float *energies;		       /* malloced */
   int num_energies;
   float *dxs;			       /* not malloced */
   float *dys;			       /* not malloced */
}
Subpix_Type;

typedef struct
{
#define NUM_FLIGHT_GRADES 256
   Subpix_Type *s[NUM_FLIGHT_GRADES];
}
Subpix_CCD_Type;

#define NUM_ACIS_CHIPS 10
struct _Marx_Subpix_Table_Type
{
   Subpix_CCD_Type *subpix_tables[NUM_ACIS_CHIPS];
   Subpix_CCD_Type *fi, *bi;
};

static Subpix_Type *allocate_subpix_type (int num_energies)
{
   Subpix_Type *s;

   if (NULL == (s = (Subpix_Type *) marx_malloc (sizeof(Subpix_Type))))
     return NULL;

   if (NULL == (s->energies = (float *)marx_malloc (3*sizeof(float)*num_energies)))
     {
	marx_free ((char *)s);
	return NULL;
     }
   s->num_energies = num_energies;
   s->dxs = s->energies + num_energies;
   s->dys = s->dxs + num_energies;
   return s;
}

static void free_subpix_type (Subpix_Type *s)
{
   if (s == NULL)
     return;

   if (s->energies != NULL)
     marx_free ((char *)s->energies);
   marx_free ((char *)s);
}

static void free_subpix_ccd_type (Subpix_CCD_Type *sccd)
{
   int i;

   if (sccd == NULL)
     return;

   for (i = 0; i < NUM_FLIGHT_GRADES; i++)
     free_subpix_type (sccd->s[i]);

   marx_free ((char *) sccd);
}

#define NUM_SUBPIX_COLUMNS	5
static char *Subpix_Columns [NUM_SUBPIX_COLUMNS] =
{
   "i:FLTGRADE", "i:NPOINTS", "f:ENERGY",
   "f:CHIPX_OFFSET", "f:CHIPY_OFFSET"
};
#define FLTGRADE_COLUMN		0
#define NPOINTS_COLUMN		1
#define ENERGY_COLUMN		2
#define CHIPX_OFFSET_COLUMN	3
#define CHIPY_OFFSET_COLUMN	4

static Subpix_CCD_Type *read_subpix_ext (char *file, char *extname)
{
   JDFits_Type *f;
   JDFits_Row_Type *r;
   unsigned int i, num_rows;
   JDFits_Col_Data_Type *c;
   unsigned int max_repeat;
   Subpix_CCD_Type *sccd = NULL;

   if (NULL == (f = jdfits_open_binary_table (file, extname)))
     {
	marx_error ("Unable to open ACIS subpix file %s[%s]", file, extname);
	return NULL;
     }
   r = jdfits_bintable_aopen_rows (f, NUM_SUBPIX_COLUMNS, Subpix_Columns);
   if (r == NULL)
     {
	marx_error ("Error processing ACIS subpix file %s[%s]", file, extname);
	jdfits_close_file (f);
	return NULL;
     }

   if (0 == (num_rows = r->num_rows))
     {
	marx_error ("ACIS subpix file %s contains NO data.", file);
	goto return_error;
     }
   c = r->col_data;

   if ((c[FLTGRADE_COLUMN].repeat != 1)
       || (c[NPOINTS_COLUMN].repeat != 1))
     {
	marx_error ("The repeat count for the fltgrade/npoints column is invalid");
	goto return_error;
     }

   max_repeat = c[ENERGY_COLUMN].repeat;
   if ((max_repeat != c[CHIPX_OFFSET_COLUMN].repeat)
       || (max_repeat != c[CHIPY_OFFSET_COLUMN].repeat))
     {
	marx_error ("The chipx/y_offset_column repeat count does not match the energy repeat count");
	goto return_error;
     }

   if (NULL == (sccd = (Subpix_CCD_Type *)marx_malloc (sizeof(Subpix_CCD_Type))))
     goto return_error;
   memset ((char *)sccd, 0, sizeof(Subpix_CCD_Type));

   for (i = 0; i < num_rows; i++)
     {
	int fltgrade, npoints;
	Subpix_Type *s;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Unexpected end of ACIS subpix table %s", file);
	     goto return_error_bad_row;
	  }
	c = r->col_data;

	fltgrade = c[FLTGRADE_COLUMN].data.i[0];
	npoints = c[NPOINTS_COLUMN].data.i[0];
	if ((fltgrade < 0) || (fltgrade >= NUM_FLIGHT_GRADES))
	  {
	     marx_error ("Illegal flight grade in subpix column");
	     goto return_error_bad_row;
	  }

	if (npoints <= 0)
	  {
	     sccd->s[fltgrade] = NULL;
	     continue;
	  }
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

	if (NULL == (s = allocate_subpix_type (npoints)))
	  goto return_error_bad_row;
	sccd->s[fltgrade] = s;

	memcpy (s->energies, (char *)c[ENERGY_COLUMN].data.f, npoints * sizeof (float));
	memcpy (s->dxs, (char *)c[CHIPX_OFFSET_COLUMN].data.f, npoints * sizeof (float));
	memcpy (s->dys, (char *)c[CHIPY_OFFSET_COLUMN].data.f, npoints * sizeof (float));
     }

   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   return sccd;

return_error_bad_row:
   marx_error ("Error encountered while processing row %d of %s\n", i+1, file);
   /* drop */

return_error:
   free_subpix_ccd_type (sccd);
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   return NULL;
}

void marx_close_acis_subpixel (Marx_Subpix_Table_Type *stt)
{
   if (stt == NULL)
     return;
   if (stt->bi != NULL)
     free_subpix_ccd_type (stt->bi);
   if (stt->fi != NULL)
     free_subpix_ccd_type (stt->fi);
   marx_free ((char *)stt);
}

static Marx_Subpix_Table_Type *read_subpixel_file (char *file)
{
   Marx_Subpix_Table_Type *stt;
   int i;

   if (NULL == (stt = (Marx_Subpix_Table_Type *)marx_malloc (sizeof(Marx_Subpix_Table_Type))))
     return NULL;
   memset ((char *)stt, 0, sizeof(Marx_Subpix_Table_Type));

   if (NULL == (stt->fi = read_subpix_ext (file, "MARX_ACIS_SUBPIX_FI")))
     return NULL;

   if (NULL == (stt->bi = read_subpix_ext (file, "MARX_ACIS_SUBPIX_BI")))
     {
	marx_close_acis_subpixel (stt);
	return NULL;
     }

   for (i = 0; i < NUM_ACIS_CHIPS; i++)
     {
	if ((i == 5) || (i == 7))
	  stt->subpix_tables[i] = stt->bi;
	else
	  stt->subpix_tables[i] = stt->fi;
     }
   return stt;
}

Marx_Subpix_Table_Type *marx_open_acis_subpix (void)
{
   char *file;
   Marx_Subpix_Table_Type *stt;

   if (NULL == (file = _marx_caldb_get_file ("ACISSUBPIX")))
     return NULL;

   marx_message ("Reading subpix file %s\n", file);
   stt = read_subpixel_file (file);
   marx_free (file);
   return stt;
}

/* In the CalDB EDSER file (in structure stt), find the energy and fltgrade
 * dependend dx, dy for a photon.
 */
int marx_compute_acis_subpix (Marx_Subpix_Table_Type *stt,
			      int ccd, float energy, int fltgrade, float *dxp, float *dyp)
{
   Subpix_Type *s;
   unsigned int i, j, n;
   double w0, w1;

   if (stt == NULL)
     return -1;

   if ((fltgrade < 0) || (fltgrade >= NUM_FLIGHT_GRADES))
     return -1;

   s = stt->subpix_tables[ccd]->s[fltgrade];
   if (s == NULL)
     {
	*dxp = *dyp = 0;
	return 0;
     }

   n = (unsigned int) s->num_energies;
   i = JDMbinary_search_f (energy, s->energies, n);
   j = i+1;
   if (j == n)
     {
	j = i;
	i--;
     }
   // Linear interpolation on the energy grid to get dxp, dyp
   w1 = ((double)energy - s->energies[i])/(s->energies[j] - s->energies[i]);
   w0 = 1.0 - w1;
   *dxp = w0 * s->dxs[i] + w1 * s->dxs[j];
   *dyp = w0 * s->dys[i] + w1 * s->dys[j];

   return 0;
}
