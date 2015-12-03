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
# include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"

/* distance from HRC-S geometric center to bottom edge */
static double Drake_X_Offset = 28.3;
static double Drake_Z_Offset = 5.0;
static double Drake_Length = 100.0;

/* Distance from geometric center to start of plates in y direction */
static double Drake_Gap_Y1 = 27.5;
static double Drake_Gap_Y2 = 27.5;

/* The Drake_Cr_Width corresponds to width of strip of chromium at the
 * edge near the gap.
 */
static double Drake_Cr_Width = 15.7;

#define MAX_PLATES 4
static double Drake_Heights[MAX_PLATES];
static double Drake_Angles[MAX_PLATES];

typedef struct /*{{{*/
{
   JDMVector_Type a;		       /* position of one corner */
   JDMVector_Type a1;		       /* position of one corner */
   JDMVector_Type e1, e2;	       /* unit vectors along side */
   JDMVector_Type normal;
   double len1, len2;		       /* length of sides */

}

/*}}}*/
Rectangle_Type;

/* Note: each flat consists of 2 sections */
static Rectangle_Type Drake_Flats [2 * MAX_PLATES];

static int Drake_N_Plates = 2;

/* Arrays of optical constants.  The first group applies to carbon and the
 * second group applies to the Cr strip.
 */
static float *Betas_C;
static float *Deltas_C;
static float *Energies_C;
static unsigned int Num_Energies_C;
static char *Drake_C_Opt_File;

static float *Betas_Cr;
static float *Deltas_Cr;
static float *Energies_Cr;
static unsigned int Num_Energies_Cr;
static char *Drake_Cr_Opt_File;

static Param_Table_Type Drake_Parm_Table [] = /*{{{*/
{
   {"HESFOffsetX",	PF_REAL_TYPE,		&Drake_X_Offset},
   {"HESFOffsetZ",	PF_REAL_TYPE,		&Drake_Z_Offset},
   {"HESFGapY1",	PF_REAL_TYPE,		&Drake_Gap_Y1},
   {"HESFGapY2",	PF_REAL_TYPE,		&Drake_Gap_Y2},
   {"HESFN",		PF_INTEGER_TYPE,	&Drake_N_Plates},
   {"HESFOptConstC",	PF_FILE_TYPE,		&Drake_C_Opt_File},
   {"HESFOptConstCr",	PF_FILE_TYPE,		&Drake_Cr_Opt_File},
   {"HESFLength",	PF_REAL_TYPE,		&Drake_Length},
   {"HESFCrWidth",	PF_REAL_TYPE, 		&Drake_Cr_Width},
   {NULL, 0, NULL}
};

/*}}}*/

static void free_optical_constants (void) /*{{{*/
{
   if (NULL != Betas_C) JDMfree_float_vector (Betas_C);
   if (NULL != Deltas_C) JDMfree_float_vector (Deltas_C);
   if (NULL != Energies_C) JDMfree_float_vector (Energies_C);
   Betas_C = Deltas_C = Energies_C = NULL;

   if (NULL != Betas_Cr) JDMfree_float_vector (Betas_Cr);
   if (NULL != Deltas_Cr) JDMfree_float_vector (Deltas_Cr);
   if (NULL != Energies_Cr) JDMfree_float_vector (Energies_Cr);
   Betas_Cr = Deltas_Cr = Energies_Cr = NULL;
}

/*}}}*/

static int read_opt_constants (char *file,
			       float **en, float **betas, float **deltas,
			       unsigned int *num_energies)
{
   unsigned int nread;

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;

   marx_message ("\t%s\n", file);

   /* The optical constant file consists of:
    *   energy (KeV), beta, delta
    */
   if (-1 == marx_f_read_bdat (file, &nread, 3, en, betas, deltas))
     {
	marx_free (file);
	return -1;
     }

   marx_free (file);

   *num_energies = (unsigned int) nread;

   return 0;
}

static int read_drake_opt_constants (void)
{
   free_optical_constants ();

   marx_message ("Reading Drake QE binary data files:\n");

   if ((Drake_C_Opt_File != NULL)
       && (*Drake_C_Opt_File != 0))
     {
	if (-1 == read_opt_constants (Drake_C_Opt_File,
				      &Energies_C, &Betas_C, &Deltas_C,
				      &Num_Energies_C))
	  return -1;
     }

   if ((Drake_Cr_Opt_File != NULL)
       && (*Drake_Cr_Opt_File != 0))
     {
	if (-1 == read_opt_constants (Drake_Cr_Opt_File,
				      &Energies_Cr, &Betas_Cr, &Deltas_Cr,
				      &Num_Energies_Cr))
	  return -1;
     }

   return 0;
}

int _marx_drake_flat_init (Param_File_Type *p) /*{{{*/
{
   char name[80];
   int i;
   double corner_x, corner_z;

   if ((-1 == pf_get_parameters (p, Drake_Parm_Table))
       || (Drake_C_Opt_File == NULL)
       || (Drake_Cr_Opt_File == NULL)
       || (Drake_N_Plates <= 0)
       || (Drake_N_Plates > MAX_PLATES))
     return -1;

   for (i = 0; i < Drake_N_Plates; i++)
     {
	double theta;

	sprintf (name, "HESFHeight%d", i + 1);
	if (-1 == pf_get_double (p, name, i + Drake_Heights))
	  return -1;

	sprintf (name, "HESFTheta%d", i + 1);
	if (-1 == pf_get_double (p, name, &theta))
	  return -1;
	Drake_Angles [i] = theta * (PI/180.0);
	/* Convert to radians */
     }

   if (-1 == read_drake_opt_constants ())
     return -1;

   /* Now massage the geometric information into normals and rectangular
    * regions.
    */

   corner_x = Drake_X_Offset;
   corner_z = Drake_Z_Offset;

   for (i = 0; i < Drake_N_Plates; i++)
     {
	Rectangle_Type *rect;
	double h, h_tan_theta;

	rect = Drake_Flats + i;
	h = Drake_Heights[i];
	h_tan_theta = h * tan (Drake_Angles[i]);

	rect->e1 = JDMv_vector (0.0, 1.0, 0.0);
	rect->len1 = Drake_Length;
	rect->e2 = JDMv_vector (h, 0.0, -h_tan_theta);
	rect->len2 = JDMv_length (rect->e2);
	rect->e2 = JDMv_unit_vector (rect->e2);
	rect->normal = JDMv_cross_prod (rect->e1, rect->e2);

	/* Corner position */
	rect->a = JDMv_vector (corner_x, Drake_Gap_Y1, corner_z);
	rect->a = JDMv_sum (_Marx_HRC_Geometric_Center, rect->a);

	rect->a1 = JDMv_ax1_bx2 (1.0, rect->a, rect->len1, rect->e1);

	/* Now do other one in pair */
	rect += Drake_N_Plates;

	rect->e1 = JDMv_vector (0.0, -1.0, 0.0);
	rect->len1 = Drake_Length;
	rect->e2 = JDMv_vector (h, 0.0, -h_tan_theta);
	rect->len2 = JDMv_length (rect->e2);
	rect->e2 = JDMv_unit_vector (rect->e2);
	/* Note that since we want the outward normal, reverse cross product */
	rect->normal = JDMv_cross_prod (rect->e2, rect->e1);

	/* corner position */
	rect->a = JDMv_vector (corner_x, -Drake_Gap_Y2, corner_z);
	rect->a = JDMv_sum (_Marx_HRC_Geometric_Center, rect->a);

	rect->a1 = JDMv_ax1_bx2 (1.0, rect->a, rect->len1, rect->e1);

	/* Now move to next corner */
	corner_z -= h_tan_theta;
	corner_x += h;
     }

   return 0;
}

/*}}}*/

static Rectangle_Type *drake_intersection (JDMVector_Type *x, JDMVector_Type *p, /*{{{*/
					   double *cos_theta, int *use_cr)
{
   unsigned int i, imax;
   Rectangle_Type *rect;
   double t, p_dot_n;
   JDMVector_Type new_x, old_x, old_p, x_prime;

   old_x = *x; old_p = *p;

   imax = 2 * Drake_N_Plates;

   for (i = 0; i < imax; i++)
     {
	double xx, yy;
	rect = Drake_Flats + i;

	if (0.0 == (p_dot_n = JDMv_pdot_prod (&old_p, &rect->normal)))
	  continue;

	t = JDMv_dot_prod (JDMv_diff (rect->a, old_x), rect->normal) / p_dot_n;

	new_x = JDMv_ax1_bx2 (1.0, old_x, t, old_p);
	x_prime = JDMv_diff (new_x, rect->a);

	xx = JDMv_pdot_prod (&x_prime, &rect->e1);
	if ((xx < 0.0) || (xx >= rect->len1))
	  continue;

	yy = JDMv_pdot_prod (&x_prime, &rect->e2);
	if ((yy < 0.0) || (yy >= rect->len2))
	  continue;

	/* We have a hit. */

	if (xx < Drake_Cr_Width)
	  *use_cr = 1;
	else
	  *use_cr = 0;

	*x = new_x;
	*cos_theta = p_dot_n;
	return rect;
     }
   return NULL;
}

/*}}}*/

int _marx_drake_reflect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int n, i, *sorted_index;

   marx_prune_photons (pt);

   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   for (i = 0; i < n; i++)
     {
	Rectangle_Type *plate;
	double rfl;
	double p_dot_n;
	int use_cr;

	at = photon_attributes + sorted_index[i];
	if (at->flags & BAD_PHOTON_MASK) continue;

	if (NULL == (plate = drake_intersection (&at->x, &at->p, &p_dot_n, &use_cr)))
	  continue;

	/* Now we have intersected at at->x.  Compute new direction
	 * but only if the reflectivity permits.
	 */
	if (use_cr)
	  {
	     if (Betas_Cr != NULL)
	       rfl = marx_interp_reflectivity (at->energy, p_dot_n, Energies_Cr, Betas_Cr, Deltas_Cr, Num_Energies_Cr);
	     else rfl = 1.0;
	  }
	else
	  {
	     if (Betas_C != NULL)
	       rfl = marx_interp_reflectivity (at->energy, p_dot_n, Energies_C, Betas_C, Deltas_C, Num_Energies_C);
	     else rfl = 1.0;
	  }

	if (rfl < JDMrandom ())
	  {
	     at->flags |= PHOTON_DRAKE_BLOCKED;
	     continue;
	  }

	/* New direction */
	at->p = JDMv_ax1_bx2 (1.0, at->p, -2.0 * p_dot_n, plate->normal);
	at->flags |= PHOTON_DRAKE_REFLECTED;
     }

   return 0;
}

/*}}}*/

