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
/* Mirror Blur module */

#include "config.h"
#include "marx-feat.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"

static int Perform_Blur = 0;

/* This structure defines theta vs photon energy and encircled energy */
typedef struct /*{{{*/
{
   unsigned int num_ee;
   unsigned int num_energies;
   float *energies;
   float *encircled_energies;
   float **thetas[MARX_NUM_MIRRORS];		       /* num_energies x num_ee */
}
/*}}}*/
Mirror_Blur_Type;

static Mirror_Blur_Type Mirr_Blur;

static void free_mirror_blurs (void) /*{{{*/
{
   unsigned int i;

   if (Mirr_Blur.energies != NULL)
     {
	JDMfree_float_vector (Mirr_Blur.energies);
	Mirr_Blur.energies = NULL;
     }

   if (Mirr_Blur.encircled_energies != NULL)
     {
	JDMfree_float_vector (Mirr_Blur.encircled_energies);
	Mirr_Blur.encircled_energies = NULL;
     }

   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	if (Mirr_Blur.thetas[i] != NULL)
	  {
	     JDMfree_float_matrix (Mirr_Blur.thetas[i], Mirr_Blur.num_energies);
	     Mirr_Blur.thetas[i] = NULL;
	  }
     }
}

/*}}}*/

int _marx_init_mirror_blur (Param_File_Type *p) /*{{{*/
{
   FILE *fp;
   unsigned int num_energies;
   unsigned int num_ee;
   unsigned int i;
   char filebuf[PF_MAX_LINE_LEN];
   char *file;

   if (-1 == pf_get_boolean (p, "MirrorBlur", &Perform_Blur))
     return -1;

   if (Perform_Blur == 0) return 0;

   if (-1 == pf_get_file (p, "MirrorBlurFile", filebuf, sizeof (filebuf)))
     {
	Perform_Blur = 0;
	return -1;
     }

   if (NULL == (file = marx_make_data_file_name (filebuf)))
     return -1;

   marx_message ("Reading blur data file:\n\t%s\n", file);

   if (NULL == (fp = fopen (file, "rb")))
     {
	marx_error ("Unable to open file: %s", file);
	marx_free (file);
	JDMath_Error = JDMATH_FILE_OPEN_ERROR;
	return -1;
     }

   if ((1 != JDMread_int32 ((int32 *)&num_energies, 1, fp))
       || (1 != JDMread_int32 ((int32 *)&num_ee, 1, fp)))
     {
	JDMath_Error = JDMATH_FILE_READ_ERROR;
	goto return_error;
     }

   marx_free (file);
   file = NULL;

   free_mirror_blurs ();

   /* Now allocate space for the blur data */
   Mirr_Blur.num_energies = num_energies;
   if (NULL == (Mirr_Blur.energies = JDMfloat_vector (num_energies)))
     goto return_error;

   if (num_energies != JDMread_float32(Mirr_Blur.energies, num_energies, fp))
     {
	JDMath_Error = JDMATH_FILE_READ_ERROR;
	goto return_error;
     }

   Mirr_Blur.num_ee = num_ee;
   if (NULL == (Mirr_Blur.encircled_energies = JDMfloat_vector (num_ee)))
     goto return_error;

   if (num_ee != JDMread_float32(Mirr_Blur.encircled_energies, num_ee, fp))
     {
	JDMath_Error = JDMATH_FILE_READ_ERROR;
	goto return_error;
     }

   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	unsigned int n;
	if (NULL == (Mirr_Blur.thetas[i] = JDMfloat_matrix (num_energies, num_ee)))
	  goto return_error;

	for (n = 0; n < num_energies; n++)
	  {
	     if (num_ee != JDMread_float32(Mirr_Blur.thetas[i][n], num_ee, fp))
	       {
		  JDMath_Error = JDMATH_FILE_READ_ERROR;
		  goto return_error;
	       }
	  }
     }
   fclose (fp);
   return 0;

   /* Get here only if something goes wrong. */
   return_error:
   if (file != NULL) marx_free (file);
   free_mirror_blurs ();
   fclose (fp);
   return -1;
}

/*}}}*/

static void blur_photon (Marx_Photon_Attr_Type *at, double blur) /*{{{*/
{
   double alpha, beta, px;
   JDMVector_Type *p;

   p = &(at->p);
   px = p->x;

   /* Note: by construction in calling routines, px is NOT 1.0 */
   beta = sin (blur) / sqrt (1.0 - px * px);
   alpha = cos (blur) - beta * px;

   p->x = alpha * px + beta;
   p->y = alpha * p->y;
   p->z = alpha * p->z;
}

/*}}}*/

int _marx_mirror_blur (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int n, i, *sorted_index;
   float *encircled_energies, *energies;
   unsigned int num_ee, num_energies;

   if (Perform_Blur == 0) return 0;

   marx_prune_photons (pt);
   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   encircled_energies = Mirr_Blur.encircled_energies;
   num_ee = Mirr_Blur.num_ee;
   energies = Mirr_Blur.energies;
   num_energies = Mirr_Blur.num_energies;

   for (i = 0; i < n; i++)
     {
	unsigned int j;
	float energy, r, theta0, theta1;
	float e0, e1;
	double blur;

	at = photon_attributes + sorted_index[i];

	energy = (float) at->energy;
	/* Now find which two energy curves to use: j and j + 1 */
	j = JDMbinary_search_f (energy, energies, num_energies);
	e0 = energies[j];
	e1 = energies[j + 1];

	/* Ideally a random number between 0 and 1 would be used.  However,
	 * the blur data close to 1.0 cannot be trusted.
	 */
	r = (float) 0.99 * JDMrandom ();

	/* Now interpolate thetas on the two energy curves */
	theta0 = JDMinterpolate_f (r, encircled_energies,
				   Mirr_Blur.thetas[at->mirror_shell][j],
				   num_ee);
	theta1 = JDMinterpolate_f (r, encircled_energies,
				   Mirr_Blur.thetas[at->mirror_shell][j + 1],
				   num_ee);
	blur = theta0 + (theta1 - theta0) * (energy - e0) / (e1 - e0);

	/* pick the sign of the blur angle */
	if (2.0 * JDMrandom () < 1.0) blur = -blur;
	blur_photon (at, blur);
     }
   return 0;
}

/*}}}*/
