/* -*- mode: C; mode: fold -*- */
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

/*{{{ #include files */

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

static double Mirror_Vignetting_Factor = 1.0;

/* This distance is distance from 1/2 H segment on HRMA to Focus. */
static double HRMA_Position = (1000.0 * (10.0692 - 0.83 / 2.0));   /* mm */
static char *HRMA_EA_File;

/* These Numbers are taken from SIN */
static double Mirror_Radii[MARX_NUM_MIRRORS] = /*{{{*/
{
   600.0, 480.0, 425.0, 310.0	       /* mm */
};

/*}}}*/

/* The mirror shutters consist of a string characters (0,1) that represent
 * a bitmapped array for the mirror shutters.  1 means the shutter is closed
 * for that quadrant.
 */
static char *Mirror_Shutters [MARX_NUM_MIRRORS];

static int Mirror_Use_Effective_Area;
static Param_Table_Type HRMA_Parm_Table [] = /*{{{*/
{
   {"MirrorVig",	PF_REAL_TYPE,		&Mirror_Vignetting_Factor},
   {"MirrorRadius1",	PF_REAL_TYPE,	 	&Mirror_Radii[0]},
   {"MirrorRadius3",	PF_REAL_TYPE,	 	&Mirror_Radii[1]},
   {"MirrorRadius4",	PF_REAL_TYPE,	 	&Mirror_Radii[2]},
   {"MirrorRadius6",	PF_REAL_TYPE,	 	&Mirror_Radii[3]},
   {"MirrorEAFile",	PF_FILE_TYPE,		&HRMA_EA_File},
   {"MirrorUseEA",	PF_BOOLEAN_TYPE,	&Mirror_Use_Effective_Area},
   {"Shutters1", 	PF_STRING_TYPE,		&Mirror_Shutters[0]},
   {"Shutters3", 	PF_STRING_TYPE,		&Mirror_Shutters[1]},
   {"Shutters4", 	PF_STRING_TYPE,		&Mirror_Shutters[2]},
   {"Shutters6", 	PF_STRING_TYPE,		&Mirror_Shutters[3]},
   {NULL, 0, NULL}
};

/*}}}*/

typedef struct /*{{{*/
{
   unsigned int n_energies;
   float *energies;
   float *cum_eff_area[MARX_NUM_MIRRORS];  /* normalized by the total geometric area */
   double radius[MARX_NUM_MIRRORS];
   unsigned int n_mirrors;
   /* It is not quite clear that these parameters should be here or elsewhere. */
   double vignetting_factor;
   unsigned int shutter_bitmap[MARX_NUM_MIRRORS];
   unsigned int num_open_shutters[MARX_NUM_MIRRORS];
}

/*}}}*/
Mirror_Type;

static Mirror_Type Mirrors;

#ifdef NOT_IMPLEMENTED
static void create_on_axis_ray (Marx_Photon_Attr_Type *at) /*{{{*/
{
   double azimuth;

   /* Assign an azimuthal angle */
   azimuth = JDMrandom () * (2.0 * PI);
   at->x.x = 0.0;
   at->x.y = cos (azimuth);
   at->x.z = sin (azimuth);
   at->p.x = 1.0;
   at->p.y = 0.0;
   at->p.z = 0.0;
}

/*}}}*/
#endif

/* The mirror file contains 7 columns:
 *    wavelength, energy, total_area, area_1, area_3, area_4, area_6
 * The first line of this format contains geometric information (zero energy limit.)
 * The following lines are effective areas.
 */
int _marx_ea_mirror_init (Param_File_Type *p) /*{{{*/
{
   double g_area;
   unsigned int nread;
   unsigned int i, j;
   unsigned int n_mirrors = MARX_NUM_MIRRORS;
   char *file;

   if ((-1 == pf_get_parameters (p, HRMA_Parm_Table))
       || (HRMA_EA_File == NULL) || (0 == *HRMA_EA_File))
     return -1;

   for (j = 0; j < n_mirrors; j++)
     {
	if (-1 == _marx_parse_shutter_string (Mirror_Shutters[j],
					      &Mirrors.shutter_bitmap[j], &Mirrors.num_open_shutters[j]))
	  return -1;
     }

   if (-1 == _marx_init_mirror_blur (p))
     {
	JDMmsg_error ("Unable to initial mirror blur.");
	return -1;
     }

   if (NULL == (file = marx_make_data_file_name (HRMA_EA_File)))
     return -1;

   /* clean up from a previous call */
   for (j = 0; j < n_mirrors; j++)
     {
	if (Mirrors.cum_eff_area[j] != NULL)
	  JDMfree_float_vector (Mirrors.cum_eff_area[j]);
	Mirrors.cum_eff_area[j] = NULL;
     }
   if (Mirrors.energies != NULL) JDMfree_float_vector (Mirrors.energies);
   Mirrors.energies = NULL;

   marx_message ("Reading binary mirror effective area file\n\t%s\n", file);

   if (-1 == marx_f_read_bdat (file, &nread, 7,
			       NULL,   /* wavelength */
			       &Mirrors.energies,
			       NULL,   /* total area */
			       &Mirrors.cum_eff_area[0],
			       &Mirrors.cum_eff_area[1],
			       &Mirrors.cum_eff_area[2],
			       &Mirrors.cum_eff_area[3]))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   /* The zeroth element of the arrays is the geometric area.  Find the total
    * geometric area taking into account of shutters.
    */
   g_area = 0.0;
   for (j = 0; j < n_mirrors; j++)
     {
	if (Mirrors.num_open_shutters[j])
	  {
	     g_area += 0.25 * Mirrors.cum_eff_area [j][0] * Mirrors.num_open_shutters[j];
	  }
     }

   if (g_area <= 0.0)
     {
	marx_error ("The mirror geometric area is 0.  Check shutters.");
	return -1;
     }

   Marx_Mirror_Geometric_Area = g_area;

   /* Now compute the cumulative effective area.  Take shutters into account. */
   for (i = 0; i < nread; i++)
     {
	double total = 0.0;
	for (j = 0; j < n_mirrors; j++)
	  {
	     if (Mirrors.num_open_shutters[j])
	       {
		  total += 0.25 * Mirrors.cum_eff_area[j][i] * Mirrors.num_open_shutters[j];
	       }
	     Mirrors.cum_eff_area[j][i] = total / g_area;
	  }
     }

   Mirrors.n_mirrors = n_mirrors;
   Mirrors.n_energies = nread;
   Mirrors.vignetting_factor = Mirror_Vignetting_Factor;

   for (j = 0; j < n_mirrors; j++)
     {
	Mirrors.radius[j] = Mirror_Radii[j];
     }
   HRMA_Position = (Marx_Focal_Length - 830.0 / 2.0);   /* mm */

   return 0;
}

/*}}}*/

/* Now we give the photon a position and direction.  We must first project
 * photons to the mirror as on-axis rays since that is all this module
 * was designed for.
 */
static void reflect_photon (Marx_Photon_Attr_Type *at, unsigned int shell) /*{{{*/
{
   JDMVector_Type *x, *p;
   double radius, theta;
   unsigned int bitmap, quad;

   at->mirror_shell = shell;
   radius = Mirrors.radius[shell];
   bitmap = Mirrors.shutter_bitmap[shell];

   do
     {
	theta = JDMrandom ();
	quad = (unsigned int) (4.0 * theta);   /* 0, 1, 2, 3 */
     }
   while (0 == (bitmap & (1 << quad)));

   theta = (2.0 * PI) * (theta - 1.0 / 8.0);

   x = &(at->x);

   x->y = radius * sin (theta);
   x->z = radius * cos (theta);
   x->x = HRMA_Position;

   p = &(at->p);
   /* For a perfectly focused ray, the direction is just -x/|x|. */
   p->x = -x->x;
   p->y = -x->y;
   p->z = -x->z;
   JDMv_normalize (p);
}

/*}}}*/

static int mirror_perfect_reflect (Marx_Photon_Type *);

int _marx_ea_mirror_reflect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   double r;
   double *photon_energies;
   unsigned int n, i, *sorted_index;
   double **interp_areas;
   unsigned int n_mirrors = Mirrors.n_mirrors;

   if (pt->history & MARX_MIRROR_SHELL_OK)
     return 0;			       /* been here already */

   pt->history |= MARX_MIRROR_SHELL_OK;

   if (Mirror_Use_Effective_Area == 0)
     return mirror_perfect_reflect (pt);

   marx_prune_photons (pt);
   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   /* First of all, apply to vignetting factor to kill a certain percentage
    * of photons.
    */
   for (i = 0; i < n; i++)
     {
	if (JDMrandom () > Mirrors.vignetting_factor)
	  {
	     (photon_attributes + sorted_index[i])->flags |= PHOTON_MIRROR_VBLOCKED;
	  }
     }

   /* I could have pruned in the previous loop but it is a better idea to
    * leave it for a function call.
    */
   marx_prune_photons (pt);
   photon_energies = pt->sorted_energies;
   n = pt->num_sorted;

   interp_areas = JDMdouble_matrix (n_mirrors, n);
   if (interp_areas == NULL) return -1;

   /* Now perform the interpolation */
   JDMinterpolate_n_dfvector (photon_energies, interp_areas, n,
			      Mirrors.energies, Mirrors.cum_eff_area, Mirrors.n_energies,
			      n_mirrors);

   for (i = 0; i < n; i++)
     {
	unsigned int shell;

	r = JDMrandom ();
	at = photon_attributes + sorted_index[i];

	for (shell = 0; shell < n_mirrors; shell++)
	  {
	     if (r <= interp_areas[shell][i])
	       {
		  reflect_photon (at, shell);
		  break;
	       }
	  }
	if (shell == n_mirrors)
	  {
#if 0
	     static FILE *fp;
	     if (fp == NULL) fp = fopen ("out", "w");
	     if (fp != NULL)
	       {
		  fprintf (fp, "%d\t%16g%16g\n", i, at->arrival_time, at->energy);
	       }
#endif
	     at->flags |= PHOTON_UNREFLECTED;
	  }
     }

   JDMfree_double_matrix (interp_areas, n_mirrors);
   return _marx_mirror_blur (pt);
}

/*}}}*/

static int mirror_perfect_reflect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   double r;
   double *photon_energies;
   unsigned int n, i, *sorted_index;
   unsigned int n_mirrors = Mirrors.n_mirrors, shell;
   double cum[MARX_NUM_MIRRORS];

   marx_prune_photons (pt);
   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   /* First of all, apply to vignetting factor to kill a certain percentage
    * of photons.
    */
   for (i = 0; i < n; i++)
     {
	if (JDMrandom () > Mirrors.vignetting_factor)
	  {
	     (photon_attributes + sorted_index[i])->flags |= PHOTON_MIRROR_VBLOCKED;
	  }
     }

   marx_prune_photons (pt);
   photon_energies = pt->sorted_energies;
   n = pt->num_sorted;

   for (shell = 0; shell < n_mirrors; shell++)
     {
	cum[shell] = Mirrors.cum_eff_area[shell][0];
     }

   for (i = 0; i < n; i++)
     {
	r = JDMrandom ();
	at = photon_attributes + sorted_index[i];

	for (shell = 0; shell < n_mirrors; shell++)
	  {
	     if (r <= cum[shell])
	       {
		  reflect_photon (at, shell);
		  break;
	       }
	  }
	if (shell == n_mirrors)
	  at->flags |= PHOTON_UNREFLECTED;
     }

   return _marx_mirror_blur (pt);
}

/*}}}*/

#ifdef NOT_IMPLEMENTED
int init_mirror_module (void) /*{{{*/
{
   return 0;
}

/*}}}*/
int dump_mirror_curves (double energy, FILE *fp) /*{{{*/
{
   int i;
   double **eff_area;
   double total;

   for (i = 0; i < Mirrors.n_mirrors; i++)
     {
	if (Mirrors.cum_eff_area[i] == NULL)
	  {
	     marx_error ("dump_mirror_curves: no data.");
	     return -1;
	  }
     }

   eff_area = JDMdouble_matrix (Mirrors.n_mirrors, 1);
   if (eff_area == NULL)
     {
	marx_error ("dump_mirror_curves: Not enough memory.");
	return -1;
     }

   JDMinterpolate_n_dvector (&energy, eff_area, 1,
			     Mirrors.energies, Mirrors.cum_eff_area,
			     Mirrors.n_energies,
			     Mirrors.n_mirrors);

   total = 0.0;
   for (i = 0; i < Mirrors.n_mirrors; i++)
     {
	double value;

	value = eff_area[i][0] - total;
	total += value;

	if (EOF == fprintf (fp, "\t%e", value))
	  {
	     marx_error ("dump_mirror_curves: write error.");
	     JDMfree_double_matrix (eff_area, Mirrors.n_mirrors);
	     return -1;
	  }
     }
#if 1
   fprintf (fp, "\t%e", total);
#endif
   JDMfree_double_matrix (eff_area, Mirrors.n_mirrors);
   return 0;
}

/*}}}*/
unsigned int dump_mirror_curves_open (unsigned int i, FILE *fp) /*{{{*/
{
   unsigned int j;

   fputs ("\n#Mirror:", fp);
   for (j = 0; j < 4; j++)
     {
	fprintf (fp, " Col_%u:Mirror-Shell-%u", i, j);
	i++;
     }

   fprintf (fp, " Col_%u:Total_Mirror", i);
   i++;
   return i;
}

/*}}}*/
int dump_mirror_variables (FILE *fp) /*{{{*/
{
   fprintf (fp, "--Mirror Module Variables--\n");
   fprintf (fp, " Mirror_Geometric_Area: %f\n", Mirror_Geometric_Area);
   fprintf (fp, " Mirror_Vignetting_Factor: %f\n", Mirror_Vignetting_Factor);
   return 0;
}

/*}}}*/
#endif
