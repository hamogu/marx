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

#define NUM_FILTER_REGIONS	4

static _Marx_HRC_QE_Type MCP_QEs [_MARX_NUM_HRC_S_CHIPS];
static _Marx_HRC_QE_Type Filter_QEs [NUM_FILTER_REGIONS];

/* The MCP ids for the HRC-S are numbered 1,2,3.  We want to map these two
 * 2,1,0 to correspond to the regions.  Using a lookup table is more general.
 */
static int Mcp_Id_Mapping[4] =
{
   -1, 2, 1, 0
};

#if MARX_HAS_DRAKE_FLAT
static int Use_Drake_Flat;
#endif

JDMVector_Type _Marx_HRC_Geometric_Center;

/*---------------------------------------------------------------------------
 *
 * HRC-S CCD Geometry parameters
 *
 *---------------------------------------------------------------------------*/
static Marx_HRC_Blur_Parm_Type *HRC_S_Blur_Parms;
static Marx_Detector_Geometry_Type *HRC_S_MCPs;

static double Shield_OffsetT;   /* mm */
static double Shield_OffsetL;   /* mm */
static double Shield_OffsetR;   /* mm */
static double Shield_OffsetX;   /* mm */
static double Shield_OffsetSR;
static double Shield_OffsetSL;
static double Shield_GapL;
static double Shield_GapR;
static double Shield_OffsetSL_Gap;
static double Shield_OffsetSR_Gap;

/* Location of detector center in marx system */
static double Shield_Y_Center;
static double Shield_Z_Center;

static Param_Table_Type HRC_Parm_Table [] = /*{{{*/
{
   {"HRC-S-QEFile0",	PF_FILE_TYPE,		&MCP_QEs[0].file},
   {"HRC-S-QEFile1",	PF_FILE_TYPE,		&MCP_QEs[1].file},
   {"HRC-S-QEFile2",	PF_FILE_TYPE,		&MCP_QEs[2].file},
   {"HRC-S-UVISFile0",PF_FILE_TYPE,		&Filter_QEs[0]},
   {"HRC-S-UVISFile1",PF_FILE_TYPE,		&Filter_QEs[1]},
   {"HRC-S-UVISFile2",PF_FILE_TYPE,		&Filter_QEs[2]},
   {"HRC-S-UVISFile3",PF_FILE_TYPE,		&Filter_QEs[3]},
#if MARX_HAS_DRAKE_FLAT
   {"HRC-HESF",		PF_BOOLEAN_TYPE,	&Use_Drake_Flat},
#endif
   {NULL, 0, NULL}
};

/*}}}*/

static int read_hrc_s_blur_parms (Param_File_Type *p)
{
   if (NULL == (HRC_S_Blur_Parms = _marx_hrc_blur_open (p, 'S')))
     return -1;

   return 0;
}


/* This routine is based on a memo dated March 14 1996 from Mike Juda with
 * title: New HRC UV/Ion Shields.   This memo shows the following figure:
 * @
 * @   +-----------------+ +-------------------+ +-------------------+
 * @   |      2(5)       | |         0(1)      | |      2(6)         |
 * @   +-----------------+ +-----+        +----+ +-------------------+  ^
 * @   |                 | |     |        |    | |                   |  +- T
 * @   |        +        | |     |   +    |    | |         +         |  v
 * @   |                 | |     |        |    | |                   |
 * @   |     3           | | 1(2)|        | 1(2) |       3(4)        |
 * @   |                 | |     |        |    | |                   |
 * @   +-----------------+ +-----+--------+----+ +-------------------+
 * @                             <-L-><-R->
 * @                       <---SL----><---SR-->
 * @                   --> <-- GapL         -->  <-- GapR
 *
 * The labeling of the regions 0-4 is mine.  The values in parenthesis are
 * what are using on the HRC web pages (mine predate theirs by about 5 years).
 * The distances L, R, S, and T
 * correspond to distances from the geometric center of the middle
 * facet.
 *
 * Apparantly, this filter lies 10mm along the positive x axis from the center
 * MCP.  It lies in a flat plane. (Shield_OffsetX)
 */
static int
get_filter_region (double y, double z, unsigned char *r) /*{{{*/
{
   /* y and z are specified in the Marx Coordinate system.  The filter is
    * offset from this.
    */
   y -= Shield_Y_Center;
   z -= Shield_Z_Center;

   if (y < 0)
     {
	y = -y;
	if (y < Shield_OffsetL)
	  *r = 0;
	else if (y < Shield_OffsetSL)
	  {
	     if (z >= Shield_OffsetT)
	       *r = 0;
	     else
	       *r = 1;
	  }
	else if (y >= Shield_OffsetSL_Gap)
	  {
	     if (z >= Shield_OffsetT)
	       *r = 2;
	     else
	       *r = 3;
	  }
	else return -1;
     }
   else
     {
	if (y < Shield_OffsetR)
	  *r = 0;
	else if (y < Shield_OffsetSR)
	  {
	     if (z >= Shield_OffsetT)
	       *r = 0;
	     else
	       *r = 1;
	  }
	else if (y >= Shield_OffsetSR_Gap)
	  {
	     if (z >= Shield_OffsetT)
	       *r = 2;
	     else
	       *r = 3;
	  }
	else return -1;
     }

   return 0;
}

/*}}}*/

static int
apply_hrc_qe (Marx_Photon_Attr_Type *at, unsigned int mcp_id)
{
   _Marx_HRC_QE_Type *mcp;
   double t, y, z;
   unsigned char region;
   double qe;
   double energy;
   JDMVector_Type *x, *p;

   energy = at->energy;
   x = &at->x;
   p = &at->p;

   mcp = MCP_QEs + Mcp_Id_Mapping[mcp_id];
   if (mcp->num_energies != 0)
     {
	qe = JDMinterpolate_f ((float) energy,
			       mcp->energies, mcp->eff, mcp->num_energies);
	if (JDMrandom () >= qe)
	  return -1;
     }

   /* Now project it back to see where it hit the filter */
   t = (Shield_OffsetX - x->x) / p->x;
   y = x->y + t * p->y;
   z = x->z + t * p->z;

   if (-1 == get_filter_region (y, z, &region))
     return -1;

   mcp = Filter_QEs + (unsigned int) region;
   if (mcp->num_energies != 0)
     {
	qe = JDMinterpolate_f ((float) energy,
			       mcp->energies, mcp->eff, mcp->num_energies);
	if (JDMrandom () >= qe)
	  return -1;
     }
   at->detector_region = region;
   at->pulse_height = _marx_hrc_compute_pha (energy);
   return 0;
}

int _marx_hrc_s_detect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *at, *attrs;
   unsigned int n_photons, i;
   unsigned int *sorted_index;

   if (pt->history & MARX_DET_NUM_OK)
     return 0;

   pt->history |= (MARX_DET_REGION_OK | MARX_PULSEHEIGHT_OK
		   | MARX_DET_PIXEL_OK | MARX_DET_NUM_OK);

   pt->history |= (MARX_DET_UV_PIXEL_OK);

#if MARX_HAS_DRAKE_FLAT
   if (Use_Drake_Flat && (-1 == _marx_drake_reflect (pt)))
     return -1;
#endif
   marx_prune_photons (pt);

   attrs = pt->attributes;
   n_photons = pt->num_sorted;
   sorted_index = pt->sorted_index;

   for (i = 0; i < n_photons; i++)
     {
	Marx_Detector_Geometry_Type *d;
	double dx, dy, u, v;

	at = attrs + sorted_index[i];

#if MARX_HAS_DITHER
	_marx_dither_detector (&at->dither_state);
#endif
	/* Transform ray into local system */
	_marx_transform_ray (&at->x, &at->p,
			    &_Marx_Det_XForm_Matrix);

	d = _marx_intersect_with_detector (at->x, at->p,
					   HRC_S_MCPs,
					   &at->x, &dx, &dy,
					  _Marx_Det_Extend_Flag);

	if (d == NULL)
	  {
	     at->flags |= PHOTON_MISSED_DETECTOR;
	     at->ccd_num = -1;
	  }
	else if (-1 == apply_hrc_qe (at, d->id))
	  {
	     at->flags |= PHOTON_UNDETECTED;
	     at->ccd_num = -1;
	  }
	else
	  {
	     _marx_hrc_blur_position (HRC_S_Blur_Parms, &dx, &dy);

	     at->ccd_num = d->id;
	     (void) _marx_hrc_s_compute_pixel (at->ccd_num, dx, dy, &dx, &dy,
					       &u, &v);

	     at->y_pixel = dx;
	     at->z_pixel = dy;
	     at->u_pixel = u;
	     at->v_pixel = v;
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

int _marx_hrc_read_efficiencies (_Marx_HRC_QE_Type *dt)
{
   char *file;
#if 0
   float *en, *en_max;
#endif
   if (dt == NULL)
     return 0;

   if (dt->energies != NULL) JDMfree_float_vector (dt->energies);
   dt->energies = NULL;

   if (dt->eff != NULL) JDMfree_float_vector (dt->eff);
   dt->eff = NULL;
   dt->num_energies = 0;

   file = dt->file;
   if (_Marx_Det_Ideal_Flag || (*file == 0))
     return 0;     /* perfect effic */

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;

   marx_message ("\t%s\n", file);

   if (-1 == marx_f_read_bdat (file, &dt->num_energies, 2, &dt->energies, &dt->eff))
     {
	marx_free (file);
	return -1;
     }
#if 0
   /* The energies are in eV, so scale them to KeV */
   en = dt->energies;
   en_max = en + dt->num_energies;
   while (en < en_max)
     {
	*en = *en * 0.001;
	en++;
     }
#endif
   marx_free (file);
   return 0;
}

static _Marx_Simple_Data_Type IonShield_Data_Table [] =
{
   {"Ion_Shield_T",	1,	&Shield_OffsetT,	1.0, 0},
   {"Ion_Shield_L",	1,	&Shield_OffsetL,	1.0, 0},
   {"Ion_Shield_R",	1,	&Shield_OffsetR,	1.0, 0},
   {"Ion_Shield_X",	1,	&Shield_OffsetX,	1.0, 0},
   {"Ion_Shield_SL",	1,	&Shield_OffsetSL,	1.0, 0},
   {"Ion_Shield_SR",	1,	&Shield_OffsetSR,	1.0, 0},
   {"Ion_Shield_GapL",	1,	&Shield_GapL,		1.0, 0},
   {"Ion_Shield_GapR",	1,	&Shield_GapR,		1.0, 0},
   {NULL, 0, NULL, 0.0, 0}
};

int _marx_hrc_s_init (Param_File_Type *p) /*{{{*/
{
   unsigned int i;
   Marx_Detector_Type *hrc;
   Marx_Detector_Geometry_Type *middle_mcp;
   char *file;

   if (-1 == pf_get_parameters (p, HRC_Parm_Table))
     return -1;

   if (-1 == _marx_hrc_s_geom_init (p))
     return -1;

   file = "hrc/hrc_s_geom.txt";

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;

   marx_message ("\t%s\n", file);

   if (-1 == _marx_read_simple_data_file (file, IonShield_Data_Table))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   /* The blur init must come after the geometry */
   if (-1 == read_hrc_s_blur_parms (p))
     return -1;

   if (NULL == (hrc = marx_get_detector_info ("HRC-S")))
     return -1;

   HRC_S_MCPs= hrc->facet_list;
   if (NULL == (middle_mcp = _marx_find_detector_facet (hrc, 2)))
     {
	marx_error ("Internal error: Unable to find HRC-S with id=2");
	return -1;
     }

   while ((middle_mcp != NULL) && (middle_mcp->id != 2))
     middle_mcp = middle_mcp->next;

   marx_message ("Reading binary HRC-S QE/UVIS data files:\n");

   for (i = 0; i < _MARX_NUM_HRC_S_CHIPS; i++)
     {
	if (-1 == _marx_hrc_read_efficiencies (MCP_QEs + i))
	  return -1;
     }

   for (i = 0; i < NUM_FILTER_REGIONS; i++)
     {
	if (-1 == _marx_hrc_read_efficiencies (Filter_QEs + i))
	  return -1;
     }

   /* Now compute geometric center position */
#if 0
   Shield_OffsetS = (middle_mcp->x_lr.y
		     - middle_mcp->x_ll.y) / 2.0;
#endif

   _Marx_HRC_Geometric_Center
     = JDMv_sum (JDMv_vector (_Marx_Det_XForm_Matrix.dx,
			      _Marx_Det_XForm_Matrix.dy,
			      _Marx_Det_XForm_Matrix.dz),
		 JDMv_ax1_bx2 (0.5, middle_mcp->x_ll,
			       0.5, middle_mcp->x_ur));

   Shield_Y_Center = _Marx_HRC_Geometric_Center.y;
   Shield_Z_Center = _Marx_HRC_Geometric_Center.z;
   Shield_OffsetSR_Gap = Shield_OffsetSR + Shield_GapR;
   Shield_OffsetSL_Gap = Shield_OffsetSL + Shield_GapL;

#if MARX_HAS_DRAKE_FLAT
   if (Use_Drake_Flat)
     {
	marx_message ("Initializing Drake flat...\n");
	if (-1 == _marx_drake_flat_init (p))
	  return -1;
     }
#endif
   return 0;
}

/*}}}*/

