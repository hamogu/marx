/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2002-2004 Massachusetts Institute of Technology

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
#include <string.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <jdmath.h>

#include "marx.h"
#include "_marx.h"

_Marx_Coord_Transform_Type _Marx_Det_XForm_Matrix;
int _Marx_Det_Ideal_Flag;

static Param_Table_Type Det_Parm_Table [] =
{
   {"DetOffsetX",	PF_REAL_TYPE,		&_Marx_Det_XForm_Matrix.dx},
   {"DetOffsetY",	PF_REAL_TYPE,		&_Marx_Det_XForm_Matrix.dy},
   {"DetOffsetZ",	PF_REAL_TYPE,		&_Marx_Det_XForm_Matrix.dz},
   {"DetIdeal",		PF_BOOLEAN_TYPE,	&_Marx_Det_Ideal_Flag},
   {NULL, 0, NULL}
};


Marx_Detector_Geometry_Type *
_marx_intersect_with_detector (JDMVector_Type x0, JDMVector_Type p, /*{{{*/
			       Marx_Detector_Geometry_Type *d, unsigned int n,
			       JDMVector_Type *x, double *dx, double *dy)
{
   JDMVector_Type normal, r;
   double pdotn;
   double rx, ry;
   Marx_Detector_Geometry_Type *dmax;
   
   for (dmax = d + n; d < dmax; d++)
     {
	normal = d->normal;
	pdotn = JDMv_pdot_prod (&p, &normal);

	if (pdotn == 0)
	  continue;
	
	r = JDMv_diff (x0, d->x_ll);
	
	r = JDMv_pax1_bx2 (1.0, &r, 
			   -1.0 * JDMv_pdot_prod (&r, &normal)/pdotn, &p);
	/* This r is defined in the plane of the detector with origin 
	 * at d->x_ll.
	 */
	
	/* Now check to see if the coordinates in the plane of the detector
	 * chip is within the boundaries of the chip.  Here we assume the 
	 * chip is rectangular.
	 */
	
	rx = JDMv_dot_prod (r, d->xhat);
	if ((rx < 0.0) || (rx >= d->xlen))
	  continue;
	
	ry = JDMv_dot_prod (r, d->yhat);
	if ((ry < 0.0) || (ry >= d->ylen))
	  continue;
	
	*x = JDMv_sum (r, d->x_ll);
	*dx = rx;
	*dy = ry;
	
	return d;
     }

   /* Missed */
   return NULL;
}

/*}}}*/

typedef struct
{
   char *name;
   int (*init) (Param_File_Type *);
   int (*detect) (Marx_Photon_Type *);
   int type;
}
Detector_Cap_Type;

static Detector_Cap_Type *Detector;

static int plane_init (Param_File_Type *p)
{
   (void) p;
   return 0;
}

static int plane_detect (Marx_Photon_Type *pt)
{
   Marx_Photon_Attr_Type *at, *attrs;
   unsigned int n_photons, i;
   unsigned int *sorted_index;

   if (pt->history & MARX_CCD_NUM_OK)
     return 0;
   
   pt->history |= (MARX_HRC_REGION_OK | MARX_PULSEHEIGHT_OK 
		   | MARX_Y_PIXEL_OK | MARX_Z_PIXEL_OK | MARX_CCD_NUM_OK);
   
   marx_prune_photons (pt);

   attrs = pt->attributes;
   n_photons = pt->num_sorted;
   sorted_index = pt->sorted_index;
   
   for (i = 0; i < n_photons; i++)
     {
	double t, pixel_size = 0.00643;/* 6.43 um */
	
	at = attrs + sorted_index[i];
	
	_marx_transform_ray (&at->x, &at->p, &_Marx_Det_XForm_Matrix);
	/* Solve: 0 = X_x + P_x t
	 * and project:
	 *       X_y = X_y + t P_y, etc... 
	 */
	t = -at->x.x / at->p.x;
	at->x.x = 0;
	at->x.y += t * at->p.y;
	at->x.z += t * at->p.z;

	at->ccd_num = 0;
	at->y_pixel = at->x.y / pixel_size;
	at->z_pixel = at->x.z / pixel_size;
	at->detector_region = 0;
	
	/* Transform detected photon back to original system */
	_marx_transform_ray_reverse (&at->x, &at->p, &_Marx_Det_XForm_Matrix);
     }
   return 0;
}

static Detector_Cap_Type Detector_Caps [] =
{
   {"ACIS-S", _marx_acis_s_init, _marx_acis_s_detect, MARX_DETECTOR_ACIS_S},
   {"ACIS-I", _marx_acis_i_init, _marx_acis_i_detect, MARX_DETECTOR_ACIS_I},
   {"HRC-S", _marx_hrc_s_init, _marx_hrc_s_detect, MARX_DETECTOR_HRC_S},
   {"HRC-I", _marx_hrc_i_init, _marx_hrc_i_detect, MARX_DETECTOR_HRC_I},
   {"PLANE", plane_init, plane_detect, MARX_DETECTOR_PLANE},
   {"NONE", NULL, NULL, 0},
   {NULL, NULL, NULL}
};

char *Marx_Supported_Detectors = "ACIS-S, ACIS-I, HRC-S, HRC-I";

static void rotate_matrix (double *matrix, double theta)
{
   double m00, m01, m10, m11;
   double c, s;

   if (theta == 0)
     {
	s = 0; 
	c = 1;
     }
   else 
     {
	s = sin (theta);
	c = cos (theta);
     }
   
   m00 = matrix[4]; m01=matrix[5];
   m10 = matrix[7]; m11=matrix[8];
   matrix [4] = m00*c - m01*s; matrix[5] = m00*s + m01*c;
   matrix [7] = m10*c - m11*s; matrix[8] = m10*s + m11*c;
}

int _marx_set_detector_angle (double theta)
{
   rotate_matrix (_Marx_Det_XForm_Matrix.matrix, theta);
   return 0;
}

      

#if MARX_HAS_DITHER
void _marx_dither_detector (Marx_Dither_Type *d)
{
   if (_Marx_Dither_Mode == 0)
     return;

   _Marx_Det_XForm_Matrix.dy += d->dy;
   _Marx_Det_XForm_Matrix.dz += d->dz;

   rotate_matrix (_Marx_Det_XForm_Matrix.matrix, d->dtheta);
}

void _marx_undither_detector (Marx_Dither_Type *d)
{
   if (_Marx_Dither_Mode == 0)
     return;

   _Marx_Det_XForm_Matrix.dy -= d->dy;
   _Marx_Det_XForm_Matrix.dz -= d->dz;
   
   rotate_matrix (_Marx_Det_XForm_Matrix.matrix, -d->dtheta);
}
#endif				       /* MARX_HAS_DITHER */

int marx_detector_init (Param_File_Type *pf) /*{{{*/
{
   char buf[PF_MAX_LINE_LEN];
   Detector_Cap_Type *d;
   double *matrix;

   /* Possible values: ACIS-S,HRC-S,NONE */
   
   Detector = NULL;

   if (-1 == pf_get_string (pf, "DetectorType", buf, sizeof (buf)))
     return -1;
   
   d = Detector_Caps;
   
   while (d->name != NULL)
     {
	if (!strcmp (d->name, buf))
	  {
	     Detector = d;
	     break;
	  }
	d++;
     }
   
   if (Detector == NULL)
     {
	marx_error ("DetectorType \"%s\" not supported.", buf);
	return -1;
     }
   
   if (Detector->detect == NULL)
     return 0;			       /* NONE */

   marx_message ("Initializing %s detector...\n", Detector->name);

   matrix = _Marx_Det_XForm_Matrix.matrix;

   /* This matrix is used to rotate the detector during aspect motion.  
    * Actually, as far as the detector goes, the AXAF aspect solution 
    * provides information about the offset of the detector origin in the 
    * (y,z) plane and a rotation about the X axis.  These should not be 
    * confused with the dither of the optical axis.  The detector offsets
    * are a result of thermal expansion, etc...
    */
   matrix[0] = 1; matrix [1] = 0; matrix [2] = 0;
   matrix[3] = 0; matrix [4] = 1; matrix [5] = 0;
   matrix[6] = 0; matrix [7] = 0; matrix [8] = 1;

   if (-1 == pf_get_parameters (pf, Det_Parm_Table))
     return -1;

   if (-1 == (*Detector->init) (pf))
     {
	marx_error ("Error initializing DetectorType %s.", buf);
	return -1;
     }
   
   return Detector->type;
}

/*}}}*/

int marx_detect (Marx_Photon_Type *pt, int verbose) /*{{{*/
{
   if (Detector == NULL)
     {
	marx_error ("Detector not initialized.");
	return -1;
     }
   
   if (Detector->detect != NULL)
     {
	if (verbose) 
	  {
	     marx_message ("Detecting with %s\n", Detector->name);
	  }
	return (*Detector->detect) (pt);
     }
   
   return 0;
}

/*}}}*/
