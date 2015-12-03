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

/*{{{ #includes */

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

static double FF_MinZ;
static double FF_MinY;
static double FF_MaxY;
static double FF_MaxZ;
static double FF_XPos;

static Param_Table_Type FlatField_Parm_Table [] =
{
   {"FF_MinY",	PF_REAL_TYPE,	&FF_MinY},
   {"FF_MinZ",	PF_REAL_TYPE,	&FF_MinZ},
   {"FF_MaxY",	PF_REAL_TYPE,	&FF_MaxY},
   {"FF_MaxZ",	PF_REAL_TYPE,	&FF_MaxZ},
   {"FF_XPos",	PF_REAL_TYPE,	&FF_XPos},

   {NULL, 0, NULL}
};

static void project_photon (Marx_Photon_Attr_Type *at, /*{{{*/
			    double source_distance)
{
   at->x.z = FF_MinZ + JDMrandom () * (FF_MaxZ - FF_MinZ);
   at->x.y = FF_MinY + JDMrandom () * (FF_MaxY - FF_MinY);
   at->x.x = FF_XPos;

   /* Source at infinity has all rays directed to
    * origin.  In this case, the photon generating
    * function sets the position to the unit vector
    * pointing toward the origin.  In other words, at->p
    * will stay the same.  So, consider only change for finite
    * source.
    */
   if (source_distance > 0.0)
     {
	at->p = JDMv_ax1_bx2 (1.0, at->x, source_distance, at->p);
	JDMv_normalize (&at->p);
     }
}

/*}}}*/

int _marx_ff_mirror_reflect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int n, i, *sorted_index;
   double source_distance;

   marx_prune_photons (pt);
   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   /* I could have pruned in the previous loop but it is a better idea to
    * leave it for a function call.
    */
   /* source_distance is in mm */
   source_distance = pt->source_distance;

   for (i = 0; i < n; i++)
     {
	at = photon_attributes + sorted_index[i];
	project_photon (at, source_distance);
     }

   pt->history |= MARX_MIRROR_SHELL_OK;

   return 0;
}

/*}}}*/

int _marx_ff_mirror_init (Param_File_Type *p) /*{{{*/
{
   if (-1 == pf_get_parameters (p, FlatField_Parm_Table))
     return -1;

   if ((FF_MaxY <= FF_MinY)
       || (FF_MaxZ <= FF_MinZ)
       || (FF_XPos < 0.0))
     {
	marx_error ("FF_* parameters not physical");
	return -1;
     }
   Marx_Mirror_Geometric_Area
     = 0.01 * (FF_MaxY - FF_MinY) * (FF_MaxZ - FF_MinZ);   /* cm^2 */
   return 0;
}

/*}}}*/
