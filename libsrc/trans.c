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

#include "marx.h"
#include "_marx.h"

static void rotate_vector (double *m, JDMVector_Type *v)
{
   double x, y, z;

   x = v->x;
   y = v->y;
   z = v->z;

   v->x = m[0] * x + m[1] * y + m[2] * z;
   v->y = m[3] * x + m[4] * y + m[5] * z;
   v->z = m[6] * x + m[7] * y + m[8] * z;
}

static void rotate_vector_inv (double *m, JDMVector_Type *v)
{
   double x, y, z;

   x = v->x;
   y = v->y;
   z = v->z;

   v->x = m[0] * x + m[3] * y + m[6] * z;
   v->y = m[1] * x + m[4] * y + m[7] * z;
   v->z = m[2] * x + m[5] * y + m[8] * z;
}

void _marx_transform_ray (JDMVector_Type *x, JDMVector_Type *p, /*{{{*/
			  _Marx_Coord_Transform_Type *a)
{
   x->x -= a->dx;
   x->y -= a->dy;
   x->z -= a->dz;

   rotate_vector (a->matrix, x);
   rotate_vector (a->matrix, p);
}

/*}}}*/

void _marx_transform_ray_reverse (JDMVector_Type *x, JDMVector_Type *p, /*{{{*/
				  _Marx_Coord_Transform_Type *a)
{
   rotate_vector_inv (a->matrix, p);
   rotate_vector_inv (a->matrix, x);

   x->x += a->dx;
   x->y += a->dy;
   x->z += a->dz;
}

/*}}}*/
