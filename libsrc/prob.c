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
/* These routines are used to randomally pick a number from a probabliity
 * distribution.
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

/* This routine returns m random values from a cumulative spectrum ycum of n
 * values.  Here ycum is assumed to be normalized such that it is a
 * monotonically increasing function with values that lie from 0 to 1.
 */
void marx_get_random_event (double *xevents, double *ycum, unsigned int n, /*{{{*/
			    double *xrnd, unsigned int m)
{
   double *xrnd_max = xrnd + m;
   double r;

   while (xrnd < xrnd_max)
     {
	r = JDMrandom ();
	*xrnd = JDMinterpolate_d (r, ycum, xevents, n);
	xrnd++;
     }
}

/*}}}*/

