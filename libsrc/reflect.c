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
#include <jdmath.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "marx.h"
#include "_marx.h"

/* This routine assumes that cos_theta is positive */
double marx_reflectivity (double cos_theta, double beta, double delta) /*{{{*/
{
   JDMComplex_Type n, root, nsqr, e_perp, e_par, num, den;
   double sin_theta;
   /* See section 7.3 of Jackson (1975).
    * We assume that n = 1 and mu = 1.  Below. Jackson's n' is represented by
    * just n.
    */

   /* Note:  If delta == 0, and beta == 0, the reflectivity will be zero.
    * However, if delta == 1, and beta == 0, the reflectivity will be unity
    */

   n.r = (1.0 - delta);
   n.i = beta;

   sin_theta = sqrt (1.0 - cos_theta * cos_theta);

   nsqr = JDMc_mul (n, n);

   root = JDMc_sqrt (JDMc_a_bz (-sin_theta*sin_theta, 1.0, nsqr));

   /* eq 7.39 */
   num = JDMc_a_bz (cos_theta, -1.0, root);
   den = JDMc_a_bz (cos_theta,  1.0, root);

   e_perp = JDMc_div (num, den);

   /* eq 7.41 */
   num = JDMc_az1_bz2 (cos_theta, nsqr, -1.0, root);
   den = JDMc_az1_bz2 (cos_theta, nsqr,  1.0, root);

   e_par = JDMc_div (num, den);

   /* Now average over polarizations */

   return 0.5 * (e_par.r * e_par.r + e_par.i * e_par.i
		 + e_perp.r * e_perp.r + e_perp.i * e_perp.i);
}

/*}}}*/

double marx_interp_reflectivity (double energy, double cos_theta, /*{{{*/
				 float *energies, float *betas, float *deltas,
				 unsigned int num_energies)
{
   float beta, delta;

   beta = JDMinterpolate_f (energy, energies, betas, num_energies);
   delta= JDMinterpolate_f (energy, energies, deltas, num_energies);

   return marx_reflectivity (fabs (cos_theta), beta, delta);
}

/*}}}*/

