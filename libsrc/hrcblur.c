/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2012 Massachusetts Institute of Technology

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

/* The HRC blur is parametrized as combination of two gaussians plus a lorentzian.
 * Since the 2d lorentzian diverges, it has to be cutoff at some
 * radius, Rmax.  Its normalized form is given by:
 *
 *    L(r) = L0/(1+r^2/g^2)
 *
 * where 1/L0 = g^2 * PI * log(1 + Rmax^2/g^2).
 * Then 1 = \int_0^Rmax dr (2*PI*r)*L(r)
 *
 * The 2-d gaussian has the normalized form
 *
 *    G(r) = exp(-r^2/2s^2) / (2*PI*s^2);    1 = \int dr (2*PI*r) * G(r)
 */

struct _Marx_HRC_Blur_Parm_Type
{
   double gauss1_sigma;
   double gauss1_xctr;
   double gauss1_yctr;
   double gauss1_wgt;

   double gauss2_sigma;
   double gauss2_xctr;
   double gauss2_yctr;
   double gauss2_wgt;

   double lorentz1_hwhm;
   double lorentz1_xctr;
   double lorentz1_yctr;
   double lorentz1_rmax;
   double lorentz1_wgt;
};

static Marx_HRC_Blur_Parm_Type HRC_I_Blur_Parms;
static Marx_HRC_Blur_Parm_Type HRC_S_Blur_Parms;

static Param_Table_Type HRC_I_Blur_Parm_Table [] =
{
   {"HRC-I-BlurG1FWHM",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss1_sigma},
   {"HRC-I-BlurG1XCTR",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss1_xctr},
   {"HRC-I-BlurG1YCTR",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss1_yctr},
   {"HRC-I-BlurG1AMP",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss1_wgt},

   {"HRC-I-BlurG2FWHM",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss2_sigma},
   {"HRC-I-BlurG2XCTR",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss2_xctr},
   {"HRC-I-BlurG2YCTR",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss2_yctr},
   {"HRC-I-BlurG2AMP",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.gauss2_wgt},

   {"HRC-I-BlurL1FWHM",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.lorentz1_hwhm},
   {"HRC-I-BlurL1XCTR",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.lorentz1_xctr},
   {"HRC-I-BlurL1YCTR",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.lorentz1_yctr},
   {"HRC-I-BlurL1AMP",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.lorentz1_wgt},
   {"HRC-I-BlurL1RMAX",	PF_REAL_TYPE,	&HRC_I_Blur_Parms.lorentz1_rmax},
   {NULL, 0, NULL}
};

static Param_Table_Type HRC_S_Blur_Parm_Table [] =
{
   {"HRC-S-BlurG1FWHM",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss1_sigma},
   {"HRC-S-BlurG1XCTR",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss1_xctr},
   {"HRC-S-BlurG1YCTR",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss1_yctr},
   {"HRC-S-BlurG1AMP",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss1_wgt},

   {"HRC-S-BlurG2FWHM",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss2_sigma},
   {"HRC-S-BlurG2XCTR",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss2_xctr},
   {"HRC-S-BlurG2YCTR",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss2_yctr},
   {"HRC-S-BlurG2AMP",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.gauss2_wgt},

   {"HRC-S-BlurL1FWHM",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.lorentz1_hwhm},
   {"HRC-S-BlurL1XCTR",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.lorentz1_xctr},
   {"HRC-S-BlurL1YCTR",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.lorentz1_yctr},
   {"HRC-S-BlurL1AMP",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.lorentz1_wgt},
   {"HRC-S-BlurL1RMAX",	PF_REAL_TYPE,	&HRC_S_Blur_Parms.lorentz1_rmax},
   {NULL, 0, NULL}
};

static int fixup_blur_parms (Marx_HRC_Blur_Parm_Type *bt, double pixelsize)
{
   double fwhm_to_sigma;
   double g2, norm;

   fwhm_to_sigma = 0.42466090014400953; /* 1/sqrt(8*log(2)) */

   if ((bt->gauss1_wgt < 0)
       || (bt->gauss2_wgt < 0)
       || (bt->lorentz1_wgt < 0)
       || (bt->gauss1_sigma <= 0)
       || (bt->gauss2_sigma <= 0)
       || (bt->lorentz1_hwhm <= 0)
       || (bt->lorentz1_rmax <= 0))
     {
	marx_error ("%s", "One or more of the HRC blur parameters has an invalid/unsupported value");
	return -1;
     }

   /* Convert units from pixels to mm */
   bt->gauss1_sigma *= pixelsize;
   bt->gauss1_xctr *= pixelsize;
   bt->gauss1_yctr *= pixelsize;
   bt->gauss2_sigma *= pixelsize;
   bt->gauss2_xctr *= pixelsize;
   bt->gauss2_yctr *= pixelsize;
   bt->lorentz1_hwhm *= pixelsize;
   bt->lorentz1_xctr *= pixelsize;
   bt->lorentz1_yctr *= pixelsize;
   bt->lorentz1_rmax *= pixelsize;

   /* Shift to origin of first gaussian */
   bt->gauss2_xctr -= bt->gauss1_xctr;
   bt->gauss2_yctr -= bt->gauss1_yctr;
   bt->lorentz1_xctr -= bt->gauss1_xctr;
   bt->lorentz1_yctr -= bt->gauss1_yctr;
   bt->gauss1_xctr = 0;
   bt->gauss1_yctr = 0;

   /* adjust the parameters */
   bt->gauss1_sigma *= fwhm_to_sigma;
   bt->gauss2_sigma *= fwhm_to_sigma;
   bt->lorentz1_hwhm *= 0.5;

   /* Compute the weights.  At this point, the wgt values are the amplitudes,
    * defined as the function values at r=0.  The functional forms are
    * normalized to have an integrated value of 1:
    *
    *    F(r) = N*f(r)
    *    1 = \int F(r) 2\pi r
    *
    * For the forms here, the shape function f(r) at r=0 is 1.
    *   ==> F(0) = N.
    *
    * it follows that the weight w is related to the desired amplitude by
    *
    *   amp = w*F(0) = w*N   ==> w = amp/N
    */

   norm = 1.0/(2*PI*bt->gauss1_sigma*bt->gauss1_sigma);
   bt->gauss1_wgt /= norm;

   norm = 1.0/(2*PI*bt->gauss2_sigma*bt->gauss2_sigma);
   bt->gauss2_wgt /= norm;

   g2 = bt->lorentz1_hwhm; g2 = g2*g2;
   norm = 1.0 / (g2 * PI * log(1 + bt->lorentz1_rmax*bt->lorentz1_rmax/g2));
   bt->lorentz1_wgt /= norm;

   norm = bt->gauss1_wgt + bt->gauss2_wgt + bt->lorentz1_wgt;
   if (norm > 0)
     {
	bt->gauss1_wgt /= norm;
	bt->gauss2_wgt /= norm;
	bt->lorentz1_wgt /= norm;
     }
   return 0;
}

void _marx_hrc_blur_close (Marx_HRC_Blur_Parm_Type *bt)
{
   (void) bt;
}

Marx_HRC_Blur_Parm_Type *_marx_hrc_blur_open (Param_File_Type *p, int det)
{
   Marx_HRC_Blur_Parm_Type *bt;
   Param_Table_Type *t;
   double dx, dy, pixelsize;

   if (det == 'I')
     {
	bt = &HRC_I_Blur_Parms;
	t = HRC_I_Blur_Parm_Table;
	if (-1 == _marx_hrc_i_get_pixel_size (&dx, &dy))
	  return NULL;
     }
   else if (det == 'S')
     {
	bt = &HRC_S_Blur_Parms;
	t = HRC_S_Blur_Parm_Table;
	if (-1 == _marx_hrc_s_get_pixel_size (&dx, &dy))
	  return NULL;
     }
   else
     {
	marx_error ("_marx_hrc_blur_open: Unsupported HRC detector : '%c'", det);
	return NULL;
     }

   /* parameter file values uses 1/4 pixels */
   dx *= 0.25;
   dy *= 0.25;
   pixelsize = 0.5 * (dx + dy);

   if (-1 == pf_get_parameters (p, t))
     return NULL;

   if (-1 == fixup_blur_parms (bt, pixelsize))
     return NULL;

   return bt;
}

static double rand_lorentz2d_radius (double hwhm, double rmax)
{
   double c;

   c = JDMrandom ();
   rmax /= hwhm;
   return hwhm * sqrt (expm1(c*log1p(rmax*rmax)));
}

static double rand_gauss2d_radius (double sigma)
{
   double c;

   do
     {
	c = JDMrandom ();
     }
   while (c == 0.0);
   return sigma * sqrt (-2*log(c));
}

void _marx_hrc_blur_position (Marx_HRC_Blur_Parm_Type *bt, double *dx, double *dy)
{
   double r, theta;
   double x_0, y_0;

   if ((_Marx_Det_Ideal_Flag) || (bt == NULL))
     return;

   r = JDMrandom ();
   if (r < bt->gauss1_wgt)
     {
	r = rand_gauss2d_radius (bt->gauss1_sigma);
	x_0 = bt->gauss1_xctr;
	y_0 = bt->gauss1_yctr;
     }
   else if (r < bt->gauss1_wgt + bt->gauss2_wgt)
     {
	r = rand_gauss2d_radius (bt->gauss2_sigma);
	x_0 = bt->gauss2_xctr;
	y_0 = bt->gauss2_yctr;
     }
   else
     {
	r = rand_lorentz2d_radius (bt->lorentz1_hwhm, bt->lorentz1_rmax);
	x_0 = bt->lorentz1_xctr;
	y_0 = bt->lorentz1_yctr;
     }

   theta = (2.0 * PI) * JDMrandom ();

   *dx += x_0 + r * cos(theta);
   *dy += y_0 + r * sin(theta);
   if (_Marx_Det_Extend_Flag == 0)
     {
	if (*dx < 0.0) *dx = 0.0;
	if (*dy < 0.0) *dy = 0.0;
     }
}

