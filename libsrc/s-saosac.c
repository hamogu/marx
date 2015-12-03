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
#include <jdfits.h>

#include "marx.h"
#include "_marx.h"
#include "source.def"

static JDFits_BTable_Read_Type *Saosac_Bin_Table;
static int File_Has_Time_Column;

#define SCALE_SAOSAC_WEIGHTS 0
#if SCALE_SAOSAC_WEIGHTS
static double Wgts_Scale_Factor = 1.0;
#endif

static int Use_Color_Rays = 0;

#if SCALE_SAOSAC_WEIGHTS
static int compute_weights_scale_factor (char *file, double *wp)
{
   JDFits_BTable_Read_Type *bt;
   char *extname;
   static char *columns [] =
     {
	  "RT_WGHT",
     };
   double w, wmax;

   marx_message ("Reading RT_WGHT column from %s\n", file);

   extname = "RAYTRACE";

   bt = jdfits_simple_aopen_btable (file,
				    extname,
				    1, columns);

   if (bt == NULL)
     {
	marx_error ("Unable to open proper SAOSAC rayfile called %s",
		    file);
	return -1.0;
     }

   wmax = 0.0;

   while (-1 != jdfits_simple_d_read_btable (bt, &w))
     {
	if (w > wmax)
	  wmax = w;
     }

   jdfits_simple_close_btable (bt);

   *wp = 1.0/wmax;
   marx_message ("Scaling SAOSAC Weights by %g\n", *wp);
   return 0;
}
#endif

static JDFits_BTable_Read_Type *open_saosac_fits_file (char *file)
{
   JDFits_BTable_Read_Type *bt;
   JDFits_Type *f;
   char *extname;
   static char *columns [] =
     {
	"RT_X",
	  "RT_Y",
	  "RT_Z",
	  "RT_COSX",
	  "RT_COSY",
	  "RT_COSZ",
	  "RT_KEV",
	  "RT_WGHT",
	  "RT_TIME"
     };

   marx_message ("Opening SAOSAC fits file %s\n", file);

   extname = "RAYTRACE";

   /* See if time column is present */
   if (NULL == (f = jdfits_open_binary_table (file, extname)))
     {
	marx_error ("Unable to open proper SAOSAC rayfile called %s with extname=%s",
		    file, extname);
	return NULL;
     }
   if (1 == jdfits_bintable_column_exists (f, "RT_TIME"))
     File_Has_Time_Column = 1;
   else
     File_Has_Time_Column = 0;

   (void) jdfits_close_file (f);

   bt = jdfits_simple_aopen_btable (file,
				    extname,
				    (File_Has_Time_Column ? 9 : 8),
				    columns);

   if (bt == NULL)
     {
	marx_error ("Unable to open proper SAOSAC rayfile called %s",
		    file);
	return NULL;
     }

   return bt;
}

static int saosac_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;

   if (NULL == Saosac_Bin_Table)
     return -1;

   return 0;
}

/*}}}*/

static int saosac_close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;

   if (Saosac_Bin_Table == NULL)
     return 0;

   jdfits_simple_close_btable (Saosac_Bin_Table);
   Saosac_Bin_Table = NULL;

   return 0;
}

/*}}}*/

static double Start_Time = 0.0;
static int saosac_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				   unsigned int num, unsigned int *num_created)
{
   unsigned int num_read;
   Marx_Photon_Attr_Type *at;
   JDFits_BTable_Read_Type *bt;
   double buf [9];
   double this_time;
   int (*efun) (Marx_Spectrum_Type *, double *);

   bt = Saosac_Bin_Table;
   efun = st->spectrum.energy_function;

   if (bt == NULL)
     return -1;

   at = pt->attributes;

   num_read = 0;
   this_time = 0.0;

   while (num_read < num)
     {
	double t, x, y, z, r, px, py, pz, w;

	if (-1 == jdfits_simple_d_read_btable (bt, buf))
	  break;

	w = buf[7];
#if SCALE_SAOSAC_WEIGHTS
	w *= Wgts_Scale_Factor;
#endif
	if (JDMrandom () >= w)
	  at->flags |= PHOTON_MIRROR_VBLOCKED;

	/* The rays are in the XRCF coord system.  This is the same as MARX
	 * except that the origin is at the fore end of the CAP.  Here I will
	 * project them to the CAP
	 */
	x = buf[0];
	px = buf[3];
	py = buf[4];
	pz = buf[5];

	/* According to Diab, the rays in the saosac file are projected to
	 * the its best focus position at 1.49 keV.  This position should
	 * correspond to the origin of the marx system.
	 */
	/* Ummm...  I am not so sure about the coordinate system. */
	t = -x / px;

	at->x.x = _Marx_HRMA_Cap_Position - 0.016750;
	at->x.y = y = buf[1] + py * t;
	at->x.z = z = buf[2] + pz * t;

	at->p.x = px;
	at->p.y = py;
	at->p.z = pz;

	/* Now deduce the mirror shell based on position of the ray.  The
	 * following values will be used:
	 *				      Radius, mm
	 *
	 * Mirror            Parabola             Intersect            Hyperbola
	 * Number       Front         Back          Plane         Front         Back
	 *
	 *    1     612.69078     600.34506     599.45019     598.52343     560.86506
	 *    3     493.40861     483.46567     482.85824     481.77216     451.45717
	 *    4     435.62935     426.85049     426.35484     425.27403     398.51324
	 *    6     323.81612     317.29022     316.96956     316.02294     296.13584
	 *
	 * In particular, the intersect plane values will be used.
	 */

	r = z * z + y * y;

	if (r < 138012.25)	       /* ((426 + 317)/2.0)^2  */
	  at->mirror_shell = 3;
	else if (r < 206116.000000)    /* (482 + 426)^2/4.0 */
	  at->mirror_shell = 2;
	else if (r < 292681.000000)    /* (599 + 483)^2/4 */
	  at->mirror_shell = 1;
	else
	  at->mirror_shell = 0;

	if (Use_Color_Rays)
	  {
	     if (-1 == (*efun) (&st->spectrum, &at->energy))
	       return -1;
	  }
	else
	  {
	     if (File_Has_Time_Column)
	       {
		  if (Start_Time == 0)
		    {
		       Start_Time = buf[8];
		       _marx_dither_set_ray_tstart (Start_Time);
		    }
		  this_time = buf[8] - Start_Time;
	       }
	     else
	       this_time += 1.0;

	     if (this_time < 0)
	       continue;

	     at->energy = buf[6];
	     at->arrival_time = this_time;
	  }

	at++;
	num_read++;
     }

   pt->history |= (MARX_ENERGY_OK
		   | MARX_X_VECTOR_OK
		   | MARX_P_VECTOR_OK);

   pt->history |= MARX_MIRROR_SHELL_OK;

   if (File_Has_Time_Column)
     {
	Start_Time += this_time;
	pt->history |= MARX_TIME_OK;
     }

   *num_created = num_read;

   return 0;
}

/*}}}*/

int marx_select_saosac_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
				char *name, unsigned int source_id)
{
   char file [PF_MAX_LINE_LEN];
   JDFits_BTable_Read_Type *bt;
#if SCALE_SAOSAC_WEIGHTS
   int scale_wgts;
#endif
   char dither_model[PF_MAX_LINE_LEN];

   (void) source_id; (void) name;
   st->open_source = saosac_open_source;
   st->create_photons = saosac_create_photons;
   st->close_source = saosac_close_source;

   if (-1 == pf_get_string (p, "SAOSACFile", file, sizeof (file)))
     return -1;
   if (-1 == pf_get_string (p, "DitherModel", dither_model, sizeof (dither_model)))
     return -1;
   if (0 == strcmp (dither_model, "FILE"))
     {
	/* Initialized from marx.c from the aspect file.  Yes this is ugly. */
	Start_Time = 0*_Marx_TStart_MJDsecs;
     }
   else Start_Time = 0;

#if 0
   if (-1 == pf_get_boolean (p, "SAOSAC_Color_Rays", &Use_Color_Rays))
     return -1;
#endif
   Use_Color_Rays = 0;		       /* nolonger supported */

#if SCALE_SAOSAC_WEIGHTS
   Wgts_Scale_Factor = 1.0;
   if (-1 == pf_get_boolean (p, "SAOSAC_Scale_Wgts", &scale_wgts))
     return -1;
   if (scale_wgts)
     {
	if (-1 == compute_weights_scale_factor (file, &Wgts_Scale_Factor))
	  return -1;
     }
#endif

   bt = open_saosac_fits_file (file);
   if (bt == NULL)
     return -1;

   Saosac_Bin_Table = bt;

   if (Use_Color_Rays)
     {
	return _marx_get_simple_specrum_parms (p, st, name);
     }
   st->spectrum.type = MARX_SAOSAC_SPECTRUM;

   return 0;
}

/*}}}*/
