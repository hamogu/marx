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
/* Dither routines */

#include "config.h"
#include "marx-feat.h"

#include <stdio.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <string.h>

#include <jdmath.h>
#include <pfile.h>
#include <jdfits.h>

#include "marx.h"
#include "_marx.h"

int _Marx_Dither_Mode = _MARX_DITHER_MODE_NONE;

static double Ra_Amp;
static double Dec_Amp;
static double Roll_Amp;
static double Ra_Period;
static double Dec_Period;
static double Roll_Period;
static double Ra_Phase;
static double Dec_Phase;
static double Roll_Phase;

static double Nominal_Ra;
static double Nominal_Dec;
static double Nominal_Roll;
static JDMVector_Type Nominal_Pointing;

static double Pointing_Ra;
static double Pointing_Dec;
static double Pointing_Roll;

static double Pointing_Offset_Y;
static double Pointing_Offset_Z;

static double Aspect_Blur;

static char *Dither_Model;
static char *Dither_File_Name;

/* The last parameter is the value of the roll to use when applying the
 * dither.  SAOSAC rays are already rolled, so this value will be 0.
 */
static int (*Get_Dither_Function)(double, Marx_Dither_Type *, double *, double *, double *);


static int (*Get_Dither_Par_Means)(double, double *, double *, double *);

static Param_Table_Type Dither_Parm_Table [] =
{
     {"DitherModel",		PF_STRING_TYPE, &Dither_Model},
     {"DitherFile",		PF_FILE_TYPE,	&Dither_File_Name},

     {"DitherAmp_RA",		PF_REAL_TYPE,	&Ra_Amp},
     {"DitherAmp_Dec",		PF_REAL_TYPE,	&Dec_Amp},
     {"DitherAmp_Roll",		PF_REAL_TYPE,	&Roll_Amp},

     {"DitherPeriod_RA",	PF_REAL_TYPE,	&Ra_Period},
     {"DitherPeriod_Dec",	PF_REAL_TYPE,	&Dec_Period},
     {"DitherPeriod_Roll",	PF_REAL_TYPE,	&Roll_Period},

     {"DitherPhase_RA",		PF_REAL_TYPE,	&Ra_Phase},
     {"DitherPhase_Dec",	PF_REAL_TYPE,	&Dec_Phase},
     {"DitherPhase_Roll",	PF_REAL_TYPE,	&Roll_Phase},

     {"AspectBlur",		PF_REAL_TYPE,	&Aspect_Blur},
     {"Ra_Nom",			PF_REAL_TYPE,	&Nominal_Ra},
     {"Dec_Nom",		PF_REAL_TYPE,	&Nominal_Dec},
     {"Roll_Nom",		PF_REAL_TYPE,	&Nominal_Roll},
     {"PointingOffsetY",	PF_REAL_TYPE,	&Pointing_Offset_Y},
     {"PointingOffsetZ",	PF_REAL_TYPE,	&Pointing_Offset_Z},

     {NULL, 0, NULL}
};

/* Convert angle to value to range [-PI,PI) */
static void normalize_angle (double *xp)
{
   double x = *xp;

   x = fmod (*xp, 2*PI);

   if (x >= PI) x -= (2*PI);
   if (x < -PI) x += (2*PI);

   *xp = x;
}

static int setup_pointing_ra_dec (void)
{
   JDMVector_Type p, nom, a, d;
   double cos_delta;
   double alpha, delta;

   alpha = Pointing_Offset_Y;
   delta = Pointing_Offset_Z;

   /* Compute (p,a,d) as orthonormal system */
   JDMv_spherical_to_triad (PI/2.0 - Nominal_Dec, Nominal_Ra, &nom, &d, &a);
   /* Flip the orientation of 'd' since it is opposite of spherical
    * coordinate theta_hat.
    */
   d.x = -d.x; d.y = -d.y; d.z = -d.z;

   /* Here nom is now in the nominal direction.  We wish to find the
    * pointing direction.  Find it by taking the roll into account.
    */
   cos_delta = cos (delta);
   p = JDMv_ax1_bx2_cx3 (cos_delta * cos(alpha), nom,
			 cos_delta * sin(alpha), a,
			 sin (delta), d);

   p = JDMv_rotate_unit_vector (p, nom, Nominal_Roll);
   JDMv_unit_vector_to_spherical (p, &Pointing_Dec, &Pointing_Ra);
   Pointing_Dec = PI*0.5 - Pointing_Dec;
   Pointing_Roll = Nominal_Roll;

   return 0;
}

static int get_dither_par_means_0(double t, double *dy, double *dz, double *dtheta)
{
  *dy = 0;
  *dz = 0;
  *dtheta = 0;

  return 1;
}


static int get_internal_dither (double t, Marx_Dither_Type *d,
				double *rap, double *decp, double *rollp)
{
   t = (2.0*PI) * t;

   /* Note that the ra/dec values are offsets */
   *rap = d->ra = Ra_Amp * sin (t/Ra_Period + Ra_Phase);
   *decp = d->dec = Dec_Amp * sin (t/Dec_Period + Dec_Phase);
   *rollp = d->roll = Nominal_Roll + Roll_Amp * sin (t/Roll_Period + Roll_Phase);

   d->dy = 0;
   d->dz = 0;
   d->dtheta = 0;

   return 1;
}

static int get_zeroamp_internal_dither (double t, Marx_Dither_Type *d,
					double *rap, double *decp, double *rollp)
{
   (void) t;
   *rap = d->ra = 0.0;
   *decp = d->dec = 0.0;
   d->roll = Nominal_Roll;	       /* roll still gets applied */
   *rollp = 0.0;

   return 1;
}

static int init_internal_dither (int dither_flag)
{
   _Marx_Dither_Mode = _MARX_DITHER_MODE_INTERNAL;

   if (dither_flag & _MARX_DITHER_ZERO_AMP)
     {
	Ra_Amp = Dec_Amp = Roll_Amp = 0;
	Get_Dither_Function = get_zeroamp_internal_dither;

	marx_message ("[Using INTERNAL dither model with 0 amplitudes]\n");
     }
   else
     {
	Ra_Amp = Ra_Amp * (PI/(180.0 * 3600));
	Dec_Amp = Dec_Amp * (PI/(180.0 * 3600));
	Roll_Amp = Roll_Amp * (PI/(180.0 * 3600));
	Get_Dither_Function = get_internal_dither;
	marx_message ("[Using INTERNAL dither model]\n");
     }
   Get_Dither_Par_Means = get_dither_par_means_0;
   if (-1 == setup_pointing_ra_dec ())
     return -1;

   return 0;
}

static int init_no_dither (int dither_flags)
{
   (void) dither_flags;
   _Marx_Dither_Mode = _MARX_DITHER_MODE_NONE;
   Get_Dither_Function = NULL;
   Get_Dither_Par_Means = get_dither_par_means_0;
   Pointing_Roll = 0.0;

   if (-1 == setup_pointing_ra_dec ())
     return -1;

   return 0;
}

static double Start_Time;
static JDMVector_Type Nominal_Pointing;

typedef struct
{
   JDFits_BTable_Read_Type *bt;
   double t0, ra0, dec0, roll0, dy0, dz0, dtheta0;
   double t1, ra1, dec1, roll1, dy1, dz1, dtheta1;
  double sum_dy, sum_dz, sum_dtheta;
  int num_steps;
}
Aspsol_Type;

static Aspsol_Type Aspsol;

static void close_aspsol (void)
{
   if (Aspsol.bt != NULL)
     {
	jdfits_simple_close_btable (Aspsol.bt);
	Aspsol.bt = NULL;
     }
}

/** Read the next point in the asol fits file and set variable Aspsol to it.
  * There is no randomn access, this will always read exactly the next line
  * in the fits file.
  * Aspsol is a module level variable that is defined when the FILE dither
  * model is in use.
  */ 
static int get_single_aspsol_point (void)
{
   double buf[7];
   JDMVector_Type p;

   Aspsol.t0 = Aspsol.t1;
   Aspsol.ra0 = Aspsol.ra1;
   Aspsol.dec0 = Aspsol.dec1;
   Aspsol.roll0 = Aspsol.roll1;
   Aspsol.dy0 = Aspsol.dy1;
   Aspsol.dz0 = Aspsol.dz1;
   Aspsol.dtheta0 = Aspsol.dtheta1;

   if (-1 == jdfits_simple_d_read_btable (Aspsol.bt, buf))
     {
	close_aspsol ();
	return -1;
     }

   Aspsol.t1 = buf[0];
   Aspsol.ra1 = buf[1] * (PI/180.0);
   Aspsol.dec1 = buf[2] * (PI/180.0);
   Aspsol.roll1 = buf[3] * (PI/180.0);
   Aspsol.dy1 = buf[4];
   Aspsol.dz1 = buf[5];
   Aspsol.dtheta1 = buf[6] * (PI/180.0);

   /* Keep cumulative sum, because we eventually need the mean of those values*/
   Aspsol.num_steps++;
   Aspsol.sum_dy += Aspsol.dy1;
   Aspsol.sum_dz += Aspsol.dz1;
   Aspsol.sum_dtheta += Aspsol.dtheta1;

   /* We need to convert the ra/dec information to the unrolled values
    * since in the unrolled frame they are equivalent to yaw/pitch.  Note
    * that ra/dec system differs from the spherical system in the definition
    * of the polar angle.
    */
   p = JDMv_spherical_to_vector (1.0, (PI/2.0) - Aspsol.dec1, Aspsol.ra1);
   p = JDMv_rotate_unit_vector (p, Nominal_Pointing, -Aspsol.roll1);
   JDMv_unit_vector_to_spherical (p, &Aspsol.dec1, &Aspsol.ra1);
   Aspsol.dec1 = (PI/2.0) - Aspsol.dec1;

   normalize_angle (&Aspsol.ra1);
   normalize_angle (&Aspsol.dec1);
   normalize_angle (&Aspsol.roll1);
   normalize_angle (&Aspsol.dtheta1);

   marx_compute_ra_dec_offsets (Nominal_Ra, Nominal_Dec,
				Aspsol.ra1, Aspsol.dec1,
				&Aspsol.ra1, &Aspsol.dec1);

   /* Aspsol.roll1 -= Nominal_Roll; */

   Aspsol.t1 -= Start_Time;
   return 0;
}

static int get_aspsol_point (double t)
{
   while (t >= Aspsol.t1)
     {
	if (-1 == get_single_aspsol_point ())
	  return -1;
     }
   return 1;
}

static int get_aspsol_dither (double t, Marx_Dither_Type *d,
			      double *rap, double *decp, double *rollp)
{
   double dt;

   if (-1 == get_aspsol_point (t))
     return -1;

   dt = Aspsol.t1 - Aspsol.t0;
   if (dt != 0)
     dt = (t - Aspsol.t0)/dt;

   *rap = d->ra = Aspsol.ra0 + dt * (Aspsol.ra1 - Aspsol.ra0);
   *decp = d->dec = Aspsol.dec0 + dt * (Aspsol.dec1 - Aspsol.dec0);
   *rollp = d->roll = Aspsol.roll0 + dt * (Aspsol.roll1 - Aspsol.roll0);
   d->dy = Aspsol.dy0 + dt * (Aspsol.dy1 - Aspsol.dy0);
   d->dz = Aspsol.dz0 + dt * (Aspsol.dz1 - Aspsol.dz0);
   d->dtheta = Aspsol.dtheta0 + dt * (Aspsol.dtheta1 - Aspsol.dtheta0);

   return 1;
}


static int get_aspsol_dither_mean(double total_time, double *dy, double *dz, double *dtheta)
{
  int iret;

  // Make sure we read the file up to the total exposure time
  iret = get_aspsol_point (total_time);
  // marx_message("get_aspsol_dither_mean(): return value of get_aspsol_point() = %d\n",iret);
  /*   if (-1 == get_aspsol_point (total_time))
     return -1;
  */
  if(iret==-1) return iret;

  marx_message("get_aspsol_dither_mean(): Aspsol.num_steps = %d\n",Aspsol.num_steps);
   if ( 0 == Aspsol.num_steps)
     {
       marx_error("Exposure is too short to calculate mean ASPSOL values.");
       return -1;
     }
   *dy = Aspsol.sum_dy / Aspsol.num_steps;
   *dz = Aspsol.sum_dz / Aspsol.num_steps;
   *dtheta = Aspsol.sum_dtheta / Aspsol.num_steps;
   return iret;
}

/* Record the dither state but do not actually dither the rays, which are
 * already assumed to be dithered with the asol file.
 */
static int get_aspsol_dither_record_only (double t, Marx_Dither_Type *d,
					  double *rap, double *decp, double *rollp)
{
   int status;

   if (1 != (status = get_aspsol_dither (t, d, rap, decp, rollp)))
     return status;

   *rap = 0;
   *decp = 0;
   *rollp = 0;

   return 1;
}

static int get_keyword_double (JDFits_Type *ft, char *name, double *v)
{
   JDFits_Keyword_Type *k;

   if ((NULL == (k = jdfits_find_keyword (ft, name)))
       || (-1 == jdfits_extract_double (k, v)))
     {
	marx_error ("Unable to find keyword %s in ASPSOL file\n", name);
	return -1;
     }

   return 0;
}

static int init_aspsol_dither (int dither_flag)
{
   JDFits_BTable_Read_Type *bt;

   marx_message ("Opening ASPSOL fits file %s\n", Dither_File_Name);

   bt = jdfits_simple_open_btable (Dither_File_Name,
				   "ASPSOL",
				   7,
				   "time",
				   "ra",
				   "dec",
				   "roll",
				   "dy",
				   "dz",
				   "dtheta");

   if (bt == NULL)
     {
	marx_error ("Unable to open a proper ASPSOL file called %s",
		    Dither_File_Name);
	return -1;
     }

   if ((-1 == get_keyword_double (bt->ft, "RA_NOM", &Nominal_Ra))
       || (-1 == get_keyword_double (bt->ft, "DEC_NOM", &Nominal_Dec))
       || (-1 == get_keyword_double (bt->ft, "ROLL_NOM", &Nominal_Roll))
       || (-1 == get_keyword_double (bt->ft, "TSTART", &Start_Time)))
     {
	jdfits_simple_close_btable (bt);
	return -1;
     }
   Aspsol.bt = bt;

   Nominal_Ra *= PI/180.0;
   Nominal_Dec *= PI/180.0;
   Nominal_Roll *= PI/180.0;
   if (Nominal_Ra >= PI) Nominal_Ra -= (2.0*PI);
   if (Nominal_Dec >= PI) Nominal_Dec -= (2.0*PI);
   if (Nominal_Roll >= PI) Nominal_Roll -= (2.0*PI);

   Nominal_Pointing = JDMv_spherical_to_vector (1.0, PI/2.0 - Nominal_Dec, Nominal_Ra);

   if (-1 == setup_pointing_ra_dec ())
     {
	close_aspsol ();
	return -1;
     }

   if (-1 == get_single_aspsol_point ())
     {
	marx_error ("ASPSOL file contains no data");
	close_aspsol ();
	return -1;
     }

   Aspsol.t1 -= Start_Time;

   _Marx_Dither_Mode = _MARX_DITHER_MODE_ASPSOL;
   if (dither_flag & _MARX_DITHER_RECORD_ONLY)
     Get_Dither_Function = get_aspsol_dither_record_only;
   else
     Get_Dither_Function = get_aspsol_dither;

   Get_Dither_Par_Means = get_aspsol_dither_mean;
   marx_message ("[Using ASPSOL dither model]\n");

   return 0;
}

void _marx_close_dither (void)
{
   if (_Marx_Dither_Mode == _MARX_DITHER_MODE_ASPSOL)
     close_aspsol ();
}

int _marx_init_dither (Param_File_Type *pf, int dither_flags, double *yoff, double *zoff)
{
   _Marx_Dither_Mode = _MARX_DITHER_MODE_NONE;
   Get_Dither_Function = NULL;
   Get_Dither_Par_Means = NULL;

   if (-1 == pf_get_parameters (pf, Dither_Parm_Table))
     {
	marx_error ("error getting dither parameters");
	return -1;
     }

   Aspect_Blur *= PI/(180.0*3600.0);
   if (Aspect_Blur < 0)
     Aspect_Blur = 0.0;

   Pointing_Offset_Y *= PI/(180.0*3600.0);
   Pointing_Offset_Z *= PI/(180.0*3600.0);

   *yoff = Pointing_Offset_Y;
   *zoff = Pointing_Offset_Z;

   Nominal_Ra *= PI/180.0;
   Nominal_Dec *= PI/180.0;
   Nominal_Roll *= PI/180.0;
   normalize_angle (&Nominal_Roll);

   if ((dither_flags & _MARX_DITHER_UNSUPPORTED)
       || (0 == strcmp (Dither_Model, "NONE")))
     return init_no_dither (dither_flags);

   if (0 == strcmp (Dither_Model, "INTERNAL"))
     return init_internal_dither (dither_flags);

   if (0 == strcmp (Dither_Model, "FILE"))
     return init_aspsol_dither (dither_flags);

   marx_error ("DitherModel = %s is not supported", Dither_Model);
   return -1;
}

#define VERY_TINY_NUMBER 1e-20

static JDMVector_Type apply_dither (double ra, double dec, double roll,
				    JDMVector_Type p)
{
   double cos_dec, sin_dec, cos_ra, sin_ra;
   double cos_theta, sin_theta;
   JDMVector_Type n;

   /* First of all, roll the spacecraft about the x axis */
   p = JDMv_rotate_unit_vector (p, JDMv_vector (1, 0, 0), -roll);

   /* Now dither in the y-z plane by ra, dec.  Here, these are aspect offsets. */
   cos_ra = cos (ra);
   sin_ra = sin (ra);
   cos_dec = cos (dec);
   sin_dec = sin (dec);

   /* Here, n is chosen to be orthogonal to the instantaneous pointing.
    * See my aspect-offsets memo (eq 4).
    */
   cos_theta = cos_dec * cos_ra;
   n = JDMv_vector (0, sin_dec, -cos_dec * sin_ra);
   sin_theta = JDMv_length (n);
   if (sin_theta <= VERY_TINY_NUMBER)
     return p;

   n.x /= sin_theta;
   n.y /= sin_theta;
   n.z /= sin_theta;

   return JDMv_rotate_unit_vector1 (p, n, cos_theta, sin_theta);
}

static JDMVector_Type unapply_dither (double ra, double dec, double roll,
				      JDMVector_Type p)
{
   double cos_dec, sin_dec, cos_ra, sin_ra;
   double cos_theta, sin_theta;
   JDMVector_Type n;

   cos_ra = cos (ra);
   sin_ra = sin (ra);
   cos_dec = cos (dec);
   sin_dec = sin (dec);

   cos_theta = cos_dec * cos_ra;
   n = JDMv_vector (0, -sin_dec, cos_dec * sin_ra);
   sin_theta = JDMv_length (n);
   if (sin_theta >= VERY_TINY_NUMBER)
     {
	n.x /= sin_theta;
	n.y /= sin_theta;
	n.z /= sin_theta;
	p = JDMv_rotate_unit_vector1 (p, n, cos_theta, sin_theta);
     }

   return JDMv_rotate_unit_vector (p, JDMv_vector (1, 0, 0), roll);
}

static int dither_ray (Marx_Photon_Attr_Type *at, double t)
{
   int status;
   Marx_Dither_Type *d;
   double delta_ra, delta_dec;
   double ra, dec, roll;

   t += at->arrival_time;
   d = &at->dither_state;

   status =  (*Get_Dither_Function)(t, d, &ra, &dec, &roll);
   if (status != 1)
     return status;

   delta_ra = Aspect_Blur * JDMgaussian_random ();
   delta_dec = Aspect_Blur * JDMgaussian_random ();

   at->p = apply_dither (ra + delta_ra, dec + delta_dec, roll, at->p);
   return 1;
}

int _marx_dither_photons (Marx_Photon_Type *pt, unsigned int *num_dithered)
{
   Marx_Photon_Attr_Type *at, *at_max;
   unsigned int num;
   double time_offset;

   if (Get_Dither_Function == NULL)
     {
	*num_dithered = pt->n_photons;
	return 0;
     }

   time_offset = pt->start_time;
   at = pt->attributes;
   at_max = at + pt->n_photons;

   num = 0;
   while (at < at_max)
     {
	if (1 != dither_ray (at, time_offset))
	  break;

	at++;
	num++;
     }

   *num_dithered = num;
   return 0;
}

#if 0
static void apply_dither_blur (JDMVector_Type *p)
{
   JDMVector_Type n;
   double rnd;

   if (Dither_Blur <= 0.0)
     return;

   /* Construct a vector normal to p.  We want:
    *
    *    0 = p_x n_x + p_y n_y + p_z n_z
    *
    * We expect p_x to be close to 1.  So, choose n_z = 0, n_y = 1,
    * and n_x = -p_y/p_x ==> p.n = -p_y + p_y = 0.
    */

   n.z = 0.0;
   n.y = 1.0;
   n.x = -(p->y/p->x);
   JDMv_normalize (&n);

   /* Now rotate n by a random angle about p */
   n = JDMv_rotate_unit_vector (n, *p, (2.0 * PI) * JDMrandom ());
#if 0
   rnd = JDMgaussian_random ();
#else
   /* We want the blur to be a gaussian in sky coordinates and NOT in the
    * direction of the ray.  So, we have to do what was done for the gaussian
    * source model.
    */
   do
     {
	rnd = JDMrandom ();
     }
   while (rnd == 0.0);

   rnd = sqrt (-log (rnd));
#endif
   *p = JDMv_rotate_unit_vector (*p, n, Dither_Blur * rnd);
}
#endif

void marx_undither_mnc (JDMVector_Type *mnc, Marx_Dither_Type *d)
{
   /* apply_dither_blur (mnc); */
   *mnc = unapply_dither (d->ra, d->dec, d->roll, *mnc);
}

void _marx_ray_to_sky_ra_dec (Marx_Photon_Attr_Type *at, double *ra, double *dec)
{
   double flen;
   JDMVector_Type p;

   /* I am going to cheat by using the xpos,ypos, and zpos coordinates, and
    * assume a focal length.  The marx origin is at the nominal focal point.
    */
   flen = Marx_Focal_Length;

   p = JDMv_diff (at->x, JDMv_vector (flen, 0, 0));
   /* Now p is the vector to the point through the mirror node.  That is,
    * it is a vector in the mirror nodal system.  Make it a unit
    * vector.
    */
   JDMv_normalize (&p);
   marx_undither_mnc (&p, &at->dither_state);

   marx_mnc_to_ra_dec (&p, ra, dec);
}

int marx_get_nominal_pointing (double *ra_nom, double *dec_nom, double *roll_nom)
{
   *ra_nom = Nominal_Ra;
   *dec_nom = Nominal_Dec;
   *roll_nom = Nominal_Roll;

   return 0;
}

int marx_get_pointing (double *ra, double *dec, double *roll)
{
   *ra = Pointing_Ra;
   *dec = Pointing_Dec;
   *roll = Pointing_Roll;

   return 0;
}

/* This function is called by the SAOSAC source model.  It gives the value
 * of the time of the first ray in the ray file.  Use it instead of the
 * TSTART value in the asol file.
 */
void _marx_dither_set_ray_tstart (double val)
{
   Start_Time = val;
}

int marx_average_dither(double t, double *dy, double *dz, double *dtheta)
{
  int status;
  status = (*Get_Dither_Par_Means)(t, dy, dz, dtheta);

  // marx_message("marx_average_dither(): status = %d\n",status);

  /* if (status != 1)
     return status;
  */
  return status;
}
