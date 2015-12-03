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

#undef MARX_HAS_HRMA_PITCH_YAW
#define MARX_HAS_HRMA_PITCH_YAW 1

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

/*{{{ HRMA_Type structure and initialization */

typedef struct
{
   unsigned int mirror_number;

   double length_p;		       /* length of parabola */
   double length_h;		       /* length of hyperbola */
   double osac_z0_p;		       /* offset of P from CAP */
   double osac_z0_h;		       /* offset of H from CAP */
   double osac_y0_p;		       /* offset of P from CAP */
   double osac_y0_h;		       /* offset of H from CAP */
   double osac_x0_p;		       /* offset of P from CAP */
   double osac_x0_h;		       /* offset of H from CAP */
   double osac_p_p;		       /* P parameter for parabola */
   double osac_p_h;
   double osac_k_p;		       /* K parameter */
   double osac_k_h;
   double osac_r_p;		       /* r parameter */
   double osac_r_h;
   double osac_az_p;		       /* azmith for P */
   double osac_az_h;		       /* azmith for H */
   double osac_el_p;		       /* elmith for P */
   double osac_el_h;		       /* elmith for H */

   double xoffset;		       /* tweak to Cap_Position to shift focus */
   double p_blur;
   double h_blur;
#if MARX_HAS_WFOLD
   double h_scat_factor;
   double p_scat_factor;
#endif
   /* These are computed */
   /* The equation of the conic is rearranged in form
    * x^2 + y^2 + z^2 = a x^2 + b x + c
    * (0,0,0) is at the center of the conic.  This is the origin of the OSAC
    * coordinate system.
    */
   double conic_a_p;		       /* 1-osac_p_p  ;was: osac_r_p ^2 */
   double conic_b_p;		       /* -2.0 * osac_k_p */
   double conic_c_p;		       /* osac_r_p^2  ;was: 1.0 - osac_p_p */
   double conic_a_h;		       /* 1-osac_p_h  ;was: osac_r_h ^2 */
   double conic_b_h;		       /* -2.0 * osac_k_h */
   double conic_c_h;		       /* osac_r_h^2  ;was: 1.0 - osac_p_h */
   double conic_xmin_h;		       /* -length_h / 2 */
   double conic_xmax_h;		       /* +length_h / 2 */
   double conic_xmin_p;		       /* -length_p / 2 */
   double conic_xmax_p;		       /* +length_p / 2 */

   double front_position;	       /* front aperature position in Marx system */
   /* These vectors translate the ray to the OSAC coordinate system.
    * The sign is such that SAOCAC = MARX + to_saosac.
    */
   JDMVector_Type to_osac_h;
   JDMVector_Type to_osac_p;

   double area_fraction;	       /* cumul fractional area of aperature */
   double min_radius, max_radius;      /* aperature limits */
   unsigned int shutter_bitmap;
   unsigned int num_open_shutters;

   /* These are the energy dependent correction factors used to match the
    * marx hrma eff-area to the caldb values.
    */
   float *correction_energies;
   float *correction_factors;
   unsigned int num_correction_factors;

#if MARX_HAS_WFOLD
   Marx_WFold_Table_Type *p_wfold_table;
   Marx_WFold_Table_Type *h_wfold_table;
#endif
#if MARX_HAS_HRMA_PITCH_YAW
   JDM_3Matrix_Type fwd_matrix_p;
   JDM_3Matrix_Type bwd_matrix_p;
   JDM_3Matrix_Type fwd_matrix_h;
   JDM_3Matrix_Type bwd_matrix_h;
#endif
}
HRMA_Type;

static HRMA_Type HRMA_Mirrors [MARX_NUM_MIRRORS] =
{
     {
	1,
	842.2150,		       /* length_p */
	842.1920,		       /* length_h */
	-426.5761,		       /* osac_z0_p */
	481.0146,		       /* osac_z0_h */
	0.2151,			       /* osac_y0_p */
	-0.2060,		       /* osac_y0_h */
	0.1239,			       /* osac_x0_p */
	-0.1154,		       /* osac_x0_h */
	0.0,			       /* osac_p_p */
	-1.7797716637950739e-03,       /* osac_p_h */
	-8.9333113530131421E+00,       /* osac_k_p */
	-26.0506034413416841,	       /* osac_k_h */
	606.86080963697918,	       /* osac_r_p */
	579.89015840093919,	       /* osac_r_h (rho0) */
	/* These come from table 6 of "Focus and Alignment of AXAF Optics"
	 * by Gaetz, et. al.
	 */
	0.0,			       /* osac_az_p = */
	2.4194219,		       /* osac_az_h = -tilt_y XRFC */
	0.0,			       /* osac_el_p */
	4.4454479		       /* osac_el_h = -tilt_z XRCF */
     },
     {
	3,
	842.2080,		       /* length_p */
	842.1970,		       /* length_h */
	-436.7098,		       /* osac_z0_p */
	480.9282,		       /* osac_z0_h */
	0.2437,			       /* osac_y0_p */
	-0.2345,		       /* osac_y0_h */
	0.08675,		       /* osac_x0_p */
	-0.08365,		       /* osac_x0_h */
	0.0,			       /* osac_p_p */
	-1.1532395834759916E-03,       /* osac_p_h */
	-5.7939624154424676E+00,       /* osac_k_p */
	-1.6875942397594130E+01,       /* osac_k_h */
	4.8846244215611011E+02,	       /* osac_r_p */
	4.6664379784205374E+02,	       /* osac_r_h */
	0.0,			       /* osac_az_p */
	1.8542174,		       /* osac_az_h */
	0.0,			       /* osac_el_p */
	4.9943249		       /* osac_el_h */
     },
     {
	4,
	842.2080,		       /* length_p */
	842.2250,		       /* length_h */
	-440.3572,		       /* osac_z0_p */
	480.8279,		       /* osac_z0_h */
	0.2168,			       /* osac_y0_p */
	-0.2065,		       /* osac_y0_h */
	0.08634,		       /* osac_x0_p */
	-0.08386,		       /* osac_x0_h */
	0.0,			       /* osac_p_p */
	-8.9864417477996457E-04,       /* osac_p_h */
	-4.5165799273846270E+00,       /* osac_k_p */
	-1.3150318066441841E+01,       /* osac_k_h */
	4.3126225933154404E+02,	       /* osac_r_p */
	4.1191935912458598E+02,	       /* osac_r_h */
	0.0,			       /* osac_az_p */
	1.8468078,		       /* osac_az_h */
	0.0,			       /* osac_el_p */
	4.4350269		       /* osac_el_h */
     },
     {
	6,
	842.2090,		       /* length_p */
	842.2000,		       /* length_h */
	-445.0821,		       /* osac_z0_p */
	479.2152,		       /* osac_z0_h */
	0.2245,			       /* osac_y0_p */
	-0.2067,		       /* osac_y0_h */
	0.08625,		       /* osac_x0_p */
	-0.1096,		       /* osac_x0_h */
	0.0,			       /* osac_p_p */
	-4.9625995845653374E-04,       /* osac_p_h */
	-2.4957050467401789E+00,       /* osac_k_p */
	-7.2620248152618760E+00,       /* osac_k_h */
	3.2056977725634789E+02,	       /* osac_r_p */
	3.0609851668776219E+02,	       /* osac_r_h */
	0.0,			       /* osac_az_p */
	2.3720568,		       /* osac_az_h */
	0.0,			       /* osac_el_p */
	4.4891913		       /* osac_el_h */
     }
};

/*}}}*/

/*{{{ Static variables and HRMA Parms */

static char *HRMA_Opt_File;

static double HRMA_Vignetting_Factor = 0.9;

double _Marx_HRMA_Cap_Position = 10079.5;

static int Use_Scale_Factors = 1;

static char *Mirror_Shutters [MARX_NUM_MIRRORS];
#if MARX_HAS_WFOLD
static int Use_Wfold_Tables = 0;
#endif
static int Use_Blur_Factors = 1;
static char *Geometry_File;
static int HRMA_Is_Ideal = 0;

#if MARX_HRMA_HAS_STRUTS
static int HRMA_Use_Struts = 1;
#endif
static Param_Table_Type HRMA_Parm_Table [] =
{
   {"HRMAOptConst",	PF_FILE_TYPE,		&HRMA_Opt_File},
   {"HRMAVig",		PF_REAL_TYPE,		&HRMA_Vignetting_Factor},
   {"HRMA_Cap_X",	PF_REAL_TYPE,		&_Marx_HRMA_Cap_Position},
   {"Shutters1", 	PF_STRING_TYPE,		&Mirror_Shutters[0]},
   {"Shutters3", 	PF_STRING_TYPE,		&Mirror_Shutters[1]},
   {"Shutters4", 	PF_STRING_TYPE,		&Mirror_Shutters[2]},
   {"Shutters6", 	PF_STRING_TYPE,		&Mirror_Shutters[3]},
   {"HRMA_P1H1_XOffset",PF_REAL_TYPE,		&HRMA_Mirrors[0].xoffset},
   {"HRMA_P3H3_XOffset",PF_REAL_TYPE,		&HRMA_Mirrors[1].xoffset},
   {"HRMA_P4H4_XOffset",PF_REAL_TYPE,		&HRMA_Mirrors[2].xoffset},
   {"HRMA_P6H6_XOffset",PF_REAL_TYPE,		&HRMA_Mirrors[3].xoffset},
   {"P1Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[0].p_blur},
   {"H1Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[0].h_blur},
   {"P3Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[1].p_blur},
   {"H3Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[1].h_blur},
   {"P4Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[2].p_blur},
   {"H4Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[2].h_blur},
   {"P6Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[3].p_blur},
   {"H6Blur",		PF_REAL_TYPE,		&HRMA_Mirrors[3].h_blur},
   {"HRMA_Use_Blur",	PF_BOOLEAN_TYPE,	&Use_Blur_Factors},
   {"HRMA_Use_Scale_Factors",	PF_BOOLEAN_TYPE,	&Use_Scale_Factors},
#if MARX_HAS_WFOLD
   {"HRMA_Use_WFold",	PF_BOOLEAN_TYPE,	&Use_Wfold_Tables},
   {"H1ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[0].h_scat_factor},
   {"P1ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[0].p_scat_factor},
   {"H3ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[1].h_scat_factor},
   {"P3ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[1].p_scat_factor},
   {"H4ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[2].h_scat_factor},
   {"P4ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[2].p_scat_factor},
   {"H6ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[3].h_scat_factor},
   {"P6ScatFactor", 	PF_REAL_TYPE,		&HRMA_Mirrors[3].p_scat_factor},
#endif
   {"HRMA_Ideal",	PF_BOOLEAN_TYPE,	&HRMA_Is_Ideal},
#if MARX_HRMA_HAS_STRUTS
   {"HRMA_Use_Struts",	PF_BOOLEAN_TYPE,	&HRMA_Use_Struts},
#endif

   {"HRMA_Geometry_File",PF_FILE_TYPE,		&Geometry_File},

   {NULL, 0, NULL}
};

/* Arrays of optical constants */
static float *Betas;
static float *Deltas;
static float *Energies;
static unsigned int Num_Energies;

/*}}}*/

static int get_rdb_value (Marx_RDB_File_Type *rdb, int row,
			  char *colname, double *value)
{
   int c;
   char *s;

   c = marx_rdb_get_col (rdb, colname);
   if (c == -1)
     return -1;

   s = marx_rdb_get_value (rdb, (unsigned int) row, (unsigned int) c);
   if (s == NULL)
     return -1;

   if (1 != sscanf (s, "%lf", value))
     {
	marx_error ("Error parsing field %s as a double", colname);
	return -1;
     }

   return 0;
}

static int read_hrma_geometry_file (char *file)
{
   Marx_RDB_File_Type *rdb;
   HRMA_Type *h, *hmax;

   marx_message ("\t%s\n", file);

   if (NULL == (rdb = marx_open_rdb_file (file)))
     return -1;

   h = HRMA_Mirrors;
   hmax = HRMA_Mirrors + MARX_NUM_MIRRORS;

   while (h < hmax)
     {
	int r;
	char mirror[3];

	/* Parabola */
	sprintf (mirror, "p%u", h->mirror_number);

	r = marx_rdb_get_row (rdb, "mirror", mirror);
	if (r == -1)
	  {
	     marx_close_rdb_file (rdb);
	     return -1;
	  }

	if ((-1 == get_rdb_value (rdb, r, "x0", &h->osac_x0_p))
	    || (-1 == get_rdb_value (rdb, r, "y0", &h->osac_y0_p))
	    || (-1 == get_rdb_value (rdb, r, "z0", &h->osac_z0_p))
	    || (-1 == get_rdb_value (rdb, r, "p", &h->osac_p_p))
	    || (-1 == get_rdb_value (rdb, r, "k", &h->osac_k_p))
	    || (-1 == get_rdb_value (rdb, r, "rho0", &h->osac_r_p))
	    || (-1 == get_rdb_value (rdb, r, "az_mis", &h->osac_az_p))
	    || (-1 == get_rdb_value (rdb, r, "el_mis", &h->osac_el_p))
	    || (-1 == get_rdb_value (rdb, r, "l", &h->length_p)))
	  {
	     marx_close_rdb_file (rdb);
	     return -1;
	  }

	/* hyperbola */

	sprintf (mirror, "h%u", h->mirror_number);

	r = marx_rdb_get_row (rdb, "mirror", mirror);
	if (r == -1)
	  {
	     marx_close_rdb_file (rdb);
	     return -1;
	  }

	if ((-1 == get_rdb_value (rdb, r, "x0", &h->osac_x0_h))
	    || (-1 == get_rdb_value (rdb, r, "y0", &h->osac_y0_h))
	    || (-1 == get_rdb_value (rdb, r, "z0", &h->osac_z0_h))
	    || (-1 == get_rdb_value (rdb, r, "p", &h->osac_p_h))
	    || (-1 == get_rdb_value (rdb, r, "k", &h->osac_k_h))
	    || (-1 == get_rdb_value (rdb, r, "rho0", &h->osac_r_h))
	    || (-1 == get_rdb_value (rdb, r, "az_mis", &h->osac_az_h))
	    || (-1 == get_rdb_value (rdb, r, "el_mis", &h->osac_el_h))
	    || (-1 == get_rdb_value (rdb, r, "l", &h->length_h)))
	  {
	     marx_close_rdb_file (rdb);
	     return -1;
	  }

	h++;
     }

   marx_close_rdb_file (rdb);
   return 0;
}

/*{{{ conic section routines */

static void blur_normal (JDMVector_Type *, double, double);

/* The conic is given by x^2 + r^2 = a x^2 + b x + c */
static double compute_conic_radius (double a, double b, double c, double x)
{
   return sqrt (c + x * (b + x * (a - 1)));
}

/* This function computes the intersection of a ray (x0, p) with the
 * portion of the surface
 *   x^2 + y^2 + z^2 = a x^2 + b x + c
 * that lies between planes x=xmin and x=xmax.
 * If the ray intersects, 0 is returned as well as the
 * intersection point (x0) and normal (via parameter list).  If no
 * intersection, -1 is returned.
 */
static int compute_conic_intersection (double a, double b, double c,
				       JDMVector_Type *x0, JDMVector_Type p,
				       JDMVector_Type *normal,
				       double xmin, double xmax)
{
   double alpha, beta, gamma, x_y, x_z;
   double t_plus, t_minus, x_plus, x_minus;

   /* project ray to x = 0 plane */
   t_plus = -x0->x / p.x;
   x_y = x0->y + t_plus * p.y;
   x_z = x0->z + t_plus * p.z;

   alpha = a * p.x * p.x - 1.0;
   beta = b * p.x - 2.0 * (p.y * x_y + p.z * x_z);
   gamma = c - x_z * x_z - x_y * x_y;

   if (alpha == 0.0)
     {
	if (beta == 0.0) return -1;
	/* beta t + gamma = 0 */
	t_plus = t_minus = -gamma / beta;
     }
   else
     if (0 >= JDMquadratic_root (alpha, beta, gamma, &t_plus, &t_minus))
       return -1;

   /* Now find out what x coordinate the ts correspond to. */
   x_plus = p.x * t_plus;
   x_minus = p.x * t_minus;

   if ((x_plus >= xmin) && (x_plus < xmax))
     {
	/* x_plus looks good.  Check x_minus. If ok, choose greatest */
	if ((x_minus >= xmin) && (x_minus < xmax)
	    && (x_minus > x_plus))
	  {
	     x0->x = x_minus;
	     x0->y = x_y + p.y * t_minus;
	     x0->z = x_z + p.z * t_minus;
	  }
	else
	  {
	     x0->x = x_plus;
	     x0->y = x_y + p.y * t_plus;
	     x0->z = x_z + p.z * t_plus;
	  }
     }
   else if ((x_minus >= xmin) && (x_minus < xmax))
     {
	x0->x = x_minus;
	x0->y = x_y + p.y * t_minus;
	x0->z = x_z + p.z * t_minus;
     }
   else
     {
	/* Out of range. */
	return -1;
     }

   /* Now compute inward normal */
   normal->x = (a - 1) * x0->x + 0.5 * b;
   normal->y = -x0->y;
   normal->z = -x0->z;
   JDMv_normalize (normal);

   return 0;
}

/* Conic has equation
 * x^2 + y^2 + z^2 = a x^2 + b x + c
 *
 * If ray missed the conic, return -2.
 * If ray absorbed, return -1.
 * Otherwise return 0.
 */

static int reflect_from_conic (double a, double b, double c,
			       JDMVector_Type *x, JDMVector_Type *p,
			       double xmin, double xmax,
#if MARX_HAS_WFOLD
			       Marx_WFold_Table_Type *wfold,
			       double scattering_factor,
#endif
			       double blur,
			       double energy, double beta, double delta,
			       double correction_factor)
{
   JDMVector_Type normal;
   double p_dot_n;
   double r, rfl;
#if MARX_HAS_WFOLD
   double sin_grazing;
   double delta_grazing;
#endif
   if (-1 == compute_conic_intersection (a, b, c,
					 x, *p, &normal,
					 xmin, xmax))
     return -2;

   if (Use_Blur_Factors)
     blur_normal (&normal, blur, energy);

   p_dot_n = JDMv_pdot_prod (p, &normal);

   if (HRMA_Is_Ideal == 0)
     {
	r = JDMrandom ();
	rfl = marx_reflectivity (fabs(p_dot_n), beta, delta);
	if (r >= rfl * correction_factor)
	  return -1;
     }

   /* If p is normaized, then this transformation will keep it normalized. */
   *p = JDMv_pax1_bx2 (1.0, p, -2.0 * p_dot_n, &normal);

#if MARX_HAS_WFOLD
   if (Use_Wfold_Tables == 0)
     return 0;

   sin_grazing = -p_dot_n;
   r = JDMrandom ();
   delta_grazing = marx_wfold_table_interp (wfold, energy, sin_grazing, r);
   delta_grazing *= scattering_factor;

   if (delta_grazing > PI/4) return -1;

   if (JDMrandom () < 0.5) delta_grazing = -delta_grazing;

   *p = JDMv_rotate_unit_vector (*p,
				 JDMv_pcross_prod (p, &normal),
				 delta_grazing);
#endif
   return 0;
}

/*}}}*/

/*{{{ init_hrma_shells */

static int get_shell_geometry (Param_File_Type *pf, int shell) /*{{{*/
{
   HRMA_Type *h;
   int id;
   double cap_pos;

   (void) pf;
   h = HRMA_Mirrors + shell;
   id = h->mirror_number;

   cap_pos = _Marx_HRMA_Cap_Position + h->xoffset;

   h->conic_c_p = h->osac_r_p * h->osac_r_p;
   h->conic_b_p = -2.0 * h->osac_k_p;
   h->conic_a_p = 1.0 - h->osac_p_p;

   h->conic_xmin_p = -0.5 * h->length_p;
   h->conic_xmax_p = 0.5 * h->length_p;

   h->conic_c_h = h->osac_r_h * h->osac_r_h;
   h->conic_b_h = -2.0 * h->osac_k_h;
   h->conic_a_h = 1.0 - h->osac_p_h;

   h->conic_xmin_h = -0.5 * h->length_h;
   h->conic_xmax_h = 0.5 * h->length_h;

   /* Note that a ray whose X coordinate is cap_pos, will have an OSAC coord
    * of -|osac_z0_p| for the parabola, and +|osac_z0_h| for the hyperbola.
    * For reference, the relations between the systems are:
    *
    *    MARX_X = -SAOSAC_Z
    *    MARX_Y = +SAOSAC_X
    *    MARX_Z = -SAOSAC_Y
    *
    * Also, a ray at the OSAC origin will have the MARX coordinate of
    * cap - osac_z0_p.
    */
   h->to_osac_p.x = -(-h->osac_z0_p);
   h->to_osac_p.y =  (-h->osac_x0_p);
   h->to_osac_p.z = -(-h->osac_y0_p);

   h->to_osac_p.x -= cap_pos;

   h->to_osac_h.x = -(-h->osac_z0_h);
   h->to_osac_h.y =  (-h->osac_x0_h);
   h->to_osac_h.z = -(-h->osac_y0_h);

   h->to_osac_h.x -= cap_pos;

   h->front_position = cap_pos - h->osac_z0_p + h->conic_xmax_p;
   return 0;
}

/*}}}*/
#if MARX_HAS_WFOLD
static Marx_WFold_Table_Type *get_fold_table (Param_File_Type *pf, /*{{{*/
					      char *fmt, int id)
{
   char pname[80];
   char *file, filebuf[PF_MAX_LINE_LEN];
   Marx_WFold_Table_Type *t;

   sprintf (pname, fmt, id);
   if (-1 == pf_get_file (pf, pname, filebuf, sizeof (filebuf)))
     return NULL;

   if (NULL == (file = marx_make_data_file_name (filebuf)))
     return NULL;

   marx_message ("\t%s\n", file);

   t = marx_read_wfold_file (file);

   marx_free (file);

   return t;
}

/*}}}*/

static int get_shell_wfold_tables (Param_File_Type *pf, unsigned int shell) /*{{{*/
{
   HRMA_Type *h;
   int id;

   h = HRMA_Mirrors + shell;
   id = h->mirror_number;

   if (h->p_wfold_table != NULL)
     {
	marx_free_wfold_table (h->p_wfold_table);
	h->p_wfold_table = NULL;
     }

   if (h->h_wfold_table != NULL)
     {
	marx_free_wfold_table (h->h_wfold_table);
	h->h_wfold_table = NULL;
     }

   if (NULL == (h->p_wfold_table = get_fold_table (pf, "WFold_P%d_File", id)))
     return -1;

   if (NULL == (h->h_wfold_table = get_fold_table (pf, "WFold_H%d_File", id)))
     {
	marx_free_wfold_table (h->p_wfold_table);
	h->p_wfold_table = NULL;
	return -1;
     }

   return 0;
}

/*}}}*/
#endif				       /* HAS_WFOLD */

static int init_hrma_shells (Param_File_Type *pf) /*{{{*/
{
   unsigned int i;
   HRMA_Type *h;
   double total_area;

   /* Note: the area_fraction is really the cumulative area fraction. */

#if MARX_HAS_WFOLD
   if (Use_Wfold_Tables)
     marx_message ("Reading scattering tables\n");
#endif

   total_area = 0.0;
   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	if (-1 == get_shell_geometry (pf, i))
	  return -1;

	if (Use_Wfold_Tables)
	  {
	     if (-1 == get_shell_wfold_tables (pf, i))
	       return -1;
	  }
	h = HRMA_Mirrors + i;

	/* For now use back of parabola as min radius and front as max */
#if 1
	h->min_radius = compute_conic_radius (h->conic_a_p, h->conic_b_p, h->conic_c_p,
					      h->conic_xmin_p);
#else
	h->min_radius = compute_conic_radius (h->conic_a_h, h->conic_b_h, h->conic_c_h,
					      h->conic_xmin_h);
#endif
	h->max_radius = compute_conic_radius (h->conic_a_p, h->conic_b_p, h->conic_c_p,
					      h->conic_xmax_p);

	/* Tweak the min radius to allow off-axis angles to hit all parts
	 * of the conic.  Assume that the maximum off-axis angle of interest
	 * is 40'.  Use: dr/length=tan(theta)
	 */
	h->min_radius -= h->length_p * tan (0.67 * PI/180);

	if (h->num_open_shutters)
	  {
	     total_area += (h->num_open_shutters / 4.0) *
	       (h->max_radius - h->min_radius) * (h->max_radius + h->min_radius);
	  }
	h->area_fraction = total_area;
     }

   if (total_area <= 0.0)
     {
	marx_error ("The mirror geometric area is 0.  Check shutters.");
	return -1;
     }

   /* This is used in the Flux function to get the time between photons.
    * cm^2
    */
   Marx_Mirror_Geometric_Area = (total_area * PI) / 100.0;

   /* Now normalize the cumulative area fraction. */
   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	h = HRMA_Mirrors + i;
	h->area_fraction = h->area_fraction / total_area;
     }
   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ optical constant routines */

static void free_optical_constants (void)
{
   if (NULL != Betas) JDMfree_float_vector (Betas);
   if (NULL != Deltas) JDMfree_float_vector (Deltas);
   if (NULL != Energies) JDMfree_float_vector (Energies);

   Betas = Deltas = Energies = NULL;
   Num_Energies = 0;
}

static int read_hrma_opt_constants (void)
{
   unsigned int nread;
   char *file;

   free_optical_constants ();

   file = HRMA_Opt_File;
   if ((file == NULL) || (*file == 0))
     return 0;

   if (NULL == (file = marx_make_data_file_name (file)))
     return -1;

   /* The optical constant file consists of:
    *   energy (KeV), beta, delta
    */
   marx_message ("Reading binary HRMA optical constants:\n\t%s\n", file);

   if (-1 == marx_f_read_bdat (file, &nread, 3, &Energies, &Betas, &Deltas))
     {
	marx_free (file);
	return -1;
     }

   marx_free (file);

   Num_Energies = nread;
   return 0;
}

static void free_correction_factors (HRMA_Type *h)
{
   marx_free ((char *) h->correction_factors);
   marx_free ((char *) h->correction_energies);
   h->correction_energies = NULL;
   h->correction_factors = NULL;
   h->num_correction_factors = 0;
}

static int read_hrma_correction_factors (void)
{
   char *file;
   char filebuf[32];
   HRMA_Type *h, *hmax;

   h = HRMA_Mirrors;
   hmax = HRMA_Mirrors + MARX_NUM_MIRRORS;

   while (h < hmax)
     {
	free_correction_factors (h);

	sprintf (filebuf, "hrma/corr_%d.dat", h->mirror_number);
	if (NULL == (file = marx_make_data_file_name (filebuf)))
	  return -1;

	marx_message ("\t%s\n", file);

	if (-1 == marx_f_read_bdat (file, &h->num_correction_factors, 2,
				    &h->correction_energies,
				    &h->correction_factors))
	  {
	     marx_free (file);
	     return -1;
	  }
	marx_free (file);

	h++;
     }
   return 0;
}

/*}}}*/

int _Marx_Gratings_Locked_To_Hrma = 0;

#if MARX_HAS_HRMA_PITCH_YAW
static int init_yaw_pitch (void)
{
   HRMA_Type *h, *hmax;

   h = HRMA_Mirrors;
   hmax = HRMA_Mirrors + MARX_NUM_MIRRORS;

   while (h < hmax)
     {
	JDMVector_Type new_y_axis;
	JDM_3Matrix_Type rz, ry;
	double az, el;

	/* According to the report 'XRCF Phase 1 Testing: Preliminary Results',
	 * from May 28 1997, osac_az specifies a positive rotation about the
	 * SAOSAC +Y axis.  This axis is the -Z marx axis.  Then osac_el
	 * is a NEGATIVE rotation about the SAOSAC X' axis.  The X' axis is
	 * the new X axis after the osac_az rotation.  The SAOSAC X axis
	 * corresponds to the MARX Y axis.  Thus:
	 *
	 * 1. Rotate by osac_az about the -Z Marx axis.  Or, equivalently,
	 *    rotate by -osac_az about the +Z axis.
	 *
	 * 2. Rotate about the new Y Marx axis by -osac_el.
	 */

	el = -h->osac_el_p * (PI/180.0);   /* - */
	az = -h->osac_az_p * (PI/180); /* - */

	JDM3m_rot_z_matrix (rz, -az);
	new_y_axis = JDM3m_vector_mul (rz, JDMv_vector (0, 1, 0));
	JDM3m_rot_matrix (ry, new_y_axis, -el);
	JDM3m_mul (h->fwd_matrix_p, ry, rz);
	/* backward */
	JDM3m_rot_matrix (ry, new_y_axis, el);
	JDM3m_rot_z_matrix (rz, az);
	JDM3m_mul (h->bwd_matrix_p, rz, ry);

	/* hyperboloid matrices */
	el = -h->osac_el_h * (PI/180.0);
	az = -h->osac_az_h * (PI/180);

	JDM3m_rot_z_matrix (rz, -az);
	new_y_axis = JDM3m_vector_mul (rz, JDMv_vector (0, 1, 0));
	JDM3m_rot_matrix (ry, new_y_axis, -el);
	JDM3m_mul (h->fwd_matrix_h, ry, rz);
	/* backward */
	JDM3m_rot_matrix (ry, new_y_axis, el);
	JDM3m_rot_z_matrix (rz, az);
	JDM3m_mul (h->bwd_matrix_h, rz, ry);

	h++;
     }

   return 0;
}
#endif

#if MARX_HRMA_HAS_STRUTS
typedef struct
{
   double xpos;			       /* offset from cap */
   double half_width;
}
HRMA_Strut_Type;

#define NUM_PRECOLIMATOR_STRUTS	2
static HRMA_Strut_Type Precolimator_Struts [NUM_PRECOLIMATOR_STRUTS] =
{
   {1492.060,	0.5*0.5*25.4},
   {942.266,	0.5*0.5*25.4}	       /* FAP strut */
};
#define NUM_CAP_STRUTS	2
static HRMA_Strut_Type Cap_Struts [NUM_CAP_STRUTS] =
{
   {0.5*1.965*25.4, 0.5*0.75*25.4},
   {-0.5*1.965*25.4, 0.5*0.75*25.4}
};

#define NUM_POSTCOLIMATOR_STRUTS 2
static HRMA_Strut_Type Postcolimator_Struts [NUM_POSTCOLIMATOR_STRUTS] =
{
   {-1050.353, 0.5*0.5*25.4},
   {-1271.333, 0.5*0.5*25.4}
};

static int intersects_struts (Marx_Photon_Attr_Type *at, HRMA_Strut_Type *s, unsigned int num)
{
   double theta = 30.0*(PI/180.0);
   double cos_theta = cos(theta);
   double sin_theta = sin(theta);
   HRMA_Strut_Type *smax = s + num;
   unsigned int i, num_rotations = 3;

   while (s < smax)
     {
	double half_width = s->half_width;
	double x, y, z;
	double px, py, pz;
	double t;

	x = at->x.x; y = at->x.y; z = at->x.z;
	px = at->p.x; py = at->p.y; pz = at->p.z;

	t = (s->xpos + _Marx_HRMA_Cap_Position - x)/px;
	y += py*t; z += pz*t;

	for (i = 0; i < num_rotations; i++)
	  {
	     if (i != 0)
	       {
		  double tmp = cos_theta*y - sin_theta*z;
		  z = sin_theta*y + cos_theta*z;
		  y = tmp;
	       }

	     if (((-half_width < y) && (y < half_width))
		 || ((-half_width < z) && (z < half_width)))
	       {
		  at->flags |= PHOTON_MIRROR_VBLOCKED;
		  return 1;
	       }
	  }
	s++;
     }
   return 0;
}

static int intersects_precolimator_struts (Marx_Photon_Attr_Type *at)
{
   return intersects_struts (at, Precolimator_Struts, NUM_PRECOLIMATOR_STRUTS);
}
static int intersects_cap_struts (Marx_Photon_Attr_Type *at)
{
   return intersects_struts (at, Cap_Struts, NUM_CAP_STRUTS);
}
static int intersects_postcolimator_struts (Marx_Photon_Attr_Type *at)
{
   return intersects_struts (at, Postcolimator_Struts, NUM_POSTCOLIMATOR_STRUTS);
}
#endif

static int project_photon_to_hrma (Marx_Photon_Attr_Type *at, /*{{{*/
				   double source_distance)
{
   double r;
   unsigned int i;
   HRMA_Type *h;

   while (1)
     {
	r = JDMrandom ();
	for (i = 0; i < MARX_NUM_MIRRORS; i++)
	  {
	     h = HRMA_Mirrors + i;

	     if (r < h->area_fraction)
	       {
		  double theta, radius;
		  unsigned int bitmap, quad;

		  at->mirror_shell = i;

		  radius = h->min_radius
		    + (h->max_radius - h->min_radius) * JDMrandom ();

		  bitmap = h->shutter_bitmap;

		  do
		    {
		       theta = JDMrandom ();
		       quad = (unsigned int) (4.0 * theta);   /* 0, 1, 2, 3 */
		    }
		  while (0 == (bitmap & (1 << quad)));

		  theta = (2.0 * PI) * (theta - 1.0 / 8.0);

		  at->x.z = radius * cos (theta);
		  at->x.y = radius * sin (theta);
		  at->x.x = h->front_position;

		  /* Account for mis-alignment of this shell */
		  at->x.z -= h->to_osac_p.z;
		  at->x.y -= h->to_osac_p.y;

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

#if MARX_HRMA_HAS_STRUTS
		  if (HRMA_Use_Struts
		      && (1 == intersects_precolimator_struts (at)))
		    return -1;
#endif
		  return 0;
	       }
	  }
     }
}

/*}}}*/

static void blur_normal (JDMVector_Type *n, double blur, double energy) /*{{{*/
{
   double phi;
   double n_y, n_z, len;
   JDMVector_Type perp;

   (void) energy;
   /* Choose a random axis perp to p and rotate by a gaussian distributed angle
    * about it.  This is accompished in two steps:
    *
    *   1.  Pick a vector perp to p.  Now rotate it by a random angle about p.
    *   2.  Rotate the normal by gaussian distributed angle about this
    *        random axis.
    */

   /* Step 1.
    * The vector n will never have both y and z components zero.  So
    * a vect perp to it is:
    */
   n_y = n->y;
   n_z = n->z;
   len = sqrt (n_y * n_y + n_z * n_z);

   perp.x = 0.0;
   perp.y = n_z / len;
   perp.z = -n_y / len;

   /* Now rotate this about n. */
   phi = (2.0 * PI) * JDMrandom ();
   perp = JDMv_rotate_unit_vector (perp, *n, phi);

   /* Step 2.
    * Rotate n about perp by gaussian distributed angle.
    */

   phi = blur * (1.0 / 3600.0 *  PI / 180.0);    /* 1 arc sec */
   phi = phi * JDMgaussian_random ();

   *n = JDMv_rotate_unit_vector (*n, perp, phi);
}

/*}}}*/

int _marx_hrma_mirror_init (Param_File_Type *p) /*{{{*/
{
   unsigned int j;
   char *file;

   if (-1 == pf_get_parameters (p, HRMA_Parm_Table))
     return -1;

   if (HRMA_Is_Ideal)
     {
	Use_Blur_Factors = 0;
	Use_Wfold_Tables = 0;
     }

   if (NULL == (file = marx_make_data_file_name (Geometry_File)))
     return -1;

   if (-1 == read_hrma_geometry_file (file))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   for (j = 0; j < MARX_NUM_MIRRORS; j++)
     {
	if (-1 == _marx_parse_shutter_string (Mirror_Shutters[j],
					      &HRMA_Mirrors[j].shutter_bitmap,
					      &HRMA_Mirrors[j].num_open_shutters))
	  return -1;
     }

   if (HRMA_Is_Ideal == 0)
     {
	if (-1 == read_hrma_opt_constants ())
	  return -1;

	if (Use_Scale_Factors)
	  {
	     if (-1 == read_hrma_correction_factors ())
	       return -1;
	  }
     }

   if (-1 == init_hrma_shells (p))
     return -1;

#if MARX_HAS_HRMA_PITCH_YAW
   /* This function needs shell information to be correct.  For that reason
    * it must follow init_hrma_shells.
    */
   if (-1 == init_yaw_pitch ())
     return -1;
#endif

   return 0;
}

/*}}}*/

int _marx_hrma_mirror_reflect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   double *photon_energies;
   unsigned int n, i, *sorted_index;
   double last_energy;
   double beta = 0.0, delta = 0.0, correction_factor = 1.0;
   double source_distance;
   HRMA_Type *last_h;

   if (pt->history & MARX_MIRROR_SHELL_OK)
     return 0;			       /* been here already */
   pt->history |= MARX_MIRROR_SHELL_OK;

   marx_prune_photons (pt);
   n = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   /* First of all, apply vignetting factor to kill a certain percentage
    * of photons.
    */
   if (HRMA_Is_Ideal == 0)
     {
	for (i = 0; i < n; i++)
	  {
	     at = photon_attributes + sorted_index[i];
	     if (JDMrandom () > HRMA_Vignetting_Factor)
	       {
		  at->flags |= PHOTON_MIRROR_VBLOCKED;
	       }
	  }

	/* I could have pruned in the previous loop but it is a better idea to
	 * leave it for a function call.
	 */
	marx_prune_photons (pt);
	photon_energies = pt->sorted_energies;
	n = pt->num_sorted;
     }

   /* source_distance is in mm */
   source_distance = pt->source_distance;

   last_energy = -1.0;
   last_h = NULL;
   for (i = 0; i < n; i++)
     {
	HRMA_Type *h;
	double energy;
	int status;

	at = photon_attributes + sorted_index[i];
	project_photon_to_hrma (at, source_distance);

	h = HRMA_Mirrors + at->mirror_shell;

	/* The conic intersection/reflection routines are expressed in a
	 * coordinate system whose origin is at the center of the conic.
	 * For that reason, we need to move our ray to that position.
	 */
	at->x = JDMv_sum (at->x, h->to_osac_p);

#if MARX_HAS_HRMA_PITCH_YAW
	at->x = JDM3m_vector_mul (h->fwd_matrix_p, at->x);
	/* We only consider the change in direction induced by the
	 * HRMA orientation and not the change in the photon position.
	 * The reason for this is that except for very
	 * near sources, the incoming rays will be parallel and the only real
	 * effect of the pitch/yaw would be to reduce the effective area
	 * by an extremely small amount.
	 *
	 * Unfortunately, this assumption is not valid after reflection.
	 */
	at->p = JDM3m_vector_mul (h->fwd_matrix_p, at->p);
#endif

	energy = at->energy;

	if (Num_Energies == 0)
	  {
	     beta = 0.0;
	     delta = 1.0;	       /* perfect reflect */
	     correction_factor = 1.0;
	  }
	else
	  {
	     if (energy != last_energy)
	       {
		  beta = JDMinterpolate_f (energy, Energies, Betas, Num_Energies);
		  delta= JDMinterpolate_f (energy, Energies, Deltas, Num_Energies);
	       }

	     if (Use_Scale_Factors)
	       {
		  if ((energy != last_energy) || (h != last_h))
		    {
		       correction_factor = JDMinterpolate_f (energy, h->correction_energies, h->correction_factors, h->num_correction_factors);
		       correction_factor = sqrt (correction_factor);
		    }
	       }
	  }

	last_energy = energy;
	last_h = h;

	status = reflect_from_conic (h->conic_a_p, h->conic_b_p, h->conic_c_p,
				     &at->x, &at->p,
				     h->conic_xmin_p, h->conic_xmax_p,
#if MARX_HAS_WFOLD
				     h->p_wfold_table,
				     h->p_scat_factor,
#endif
				     h->p_blur,
				     energy, beta, delta, correction_factor);
	if (status == -1)
	  {
	     at->flags |= PHOTON_UNREFLECTED;
	     continue;
	  }
	if (status == -2)
	  {
	     at->flags |= PHOTON_UNREFLECTED;
	     continue;
	  }

#if MARX_HAS_HRMA_PITCH_YAW
	/* The backward transformation is more complicated since the position
	 * must be taken into account.
	 */
	at->p = JDM3m_vector_mul (h->bwd_matrix_p, at->p);
	at->x = JDM3m_vector_mul (h->bwd_matrix_p, at->x);
#endif

	/* Now go back to our coordinate system and then into OSAC for hyperbola */
	at->x = JDMv_diff (at->x, h->to_osac_p);

#if MARX_HRMA_HAS_STRUTS
	if (HRMA_Use_Struts
	    && (1 == intersects_cap_struts (at)))
	  continue;
#endif
	at->x = JDMv_sum (at->x, h->to_osac_h);

#if MARX_HAS_HRMA_PITCH_YAW
	/* Now rotate into the coordinate system of the hyperbola */

	at->p = JDM3m_vector_mul (h->fwd_matrix_h, at->p);
	at->x = JDM3m_vector_mul (h->fwd_matrix_h, at->x);
#endif

	if (0 != reflect_from_conic (h->conic_a_h, h->conic_b_h, h->conic_c_h,
				     &at->x, &at->p,
				     h->conic_xmin_h, h->conic_xmax_h,
#if MARX_HAS_WFOLD
				     h->h_wfold_table,
				     h->h_scat_factor,
#endif
				     h->h_blur,
				     energy, beta, delta, correction_factor))
	  {
	     at->flags |= PHOTON_UNREFLECTED;
	     continue;
	  }

#if MARX_HAS_HRMA_PITCH_YAW
	at->p = JDM3m_vector_mul (h->bwd_matrix_h, at->p);
	at->x = JDM3m_vector_mul (h->bwd_matrix_h, at->x);
#endif

	at->x = JDMv_diff (at->x, h->to_osac_h);

#if MARX_HRMA_HAS_STRUTS
	if (HRMA_Use_Struts
	    && (1 == intersects_postcolimator_struts (at)))
	  continue;
#endif

     }
   return 0;
}

/*}}}*/

