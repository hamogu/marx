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

/*{{{ #include files */

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

#define USE_GEFF_CALDB_FILE	1
/*{{{ Grating Structure */

/* Geometric parameters, etc... */
typedef struct
{
   unsigned int num_refs;	       /* number of objects referencing this */
   unsigned int num_orders;	       /* 2n+1 */
   int *order_list;		       /* malloced */
   /* The energies pointer will only be malloced if the energies are
    * read in from a file.  Otherwise, it will point to space allocated to
    * the global variable Energies.
    */
   float *energies;
   unsigned int num_energies;

   float **cum_efficiencies;	       /* malloced quantity.  In MARX, these
					* are cumulative efficiencies.
					* num_orders x num_energies
					*/
   double dispersion_angle;
   double period;
   double dp_over_p;
   double theta_blur;
   double vig;			       /* vignetting */

   int flags;
#define GRATING_EFF_CORRECTED		0x1
#define GRATING_DTHETA_CORRECTED	0x2
}
Grating_Type;

/*}}}*/

typedef struct
{
   unsigned int num_sectors;
   double *min_angle;
   double *max_angle;
   double *dtheta;		       /* rotation of facet from nominal */
   double *dtheta_blur;		       /* rotation of facet from nominal */
   double *dp_over_p;		       /* deviation from nominal period */
   double *dp_over_p_blur;	       /* deviation from nominal period */
}
Grating_Sector_Type;
/*{{{ static global variables */

static Grating_Sector_Type *Grating_Sectors [MARX_NUM_MIRRORS];
static Grating_Sector_Type *read_sector_info (char *);
static int Use_Letg_Sector_Files;
static int Use_Hetg_Sector_Files;
static int Use_Unit_Efficiencies = 0;
static int Use_This_Order = 0;

static int Use_File_Efficiencies = 0;
/* if 0, perform computation of efficiencies.  Otherwise, read efficiencies
 * from a file.
 */

/* These variables are used when marx computes the grating efficiencies. */
static float *Energies;
static unsigned int Num_Energies = 1024;
static double Min_Energy = 0.03;
static double Max_Energy = 10.0;
static char *Optical_Constants;
/* Energy dependent optical constants.  There are Num_Energies long and
 * correspond to the Energies array.
 */
static float *Beta_Factor;
static float *Delta_Factor;
static float *Poly_Factor;
static float *Gold_Factor;
static float *Nickel_Factor;
static float *Chromium_Factor;

/*}}}*/

/*{{{ free_grating */

static void free_grating (Grating_Type *g)
{
   if (g == NULL) return;

   if (g->num_refs > 1)
     {
	g->num_refs -= 1;
	return;
     }

   if (g->cum_efficiencies != NULL)
     {
	JDMfree_float_matrix (g->cum_efficiencies, g->num_orders);
	g->cum_efficiencies = NULL;
     }
   if (g->order_list != NULL)
     {
	JDMfree_integer_vector (g->order_list);
	g->order_list = NULL;
     }
   if (g->energies != NULL)
     {
	if (g->energies != Energies)
	  JDMfree_float_vector (g->energies);
	g->energies = NULL;
     }
   marx_free ((char *) g);
}

/*}}}*/
/*{{{ free_grating_factors */

static void free_grating_factors (void)
{
   if (Beta_Factor != NULL) JDMfree_float_vector (Beta_Factor);
   if (Delta_Factor != NULL) JDMfree_float_vector (Delta_Factor);
   if (Poly_Factor != NULL) JDMfree_float_vector (Poly_Factor);
   if (Gold_Factor != NULL) JDMfree_float_vector (Gold_Factor);
   if (Chromium_Factor != NULL) JDMfree_float_vector (Chromium_Factor);
   if (Nickel_Factor != NULL) JDMfree_float_vector (Nickel_Factor);

   Beta_Factor = Delta_Factor = Poly_Factor = NULL;
   Gold_Factor = Chromium_Factor = Nickel_Factor = NULL;
}

/*}}}*/
/*{{{ read_grating_optical_consts */

static int read_grating_optical_consts (char *optfile)
{
   float *data[7];
   float *new_data[7 - 1];
   unsigned int nread;
   unsigned int i;

   if (NULL == optfile) return -1;

   free_grating_factors ();

   if ((NULL == (Beta_Factor = JDMfloat_vector (Num_Energies)))
       || (NULL == (Delta_Factor = JDMfloat_vector (Num_Energies)))
       || (NULL == (Gold_Factor = JDMfloat_vector (Num_Energies)))
       || (NULL == (Nickel_Factor = JDMfloat_vector (Num_Energies)))
       || (NULL == (Chromium_Factor = JDMfloat_vector (Num_Energies)))
       || (NULL == (Poly_Factor = JDMfloat_vector (Num_Energies))))
     {
	free_grating_factors ();
	return -1;
     }

   /* The optical constant file consists of:
    *   energy, beta, delta, t-gold, t-poly, t-cr, t-ni
    * The t-?? correspond to 0.001 thickness of matter.
    */
   marx_message ("\t%s\n", optfile);

   if (-1 == marx_f_read_bdat (optfile, &nread, 7,
			       &data[0],
			       &data[1],
			       &data[2],
			       &data[3],
			       &data[4],
			       &data[5],
			       &data[6]))
     {
	free_grating_factors ();
	return -1;
     }

   /* interpolate to the new grid */
   new_data[0] = Beta_Factor;
   new_data[1] = Delta_Factor;
   new_data[2] = Gold_Factor;
   new_data[3] = Poly_Factor;
   new_data[4] = Chromium_Factor;
   new_data[5] = Nickel_Factor;

   JDMinterpolate_n_fvector (Energies, new_data, Num_Energies,
			     data[0], &data[1], (unsigned int) nread, 7 - 1);
   for (i = 0; i < 7; i++)
     {
	JDMfree_float_vector (data[i]);
     }

   /* Compute logs now so we do not have to do it later.  Note the
    * factor of 1000 derives from the fact that the transmission constants
    * in the file correspond to a distance 0f 0.001 microns
    */
   for (i = 0; i < Num_Energies; i++)
     {
	double factor;

	factor = Gold_Factor[i];
	if (factor <= 0.0) factor = 0.0;
	else factor = log (factor);
	Gold_Factor[i] = 1000.0 * factor;

	factor = Poly_Factor[i];
	if (factor <= 0.0) factor = 0.0;
	else factor = log (factor);
	Poly_Factor[i] = 1000.0 * factor;

	factor = Chromium_Factor[i];
	if (factor <= 0.0) factor = 0.0;
	else factor = log (factor);
	Chromium_Factor[i] = 1000.0 * factor;

	factor = Nickel_Factor[i];
	if (factor <= 0.0) factor = 0.0;
	else factor = log (factor);
	Nickel_Factor[i] = 1000.0 * factor;
     }
   return 0;
}

/*}}}*/

/* The grating efficiency function is only used for the LETG support
 * structures as well as the efficiency computation when Use_File_Efficiencies
 * is 0.
 */
/*{{{ compute_efficiency functions */

static double grating_efficiency (int order, double energy,
				  double beta, double delta,
				  double poly, double gold,
				  double chromium, double nickel,
				  Marx_Grating_Info_Type *g)
{
   double tfactor;
   double kh, khb, k;
   double eff;
   double ekhb;
   double a;

   a = 1.0 - g->bar_width / g->period;
   k = energy / HBAR_C;
   kh = g->bar_height * k;
   khb = kh * beta;
   ekhb = exp (-khb);

   if (order != 0)
     {
	double mpi = order * PI;
	a = sin (mpi * a) / mpi;
	eff = a * a * (1 + ekhb * (ekhb -
				   2.0 * cos (kh * delta)));
     }
   else
     {
	double b = 1 - a;
	eff = a * a + b * ekhb * (b * ekhb + 2 * a * cos (kh * delta));
     }

   /* Compute the transmittance--- the logs have been computed earlier */
   tfactor = exp (g->polyimide * poly
		  + g->gold * gold
		  + g->chromium * chromium
		  + g->nickel * nickel);

   return eff * tfactor;
}

static double compute_efficiency (int order, unsigned int nth_energy,
				  Marx_Grating_Info_Type *g)
{
   return grating_efficiency (order,
			      Energies [nth_energy],
			      Beta_Factor [nth_energy],
			      Delta_Factor [nth_energy],
			      Poly_Factor [nth_energy],
			      Gold_Factor [nth_energy],
			      Chromium_Factor [nth_energy],
			      Nickel_Factor [nth_energy],
			      g);
}

/* Here we compute the efficiency of some energy at some order.  What we do
 * is compute two efficiencies on the grid defined by the Energies vector
 * and then interpolate the result.
 */
double marx_compute_grating_efficiency (double energy, int order,
					Marx_Grating_Info_Type *geom)
{
   unsigned int n0, n1;
   double eff0, eff1, e0, e1;

   if (Energies == NULL)
     {
	marx_error ("marx_compute_grating_efficiency: not initialized.");
	return -1.0;
     }

   /* Find where energy lies on the Energies grid. */
   n0 = JDMbinary_search_f (energy, Energies, Num_Energies);
   n1 = n0 + 1;

   if (n1 == Num_Energies)
     {
	/* Oops.  energy is out of range.
	 * Use last two points to extrapolate. */
	n0--;
	n1--;
     }

   e0 = Energies[n0];
   eff0 = compute_efficiency (order, n0, geom);

   if (energy == e0)
     return eff0;

   e1 = Energies[n1];
   eff1 = compute_efficiency (order, n1, geom);

   return eff0 + (eff1 - eff0)/(e1 - e0) * (energy - e0);
}

/*}}}*/

/* This routines only initializes a grating whose efficiencies are computed.
 * A separate routine is used to initialize the grating structure when the
 * efficiencies are read from a file.
 */
static Grating_Type *init_one_grating (Marx_Grating_Info_Type *info, int use_this_order, int use_unit_eff)    /*{{{*/
{
   unsigned int i, order_number, j;
   unsigned int num_orders, num_energies;
   Grating_Type *g;

   if (NULL == (g = (Grating_Type *) marx_malloc (sizeof (Grating_Type))))
     return NULL;
   memset ((char *) g, 0, sizeof (Grating_Type));

   num_orders = info->num_orders;
   if (use_this_order) num_orders = 2;
   g->num_orders = num_orders;

   g->num_energies = num_energies = Num_Energies;
   g->energies = Energies;

   g->dispersion_angle = info->dispersion_angle;
   g->theta_blur = info->theta_blur;
   g->period = info->period;
   g->vig = info->vig;
   g->dp_over_p = info->dp_over_p;

   if (NULL == (g->cum_efficiencies = JDMfloat_matrix (num_orders,
						       num_energies)))
     {
	marx_free ((char *) g);
	return NULL;
     }

   if (NULL == (g->order_list = JDMinteger_vector (num_orders)))
     {
	JDMfree_float_matrix (g->cum_efficiencies, num_orders);
	g->cum_efficiencies = NULL;
	JDMfree_float_vector (g->energies);
	g->energies = NULL;
	marx_free ((char *) g);
	return NULL;
     }

   /* Initialize the order list.  The orders are arranged as 0, +1, -1, +2...
    * This is approximately in the direction of decreasing efficiency.  I do
    * this to minimize the amount of time spent determining the order that
    * a photon will go into.
    */

   if (use_this_order)
     {
	g->order_list[0] = use_this_order;
	g->order_list[1] = -use_this_order;
     }
   else
     {
	g->order_list[0] = 0;
	j = 1;
	for (i = 1; i < num_orders - 1; i += 2)
	  {
	     g->order_list[i] = j;
	     g->order_list[i + 1] = -j;
	     j++;
	  }
     }

   /* Now compute efficiencies */
   for (i = 0; i < num_energies; i++)
     {
	double sum = 0.0;

	for (order_number = 0; order_number < num_orders; order_number++)
	  {
	     if (use_unit_eff)
	       sum += (1.0 / num_orders);
	     else
	       sum += compute_efficiency (g->order_list [order_number],
					  i,
					  info);

	     g->cum_efficiencies[order_number][i] = sum;
	  }
     }

   g->num_refs = 0;

   return g;
}

/*}}}*/

/*{{{ Grating_Type Structure initialization */

static int Sim_Use_LETG = 0;
static double MEG_Rowland_Diameter = 8634.0;	       /* mm */
static double HEG_Rowland_Diameter = 8634.0;	       /* mm */
static double LEG_Rowland_Diameter = 8634.0;	       /* mm */

#define REST_OF_STRUCT_FIELDS 0, 0, 0
static Marx_Grating_Info_Type HEG_Grating_Info =
{
   0.2,				       /* period (um)*/
   0.7,				       /* bar_height (um) */
   0.11,			       /* bar_width */
   1.0,				       /* double polyimide */
   1.0,				       /* gold */
   0.0,				       /* chromium */
   0.0,				       /* nickel */
   0.0,				       /* dp_over_p */
   -5.0,				       /* dispersion angle */
   REST_OF_STRUCT_FIELDS
};

static Marx_Grating_Info_Type MEG_Grating_Info =
{
   0.4,				       /* period (um)*/
   0.4,				       /* bar_height (um) */
   0.22,			       /* bar_width */
   0.5,				       /* double polyimide */
   1.0,				       /* gold*/
   0.0,				       /* chromium */
   0.0,				       /* nickel */
   0.0,				       /* dp_over_p */
   5.0,				       /* dispersion angle */
   REST_OF_STRUCT_FIELDS
};

static Marx_Grating_Info_Type LEG_Grating_Info =
{
   0.9921,			       /* period (um)*/
   0.5,				       /* bar_height (um) */
   0.4712,			       /* bar_width */
   0.0,				       /* double polyimide */
   0.0,				       /* gold*/
   0.0,				       /* chromium */
   0.0,				       /* nickel */
   0.0,				       /* dp_over_p */
   0.0,				       /* dispersion angle */
   REST_OF_STRUCT_FIELDS
};

static Marx_Grating_Info_Type LEG_Fine_Grating_Info =
{
   25.4,			       /* period */
   2.5,				       /* height */
   2.5,				       /* bar width */
   0.0,				       /* double polyimide */
   0.0,				       /* gold*/
   0.0,				       /* chromium */
   0.0,				       /* nickel */
   0.0,				       /* dp_over_p */
   0.0,				       /* dispersion angle */
   REST_OF_STRUCT_FIELDS
};

static Marx_Grating_Info_Type LEG_Coarse_Grating_Info =
{
   2000.0,			       /* period */
   26.0,			       /* height */
   68.0,			       /* bar width */
   0.0,				       /* double polyimide */
   0.0,				       /* gold*/
   0.0,				       /* chromium */
   0.0,				       /* nickel */
   0.0,				       /* dp_over_p */
   0.0,				       /* dispersion angle */
   REST_OF_STRUCT_FIELDS
};

static Grating_Type *Gratings[MARX_NUM_MIRRORS];
static Grating_Type *LEG_Fine_Grating;
static Grating_Type *LEG_Coarse_Grating;

/*}}}*/

/*{{{ Grating Parameter table definition */

static double LEG_Eff_Scale_Factor;
static Param_Table_Type Grating_Parm_Table [] =
{
     {"UseGratingEffFiles",	PF_BOOLEAN_TYPE,&Use_File_Efficiencies},
     {"Use_HETG_Sector_Files", 	PF_BOOLEAN_TYPE,&Use_Hetg_Sector_Files},
     {"Use_LETG_Sector_Files", 	PF_BOOLEAN_TYPE,&Use_Letg_Sector_Files},
     {"Use_Unit_Efficiencies",	PF_BOOLEAN_TYPE,&Use_Unit_Efficiencies},
     {"Use_This_Order",		PF_INTEGER_TYPE,&Use_This_Order},

     {"GratingOptConsts",	PF_FILE_TYPE,	&Optical_Constants},
     {"MEGRowlandDiameter",	PF_REAL_TYPE,	&MEG_Rowland_Diameter},
     {"HEGRowlandDiameter",	PF_REAL_TYPE,	&HEG_Rowland_Diameter},
     {"LEGRowlandDiameter",	PF_REAL_TYPE,	&LEG_Rowland_Diameter},
     {"HEGVig",			PF_REAL_TYPE,	&HEG_Grating_Info.vig},
     {"MEGVig",			PF_REAL_TYPE,	&MEG_Grating_Info.vig},
     {"LEGVig",			PF_REAL_TYPE,	&LEG_Grating_Info.vig},

     {"LETG_Eff_Scale_Factor",	PF_REAL_TYPE,	&LEG_Eff_Scale_Factor},

     {"hegNumOrders",		PF_INTEGER_TYPE,&HEG_Grating_Info.num_orders},
     {"hegGold",		PF_REAL_TYPE,	&HEG_Grating_Info.gold},
     {"hegNickel",		PF_REAL_TYPE,	&HEG_Grating_Info.nickel},
     {"hegChromium",		PF_REAL_TYPE,	&HEG_Grating_Info.chromium},
     {"hegPolyimide",		PF_REAL_TYPE,	&HEG_Grating_Info.polyimide},
     {"hegBarHeight",		PF_REAL_TYPE,	&HEG_Grating_Info.bar_height},
     {"hegBarWidth",		PF_REAL_TYPE,	&HEG_Grating_Info.bar_width},
     {"hegPeriod",		PF_REAL_TYPE,	&HEG_Grating_Info.period},
     {"hegdPoverP",		PF_REAL_TYPE,	&HEG_Grating_Info.dp_over_p},
     {"hegTheta",		PF_REAL_TYPE, 	&HEG_Grating_Info.dispersion_angle},
     {"hegdTheta",		PF_REAL_TYPE, 	&HEG_Grating_Info.theta_blur},

     {"legNumOrders",		PF_INTEGER_TYPE,&LEG_Grating_Info.num_orders},
     {"legGold",		PF_REAL_TYPE,	&LEG_Grating_Info.gold},

     {"legBarHeight",		PF_REAL_TYPE,	&LEG_Grating_Info.bar_height},
     {"legBarWidth",		PF_REAL_TYPE,	&LEG_Grating_Info.bar_width},
     {"legPeriod",		PF_REAL_TYPE,	&LEG_Grating_Info.period},
     {"legdPoverP",		PF_REAL_TYPE,	&LEG_Grating_Info.dp_over_p},
     {"legTheta",		PF_REAL_TYPE, 	&LEG_Grating_Info.dispersion_angle},
     {"legdTheta",		PF_REAL_TYPE, 	&LEG_Grating_Info.theta_blur},

     {"megNumOrders",		PF_INTEGER_TYPE,&MEG_Grating_Info.num_orders},
     {"megGold",		PF_REAL_TYPE,	&MEG_Grating_Info.gold},
     {"megNickel",		PF_REAL_TYPE,	&MEG_Grating_Info.nickel},
     {"megChromium",		PF_REAL_TYPE,	&MEG_Grating_Info.chromium},
     {"megPolyimide",		PF_REAL_TYPE,	&MEG_Grating_Info.polyimide},
     {"megBarHeight",		PF_REAL_TYPE,	&MEG_Grating_Info.bar_height},
     {"megBarWidth",		PF_REAL_TYPE,	&MEG_Grating_Info.bar_width},
     {"megPeriod",		PF_REAL_TYPE,	&MEG_Grating_Info.period},
     {"megdPoverP",		PF_REAL_TYPE,	&MEG_Grating_Info.dp_over_p},
     {"megTheta",		PF_REAL_TYPE, 	&MEG_Grating_Info.dispersion_angle},
     {"megdTheta",		PF_REAL_TYPE, 	&MEG_Grating_Info.theta_blur},

     {"legFineNumOrders",	PF_INTEGER_TYPE,&LEG_Fine_Grating_Info.num_orders},
     {"legCoarseNumOrders",	PF_INTEGER_TYPE,&LEG_Coarse_Grating_Info.num_orders},

     {"legFineBarWidth",	PF_REAL_TYPE,&LEG_Fine_Grating_Info.bar_width},
     {"legFineBarHeight",	PF_REAL_TYPE,&LEG_Fine_Grating_Info.bar_height},
     {"legFinePeriod",		PF_REAL_TYPE,&LEG_Fine_Grating_Info.period},

     {"legCoarseBarWidth",	PF_REAL_TYPE,&LEG_Coarse_Grating_Info.bar_width},
     {"legCoarseBarHeight",	PF_REAL_TYPE,&LEG_Coarse_Grating_Info.bar_height},
     {"legCoarsePeriod",	PF_REAL_TYPE,&LEG_Coarse_Grating_Info.period},

     {NULL, 0, NULL}
};

/*}}}*/

/* The routines here compute the intersection point of the incoming ray
 * with the Rowland torus.  Newton's method is used to solve the quartic.
 */
/*{{{ newtons_quartic, intersect, rotate */
static int newtons_quartic (double a, double b, double c, double d,
			    double t0, double *tp)
{
   unsigned int max_it = 10;
   double a2, a3, b2;
   double t2, t;
   double eps = 1.0e-4;
   double num, den;

   a2 = 2.0 * a;
   a3 = 3.0 * a;
   b2 = 2.0 * b;

   while (1)
     {
	t2 = t0 * t0;
	num = t2 * (3.0 * t2 + a2 * t0 + b) - d;
	den = t2 * (4.0 * t0 + a3) + b2 * t0 + c;
	t = num / den;

	if (fabs (t - t0) < eps) break;

	max_it--;
	if (max_it == 0) return -1;

	t0 = t;
     }

   *tp = t;
   return 0;
}

static int intersect (JDMVector_Type *x0, JDMVector_Type p, double rowland)
{
   double a, b, c, d, t0, t;
   double x2, pdotx;
   double pxpz_len;
   double r2;

   /* Here we assume that x is projected to the x = 0 plane */
   t = -x0->x / p.x;
   x0->x = 0.0;
   x0->y = x0->y + p.y * t;
   x0->z = x0->z + p.z * t;

   pxpz_len = p.x * p.x + p.z * p.z;
   x2 = x0->z * x0->z + x0->y * x0->y;	       /* x.x = 0 */
   r2 = rowland * rowland;
   pdotx = JDMv_pdot_prod (&p, x0);

   a = 4.0 * pdotx;
   b = 2.0 * x2 + a * pdotx - r2 * pxpz_len;
   c = a * x2 - 2.0 * r2 * p.z * x0->z;
   d = x2 * x2 - r2 * x0->z * x0->z;

   t0 = -rowland * sqrt (pxpz_len);

   if (-1 == newtons_quartic (a, b, c, d, t0, &t))
     {
	fprintf (stderr, "Quartic failed.\n");
	return -1;
     }

   *x0 = JDMv_pax1_bx2 (1.0, x0, t, &p);
   return 0;
}

/* Rotate about the x axis by an angle theta */
static JDMVector_Type rotate (JDMVector_Type a, double theta)
{
   double c, s, ay, az;
   c = cos(theta);
   s = sin(theta);
   ay = a.y;
   az = a.z;
   a.y = c * ay - s * az;
   a.z = s * ay + c * az;
   return a;
}

/*}}}*/

/*{{{ diffract_photon */

static int diffract_photon (Grating_Type *g, double theta,
			    Marx_Photon_Attr_Type *at,
			    int order, Grating_Sector_Type *gs)
{
   JDMVector_Type p, x;
   JDMVector_Type n, l, d, dp;
   double factor;
   double n_lambda_over_d, dp_over_p, dtheta;
   double p_n, p_d, p_l;

   n_lambda_over_d = order * (2.0 * PI * HBAR_C) / g->period / at->energy;

   x = at->x;
   p = at->p;

   /* Now find new direction of the ray.  We assume that the facets are
    * perfect in the sense that the facet normal lies in the direction of
    * the origin.  In the future, we will want to get this information from
    * the facet database as a function of x.
    *
    * n, l, and d are the basis vectors of the facet.  n is normal to the
    * facet, l is in the direction of the lines, and d is in the is normal
    * to the lines.
    */
   /* n = -x/|x| */
   n = JDMv_unit_vector (x);
   n.x = -n.x; n.y = -n.y; n.z = -n.z;

   l.x = n.z;
   l.y = 0.0;
   l.z = -n.x;

   l = JDMv_unit_vector (l);
   d = JDMv_cross_prod (n, l);

   if (gs != NULL)
     {
	unsigned int sector_num;
	double sector;
	/* The sector is determined by the y and z coordinates of the ray. */

	sector = atan2 (x.y, x.z);
	if (sector < 0) sector = 2*PI + sector;

	sector_num = JDMbinary_search_d (sector, gs->min_angle, gs->num_sectors);
	if ((gs->max_angle[sector_num] <= sector)
	    || (gs->min_angle[sector_num] > sector))
	  return -1;

	dtheta = gs->dtheta[sector_num]
	  + gs->dtheta_blur[sector_num] * JDMgaussian_random ();

	dp_over_p = gs->dp_over_p[sector_num]
	  + gs->dp_over_p_blur[sector_num] * JDMgaussian_random ();
     }
   else
     {
	/* Until more information is obtained upon a facet by facet basis involving
	 * the misalignment of the grating, apply a statistical blur.
	 */
	dtheta = g->theta_blur * JDMgaussian_random ();
	dp_over_p = g->dp_over_p * JDMgaussian_random ();
     }

   theta -= dtheta;

   /* Now rotate about n axis by theta.  This handles the support structure
    * for the LETG */

   if (theta != 0.0)
     {
	JDMVector_Type l_tmp, d_tmp;
	double c, s;

	c = cos(theta);
	s = sin(theta);
	l_tmp = l;
	d_tmp = d;
	l = JDMv_ax1_bx2 (c, l_tmp, s, d_tmp);
	d = JDMv_ax1_bx2 (-s, l_tmp, c, d_tmp);
     }

   p_d = n_lambda_over_d + JDMv_pdot_prod (&p, &d);
   p_l = JDMv_pdot_prod (&p, &l);
   p_n = 1.0 - p_l * p_l - p_d * p_d;

   if (p_n < 0.0)
     return -1;

   p_n = sqrt (p_n);

   p = JDMv_ax1_bx2 (p_d, d, p_n, n);
   at->p = JDMv_ax1_bx2 (1.0, p, p_l, l);

   /* Now apply dp/p blur */
   if (dp_over_p == 0)
     return 0;

   factor = n_lambda_over_d * dp_over_p;

   dp = JDMv_ax1_bx2 (-factor, d, factor * (p_d / p_n), n);

   at->p.x += dp.x;
   at->p.y += dp.y;
   at->p.z += dp.z;

   at->p = JDMv_unit_vector (at->p);
   return 0;
}

/*}}}*/
/*{{{ diffract_photon_from_grating */

static int diffract_photon_from_grating (Grating_Type *g, double theta,
					 Marx_Photon_Attr_Type *at,
					 float **cum_eff, unsigned int photon_number,
					 SIGNED_CHAR *order,
					 Grating_Sector_Type *gs)
{
   double r;
   unsigned int num_orders, k;

   r = JDMrandom ();
   num_orders = g->num_orders;

   for (k = 0; k < num_orders; k++)
     {
	if (r <= cum_eff[k][photon_number])
	  {
	     *order = g->order_list[k];
	     return diffract_photon (g, theta, at, *order, gs);
	  }
     }
   return -1;
}

/*}}}*/
/*{{{ rotate_photons */

static void rotate_photons (Marx_Photon_Type *pt, int dir)
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int n, i;
   unsigned int *sorted_index;
   Grating_Type *g;

   n = pt->num_sorted;
   sorted_index = pt->sorted_index;
   photon_attributes = pt->attributes;

   for (i = 0; i < n; i++)
     {
	double theta;
	at = photon_attributes + sorted_index[i];
	g = Gratings[at->mirror_shell];
	if (g == NULL) continue;       /* should not happen */
	theta = dir * g->dispersion_angle;
	at->x = rotate (at->x, theta);
	at->p = rotate (at->p, theta);
     }
}

/*}}}*/
/*{{{ intersect_photons */

static void intersect_photons (Marx_Photon_Type *pt)
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int n, i;
   unsigned int *sorted_index;
   double rowland[4];

   n = pt->num_sorted;
   sorted_index = pt->sorted_index;
   photon_attributes = pt->attributes;

   if (Sim_Use_LETG)
     {
	rowland[0] = LEG_Rowland_Diameter;
	rowland[1] = LEG_Rowland_Diameter;
	rowland[2] = LEG_Rowland_Diameter;
	rowland[3] = LEG_Rowland_Diameter;
     }
   else
     {
	rowland[0] = MEG_Rowland_Diameter;
	rowland[1] = MEG_Rowland_Diameter;
	rowland[2] = HEG_Rowland_Diameter;
	rowland[3] = HEG_Rowland_Diameter;
     }

   for (i = 0; i < n; i++)
     {
	at = photon_attributes + sorted_index[i];
	if (-1 == intersect (&(at->x), at->p, rowland[(unsigned int)at->mirror_shell]))
	  at->flags |= PHOTON_UNDIFFRACTED;
     }
}

/*}}}*/
/*{{{ diffract_from_support_grating */

static void diffract_from_support_grating (Grating_Type *g, double theta,
					   Marx_Photon_Type *pt,
					   float *tmp_energies,
					   float **tmp_cum_efficiencies,
					   int support_index)
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int j;
   unsigned int photon_number = 0;
   unsigned int num_orders;
   unsigned int n;
   unsigned int *sorted_index;

   n = pt->num_sorted;
   sorted_index = pt->sorted_index;
   photon_attributes = pt->attributes;

   for (j = 0; j < n; j++)
     {
	at = photon_attributes + sorted_index[j];
	if (at->flags & BAD_PHOTON_MASK) continue;       /* already absorbed */
	tmp_energies[photon_number] = at->energy;
	photon_number++;
     }
   if (photon_number == 0) return;

   num_orders = g->num_orders;

   /* Now perform the interpolation */
   JDMinterpolate_n_fvector (tmp_energies, tmp_cum_efficiencies, photon_number,
			     Energies, g->cum_efficiencies, Num_Energies,
			     num_orders);

   /* and loop through photon list finding diffraction order */

   photon_number = 0;
   for (j = 0; j < n; j++)
     {
	Grating_Sector_Type *gs;

	at = photon_attributes + sorted_index[j];
	if (at->flags & BAD_PHOTON_MASK) continue;

	gs = NULL;

	if (-1 == diffract_photon_from_grating (g, theta, at, tmp_cum_efficiencies,
						photon_number,
						at->support_orders + support_index,
						gs))
	  at->flags |= PHOTON_UNDIFFRACTED;

	photon_number++;
     }
}

/*}}}*/

static int diffract (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   float *tmp_energies, **tmp_cum_efficiencies;
   double vig;
   unsigned int num_sorted, i, *sorted_index;
   unsigned int num_orders;

   if (pt->history & MARX_ORDER_OK)
     return 0;
   pt->history |= MARX_ORDER_OK;

   marx_prune_photons (pt);
   num_sorted = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   /* Apply the vignetting factor to kill a certain percentage
    * of photons.
    */
   for (i = 0; i < num_sorted; i++)
     {
	at = photon_attributes + sorted_index[i];
	vig = Gratings [at->mirror_shell]->vig;

	if (JDMrandom () > vig)
	  at->flags |= PHOTON_MIRROR_VBLOCKED;
     }

   /* I could have pruned in the previous loop but it is a better idea to
    * leave it for a function call.
    */
   marx_prune_photons (pt);
   num_sorted = pt->num_sorted;

   /* Rotate the coord system in a manner that is appropriate for the
    * gratings, e.g, MEG -5 degrees, etc...  Then intersect with the torus.
    * At end of routine, we rotate back.
    */
   rotate_photons (pt, -1);
   intersect_photons (pt);

   tmp_energies = JDMfloat_vector (num_sorted);
   if (tmp_energies == NULL)
     return -1;

   /* Loop through all valid photons and diffract off the appropriate
    * grating.  The grating is determined by the Mirror shell the photon
    * diffracted from.
    */
   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	Grating_Type *g;
	unsigned int j;
	unsigned int photon_number;
	Grating_Sector_Type *gs;

	g = Gratings[i];
	if (g == NULL) continue;
	num_orders = g->num_orders;

	tmp_cum_efficiencies = JDMfloat_matrix (num_orders, num_sorted);
	if (tmp_cum_efficiencies == NULL)
	  {
	     JDMfree_float_vector (tmp_energies);
	     return -1;
	  }

	photon_number = 0;
	for (j = 0; j < num_sorted; j++)
	  {
	     at = photon_attributes + sorted_index[j];
	     if ((unsigned int)at->mirror_shell != i) continue;
	     tmp_energies[photon_number] = at->energy;
	     photon_number++;
	  }

	if (photon_number == 0)
	  {
	     JDMfree_float_matrix (tmp_cum_efficiencies, num_orders);
	     continue;
	  }

	/* Now perform the interpolation */
	JDMinterpolate_n_fvector (tmp_energies, tmp_cum_efficiencies, photon_number,
				  g->energies, g->cum_efficiencies, g->num_energies,
				  num_orders);

	/* and loop through photon list finding diffraction order */

	gs = Grating_Sectors[i];
	photon_number = 0;
	for (j = 0; j < num_sorted; j++)
	  {
	     at = photon_attributes + sorted_index[j];
	     if ((unsigned int) at->mirror_shell != i) continue;

	     if (-1 == diffract_photon_from_grating (g, 0.0, at, tmp_cum_efficiencies,
						     photon_number, &at->order, gs))
	       at->flags |= PHOTON_UNDIFFRACTED;

	     photon_number++;
	  }

	JDMfree_float_matrix (tmp_cum_efficiencies, num_orders);
     }

   if (Sim_Use_LETG)
     {
	num_orders = LEG_Fine_Grating->num_orders;
	if (num_orders < LEG_Coarse_Grating->num_orders)
	  num_orders = LEG_Coarse_Grating->num_orders;

	if (num_orders)
	  {
	     tmp_cum_efficiencies = JDMfloat_matrix (num_orders, num_sorted);
	     if (tmp_cum_efficiencies == NULL)
	       {
		  JDMfree_float_vector (tmp_energies);
		  return -1;
	       }
	  }
	else tmp_cum_efficiencies = NULL;/* GCC complains otherwise */

	if (LEG_Fine_Grating->num_orders)
	  {

	     diffract_from_support_grating (LEG_Fine_Grating, PI/2.0, pt,
					    tmp_energies, tmp_cum_efficiencies, 0);
	     pt->history |= MARX_ORDER1_OK;
	  }

	if (LEG_Coarse_Grating->num_orders)
	  {
	     diffract_from_support_grating (LEG_Coarse_Grating, PI/3.0, pt,
					    tmp_energies, tmp_cum_efficiencies, 1);
	     diffract_from_support_grating (LEG_Coarse_Grating, 2.0 * PI/3.0, pt,
					    tmp_energies, tmp_cum_efficiencies, 2);
	     diffract_from_support_grating (LEG_Coarse_Grating, 0.0, pt,
					    tmp_energies, tmp_cum_efficiencies, 3);

	     pt->history |= (MARX_ORDER2_OK
			     | MARX_ORDER3_OK
			     | MARX_ORDER4_OK);
	  }

	if (num_orders)
	  JDMfree_float_matrix (tmp_cum_efficiencies, num_orders);
     }

   JDMfree_float_vector (tmp_energies);

   /* Now rotate photons back. */
   rotate_photons (pt, 1);

   return 0;
}

/*}}}*/

int _marx_letg_diffract (Marx_Photon_Type *pt)
{
   if (Sim_Use_LETG == 0) return -1;
   return diffract (pt);
}

int _marx_hetg_diffract (Marx_Photon_Type *pt)
{
   if (Sim_Use_LETG) return -1;
   return diffract (pt);
}

/*{{{ Initialization routines */

static int get_pfile_parms (Param_File_Type *pfile)
{
   if (pfile == NULL) return -1;
   if (-1 == pf_get_parameters (pfile, Grating_Parm_Table))
     return -1;

   HEG_Grating_Info.dispersion_angle *= PI / 180.0;
   LEG_Grating_Info.dispersion_angle *= PI / 180.0;
   MEG_Grating_Info.dispersion_angle *= PI / 180.0;

   HEG_Grating_Info.theta_blur *= (PI / 180.0) / 60.0;
   LEG_Grating_Info.theta_blur *= (PI / 180.0) / 60.0;
   MEG_Grating_Info.theta_blur *= (PI / 180.0) / 60.0;

   if (Use_Unit_Efficiencies) Use_File_Efficiencies = 0;

   return 0;
}

int marx_create_grating_opt_const_tables (char *file)
{
   unsigned int i;
   double energy, d_energy;

   if (Energies != NULL) JDMfree_float_vector (Energies);

   if (NULL == (Energies = JDMfloat_vector (Num_Energies)))
     return -1;

   d_energy = (Max_Energy - Min_Energy) / Num_Energies;
   energy = Min_Energy;
   for (i = 0; i < Num_Energies; i++)
     {
	Energies[i] = energy;
	energy += d_energy;
     }

   if (-1 == read_grating_optical_consts (file))
     {
	JDMfree_float_vector (Energies);
	Energies = NULL;
	return -1;
     }
   return 0;
}

static int create_opt_const_tables (void)
{
   int ret;
   char *file;

   if (NULL == (file = marx_make_data_file_name (Optical_Constants)))
     return -1;

   ret = marx_create_grating_opt_const_tables (file);
   marx_free (file);
   return ret;
}

static int grating_pre_init (Param_File_Type *pf)
{
   unsigned int i;

   if (-1 == get_pfile_parms (pf))
     return -1;

   if ((Use_File_Efficiencies == 0)
       || Sim_Use_LETG)
     {
	if (-1 == create_opt_const_tables ())
	  return -1;
     }

   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	free_grating (Gratings[i]);
	Gratings[i] = NULL;
     }

   free_grating (LEG_Fine_Grating);
   free_grating (LEG_Coarse_Grating);

   return 0;
}

/* Re-arrange orders from -N .. 0 .. N to 0 +1 -1 ... +N -N */
static int create_order_list (Grating_Type *g)
{
   int *order_list;
   unsigned int num_orders;
   int order, max_order;

   num_orders = g->num_orders;
   if (NULL == (order_list = JDMinteger_vector (num_orders)))
     return -1;

   max_order = (num_orders - 1) / 2;

   for (order = -max_order; order <= max_order; order++)
     order_list[order + max_order] = order;

   g->order_list = order_list;
   return 0;
}

/* Compute cumulative efficiencies from those read in from file. */
static int sum_file_efficiencies (Grating_Type *g, double scale)
{
   float **cum_efficiencies;
   unsigned int num_energies, num_orders;
   unsigned int i, j;
   double sum;

   cum_efficiencies = g->cum_efficiencies;
   num_orders = g->num_orders;
   num_energies = g->num_energies;

   if (Use_This_Order > 0)
     {
	int use_this_order = Use_This_Order;
	for (i = 0; i < num_energies; i++)
	  {
	     sum = 0.0;
	     for (j = 0; j < num_orders; j++)
	       {
		  if (abs(g->order_list[j]) == use_this_order)
		    sum += scale * (double) cum_efficiencies [j][i];
		  cum_efficiencies [j][i] = (float) sum;
	       }
	     if (sum > 0)
	       for (j = 0; j < num_orders; j++)
		 cum_efficiencies [j][i] /= sum;
	  }
	return 0;
     }

   for (i = 0; i < num_energies; i++)
     {
	sum = 0.0;
	for (j = 0; j < num_orders; j++)
	  {
	     sum += scale * (double) cum_efficiencies [j][i];
	     cum_efficiencies [j][i] = (float) sum;
	  }
	if (sum > 1.0)
	  {
	     for (j = 0; j < num_orders; j++)
	       cum_efficiencies [j][i] /= sum;
	  }
     }

   return 0;
}

static int create_efficiency_arrays (unsigned int num_energies,
				     unsigned int num_orders,
				     float ***efficiencies_p,
				     float **energies_p,
				     int *max_order_p)
{
   float **efficiencies;
   float *energies;

   if (num_energies == 0)
     {
	marx_error ("Efficiency file contains no rows");
	return -1;
     }

   if ((num_orders < 3)		       /* -1, 0, 1 */
       || ((num_orders % 2) != 1))
     {
	marx_error ("Grating efficiency file does not contain enough columns or the wrong number");
	return -1;
     }

   /* In this file, the rows of the efficiency matrix correspond to order. */
   efficiencies = JDMfloat_matrix (num_orders, num_energies);
   if (efficiencies == NULL)
     return -1;

   energies = JDMfloat_vector (num_energies);
   if (energies == NULL)
     {
	JDMfree_float_matrix (efficiencies, num_orders);
	return -1;
     }

   *energies_p = energies;
   *efficiencies_p = efficiencies;
   *max_order_p = num_orders / 2;
   return 0;
}

#if USE_GEFF_CALDB_FILE
static int read_geff_caldb_file (char *file, int shell,
				 float ***efficiencies_p,
				 float **energies_p,
				 unsigned int *num_energies_p,
				 int *max_order_p)
{
   int max_order;
   unsigned int num_energies, num_orders = 0;
   float **efficiencies = NULL;
   float *energies = NULL;
   JDFits_Type *f;
   char hduname[16];
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   unsigned int i;

   sprintf (hduname, "AXAF_GREFF%d", shell+1);
   marx_message ("\t%s[%s]\n", file, hduname);

   if (NULL == (f = _marx_open_binary_hdu (file, hduname)))
     return -1;

   r = jdfits_bintable_open_rows (f, 2, "f:ENERGY", "f:EFF");
   if (r == NULL)
     goto return_error;

   num_energies = r->num_rows;
   c = r->col_data;
   if (c[0].repeat != 1)
     {
	marx_error ("The ENERGY column must be a scalar column");
	goto return_error;
     }
   num_orders = c[1].repeat;

   if (-1 == create_efficiency_arrays (num_energies, num_orders,
				       &efficiencies, &energies, &max_order))
     goto return_error;

   for (i = 0; i < num_energies; i++)
     {
	unsigned int j;
	float *eff;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Error reading row %u", i+1);
	     goto return_error;
	  }

	energies[i] = c[0].data.f[0];
	eff = c[1].data.f;
	for (j = 0; j < num_orders; j++)
	  efficiencies[j][i] = eff[j];
     }
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   *efficiencies_p = efficiencies;
   *energies_p = energies;
   *num_energies_p = num_energies;
   *max_order_p = max_order;
   return 0;

   return_error:
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   if (efficiencies != NULL)
     JDMfree_float_matrix (efficiencies, num_orders);
   if (energies != NULL)
     JDMfree_float_vector (energies);
   *efficiencies_p = NULL;
   *energies_p = NULL;

   return -1;
}

#else
static int read_geff_bdat_file (char *file, int shell,
				float ***efficiencies_p,
				float **energies_p,
				unsigned int *num_energies_p,
				int *max_order_p)
{
   int max_order;
   unsigned int num_energies;
   float **efficiencies = NULL;
   float *energies = NULL;
   int order;
   unsigned int num;
   JDMBData_File_Type *bf;
   FILE *fp;

   (void) shell;

   marx_message ("\t%s\n", file);

   bf = JDMbdata_open_file (file);
   if (bf == NULL)
     {
	marx_error ("Unable to open %s", file);
	goto return_error;
     }

   /* The data type for these files should be 'E' float32 */
   if (bf->data_type != 'E')
     {
	marx_error ("%s has wrong data type (found '%c')", file, bf->data_type);
	goto return_error;
     }

   /* We expect columns for the energy, -max_order ... max_order.  That is,
    * 2 max_order + 2 columns.
    */
   num_energies = bf->nrows;

   if (-1 == create_efficiency_arrays (num_energies, bf->ncols-1,
				       &efficiencies, &energies, &max_order))
     goto return_error;

   fp = bf->fp;

   for (num = 0; num < (unsigned int) num_energies; num++)
     {
	if (1 != JDMread_f_float32 (energies + num, 1, fp))
	  {
	     marx_error ("%s: read error", file);
	     goto return_error;
	  }

	for (order = -max_order; order <= max_order; order++)
	  {
	     if (1 != JDMread_f_float32 (&efficiencies[order + max_order][num], 1, fp))
	       {
		  marx_error ("%s: read error", file);
		  goto return_error;
	       }
	  }
     }

   JDMbdata_close_file (bf);
   *efficiencies_p = efficiencies;
   *energies_p = energies;
   *num_energies_p = num_energies;
   *max_order_p = max_order;
   return 0;

   return_error:
   if (bf != NULL)
     JDMbdata_close_file (bf);
   if (efficiencies != NULL)
     JDMfree_float_matrix (efficiencies, max_order * 2 + 1);
   if (energies != NULL)
     JDMfree_float_vector (energies);
   *efficiencies_p = NULL;
   *energies_p = NULL;

   return -1;
}
#endif				       /* USE_GEFF_CALDB_FILE */

static Grating_Type *init_one_file_grating (char *file, int shell, double scale)
{
   int max_order;
   unsigned int num_energies;
   float **efficiencies;
   float *energies;
   Grating_Type *g;
   int status;

   g = (Grating_Type *) marx_malloc (sizeof(Grating_Type));
   if (g == NULL)
     return NULL;

   memset ((char *) g, 0, sizeof (Grating_Type));

#if USE_GEFF_CALDB_FILE
   if (NULL == (file = _marx_caldb_get_file (file)))
     goto return_error;
   status = read_geff_caldb_file (file, shell, &efficiencies, &energies,
				  &num_energies, &max_order);
#else
   if (NULL == (file = marx_make_data_file_name (file)))
     goto return_error;
   status = read_geff_bdat_file (file, shell, &efficiencies, &energies,
				 &num_energies, &max_order);
#endif
   marx_free (file);

   if (status == -1)
     goto return_error;

   g->num_orders = 2 * max_order + 1;
   g->cum_efficiencies = efficiencies;
   g->energies = energies;
   g->num_energies = num_energies;

   if ((0 == create_order_list (g))
       && (0 == sum_file_efficiencies (g, scale)))
   return g;

   /* Get below here only if an error occurs */

   return_error:
   if (g != NULL) free_grating (g);
   return NULL;
}

static int init_file_efficiencies (Param_File_Type *pf, char *prefix_name, double scale)
{
   char *shells = "1346";
   unsigned int i, imax;
   char pname[PF_MAX_LINE_LEN];
   char file[PF_MAX_LINE_LEN];

#if USE_GEFF_CALDB_FILE
   sprintf (file, "%sEFF", prefix_name);
#endif

   imax = strlen (shells);

   if (imax != MARX_NUM_MIRRORS)
     {
	marx_error ("Application error.  init_file_efficiencies: MARX_NUM_MIRRORS mismatch");
	return -1;
     }

   for (i = 0; i < imax; i++)
     {
	unsigned int shell;
	double vig;
	double theta;
	double dtheta;
	double period;
	double dp_over_p;

	Grating_Type *g;

	shell = shells[i] - '0';
#if !USE_GEFF_CALDB_FILE
	sprintf (pname, "%s_Shell%d_File", prefix_name, shell);
	if (-1 == pf_get_file (pf, pname, file, sizeof (file)))
	  return -1;
#endif
	sprintf (pname, "%s_Shell%d_Vig", prefix_name, shell);
	if (-1 == pf_get_double (pf, pname, &vig))
	  return -1;

	sprintf (pname, "%s_Shell%d_Theta", prefix_name, shell);
	if (-1 == pf_get_double (pf, pname, &theta))
	  return -1;

	sprintf (pname, "%s_Shell%d_dTheta", prefix_name, shell);
	if (-1 == pf_get_double (pf, pname, &dtheta))
	  return -1;

	sprintf (pname, "%s_Shell%d_Period", prefix_name, shell);
	if (-1 == pf_get_double (pf, pname, &period))
	  return -1;

	sprintf (pname, "%s_Shell%d_dPoverP", prefix_name, shell);
	if (-1 == pf_get_double (pf, pname, &dp_over_p))
	  return -1;

	g = init_one_file_grating (file, i, scale);

	if (g == NULL)
	  return -1;

	Gratings [i] = g;
	g->dispersion_angle = (theta * (PI / 180.0));
	g->theta_blur = (dtheta * (PI / 180.0) / 60.0);
	g->period = period;
	g->dp_over_p = dp_over_p;
	g->vig = vig;
     }

   return 0;
}

static void free_sector_info (void)
{
   unsigned int i;

   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	Grating_Sector_Type *gs = Grating_Sectors[i];
	if (gs == NULL) continue;

	marx_free ((char *)gs->min_angle);
	marx_free ((char *)gs->max_angle);
	marx_free ((char *)gs->dtheta);
	marx_free ((char *)gs->dp_over_p);
	marx_free ((char *) gs);

	Grating_Sectors[i] = NULL;
     }
}

static int read_sector_files (Param_File_Type *pf, char *name)
{
   char *shells = "1346";
   char filename[1024];
   char *file;
   unsigned int i;

   free_sector_info ();

   i = 0;
   while (shells[i] != 0)
     {
	Grating_Sector_Type *gs;
	char pname[PF_MAX_LINE_LEN];

	sprintf (pname, "%s_Sector%c_File", name, shells[i]);
	if (-1 == pf_get_file (pf, pname, filename, sizeof (filename)))
	  return -1;

	if (NULL == (file = marx_make_data_file_name (filename)))
	  {
	     free_sector_info ();
	     return -1;
	  }

	marx_message ("\t%s\n", file);

	if (NULL == (gs = read_sector_info (file)))
	  {
	     marx_free (file);
	     free_sector_info ();
	     return -1;
	  }

	marx_free (file);
	Grating_Sectors [i] = gs;
	i++;
     }

   return 0;
}

int _marx_hetg_init (Param_File_Type *pf)
{
   Grating_Type *heg, *meg;

   Sim_Use_LETG = 0;

   marx_message ("Initializing HETG...\n");

   if (-1 == grating_pre_init (pf))
     return -1;

   if (Use_Hetg_Sector_Files)
     {
	if (-1 == read_sector_files (pf, "HETG"))
	  return -1;
     }

   if (Use_File_Efficiencies)
     {
	return init_file_efficiencies (pf, "HETG", 1.0);
     }

   if (NULL == (heg = init_one_grating (&HEG_Grating_Info, Use_This_Order, Use_Unit_Efficiencies)))
     return -1;

   if (NULL == (meg = init_one_grating (&MEG_Grating_Info, Use_This_Order, Use_Unit_Efficiencies)))
     {
	free_grating (heg);
	return -1;
     }

   Gratings [0] = Gratings [1] = meg;
   meg->num_refs = 2;

   Gratings [2] = Gratings [3] = heg;
   heg->num_refs = 2;

   return 0;
}

int _marx_letg_init (Param_File_Type *pf)
{
   Grating_Type *g;
   unsigned int i;
   char buf[PF_MAX_LINE_LEN];

   Sim_Use_LETG = 1;

   marx_message ("Initializing LETG...\n");

   if (-1 == grating_pre_init (pf))
     return -1;

   if (Use_Letg_Sector_Files)
     {
	if (-1 == read_sector_files (pf, "LETG"))
	  return -1;
     }

   if (Use_File_Efficiencies)
     {
	LEG_Eff_Scale_Factor = 1;
	if (-1 == init_file_efficiencies (pf, "LETG", LEG_Eff_Scale_Factor))
	  return -1;
     }
   else
     {
	if (NULL == (g = init_one_grating (&LEG_Grating_Info, Use_This_Order, Use_Unit_Efficiencies)))
	  return -1;

	for (i = 0; i < MARX_NUM_MIRRORS; i++)
	  {
	     Gratings [i] = g;
	     g->num_refs += 1;
	  }
     }

   if (-1 == pf_get_string (pf, "DetectorType", buf, sizeof (buf)))
     return -1;
   if (0 == strncmp (buf, "ACIS", 4))
     {
	double delta_theta;
	/* Ugly hack to account for aspect-focal plane misalignmnent problem */
	if (-1 == pf_get_double (pf, "LETG_ACIS_dTheta", &delta_theta))
	  return -1;
	delta_theta *= PI/180.0;
	for (i = 0; i < MARX_NUM_MIRRORS; i++)
	  {
	     if (Gratings[i]->flags & GRATING_DTHETA_CORRECTED)
	       continue;
	     Gratings[i]->dispersion_angle += delta_theta;
	     Gratings[i]->flags |= GRATING_DTHETA_CORRECTED;
	  }
     }

   if ((LEG_Fine_Grating_Info.num_orders)
       && (NULL == (LEG_Fine_Grating = init_one_grating (&LEG_Fine_Grating_Info, 0, 0))))
     return -1;

   if ((LEG_Coarse_Grating_Info.num_orders)
       && (NULL == (LEG_Coarse_Grating = init_one_grating (&LEG_Coarse_Grating_Info, 0, 0))))
     return -1;

   if ((Use_File_Efficiencies == 0) && (Use_Unit_Efficiencies == 0))
     return 0;

   /* Now correct the LEG primary efficiencies to account for loss by the
    * support orders.  The idea here is that the file-based efficiencies represents
    * the efficiency for the primary order without regard for support structures.
    * In other words, let P(1,*) be the tabulated efficiency for first order.
    * The "*" here means any support order.  Then:
    *   P(1,*) = P(1) \sum_i p_i
    * where P(1) is efficiency of the primary grating, and p_i are the efficiencies
    * of the others.  Note that \sum_i p_i may not add to 1 to indicate absorbtion.
    * Hence, we need to use
    *   P(1) = P(1,*)/\sum_i p_i
    */

   /* We diffract from the fine grating once, and the coarse gratings 3 times. */
   for (i = 0; i < MARX_NUM_MIRRORS; i++)
     {
	unsigned int num_energies, num_orders, j, k;
	float **cum_efficiencies, *energies;

	g = Gratings[i];

	if (g->flags & GRATING_EFF_CORRECTED)
	  continue;
	g->flags |= GRATING_EFF_CORRECTED;

	num_orders = g->num_orders;
	num_energies = g->num_energies;
	energies = g->energies;
	cum_efficiencies = g->cum_efficiencies;

	if (LEG_Fine_Grating->num_orders)
	  {
	     float *support_effs, *support_ens;
	     unsigned num_support_ens;

	     support_effs = LEG_Fine_Grating->cum_efficiencies[LEG_Fine_Grating->num_orders-1];
	     support_ens = LEG_Fine_Grating->energies;
	     num_support_ens = LEG_Fine_Grating->num_energies;

	     for (j = 0; j < num_energies; j++)
	       {
		  float factor = JDMinterpolate_f (energies[j], support_ens, support_effs, num_support_ens);
		  factor = 1.0/factor;
		  if (factor <= 0)
		    continue;
		  for (k = 0; k < num_orders; k++)
		    cum_efficiencies[k][j] *= factor;
	       }
	  }
	if (LEG_Coarse_Grating->num_orders)
	  {
	     float *support_effs, *support_ens;
	     unsigned num_support_ens;

	     support_effs = LEG_Coarse_Grating->cum_efficiencies[LEG_Coarse_Grating->num_orders-1];
	     support_ens = LEG_Coarse_Grating->energies;
	     num_support_ens = LEG_Coarse_Grating->num_energies;

	     for (j = 0; j < num_energies; j++)
	       {
		  float factor = JDMinterpolate_f (energies[j], support_ens, support_effs, num_support_ens);
		  factor = factor*factor*factor;   /* 3 supp gratings */
		  if (factor <= 0)
		    continue;
		  factor = 1.0/factor;
		  for (k = 0; k < num_orders; k++)
		    cum_efficiencies[k][j] *= factor;
	       }
	  }
     }
   return 0;
}

/*}}}*/

/*{{{ marx_set_grating_table_parms */
int marx_set_grating_table_parms (double min_energy, double max_energy, unsigned int num_energies)
{
   if (Energies != NULL)
     return -1;

   if (min_energy <= 0.0)
     return -1;
   if (max_energy >= 10.0)
     return -1;
   if (min_energy >= max_energy)
     return -1;
   if (num_energies < 2)
     return -1;

   Min_Energy = min_energy;
   Max_Energy = max_energy;
   Num_Energies = num_energies;
   return 0;
}

/*}}}*/

/*
 * The sector data base file is assumed to be as ASCII file with columns of
 * the form:
 *
 *    MIN_ANGLE MAX_ANGLE DELTA_THETA DELTA_THETA_BLUR DP_OVER_P DP_OVER_P_BLUR...
 */
#define NUM_SECTOR_COLUMNS 6
static Grating_Sector_Type *read_sector_info (char *file)
{
   double *data[NUM_SECTOR_COLUMNS];
   int cindex[NUM_SECTOR_COLUMNS];
   unsigned int num_read;
   Grating_Sector_Type *gs;
   unsigned int i;

   for (i = 0; i < NUM_SECTOR_COLUMNS; i++)
     {
	cindex[i] = 1;
	data[i] = NULL;
     }

   if (-1 == JDMread_column_ddata (file, data, cindex, NUM_SECTOR_COLUMNS, &num_read))
     return NULL;

   if (NULL == (gs = (Grating_Sector_Type *) marx_malloc (sizeof (Grating_Sector_Type))))
     goto return_error;

   memset ((char *) gs, 0, sizeof (Grating_Sector_Type));

   gs->num_sectors = num_read;

   gs->min_angle = data[0];
   gs->max_angle = data[1];
   gs->dtheta = data[2];
   gs->dtheta_blur = data[3];
   gs->dp_over_p = data[4];
   gs->dp_over_p_blur = data[5];

   for (i = 0; i < num_read; i++)
     {
	/* Convert from degrees to radians */
	gs->min_angle[i] *= PI/180.0;
	gs->max_angle[i] *= PI/180.0;

	/* dtheta is in ARC-minutes */
	/* dtheta_blur is in ARC-minutes */
	gs->dtheta[i] *= PI/(180.0 * 60.0);
	gs->dtheta_blur[i] *= PI/(180.0 * 60.0);
     }

   return gs;

   return_error:

   for (i = 0; i < NUM_SECTOR_COLUMNS; i++)
     marx_free ((char *) data[i]);

   return NULL;
}

