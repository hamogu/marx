/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2011-2013 Massachusetts Institute of Technology

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

/*{{{ Wolter1_Type structure and initialization */

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

   /* These are computed */
   /* The equation of the conic is rearranged in form
    * x^2 + y^2 + z^2 = a x^2 + b x + c
    * (0,0,0) is at the center of the conic.  This is the origin of the OSAC
    * coordinate system.
    */
   double conic_a_p;		       /* osac_r_p ^2 */
   double conic_b_p;		       /* -2.0 * osac_k_p */
   double conic_c_p;		       /* 1.0 - osac_p_p */
   double conic_a_h;		       /* osac_r_h ^2 */
   double conic_b_h;		       /* -2.0 * osac_k_h */
   double conic_c_h;		       /* 1.0 - osac_p_h */
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

#if MARX_HAS_HRMA_PITCH_YAW
   JDM_3Matrix_Type fwd_matrix_p;
   JDM_3Matrix_Type bwd_matrix_p;
   JDM_3Matrix_Type fwd_matrix_h;
   JDM_3Matrix_Type bwd_matrix_h;
#endif
}
Wolter1_Type;

typedef struct
{
   unsigned int num_shells;
   Wolter1_Type *shell_info;
}
IXO_Mirror_Type;

static IXO_Mirror_Type *IXO_Mirrors;

/*}}}*/

/*{{{ Static variables and HRMA Parms */

static char *IXO_Opt_File;

static double IXO_Vignetting_Factor = 0.9;
static double IXO_Cap_Position = 20000.0;    /* mm */

static char *Geometry_File;
static int Mirror_Is_Ideal = 0;
static int Use_Blur_Factors = 1;
static double IXO_Mirror_Blur = 0.0;
static double Az_Blur_Sigma = 0.0;
static double El_Blur_Sigma = 0.0;
static double Lateral_Disp_Y_Blur_Sigma = 0;   /* marx Y, osac X */
static double Lateral_Disp_Z_Blur_Sigma = 0;   /* marx Z, osac Y */
static double Defocus_Blur_Sigma = 0;
static char *Mirror_Sector_String;
static Marx_Range_Type *Mirror_Sector_Ranges;
static double Sector_Fraction;

static Param_Table_Type IXOMirror_Parm_Table [] =
{
   {"IXOMirrorOptConst",	PF_FILE_TYPE,		&IXO_Opt_File},
   {"IXOMirrorVig",		PF_REAL_TYPE,		&IXO_Vignetting_Factor},
   {"IXOMirror_Cap_X",		PF_REAL_TYPE,		&IXO_Cap_Position},
   {"IXOMirror_Use_Blur",	PF_BOOLEAN_TYPE,	&Use_Blur_Factors},
   {"IXOMirror_Blur",		PF_REAL_TYPE,		&IXO_Mirror_Blur},
   {"IXOMirror_Ideal",		PF_BOOLEAN_TYPE,	&Mirror_Is_Ideal},
   {"IXOMirror_Geometry_File",	PF_FILE_TYPE,		&Geometry_File},
   {"IXOMirror_Az_Blur",	PF_REAL_TYPE,		&Az_Blur_Sigma},
   {"IXOMirror_El_Blur",	PF_REAL_TYPE,		&El_Blur_Sigma},
   {"IXOMirror_Lateral_Y_Blur",	PF_REAL_TYPE,		&Lateral_Disp_Y_Blur_Sigma},
   {"IXOMirror_Lateral_Z_Blur",	PF_REAL_TYPE,		&Lateral_Disp_Z_Blur_Sigma},
   {"IXOMirror_Defocus_Blur",	PF_REAL_TYPE,		&Defocus_Blur_Sigma},
   {"IXOMirror_Sectors",	PF_STRING_TYPE,		&Mirror_Sector_String},
   {NULL, 0, NULL}
};

/* Arrays of optical constants */
static float *Betas;
static float *Deltas;
static float *Energies;
static unsigned int Num_Energies;

/*}}}*/

static void deallocate_ixo_mirror (IXO_Mirror_Type *m)
{
   if (m == NULL)
     return;
   if (m->shell_info == NULL)
     marx_free ((char *) m->shell_info);
   marx_free ((char *)m);
}

static IXO_Mirror_Type *allocate_ixo_mirror (unsigned int num_shells)
{
   IXO_Mirror_Type *m;

   if (NULL == (m = (IXO_Mirror_Type *)marx_malloc (sizeof(IXO_Mirror_Type))))
     return NULL;

   m->num_shells = num_shells;
   if (NULL == (m->shell_info = (Wolter1_Type *)marx_calloc (num_shells, sizeof(Wolter1_Type))))
     {
	marx_free ((char *)m);
	return NULL;
     }
   return m;
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
			       double blur,
			       double energy, double beta, double delta,
			       double correction_factor)
{
   JDMVector_Type normal;
   double p_dot_n;
   double r, rfl;

   if (-1 == compute_conic_intersection (a, b, c,
					 x, *p, &normal,
					 xmin, xmax))
     return -2;

   if (Use_Blur_Factors)
     blur_normal (&normal, blur, energy);

   p_dot_n = JDMv_pdot_prod (p, &normal);

   if (Mirror_Is_Ideal == 0)
     {
	r = JDMrandom ();
	rfl = marx_reflectivity (fabs(p_dot_n), beta, delta);
	if (r >= rfl * correction_factor)
	  return -1;
     }

   /* If p is normaized, then this transformation will keep it normalized. */
   *p = JDMv_pax1_bx2 (1.0, p, -2.0 * p_dot_n, &normal);

   return 0;
}

/*}}}*/

/*{{{ init_hrma_shells */

static int get_shell_geometry (Param_File_Type *pf, int shell) /*{{{*/
{
   Wolter1_Type *h;
   int id;
   double cap_pos;

   (void) pf;
   h = IXO_Mirrors->shell_info + shell;
   id = h->mirror_number;

   cap_pos = IXO_Cap_Position + h->xoffset;

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

static int init_mirror_shells (Param_File_Type *pf) /*{{{*/
{
   unsigned int i;
   Wolter1_Type *h, *shell_info;
   double total_area;
   unsigned int num_shells;

   /* Note: the area_fraction is really the cumulative area fraction. */

   num_shells = IXO_Mirrors->num_shells;
   shell_info = IXO_Mirrors->shell_info;

   total_area = 0.0;
   for (i = 0; i < num_shells; i++)
     {
	if (-1 == get_shell_geometry (pf, i))
	  return -1;

	h = shell_info + i;

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
	/* h->min_radius -= h->length_p * tan (0.67 * PI/180); */
	total_area += (h->max_radius - h->min_radius) * (h->max_radius + h->min_radius);

	h->area_fraction = total_area;
     }

   if (total_area <= 0.0)
     {
	marx_error ("The mirror geometric area <= 0.");
	return -1;
     }

   /* This is used in the Flux function to get the time between photons.
    * Need to multiply by PI since it was omitted above and divide by 100
    * to get cm^2.  The multiply by the fraction of open sectors.
    */
   Marx_Mirror_Geometric_Area = Sector_Fraction * (total_area * PI) / 100.0;

   /* Now normalize the cumulative area fraction. */
   for (i = 0; i < num_shells; i++)
     {
	h = shell_info + i;
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

static int read_opt_constants (void)
{
   unsigned int nread;
   char *file;

   free_optical_constants ();

   file = IXO_Opt_File;
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

/*}}}*/

#if MARX_HAS_HRMA_PITCH_YAW
static int init_yaw_pitch (void)
{
   Wolter1_Type *h, *hmax;

   h = IXO_Mirrors->shell_info;
   hmax = h + IXO_Mirrors->num_shells;

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

static int project_photon_to_mirror (Marx_Photon_Attr_Type *at, /*{{{*/
				     double source_distance)
{
   double r;
   unsigned int i;
   Wolter1_Type *shell_info;
   unsigned int num_shells;

   num_shells = IXO_Mirrors->num_shells;
   shell_info = IXO_Mirrors->shell_info;

   /* FIXME: Use a lookup table!!! */
   while (1)
     {
	r = JDMrandom ();
	for (i = 0; i < num_shells; i++)
	  {
	     Wolter1_Type *h = shell_info + i;

	     if (r < h->area_fraction)
	       {
		  double theta, radius;
		  double dr, r0, r1;

		  at->mirror_shell = i;

		  /* Use radius^2 = r0^2 + (r1^2-r0^2)*r */
		  r0 = h->min_radius; r1 = h->max_radius; dr = r1-r0;
		  radius = r0*sqrt(1.0 + ((r0+r1)/(r0*r0))*(r1-r0)*JDMrandom());

		  theta = 360.0*JDMrandom();
		  if (Mirror_Sector_Ranges != NULL) while (1)
		    {
		       theta = JDMrandom () * 360.0;
		       if (_marx_is_in_range (Mirror_Sector_Ranges, theta))
			 break;
		    }
		  theta *= PI/180.0;

		  at->x.x = h->front_position;
		  at->x.y = radius * cos (theta);
		  at->x.z = radius * sin (theta);

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

   blur = blur * (1.0 / 3600.0 *  PI / 180.0);    /* 1 arc sec */
#if 1
   phi = blur * JDMgaussian_random ();
#else
   do
     {
	phi = JDMrandom ();
     }
   while (phi == 0.0);

   phi = blur * sqrt (-log (phi));
#endif
   *n = JDMv_rotate_unit_vector (*n, perp, phi);
}

/*}}}*/

static int read_mirror_geometry_file (char *file)
{
   JDFits_Row_Type *r;
   JDFits_Type *f;
   unsigned int num_cols = 8;
   static char *columns[8] =
     {
	"i:shell", "d:d", "d:e", "d:a",
	"d:zmin", "d:zmax", "d:h_zmin", "d:h_zmax"
     };
#define GEOM_COLUMN_shell	0
#define GEOM_COLUMN_d		1
#define GEOM_COLUMN_e		2
#define GEOM_COLUMN_a		3
#define GEOM_COLUMN_zmin	4
#define GEOM_COLUMN_zmax	5
#define GEOM_COLUMN_h_zmin	6
#define GEOM_COLUMN_h_zmax	7

   char *extname = "IXO_MIRROR_GEOM";
   unsigned int i, num_rows;
   Wolter1_Type *shell_info;
   IXO_Mirror_Type *mirror;

   if (IXO_Mirrors != NULL)
     {
	deallocate_ixo_mirror (IXO_Mirrors);
	IXO_Mirrors = NULL;
     }

   marx_message ("Opening IXO Mirror Geometry fits file %s\n", file);
   if (NULL == (f = jdfits_open_binary_table (file, extname)))
     {
	marx_error ("Unable to find a binary table called %s in %s\n",
		    extname, file);
	return -1;
     }

   if (NULL == (r = jdfits_bintable_aopen_rows (f, num_cols, columns)))
     {
	marx_error ("Unable to open a %s as an IXO Mirror Geometry file\n", file);
	jdfits_close_file (f);
	return -1;
     }

   num_rows = r->num_rows;

   if (NULL == (mirror = allocate_ixo_mirror (num_rows)))
     goto return_error;
   shell_info = mirror->shell_info;

   for (i = 0; i < num_rows; i++)
     {
	JDFits_Col_Data_Type *c;
	double d, e, a, zmin, zmax, h_zmin, h_zmax;
	double z0;
	int shell;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Unexpected end of IXO Mirror Geom table %s", file);
	     goto return_error;
	  }

	c = r->col_data;
	shell =c[GEOM_COLUMN_shell].data.i[0];
	d = c[GEOM_COLUMN_d].data.d[0];
	e = c[GEOM_COLUMN_e].data.d[0];
	a = c[GEOM_COLUMN_a].data.d[0];
	zmin = c[GEOM_COLUMN_zmin].data.d[0];
	zmax = c[GEOM_COLUMN_zmax].data.d[0];
	h_zmin = c[GEOM_COLUMN_h_zmin].data.d[0];
	h_zmax = c[GEOM_COLUMN_h_zmax].data.d[0];

	shell_info->mirror_number = shell;
	shell_info->length_p = zmax - zmin;
	shell_info->length_h = h_zmax - h_zmin;

	/* For the paraboloid, The IXO mirror prescription uses:
	 *   r^2 = (d+a+z)^2 - (a+z)^2
	 * Where z is measured from the focus.
	 * In contrast, SAOSAC uses
	 *   r^2 = r0^2 + 2*k*z' - p*z'^2
	 * where z' is measured from the body center of the optic and increases
	 * towards the narrow end of the optic, i.e., to the focus.  Let
	 * z0=(zmin+zmax)/2 in the IXO system, which corresponds to z'=0 in
	 * the SAOSAC system.  Then z = z0 - z'.
	 * ==>
	 *   r^2 = (d+a+z0-z')^2 - (a+z0-z')^2
	 *       = d^2 + 2*d*(a+z0-z')
	 *       = d^2 + 2*d*(a+z0) - 2*d*z'
	 * ==> r0^2 = d*(d+2*(a+z0)), k = -d, p = 0
	 */
	z0 = 0.5*(zmin+zmax);
	shell_info->osac_p_p = 0.0;
	shell_info->osac_k_p = -d;
	shell_info->osac_r_p = sqrt(d*(d+2.0*(a+z0)));
	shell_info->osac_z0_p = -(z0 - IXO_Cap_Position);
	shell_info->osac_x0_p = 0.0;
	shell_info->osac_y0_p = 0.0;

	/* For the hyperbolid, the IXO mirror prescription uses:
	 *  r^2 = e^2*(d+z)^2 - z^2
	 *      = e^2*(d+z0-z')^2 - (z0-z')^2
	 *      = e^2*(d+z0)^2 - 2*e^2*(d+z0)*z' + e^2*z'^2 - (z0^2-2*z0*z'+z'^2)
	 *      = e^2*(d+z0)^2 - z0^2 - 2*z'*(e^2*(d+z0)-z0) - (1-e^2)z'^2
	 * ==> r0^2 = e^2*(d+z0)^2 - z0^2
	 *          = e^2*(d^2 + 2*d*z0 + z0^2) - z0^2
	 *          = e^2*d*(d+2*z0) + z0^2*(e^2-1)
	 *        k = z0 - e^2*(d+z0) = z0*(1-e^2) - e^2*d
	 *        p = (1-e^2)
	 */
	z0 = 0.5*(h_zmin + h_zmax);
	shell_info->osac_p_h = (1.0+e)*(1.0-e);
	shell_info->osac_k_h = z0*(1+e)*(1-e) - e*e*d;
	shell_info->osac_r_h = sqrt(e*e*d*(d+2.0*z0) + z0*z0*(e-1)*(e+1));
	shell_info->osac_z0_h = -(z0 - IXO_Cap_Position);
	shell_info->osac_x0_h = 0.0;
	shell_info->osac_y0_h = 0.0;

	shell_info->p_blur = IXO_Mirror_Blur;
	shell_info->h_blur = IXO_Mirror_Blur;
	shell_info->osac_az_p = Az_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_el_p = El_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_az_h = Az_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_el_h = El_Blur_Sigma * JDMgaussian_random();

	shell_info->osac_x0_p += Lateral_Disp_Y_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_y0_p += Lateral_Disp_Z_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_z0_p += Defocus_Blur_Sigma * JDMgaussian_random();

	shell_info->osac_x0_h += Lateral_Disp_Y_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_y0_h += Lateral_Disp_Z_Blur_Sigma * JDMgaussian_random();
	shell_info->osac_z0_h += Defocus_Blur_Sigma * JDMgaussian_random();
	shell_info++;
     }
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   IXO_Mirrors = mirror;
   return 0;

return_error:
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   deallocate_ixo_mirror (mirror);
   return -1;
}


int _marx_ixo_mirror_init (Param_File_Type *p) /*{{{*/
{
   char *file;

   if (-1 == pf_get_parameters (p, IXOMirror_Parm_Table))
     return -1;

   if (*Mirror_Sector_String != 0)
     {
	double len;

	Mirror_Sector_Ranges = _marx_parse_range_string (Mirror_Sector_String);
	if (Mirror_Sector_Ranges == NULL)
	  return -1;
	len = _marx_compute_range_length (Mirror_Sector_Ranges);
	if (len < 0)
	  Sector_Fraction = 1.0;
	else if (len == 0.0)
	  {
	     marx_error ("Mirror sector ranges are empty/invalid");
	     return -1;
	  }
	else
	  {
	     Sector_Fraction = len/360.0;
	     if (Sector_Fraction > 1.0) Sector_Fraction = 1.0;
	  }
     }
   else
     {
	Mirror_Sector_Ranges = NULL;
	Sector_Fraction = 1.0;
     }

   if (Mirror_Is_Ideal)
     {
	/* Use_Blur_Factors = 0; */
     }

   if (NULL == (file = marx_make_data_file_name (Geometry_File)))
     return -1;

   if (-1 == read_mirror_geometry_file (file))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);

   if (Mirror_Is_Ideal == 0)
     {
	if (-1 == read_opt_constants ())
	  return -1;
     }

   if (-1 == init_mirror_shells (p))
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

int _marx_ixo_mirror_reflect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   double *photon_energies;
   unsigned int n, i, *sorted_index;
   double last_energy;
   double beta = 0.0, delta = 0.0, correction_factor = 1.0;
   double source_distance;
   Wolter1_Type *last_h;
   Wolter1_Type *shell_info;

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
   if (Mirror_Is_Ideal == 0)
     {
	for (i = 0; i < n; i++)
	  {
	     at = photon_attributes + sorted_index[i];
	     if (JDMrandom () > IXO_Vignetting_Factor)
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
   shell_info = IXO_Mirrors->shell_info;
   for (i = 0; i < n; i++)
     {
	Wolter1_Type *h;
	double energy;
	int status;

	at = photon_attributes + sorted_index[i];
	project_photon_to_mirror (at, source_distance);

	h = shell_info + at->mirror_shell;

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
	  }

	last_energy = energy;
	last_h = h;

	status = reflect_from_conic (h->conic_a_p, h->conic_b_p, h->conic_c_p,
				     &at->x, &at->p,
				     h->conic_xmin_p, h->conic_xmax_p,
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

	at->x = JDMv_sum (at->x, h->to_osac_h);

#if MARX_HAS_HRMA_PITCH_YAW
	/* Now rotate into the coordinate system of the hyperbola */

	at->p = JDM3m_vector_mul (h->fwd_matrix_h, at->p);
	at->x = JDM3m_vector_mul (h->fwd_matrix_h, at->x);
#endif

	if (0 != reflect_from_conic (h->conic_a_h, h->conic_b_h, h->conic_c_h,
				     &at->x, &at->p,
				     h->conic_xmin_h, h->conic_xmax_h,
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
     }
   return 0;
}

/*}}}*/

