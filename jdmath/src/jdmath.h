/*
 Copyright (c) 2002,2013 John E. Davis

 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the Free
 Software Foundation; either version 2 of the License, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc., 675
 Mass Ave, Cambridge, MA 02139, USA.
*/
#ifndef JDMATH_H_INCLUDED
#define JDMATH_H_INCLUDED

#define JDMATH_VERSION "1.83"
#include <stdio.h>
#include <math.h>

#ifndef HAS_BASIC_TYPEDEFS_DEFINED
# if defined(__alpha__) || defined(__ALPHA__) || defined(__alpha)
#  define INT16_BASIC_TYPE	short
#  define INT32_BASIC_TYPE	int
#  define FLOAT32_BASIC_TYPE	float
#  define FLOAT64_BASIC_TYPE	double
# else
#  include <limits.h>
/* These provide defaults.  They break on a 64 bit system. */
#  if defined(INT_MAX) && defined(LONG_MAX)
#   if (LONG_MAX != 2147483647) && (INT_MAX == 2147483647)
#    define INT32_BASIC_TYPE	int
#   endif
#  endif
#  ifndef INT32_BASIC_TYPE
#   define INT32_BASIC_TYPE	long
#  endif
#  define INT16_BASIC_TYPE	short
#  define FLOAT32_BASIC_TYPE	float
#  define FLOAT64_BASIC_TYPE	double
# endif
# define HAS_BASIC_TYPEDEFS_DEFINED
#endif

#ifdef SIZEOF_INT
# if SIZEOF_INT == 2
#  undef INT16_BASIC_TYPE
#  define INT16_BASIC_TYPE int
# else
#  if SIZEOF_INT == 4
#   undef INT32_BASIC_TYPE
#   define INT32_BASIC_TYPE int
#  endif
# endif
#endif

#ifdef SIZEOF_SHORT
# if SIZEOF_SHORT == 2
#  undef INT16_BASIC_TYPE
#  define INT16_BASIC_TYPE short
# else
#  if SIZEOF_SHORT == 4
#   undef INT32_BASIC_TYPE
#   define INT32_BASIC_TYPE short
#  endif
# endif
#endif

#ifdef SIZEOF_FLOAT
# if SIZEOF_FLOAT == 4
#  undef FLOAT32_BASIC_TYPE
#  define FLOAT32_BASIC_TYPE float
# else
#  if SIZEOF_FLOAT == 8
#   undef FLOAT64_BASIC_TYPE
#   define FLOAT64_BASIC_TYPE float
#  endif
# endif
#endif

#ifdef SIZEOF_DOUBLE
# if SIZEOF_DOUBLE == 4
#  undef FLOAT32_BASIC_TYPE
#  define FLOAT32_BASIC_TYPE double
# else
#  if SIZEOF_DOUBLE == 8
#   undef FLOAT64_BASIC_TYPE
#   define FLOAT64_BASIC_TYPE double
#  endif
# endif
#endif

#ifndef INT16_TYPEDEFED
# define INT16_TYPEDEFED
  typedef INT16_BASIC_TYPE int16;
  typedef unsigned INT16_BASIC_TYPE uint16;
#endif

#ifndef INT32_TYPEDEFED
# define INT32_TYPEDEFED
  typedef INT32_BASIC_TYPE int32;
  typedef unsigned INT32_BASIC_TYPE uint32;
#endif

#ifndef FLOAT32_TYPEDEFED
# define FLOAT32_TYPEDEFED
  typedef FLOAT32_BASIC_TYPE float32;
#endif

#ifndef FLOAT64_TYPEDEFED
# define FLOAT64_TYPEDEFED
  typedef FLOAT64_BASIC_TYPE float64;
#endif

extern int JDMcheck_types (unsigned long, unsigned long, unsigned long, unsigned long);
#define JDMATH_INIT JDMcheck_types(sizeof(int16), sizeof(int32), sizeof(float32), sizeof(float64))

/* Error codes */
#define JDMATH_INVALID_PARAMETER	1
#define JDMATH_MALLOC_ERROR		2
#define JDMATH_FILE_OPEN_ERROR		3
#define JDMATH_FILE_READ_ERROR		4
#define JDMATH_FILE_WRITE_ERROR		5
#define JDMATH_FILE_CLOSE_ERROR		6
#define JDMATH_CORRUPT_FILE_ERROR	7
#define JDMATH_DIVIDE_ZERO_ERROR	8

extern int JDMath_Error;
extern int JDMUser_Break;

extern void JDMmsg_error (char *s);
extern void JDMmsg_error2 (char *, char *);

extern unsigned int JDMbinary_search_d (double, double *, unsigned int);
extern unsigned int JDMbinary_search_f (float, float *, unsigned int);

extern double JDMinterpolate_d (double, double *, double *, unsigned int);
extern float JDMinterpolate_f (float, float *, float *, unsigned int);
extern float JDMlog_interpolate_f (float, float *, float *, unsigned int);

extern int JDMinterpolate_dvector (double *, double *, unsigned int,
				   double *, double *, unsigned int);

extern int JDMinterpolate_fvector (float *, float *, unsigned int,
				   float *, float *, unsigned int);

extern int JDMinterpolate_dfvector (double *, double *, unsigned int,
					float *, float *, unsigned int);

extern int JDMlog_interpolate_fvector (float *, float *, unsigned int,
				       float *, float *, unsigned int);

extern int JDMinterpolate_n_dvector (double *, double **, unsigned int,
				     double *, double **, unsigned int,
				     unsigned int);
extern int JDMinterpolate_n_dfvector (double *, double **, unsigned int,
				      float *, float **, unsigned int,
				      unsigned int);
extern int JDMinterpolate_n_fvector (float *, float **, unsigned int,
				     float *, float **, unsigned int,
				     unsigned int);
extern int JDMlog_interpolate_n_fvector (float *, float **, unsigned int,
					 float *, float **, unsigned int,
					 unsigned int);

typedef struct _JDM_Bilinear_Interp_Type JDM_Bilinear_Interp_Type;

extern JDM_Bilinear_Interp_Type *JDM_bilinear_open (
			      double, double,   /* lower left corner coords */
			      double, double,   /* length of sides */
			      unsigned int,   /* num values */
			      double *, /* Values at LL corner */
			      double *, /* Values at LR corner */
			      double *, /* Values at UL corner */
			      double * /* Values at UR corner */);
extern void JDM_bilinear_close (JDM_Bilinear_Interp_Type *);

extern int JDM_bilinear_interp (JDM_Bilinear_Interp_Type *,
				double, double,   /* coords of point */
				double *);   /* result(s) */

extern int *JDMinteger_vector (unsigned int);
extern void JDMfree_integer_vector (int *);
extern double *JDMdouble_vector (unsigned int);
extern void JDMfree_double_vector (double *);
extern double **JDMdouble_matrix (unsigned int, unsigned int);
extern void JDMfree_double_matrix (double **, unsigned int);
extern float *JDMfloat_vector (unsigned int);
extern void JDMfree_float_vector (float *);
extern float **JDMfloat_matrix (unsigned int, unsigned int);
extern void JDMfree_float_matrix (float **, unsigned int);

/* Matrix routines */
extern int JDMgauss_elimin (double **, double *, unsigned int, double, int *);
extern int JDMgauss_elimin_n (double **, unsigned int,
			      double **, unsigned int,
			      double, int *);
extern int JDMback_subst (double **, double *, unsigned int);
extern int JDMback_subst_n (double **, unsigned int,
			    double **, unsigned int);
extern int JDMmatrix_inverse (double **, unsigned int);

extern int JDM_lu_decomp (double **, unsigned int, unsigned int *, double epsilon, int *);
extern int JDM_lu_backsubst (double **, unsigned int, unsigned int *, double *);
extern int JDM_ludecomp_inverse (double **, unsigned int);

extern int JDM_bisection (double (*f)(double, void *), double a, double b, void *cd, double *xp);

extern unsigned int *JDMsort_doubles (double *, unsigned int);
extern unsigned int *JDMsort_floats (float *, unsigned int);

extern int JDMread_column_ddata (char *, double **, int *, unsigned int, unsigned int *);
extern int JDMread_column_fdata (char *, float **, int *, unsigned int, unsigned int *);
extern int JDMread_float_xy (char *, float **, float **, int, int, unsigned int *);
extern int JDMread_double_xy (char *, double **, double **, int, int, unsigned int *);

/* Histogram routines */
extern int JDMhistogram_d (double *, unsigned int,
			   double *, unsigned int,
			   unsigned int *, int *);
extern int JDMhistogram_f (float *, unsigned int,
			   float *, unsigned int,
			   unsigned int *, int *);
#if 0
unsigned int JDMhist_examine_ddata (double *, unsigned int, double,
				    double *, double *);
unsigned int JDMhist_examine_fdata (float *, unsigned int, float,
				    float *, float *);
#endif

/* random number routines */
typedef struct _JDMRandom_Type JDMRandom_Type;
extern JDMRandom_Type *JDMcreate_random(void);
void JDMfree_random (JDMRandom_Type *);

extern uint32 JDMgenerate_uint32_random (JDMRandom_Type *);
extern double JDMgenerate_random (JDMRandom_Type *);
extern int JDMseed_random (JDMRandom_Type *, unsigned long);

extern uint32 JDMuint32_random (void);
extern double JDMrandom (void);
extern int JDMsrandom (unsigned long);

extern double JDMgaussian_random (void);
extern double JDMexpn_random (void);

/* This generator is not the best but it is fast.  For a good generator,
 * use JDMrandom.  To seed this, call it some number of times.
 */
extern uint32 JDMfast_uint32_random (void);
extern void JDMseed_fast_random (unsigned long);
extern double JDMfast_random (void);

/* Binary read/write operations */
extern unsigned int JDMread_float64 (float64 *, unsigned int, FILE *);
extern unsigned int JDMread_float32 (float32 *, unsigned int, FILE *);
extern unsigned int JDMread_int32 (int32 *, unsigned int, FILE *);
extern unsigned int JDMread_int16 (int16 *, unsigned int, FILE *);
extern unsigned int JDMwrite_float64 (float64 *, unsigned int, FILE *);
extern unsigned int JDMwrite_float32 (float32 *, unsigned int, FILE *);
extern unsigned int JDMwrite_int32 (int32 *, unsigned int, FILE *);
extern unsigned int JDMwrite_int16 (int16 *, unsigned int, FILE *);

/* The JDMstr_ functions return the position in the string where the next
 * read/write would occur
 */
extern unsigned char *JDMstr_read_int32 (int32 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_read_int16 (int16 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_read_float64 (float64 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_read_float32 (float32 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_write_int32 (int32 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_write_int16 (int16 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_write_float64 (float64 *, unsigned int, unsigned char *);
extern unsigned char *JDMstr_write_float32 (float32 *, unsigned int, unsigned char *);

unsigned char *JDMstr_read_s_int16 (short *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_s_int32 (short *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_i_int16 (int *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_i_int32 (int *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_l_int16 (long *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_l_int32 (long *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_f_float32 (float *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_f_float64 (float *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_d_float32 (double *, unsigned int, unsigned char *);
unsigned char *JDMstr_read_d_float64 (double *, unsigned int, unsigned char *);

extern unsigned int JDMread_s_int16 (short *, unsigned int, FILE *);
extern unsigned int JDMread_s_int32 (short *, unsigned int, FILE *);
extern unsigned int JDMread_i_int16 (int *, unsigned int, FILE *);
extern unsigned int JDMread_i_int32 (int *, unsigned int, FILE *);
extern unsigned int JDMread_l_int16 (long *, unsigned int, FILE *);
extern unsigned int JDMread_l_int32 (long *, unsigned int, FILE *);
extern unsigned int JDMread_f_float32 (float *, unsigned int, FILE *);
extern unsigned int JDMread_f_float64 (float *, unsigned int, FILE *);
extern unsigned int JDMread_d_float32 (double *, unsigned int, FILE *);
extern unsigned int JDMread_d_float64 (double *, unsigned int, FILE *);

extern unsigned int JDMwrite_f_float32 (float *, unsigned int, FILE *);
extern unsigned int JDMwrite_d_float32 (double *, unsigned int, FILE *);
extern unsigned int JDMwrite_f_float64 (float *, unsigned int, FILE *);
extern unsigned int JDMwrite_d_float64 (double *, unsigned int, FILE *);
extern unsigned int JDMwrite_s_int32 (short *, unsigned int, FILE *);
extern unsigned int JDMwrite_i_int32 (int *, unsigned int, FILE *);
extern unsigned int JDMwrite_l_int32 (long *, unsigned int, FILE *);
extern unsigned int JDMwrite_i_int16 (int *, unsigned int, FILE *);
extern unsigned int JDMwrite_s_int16 (short *, unsigned int, FILE *);

/* Simplex fitting routines */
typedef long *JDM_User_Type;
typedef double (*JDMFit_Funct_Type)
(
 unsigned int,			       /* ith point to be evaluated */
 double, 			       /* x[ith] */
 double *, 			       /* parm list */
 unsigned int, 			       /* num parms */
 JDM_User_Type *		       /* user defined */
 );

/* This function is a chisqr hook. */
extern double (*JDMSimplex_Chisqr)
(
 double *,			       /* x array */
 double *,			       /* y array */
 double *,			       /* dy array */
 unsigned int,			       /* num x points */
 double *,			       /* parm array */
 unsigned int,			       /* num parms */
 JDM_User_Type *,		       /* user defined */
 JDMFit_Funct_Type		       /* function to compute fit */
 );

extern int JDMsimplex_best_fit (double * /* xp */,
				double * /* yp */,
				double * /* dyp */,
				unsigned int /* npts */,
				double * /* parm_list */,
				unsigned int /* n_parms */,
				unsigned int * /* vary_list */,
				double /* scale */,
				double * /* chisqr */,
				JDMFit_Funct_Type /* f */,
				JDM_User_Type * /* user */
				);

typedef struct
{
   double reflect;
   double contract;
   double stretch;
   unsigned int max_iterations;
   double precision;
   void (*report_fun) (int, double, double, unsigned int, unsigned int);
   int (*check_validity_fun) (double *, unsigned int);
   unsigned int report_freq;
}
JDMSimplex_Type;

extern void JDMset_simplex_parameters (JDMSimplex_Type *, int);

typedef struct
{
   double r, i;
}
JDMComplex_Type;

extern JDMComplex_Type JDMComplex (double, double);
extern int JDMc_eqs (JDMComplex_Type, JDMComplex_Type);
extern int JDMc_neqs (JDMComplex_Type, JDMComplex_Type);
extern double JDMc_real (JDMComplex_Type);
extern double JDMc_imag (JDMComplex_Type);
extern JDMComplex_Type JDMc_add (JDMComplex_Type, JDMComplex_Type); /* z1 + z2 */
extern JDMComplex_Type JDMc_mul (JDMComplex_Type, JDMComplex_Type); /* z1 * z2 */
extern JDMComplex_Type JDMc_sub (JDMComplex_Type, JDMComplex_Type); /* z1 - z2 */
extern JDMComplex_Type JDMc_div (JDMComplex_Type, JDMComplex_Type); /* z1 / z2 */
extern JDMComplex_Type JDMc_exp (JDMComplex_Type);   /* exp (z) */
extern JDMComplex_Type JDMc_sqrt (JDMComplex_Type a);   /* sqrt(z) */

extern JDMComplex_Type JDMc_iexp (double, JDMComplex_Type); /* exp(iaz) */
extern JDMComplex_Type JDMc_smul (double, JDMComplex_Type); /* az */
extern JDMComplex_Type JDMc_imul (double, JDMComplex_Type); /* iaz */
extern JDMComplex_Type JDMc_sadd (double, JDMComplex_Type); /* a + z */
extern JDMComplex_Type JDMc_a_bz (double, double, JDMComplex_Type);/* a + bz */
extern JDMComplex_Type JDMc_az1_bz2 (double, JDMComplex_Type,
				     double, JDMComplex_Type);/* az_1 + bz_2 */

extern void JDMc_inc (JDMComplex_Type *, JDMComplex_Type);   /* *z += z1 */

extern double JDMc_abs (JDMComplex_Type);     /* |z| */

/* Vectors */
typedef struct
{
   double x, y, z;
}
JDMVector_Type;

extern double JDMv_length (JDMVector_Type);
extern JDMVector_Type JDMv_vector (double, double, double);
extern JDMVector_Type JDMv_smul (double, JDMVector_Type);
extern JDMVector_Type JDMv_cross_prod (JDMVector_Type, JDMVector_Type);
extern JDMVector_Type JDMv_pcross_prod (JDMVector_Type *, JDMVector_Type *);
extern JDMVector_Type JDMv_ax1_bx2 (double, JDMVector_Type, double, JDMVector_Type);
extern JDMVector_Type JDMv_ax1_bx2_cx3 (double, JDMVector_Type, double, JDMVector_Type, double, JDMVector_Type);
extern JDMVector_Type JDMv_pax1_bx2 (double, JDMVector_Type *, double, JDMVector_Type *);
extern double JDMv_dot_prod (JDMVector_Type, JDMVector_Type);
extern double JDMv_pdot_prod (JDMVector_Type *, JDMVector_Type *);
extern void JDMv_normalize (JDMVector_Type *);
extern JDMVector_Type JDMv_unit_vector (JDMVector_Type);
extern JDMVector_Type JDMv_sum (JDMVector_Type, JDMVector_Type);
extern JDMVector_Type JDMv_diff (JDMVector_Type, JDMVector_Type);
extern JDMVector_Type JDMv_rotate_unit_vector (JDMVector_Type, JDMVector_Type, double);
extern JDMVector_Type JDMv_rotate_vector1 (JDMVector_Type, JDMVector_Type, double, double);
extern JDMVector_Type JDMv_rotate_unit_vector1 (JDMVector_Type, JDMVector_Type, double, double);
extern double JDMv_distance (JDMVector_Type, JDMVector_Type);
extern void JDMv_vector_to_spherical (JDMVector_Type, double *, double *, double *);
extern void JDMv_unit_vector_to_spherical (JDMVector_Type, double *, double *);
extern JDMVector_Type JDMv_spherical_to_vector (double, double, double);
extern double JDMv_find_rotation_axis (JDMVector_Type, JDMVector_Type, JDMVector_Type *);
extern void JDMv_spherical_to_triad (double theta, double phi,
			      JDMVector_Type *r_hat, JDMVector_Type *theta_hat,
			      JDMVector_Type *phi_hat);

typedef double JDM_3Matrix_Type [3][3];
/* These define ACTIVE rotations of the coordinate system.  */
extern void JDM3m_rot_x_matrix (JDM_3Matrix_Type, double);
extern void JDM3m_rot_y_matrix (JDM_3Matrix_Type, double);
extern void JDM3m_rot_z_matrix (JDM_3Matrix_Type, double);
extern void JDM3m_rot_matrix (JDM_3Matrix_Type, JDMVector_Type, double);

extern void JDM3m_mul (JDM_3Matrix_Type, JDM_3Matrix_Type, JDM_3Matrix_Type);
extern JDMVector_Type JDM3m_vector_mul (JDM_3Matrix_Type, JDMVector_Type);

/* returns +1 if roots are real
 * returns 0 if roots are complex.  Then roots are:
 *    x1 = rplus + i rminus;
 *    x2 = rplus - i rminus;
 */
extern int JDMquadratic_root (double /* a */, double /* b */, double /* c */,
			      double * /* rplus */,
			      double * /* rminus */);

extern void JDMint_trapezoid (double (*)(double, JDM_User_Type *),
			      JDM_User_Type *, double, double,
			      unsigned int *, double *);

extern void JDMint_midpoint (double (*)(double, JDM_User_Type *),
			     JDM_User_Type *, double, double,
			     unsigned int *, double *);

/* This routine is like the JDMint_midpoint except it makes a change of
 * variable x->1/x.
 */
extern void JDMint_midpoint_inf (double (*)(double, JDM_User_Type *),
				 JDM_User_Type *, double, double,
				 unsigned int *, double *);

extern int JDMint_romberg (double (*f)(double, JDM_User_Type *),
			   JDM_User_Type *, double, double,
			   double *);

extern int JDMint_oromberg (double (*f)(double, JDM_User_Type *),
			    JDM_User_Type *, double, double,
			    int,
			    double *);

extern int JDMpoly_interp (double *, double *, unsigned int,
			   double, double *, double *);

/* Fast Fourier Transforms */

extern int JDMfftn (int /* ndim */,
		    int * /* dims */,
		    double * /* Re */,
		    double * /* Im */,
		    int /* iSign, */,
		    double /* scaling */);
extern int JDMfftnf (int /* ndim */,
		     int * /* dims */,
		     float * /* Re */,
		     float * /* Im */,
		     int /* iSign, */,
		     double /* scaling */);

extern void JDMfft_free (void);

/* Interface to GNUPLOT */

#define JDM_MAX_GNUPLOTS 10
extern char *JDMgnuplot_command;
extern int JDMgnuplot_open (unsigned int, char *);
extern int JDMgnuplot_close (unsigned int);
extern int JDMgnuplot_cmd (unsigned int, char *, ...);
extern FILE *JDMgnuplot_get_fp (unsigned int);

extern int JDMlog_grid (double *, unsigned int, double, double);
extern int JDMlog_grid_f (float *, unsigned int, float, float);

/* Binary data interface */
typedef struct
{
   FILE *fp;
#define JDMBDATA_READ_MODE	1
#define JDMBDATA_WRITE_MODE	2
   unsigned int flags;
   int data_type;
   unsigned int nrows;
   unsigned int ncols;
   char *comment;
   unsigned int comment_len;
}
JDMBData_File_Type;

extern JDMBData_File_Type *JDMbdata_open_file (char *);
extern JDMBData_File_Type *JDMbdata_create_file (char *, int,   /* data type */
						 unsigned int,   /* cols */
						 unsigned int,   /* rows */
						 char *);   /* comment */

extern int JDMbdata_close_file (JDMBData_File_Type *);

extern unsigned int JDMcheck_finite (double *, unsigned int);
extern unsigned int JDMcheck_nan (double *, unsigned int);
extern unsigned int JDMcheck_inf (double *, unsigned int, int *);
extern int JDMfinite (double);
extern int JDMisnan (double);
extern int JDMisinf (double);

extern double JDMerf (double);
extern double JDMlog_gamma (double);
extern double JDMincomplete_gamma (double a, double x);

#ifdef PI
#undef PI
#endif
#define PI 3.14159265358979323846264338327950288

extern double JDM_Mach_Eps;

#endif /* JDMATH_H_INCLUDED */
