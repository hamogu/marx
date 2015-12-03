/* -*- mode: c; mode: fold; -*- */
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
#ifndef _MARX_H_
#define _MARX_H_
#include <jdmath.h>
#include <pfile.h>

#define MARX_VERSION 50400
#define MARX_VERSION_STRING "5.4.0"

#ifndef SIGNED_CHAR
# define SIGNED_CHAR signed char
#endif

/*{{{ Marx_Photon_Type and functions to create/destroy it */

/* These 6 quantities are used to specify the orientation of the optical
 * axis, as well as the detector offset and angle.  They are derived from
 * the aspect solution.
 */
typedef struct
{
   float ra, dec, roll;		       /* the ra,dec values are with respect to nominal pointing */
   float dy, dz, dtheta;
}
Marx_Dither_Type;

/* This structure defines the attributes of photon. */
/*{{{ Marx_Photon_Attr_Type */
typedef struct
{
   double energy;		       /* energy of the photon */

   JDMVector_Type x;		       /* position in space */
   JDMVector_Type p;		       /* direction */

   /* Misc attributes that describe the history of the photon */

   double arrival_time;		       /* time offset from Marx_Photon_Type.start_time */
   unsigned int flags;		       /* if non-zero, this photon is
					* invalid.  In that case, the individual
					* bits carry more information.
					*/
#define PHOTON_UNDETECTED	0x01
#define PHOTON_UNREFLECTED	0x02
#define PHOTON_UNDIFFRACTED	0x04
#define PHOTON_MISSED_DETECTOR	0x08
#define PHOTON_MIRROR_VBLOCKED	0x10   /* blocked by vignetting */
#define PHOTON_DRAKE_BLOCKED	0x20   /* blocked by drake */
#define PHOTON_GRATING_VBLOCKED	0x40   /* blocked by grating vignetting */
#define BAD_PHOTON_MASK		0xFF

#define PHOTON_DRAKE_REFLECTED	0x100
#define PHOTON_ACIS_STREAKED	0x200
   /* CCD and HRC parameters */
   /* Note: MARX y,z pixels are 0 based. */
   float y_pixel;		       /* photon in-dispersion index (Chip X)*/
   float z_pixel;	       /* photon cross-dispersion index (Chip Y) */
   float u_pixel;
   float v_pixel;

   Marx_Dither_Type dither_state;

   float pi;
   short pulse_height;		       /* channel num */
   unsigned int mirror_shell;	       /* describes which mirror the
					* reflected off.  This same number
					* determines which grating the photon
					* will hit.
					*/
   SIGNED_CHAR ccd_num;		       /* ccd number photon hit */

   SIGNED_CHAR detector_region;      /* region of a (HRC) detector */
   SIGNED_CHAR order;			       /* diffraction order */
   SIGNED_CHAR support_orders[4];	       /* diffraction orders from support */

   unsigned int tag;		       /* unique tag for each photon */
   // short grade;                         /* ACIS grade 0-255 */// AML April 12, 2018
   int grade;                         // AML June 10, 2019
}
Marx_Photon_Attr_Type;

/*}}}*/
/*{{{ Marx_Photon_Type */

typedef struct
{
   double source_distance;	       /* (mm) if 0, source is at infinity */
   Marx_Photon_Attr_Type *attributes;
   unsigned int n_photons;
   unsigned int max_n_photons;	       /* number of photons attribute list can
					* accomodate
					*/
   double total_time;		       /* total time to accumulate THESE photons.
					* It is really the maximum value of the
					* time element of the attribute structure.
					*/
   double start_time;		       /* start time for this group */
   unsigned int *sorted_index;	       /* indices of sorted energies */
   double *sorted_energies;	       /* sorted energies */
   unsigned int num_sorted;

   /* This is a bitmapped value that indicates which entries in the attributes
    * structure are meaningful.  In some sense, it serves as a history of the
    * ray.
    */
   unsigned long history;
/* These are always valid but for completeness... */
#define MARX_ENERGY_OK		0x00000001
#define MARX_TIME_OK		0x00000002
#define MARX_X_VECTOR_OK	0x00000004
#define MARX_P_VECTOR_OK	0x00000008
#define MARX_TAG_OK		0x00000010
/* These vary */
#define MARX_PULSEHEIGHT_OK	0x00000020
#define MARX_PI_OK		0x00000040
#define MARX_DET_PIXEL_OK	0x00000080
#define MARX_DET_NUM_OK		0x00000100
#define MARX_DET_REGION_OK	0x00000200
#define MARX_DET_UV_PIXEL_OK	0x00000400
#define MARX_MIRROR_SHELL_OK	0x00000800
#define MARX_SKY_DITHER_OK	0x00001000
#define MARX_DET_DITHER_OK	0x00002000
#define MARX_ORDER_OK		0x00100000
#define MARX_ORDER1_OK		0x00200000
#define MARX_ORDER2_OK		0x00400000
#define MARX_ORDER3_OK		0x00800000
#define MARX_ORDER4_OK		0x01000000
#define MARX_DET_GRADE_OK       0x02000000    // AML April 12, 2018

   unsigned int tag_start;	       /* start tag for this group */
}
Marx_Photon_Type;

/*}}}*/

extern void marx_prune_photons (Marx_Photon_Type *);
extern Marx_Photon_Type *marx_alloc_photon_type (unsigned int);
extern int marx_dealloc_photon_type (Marx_Photon_Type *);

/*}}}*/

/*{{{ Marx_Dump_File_Type */

typedef struct Marx_Dump_File_Type
{
   FILE *fp;
   char type;		       /* A, I, J, E, D */
#define MARX_MAX_COLUMN_NAME_LEN 16
   char colname[MARX_MAX_COLUMN_NAME_LEN];
   int32 num_rows;
   int32 num_cols;		       /* 0 ==> vector */
   struct Marx_Dump_File_Type *next;
}
Marx_Dump_File_Type;

extern Marx_Dump_File_Type *marx_open_read_dump_file (char *);
extern void marx_close_read_dump_file (Marx_Dump_File_Type *);

extern FILE *marx_create_write_dump_file (char *, Marx_Dump_File_Type *);
extern int marx_close_write_dump_file (FILE *, unsigned long);

/*}}}*/

/*{{{ Marx_Source_Type and functions to manipulate it */

/*{{{ Important Remarks about the source */

/* This simulator makes the assumption that the source is distant enough that
 * one can assume that the total integrated flux (num photons/sec/cm^2)
 * from from the source is the same at the origin as it is at the
 * front of the telescope.  Certainly this assumption is valid for
 * astrononical sources such as distant galaxies. However, it may not
 * be valid at test facilities such as XRCF. Nevertheless, the error
 * introduced by this assumption should be small and is greatest for
 * the point source.  For a point source about 200 meters from the
 * origin, the flux at the origin is about 90% of what it is at the opening
 * of the telescope.  Fortunately, the physical quantity that is affected by
 * an error in the flux is the photon rate.
 *
 * When a source module is called to initialize the photon list, it is required
 * to set the direction part (p) of the Marx_Photon_Attr_Type structure to
 * correspond to a ray that will intercept the origin.  This means that the
 * point of origin of the array can be determined by tracing the line from
 * the origin of the coordinate system back towards the source along this
 * direction.  Assuming the validity of the assumption of the previous
 * paragraph, this means that the mirror module can construct a ray at any
 * point on the aperature of the mirror and from its inferred point of origin
 * a direction can be assigned for the ray.  For sources at infinity, the
 * direction given to the ray will simply be the one given to it by the source.
 * However, for sources at finite distances, a new direction will result whose
 * justification is contingent upon the assumption described above.
 */

/*}}}*/

/*{{{ Marx_Source_Type and Marx_Spectrum_Type */

typedef struct
{
   double emin, emax, flux;
}
Marx_Flat_Spectrum_Type;

typedef struct
{
   double *energies;
   double *cum_flux;		       /* normalized to 1 */
   unsigned int num;
}
Marx_File_Spectrum_Type;

typedef struct
{
   FILE *fp;
   unsigned long history;
   double last_time;
}
Marx_Ray_Spectrum_Type;

typedef struct _Marx_Spectrum_Type
{
   unsigned int type;
#define MARX_FLAT_SPECTRUM 1
#define MARX_FILE_SPECTRUM 2
#define MARX_RAY_SPECTRUM 3
#define MARX_SAOSAC_SPECTRUM 4
   union
     {
	Marx_File_Spectrum_Type file;
	Marx_Flat_Spectrum_Type flat;
	Marx_Ray_Spectrum_Type ray;
	void *saosac;
     }
   s;
   int (*energy_function) (struct _Marx_Spectrum_Type *, double *);
   double total_flux;		       /* photons/sec/cm^2 at origin */
   void (*close_spectrum) (struct _Marx_Spectrum_Type *);
}
Marx_Spectrum_Type;

typedef struct _Marx_Source_Type
{
   unsigned int source_id;
   int (*open_source)(struct _Marx_Source_Type *);
   int (*create_photons) (struct _Marx_Source_Type *, Marx_Photon_Type *, unsigned int, unsigned int *);
   int (*close_source)(struct _Marx_Source_Type *);

   /*
    * If distance is less than or equal to zero, the source is
    * considered to be at infinity.  The vector p is a unit vector
    * FROM the source TO the origin.
    */
   double distance;
   JDMVector_Type p;
   JDMVector_Type p_normal;	       /* normal to p in approx y direction */

   /* cum spatial probability distribution */
   double *cum_radial_dist;
   double *cum_radial_rvalues;
   unsigned int num_radial_points;

   Marx_Spectrum_Type spectrum;
}
Marx_Source_Type;

/*}}}*/

extern Marx_Source_Type *marx_create_source (Param_File_Type *);
extern int marx_open_source (Marx_Source_Type *);
extern int marx_close_source (Marx_Source_Type *);
extern int marx_create_photons (Marx_Source_Type *,
				Marx_Photon_Type *,
				unsigned int,
				unsigned int *, double *);

extern int marx_average_dither(double, double *, double *, double *);

extern char *marx_caldb_get_filename (char *object);

/*}}}*/

/*{{{ Marx_Grating_Info_Type and related functions */

typedef struct
{
   double period;
   double bar_height;
   double bar_width;
   double polyimide;
   double gold;
   double chromium;
   double nickel;
   double dp_over_p;
   double dispersion_angle;
   double theta_blur;		       /* alignment error of grating (HWHM) */
   unsigned int num_orders;
   double vig;			       /* vignetting */
}
Marx_Grating_Info_Type;

extern double marx_compute_grating_efficiency (double, int,
					       Marx_Grating_Info_Type *);

/*}}}*/

extern void marx_get_random_event (double *, double *, unsigned int,
				   double *, unsigned int);

#ifdef __GNUC__
# define MARX_ATTRIBUTE_(x) __attribute__ (x)
#else
# define MARX_ATTRIBUTE_(x)
#endif
#define MARX_ATTRIBUTE_PRINTF(a,b) MARX_ATTRIBUTE_((format(printf,a,b)))

extern void marx_error (char *,...) MARX_ATTRIBUTE_PRINTF(1,2);
extern void marx_message (char *,...) MARX_ATTRIBUTE_PRINTF(1,2);
extern char *marx_make_version_string (void);
extern void *marx_realloc (void *, unsigned int);
extern void *marx_malloc (unsigned int);
extern void *marx_calloc (unsigned int, unsigned int);
extern int marx_free (void *);
extern char *marx_dircat (char *, char *);
extern char *marx_find_file_in_path (char *, char *, char);
extern int marx_file_exists (char *);
extern char *marx_make_data_file_name (char *);

extern Param_File_Type *marx_pf_parse_cmd_line (char *file, char *mode, int argc, char **argv);

extern double marx_reflectivity (double, double, double);
extern double marx_interp_reflectivity (double, double,
					float *, float *, float *,
					unsigned int);

extern int marx_detector_init (Param_File_Type *);
extern int marx_detect (Marx_Photon_Type *, int);

extern int marx_mirror_init (Param_File_Type *);
extern int marx_mirror_reflect (Marx_Photon_Type *, int);
extern double Marx_Mirror_Geometric_Area;

extern int marx_grating_init (Param_File_Type *);
extern int marx_grating_diffract (Marx_Photon_Type *, int);

#define MARX_DETECTOR_HRC_S	1
#define MARX_DETECTOR_HRC_I	2
#define MARX_DETECTOR_ACIS_S	3
#define MARX_DETECTOR_ACIS_I	4
#define MARX_DETECTOR_PLANE	5
#define MARX_DETECTOR_IXO_CATGS_CCD	6
#define MARX_DETECTOR_IXO_XMS	7
#define MARX_MAX_DETECTORS	10

#define MARX_MIRROR_EA		1
#define MARX_MIRROR_HRMA	2
#define MARX_MIRROR_FFIELD	3
#define MARX_MIRROR_IXO		4

#define MARX_GRATING_HETG	1
#define MARX_GRATING_LETG	2
#define MARX_GRATING_CATGS	3

typedef struct
{
   /* Focal Plane Coordinate system info.
    * This information is derived from section 7.2 of Rev 4.2 of the coord
    * memo.
    */
   char *fp_system_name;
   double fp_delta_s0;		       /* pixel size, radians */
   double fp_x0;		       /* fp x pixel of on-axis ray */
   double fp_y0;		       /* fp y pixel of on-axis ray */
}
Marx_FP_Coord_Type;
extern double Marx_Focal_Length;
extern Marx_FP_Coord_Type *marx_get_fp_system_info (char *);
extern int marx_mnc_to_fpc (Marx_FP_Coord_Type *,
			    JDMVector_Type *,
			    double *, double *);

extern void marx_undither_mnc (JDMVector_Type *, Marx_Dither_Type *);

typedef struct _Marx_Dectector_Geometry_Type
{
   struct _Marx_Dectector_Geometry_Type *next;
   int id;
   double x_pixel_size, y_pixel_size;  /* mm */
   unsigned int num_x_pixels;
   unsigned int num_y_pixels;

   /* Pixel values at LL corner */
   double xpixel_offset;
   double ypixel_offset;

   JDMVector_Type x_ll, x_lr, x_ur, x_ul;   /* corners */

   /* xhat, yhat, normal form an orthonormal coord system. */
   JDMVector_Type normal;
   JDMVector_Type xhat, yhat;
   double xlen, ylen;		       /* mm */

   /* methods */
   int (*pha_fun) (struct _Marx_Dectector_Geometry_Type *, double, double, double, float *);
   int (*qe_fun) (struct _Marx_Dectector_Geometry_Type *, double en, double cx, double cy, double *qep);

#ifdef MARX_DET_FACET_PRIVATE_DATA
   MARX_DET_FACET_PRIVATE_DATA
#endif
}
Marx_Detector_Geometry_Type;

typedef struct _Marx_Detector_Type
{
   int detector_type;
   int is_initialized;

   int (*tiled_pixel_map_fun) (struct _Marx_Detector_Type *,
			       Marx_Detector_Geometry_Type *,
			       int, unsigned int, unsigned int,
			       unsigned int *, unsigned int *);

   /* linked list of facets */
   Marx_Detector_Geometry_Type *facet_list;
   int first_facet_id, last_facet_id;

   char *fp_system_name;	       /* focal plane system for detector */
   Marx_FP_Coord_Type *fp_coord_info;

   /* Detector specific offset (stf_stt_offset for chandra) */
   JDMVector_Type aimpoint_offset;

   int (*print_info)(struct _Marx_Detector_Type *, FILE *);

   /* After a call to marx_get_detector_info, this structure will be
    * expressed in the STF system at the nominal aimpoint.
    */
#ifdef MARX_DETECTOR_TYPE_PRIVATE_DATA
   MARX_DETECTOR_TYPE_PRIVATE_DATA
#endif
}
Marx_Detector_Type;

extern Marx_Detector_Type *marx_get_detector_info (char *);
extern int marx_compute_tiled_pixel (Marx_Detector_Type *,
				     int, unsigned int, unsigned int,
				     unsigned int *, unsigned int *);

extern int marx_set_data_directory (char *);

extern int marx_write_photons (Marx_Photon_Type *, unsigned long,
			       char *, int, double);
extern int marx_dump (int, char **);

extern char *marx_dircat (char *, char *);

/* #define HBAR_C 1.973271e-04 keV-Microns (1986) */
 /* WAS: 1.9732858e-4
  * The new value is the product of
  * 6.5821220e-19 KeV sec * 2.99792458e14 um
  * taken from the 1986 values
  */
#define HBAR_C 1.973269631e-4 /* keV-Microns (2006 NIST value) */

#define MARX_RAYFILE_MAGIC 0xF8161961UL
extern int marx_dump_to_rayfile (char *, int,
				 Marx_Photon_Type *,
				 double);
extern int marx_dump_rayfile (char *);

extern int marx_set_grating_table_parms (double, double, unsigned int);
extern int marx_create_grating_opt_const_tables (char *);

typedef struct _Marx_WFold_Table_Type Marx_WFold_Table_Type;

extern void marx_free_wfold_table (Marx_WFold_Table_Type *);
extern Marx_WFold_Table_Type *marx_read_wfold_file (char *);

extern double marx_wfold_table_interp (Marx_WFold_Table_Type *,
				       double, double, double);
extern int marx_wfold_dump_file (char *);

extern int marx_f_read_bdat (char *, unsigned int *,
			     unsigned int, float **, ...);

/* Supported Features */
extern char *Marx_Supported_Detectors;
extern char *Marx_Supported_Sources;
extern int Marx_Verbose;

/* Marx pixlib routines */
typedef struct _Marx_Chip_To_MNC_Type Marx_Chip_To_MNC_Type;

extern Marx_Chip_To_MNC_Type *marx_allocate_chip_to_mnc (char *, int);

extern int
marx_init_chip_to_mnc (Marx_Chip_To_MNC_Type *,
		       double, /* focal_length */
		       double, double, double,   /* xyz offset of sts from fc */
		       double);	       /* angle of stf with respect to FC */

extern int
marx_chip_to_mnc (Marx_Chip_To_MNC_Type *,
		  double, double,   /* x,y pixel */
		  JDMVector_Type *);   /* returned unit vector */

extern int
marx_mnc_to_chip (Marx_Chip_To_MNC_Type *,
		  JDMVector_Type *,
		  double *, double *);

extern void
marx_mnc_to_ra_dec (JDMVector_Type *, double *, double *);

extern void marx_free_chip_to_mnc (Marx_Chip_To_MNC_Type *);
extern int marx_vector_to_tan_plane (JDMVector_Type *,
				     JDMVector_Type *,
				     double *, double *);
extern int
marx_tan_plane_to_vector (double, double, JDMVector_Type *, JDMVector_Type *);

extern void
marx_compute_ra_dec_offsets (double ra_0, double dec_0,
			     double ra, double dec,
			     double *delta_ra, double *delta_dec);

#define MARX_GRATING_HEG	1
#define MARX_GRATING_MEG	2
#define MARX_GRATING_LEG	4

typedef struct _Marx_Grating_Xform_Type Marx_Grating_Xform_Type;
extern int marx_grating_to_chip (Marx_Grating_Xform_Type *,
				 Marx_Chip_To_MNC_Type *,
				 double, int,   /* lambda (um), order */
				 double *, double *);   /* x,y pixel */

extern int marx_init_grating_xform (Marx_Grating_Xform_Type *,
				    JDMVector_Type *);

extern Marx_Grating_Xform_Type *
marx_allocate_grating_xform (Marx_Chip_To_MNC_Type *, int);
extern void marx_free_grating_xform (Marx_Grating_Xform_Type *);

/* RDB routines : See comment in marx/libsrc/rdb.c */
typedef struct _Marx_RDB_File_Type Marx_RDB_File_Type;
extern void marx_close_rdb_file (Marx_RDB_File_Type *);
extern Marx_RDB_File_Type *marx_open_rdb_file (char *);
extern int marx_rdb_get_row (Marx_RDB_File_Type *, char *, char *);
extern int marx_rdb_get_col (Marx_RDB_File_Type *, char *);
extern char *marx_rdb_get_value (Marx_RDB_File_Type *,
				 unsigned int, unsigned int);

extern int marx_get_nominal_pointing (double *ra_nom, double *dec_nom, double *roll_nom);
extern int marx_get_pointing (double *ra_pnt, double *dec_pnt, double *roll_pnt);

extern int marx_apply_acis_rmf (int ccd_id, float x, float y,
				double energy, float *pip, short *phap);

extern int marx_map_energy_to_acis_pha (int ccd_id, int x, int y, double energy, short *phap);

extern int marx_init_acis_s_rmf (Param_File_Type *p);
extern int marx_init_acis_i_rmf (Param_File_Type *p);
extern int marx_init_ixo_ccd_rmf (Param_File_Type *p);

extern int marx_set_time (double, double);

typedef struct _Marx_Subpix_Table_Type Marx_Subpix_Table_Type;
extern Marx_Subpix_Table_Type *marx_open_acis_subpix (void);
extern void marx_close_acis_subpixel (Marx_Subpix_Table_Type *);
extern int marx_compute_acis_subpix (Marx_Subpix_Table_Type *stt,
			      int ccd, float energy, int fltgrade, float *dxp, float *dyp);

/* This function computes the azimuth and elevation of the source with
 * respect to the spacecraft pointing axis, which is generally different
 * from the optical axis.
 */
extern int marx_compute_elaz (double Ra, double Dec, double *azp, double *elp);
#endif				       /* _MARX_H_ */
