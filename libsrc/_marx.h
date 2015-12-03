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
#ifndef _M_MARX_H_
#define _M_MARX_H_

#include <jdfits.h>


#ifdef HAVE_FSEEKO
# define FSEEK(a,b,c) fseeko(a,b,c)
# define FTELL(a) ftello(a)
#else
# define FSEEK(a,b,c) fseek(a,b,c)
# define FTELL(a) ftell(a)
#endif

extern JDMVector_Type _Marx_HRC_Geometric_Center;

typedef struct
{
   double dx, dy, dz;
   double matrix[9];
}
_Marx_Coord_Transform_Type;
extern void _marx_transform_ray (JDMVector_Type *, JDMVector_Type *,
				_Marx_Coord_Transform_Type *);
extern void _marx_transform_ray_reverse (JDMVector_Type *, JDMVector_Type *,
					 _Marx_Coord_Transform_Type *);

#define _MARX_NUM_ACIS_S_CHIPS	6
#define _MARX_NUM_ACIS_I_CHIPS	4
#define _MARX_NUM_HRC_S_CHIPS	3
#define _MARX_NUM_HRC_I_CHIPS	1

#define MARX_NUM_MIRRORS 4
extern double _Marx_HRMA_Cap_Position;

extern JDMVector_Type _Marx_HRC_Geometric_Center;

extern _Marx_Coord_Transform_Type _Marx_Det_XForm_Matrix;

typedef struct _Marx_Acis_Chip_Type
{
   int ccd_id;

   unsigned int qe_num_energies;
   float *qe_energies;
   float *qe;

   unsigned int filter_num_energies;
   float *filter_energies;
   float *filter_qe;

   char *qe_file;
   char *filter_file;

   short (*pha_fun) (struct _Marx_Acis_Chip_Type *, float, float, double, float *);
#if !MARX_HAS_ACIS_GAIN_MAP
   double read_noise;
   double energy_gain;
   double fano_factor;

   /* These are used to convert blurred energy to pha */
   double ccd_gain;
   double ccd_offset;
#endif
   double (*contam_fun)(struct _Marx_Acis_Chip_Type *, double, double, double);
}
_Marx_Acis_Chip_Type;

typedef struct
{
   char *file;
   unsigned int num_energies;
   float *energies;
   float *eff;
}
_Marx_HRC_QE_Type;

extern int _marx_hrc_read_efficiencies (_Marx_HRC_QE_Type *);
extern int _marx_hrc_s_geom_init (Param_File_Type *);
extern int _marx_hrc_i_geom_init (Param_File_Type *);
extern int _marx_hrc_i_get_pixel_size (double *dx, double *dy);
extern int _marx_hrc_s_get_pixel_size (double *dx, double *dy);

typedef struct _Marx_HRC_Blur_Parm_Type Marx_HRC_Blur_Parm_Type;
extern Marx_HRC_Blur_Parm_Type *_marx_hrc_blur_open (Param_File_Type *p, int det);
extern void _marx_hrc_blur_close (Marx_HRC_Blur_Parm_Type *bt);
extern void _marx_hrc_blur_position (Marx_HRC_Blur_Parm_Type *bt, double *dx, double *dy);

extern short _marx_hrc_compute_pha (double);

extern int _marx_acis_apply_qe_and_pha (_Marx_Acis_Chip_Type *, Marx_Photon_Attr_Type *);
extern void _marx_free_acis_chip_type (_Marx_Acis_Chip_Type *);
extern int _marx_acis_read_chip_efficiencies (_Marx_Acis_Chip_Type *);
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
extern short _marx_acis_compute_fs_pha (_Marx_Acis_Chip_Type *, float, float, double, float *);
extern short _marx_acis_compute_bs_pha (_Marx_Acis_Chip_Type *, float, float, double, float *);
#endif

#if MARX_HAS_ACIS_FEF
/* extern int _marx_init_acis_i_rmf (Param_File_Type *p); */
/* extern int _marx_init_acis_s_rmf (Param_File_Type *p); */
extern int _marx_apply_acis_rmf (_Marx_Acis_Chip_Type *c, float x, float y,
				 double energy, float *pi, short *phap);
#endif

extern void _marx_acis_apply_streak (double, Marx_Photon_Attr_Type *, Marx_Detector_Geometry_Type *);
extern int _marx_acis_get_generic_parms (Param_File_Type *);
extern int _marx_acis_contam_init (Param_File_Type *, _Marx_Acis_Chip_Type *);

extern int _Marx_Det_Ideal_Flag;
extern int _marx_set_detector_angle (double);
extern int _marx_parse_shutter_string (char *, unsigned int *, unsigned int *);

extern int _marx_init_flat_spectrum (Param_File_Type *, Marx_Spectrum_Type *);
extern int _marx_init_file_spectrum (Param_File_Type *, Marx_Spectrum_Type *);

extern int _marx_get_simple_specrum_parms (Param_File_Type *, Marx_Source_Type *, char *);

extern int _Marx_Det_Extend_Flag;
extern Marx_Detector_Geometry_Type
 *_marx_intersect_with_detector (JDMVector_Type, JDMVector_Type,
				 Marx_Detector_Geometry_Type *,
				 JDMVector_Type *, double *, double *,
				 int extend_flag);

#if MARX_HAS_DITHER
extern int _marx_init_dither (Param_File_Type *, int, double *, double *);
extern void _marx_close_dither (void);
extern int _marx_dither_photons (Marx_Photon_Type *, unsigned int *);
extern void _marx_ray_to_sky_ra_dec (Marx_Photon_Attr_Type *, double *, double *);
extern int _Marx_Dither_Mode;
# define _MARX_DITHER_MODE_NONE 		0
# define _MARX_DITHER_MODE_INTERNAL		1
# define _MARX_DITHER_MODE_ASPSOL		2

/* source.def flags */
#define _MARX_DITHER_UNSUPPORTED		1
#define _MARX_DITHER_RECORD_ONLY		2
#define _MARX_DITHER_ZERO_AMP			4

void _marx_undither_detector (Marx_Dither_Type *);
void _marx_dither_detector (Marx_Dither_Type *);
void _marx_dither_set_ray_tstart (double);
#endif

extern int _marx_init_mirror_blur (Param_File_Type *);
extern int _marx_mirror_blur (Marx_Photon_Type *);

#define GR5_MAX_DEL 0.001    // maximum deviation of subtotal of grade probabilities from unity
extern int _marx_hrc_s_detect (Marx_Photon_Type *);
extern int _marx_hrc_i_detect (Marx_Photon_Type *);
extern int _marx_acis_s_detect (Marx_Photon_Type *);
extern int _marx_acis_i_detect (Marx_Photon_Type *);
extern int _marx_hrc_s_init (Param_File_Type *);
extern int _marx_hrc_i_init (Param_File_Type *);
extern int _marx_acis_s_init (Param_File_Type *);
extern int _marx_acis_i_init (Param_File_Type *);

extern int _marx_acis_enyz_to_grade (Marx_Photon_Attr_Type *at);  // AML Oct. 5, 2018
extern void read_acis_grades_files(int ireturn[2]);               // AML Oct. 22, 2018
extern void set_print_grfits_diag();
extern void unset_print_grfits_diag();
extern void set_extra_stderr_diag();
extern void unset_extra_stderr_diag();

#if MARX_HAS_IXO_SUPPORT
extern int _marx_ixoccd_init (Param_File_Type *);
extern int _marx_ixoccd_detect (Marx_Photon_Type *);
extern int _marx_ixoxms_detect (Marx_Photon_Type *);
#endif

extern int _marx_hrc_s_compute_pixel (int id, double, double,
				      double *, double *, double *, double *);
extern int _marx_patch_hrc_s_geom (Marx_Detector_Type *);
extern int _marx_patch_hrc_i_geom (Marx_Detector_Type *);
extern int _marx_hrc_i_compute_pixel (double, double,
				      double *, double *);
extern int _marx_patch_acis_s_geom (Marx_Detector_Type *);
extern int _marx_patch_acis_i_geom (Marx_Detector_Type *);

extern int _marx_drake_reflect (Marx_Photon_Type *);
extern int _marx_drake_flat_init (Param_File_Type *);

extern int _marx_hrma_mirror_init (Param_File_Type *);
extern int _marx_hrma_mirror_reflect (Marx_Photon_Type *);
extern int _marx_ea_mirror_init (Param_File_Type *);
extern int _marx_ea_mirror_reflect (Marx_Photon_Type *);
extern int _marx_ff_mirror_init (Param_File_Type *);
extern int _marx_ff_mirror_reflect (Marx_Photon_Type *);

extern int _marx_hetg_init (Param_File_Type *);
extern int _marx_letg_init (Param_File_Type *);
extern int _marx_letg_diffract (Marx_Photon_Type *);
extern int _marx_hetg_diffract (Marx_Photon_Type *);

#if MARX_HAS_IXO_SUPPORT
extern int _marx_ixo_mirror_init (Param_File_Type *);
extern int _marx_ixo_mirror_reflect (Marx_Photon_Type *);
extern int _marx_catgs_init (Param_File_Type *);
extern int _marx_catgs_diffract (Marx_Photon_Type *);
extern int _marx_ixoxms_init (Param_File_Type *p);
#endif

extern int _marx_init_acis_i_gain_map (Param_File_Type *);
extern int _marx_init_acis_s_gain_map (Param_File_Type *);
extern short _marx_apply_acis_gain_map (_Marx_Acis_Chip_Type *, float, float, double, float *);

extern char *_marx_caldb_get_file (char *object);
extern int _marx_caldb_patch_acis_geom (Marx_Detector_Type *d);
extern int _marx_caldb_patch_aimpoint (Marx_Detector_Type *d);
extern int _marx_caldb_patch_hrc_s_geom (Marx_Detector_Type *d);
extern int _marx_read_acis_qe (int ccdid, float **enp, float **qep, unsigned int *np);
extern JDFits_Type *_marx_open_binary_hdu (char *file, char *hduname);

typedef struct
{
   char *name;
   unsigned int num_elements;
   double *data;
   double scale;

   /* private */
   int processed;
}
_Marx_Simple_Data_Type;

extern int _marx_read_simple_data_file (char *, _Marx_Simple_Data_Type *);

extern int _marx_get_vector_parm (Param_File_Type *, char *, JDMVector_Type *);

extern double _Marx_TStart_Yrs;
extern double _Marx_TStart_MJDsecs;

extern int _marx_check_monotonicity_f (float *, unsigned int);
extern int _marx_check_monotonicity_d (double *, unsigned int);

extern int _marx_compute_detector_basis (Marx_Detector_Type *);
extern Marx_Detector_Geometry_Type *_marx_find_detector_facet (Marx_Detector_Type *, int);
extern Marx_Detector_Geometry_Type *_marx_link_detector_facet_list (Marx_Detector_Geometry_Type *, unsigned int, unsigned int);

extern Marx_Detector_Type *_marx_get_acis_s_detector (void);
extern Marx_Detector_Type *_marx_get_acis_i_detector (void);
extern Marx_Detector_Type *_marx_get_hrc_s_detector (void);
extern Marx_Detector_Type *_marx_get_hrc_i_detector (void);
#if MARX_HAS_IXO_SUPPORT
extern Marx_Detector_Type *_marx_get_ixo_ccd_detector (void);
extern Marx_Detector_Type *_marx_get_ixo_xms_detector (void);
#endif

typedef struct _Marx_QE_Type Marx_QE_Type;
extern void _marx_qe_free (Marx_QE_Type *);
extern double _marx_qe_compute (Marx_QE_Type *qeinfo, double energy);
extern void _marx_qe_inc_ref (Marx_QE_Type *qeinfo);
extern Marx_QE_Type *_marx_qe_read_file (char *file, char *ext, char *encol, char *qecol, char *filtercol);

typedef struct _Marx_Range_Type Marx_Range_Type;
extern Marx_Range_Type *_marx_parse_range_string (char *str);
extern void _marx_free_range_type (Marx_Range_Type *r);
extern int _marx_is_in_range (Marx_Range_Type *r, double x);
extern double _marx_compute_range_length (Marx_Range_Type *r);

#endif /* _M_MARX_H_ */
