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
/* The code in this file is very-VERY ugly.  It is a result of trying to keep
 * up with the ICDs and issues such as should marx produce a level 0 or a
 * level 1 file.  If you think about it, this file implements much of the level
 * 1 pipeline.
 */

#include "config.h"
#include "marx-feat.h"

/*{{{ Include Files */

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <sys/types.h>
#include <time.h>
#include <ctype.h>

#include <marx.h>
#include <jdfits.h>

/*}}}*/

#define ACIS_FITS_FILE_SPEC		"ACIS_L1_2.0"
#define HRC_FITS_FILE_SPEC		"HRC_L1_0.0"
static char *Program_Name = "marx2fits";
static char Marx2fits_Pgm[80];
#define MARX2FITS_PATCHLVL "-0"

static Marx_Detector_Type *The_Detector;
static int Simulation_Grating_Type;    /* 0==>NONE, 1==>HETG, 2==>LETG */

#define PIX_ADJ_NONE		0
#define PIX_ADJ_RANDOMIZE	1
#define PIX_ADJ_EDSER		2
#define PIX_ADJ_EXACT		3
static int Pixel_Adjust = PIX_ADJ_EDSER;
static char Pix_Adj_Str[10] = "EDSER";
static char *Pix_Adj = Pix_Adj_Str;
static double Rand_Sky = 0.0;

static Marx_Subpix_Table_Type *Acis_Subpixel_Object;

static int Simulation_Used_ACIS;
static int Simulation_Detector_Type;   /* bitmapped */
#define DETECTOR_NONE	0x00
#define DETECTOR_ACIS_S	0x01
#define DETECTOR_ACIS_I 0x02
#define DETECTOR_ACIS	(DETECTOR_ACIS_I|DETECTOR_ACIS_S)
#define DETECTOR_HRC_S	0x04
#define DETECTOR_HRC_I	0x08
#define DETECTOR_HRC	(DETECTOR_HRC_I|DETECTOR_HRC_S)

static int Pileup_Mode;
static int Simulation_Used_Dither;
static int Simulation_Used_No_Mirror;

#define MAX_ACIS_CCDID	9
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
static double Acis_PHA_Offsets[MAX_ACIS_CCDID + 1];
static double Acis_PHA_Gains[MAX_ACIS_CCDID + 1];
#endif
static double Acis_PI_Factor;
static double ACIS_Exposure_Time;
static double ACIS_Frame_Transfer_Time;

static double Nominal_RA;
static double Nominal_Dec;
static double Nominal_Roll;

static double Target_RA;
static double Target_Dec;

static double Pointing_RA;
static double Pointing_Dec;
static double Pointing_Roll;

static double Exposure_Time;
static double DT_Corr;

/*{{{ Data_Table_Type and Data_Def_Type definitions */

typedef struct /*{{{*/
{
   float64 dtt_time;

   /* marx xpixel, ypixel */
   float32 dtt_chipx;		       /* 1-based, with X.5 is at center of pixel */
   float32 dtt_chipy;

   int32 dtt_hrc_u;
   int32 dtt_hrc_v;

   /* tiled detector coords */
   int32 dtt_tdetx;
   int32 dtt_tdety;

   /* focal plane detector coords */
   float64 dtt_detx;
   float64 dtt_dety;

   /* tangent plane coords */
   float64 dtt_xsky;
   float64 dtt_ysky;

   int16 dtt_pi;
   int32 dtt_pha;
   int16 dtt_ccdid;
   int32 dtt_expno;
   float32 dtt_energy;
   float32 dtt_benergy;
   float32 dtt_marx_energy;
   float32 dtt_xpos;
   float32 dtt_ypos;
   float32 dtt_zpos;
   float32 dtt_xcos;
   float32 dtt_ycos;
   float32 dtt_zcos;
   int16 dtt_mirror;
   int16 dtt_order;
   int16 dtt_pha_island [9];
   int16 dtt_fltgrade;
   int16 dtt_grade;
   int16 dtt_node_id;
   int16 dtt_status;
   int16 dtt_nphotons;
   int32 dtt_tag;

   Marx_Dither_Type dtt_dither;
   JDMVector_Type dtt_mnc;
   int dtt_update_dither;
}

/*}}}*/
Data_Table_Type;

static Data_Table_Type Data_Table;

typedef struct _Data_Def_Type /*{{{*/
{
   int ddt_fits_type;
   void *ddt_value_ptr;
   char *ddt_filename;
   unsigned int ddt_flags;
#define DDT_REQUIRED		0x001
#define DDT_NEED_GRATING	0x002
#define DDT_NEED_ACIS		0x004
#define DDT_NEED_PILEUP		0x008
#define DDT_NOT_FOR_PILEUP	0x010
#define DDT_NEED_HRC_S		0x020
#define DDT_NEED_HRC		0x040
#define DDT_NEED_MIRROR		0x080

#define DDT_INVALID		0x100
   char *ddt_ttype;		       /* name of column */
   char *ddt_ttype_comment;
   char *ddt_tform;		       /* data type */
   char *ddt_tunit;
   char *ddt_ctype;
   int ddt_column_number;		       /* prefered column number */
   int (*ddt_compute_value)(struct _Data_Def_Type *);
   int (*ddt_write_value)(struct _Data_Def_Type *, JDFits_Type *);
   int (*ddt_open)(struct _Data_Def_Type *);
   int (*ddt_close)(struct _Data_Def_Type *);
   Marx_Dump_File_Type *ddt_dft;

   int ddt_min_max_type;
   double ddt_min_float_value;
   double ddt_max_float_value;
   long ddt_min_int_value;
   long ddt_max_int_value;
}

/*}}}*/
Data_Def_Type;

/*}}}*/

/*{{{ static function declarations */

static int write_int32 (Data_Def_Type *, JDFits_Type *);
static int write_int32_as_int16 (Data_Def_Type *, JDFits_Type *);
static int write_int16_as_int32 (Data_Def_Type *, JDFits_Type *);
static int write_int16 (Data_Def_Type *, JDFits_Type *);
static int write_float32 (Data_Def_Type *, JDFits_Type *);
static int write_float32_as_int16 (Data_Def_Type *, JDFits_Type *);
static int write_float64 (Data_Def_Type *, JDFits_Type *);
static int write_float64_as_float32 (Data_Def_Type *, JDFits_Type *);
static int write_time (Data_Def_Type *, JDFits_Type *);
static int read_int16 (Data_Def_Type *);
static int read_int32 (Data_Def_Type *);
static int read_byte_to_int16 (Data_Def_Type *);
static int read_float32 (Data_Def_Type *);
static int read_float32_to_float64 (Data_Def_Type *);
static int read_float32_add_1 (Data_Def_Type *);
static int read_float32_to_int32 (Data_Def_Type *);
static int read_int16_to_int32 (Data_Def_Type *);
#ifdef OBSOLETE_FEATURE
static int read_pha_island (Data_Def_Type *);
#endif

static int compute_xy_sky (Data_Def_Type *);

static int compute_acis_energy (Data_Def_Type *);

static int compute_tdetxy (Data_Def_Type *);
static int compute_detxy (Data_Def_Type *);
static int compute_pi (Data_Def_Type *);
#if 0
static int compute_pileup_pha (Data_Def_Type *);
#endif
static int compute_grade (Data_Def_Type *);
static int compute_fltgrade (Data_Def_Type *);
static int compute_node_id (Data_Def_Type *);
static int compute_status (Data_Def_Type *);

static int compute_expno (Data_Def_Type *);
static int read_dither_value (Data_Def_Type *);
static int read_expno_value (Data_Def_Type *);

static int open_marx_int32_file (Data_Def_Type *);
static int open_marx_int16_file (Data_Def_Type *);
static int open_marx_f32_file (Data_Def_Type *);
static int open_marx_byte_file (Data_Def_Type *);
static int close_marx_file (Data_Def_Type *);

static int open_marx_dither_file (Data_Def_Type *ddt);

static int open_marx_chip_file (Data_Def_Type *);
static int open_detxy (Data_Def_Type *);
static int close_detxy (Data_Def_Type *);

static int marx_fit_long_filename_in_fits (char *a);
/*}}}*/

#define MAX_BTABLE_COLUMNS 64

/* This table is arranged in write order.  It is dynamically constructed
 * from Data_Def_Table using the ddt_column_number field.  The elements
 * are pointers into Data_Def_Table.
 */
static Data_Def_Type *Data_Def_Write_Table [MAX_BTABLE_COLUMNS + 1];

/* The Data_Def_Table is arranged in the order necessary to compute
 * values which depend upon other values.
 */
static Data_Def_Type Data_Def_Table [] = /*{{{*/

{
   {
      'J',
      &Data_Table.dtt_tag,		       /* pointer to value */
      "tag.dat",				       /* filename */
      0,		       /* flags */
      "TAG",			       /* colname */
      "Photon Identification Tag",  /* comment */
      "J",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      1,			       /* column_number */
      read_int32,		       /* compute_value */
      write_int32,		       /* write_value */
      open_marx_int32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0,			       /* ddt_min_float_value */
      0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'D',
      &Data_Table.dtt_time,		       /* pointer to value */
      "time.dat",				       /* filename */
      DDT_REQUIRED,		       /* flags */
      "TIME",			       /* colname */
      "time since observation start",  /* comment */
      "D",			       /* type */
      "s",			       /* units */
      NULL,			       /* WCS CTYPE */
      1,			       /* column_number */
      read_float32_to_float64,	       /* compute_value */
      write_time,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_ccdid,	       /* pointer to value */
      "detector.dat",		       /* filename */
      DDT_NEED_ACIS|DDT_REQUIRED, /* flags */
      "CCD_ID",			       /* colname */
      "CCD id number",		       /* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      2,			       /* column_number */
      read_byte_to_int16,	       /* compute_value */
      write_int16,		       /* write_value */
      open_marx_byte_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'I',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      9				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_ccdid,	       /* pointer to value */
      "detector.dat",		       /* filename */
      DDT_NEED_HRC|DDT_REQUIRED,       /* flags */
      "CHIP_ID",			       /* colname */
      "MCP id number",		       /* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      2,			       /* column_number */
      read_byte_to_int16,	       /* compute_value */
      write_int16,		       /* write_value */
      open_marx_byte_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0,			       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_marx_energy,	       /* pointer to value */
      "energy.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "MARX_ENERGY",		       /* colname */
      "Energy of ray",		       /* comment */
      "E",			       /* type */
      "keV",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_benergy,	       /* pointer to value */
      "b_energy.dat",		       /* filename */
      DDT_NEED_ACIS,		       /* flags */
      "B_ENERGY",		       /* colname */
      "Energy of ray",		       /* comment */
      "E",			       /* type */
      "keV",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_xpos,	       /* pointer to value */
      "xpos.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "XPOS",			       /* colname */
      "X position in MARX coord system",   /* comment */
      "E",			       /* type */
      "mm",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_ypos,	       /* pointer to value */
      "ypos.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "YPOS",			       /* colname */
      "Y position in MARX coord system",   /* comment */
      "E",			       /* type */
      "mm",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_zpos,	       /* pointer to value */
      "zpos.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "ZPOS",			       /* colname */
      "Z position in MARX coords",     /* comment */
      "E",			       /* type */
      "mm",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_xcos,	       /* pointer to value */
      "xcos.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "XCOS",			       /* colname */
      "X direction cosine ",	       /* comment */
      "E",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_ycos,	       /* pointer to value */
      "ycos.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "YCOS",			       /* colname */
      "Y direction cosine",	       /* comment */
      "E",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_zcos,	       /* pointer to value */
      "zcos.dat",		       /* filename */
      DDT_NOT_FOR_PILEUP,	       /* flags */
      "ZCOS",			       /* colname */
      "Z direction cosine",	       /* comment */
      "E",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32,		       /* compute_value */
      write_float32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
#if 0
   {
      'I',			       /* type */
      &Data_Table.dtt_chipx,	       /* pointer to value */
      "chipx.dat",		       /* filename */
      DDT_NEED_PILEUP|DDT_REQUIRED, /* flags */
      "CHIPX",			       /* colname */
      "CHIP X",			       /* comment */
      "I",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      4,			       /* column_number */
      read_float32_add_1,	       /* compute_value */
      write_float32_as_int16,	       /* write_value */
      open_marx_chip_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'I',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      2,			       /* ddt_min_int_value */
      1023				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_chipy,	       /* pointer to value */
      "chipy.dat",		       /* filename */
      DDT_NEED_PILEUP|DDT_REQUIRED, /* flags */
      "CHIPY",			       /* colname */
      "CHIP Y",			       /* comment */
      "I",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      5,			       /* column_number */
      read_float32_add_1,		       /* compute_value */
      write_float32_as_int16,		       /* write_value */
      open_marx_chip_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'I',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      2,			       /* ddt_min_int_value */
	1023				       /* ddt_max_int_value */
   },
#endif
   {
      'I',			       /* type */
      &Data_Table.dtt_chipx,	       /* pointer to value */
      "xpixel.dat",		       /* filename */
       DDT_REQUIRED, /* flags */
      "CHIPX",			       /* colname */
      "CHIP X",			       /* comment */
      "I",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      4,			       /* column_number */
      read_float32_add_1,	       /* compute_value */
      write_float32_as_int16,		       /* write_value */
      open_marx_chip_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_chipy,	       /* pointer to value */
      "ypixel.dat",		       /* filename */
      DDT_REQUIRED, /* flags */
      "CHIPY",			       /* colname */
      "CHIP Y",			       /* comment */
      "I",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      5,			       /* column_number */
      read_float32_add_1,      /* compute_value */
      write_float32_as_int16,	       /* write_value */
      open_marx_chip_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },

   {
      'J',			       /* type */
      &Data_Table.dtt_hrc_u,	       /* pointer to value */
      "hrc_u.dat",		       /* filename */
      DDT_NEED_HRC_S,		       /* flags */
      "U",			       /* colname */
      "HRC U",			       /* comment */
      "J",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32_to_int32,	       /* compute_value */
      write_int32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'J',			       /* type */
      &Data_Table.dtt_hrc_v,	       /* pointer to value */
      "hrc_v.dat",		       /* filename */
      DDT_NEED_HRC_S,		       /* flags */
      "V",			       /* colname */
      "HRC V",			       /* comment */
      "J",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_float32_to_int32,	       /* compute_value */
      write_int32,		       /* write_value */
      open_marx_f32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_pha,	       /* pointer to value */
      "pha.dat",		       /* filename */
#if 1
      DDT_REQUIRED,		       /* flags */
#else
      DDT_REQUIRED|DDT_NOT_FOR_PILEUP, /* flags */
#endif
      "PHA",			       /* colname */
      "Total PHA for event",	       /* comment */
      "J",			       /* type */
      "adu",			       /* units */
      NULL,			       /* WCS CTYPE */
      13,			       /* column_number */
      read_int16_to_int32,	       /* compute_value */
      write_int32,		       /* write_value */
      open_marx_int16_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      36855			       /* ddt_max_int_value */
   },
   {
      'J',			       /* type */
      &Data_Table.dtt_expno,	       /* pointer to value */
      "frame.dat",		       /* filename */
      DDT_NEED_ACIS|DDT_NEED_PILEUP,   /* flags */
      "EXPNO",			       /* colname */
      "Exposure number",	       /* comment */
      "J",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      3,			       /* column_number */
      read_expno_value,		       /* compute_value */
      write_int32,		       /* write_value */
      open_marx_int32_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      'J',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0x7FFFFFFF		       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_order,	       /* pointer to value */
      "order.dat",		       /* filename */
      DDT_NEED_GRATING|DDT_NOT_FOR_PILEUP,/* flags */
      "ORDER",			       /* colname */
      "Diffraction Order",	       /* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_byte_to_int16,	       /* compute_value */
      write_int16,		       /* write_value */
      open_marx_byte_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_mirror,	       /* pointer to value */
      "mirror.dat",		       /* filename */
      DDT_NEED_MIRROR,		       /* flags */
      "SHELL",			       /* colname */
      "Mirror Shell (0=1,1=3,2=4,3=6)",/* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      0,			       /* column_number */
      read_int16,	       /* compute_value */
      write_int16,		       /* write_value */
      open_marx_int16_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },

   /* The rest have no files associated with them and their values
    * depend upon the previous.
    */
   {
      'J',			       /* type */
      &Data_Table.dtt_expno,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS|DDT_NOT_FOR_PILEUP,/* flags */
      "EXPNO",			       /* colname */
      "Exposure number",	       /* comment */
      "J",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      3,			       /* column_number */
      compute_expno,		       /* compute_value */
      write_int32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'J',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0x7FFFFFFF		       /* ddt_max_int_value */
   },
   /* Note: compute_tdetxy assumes that the next two entries are sequential */
   {
      'J',			       /* type */
      &Data_Table.dtt_tdetx,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "TDETX",			       /* colname */
      "Detector X",		       /* comment */
      "J",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      6,			       /* column_number */
      compute_tdetxy,		       /* compute_value */
      write_int32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'J',			       /* type */
      &Data_Table.dtt_tdety,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "TDETY",			       /* colname */
      "Detector Y",		       /* comment */
      "J",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      7,			       /* column_number */
      NULL,			       /* compute_value (already computed) */
      write_int32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },

   /* These dither related ones must be sequential, and must occur after
    * expno has been computed.
    */
   {
      'E',			       /* type */
      &Data_Table.dtt_dither.ra,   /* pointer to value */
      "sky_ra.dat",			       /* filename */
      DDT_REQUIRED,		       /* flags */
      NULL,			       /* colname */
      NULL,			       /* comment */
      NULL,			       /* type */
      NULL,			       /* units */
      NULL,			       /* WCS CTYPE */
      -1,			       /* column_number */
      read_dither_value,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_dither_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_dither.dec,   /* pointer to value */
      "sky_dec.dat",			       /* filename */
      DDT_REQUIRED,		       /* flags */
      NULL,			       /* colname */
      NULL,			       /* comment */
      NULL,			       /* type */
      NULL,			       /* units */
      NULL,			       /* WCS CTYPE */
      -1,			       /* column_number */
      read_dither_value,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_dither_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_dither.roll,   /* pointer to value */
      "sky_roll.dat",			       /* filename */
      DDT_REQUIRED,		       /* flags */
      NULL,			       /* colname */
      NULL,			       /* comment */
      NULL,			       /* type */
      NULL,			       /* units */
      NULL,			       /* WCS CTYPE */
      -1,			       /* column_number */
      read_dither_value,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_dither_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_dither.dy,   /* pointer to value */
      "det_dy.dat",			       /* filename */
      DDT_REQUIRED,		       /* flags */
      NULL,			       /* colname */
      NULL,			       /* comment */
      NULL,			       /* type */
      NULL,			       /* units */
      NULL,			       /* WCS CTYPE */
      -1,			       /* column_number */
      read_dither_value,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_dither_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'D',			       /* type */
      &Data_Table.dtt_dither.dz,   /* pointer to value */
      "det_dz.dat",			       /* filename */
      DDT_REQUIRED,		       /* flags */
      NULL,			       /* colname */
      NULL,			       /* comment */
      NULL,			       /* type */
      NULL,			       /* units */
      NULL,			       /* WCS CTYPE */
      -1,			       /* column_number */
      read_dither_value,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_dither_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_dither.dtheta,   /* pointer to value */
      "det_theta.dat",			       /* filename */
      DDT_REQUIRED,		       /* flags */
      NULL,			       /* colname */
      NULL,			       /* comment */
      NULL,			       /* type */
      NULL,			       /* units */
      NULL,			       /* WCS CTYPE */
      -1,			       /* column_number */
      read_dither_value,		       /* compute_value */
      NULL,			       /* write_value */
      open_marx_dither_file,	       /* open */
      close_marx_file,		       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },

   {
      'I',			       /* type */
      &Data_Table.dtt_fltgrade,	       /* pointer to value */
      "grade.dat",		       /* filename */
      // NULL,			       /* filename */
      DDT_NEED_ACIS,		       /* flags */
      "FLTGRADE",			       /* colname */
      "Event Grade Code",	       /* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      16,			       /* column_number */
      read_int16,	               /* compute_value */
      write_int16,		       /* write_value */
      open_marx_int16_file,	       /* open */
      close_marx_file,		       /* close */
      // compute_fltgrade,		       /* compute_value */
      // write_int16,		       /* write_value */
      // NULL,			       /* open */
      // NULL,			       /* close */
      NULL,			       /* cdt */
      'I',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      255			       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_grade,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS,		       /* flags */
      "GRADE",			       /* colname */
      "ACIS grade code",	       /* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      17,			       /* column_number */
      compute_grade,		       /* compute_value */
      write_int16,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   /* compute_detxy assumes next two entries are sequential */
   {
      'E',			       /* type */
      &Data_Table.dtt_detx,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "DETX",			       /* colname */
      "Focal Plane X",		       /* comment */
      "E",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      8,			       /* column_number */
      compute_detxy,		       /* compute_value */
      write_float64_as_float32,	       /* write_value */
      open_detxy,			       /* open */
      close_detxy,			       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_dety,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "DETY",			       /* colname */
      "Focal Plane Y",		       /* comment */
      "E",			       /* type */
      "pixel",			       /* units */
      NULL,			       /* WCS CTYPE */
      9,			       /* column_number */
      NULL,			       /* compute_value */
      write_float64_as_float32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   /* compute_xy_sky assumes next two are sequential */
   {
      'E',			       /* type */
      &Data_Table.dtt_xsky,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "X",			       /* colname */
      "sky X pixel",		       /* comment */
      "E",			       /* type */
      "pixel",			       /* units */
      "RA---TAN",		       /* WCS CTYPE */
      10,			       /* column_number */
      compute_xy_sky,		       /* compute_value */
      write_float64_as_float32,		       /* write_value */
      NULL,			       /* open */
	NULL,			       /* close */
	NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'E',			       /* type */
      &Data_Table.dtt_ysky,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "Y",			       /* colname */
      "sky Y pixel",		       /* comment */
      "E",			       /* type */
      "pixel",			       /* units */
      "DEC--TAN",			       /* WCS CTYPE */
      11,			       /* column_number */
      NULL,			       /* compute_value */
      write_float64_as_float32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
#ifdef OBSOLETE_FEATURE
   {
      'I',			       /* type */
      &Data_Table.dtt_pha_island,      /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS|DDT_NOT_FOR_PILEUP,/* flags */
      "PHAS",			       /* colname */
      "Event island PHAs",	       /* comment */
      "9I",			       /* type */
      "adu",			       /* units */
      NULL,			       /* WCS CTYPE */
      12,			       /* column_number */
      compute_pha_island,	       /* compute_value */
      write_pha_island,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'I',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      -4096,			       /* ddt_min_int_value */
      4095			       /* ddt_max_int_value */
   },
#endif
   {
      'E',			       /* type */
      &Data_Table.dtt_energy,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS,		       /* flags */
      "ENERGY",			       /* colname */
      "Nominal energy of event",       /* comment */
      "E",			       /* type */
      "eV",			       /* units */
      NULL,			       /* WCS CTYPE */
      14,			       /* column_number */
      compute_acis_energy,	       /* compute_value */
      write_float32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'E',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      1e6,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
#if 0
   {
      'J',			       /* type */
      &Data_Table.dtt_pha,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS|DDT_NEED_PILEUP,   /* flags */
      "PHA",			       /* colname */
      "Total PHA for event",	       /* comment */
      "J",			       /* type */
      "adu",			       /* units */
      NULL,			       /* WCS CTYPE */
      13,			       /* column_number */
      compute_pileup_pha,	       /* compute_value */
      write_int32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      'X',			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      36855			       /* ddt_max_int_value */
   },
#endif
   {
      'J',			       /* type */
      &Data_Table.dtt_pi,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS,		       /* flags */
      "PI",			       /* colname */
      "pulse invariant energy of event",/* comment */
      "J",			       /* type */
      "Chan",			       /* units */
      NULL,			       /* WCS CTYPE */
      15,			       /* column_number */
      compute_pi,		       /* compute_value */
      write_int16_as_int32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
      &Data_Table.dtt_node_id,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_NEED_ACIS|DDT_REQUIRED,       /* flags */
      "NODE_ID",		       /* colname */
      "0-4",			       /* comment */
      "I",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
      3,			       /* column_number */
      compute_node_id,		       /* compute_value */
      write_int16,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0,			       /* ddt_max_int_value */
   },
   {
      'I',			       /* type */
	&Data_Table.dtt_status,	       /* pointer to value */
      NULL,			       /* filename */
      DDT_REQUIRED,		       /* flags */
      "STATUS",			       /* colname */
      "status flags",		       /* comment */
      "32X",			       /* type */
      "",			       /* units */
      NULL,			       /* WCS CTYPE */
	19,			       /* column_number */
      compute_status,			       /* compute_value */
      write_int16_as_int32,		       /* write_value */
      NULL,			       /* open */
      NULL,			       /* close */
      NULL,			       /* cdt */
      0,			       /* ddt_min_max_type */
      0.0,			       /* ddt_min_float_value */
      0.0,			       /* ddt_max_float_value */
      0,			       /* ddt_min_int_value */
      0				       /* ddt_max_int_value */
   },
   {
      0,
      NULL
   }
};

/*}}}*/

static JDFits_BTable_Keyword_Type BTable_Keywords [MAX_BTABLE_COLUMNS];

static char *Marx_Dir;
static int32 Num_Marx_Data_Values;

/* This function returns a pointer to static area!! */
static char *make_marx_filename (char *f) /*{{{*/
{
   static char file [1024];
   unsigned int len;

   strcpy (file, Marx_Dir);
   if (0 != (len = strlen (file)))
     {
	if (file[len - 1] != '/')
	  {
	     file[len++] = '/';
	     file[len] = 0;
	  }
     }
   strcpy (file + len, f);
   return file;
}

/*}}}*/
static int Num_Marx_File_Rows;

static int open_marx_file_internal (Data_Def_Type *ddt, int type) /*{{{*/
{
   Marx_Dump_File_Type *dft;
   char *file;
   int is_required;

   file = make_marx_filename (ddt->ddt_filename);

   is_required = (ddt->ddt_flags & DDT_REQUIRED);
   if (0 == is_required)
     {
	if (0 == marx_file_exists (file))
	  {
	     ddt->ddt_flags |= DDT_INVALID;
	     return 0;
	  }
     }

   fprintf (stdout, "Examining %s\n", file);

   if (NULL == (dft = marx_open_read_dump_file (file)))
     {
	marx_error ("*** Unable to open %s.", file);
	return -1;
     }

   if (Num_Marx_File_Rows != 0)
     {
	if (Num_Marx_File_Rows != dft->num_rows)
	  {
	     marx_close_read_dump_file (dft);
	     marx_error ("*** File %s has different number of elements than expected",
			 file);
	     return -1;
	  }
     }
   else Num_Marx_File_Rows = dft->num_rows;

   if ((int) dft->type != type)
     {
	marx_error ("*** %s is not of type '%c' as expected.", file, type);
	marx_close_read_dump_file (dft);
	return -1;
     }

   ddt->ddt_dft = dft;

   return 1;
}

/*}}}*/
static int open_marx_f32_file (Data_Def_Type *ddt)
{
   return open_marx_file_internal (ddt, 'E');
}

static int open_marx_int16_file (Data_Def_Type *ddt)
{
   return open_marx_file_internal (ddt, 'I');
}

static int open_marx_int32_file (Data_Def_Type *ddt)
{
   return open_marx_file_internal (ddt, 'J');
}

static int open_marx_byte_file (Data_Def_Type *ddt)
{
   return open_marx_file_internal (ddt, 'A');
}

static int open_marx_chip_file (Data_Def_Type *ddt)
{
   int status;

   status = open_marx_file_internal (ddt, 'E');
   if (status == -1)
     return -1;

   if (Simulation_Used_ACIS)
     {
	ddt->ddt_min_max_type = 'I';
	ddt->ddt_min_int_value = 2;
	ddt->ddt_max_int_value = 1023;
     }

   return 0;
}

static int open_marx_dither_file (Data_Def_Type *ddt)
{
   if (Simulation_Used_Dither == 0)
     return 0;

   return open_marx_f32_file (ddt);
}

static int close_marx_file (Data_Def_Type *ddt) /*{{{*/
{
   if (ddt->ddt_dft != NULL)
     {
	marx_close_read_dump_file (ddt->ddt_dft);
	ddt->ddt_dft = NULL;
     }

   return 0;
}

/*}}}*/

static int close_data_def_table (void) /*{{{*/
{
   Data_Def_Type *ddt;
   int ret = 0;

   ddt = Data_Def_Table;

   while (ddt->ddt_value_ptr != NULL)
     {
	if ((ddt->ddt_close != NULL)
	    && (-1 == (*ddt->ddt_close) (ddt)))
	  ret = -1;

	ddt++;
     }
   return ret;
}
/*}}}*/
static int open_data_def_table (void) /*{{{*/
{
   Data_Def_Type *ddt;

   ddt = Data_Def_Table;

   while (ddt->ddt_value_ptr != NULL)
     {
	unsigned int flags = ddt->ddt_flags;

	if (((flags & DDT_NEED_GRATING) && (Simulation_Grating_Type == 0))
	    || ((flags & DDT_NOT_FOR_PILEUP) && Pileup_Mode)
	    || ((flags & DDT_NEED_PILEUP) && (Pileup_Mode == 0))

	    || ((flags & DDT_NEED_ACIS)
		&& (0 == Simulation_Used_ACIS))

	    || ((flags & DDT_NEED_HRC)
		&& (0 == (Simulation_Detector_Type & DETECTOR_HRC)))

	    || ((flags & DDT_NEED_MIRROR) && Simulation_Used_No_Mirror)
	    || ((flags & DDT_NEED_HRC_S)
		&& (0 == (Simulation_Detector_Type & DETECTOR_HRC_S))))
	  {
	     ddt->ddt_flags |= DDT_INVALID;
	     ddt++;
	     continue;
	  }

	if ((ddt->ddt_open != NULL)
	    && (-1 == (*ddt->ddt_open) (ddt)))
	  return -1;

	ddt++;
     }
   return 0;
}
/*}}}*/

static int compute_table_values (void) /*{{{*/
{
   Data_Def_Type *ddt;

   ddt = Data_Def_Table;

   while (ddt->ddt_value_ptr != NULL)
     {
	if ((0 == (ddt->ddt_flags & DDT_INVALID))
	    && (NULL != ddt->ddt_compute_value)
	    && (-1 == (*ddt->ddt_compute_value) (ddt)))
	  {
	     marx_error ("*** Error computing value for %s\n",
			 ddt->ddt_ttype);
	     return -1;
	  }

	ddt++;
     }
   return 0;
}

/*}}}*/
static int write_table_values (JDFits_Type *ft) /*{{{*/
{
   Data_Def_Type **ddtp, *ddt;

   ddtp = Data_Def_Write_Table;
   while (NULL != (ddt = *ddtp))
     {
	if (-1 == (*ddt->ddt_write_value) (ddt, ft))
	  {
	     marx_error ("*** Error writing value for %s\n",
		      ddt->ddt_ttype);
	     return -1;
	  }
	ddtp++;
     }
   return 0;
}
/*}}}*/

static int data_def_table_cmp (Data_Def_Type **ap, Data_Def_Type **bp)
{
   Data_Def_Type *a, *b;

   a = *ap;
   b = *bp;

   if (a->ddt_column_number == 0)
     {
	if (b->ddt_column_number > 0)
	  return 1;

	if (a - Data_Def_Table < b - Data_Def_Table) return 1;
	if (a == b) return 0;
	return -1;
     }

   if (b->ddt_column_number == 0)
     return -1;

   if (b->ddt_column_number > a->ddt_column_number)
     return -1;
   if (b->ddt_column_number < a->ddt_column_number)
     return 1;

   if (a - Data_Def_Table < b - Data_Def_Table)
     return 1;
   if (a == b) return 0;
   return -1;
}

static int init_data_def_write_table (void) /*{{{*/
{
   Data_Def_Type *ddt, **ddtp;
   void (*qsort_fun) (Data_Def_Type **, unsigned int, unsigned int,
		      int (*)(Data_Def_Type**, Data_Def_Type **));
   unsigned int num_columns;

   ddtp = Data_Def_Write_Table;
   ddt = Data_Def_Table;

   num_columns = 0;
   while (ddt->ddt_value_ptr != NULL)
     {
	if ((0 == (ddt->ddt_flags & DDT_INVALID))
	    && (ddt->ddt_write_value != NULL))
	  {
	     num_columns++;
	     *ddtp++ = ddt;
	  }

	ddt++;
     }

   *ddtp = NULL;

   qsort_fun = (void (*)(Data_Def_Type **, unsigned int, unsigned int,
			 int (*)(Data_Def_Type**, Data_Def_Type **))) qsort;

   if (num_columns > 1)
     (*qsort_fun) (Data_Def_Write_Table, num_columns, sizeof (Data_Def_Type *),
		   data_def_table_cmp);

   return 0;
}
/*}}}*/

static void patch_min_max_values (Data_Def_Type *ddt)
{
   char *ttype;

   ttype = ddt->ddt_ttype;

   /* This would probably be more elegantly expressed by a table lookup. */

   if (0 == strcmp (ttype, "CHIP_ID")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'I';
	if (Simulation_Detector_Type & DETECTOR_HRC_S)
	  {
	     ddt->ddt_min_int_value = 1;
	     ddt->ddt_max_int_value = 3;
	  }
	else
	  {
	     ddt->ddt_min_int_value = 0;
	     ddt->ddt_max_int_value = 0;
	  }

	return;
     }

/*}}}*/

   if (0 == strcmp (ttype, "TDETX")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'J';
	if (Simulation_Used_ACIS)
	  {
	     ddt->ddt_min_int_value = 2;
	     ddt->ddt_max_int_value = 8191;
	  }
	else
	  {
	     /* HRC */
	     ddt->ddt_min_int_value = 1;
	     if (Simulation_Detector_Type & DETECTOR_HRC_S)
	       ddt->ddt_max_int_value = 49368;
	     else
	       ddt->ddt_max_int_value = 16384;
	  }

	return;
     }

/*}}}*/

   if (0 == strcmp (ttype, "TDETY")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'J';
	if (Simulation_Used_ACIS)
	  {
	     ddt->ddt_min_int_value = 2;
	     ddt->ddt_max_int_value = 8191;
	  }
	else
	  {
	     /* HRC */
	     ddt->ddt_min_int_value = 1;
	     if (Simulation_Detector_Type & DETECTOR_HRC_S)
	       ddt->ddt_max_int_value = 4096;
	     else
	       ddt->ddt_max_int_value = 16384;
	  }

	return;
     }

/*}}}*/

   if ((0 == strcmp (ttype, "DETY"))
       || (0 == strcmp (ttype, "DETX"))
       || (0 == strcmp (ttype, "X"))
       || (0 == strcmp (ttype, "Y"))) /*{{{*/
     {
	ddt->ddt_min_max_type = 'D';
	ddt->ddt_min_float_value = 0.5;

	if (Simulation_Used_ACIS)
	  ddt->ddt_max_float_value = 8192.5;
	else
	  {
	     /* HRC */
	     if (Simulation_Detector_Type & DETECTOR_HRC_S)
	       ddt->ddt_max_float_value = 65536.5;
	     else
	       ddt->ddt_max_float_value = 32768.5;
	  }

	return;
     }

/*}}}*/

   if (0 == strcmp (ttype, "CHIPX")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'I';
	if (Simulation_Used_ACIS)
	  {
	     ddt->ddt_min_int_value = 2;
	     ddt->ddt_max_int_value = 1023;
	  }
	else
	  {
	     ddt->ddt_min_int_value = 1;
	     if (Simulation_Detector_Type & DETECTOR_HRC_S)
	       ddt->ddt_max_int_value = 4096;
	     else
	       ddt->ddt_max_int_value = 16384;
	  }

	return;
     }

/*}}}*/

   if (0 == strcmp (ttype, "CHIPY")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'I';
	if (Simulation_Used_ACIS)
	  {
	     ddt->ddt_min_int_value = 2;
	     ddt->ddt_max_int_value = 1023;
	  }
	else
	  {
	     ddt->ddt_min_int_value = 1;
	     ddt->ddt_max_int_value = 16384;
	  }

	return;
     }

/*}}}*/

   if (0 == strcmp (ttype, "PHA")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'I';
	if (Simulation_Used_ACIS)
	  {
	     ddt->ddt_min_int_value = 0;
	     ddt->ddt_max_int_value = 36855;
	  }
	else
	  {
	     ddt->ddt_min_int_value = 0;
	     ddt->ddt_max_int_value = 255;
	  }

	return;
     }
/*}}}*/

   if (0 == strcmp (ttype, "PI")) /*{{{*/
     {
	ddt->ddt_min_max_type = 'I';
	if (Simulation_Used_ACIS)
	  {
	     ddt->ddt_min_int_value = 1;
	     ddt->ddt_max_int_value = 1024;
	  }
	else
	  {
	     ddt->ddt_min_int_value = 1;
	     ddt->ddt_max_int_value = 1024;
	  }

	return;
     }

/*}}}*/

   fprintf (stderr, "***WARNING: Column %s does not have TLMIN/TLMAX specified.\n",
	    ttype);

}

static int get_column_wcs_info (Data_Def_Type *ddt,
				double *crval, double *crpix, double *cdelt)
{
   char *ttype;
   ttype = ddt->ddt_ttype;

   if ((0 == strcmp (ttype, "X"))
	|| (0 == strcmp (ttype, "Y")))
     {
	if ((The_Detector == NULL)
	    || (The_Detector->fp_coord_info == NULL))
	  return -1;

	if (*ttype == 'X')
	  {
	     *crpix = The_Detector->fp_coord_info->fp_x0;
	     *cdelt = -The_Detector->fp_coord_info->fp_delta_s0 * (180.0/PI);
	     *crval = Nominal_RA;
	  }
	else
	  {
	     *crpix = The_Detector->fp_coord_info->fp_y0;
	     *cdelt = The_Detector->fp_coord_info->fp_delta_s0 * (180.0/PI);
	     *crval = Nominal_Dec;
	  }

	return 0;
     }

   return -1;
}

static int create_btable_keywords (void) /*{{{*/
{
   JDFits_BTable_Keyword_Type *bkw;
   Data_Def_Type **ddtp, *ddt;

   bkw = BTable_Keywords;
   ddtp = Data_Def_Write_Table;
   while ((ddt = *ddtp) != NULL)
     {
	memset ((char *) bkw, 0, sizeof (JDFits_BTable_Keyword_Type));

	bkw->ttype = ddt->ddt_ttype;
	bkw->tform = ddt->ddt_tform;
	bkw->tunit = ddt->ddt_tunit;

	bkw->tunit_comment = NULL;
	bkw->tform_comment = NULL;
	bkw->ttype_comment = ddt->ddt_ttype_comment;

	bkw->min_max_type = ddt->ddt_fits_type;

	if ('X' == ddt->ddt_min_max_type)
	  patch_min_max_values (ddt);

	switch (ddt->ddt_min_max_type)
	  {
	   case 'A':
	   case 0:
	   default:
	     bkw->min_max_type = 0;
	     break;

	   case 'I':
	   case 'J':
	     bkw->min_value.j_val = ddt->ddt_min_int_value;
	     bkw->max_value.j_val = ddt->ddt_max_int_value;
	     break;

	   case 'E':
	   case 'D':
	     bkw->min_value.d_val = ddt->ddt_min_float_value;
	     bkw->max_value.d_val = ddt->ddt_max_float_value;
	     break;
	  }

	if ((ddt->ddt_ctype != NULL)
	    && (-1 != get_column_wcs_info (ddt, &bkw->crval, &bkw->crpix, &bkw->cdelt)))
	  bkw->ctype = ddt->ddt_ctype;

	ddtp++;
	bkw++;
     }

   return 0;
}

/*}}}*/

typedef struct /*{{{*/
{
   unsigned int location;	       /* a bitmapped quantity */
#define FULL_COMPONENT 1
#define SHORT_COMPONENT        2

   char *keyword;
   int type;
#define H_PINT	1
#define H_PFLT	2
#define H_STR	3
#define H_PSTR	4
#define H_COM	5
#define H_ENV	6
#define H_SMARX	7		       /* read from parameter file */
#define H_DMARX	8		       /* read from parameter file */
#define H_FILE	9   		       /* read from obs.par file */
#define H_LOG	10		       /* logical */
   void *value;
   char *comment;
}

/*}}}*/
Fits_Header_Table_Type;

static Param_File_Type *Obs_Par_Parms;

static char *GratingType;
static char *DetectorType;
static char *Instrum_Name;
static char *SourceType;
static char *HDU_Class;
static char *HDU_Class1;
static char *HDU_Class2;
static char *Content_Hdr;
static char *HDU_Name_Hdr;
#if 0
static int Naxlen_Int = 2;
static int Axlen1_Int = 32768;
static int Axlen2_Int = 32768;
#endif

static char Todays_Date [64];	       /*  */
static double Sim_X;
static double Sim_Y;
static double Sim_Z;

static double Focal_Length = 10079.77;

static double Frame_Exposure_Time;/* = 0.0000001; */
static double TimeDel;	       /* time between exposures */
static double Time_Start;
static double Int_0 = 0;
static double Int_1 = 1;
/* static double Int_2 = 2; */
static double Int_1024 = 1024;
static double FP_Temp = 153.0;

static Fits_Header_Table_Type CC_NULL_Component [] =
{
     {3,"COMMENT",	H_COM,	NULL,	"\n------- Configuration Control Component -------\n\n"},
     {3,"ORIGIN",	H_STR,	"ASC",		NULL},
     {3,"CREATOR",	H_STR,	Marx2fits_Pgm, "Program creating this file"},
#if 0
     {3, "CHECKSUM",	H_STR,	"0000000000000000", NULL},
     {3, "DATASUM",	H_STR,	"0", NULL},
#endif
     {3,"HDUNAME",	H_PSTR,	&HDU_Name_Hdr,		NULL},
     {3,"HDUDOC",	H_STR,	"ASC-FITS-2.0",		NULL},
     {3,"HDUVERS",	H_STR,	"1.0.0",		NULL},

     {3,"HDUCLASS",	H_PSTR,	&HDU_Class,		NULL},
     {3,"HDUCLAS1",	H_PSTR,	&HDU_Class1,		NULL},
     {3,"HDUCLAS2",	H_PSTR,	&HDU_Class2,		NULL},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type CC_Component [] =
{
     {3,"COMMENT",	H_COM,	NULL,	"\n------- Configuration Control Component -------\n\n"},
     {3,"ORIGIN",	H_STR,	"ASC",		NULL},
     {3,"CREATOR",	H_STR,	Marx2fits_Pgm, "Program creating this file"},
     {1,"REVISION",	H_FILE,	NULL,		"Processing revision"},
#if 0
     {3, "CHECKSUM",	H_STR,	"0000000000000000", NULL},
     {3, "DATASUM",	H_STR,	"0", NULL},
#endif
     {3,"CONTENT",	H_PSTR,	&Content_Hdr,	NULL},

     {3,"HDUNAME",	H_PSTR,	&HDU_Name_Hdr,		NULL},
     {1,"HDUSPEC",	H_STR,	"Level 1 Data Products ICD",		NULL},
     {3,"HDUDOC",	H_STR,	"ASC-FITS-2.0",		NULL},
     {3,"HDUVERS",	H_STR,	"1.0.0",		NULL},

     {3,"HDUCLASS",	H_PSTR,	&HDU_Class,		NULL},
     {3,"HDUCLAS1",	H_PSTR,	&HDU_Class1,		NULL},
     {3,"HDUCLAS2",	H_PSTR,	&HDU_Class2,		NULL},

     {1,"LONGSTRN",	H_STR, "OGIP 1.0",	"Unofficial Convention"},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Timing_Component [] =
{
     {3,"COMMENT",	H_COM,	NULL,	"\n-------- Timing Component -------\n\n"},
     {3,"DATE",		H_STR,	Todays_Date,	"Date and time of file creation (UTC)"},
     {3,"DATE-OBS",	H_FILE,	NULL,		"TT, with clock correction if CLOCKAPP"},
     {3,"DATE-END",	H_FILE,	NULL,		"TT, with clock correction if CLOCKAPP"},
     {3,"TIMESYS",	H_STR,	"TT",		"AXAF time will be Terrestrial Time"},
     {3,"MJDREF",	H_FILE,	NULL,		"MJD of clock start"},

     {3,"TIMEZERO",	H_FILE,	NULL,	"Clock Correction"},
     {3,"TIMEUNIT",	H_STR,	"s",	"seconds"},
     {1,"BTIMNULL",	H_FILE,	NULL,	"Basic Time offset (s)"},
     {1,"BTIMRATE",	H_FILE,	NULL,	"Basic Time clock rate (s/VCDUcount)"},
     {1,"BTIMDRFT",	H_FILE,	NULL,	"Basic Time clock drift (s/VCDUcount^2)"},
     {1,"BTIMCORR",	H_FILE,	NULL,	"Correction applied to Basic Time rate (s)"},

     {1,"TIMEREF",	H_STR,	"LOCAL",	"Time is local for data"},
     {1,"TASSIGN",	H_STR,	"SATELLITE",	"Source of time assignment"},
     {3,"CLOCKAPP",	H_LOG,	(void *)1,	"Clock correction applied"},

     {1,"TIERRELA",	H_FILE,	NULL,	"Short term clock stability"},
     {1,"TIERABSO",	H_FILE,	NULL,	"Absolute precision of clock correction"},

     {1,"TIMVERSN",	H_STR,	"ASC-FITS-1.1", "AXAF Fits design document"},

     {3,"TSTART",	H_FILE,	NULL,		"As in the TIME column: raw space craft clock"},
     {3,"TSTOP",	H_FILE,	NULL,		"  add TIMEZERO and MJDREF for absolute TT"},

     {1,"TIMEPIXR",	H_PINT,	(void *)&Int_0,	"Time stamp refers to start of bin"},
     {1,"TIMEDEL",	H_PFLT,	(void *)&TimeDel,	"Time resolution of data in seconds"},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_Timing_Component [] =
{
     {1,"STARTMJF", H_PINT, (void *)&Int_0,	"Major frame count at start"},
     {1,"STARTMNF", H_PINT, (void *)&Int_0,	"Minor frame count at start"},
     {1,"STARTOBT", H_PINT, (void *)&Int_0,	"Onboard MET close to STARTMJF and STARTMNF"},
     {1,"STOPMJF", H_PINT, (void *)&Int_0,	"Major frame count at stop"},
     {1,"STOPMNF", H_PINT, (void *)&Int_0,	"Minor frame count at stop"},
     {0,NULL, 0, NULL, NULL}
};

static double Defocus = 0.0;
static char *Canonical_Detname;
static char *Dither_File;
static char *Geom_File;
static char *Aimpts_File;
static char *Subpix_File;

static Fits_Header_Table_Type Obs_Info_Component [] =
{
     {3,"COMMENT",	H_COM, NULL, "\n------- Observation Information -------\n\n"},
     {1,"OBSERVER",	H_ENV,	"USER",		"Observer or PI"},
     {1,"TITLE",	H_FILE,	NULL,	"Title of Observation"},
     {3,"OBS_ID",	H_STR,	"0",		"Observation ID (*)"},
     {3,"SEQ_NUM",      H_STR,  "0", "Sequence number"},
     {3,"MISSION",	H_STR,	"AXAF",	"Advanced X-ray Astrophysics Facility"},
     {3,"TELESCOP",	H_STR,	"CHANDRA",	"Telescope used"},
     {3,"INSTRUME",	H_PSTR,	&Instrum_Name,	NULL},
     {1,"DETNAM",	H_PSTR,	&Canonical_Detname,	NULL},
     {1,"GRATING",	H_PSTR,	&GratingType,	"Grating type used"},
     {1,"OBS_MODE",	H_STR,	"POINTING",	"Observation mode"},

     {1,"SIM_X",	H_PFLT,	&Sim_X,		"SIM offset, mm"},
     {1,"SIM_Y",	H_PFLT,	&Sim_Y,		"SIM offset, mm"},
     {1,"SIM_Z",	H_PFLT,	&Sim_Z,		"SIM offset, mm"},
     {1,"DEFOCUS",	H_PFLT,	&Defocus,	"Needs clarification"},
     {1,"FOC_LEN",	H_PFLT,	&Focal_Length,	NULL},
     {1,"ONTIME",	H_FILE,	NULL,		"Ontime in seconds"},
     {1,"LIVETIME",	H_FILE,	NULL,		"Livetime in seconds"},
     {1,"EXPOSURE",	H_PFLT,	&Exposure_Time,	"Includes all corrections"},
     {1,"DTCOR",	H_PFLT,	&DT_Corr,	"Dead time correction factor"},
     {1,"OBJECT",	H_FILE, NULL,	"Source Name"},
     {1,"RA_TARG",	H_PFLT,	&Target_RA,		"Target RA in degrees"},
     {1,"DEC_TARG",	H_PFLT,	&Target_Dec,		"Target DEC in degrees"},
     {1,"RA_NOM",	H_PFLT,	&Nominal_RA,		"Nominal RA in degrees"},
     {1,"DEC_NOM",	H_PFLT,	&Nominal_Dec,		"Nominal Dec in degrees"},
     {1,"ROLL_NOM",	H_PFLT,	&Nominal_Roll,		"Nominal Roll in degrees"},
     {1,"RA_PNT",	H_PFLT,	&Pointing_RA,		"Pointing RA in degrees"},
     {1,"DEC_PNT",	H_PFLT,	&Pointing_Dec,		"Pointing Dec in degrees"},
     {1,"ROLL_PNT",	H_PFLT,	&Pointing_Roll,		"Pointing Roll in degrees"},
     {1,"EQUINOX",	H_FILE,	NULL,		"Equinox"},
     {1,"RADECSYS",	H_STR,	"ICRS",		"WCS system"},
     {1,"DATACLAS",	H_STR,	"SIMULATED",	"File contains simulated data produced by MARX"},
     {1,"DY_AVG",	H_FILE,	NULL,		"[mm] Mean DY during observation"},
     {1,"DZ_AVG",	H_FILE,	NULL,		"[mm] Mean DZ during observation"},
     {1,"DTH_AVG",	H_FILE,	NULL,		"[deg] Mean DTHETA during observation"},
     {1, "ASOLFILE",    H_PSTR, &Dither_File, NULL},
   {1,"AIMPFILE",       H_PSTR, &Aimpts_File, NULL},
   {1,"GEOMFILE",       H_PSTR, &Geom_File, NULL},
   {1,"SKYFILE",       H_STR, "/dev/null", NULL},
   {1,"TDETFILE",       H_STR, "/dev/null", NULL},
   {1,"SHELLFIL",       H_STR, "mirror.dat (MARX)", NULL}, 
     {1,"RAND_PI",      H_STR, "1.0000000000000E+00", NULL},
     {1,"RUN_ID",       H_STR, "1", "Science run index"},
     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_Obs_Info_Component [] =
{
   {1,"DATAMODE",	H_FILE,	NULL,	"telemetry mode"},
   {1,"CYCLE",		H_STR,	"P",	"Events from Primary exposures"},
   {1,"CTI_CORR",	H_LOG,	(void*)0, "CTI correction not applied"},
   {1,"CTI_APP",	H_STR,	"NNNNNNNNNN", "To ANY chips"},
   {1,"CTIFILE",	H_STR,	"NONE", NULL},
   {1,"FP_TEMP",	H_PFLT,	&FP_Temp, "[K] focal plane temperature"},
   /* Added in CIAO 4.6, used to be part of pbk file */
   {1,"OCLKPAIR",       H_STR, "8", "# of pairs of overclock pixels per output"},
   {1,"ORC_MODE",       H_STR, "0", "Output register clocking mode"},
   {1,"SUM_2X2",        H_STR, "0", "On-chip summing. 0:None; 1:Sum 2x2"},
   {1,"MASKFILE",       H_STR, "/dev/null", NULL},
   {1,"PBKFILE",        H_STR, "/dev/null", NULL},
   {1,"PIX_ADJ",        H_PSTR, &Pix_Adj, "Subpixel adjustment algorithm"},
   {1,"RAND_SKY",       H_PFLT, &Rand_Sky, NULL},
   {1,"SUBPIXFL",       H_PSTR, &Subpix_File, NULL},
   {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_S_Obs_Info_Component [] =
{
   {1,"EXPOSUR4",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR5",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR6",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR7",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR8",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR9",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"ONTIME4",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME5",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME6",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME7",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME8",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME9",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME4",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME5",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME6",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME7",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME8",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME9",	H_FILE,	NULL,	"[secs]"},

   /* These keywords are necessary for dmcopy to copy the GTI extensions.
    * I have no idea what they mean.
    */
   {1,"DSTYP1",		H_STR,	"ccd_id",	NULL},
   {1,"DSVAL1",		H_STR,	"7:7",		NULL},
   {1,"2DSVAL1",	H_STR,	"4:4",		NULL},
   {1,"3DSVAL1",	H_STR,	"5:5",		NULL},
   {1,"4DSVAL1",	H_STR,	"6:6",		NULL},
   {1,"5DSVAL1",	H_STR,	"8:8",		NULL},
   {1,"6DSVAL1",	H_STR,	"9:9",		NULL},

   {1,"DSTYP3",		H_STR,	"time",		NULL},
   {1,"DSVAL3",		H_STR,	"TABLE", 	NULL},
   {1,"DSFORM3",	H_STR,	"D", 		NULL},
   {1,"DSUNIT3",	H_STR,	"s",		NULL},
   {1,"DSREF3",		H_STR,	":GTI7",	NULL},
   {1,"2DSREF3",	H_STR,	":GTI4",	NULL},
   {1,"3DSREF3",	H_STR,	":GTI5",	NULL},
   {1,"4DSREF3",	H_STR,	":GTI6",	NULL},
   {1,"5DSREF3",	H_STR,	":GTI8",	NULL},
   {1,"6DSREF3",	H_STR,	":GTI9",	NULL},
   /* added in CIAO 4.6. Used to be part of asol file. 
    * Used for dead area calibration, so technically this shoul not be used
    * because MARX does not simulate the dead area, but it is required for 
    * many CIAO tools.
    */
   {1,"FEP_CCD",        H_STR,   "475689", "CCD to FEPID mapping, fep0 is left most digit"},

   {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_I_Obs_Info_Component [] =
{
   {1,"EXPOSUR0",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR1",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR2",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"EXPOSUR3",	H_PFLT, &Exposure_Time,	"[secs]"},
   {1,"ONTIME0",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME1",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME2",	H_FILE,	NULL,	"[secs]"},
   {1,"ONTIME3",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME0",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME1",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME2",	H_FILE,	NULL,	"[secs]"},
   {1,"LIVTIME3",	H_FILE,	NULL,	"[secs]"},

   {1,"DSTYP1",		H_STR,	"ccd_id",	NULL},
   {1,"DSVAL1",		H_STR,	"3:3",		NULL},
   {1,"2DSVAL1",	H_STR,	"0:0",		NULL},
   {1,"3DSVAL1",	H_STR,	"1:1",		NULL},
   {1,"4DSVAL1",	H_STR,	"2:2",		NULL},
   {1,"DSTYP3",		H_STR,	"time",		NULL},
   {1,"DSVAL3",		H_STR,	"TABLE", 	NULL},
   {1,"DSFORM3",	H_STR,	"D", 		NULL},
   {1,"DSUNIT3",	H_STR,	"s",		NULL},
   {1,"DSREF3",		H_STR,	":GTI3",	NULL},
   {1,"2DSREF3",	H_STR,	":GTI0",	NULL},
   {1,"3DSREF3",	H_STR,	":GTI1",	NULL},
   {1,"4DSREF3",	H_STR,	":GTI2",	NULL},
   /* added in CIAO 4.6. Used to be part of asol file. 
    * Used for dead area calibration, so technically this should not be used
    * because MARX does not simulate the dead area, but it is required for 
    * many CIAO tools.
    */
   {1,"FEP_CCD",        H_STR,   "0123xx", "CCD to FEPID mapping, fep0 is left most digit"},

   {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_Coord_Sys_Component [] =
{
   {3,"ACSYS1",	H_STR,	"CHIP:AXAF-ACIS-1.0","Ref for chip cood system"},
   {3,"ACSYS2",	H_STR,	"TDET:AXAF-ACIS-2.2","Ref for tiled det cood system"},
   {3,"ACSYS3",	H_STR,	"DET:ASC-FP-1.1","Ref for focal plane coord system"},
   {3,"ACSYS4",	H_STR,	"SKY:ASC-FP-1.1","Ref for sky coord system"},
   {3,"ACSYS5", H_STR,  "GDP:ASC-GDP-1.1","Grating coordinate system"},
   {0, NULL, 0, NULL, NULL},
};

static Fits_Header_Table_Type HRC_S_Coord_Sys_Component [] =
{
   {3,"ACSYS1",	H_STR,	"RAW:AXAF-HRC-1.1","Ref for raw detector coods"},
   {3,"ACSYS2",	H_STR,	"CHIP:AXAF-HRC-1.1","Ref for chip cood system"},
   {3,"ACSYS3",	H_STR,	"TDET:AXAF-HRC-2.7S","Ref for tiled det cood system"},
   {3,"ACSYS4",	H_STR,	"DET:ASC-FP-2.3","Ref for focal plane coord system"},
   {3,"ACSYS5",	H_STR,	"SKY:ASC-FP-2.3","Ref for sky coord system"},
   {0, NULL, 0, NULL, NULL},
};

static Fits_Header_Table_Type HRC_I_Coord_Sys_Component [] =
{
   {3,"ACSYS1",	H_STR,	"RAW:AXAF-HRC-1.1","Ref for raw detector coods"},
   {3,"ACSYS2",	H_STR,	"CHIP:AXAF-HRC-1.1","Ref for chip cood system"},
   {3,"ACSYS3",	H_STR,	"TDET:AXAF-HRC-2.3I","Ref for tiled det cood system"},
   {3,"ACSYS4",	H_STR,	"DET:ASC-FP-2.1","Ref for focal plane coord system"},
   {3,"ACSYS5",	H_STR,	"SKY:ASC-FP-2.1","Ref for sky coord system"},
   {0, NULL, 0, NULL, NULL},
};

static Fits_Header_Table_Type Acis_Faint_Header_Keywords [] =
{
     {3,"COMMENT",	H_COM,	NULL, "\nAXAF FITS Event File: ACIS Level 1\n\n"},
     {3,"READMODE",	H_STR,	"TIMED",	"CCD exposure mode"},

     {3,"FIRSTROW",	H_PINT,	(void *)&Int_1,	"Index of first row of CCD readout"},
     {3,"NROWS",	H_PINT,	(void *)&Int_1024,	"Number of rows in readout"},
     {3, "EXPTIME",	H_PFLT,	(void *)&Frame_Exposure_Time,	"Commanded exposure time in secs"},

     {3,"COMMENT",	H_COM,	NULL, "\nApplied event correction/flagging reference files\n\n"},

     {3,"BIASFIL0",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL1",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL2",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL3",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL4",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL5",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL6",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL7",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL8",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL9",	H_STR,	"/dev/null",	NULL},
     {3,"BPIXFILE",	H_STR,	"/dev/null",	NULL},

     {3,"COMMENT",	H_COM,	NULL, "\nApplied event calibration/transform reference files/systems\n\n"},

     {3,"GAINFILE",	H_STR,	"/dev/null",	NULL},
     {3,"GRD_FILE",	H_STR,	"/dev/null",	NULL},
     {3,"GRD_SCHM",	H_STR,	"ACIS",	NULL},

     {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type DM_Keywords [] =
{
     {2, "MTYPE1",	H_STR,	"chip", "DM Keyword: Descriptor name"},
     {2, "MFORM1",	H_STR,	"chipx,chipy", "DM Keyword: Descriptor value"},
     {2, "MTYPE2",	H_STR,	"tdet", "DM Keyword: Descriptor name"},
     {2, "MFORM2",	H_STR,	"tdetx,tdety", "DM Keyword: Descriptor value"},
     {2, "MTYPE3",	H_STR,	"det", "DM Keyword: Descriptor name"},
     {2, "MFORM3",	H_STR,	"detx,dety", "DM Keyword: Descriptor value"},
     {2, "MTYPE4",	H_STR,	"sky", "DM Keyword: Descriptor name"},
     {2, "MFORM4",	H_STR,	"x,y", "DM Keyword: Descriptor value"},
     {2, "MFORM5",	H_STR,	"RA,DEC", "[deg"},
     {2, "MTYPE5",	H_STR,	"EQPOS", "DM Keyword: Descriptor name"},
#if 0
     {2, "MTYPE5",	H_STR,	"CPC", "DM Keyword: Descriptor name"},
     {2, "MFORM5",	H_STR,	"CPCX,CPCY", "[mm"},
#endif
     {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Hrc_Header_Keywords [] =
{
     {0,NULL, 0, NULL, NULL}
};

static int Extver = 1;
static Fits_Header_Table_Type GoodTime_Header_Keywords [] =
{
     {2,"EXTVER",	H_PINT,	&Extver,	"Convention for HDUCLASS"},
#if 0
     {2,"HDUCLASS",	H_STR,	"OGIP",	"Convention for HDUCLASS"},
     {2,"HDUCLAS1",	H_STR,	"GTI",	"File contains good time intervals"},
     {2,"HDUCLAS2",	H_STR,	"STANDARD","File contains standard good time intervals"},
#endif
#if 0
     {2,"COMMENT",	H_COM,	NULL,		"\nData model support keywords"},
     {2,"CFIELDS",	H_PINT,	(void *)&Int_1,"Number of ASC Table Columns"},
     {2,"CNAM1",	H_STR,	"TIME",	"Time"},

     {2,"CNC1",		H_PINT,	(void *)&Int_2,"Number of FITS cols for ASC col 1"},
     {2,"CETYP1",	H_STR,	"R",	"Data is an interval"},
     {2,"CITYP1",	H_STR,	"[)",	"Interval is semi-open"},
     {2,"TDISP1",	H_STR,	"F20.6","Display format for FITS col 1"},
     {2,"TDISP2",	H_STR,	"F20.6","Display format for FITS col 2"},
#endif
#if 0
     {2,"MTYPE1",	H_STR,	"TIME", "Data Model Keyword"},
     {2,"MFORM1",	H_STR,	"START,STOP", "Data Model Keyword"},
     {2,"METYP1",	H_STR,	"R", "Data Model Keyword"},
#endif
     {0,NULL, 0, NULL, NULL}
};

static char *Dither_Model;
static char *Mirror_Type;

static double DetOffset_X;
static double DetOffset_Y;
static double DetOffset_Z;

static Param_Table_Type Parm_Table [] = /*{{{*/
{
     {"GratingType",	PF_STRING_TYPE,		&GratingType},
     {"DetectorType",	PF_STRING_TYPE,		&DetectorType},
     {"SourceType",	PF_STRING_TYPE,		&SourceType},
     {"DetOffsetX",	PF_DOUBLE_TYPE,		&DetOffset_X},
     {"DetOffsetY",	PF_DOUBLE_TYPE,		&DetOffset_Y},
     {"DetOffsetZ",	PF_DOUBLE_TYPE,		&DetOffset_Z},
     {"DitherModel",	PF_STRING_TYPE,		&Dither_Model},
     {"DitherFile",     PF_STRING_TYPE,         &Dither_File},
     {"MirrorType",	PF_STRING_TYPE,		&Mirror_Type},
     {"FocalLength",	PF_DOUBLE_TYPE,		&Focal_Length},

     {"SourceRA",	PF_DOUBLE_TYPE,		&Target_RA},
     {"SourceDec",	PF_DOUBLE_TYPE,		&Target_Dec},

     {"ACIS_Exposure_Time",	PF_REAL_TYPE,	&ACIS_Exposure_Time},
     {"ACIS_Frame_Transfer_Time",PF_REAL_TYPE,	&ACIS_Frame_Transfer_Time},

     {NULL, 0, NULL}
};

/*}}}*/

static int get_date (char *str) /*{{{*/
{
   time_t tloc;
   struct tm *tms;

   time (&tloc);
   tms = gmtime (&tloc);

   sprintf (str, "%4d-%02d-%02dT%02d:%02d:%02d",
	    1900 + tms->tm_year, tms->tm_mon + 1, tms->tm_mday,
	    tms->tm_hour, tms->tm_min, tms->tm_sec);

   return 0;
}

/*}}}*/

static Param_Table_Type Pileup_Parm_Table [] =
{
     {"FrameTime",		PF_REAL_TYPE,	&ACIS_Exposure_Time},
     {"FrameTransferTime",	PF_REAL_TYPE,	&ACIS_Frame_Transfer_Time},
     {NULL, 0, NULL}
};

static int read_pileup_parms (void)
{
   char *file;
   Param_File_Type *pf;

   Simulation_Used_No_Mirror = 1;

   file = make_marx_filename ("marxpileup.par");
   pf = pf_open_parameter_file (file, "r");
   if (pf == NULL)
     {
	marx_error ("*** Unable to open %s\n", file);
	return -1;
     }

   if (-1 == pf_get_parameters (pf, Pileup_Parm_Table))
     {
	pf_close_parameter_file (pf);
	return -1;
     }

   pf_close_parameter_file (pf);

   return 0;
}

/** Shorten string so that it fits into the value field of a fits header
 */
static int marx_fit_long_filename_in_fits (char *a)
{
  char buf[67];
  if (strlen(a) > 67)
    {
      strncpy(buf, a, 63);
      strncat(buf, "...", 63);
      strcpy(a, buf);
    }
  return 0;
}

#if 0
static
int _marx_strcasecmp (char *a, char *b)
{
   while (1)
     {
	char cha, chb;

	cha = *a;
	chb = *b;
	if (toupper(cha) != toupper(chb))
	  return (int)cha - (int)chb;

	if (cha == 0)
	  return 0;

	a++;
	b++;
     }
}
#endif
static int get_marx_pfile_info (void) /*{{{*/
{
   Param_File_Type *pf;
   char *file;
   file = make_marx_filename ("marx.par");
   pf = pf_open_parameter_file (file, "r");
   if (pf == NULL)
     {
	marx_error ("*** Unable to open %s\n", file);
	return -1;
     }

   if (-1 == pf_get_parameters (pf, Parm_Table))
     {
	pf_close_parameter_file (pf);
	return -1;
     }

   if (-1 == marx_set_data_directory ("$MARX_DATA_DIR"))
     return -1;

   if (0 == strcmp (Mirror_Type, "FLATFIELD"))
     Simulation_Used_No_Mirror = 1;

   Canonical_Detname = DetectorType;

   if (0 == strcmp (DetectorType, "ACIS-S"))
     {
	Simulation_Detector_Type = DETECTOR_ACIS_S;
	Canonical_Detname = "ACIS-456789";
     }
   else if (0 == strcmp (DetectorType, "ACIS-I"))
     {
	Simulation_Detector_Type = DETECTOR_ACIS_I;
	Canonical_Detname = "ACIS-0123";
     }
   else if (0 == strcmp (DetectorType, "HRC-S"))
     Simulation_Detector_Type = DETECTOR_HRC_S;
   else if (0 == strcmp (DetectorType, "HRC-I"))
     Simulation_Detector_Type = DETECTOR_HRC_I;

   Simulation_Used_ACIS = Simulation_Detector_Type & DETECTOR_ACIS;

  if (0 == strcmp (GratingType, "HETG"))
     Simulation_Grating_Type = MARX_GRATING_HETG;
   else if (0 == strcmp (GratingType, "LETG"))
     Simulation_Grating_Type = MARX_GRATING_LETG;
   else Simulation_Grating_Type = 0;

   The_Detector = marx_get_detector_info (DetectorType);
   if (The_Detector == NULL)
     {
	pf_close_parameter_file (pf);
	marx_error ("*** DetectorType %s not supported.\n", DetectorType);
	return -1;
     }
#if 1
   Sim_X = DetOffset_X + The_Detector->aimpoint_offset.x;
   Sim_Y = DetOffset_Y + The_Detector->aimpoint_offset.y;
   Sim_Z = DetOffset_Z + The_Detector->aimpoint_offset.z;
#else
   Sim_X = -DetOffset_X + The_Detector->aimpoint_offset.x;
   Sim_Y = -DetOffset_Y + The_Detector->aimpoint_offset.y;
   Sim_Z = -DetOffset_Z + The_Detector->aimpoint_offset.z;
#endif
   Instrum_Name = "ACIS";

   if (Simulation_Detector_Type & DETECTOR_HRC)
     Instrum_Name = "HRC";

   if (0 != strcmp (Dither_Model, "NONE"))
     {
	  Simulation_Used_Dither = 1;
     }

   if (0 == strcmp (Dither_Model, "NONE"))
     {
       Dither_File = "/dev/null";
     }
   else if (0 == strcmp (Dither_Model, "INTERNAL"))
     {
       Dither_File = "MARX Internal Model";
     }
   else if (0 == strcmp(Dither_Model,"FILE"))
     {
       // Dither_File keeps the value in marx.par
       // But if some user has Dither_File names too long to fit into fits:
       marx_fit_long_filename_in_fits(Dither_File);
     }
   else
     {
       marx_error("*** Dithermodel %sa not supported.\n", Dither_Model);
     }



   if (Simulation_Used_ACIS)
     {
#if !MARX_HAS_ACIS_GAIN_MAP && !MARX_HAS_ACIS_FEF
	unsigned int i;
	for (i = 0; i <= MAX_ACIS_CCDID; i++)
	  {
	     char buf[32];

	     sprintf (buf, "ACIS_CCD%d_Gain", i);
	     if (-1 == pf_get_double (pf, buf, Acis_PHA_Gains + i))
	       {
		  pf_close_parameter_file (pf);
		  return -1;
	       }

	     sprintf (buf, "ACIS_CCD%d_Offset", i);
	     if (-1 == pf_get_double (pf, buf, Acis_PHA_Offsets + i))
	       {
		  pf_close_parameter_file (pf);
		  return -1;
	       }
	  }
#endif
	if (-1 == pf_get_double (pf, "ACIS_eV_Per_PI", &Acis_PI_Factor))
	  {
	     pf_close_parameter_file (pf);
	     return -1;
	  }
	if (Acis_PI_Factor <= 0.0)
	  Acis_PI_Factor = 14.6;

	Acis_PI_Factor = 1.0 / Acis_PI_Factor;

	if ((-1 == pf_get_double (pf, "ACIS_Exposure_Time", &ACIS_Exposure_Time))
	    || (-1 == pf_get_double (pf, "ACIS_Frame_Transfer_Time", &ACIS_Frame_Transfer_Time)))
	  {
	     pf_close_parameter_file (pf);
	     return -1;
	  }
#if MARX_HAS_ACIS_FEF
#if 0
	if (Pileup_Mode)
	  {
	     if (Simulation_Detector_Type == DETECTOR_ACIS_I)
	       {
		  if (-1 == marx_init_acis_i_rmf (pf))
		    {
		       pf_close_parameter_file (pf);
		       return -1;
		    }
	       }
	     else if (-1 == marx_init_acis_s_rmf (pf))
	       {
		  pf_close_parameter_file (pf);
		  return -1;
	       }
	  }
#endif
#endif

     }

   pf_close_parameter_file (pf);

   if (Pileup_Mode)
     {
	if (-1 == read_pileup_parms ())
	  return -1;
     }

   if (Simulation_Used_ACIS)
     {
	Frame_Exposure_Time = TimeDel = ACIS_Exposure_Time;
	if ((ACIS_Exposure_Time > 0.0)
	    && (ACIS_Frame_Transfer_Time > 0.0))
	  TimeDel += ACIS_Frame_Transfer_Time;

	if (ACIS_Exposure_Time > 0)
	  DT_Corr = ACIS_Exposure_Time / TimeDel;

        if (NULL == (Subpix_File = marx_caldb_get_filename ("ACISSUBPIX")))
          return -1;

     }


   if (NULL == (Geom_File = marx_caldb_get_filename ("GEOM")))
     return -1;
   marx_fit_long_filename_in_fits(Geom_File);

   if (NULL == (Aimpts_File = marx_caldb_get_filename ("AIMPTS")))
     return -1;
   marx_fit_long_filename_in_fits(Aimpts_File);

   return 0;
}

/*}}}*/

static int write_parfile_value (JDFits_Type *f, Param_File_Type *p,
				Fits_Header_Table_Type *h)
{
   char buf[1024];
   char *name;
   int ret, type;
   int i;
   double d;

   name = h->keyword;

   if ((p == NULL)
       || (-1 == (type = pf_get_type (p, name))))
     {
	fprintf (stderr, "**Warning: %s not found obs.par file.\n", name);
	return 0;
     }

   switch (type)
     {
      default:
	fprintf (stderr, "**Warning: obs.par parameter %s has unsupport type.\n", name);
	return 0;

      case PF_FILE_TYPE:
      case PF_STRING_TYPE:
	if (-1 == pf_get_string (p, name, buf, sizeof (buf)))
	  return -1;
	ret = jdfits_write_header_string (f, name, buf, h->comment);
	break;

      case PF_BOOLEAN_TYPE:
	if (-1 == pf_get_boolean (p, name, &i))
	  return -1;
	ret = jdfits_write_header_logical (f, name, i, h->comment);
	break;

      case PF_INTEGER_TYPE:
	if (-1 == pf_get_integer (p, name, &i))
	  return -1;
	ret = jdfits_write_header_integer (f, name, i, h->comment);
	break;

      case PF_REAL_TYPE:
      case PF_DOUBLE_TYPE:
	if (-1 == pf_get_double (p, name, &d))
	  return -1;
	ret = jdfits_write_header_double (f, name, d, h->comment);
	break;
     }

   return ret;
}

static int write_extra_headers (JDFits_Type *ft,
				Fits_Header_Table_Type *h,
				unsigned int mask) /*{{{*/
{
   get_date (Todays_Date);

   while (h->keyword != NULL)
     {
	int ret;
	char *str;
	double d;
	int i;

	if ((h->location & mask) == 0)
	  {
	     h++;
	     continue;
	  }

	switch (h->type)
	  {
	   case H_FILE:
	     ret = write_parfile_value (ft, Obs_Par_Parms, h);
	     break;

	   case H_ENV:
	     str = getenv ((char *) h->value);
	     if (str == NULL)
	       {
		  fprintf (stderr, "Warning: Environment variable %s not set.\n", (char *) h->value);
		  str = "Unknown";
	       }
	     ret = jdfits_write_header_string (ft, h->keyword, str, h->comment);
	     break;

	   case H_STR:
	     str = (char *) h->value;
	     if (str == NULL)
	       str = "Unknown";

	     ret = jdfits_write_header_string (ft, h->keyword, str, h->comment);
	     break;

	   case H_PSTR:
	     if (h->value == NULL) str = "Unknown";
	     else str = *(char **)h->value;
	     ret = jdfits_write_header_string (ft, h->keyword, str, h->comment);
	     break;

	   case H_PFLT:
	     if (h->value == NULL) d = 0.0;
	     else d = *(double *) h->value;
	     ret = jdfits_write_header_double (ft, h->keyword, d, h->comment);
	     break;

	   case H_LOG:
	     i = (h->value != NULL);
	     ret = jdfits_write_header_logical (ft, h->keyword, i, h->comment);
	     break;

	   case H_PINT:
	     if (h->value == NULL) i = 0;
	     else i = *(int *) h->value;

	     ret = jdfits_write_header_integer (ft, h->keyword, i, h->comment);
	     break;

	   case H_COM:
	     ret = jdfits_write_header_comment (ft, h->keyword, h->comment);
	     break;

	   default:
	     fprintf (stderr, "write_extra_headers: %s: type %d not supported.\n",
		      h->keyword, h->type);
	     ret = 0;
	  }

	if (ret == -1)
	  return -1;

	if (h->value == NULL)
	  {
	     if ((h->type != H_COM) && (h->type != H_FILE) && (h->type != H_LOG))
	       fprintf (stderr, "Warning: Fits keyword %s has no value\n",
			h->keyword);
	  }

	h++;
     }
   return 0;
}

/*}}}*/

static int write_marx_par_comments (JDFits_Type *ft)
{
   char *file = make_marx_filename ("marx.par");
   return jdfits_add_comments_from_file (ft, file,
					"COMMENT", "#@#", 0);
}

static int init_marx_fits_file (JDFits_Type *ft) /*{{{*/
{
   if (-1 == jdfits_init_null_primary_header (ft))
     return -1;

   HDU_Name_Hdr = "NULL";
   HDU_Class = "ASC";
   HDU_Class1 = "";
   HDU_Class2 = "";
   Content_Hdr = "";

   if (-1 == write_extra_headers (ft, CC_NULL_Component, 1))
     return -1;
   if (-1 == write_extra_headers (ft, Timing_Component, 2))
     return -1;

   if (Simulation_Used_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Timing_Component, 2))
	  return -1;
     }

   if (-1 == write_extra_headers (ft, Obs_Info_Component, 2))
     return -1;

   if (Simulation_Used_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Obs_Info_Component, 2))
	  return -1;
     }

   if (-1 == write_marx_par_comments (ft))
     return -1;

   if (-1 == jdfits_end_header (ft))
     return -1;

   return 0;
}

/*}}}*/

static int marx2fits (JDFits_Type *ft) /*{{{*/
{
   int32 i;

   if (-1 == create_btable_keywords ())
     return -1;

   if (-1 == jdfits_create_btable_extension (ft,
					     BTable_Keywords,
					     Num_Marx_Data_Values,
					     0, 1,
					     "EVENTS"))
     return -1;

   HDU_Name_Hdr = "EVENTS";
   HDU_Class = "OGIP";
   HDU_Class1 = "EVENTS";
   HDU_Class2 = "ALL";
   Content_Hdr = "EVT1";

   if (-1 == write_extra_headers (ft, CC_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, Timing_Component, 3))
     return -1;

   if (Simulation_Used_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Timing_Component, 3))
	  return -1;
     }

   if (-1 == write_extra_headers (ft, Obs_Info_Component, 3))
     return -1;

   if (Simulation_Used_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Obs_Info_Component, 3))
	  return -1;
     }

   if (Simulation_Detector_Type & DETECTOR_ACIS_I)
     {
	if (-1 == write_extra_headers (ft, Acis_I_Obs_Info_Component, 3))
	  return -1;
     }

   if (Simulation_Detector_Type & DETECTOR_ACIS_S)
     {
	if (-1 == write_extra_headers (ft, Acis_S_Obs_Info_Component, 3))
	  return -1;
     }

   if (Simulation_Used_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Faint_Header_Keywords, 3))
	  return -1;

	if (-1 == write_extra_headers (ft, Acis_Coord_Sys_Component, 3))
	  return -1;
     }
   else
     {
	if (-1 == write_extra_headers (ft, Hrc_Header_Keywords, 3))
	  return -1;

	if (Simulation_Detector_Type & DETECTOR_HRC_I)
	  {
	     if (-1 == write_extra_headers (ft, HRC_I_Coord_Sys_Component, 3))
	       return -1;
	  }
	if (Simulation_Detector_Type & DETECTOR_HRC_S)
	  {
	     if (-1 == write_extra_headers (ft, HRC_S_Coord_Sys_Component, 3))
	       return -1;
	  }
     }

   if (-1 == write_extra_headers (ft, DM_Keywords, 2))
     return -1;

   if (-1 == write_marx_par_comments (ft))
     return -1;

   if (-1 == jdfits_end_header (ft))
     return -1;

   if (0 == Simulation_Used_ACIS)
     Data_Table.dtt_update_dither = 1;

   i = Num_Marx_File_Rows;
   while (i > 0)
     {
	i--;
	if (-1 == compute_table_values ())
	  return -1;

	if (Data_Table.dtt_pha == -1)
	  continue;

	if (-1 == write_table_values (ft))
	  return -1;

     }

   return jdfits_end_data (ft);
}

/*}}}*/

static int get_simulation_info (void) /*{{{*/
{
   Marx_Dump_File_Type *dft;
   char *file;
   int ret = 0;
   char type;
   int is_pha;

   if (-1 == get_marx_pfile_info ())
     return -1;

   /* Compute number of rows.  This is necessary because the jdfits library
    * does not have the capability to update header values.
    */
   file = make_marx_filename ("pha.dat");
   if (1 == marx_file_exists (file))
     {
	is_pha = 1;
	type = 'I';
     }
   else
     {
	is_pha = 0;
	file = make_marx_filename ("time.dat");
	type = 'E';
     }

   if (NULL == (dft = marx_open_read_dump_file (file)))
     {
	marx_error ("*** Unable to open %s.\n", file);
	return -1;
     }
   if ((int) dft->type != type)
     {
	marx_error ("*** %s is not of type '%c'", file, type);
	marx_close_read_dump_file (dft);
	return -1;
     }

   Num_Marx_Data_Values = (int) dft->num_rows;

   if (is_pha)
     {
	FILE *fp = dft->fp;
	unsigned int row = (unsigned int) Num_Marx_Data_Values;

	while (row != 0)
	  {
	     int16 phas[1024];
	     unsigned int nread = 1024, i;
	     if (row < nread)
	       nread = row;

	     nread = JDMread_int16 (phas, nread, fp);
	     if (nread == 0)
	       {
		  marx_error ("*** Error reading %s", file);
		  ret = -1;
		  break;
	       }

	     for (i = 0; i < nread; i++)
	       {
		  if (phas[i] == -1) Num_Marx_Data_Values--;
	       }
	     row -= nread;
	  }
     }

   (void) marx_close_read_dump_file (dft);
   return ret;
}

/*}}}*/
static Param_Table_Type Obspar_Parm_Table [] =
{
     {"RA_Nom",		PF_DOUBLE_TYPE,		&Nominal_RA},
     {"Dec_Nom",	PF_DOUBLE_TYPE,		&Nominal_Dec},
     {"Roll_Nom",	PF_DOUBLE_TYPE,		&Nominal_Roll},
     {"RA_Pnt",		PF_DOUBLE_TYPE,		&Pointing_RA},
     {"Dec_Pnt",	PF_DOUBLE_TYPE,		&Pointing_Dec},
     {"Roll_Pnt",	PF_DOUBLE_TYPE,		&Pointing_Roll},
     {"TSTART",		PF_DOUBLE_TYPE,		&Time_Start},
     {"EXPOSURE",	PF_DOUBLE_TYPE,		&Exposure_Time},

     {NULL, 0, NULL}
};

static Param_File_Type *read_obspar_file (void)
{
   char *file;
   Param_File_Type *pf;

   file = make_marx_filename ("obs.par");
   pf = pf_open_parameter_file (file, "r");
   if (pf == NULL)
     {
	fprintf (stderr, "Unable to open obs.par file %s\n", file);
	return NULL;
     }

   if (-1 == pf_get_parameters (pf, Obspar_Parm_Table))
     {
	pf_close_parameter_file (pf);
	return NULL;
     }

   if (DT_Corr == 0)
     DT_Corr = 1.0;

   Exposure_Time *= DT_Corr;
   return pf;
}

static int add_marx_par_to_file (JDFits_Type *ft)
{
   char *file;
   FILE *fp;
   char line[512];
   unsigned int nrows, ncols;
   JDFits_BTable_Keyword_Type keywords[2];
   char tform [32];

   file = make_marx_filename ("marx.par");
   fp = fopen (file, "r");

   if (fp == NULL)
     {
	fprintf (stderr, "***Warning: %s not opened-- not added to fits file.\n",
		 file);
	return 0;
     }

   nrows = 0;
   ncols = 0;
   while (NULL != fgets (line, sizeof(line), fp))
     {
	unsigned int len;

	len = strlen (line);
	if (len > ncols)
	  ncols = len;

	nrows++;
     }

   rewind (fp);

   memset ((char *) keywords, 0, sizeof (keywords));

   sprintf (tform, "%uA", ncols);
   keywords[0].tform = tform;
   keywords[0].ttype = "MARXPAR_LINE";

   if (-1 == jdfits_create_btable_extension (ft, keywords,
					     nrows, 0, 1,
					     "MARX_PAR"))
     {
	fclose (fp);
	return -1;
     }

   if (-1 == jdfits_end_header (ft))
     return -1;

   while (nrows--)
     {
	unsigned int len;
	memset (line, 0, sizeof (line));

	if (NULL == fgets (line, sizeof (line), fp))
	  {
	     fprintf (stderr, "Error reading %s.\n", file);
	     fclose (fp);
	     return -1;
	  }

	len = strlen (line);
	if (len && (line[len - 1] == '\n'))
	  line[len - 1] = 0;

	if (-1 == jdfits_write (ft, (unsigned char *) line, ncols))
	  {
	     fclose (fp);
	     return -1;
	  }
     }

   fclose (fp);
   return jdfits_end_data (ft);
}

static int add_goodtime_extension_n (JDFits_Type *ft, int ccdid)
{
   char *file;
   FILE *fp;
   JDFits_BTable_Keyword_Type keywords[3];
   Marx_Dump_File_Type *dft;
   float64 tmin, tmax;
   float32 t;
   char gti_n[12];

   file = make_marx_filename ("time.dat");
   dft = marx_open_read_dump_file (file);

   if (dft == NULL)
     {
	fprintf (stderr, "Unable to open time.dat.\n");
	return -1;
     }
   fp = dft->fp;

   tmin = 0.0;
   tmax = 0.0;
   while (1 == JDMread_float32 (&t, 1, fp))
     tmax = t;

   tmin += Time_Start;
   tmax += Time_Start;

   marx_close_read_dump_file (dft);

   memset ((char *) keywords, 0, sizeof (keywords));

   keywords[0].tform = "1D";
   keywords[0].ttype = "START";
   keywords[0].tunit = "s";

   keywords[1].tform = "1D";
   keywords[1].ttype = "STOP";
   keywords[1].tunit = "s";

   if (-1 == jdfits_create_btable_extension (ft, keywords,
					     1, 0, 1,
					     "GTI"))
     {
	fclose (fp);
	return -1;
     }

   if (ccdid != -1)
     {
	HDU_Name_Hdr = gti_n;
	sprintf (gti_n, "GTI%d", ccdid);
	Extver = ccdid;
     }
   else
     {
	HDU_Name_Hdr = "GTI";
	Extver = 1;
     }

   HDU_Class = "OGIP";
   HDU_Class1 = "GTI";
   HDU_Class2 = "ALL";
   Content_Hdr = "GTI";

   if (-1 == write_extra_headers (ft, CC_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, Timing_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, Obs_Info_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, GoodTime_Header_Keywords, 0xFF))
     return -1;

   if (ccdid != -1)
     {
	if ((-1 == jdfits_write_header_integer (ft, "CCD_ID", ccdid, NULL))
	    || (-1 == jdfits_write_header_integer (ft, "FEP_ID", ccdid, NULL)))
	  return -1;
     }

   if (-1 == jdfits_end_header (ft))
     return -1;

   (void) jdfits_write_float64 (ft, &tmin, 1);
   (void) jdfits_write_float64 (ft, &tmax, 1);

   return jdfits_end_data (ft);
}

static int add_goodtime_extensions (JDFits_Type *f)
{
   int imin, imax, i;

   if (0 == Simulation_Used_ACIS)
     return add_goodtime_extension_n (f, -1);

   if (Simulation_Detector_Type & DETECTOR_ACIS_I)
     {
	imin = 0;
	imax = 3;
     }
   else
     {
	imin = 4;
	imax = 9;
     }

   for (i = imin; i <= imax; i++)
     {
	if (-1 == add_goodtime_extension_n (f, i))
	  return -1;
     }
   return 0;
}

static int usage (void) /*{{{*/
{
   fprintf (stderr, "%s:\n", Marx2fits_Pgm);
   fprintf (stderr, "Usage: %s [options] marxdir outfile\n", Program_Name);
   fprintf (stderr, "Options:\n");
   fprintf (stderr, "  --pileup             Process a marxpileup simulation\n");
   fprintf (stderr, "  --pixadj=EDSER       Use a subpixel algorithm (default)\n");
   fprintf (stderr, "  --pixadj=RANDOMIZE   Randomize within a detector pixel\n");
   fprintf (stderr, "  --pixadj=NONE        Do not randomize within a detector pixel\n");
   fprintf (stderr, "  --pixadj=EXACT       Use exact chip coordinates\n");

   return 1;
}

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   char *fits_file;
   JDFits_Type *ft;

   sprintf (Marx2fits_Pgm, "marx2fits v%s", MARX_VERSION MARX2FITS_PATCHLVL);

   while (1)
     {
	char *arg;

	if (argc < 3)
	  return usage ();

	arg = argv[1];
	if ((argc == 3) && (*arg != '-'))
	  {
	     Marx_Dir = argv[1];
	     fits_file = argv[2];
	     break;
	  }

	if (0 == strcmp (arg, "--pileup"))
	  {
	     Pileup_Mode = 1;
	     argv++;
	     argc--;
	     continue;
	  }
	if (0 == strncmp (arg, "--pixadj", 8))
	  {
	     if (arg[8] == '=')
	       arg += 9;
	     else if (arg[8] == 0)
	       {
		  /* Note: argc > 2 here */
		  argv++;
		  argc--;
		  arg = argv[1];
	       }
	     else return usage ();

	     argc--; argv++;
	     if ((0 == strcmp (arg, "randomize")
		  || (0 == strcmp (arg, "RANDOMIZE"))))
	       {
		  Pixel_Adjust = PIX_ADJ_RANDOMIZE;
		  Pix_Adj = "RANDOMIZE";
		  Rand_Sky = 0.5;
		  continue;
	       }
	     if ((0 == strcmp (arg, "NONE")
		  || (0 == strcmp (arg, "none"))))
	       {
		  Pixel_Adjust = PIX_ADJ_NONE;
		  Pix_Adj = "NONE";
		  Rand_Sky = 0.;
		  continue;
	       }
	     if ((0 == strcmp (arg, "EDSER")
		  || (0 == strcmp (arg, "edser"))))
	       {
		  Pixel_Adjust = PIX_ADJ_EDSER;
		  Pix_Adj = "EDSER";
		  Rand_Sky = 0.;
		  continue;
	       }
	     if ((0 == strcmp (arg, "EXACT")
		  || (0 == strcmp (arg, "exact"))))
	       {
		  Pixel_Adjust = PIX_ADJ_EXACT;
		  Pix_Adj = "EXACT";
		  Rand_Sky = 0.;
		  continue;
	       }
	     fprintf (stderr, "***** Unsupported --pixadj option: %s\n", arg);
	     return usage ();
	  }

	if (0 == strcmp (arg, "--help"))
	  return usage ();

	fprintf (stderr, "***** Unsupported option: %s\n", arg);
	return usage ();
     }

   if (-1 == get_simulation_info ())
     return 1;

   if (Pileup_Mode
       && (0 == Simulation_Used_ACIS))
     {
	fprintf (stderr, "The simulation does not appear to be an ACIS simulation\n");
	return 1;
     }

   if (Simulation_Used_ACIS)
     {
	if (Pixel_Adjust == PIX_ADJ_EDSER)
	    {
	       if (NULL == (Acis_Subpixel_Object = marx_open_acis_subpix ()))
		 {
		    fprintf (stderr, "Error opening the subpixel file\n");
		    return 1;
		 }
	    }
     }
   else if (Pixel_Adjust == PIX_ADJ_EDSER)
     Pixel_Adjust = PIX_ADJ_RANDOMIZE;

   if (NULL == (Obs_Par_Parms = read_obspar_file ()))
     {
     }

   if (-1 == open_data_def_table ())
     return 1;

   if (-1 == init_data_def_write_table ())
     {
	(void) close_data_def_table ();
	return 1;
     }

   if (NULL == (ft = jdfits_open_file (fits_file, JDFITS_WRITE_MODE)))
     {
	marx_error ("*** Unable to open output file %s\n", fits_file);
	(void) close_data_def_table ();
	return 1;
     }

   if (-1 == init_marx_fits_file (ft))
     {
	(void) jdfits_close_file (ft);
	(void) close_data_def_table ();
	return 1;
     }

   if (-1 == marx2fits (ft))
     {
	(void) jdfits_close_file (ft);
	(void) close_data_def_table ();
	return 1;
     }

   if (-1 == add_goodtime_extensions (ft))
     {
	(void) jdfits_close_file (ft);
	(void) close_data_def_table ();
	return 1;
     }

   if (-1 == add_marx_par_to_file (ft))
     return -1;

   if (-1 == jdfits_close_file (ft))
     {
	(void) close_data_def_table ();
	return 1;
     }

   if (-1 == close_data_def_table ())
     return 1;

   return 0;
}

/*}}}*/

/*{{{ ddt_write_value functions */

static int write_int16 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   return jdfits_write_int16 (ft, (int16 *) ddt->ddt_value_ptr, 1);
}

/*}}}*/

static int write_float32_as_int16 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   int16 x = *(float32 *)ddt->ddt_value_ptr;
   return jdfits_write_int16 (ft, &x, 1);
}
/*}}}*/

static int write_int32_as_int16 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   int16 x = *(int32 *)ddt->ddt_value_ptr;
   return jdfits_write_int16 (ft, &x, 1);
}
/*}}}*/

static int write_int16_as_int32 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   int32 x = *(int16 *)ddt->ddt_value_ptr;
   return jdfits_write_int32 (ft, &x, 1);
}
/*}}}*/


static int write_float64_as_float32 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   float32 x = *(float64 *)ddt->ddt_value_ptr;
   return jdfits_write_float32 (ft, &x, 1);
}
/*}}}*/

static int write_int32 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   return jdfits_write_int32 (ft, (int32 *) ddt->ddt_value_ptr, 1);
}

/*}}}*/

static int write_float32 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   return jdfits_write_float32 (ft, (float32 *) ddt->ddt_value_ptr, 1);
}

/*}}}*/

static int write_float64 (Data_Def_Type *ddt, JDFits_Type *ft) /*{{{*/
{
   return jdfits_write_float64 (ft, (float64 *) ddt->ddt_value_ptr, 1);
}
/*}}}*/

#ifdef OBSOLETE_FEATURE
static int write_pha_island (Data_Def_Type *ddt, JDFits_Type *ft)
{
   (void) ddt;
   return jdfits_write_int16 (ft, Data_Table.dtt_pha_island, 9);
}
#endif
static int write_time (Data_Def_Type *ddt, JDFits_Type *ft)
{
   float64 t;

   t = *(float64 *)ddt->ddt_value_ptr;
   if (TimeDel > 0.0)
     t = TimeDel * Data_Table.dtt_expno;

   t += Time_Start;

   return jdfits_write_float64 (ft,  &t, 1);
}

/*}}}*/

/*{{{ ddt_compute_value functions */

/*{{{ Read data from marx functions */
static int read_int16 (Data_Def_Type *ddt) /*{{{*/
{
   if (1 != JDMread_int16 ((int16 *) ddt->ddt_value_ptr, 1, ddt->ddt_dft->fp))
     return -1;

   return 0;
}
/*}}}*/

static int read_float32_add_1 (Data_Def_Type *ddt) /*{{{*/
{
   float32 x;

   if (1 != JDMread_float32 (&x, 1, ddt->ddt_dft->fp))
     return -1;

   *(float32 *) ddt->ddt_value_ptr = x + 1;

   return 0;
}
/*}}}*/

static int read_int32 (Data_Def_Type *ddt) /*{{{*/
{
   if (1 != JDMread_int32 ((int32 *) ddt->ddt_value_ptr, 1, ddt->ddt_dft->fp))
     return -1;

   return 0;
}

/*}}}*/

static int read_float32_to_int32 (Data_Def_Type *ddt) /*{{{*/
{
   float32 f32;

   if (1 != JDMread_float32 (&f32, 1, ddt->ddt_dft->fp))
     return -1;

   *(int32 *) ddt->ddt_value_ptr = (int32) f32;

   return 0;
}

/*}}}*/

static int read_int16_to_int32 (Data_Def_Type *ddt) /*{{{*/
{
   int16 i16;

   if (1 != JDMread_int16 (&i16, 1, ddt->ddt_dft->fp))
     return -1;

   *(int32 *) ddt->ddt_value_ptr = (int32) i16;

   return 0;
}

/*}}}*/

static int read_byte_to_int16 (Data_Def_Type *ddt) /*{{{*/
{
   signed char ch;

   if (1 != fread (&ch, 1, 1, ddt->ddt_dft->fp))
     return -1;

   *(int16 *) ddt->ddt_value_ptr = (int16) ch;

   return 0;
}

/*}}}*/

static int read_float32 (Data_Def_Type *ddt) /*{{{*/
{
   if (1 != JDMread_float32 ((float32 *) ddt->ddt_value_ptr, 1, ddt->ddt_dft->fp))
     return -1;

   return 0;
}

/*}}}*/

static int read_float32_to_float64 (Data_Def_Type *ddt) /*{{{*/
{
   float32 f32;

   if (1 != JDMread_float32 (&f32, 1, ddt->ddt_dft->fp))
     return -1;

   *(float64 *) ddt->ddt_value_ptr = (float64) f32;
   return 0;
}

/*}}}*/

/*}}}*/

static int read_expno_value (Data_Def_Type *ddt) /*{{{*/
{
   static int32 last_expno = -1;
   int32 val;

   if (1 != JDMread_int32 (&val, 1, ddt->ddt_dft->fp))
     return -1;

   Data_Table.dtt_update_dither = (last_expno != val);
   last_expno = val;
   *(int32 *) ddt->ddt_value_ptr = val;

   return 0;
}

/*}}}*/

static int read_dither_value (Data_Def_Type *ddt) /*{{{*/
{
   float32 val;

   if (Simulation_Used_Dither == 0)
     val = 0;
   else if (1 != JDMread_float32 (&val, 1, ddt->ddt_dft->fp))
     return -1;

   if ((Data_Table.dtt_update_dither) || (Pixel_Adjust == PIX_ADJ_EXACT))
     *(float32 *) ddt->ddt_value_ptr = val;

   return 0;
}

/*}}}*/

static int compute_tdetxy (Data_Def_Type *ddt) /*{{{*/
{
   unsigned int x, y;

   (void) ddt;

   if (-1 == marx_compute_tiled_pixel (The_Detector, Data_Table.dtt_ccdid,
				       Data_Table.dtt_chipx, Data_Table.dtt_chipy,
				       &x, &y))
     return -1;

   Data_Table.dtt_tdetx = x + 1;
   Data_Table.dtt_tdety = y + 1;
   return 0;
}

/*}}}*/

static Marx_Chip_To_MNC_Type **Chip_To_Mncs;
static int First_Chip_Id;
static int Last_Chip_Id;

static void delete_chip_to_mncs (void)
{
   Marx_Chip_To_MNC_Type **m, **mmax;

   if (Chip_To_Mncs == NULL)
     return;

   m = Chip_To_Mncs;
   mmax = m + (Last_Chip_Id - First_Chip_Id + 1);

   while (m < mmax)
     {
	if (*m != NULL)
	  {
	     marx_free_chip_to_mnc (*m);
	     *m = NULL;
	  }
	m++;
     }
   free ((char *) Chip_To_Mncs);
   Chip_To_Mncs = NULL;
}

static int open_detxy (Data_Def_Type *ddt)
{
   int num, i;

   (void) ddt;
   First_Chip_Id = The_Detector->first_facet_id;
   Last_Chip_Id = The_Detector->last_facet_id;

   num = Last_Chip_Id - First_Chip_Id + 1;
   if (num <= 0)
     {
	fprintf (stderr, "open_detxy: internal error.\n");
	return -1;
     }

   Chip_To_Mncs = (Marx_Chip_To_MNC_Type **) malloc (sizeof (Marx_Chip_To_MNC_Type *) * num);
   if (NULL == Chip_To_Mncs)
     {
	fprintf (stderr, "Not enough memory.\n");
	return -1;
     }
   memset ((char *) Chip_To_Mncs, 0, num * sizeof(Marx_Chip_To_MNC_Type *));

   for (i = First_Chip_Id; i <= Last_Chip_Id; i++)
     {
	Marx_Chip_To_MNC_Type *m;

	if (NULL == (m = marx_allocate_chip_to_mnc (DetectorType, i)))
	  {
	     delete_chip_to_mncs ();
	     return -1;
	  }

	Chip_To_Mncs[i - First_Chip_Id] = m;
     }

   return 0;
}

static int close_detxy (Data_Def_Type *ddt)
{
   (void) ddt;

   delete_chip_to_mncs ();
   return 0;
}

static int compute_detxy (Data_Def_Type *ddt) /*{{{*/
{
   double x, y;
   float dx, dy;
   Marx_Chip_To_MNC_Type *chip2mnc;

   (void) ddt;

   chip2mnc = Chip_To_Mncs[Data_Table.dtt_ccdid - First_Chip_Id];

   if (-1 == marx_init_chip_to_mnc (chip2mnc, Focal_Length,
				    DetOffset_X,
				    DetOffset_Y + Data_Table.dtt_dither.dy,
				    DetOffset_Z + Data_Table.dtt_dither.dz,
				    Data_Table.dtt_dither.dtheta))
     return -1;

   /* Go from a 1-based to a 0-based system */
   x = Data_Table.dtt_chipx - 1.0;
   y = Data_Table.dtt_chipy - 1.0;
   switch (Pixel_Adjust)
     {
      case PIX_ADJ_EXACT:
	break;

      case PIX_ADJ_NONE:
	x = floor(x) + 0.5;
	y = floor(y) + 0.5;
	break;

      case PIX_ADJ_RANDOMIZE:
	x = floor(x) + JDMrandom ();
	y = floor(y) + JDMrandom ();
	break;

      case PIX_ADJ_EDSER:
	if (-1 == marx_compute_acis_subpix (Acis_Subpixel_Object, Data_Table.dtt_ccdid,
					    Data_Table.dtt_benergy, Data_Table.dtt_fltgrade, &dx, &dy))
	  return -1;

	x = floor(x) + 0.5 + dx;
	y = floor(y) + 0.5 + dy;
	break;
     }

   if (-1 == marx_chip_to_mnc (chip2mnc, x, y, &Data_Table.dtt_mnc))
     return -1;

   /* (void) marx_mnc_to_chip (chip2mnc, &mnc, &x, &y); */

   if (-1 == marx_mnc_to_fpc (The_Detector->fp_coord_info, &Data_Table.dtt_mnc,
			      &x, &y))
     return -1;

   /* Note that there is no need to add 0.5 because marx_mnc_to_fpc returns
    * the value already in the FPC system
    */
   Data_Table.dtt_detx = x;
   Data_Table.dtt_dety = y;

   return 0;
}

/*}}}*/

static int compute_expno (Data_Def_Type *ddt) /*{{{*/
{
   static long last_expno = -1;
   long expno;

   (void) ddt;

   if (TimeDel <= 0.0)
     {
	Data_Table.dtt_expno = last_expno;
	last_expno++;
	Data_Table.dtt_update_dither = 1;
	return 0;
     }
   expno = (long)(Data_Table.dtt_time / TimeDel);
   Data_Table.dtt_expno = expno;
   Data_Table.dtt_update_dither = (expno != last_expno);
   last_expno = expno;

   return 0;
}

/*}}}*/

static int compute_node_id (Data_Def_Type *ddt) /*{{{*/
{
   (void) ddt;

   Data_Table.dtt_node_id = (int16) ((Data_Table.dtt_chipx-1)/256);
   return 0;
}

/*}}}*/

static int compute_status (Data_Def_Type *ddt) /*{{{*/
{
   (void) ddt;

   Data_Table.dtt_status = 0;

   if (Simulation_Used_ACIS)
     {
	if (((Data_Table.dtt_chipx<2.0)||(Data_Table.dtt_chipx>=1024.0))
	    || ((Data_Table.dtt_chipy<2.0)||(Data_Table.dtt_chipy>=1024.0)))
	  Data_Table.dtt_status |= 0x0001;
     }

   return 0;
}

/*}}}*/
// see $CALDB/data/chandra/acis/grade/acisD2009-11-01gradeN0005.fits
static char Grade_Map [256] =
{
   0,   1,   2,   5,   1,   1,   5,   7,   3,   5,   6,   6,   3,   5,   7,   7,
   4,   4,   6,   7,   5,   5,   6,   7,  7,   7,   7,   7,   7,   7,   7,   7,
   1,   1,   2,   5,   1,   1,   5,   7,   5,   7,   7,   7,   5,   7,   7,   7,
   4,   4,   6,   7,   5,   5,   6,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   2,   2,  7,   7,   2,   2,   7,   7,   6,   7,   7,   7,   6,   7,   7,   7,
   6,   6,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   5,   5,   7,   7,   5,   5,   7,   7,   6,   7,   7,  7,   6,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   1,   1,   2,   5,   1,   1,   5,   7,   3,   5,   6,   6,   3,   5,   7,   7,
   5,   5,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   1,   1,   2,   5,   1,   1,   5,   7,   5,   7,   7,   7,   5,   7,   7,   7,
   5,   5,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   5,   5,   7,   7,   5,   5,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   6,   6,   7,   7,   7,   7,  7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,  7,
};

static int compute_grade (Data_Def_Type *ddt) /*{{{*/
{
   (void) ddt;
   Data_Table.dtt_grade = Grade_Map[(unsigned char)Data_Table.dtt_fltgrade];
   return 0;
}

/*}}}*/

/* fltgrades assigned to the 3*3 subpixels. Each subpixel had four fltgradea
 * listed here. compute_fltgrade shoses one of them randomly.
 * This array is 9 * 4 instead of 3 * 3 * 4, just to reduce the number of
 * dimensions and simplify typing.
 */
static int Flight_Grade_Table [9][4] =
{
   /* -- */ { 10,  11,  138,  139 },
   /* 0- */ {  2,  34,  130,  162 },
   /* +- */ { 18,  22,   50,   54 },
   /* -0 */ {  8,  12,  136,  140 },
   /* 00 */ {  0,   0,    0,    0 },
   /* +0 */ { 16,  17,   48,   49 },
   /* -+ */ { 72,  76,  104,  108 },
   /* 0+ */ { 64,  65,   68,   69 },
   /* ++ */ { 80,  81,  208,  209 }
};

/* Each pixel is split into 3*3 sub-pixels. Photons hitting the middle one get
 * assigned fltgrade 0 = grade 0, those hitting one of the 4 corners get
 * randomly assigned one of four appropriate corner-split grades (e.g. 10,
 * 11, 138, 139 for the bottom left corner), and those hitting a side one
 * of four vertical or horizontal split grades (e.g. 2, 34, 130, 162 for
 * the middle bottom sub-pixel). That makes almost all fltgrades equally
 * likely, and grade 0 four time more likely than the rest.but, in the words of
 * Cathrine Grant from the ACIS team, any scheme
 * that's not energy dependent and does not distinguish between FI and BI
 * "is just pure fiction". For example, for BI chips that low energies,
 * photons will generate only small charge clouds and even events that hit
 * fairly close to the pixel boundary do not generate split events (see
 * e.g. https://iopscience.iop.org/article/10.1086/421866/fulltext/ - Fig. 1).
 */
static int compute_fltgrade (Data_Def_Type *ddt) /*{{{*/
{
   int dx, dy;
   int *fltgrades;

   (void) ddt;

   dx = (int)(3.0*(Data_Table.dtt_chipx - floor(Data_Table.dtt_chipx)));
   dy = (int)(3.0*(Data_Table.dtt_chipy - floor(Data_Table.dtt_chipy)));
   fltgrades = Flight_Grade_Table[3*dy + dx];
   Data_Table.dtt_fltgrade = fltgrades[(int)(4*JDMrandom())];
   return 0;
}
/*}}}*/

static int compute_xy_sky (Data_Def_Type *ddt)
{
   double x, y;
   double pixel_size;
   Marx_FP_Coord_Type *f;

   (void) ddt;

   f = The_Detector->fp_coord_info;

   if (Simulation_Used_Dither == 0)
     {
	double theta = Nominal_Roll * PI/180.0;
	double c = cos (theta);
	double s = sin (theta);

	x = Data_Table.dtt_detx - f->fp_x0;
	y = Data_Table.dtt_dety - f->fp_y0;

	Data_Table.dtt_xsky = f->fp_x0 + c*x + s*y;
	Data_Table.dtt_ysky = f->fp_y0 - s*x + c*y;
	return 0;
     }

   marx_undither_mnc (&Data_Table.dtt_mnc, &Data_Table.dtt_dither);
   marx_mnc_to_ra_dec (&Data_Table.dtt_mnc, &x, &y);
   /* x, y are in radians.  Now convert them to aspect offsets */
   marx_compute_ra_dec_offsets (0, 0, x, y, &x, &y);

   f = The_Detector->fp_coord_info;
   pixel_size = f->fp_delta_s0;

   /* Note that the ra offset needs to be reversed since RA goes to the
    * left in the sky, and we want X to go to the right.
    */
   x = -x;

   x = x/pixel_size + f->fp_x0;
   y = y/pixel_size + f->fp_y0;

   Data_Table.dtt_xsky = x;
   Data_Table.dtt_ysky = y;
   return 0;
}
#ifdef OBSOLETE_FEATURE
static int compute_pha_island (Data_Def_Type *ddt)
{
   (void) ddt;

   Data_Table.dtt_pha_island [4] = Data_Table.dtt_pha;
   return 0;
}
#endif

static int compute_acis_energy (Data_Def_Type *dtt)
{
   (void) dtt;

   Data_Table.dtt_energy = Data_Table.dtt_benergy * 1e3;
   return 0;
}

static int compute_pi (Data_Def_Type *dtt)
{
   (void) dtt;
#if 1
   Data_Table.dtt_pi = (int16) (Data_Table.dtt_benergy * Acis_PI_Factor * 1e3 + 1.0);
#else
   Data_Table.dtt_pi = (int16) (Data_Table.dtt_benergy * Acis_PI_Factor * 1e3 + 0.5);
#endif
   return 0;
}

#if 0
static int compute_pileup_pha (Data_Def_Type *dtt)
{
   short pha;
   (void) dtt;

   if (-1 == marx_map_energy_to_acis_pha (Data_Table.dtt_ccdid, Data_Table.dtt_chipx, Data_Table.dtt_chipy, Data_Table.dtt_benergy, &pha))
     return -1;

   /* Data_Table.dtt_pha = (int16) (Data_Table.dtt_benergy * Acis_PI_Factor *1e3 + 0.5); */
   Data_Table.dtt_pha = (int16) (pha);
   return 0;
}
#endif

/*}}}*/
