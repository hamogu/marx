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
#include <string.h>
#include <jdmath.h>
#include <jdfits.h>

#include "chandra.h"

#include "marx.h"
#include "_marx.h"

static char *Caldb_Index_File = "marxcaldb.par";

static char *Geom_File;
static char *Aimp_File;
#if MARX_HAS_ACIS_GAIN_MAP
static char *Gain_File;
#endif
static char *ACIS_QE_File;
static char *ACIS_Fef_File;
static char *Hetg_Geff_File;
static char *Letg_Geff_File;
static char *ACIS_Contam_File;
static char *ACIS_Subpix_File;

static Param_Table_Type Caldb_Parm_Table [] =
{
   {"GEOM",	PF_FILE_TYPE,		&Geom_File},
   {"AIMPTS",	PF_FILE_TYPE,		&Aimp_File},
#if MARX_HAS_ACIS_GAIN_MAP
   {"GAIN",	PF_FILE_TYPE,		&Gain_File},
#endif
   {"ACISQE",	PF_FILE_TYPE,		&ACIS_QE_File},
   {"ACISFEF",	PF_FILE_TYPE,		&ACIS_Fef_File},
   {"HETGEFF",	PF_FILE_TYPE,		&Hetg_Geff_File},
   {"LETGEFF",	PF_FILE_TYPE,		&Letg_Geff_File},
   {"ACISCONTAM",PF_FILE_TYPE,		&ACIS_Contam_File},
   {"ACISSUBPIX",PF_FILE_TYPE,		&ACIS_Subpix_File},
   {NULL, 0, NULL}
};

static int Caldb_Inited = 0;

static int init_caldb (void)
{
   char *file;
   Param_File_Type *pf;

   if (Caldb_Inited)
     return 0;

   if (NULL == (file = marx_make_data_file_name (Caldb_Index_File)))
     return -1;

   pf = pf_open_parameter_file (file, "rQ");
   if (pf == NULL)
     {
	marx_error ("Unable to open %s", file);
	marx_free (file);
	return -1;
     }
   if (-1 == pf_get_parameters (pf, Caldb_Parm_Table))
     {
	marx_error ("Error processing %s", file);
	marx_free (file);
	(void) pf_close_parameter_file (pf);

	return -1;
     }
   marx_free (file);
   Caldb_Inited = 1;
   return 0;
}

/** Return the filename and path of a caldb data product from marxcaldb.par
  * MARX includes a few CalDB files. The filename to be used is defined 
  * in marxcaldb.par. This function looks for a key in that parameter file
  * and returns the filename (e.g. _marx_caldb_get_file("GEOM")
  * The return value includes the full, absolute path to the file.
  */
char *_marx_caldb_get_file (char *object)
{
   Param_Table_Type *t;

   if (-1 == init_caldb ())
     return NULL;

   t = Caldb_Parm_Table;
   while (t->name != NULL)
     {
	char *file;

	if (0 != strcmp (t->name, object))
	  {
	     t++;
	     continue;
	  }

	file = *(char **) t->value;
	if (NULL == (file = marx_make_data_file_name (file)))
	  return NULL;

	return file;
     }

   marx_error ("Unable to find %s caldb file", object);
   return NULL;
}

/** Return the filename of a caldb data product from marxcaldb.par
  * MARX includes a few CalDB files. The filename to be used is defined 
  * in marxcaldb.par. This function looks for a key in that parameter file
  * and returns the filename, e.g. marx_caldb_get_filename("GEOM") .
  */
char *marx_caldb_get_filename (char *object)
{
   Param_Table_Type *t;

   if (-1 == init_caldb ())
     return NULL;

   t = Caldb_Parm_Table;
   while (t->name != NULL)
     {
	char *file;

	if (0 != strcmp (t->name, object))
	  {
	     t++;
	     continue;
	  }

	file = *(char **) t->value;
	return file;
     }

   marx_error ("Unable to find %s caldb file", object);
   return NULL;
}


static int is_keyword_str (JDFits_Type *ft, char *key, char *val)
{
   char *s;
   int status;

   if (-1 == jdfits_read_keyword_string (ft, key, &s))
     return 0;

   status = (0 == jdfits_strcasecmp (s, val));
   marx_free (s);
   return status;
}

static int is_keyword_int (JDFits_Type *ft, char *key, int val)
{
   int val1;

   if (-1 == jdfits_read_keyword_int (ft, key, &val1))
     return 0;

   return (val == val1);
}

static int check_repeat (JDFits_Col_Data_Type *c, unsigned int repeat, int n)
{
   c += n;
   n++;				       /* columns are 1-based */
   if (c->repeat != repeat)
     {
	marx_error ("Expecting a repeat count of 1 for column %d", n);
	return -1;
     }
   return 0;
}

#if 0
static void tweak_chip (Marx_Detector_Geometry_Type *g)
{
   JDMVector_Type axis;
   JDMVector_Type v;
   double theta = -3 * PI/180.0;
   double cos_theta, sin_theta;

   if (g->id < 4)
     return;

   cos_theta = cos(theta);
   sin_theta = sin(theta);

   axis = JDMv_unit_vector (JDMv_diff(g->x_ul, g->x_ll));

   v = JDMv_rotate_vector1 (JDMv_diff (g->x_lr, g->x_ll), axis, cos_theta, sin_theta);
   g->x_lr = JDMv_sum (g->x_ll, v);

   v = JDMv_rotate_vector1 (JDMv_diff (g->x_ur, g->x_ul), axis, cos_theta, sin_theta);
   g->x_ur = JDMv_sum (g->x_ul, v);
}
#endif
static void set_vector (JDMVector_Type *v, double *xyz)
{
   v->x = xyz[0];
   v->y = xyz[1];
   v->z = xyz[2];
}

static int set_corner_data (Marx_Detector_Geometry_Type *g,
			    JDFits_Col_Data_Type *c)
{
   int chip_id;

   chip_id = c[0].data.i[0];

   while (g != NULL)
     {
	if (g->id == chip_id)
	  {
	     set_vector (&g->x_ll, c[1].data.d);
	     set_vector (&g->x_lr, c[2].data.d);
	     set_vector (&g->x_ul, c[3].data.d);
	     set_vector (&g->x_ur, c[4].data.d);
	     /* tweak_chip (g); */
	     return 0;
	  }
	g = g->next;
     }

   return 0;
}

static int patch_detector_geom (Marx_Detector_Type *d,
				int (*ext_fun) (void *, JDFits_Type *))
{
   char *file;
   JDFits_Type *ft;
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   Marx_Detector_Geometry_Type *g;

   if (NULL == (file = _marx_caldb_get_file ("GEOM")))
     return -1;

   ft = jdfits_find_binary_table (file, ext_fun, NULL);
   if (ft == NULL)
     {
	marx_error ("Unable to find the proper INSTGEOM extension in %s", file);
	marx_free (file);
	return -1;
     }
   r = jdfits_bintable_open_rows (ft, 5, "i:CHIP_ID", "d:LL", "d:LR", "d:UL", "d:UR");
   if (r == NULL)
     goto return_error;

   c = r->col_data;
   if ((-1 == check_repeat (c, 1, 0))
       || (-1 == check_repeat (c, 3, 1))
       || (-1 == check_repeat (c, 3, 2))
       || (-1 == check_repeat (c, 3, 3))
       || (-1 == check_repeat (c, 3, 4)))
     goto return_error;

   g = d->facet_list;
   while (1 == jdfits_read_next_row (ft, r))
     {
	if (-1 == set_corner_data (g, c))
	  goto return_error;
     }

   marx_free (file);
   jdfits_bintable_close_rows (r);
   jdfits_close_file (ft);
   return 0;

   return_error:
   marx_error ("Error processing %s", file);
   marx_free (file);
   jdfits_bintable_close_rows (r);
   (void) jdfits_close_file (ft);
   return -1;
}

static int ext_instgeom_acis_fun (void *unused, JDFits_Type *ft)
{
   (void) unused;
   if (is_keyword_str(ft, "EXTNAME", "INSTGEOM")
       && is_keyword_str (ft, "INSTRUME", "ACIS"))
     return 1;

   return 0;
}
int _marx_caldb_patch_acis_geom (Marx_Detector_Type *d)
{
   return patch_detector_geom (d, ext_instgeom_acis_fun);
}

static int ext_instgeom_hrc_s_fun (void *unused, JDFits_Type *ft)
{
   (void) unused;
   if (is_keyword_str(ft, "EXTNAME", "INSTGEOM")
       && is_keyword_str (ft, "INSTRUME", "HRC-S"))
     return 1;

   return 0;
}

int _marx_caldb_patch_hrc_s_geom (Marx_Detector_Type *d)
{
   return patch_detector_geom (d, ext_instgeom_hrc_s_fun);
}

/* Unfortunately character strings read from the caldb files have trailing
 * whitespace.  Deal with that non-sense here. Sigh.
 */
static void remove_trailing_whitespace (char *s)
{
   char *e;

   e = s + strlen (s);
   while (e > s)
     {
	e--;
	if (*e != ' ')
	  break;
	*e = 0;
     }
}

static int read_3d_from_named_row (char *object, char *extnam, char *columns[2],
				   char *name, JDMVector_Type *v)
{
   JDFits_Type *ft;
   JDFits_Col_Data_Type *c;
   JDFits_Row_Type *r;
   char *file;

   file = _marx_caldb_get_file (object);

   if (file == NULL)
     return -1;

   ft = jdfits_open_binary_table (file, extnam);
   if (ft == NULL)
     {
	marx_free (file);
	return -1;
     }

   if (NULL == (r = jdfits_bintable_aopen_rows (ft, 2, columns)))
     goto return_error;

   c = r->col_data;
   if (-1 == check_repeat (c, 3, 1))
     goto return_error;

   while (1 == jdfits_read_next_row (ft, r))
     {
	remove_trailing_whitespace (c[0].data.a);

	if (0 != jdfits_strcasecmp (c[0].data.a, name))
	  continue;

	set_vector (v, c[1].data.d);
	marx_free (file);
	jdfits_bintable_close_rows (r);
	jdfits_close_file (ft);
	return 0;
     }

   /* drop */

   return_error:
   marx_error ("Unable to find %s %s data in %s", name, columns[0]+2, file);
   marx_free (file);
   jdfits_bintable_close_rows (r);
   jdfits_close_file (ft);
   return -1;
}

static int patch_aimpoint (Marx_Detector_Type *d)
{
   char *aimpt;
   static char *columns[2] =
     {
	"a:AIMPOINT_NAME", "d:AIMPOINT"
     };

   switch (d->detector_type)
     {
      case MARX_DETECTOR_HRC_S:
	aimpt = "HS1";
	break;

      case MARX_DETECTOR_HRC_I:
	aimpt = "HI1";
	break;

      case MARX_DETECTOR_ACIS_S:
	aimpt = "AS1";
	break;

      case MARX_DETECTOR_ACIS_I:
	aimpt = "AI2";
	break;

      default:
	return 0;
     }

   return read_3d_from_named_row ("AIMPTS", "AIMPOINTS", columns, aimpt, &d->stf_stt_offset);
}

static int patch_origin (Marx_Detector_Type *d)
{
   char *inst;
   static char *columns[2] =
     {
	"a:INSTRUMENT", "d:ORIGIN"
     };

   switch (d->detector_type)
     {
      case MARX_DETECTOR_HRC_S:
	inst = "HRC-S";
	break;

      case MARX_DETECTOR_HRC_I:
	inst = "HRC-I";
	break;

      case MARX_DETECTOR_ACIS_I:
      case MARX_DETECTOR_ACIS_S:
	inst = "ACIS";
	break;

      default:
	return 0;
     }

   return read_3d_from_named_row ("GEOM", "INSTRUMENTS", columns, inst,
				  &d->stt_lsi_offset);
}

int _marx_caldb_patch_aimpoint (Marx_Detector_Type *d)
{
   Marx_Detector_Geometry_Type *g;
   JDMVector_Type ofs;

   if ((-1 == patch_origin (d))
       || (-1 == patch_aimpoint (d)))
     return -1;

   d->aimpoint_offset = d->stf_stt_offset;
   ofs = JDMv_sum (d->stf_stt_offset, d->stt_lsi_offset);

   g = d->facet_list;
   while (g != NULL)
     {
	g->x_lr = JDMv_sum (g->x_lr, ofs);
	g->x_ll = JDMv_sum (g->x_ll, ofs);
	g->x_ur = JDMv_sum (g->x_ur, ofs);
	g->x_ul = JDMv_sum (g->x_ul, ofs);

	g = g->next;
     }
   return 0;
}

static int ext_acis_qe_fun (void *vccdid, JDFits_Type *ft)
{
   int ccdid = *(int *)vccdid;

   if (is_keyword_str(ft, "EXTNAME", "AXAF_QE")
       && is_keyword_int (ft, "CCD_ID", ccdid))
     return 1;

   return 0;
}

int _marx_read_acis_qe (int ccdid, float **enp, float **qep, unsigned int *np)
{
   JDFits_Type *ft;
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   char *file;
   unsigned int num_rows;
   unsigned int i;
   float *en, *qe;

   en = qe = NULL;

   if (NULL == (file = _marx_caldb_get_file ("ACISQE")))
     return -1;

   marx_message ("\t%s for [CCDID = %d]\n", file, ccdid);
   ft = jdfits_find_binary_table (file, ext_acis_qe_fun, (void *)&ccdid);
   if (ft == NULL)
     {
	marx_error ("Unable to find ACIS-%d QE extension in %s", ccdid, file);
	marx_free (file);
	return -1;
     }
   r = jdfits_bintable_open_rows (ft, 2, "f:ENERGY", "f:QE");
   if (r == NULL)
     goto return_error;

   c = r->col_data;
   if (-1 == check_repeat (c, 1, 0)
       || (-1 == check_repeat (c, 1, 1)))
     goto return_error;

   num_rows = r->num_rows;
   if (num_rows < 2)
     {
	marx_error ("Expecting more than 2 rows in %s for CCD_ID=%d", file, ccdid);
	goto return_error;
     }

   if ((NULL == (en = (float *) marx_malloc (num_rows * sizeof (float))))
       || (NULL == (qe = (float *) marx_malloc (num_rows * sizeof (float)))))
     goto return_error;

   if (en == NULL)
     goto return_error;

   for (i = 0; i < num_rows; i++)
     {
	if (1 != jdfits_read_next_row (ft, r))
	  {
	     marx_error ("File %s appears to be corrupt", file);
	     goto return_error;
	  }
	en[i] = c[0].data.f[0];
	qe[i] = c[1].data.f[0];
     }
   *enp = en;
   *qep = qe;
   *np = num_rows;
   marx_free (file);
   jdfits_bintable_close_rows (r);
   (void) jdfits_close_file (ft);
   return 0;

   return_error:
   marx_free ((char *) qe);
   marx_free ((char *) en);
   marx_error ("Error processing %s", file);
   marx_free (file);
   jdfits_bintable_close_rows (r);
   (void) jdfits_close_file (ft);
   return -1;
}

static int ext_hduname_fun (void *vhdunam, JDFits_Type *f)
{
   char *hduname = (char *) vhdunam;
   if (is_keyword_str (f, "HDUNAME", hduname))
     return 1;

   return 0;
}

JDFits_Type *_marx_open_binary_hdu (char *file, char *hduname)
{
   JDFits_Type *f;

   f = jdfits_find_binary_table (file, ext_hduname_fun, (void *) hduname);
   if (f == NULL)
     marx_error ("Unable to find an extension with HDUNAME=%s in %s", hduname, file);
   return f;
}
