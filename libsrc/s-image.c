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

static float *Image;
unsigned int Image_Size;
double Rad_Per_YPixel;
double Rad_Per_XPixel;
unsigned int X_Image_Size;
unsigned int Y_Image_Size;

static int read_image (JDFits_Type *ft, int bitpix, unsigned int n1, unsigned int n2)
{
   double sum;
   unsigned int i;
   unsigned char *read_buf;

   if (-1 == jdfits_read_open_data (ft))
     return -1;

   marx_message ("Allocating space for %u x %u image\n", n1, n2);

   Image_Size = n1 * n2;
   X_Image_Size = n1;
   Y_Image_Size = n2;

   Image = (float *) marx_malloc (Image_Size * sizeof (float));
   if (Image == NULL)
     return -1;

   marx_message ("Reading image data...\n");
   read_buf = NULL;

   switch (bitpix)
     {
      case -32:			       /* float32 */
	for (i = 0; i < Image_Size; i++)
	  {
	     float32 f32;
	     if (1 != jdfits_read_float32 (ft, &f32, 1))
	       goto read_error;
	     Image [i] = (float) f32;
	  }
	break;

      case -64:			       /* double */
	for (i = 0; i < Image_Size; i++)
	  {
	     float64 f64;
	     if (1 != jdfits_read_float64 (ft, &f64, 1))
	       goto read_error;
	     Image [i] = (float) f64;
	  }
	break;

      case 8:
	if (NULL == (read_buf = (unsigned char *) marx_malloc (X_Image_Size)))
	  goto return_error;
	for (i = 0; i < Y_Image_Size; i++)
	  {
	     unsigned int j;
	     float *im;

	     if (X_Image_Size != jdfits_read_bytes (ft, read_buf, X_Image_Size))
	       goto read_error;

	     im = Image + i * X_Image_Size;
	     for (j = 0; j < X_Image_Size; j++)
	       im [j] = (float) read_buf[j];
	  }
	marx_free ((char *)read_buf);
	read_buf = NULL;
	break;

      case 16:
	for (i = 0; i < Image_Size; i++)
	  {
	     int16 i16;
	     if (1 != jdfits_read_int16 (ft, &i16, 1))
	       goto read_error;
	     Image [i] = (float) i16;
	  }
	break;

      case 32:
	for (i = 0; i < Image_Size; i++)
	  {
	     int32 i32;
	     if (1 != jdfits_read_int32 (ft, &i32, 1))
	       goto read_error;
	     Image [i] = (float) i32;
	  }
	break;
      default:
	marx_error ("BITPIX = %d is not supported", bitpix);
	goto return_error;
     }

   /* compute a cumulative distribution */
   sum = 0;
   for (i = 0; i < Image_Size; i++)
     {
        sum += (double) Image [i];
        Image [i] = sum;
     }

   /* Normalize it. */
   for (i = 0; i < Image_Size; i++)
     Image [i] = Image[i] / sum;

   return 0;

   read_error:
   marx_error ("Error reading image file");
   /* drop */
   return_error:
   if (read_buf != NULL) marx_free ((char *)read_buf);
   marx_free ((char *) Image);
   return -1;
}

static int read_keyword_int (JDFits_Type *ft, char *kw, int *i, int is_error)
{
   JDFits_Keyword_Type *k;

   if (NULL == (k = jdfits_find_keyword (ft, kw)))
     {
	if (is_error)
	  marx_error ("Error locating keyword %s", kw);
	return -1;
     }

   if (0 == jdfits_extract_integer (k, i))
     return 0;

   if (is_error)
     marx_error ("Keyword %s is not an integer", kw);

   return -1;
}

static int read_keyword_double (JDFits_Type *ft, char *kw, double *d, int is_error)
{
   JDFits_Keyword_Type *k;

   if (NULL == (k = jdfits_find_keyword (ft, kw)))
     {
	if (is_error)
	  marx_error ("Error locating keyword %s", kw);
	return -1;
     }

   if (0 == jdfits_extract_double (k, d))
     return 0;

   if (is_error)
     marx_error ("Keyword %s is not a double", kw);

   return -1;
}

static int open_image_fits_file (char *file)
{
   JDFits_Type *ft;
   int naxis1, naxis2;
   double scale1, scale2;
   int bitpix;

   marx_message ("Opening IMAGE fits file %s\n", file);
   if (NULL == (ft = jdfits_open_file (file, JDFITS_READ_MODE)))
     {
	marx_error ("Failed to open %s", file);
	return -1;
     }

   if ((-1 == read_keyword_int (ft, "NAXIS1", &naxis1, 1))
       || (-1 == read_keyword_int (ft, "NAXIS2", &naxis2, 1))
       || (-1 == read_keyword_int (ft, "BITPIX", &bitpix, 1)))
     {
	jdfits_close_file (ft);
	return -1;
     }

   /* Look for scale factors */
   if (0 == read_keyword_double (ft, "CDELT1", &scale1, 0))
     {
	if (-1 == read_keyword_double (ft, "CDELT2", &scale2, 1))
	  {
	     jdfits_close_file (ft);
	     return -1;
	  }
	scale1 = fabs (scale1);
	scale2 = fabs (scale2);
     }
   else if (0 == read_keyword_double (ft, "CD1_1", &scale1, 0))
     {
	double cd1_1, cd1_2, cd2_1, cd2_2;

	cd1_1 = scale1;

	if ((-1 == read_keyword_double (ft, "CD1_2", &cd1_2, 1))
	    || (-1 == read_keyword_double (ft, "CD2_1", &cd2_1, 1))
	    || (-1 == read_keyword_double (ft, "CD2_2", &cd2_2, 1)))
	  {
	     jdfits_close_file (ft);
	     return -1;
	  }
	scale1 = sqrt (cd1_1 * cd1_1 + cd1_2 * cd1_2);
	scale2 = sqrt (cd2_1 * cd2_1 + cd2_2 * cd2_2);
     }
   else scale1 = scale2 = 0.0;

   if ((scale1 == 0.0) || (scale2 == 0.0))
     {
	marx_error ("Cannot determine scale of image.  Assuming 0.5 arc-sec per pixel");
	scale1 = scale2 = 0.5 / 3600.0;
     }

   Rad_Per_XPixel = scale1 * PI / 180.0;
   Rad_Per_YPixel = scale2 * PI / 180.0;

   if (-1 == read_image (ft, bitpix, (unsigned int) naxis1, (unsigned int) naxis2))
     {
	jdfits_close_file (ft);
	return -1;
     }

   jdfits_close_file (ft);
   return 0;
}

static int image_open_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int image_close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;
   return 0;
}

/*}}}*/

static int image_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				 unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   Marx_Photon_Attr_Type *at;
   int (*efun) (Marx_Spectrum_Type *, double *);
   JDMVector_Type p0, normal;
   double theta;

   at = pt->attributes;
   efun = st->spectrum.energy_function;

   /* p0 points from the desired center of the image to the origin.  The
    * image center is by definition (-1,0,0).  So, I need to rotate about
    * an axis that rotates (-1,0,0) into p0.
    */

   p0 = st->p;
   theta = JDMv_dot_prod (p0, JDMv_vector (-1, 0, 0));
   /* handle roundoff-error */
   if (fabs (theta) > 1.0)
     {
	if (theta < 0) theta = -1.0;
	else theta = 1.0;
     }
   theta = acos (theta);

   /* Note: normal could be the NULL vector if the vectors coincide.
    * This is ok since theta would be 0, which implies no rotation.
    */
   normal = JDMv_cross_prod (JDMv_vector (-1, 0, 0), p0);
   if (theta != 0.0)
     JDMv_normalize (&normal);

   for (i = 0; i < num; i++)
     {
	unsigned int ofs;
	double x, y, cos_y;
	JDMVector_Type p1;

	if (-1 == (*efun) (&st->spectrum, &at->energy))
	  return -1;

	ofs = JDMbinary_search_f (JDMrandom (), Image, Image_Size);

	y = (double) (ofs / X_Image_Size);
	x = (double) (ofs % X_Image_Size);

	/* Center the image and randomize within the pixel.  Note that the
	 * average of the RHS of next expression is 0.5*YSize
	 */
	y -= 0.5 * (Y_Image_Size - 1) + JDMrandom ();
	x -= 0.5 * (X_Image_Size - 1) + JDMrandom ();

	y = y * Rad_Per_YPixel;
	x = x * Rad_Per_XPixel;

	cos_y = cos (y);

	/* Note that an image on the sky has +X running in -Marx_Y, and
	 * +Y in +Marx_Z.  So, when flipping the vector to point to the
	 * origin, take these into account.
	 */
	p1.x = -cos_y * cos (x);
	p1.y = cos_y * sin (x);	       /* <-- no flip here */
	p1.z = -sin (y);

	at->p = JDMv_rotate_unit_vector (p1, normal, theta);

	at++;
     }

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_image_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   char buf [PF_MAX_LINE_LEN];

   (void) source_id;

   st->open_source = image_open_source;
   st->create_photons = image_create_photons;
   st->close_source = image_close_source;

   if (-1 == pf_get_file (p, "S-ImageFile", buf, sizeof (buf)))
     return -1;

   if (-1 == open_image_fits_file (buf))
     return -1;

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

