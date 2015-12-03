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

/*{{{ Include Files */

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <sys/types.h>
#include <time.h>

#include <marx.h>
#include <jdfits.h>

/*}}}*/

#include "argcargv.h"
#include "argcargv.c"

#define MAX_PATH_LEN	1024

static char *The_Marx_Dir;
static char *The_Fits_File;
static char *Pgm_Name = "marx2img";

static Marx_Detector_Type *Detector;

static unsigned int First_Chip_Id;
static unsigned int Last_Chip_Id;
static unsigned int Num_X_Pixels;      /* per chip */
static unsigned int Num_Y_Pixels;      /* per chip */

static unsigned int Min_X_Pixel;
static unsigned int Min_Y_Pixel;
static unsigned int Max_X_Pixel;
static unsigned int Max_Y_Pixel;

static unsigned int Pixel_Scale_Factor = 0;

/* If this >= 0, only that chip will be used. */
static int Use_This_Chip = -1;

static float Min_Energy;
static float Max_Energy;
static float Min_Time;
static float Max_Time;
static float Min_BEnergy;
static float Max_BEnergy;

static unsigned int N_Axis_1;	       /* number of x pixels in image */
static unsigned int N_Axis_2;	       /* number of y pixels in image */

static int HRC_Used;

static char *make_marx_filename (char *name)
{
   static char file [MAX_PATH_LEN + 1];
   unsigned int len;

   len = strlen (The_Marx_Dir);

   if (len + 1 + strlen (name) > MAX_PATH_LEN)
     {
	fprintf (stderr, "Filename too long.\n");
	exit (1);
     }

   strcpy (file, The_Marx_Dir);

   if (len && (file[len - 1] != '/'))
     file [len++] = '/';

   strcpy (file + len, name);

   return file;
}

static int compute_tiled_pixel (int chip,
				unsigned int x, unsigned int y,
				unsigned int *xp, unsigned int *yp)
{
   unsigned int new_x, new_y;

   if (-1 == marx_compute_tiled_pixel (Detector, chip, x, y, &new_x, &new_y))
     return -1;

   *xp = new_x / Pixel_Scale_Factor;
   *yp = new_y / Pixel_Scale_Factor;
   return 0;
}

static void get_min_max_pixels (unsigned int *xmin, unsigned int *xmax,
				unsigned int *ymin, unsigned int *ymax)
{
   unsigned int chip;
   unsigned int max_x, max_y, min_x, min_y;
   unsigned int x, y, x_1, y_1;

   /* Since chips may be rotated and oriented in strange configurations,
    * it is not at all clear what the min/max pixel values will be.
    * So, I do a brute force calculation.  Sigh.
    */

   compute_tiled_pixel (First_Chip_Id, 0, 0, &x, &y);
   min_x = max_x = x;
   min_y = max_y = y;

   x_1 = Num_X_Pixels - 1;
   y_1 = Num_Y_Pixels - 1;

   for (chip = First_Chip_Id; chip <= Last_Chip_Id; chip++)
     {
	compute_tiled_pixel (chip, 0, 0, &x, &y);
	if (x > max_x) max_x = x;
	if (x < min_x) min_x = x;
	if (y > max_y) max_y = y;
	if (y < min_y) min_y = y;

	compute_tiled_pixel (chip, x_1, 0, &x, &y);
	if (x > max_x) max_x = x;
	if (x < min_x) min_x = x;
	if (y > max_y) max_y = y;
	if (y < min_y) min_y = y;

	compute_tiled_pixel (chip, 0, y_1, &x, &y);
	if (x > max_x) max_x = x;
	if (x < min_x) min_x = x;
	if (y > max_y) max_y = y;
	if (y < min_y) min_y = y;

	compute_tiled_pixel (chip, x_1, y_1, &x, &y);
	if (x > max_x) max_x = x;
	if (x < min_x) min_x = x;
	if (y > max_y) max_y = y;
	if (y < min_y) min_y = y;
     }

   *xmin = min_x;
   *ymin = min_y;
   *xmax = max_x;
   *ymax = max_y;
}

static int read_marx_par_file (void)
{
   char *marx_par;
   char buf [PF_MAX_LINE_LEN];
   Param_File_Type *pf;
   int detector_type;

   marx_par = make_marx_filename ("marx.par");

   pf = pf_open_parameter_file (marx_par, "r");
   if (pf == NULL)
     {
	fprintf (stderr, "Unable to open parameter file %s.\n",
		 marx_par);
	return -1;
     }

   if (-1 == pf_get_string (pf, "DetectorType", buf, sizeof (buf)))
     {
	pf_close_parameter_file (pf);
	fprintf (stderr, "Perhaps %s is corrupt.\n", marx_par);
	return -1;
     }
   pf_close_parameter_file (pf);

   if (!strcmp (buf, "NONE"))
     {
	pf_close_parameter_file (pf);
	fprintf (stderr, "The simulation used no detector.\n");
	return -1;
     }

   Detector = marx_get_detector_info (buf);
   detector_type = (Detector == NULL) ? -1 : Detector->detector_type;
   switch (detector_type)
     {
      case MARX_DETECTOR_HRC_S:
	Num_X_Pixels = 4096;
	Num_Y_Pixels = 16384;
	break;

      case MARX_DETECTOR_HRC_I:
	Num_X_Pixels = 16384;
	Num_Y_Pixels = 16384;
	break;

      case MARX_DETECTOR_ACIS_S:
      case MARX_DETECTOR_ACIS_I:
	Num_Y_Pixels = 1024;
	Num_X_Pixels = 1024;
	break;

      default:
	  fprintf (stderr, "DetectorType %s not supported.\n", buf);
	return -1;
     }

   First_Chip_Id = Detector->first_facet_id;
   Last_Chip_Id = Detector->last_facet_id;

   if (Use_This_Chip >= 0)
     {
	if (((unsigned int) Use_This_Chip < First_Chip_Id)
	    || ((unsigned int) Use_This_Chip > Last_Chip_Id))
	  {
	     fprintf (stderr, "Chip-id must be between %d and %d",
		      First_Chip_Id, Last_Chip_Id);
	     return -1;
	  }
	First_Chip_Id = Last_Chip_Id = Use_This_Chip;
     }

   fprintf (stdout, "DetectorType: %s\n", buf);
   return 0;
}

static int write_marx_par_comments (JDFits_Type *ft)
{
   char *file = make_marx_filename ("marx.par");
   return jdfits_add_comments_from_file (ft, file,
					"COMMENT", "#@#", 0);
}

static int write_extra_fits_headers (JDFits_Type *ft)
{
   char buf [512];

   if ((-1 == jdfits_write_header_integer (ft, "X-MINX", Min_X_Pixel, NULL))
       || (-1 == jdfits_write_header_integer (ft, "X-MAXX", Max_X_Pixel, NULL))
       || (-1 == jdfits_write_header_integer (ft, "X-MINY", Min_Y_Pixel, NULL))
       || (-1 == jdfits_write_header_integer (ft, "X-MAXY", Max_Y_Pixel, NULL))
       || (-1 == jdfits_write_header_string (ft, "X-MARX", The_Marx_Dir, "MARX directory")))
     return -1;

   sprintf (buf, "This file was generated by %s.", Pgm_Name);
   if (-1 == jdfits_write_header_comment (ft, "COMMENT", buf))
     return -1;
   if (-1 == jdfits_write_header_comment (ft, NULL, "Contact davis@space.mit.edu for more information about the program."))
     return -1;
   if (-1 == write_marx_par_comments (ft))
     return -1;
   return 0;
}

static JDFits_Type *open_fits_file (void)
{
   JDFits_Type *ft;

   if (NULL == (ft = jdfits_open_file (The_Fits_File, JDFITS_WRITE_MODE)))
     return NULL;

   if ((-1 == jdfits_write_header_logical (ft, "SIMPLE", 1, "FITS STANDARD"))
       || (-1 == jdfits_write_header_integer (ft, "BITPIX", 32, "32 bit integer data"))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS", 2, "2d array"))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS1", (int) N_Axis_1, NULL))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS2", (int) N_Axis_2, NULL))
       || (-1 == jdfits_write_header_logical (ft, "EXTEND", 0, "No extensions")))
     {
	jdfits_close_file (ft);
	return NULL;
     }

   if (-1 == write_extra_fits_headers (ft))
     {
	jdfits_close_file (ft);
	return NULL;
     }

   if (-1 == jdfits_end_header (ft))
     {
	jdfits_close_file (ft);
	return NULL;
     }

   return ft;
}

static Marx_Dump_File_Type *open_marx_data_file (char *file)
{
   Marx_Dump_File_Type *dft;

   file = make_marx_filename (file);

   if (NULL == (dft = marx_open_read_dump_file (file)))
     {
	marx_error ("Unable to open %s.", file);
	return NULL;
     }

   return dft;
}

Marx_Dump_File_Type *XPixel_Dft, *YPixel_Dft, *Detector_Dft;
Marx_Dump_File_Type *Energy_Dft, *Time_Dft, *BEnergy_Dft;

static void close_marx_data_files (void)
{
   if (XPixel_Dft != NULL) marx_close_read_dump_file (XPixel_Dft);
   if (YPixel_Dft != NULL) marx_close_read_dump_file (YPixel_Dft);
   if (Detector_Dft != NULL) marx_close_read_dump_file (Detector_Dft);
   if (Energy_Dft != NULL) marx_close_read_dump_file (Energy_Dft);
   if (Time_Dft != NULL) marx_close_read_dump_file (Time_Dft);
   if (BEnergy_Dft != NULL) marx_close_read_dump_file (BEnergy_Dft);

   XPixel_Dft = YPixel_Dft = Detector_Dft = NULL;
   BEnergy_Dft = Energy_Dft = Time_Dft = NULL;
}

static int open_marx_data_files (void)
{
   close_marx_data_files ();

   if ((NULL == (XPixel_Dft = open_marx_data_file ("xpixel.dat")))
       || (NULL == (YPixel_Dft = open_marx_data_file ("ypixel.dat")))
       || (NULL == (Detector_Dft = open_marx_data_file ("detector.dat")))
       || (NULL == (Energy_Dft = open_marx_data_file ("energy.dat")))
       || ((HRC_Used == 0)
	   && (NULL == (BEnergy_Dft = open_marx_data_file ("b_energy.dat"))))
       || (NULL == (Time_Dft = open_marx_data_file ("time.dat"))))
     {
	close_marx_data_files ();
	return -1;
     }

   return 0;
}

static int read_marx_xy (unsigned int *xp, unsigned int *yp, float *energy, float *t, float *b)
{
   unsigned char id;
   float32 e32, t32, x32, y32, b32;

   do
     {
	if (1 != JDMread_float32 (&x32, 1, XPixel_Dft->fp))
	  return 0;
	if (1 != JDMread_float32 (&y32, 1, YPixel_Dft->fp))
	  return -1;
	if (1 != fread ((unsigned char *) &id, 1, 1, Detector_Dft->fp))
	  return -1;

	if (energy != NULL)
	  {
	     if (1 != JDMread_float32 (&e32, 1, Energy_Dft->fp))
	       return -1;
	     *energy = (float) e32;
	  }

	if (t != NULL)
	  {
	     if (1 != JDMread_float32 (&t32, 1, Time_Dft->fp))
	       return -1;
	     *t = (float) t32;
	  }
	if (b != NULL)
	  {
	     if (1 != JDMread_float32 (&b32, 1, BEnergy_Dft->fp))
	       return -1;
	     *b = (float) b32;
	  }
     }
   while (((unsigned int) id < First_Chip_Id)
	  || ((unsigned int) id > Last_Chip_Id));

   if (-1 == compute_tiled_pixel (id,
				  (unsigned int) x32, (unsigned int) y32,
				  xp, yp))
     return -1;

   return 1;
}

#define MAX_BUF_SIZE (1024 * 1024 * 8)

static int write_image (JDFits_Type *ft)
{
   int32 *buf;
   unsigned int max_buffer_size;
   unsigned int dy, dx, dxdy;
   unsigned int y_0, y_1, x_0, x_1, x, y;
   float *ep, *tp, e, t, *bp, b;
   unsigned int count;

   /* Compute buffer size */
   buf = NULL;
   dy = N_Axis_2;
   dx = N_Axis_1;

   while (1)
     {
	if (dy == 0)
	  {
	     fprintf (stderr, "Not enough memory.\n");
	     return -1;
	  }

	dxdy = dx * dy;
	max_buffer_size = 4 * dxdy;
	if ((max_buffer_size <= MAX_BUF_SIZE)
	    && (NULL != (buf = (int32 *) malloc (max_buffer_size))))
	  break;

	dy = dy / 2;
     }

   if (-1 == open_marx_data_files ())
     {
	free ((char *) buf);
	return -1;
     }

   memset ((char *) buf, 0, max_buffer_size);

   y_0 = Min_Y_Pixel; y_1 = y_0 + dy;
   if (y_1 > Max_Y_Pixel)
     y_1 = Max_Y_Pixel;

   x_0 = Min_X_Pixel; x_1 = x_0 + dx;
   if (x_1 > Max_X_Pixel)
     x_1 = Max_X_Pixel;

   if ((Max_Energy < 0.0) && (Min_Energy < 0.0))
     ep = NULL;
   else
     ep = &e;

   if ((Min_Time < 0.0) && (Max_Time < 0.0))
     tp = NULL;
   else
     tp = &t;

   if ((Min_BEnergy < 0) && (Max_BEnergy < 0))
     bp = NULL;
   else
     bp = &b;

   count = 0;
   while (1)
     {
	int status;

	status = read_marx_xy (&x, &y, ep, tp, bp);
	if (status != 1)
	  {
	     if (status == -1)
	       {
		  free ((char *) buf);
		  return -1;
	       }
	     if (-1 == jdfits_write_int32 (ft, buf, (int) dxdy))
	       {
		  free ((char *) buf);
		  fprintf (stderr, "Write error. Check disk space.\n");
		  return -1;
	       }

	     y_0 += dy;
	     if (y_0 >= Max_Y_Pixel)
	       break;

	     if (y_1 + dy > Max_Y_Pixel)
	       {
		  dy = Max_Y_Pixel - y_1;
		  dxdy = dy * dx;
	       }
	     y_1 += dy;

	     if (-1 == open_marx_data_files ())
	       {
		  free ((char *) buf);
		  fprintf (stderr, "Unable to REOPEN marx data files.\n");
		  return -1;
	       }

	     memset ((char *) buf, 0, max_buffer_size);
	     continue;
	  }

	if ((y < y_0) || (y >= y_1)
	    || (x < x_0) || (x >= x_1))
	  continue;

	if ((ep != NULL)
	     && ((e < Min_Energy)
		 || ((Max_Energy >= 0.0) && (e > Max_Energy))))
	  continue;

	if ((tp != NULL)
	     && ((t < Min_Time)
		 || ((Max_Time >= 0.0) && (t > Max_Time))))
	  continue;

	if ((bp != NULL)
	     && ((b < Min_BEnergy)
		 || ((Max_BEnergy >= 0) && (b > Max_BEnergy))))
	  continue;

	*(buf + (y - y_0) * dx + (x - x_0)) += 1;
	count++;
     }

   fprintf (stdout, "%u counts written to the fits file.\n", count);

   free ((char *) buf);
   close_marx_data_files ();
   return jdfits_end_data (ft);
}

static int find_bounding_box (void)
{
   unsigned int x, y;
   unsigned int x_0, x_1, y_0, y_1;
   float e, t;
   float e0, e1, t0, t1;

   fprintf (stdout, "Finding data limits...\n");
   fflush (stdout);

   if (-1 == open_marx_data_files ())
     return -1;

   x_0 = Max_X_Pixel;
   x_1 = Min_X_Pixel;
   y_1 = Min_Y_Pixel;
   y_0 = Max_Y_Pixel;
   t0 = e0 = 1e30;
   t1 = e1 = -1e30;

   while (1 == read_marx_xy (&x, &y, &e, &t, NULL))
     {
	if (x > x_1) x_1 = x;
	if (x < x_0) x_0 = x;
	if (y > y_1) y_1 = y;
	if (y < y_0) y_0 = y;
	if (e > e1) e1 = e;
	if (e < e0) e0 = e;
	if (t > t1) t1 = t;
	if (t < t0) t0 = t;
     }

   close_marx_data_files ();

   if ((x_0 > x_1) || (y_0 > y_1) || (e0 > e1) || (t0 > t1))
     {
	fprintf (stderr, "There are no counts in requested range.\n");
	return -1;
     }

   /*   fprintf (stdout, "\tDetector Energy limits: %f - %f\n", b0, b1); */
   fprintf (stdout, "\tEnergy limits: %f - %f\n", e0, e1);
   fprintf (stdout, "\tTime limits: %f - %f\n\n", t0, t1);
   fflush (stdout);

   Min_X_Pixel = x_0;
   Max_X_Pixel = x_1 + 1;	       /* since min <= x < max  */
   Min_Y_Pixel = y_0;
   Max_Y_Pixel = y_1 + 1;

   return 0;
}

static int Find_BBox;

static ArgcArgv_Type Cmd_Arg_Table [] =
{
   {"--bbox", ARGCARGV_BOOLEAN, (long)&Find_BBox, "Find bounding box."},
   {"--minX", ARGCARGV_INTEGER, (long)&Min_X_Pixel, "<MIN_X_PIXEL>"},
   {"--maxX", ARGCARGV_INTEGER, (long)&Max_X_Pixel, "<MAX_X_PIXEL>"},
   {"--minY", ARGCARGV_INTEGER, (long)&Min_Y_Pixel, "<MIN_Y_PIXEL>"},
   {"--maxY", ARGCARGV_INTEGER, (long)&Max_Y_Pixel, "<MAX_Y_PIXEL>"},
   {"--minE", ARGCARGV_FLOAT, (long)&Min_Energy, "<MIN_ENERGY> (KeV)"},
   {"--maxE", ARGCARGV_FLOAT, (long)&Max_Energy, "<MAX_ENERGY> (KeV)"},
   {"--minT", ARGCARGV_FLOAT, (long)&Min_Time, "<MIN_TIME> (sec)"},
   {"--maxT", ARGCARGV_FLOAT, (long)&Max_Time, "<MAX_TIME> (sec)"},
   {"--minD", ARGCARGV_FLOAT, (long)&Min_BEnergy, "<MIN_DETECTOR_ENERGY>"},
   {"--maxD", ARGCARGV_FLOAT, (long)&Max_BEnergy, "<MAX_DETECTOR_ENERGY>"},
   {"--chip", ARGCARGV_INTEGER, (long) &Use_This_Chip, "<CHIP_NUMBER>"},
   {"--scale", ARGCARGV_INTEGER, (long) &Pixel_Scale_Factor, "<PIXEL_SCALE_FACTOR>"},
   {NULL}
};

static void usage (void)
{
   ArgcArgv_Type *a;

   fprintf (stderr, "Usage: %s [options] <MARXDIR> <FITSFILE>\n",
	    Pgm_Name);
   fprintf (stderr, "Options include:\n");
   a = Cmd_Arg_Table;
   while (a->name != NULL)
     {
	fprintf (stderr, "\t%s\t\t%s\n", a->name, a->help);
	a++;
     }

   exit (1);
}

int main (int argc, char **argv)
{
   JDFits_Type *ft;
   unsigned int xmin, ymin, xmax, ymax;

   Max_X_Pixel = Min_X_Pixel = (unsigned int)-1;
   Max_Y_Pixel = Min_Y_Pixel = (unsigned int)-1;
   Min_Time = -1.0;
   Max_Time = -1.0;
   Min_Energy = -1.0;
   Max_Energy = -1.0;
   Min_BEnergy = -1.0;
   Max_BEnergy = -1.0;

   argc--; argv++;
   if (-1 == argcargv (&argc, &argv, Cmd_Arg_Table))
     usage ();

   if (argc != 2) usage ();
   The_Marx_Dir = argv[0];
   The_Fits_File = argv[1];

   if (-1 == read_marx_par_file ())
     return 1;

   switch (Detector->detector_type)
     {
      case MARX_DETECTOR_HRC_S:
      case MARX_DETECTOR_HRC_I:
	HRC_Used = 1;
	Min_BEnergy = Max_BEnergy = -1.0;
	break;

      default:
	HRC_Used = 0;
	break;
     }

   if (Pixel_Scale_Factor <= 0)
     {
	Pixel_Scale_Factor = 1;
	if ((Detector->detector_type == MARX_DETECTOR_HRC_S)
	    || (Detector->detector_type == MARX_DETECTOR_HRC_I))
	  {
	     Pixel_Scale_Factor = 8;
	     fprintf (stdout, "\n\
***Note: Pixel scale factor not specified.  Since the detector is HRC,\n\
a scale factor of %u will be used.  Set it to 1 via --scale for\n\
full resolution\n\
\n\
", Pixel_Scale_Factor);
	  }
     }

   get_min_max_pixels (&xmin, &xmax, &ymin, &ymax);

   if (Max_X_Pixel == (unsigned int)-1)
     Max_X_Pixel = xmax;

   if (Max_Y_Pixel == (unsigned int)-1)
     Max_Y_Pixel = ymax;

   if (Min_X_Pixel == (unsigned int)-1)
     Min_X_Pixel = xmin;

   if (Min_Y_Pixel == (unsigned int)-1)
     Min_Y_Pixel = ymin;

   if (Min_X_Pixel > Max_X_Pixel)
     {
	xmin = Min_X_Pixel;
	Min_X_Pixel = Max_X_Pixel;
	Max_X_Pixel = xmin;
     }
   if (Min_Y_Pixel > Max_Y_Pixel)
     {
	ymin = Min_Y_Pixel;
	Min_Y_Pixel = Max_Y_Pixel;
	Max_Y_Pixel = ymin;
     }

   if (Find_BBox)
     {
	if (-1 == find_bounding_box ())
	  return 1;
     }

   N_Axis_1 = Max_X_Pixel - Min_X_Pixel;
   N_Axis_2 = Max_Y_Pixel - Min_Y_Pixel;

   fprintf (stdout, "Pixel Size Scale Factor: %u\n", Pixel_Scale_Factor);
   fprintf (stdout, "Num X Pixels: %u (per chip)\n", Num_X_Pixels / Pixel_Scale_Factor);
   fprintf (stdout, "Num Y Pixels: %u (per chip)\n", Num_Y_Pixels / Pixel_Scale_Factor);
   fprintf (stdout, "First Chip Id: %u\n", First_Chip_Id);
   fprintf (stdout, "Last  Chip Id: %u\n", Last_Chip_Id);

   fprintf (stdout, "Min X Pixel: %u\n", Min_X_Pixel);
   fprintf (stdout, "Max X Pixel: %u\n", Max_X_Pixel);
   fprintf (stdout, "Min Y Pixel: %u\n", Min_Y_Pixel);
   fprintf (stdout, "Max Y Pixel: %u\n", Max_Y_Pixel);

   if (Min_Energy >= 0.0) fprintf (stdout, "Min Energy: %f\n", Min_Energy);
   if (Max_Energy >= 0.0) fprintf (stdout, "Max Energy: %f\n", Max_Energy);
   if (Min_Time >= 0.0) fprintf (stdout, "Min Time: %f\n", Min_Time);
   if (Max_Time >= 0.0) fprintf (stdout, "Max Time: %f\n", Max_Time);
   if (Min_BEnergy >= 0) fprintf (stdout, "Min Detector Energy: %f\n", Min_BEnergy);
   if (Max_BEnergy >= 0) fprintf (stdout, "Max Detector Energy: %f\n", Max_BEnergy);

   fprintf (stdout, "Estimated Image size: %lu\n",
	    4 * (unsigned long) N_Axis_1 * (unsigned long) N_Axis_2);

   fflush (stdout);

   if (NULL == (ft = open_fits_file ()))
     return 1;

   if (-1 == write_image (ft))
     {
	jdfits_close_file (ft);
	return 1;
     }

   if (-1 == jdfits_close_file (ft))
     return 1;

   return 0;
}
