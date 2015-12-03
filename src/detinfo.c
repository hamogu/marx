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

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <marx.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

static char *Pgm_Name = "detinfo";
static double Focal_Length = 10065.5;

static void usage (void)
{
   fprintf (stderr, "Usage: %s [--sky] ACIS-S|ACIS-I|HRC-S|HRC-I\n", Pgm_Name);
   exit (1);
}

static void dump_physical_geometry (Marx_Detector_Type *);
static void dump_sky_geometry (Marx_Detector_Type *, char *);

int main (int argc, char **argv)
{
   char *name = NULL;
   Marx_Detector_Type *det;
   int sky = 0;

   switch (argc)
     {
      case 2:
	name = argv[1];
	break;

      case 3:
	if (0 == strcmp ("--sky", argv[1]))
	  {
	     name = argv[2];
	     sky = 1;
	     break;
	  }
	/* drop */
      default:
	usage ();
     }

   if (NULL != strstr(name, ".par"))
     {
	Param_File_Type *p;
	static char namebuf[1024];

	if (NULL == (p = pf_open_parameter_file (name, "r")))
	  return 1;

	if (-1 == marx_mirror_init (p))
	  return 1;

	if (-1 == marx_grating_init (p))
	  return 1;

	if (-1 == marx_detector_init (p))
	  return 1;

	if ((-1 == pf_get_string (p, "DetectorType", namebuf, sizeof(namebuf)))
	    || (-1 == pf_get_double (p, "FocalLength", &Focal_Length)))
	  {
	     pf_close_parameter_file (p);
	     return 1;
	  }

	pf_close_parameter_file (p);
	name = namebuf;
     }

   if (NULL == (det = marx_get_detector_info (name)))
     usage ();

   fprintf (stdout, " Detector Name: %s\n", name);
   fprintf (stdout, " First Chip id: %d\n", det->first_facet_id);
   fprintf (stdout, "  Last Chip id: %d\n", det->last_facet_id);

   if (det->print_info != NULL)
     (void) (*det->print_info)(det, stdout);

   if (sky == 0)
     dump_physical_geometry (det);
   else
     dump_sky_geometry (det, name);

   return 0;
}

static void compute_ra_dec (Marx_Chip_To_MNC_Type *chip_mnc,
			    unsigned int x, unsigned int y,
			    double *ra, double *dec,
			    double *tx, double *ty)
{
   JDMVector_Type p;
   JDMVector_Type p0;

   (void) marx_chip_to_mnc (chip_mnc, x, y, &p);
   marx_mnc_to_ra_dec (&p, ra, dec);

   /* Convert to arc-sec from radians */
   *ra = *ra * (60.0 * 180.0/PI);
   *dec = *dec * (60.0 * 180.0/PI);

   p0.x = 1.0;
   p0.y = p0.z = 0;

   /* Flip sign of MNC */
   p.x = -p.x;
   p.y = -p.y;
   p.z = -p.z;

   (void) marx_vector_to_tan_plane (&p, &p0, tx, ty);
   *tx = *tx * (60.0 * 180.0/PI);
   *ty = *ty * (60.0 * 180.0/PI);
}

static void dump_sky_geometry (Marx_Detector_Type *det, char *name)
{
   int i;
   unsigned int min_x_pixel, max_x_pixel;
   unsigned int min_y_pixel, max_y_pixel;
   Marx_Detector_Geometry_Type *g;

   fprintf (stdout, "\nSky RA/Dec Geometry Follows (with respect to nominal)");

   min_x_pixel = 0;
   min_y_pixel = 0;

   g = det->facet_list;

   while (g != NULL)
     {
	Marx_Chip_To_MNC_Type *chip_mnc;
	double ra, dec;
	double tx, ty;

	i = g->id;
	max_x_pixel = g->num_x_pixels;
	max_y_pixel = g->num_y_pixels;

	fprintf (stdout, "\nChip %d:\n", i);
	fprintf (stdout, " Location of corners: (RA, Dec, RA--TAN, Dec--TAN (arc-min))\n");

	if (NULL == (chip_mnc = marx_allocate_chip_to_mnc (name, i)))
	  exit (1);

	if (-1 == marx_init_chip_to_mnc (chip_mnc, Focal_Length, 0, 0, 0, 0))
	  exit (1);

	compute_ra_dec (chip_mnc, min_x_pixel, min_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e % 16.10e\t(LL)\n", ra, dec, tx, ty);

	compute_ra_dec (chip_mnc, max_x_pixel, min_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e % 16.10e\t(LR)\n", ra, dec, tx, ty);

	compute_ra_dec (chip_mnc, max_x_pixel, max_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e % 16.10e\t(UR)\n", ra, dec, tx, ty);

	compute_ra_dec (chip_mnc, min_x_pixel, max_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e % 16.10e\t(UL)\n", ra, dec, tx, ty);

	marx_free_chip_to_mnc (chip_mnc);
	g = g->next;
     }
}

static void dump_physical_geometry (Marx_Detector_Type *det)
{
   Marx_Detector_Geometry_Type *g;

   fprintf (stdout, "\nPhysical Geometry Follows (STF coords at nominal aimpoint, units in mm)\n");

   g = det->facet_list;

   while (g != NULL)
     {
	fprintf (stdout, "\nChip %d:\n", g->id);
	fprintf (stdout, " Num X Pixels: %u (per chip)\n", g->num_x_pixels);
	fprintf (stdout, " Num Y Pixels: %u (per chip)\n", g->num_y_pixels);
	fprintf (stdout, " X Pixel Size: % 16.10e (mm)\n", g->x_pixel_size);
	fprintf (stdout, " Y Pixel Size: % 16.10e (mm)\n", g->y_pixel_size);
	fprintf (stdout, " X Length: % 16.10e (pixel-size: % 16.10e)\n", g->xlen, g->x_pixel_size);
	fprintf (stdout, " Y Length: % 16.10e (pixel-size: % 16.10e)\n", g->ylen, g->y_pixel_size);
	fprintf (stdout, " Location of corners:\n");
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e\t(LL)\n",
		 g->x_ll.x, g->x_ll.y, g->x_ll.z);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e\t(LR)\n",
		 g->x_lr.x, g->x_lr.y, g->x_lr.z);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e\t(UR)\n",
		 g->x_ur.x, g->x_ur.y, g->x_ur.z);
	fprintf (stdout, "\t% 16.10e % 16.10e % 16.10e\t(UL)\n",
		 g->x_ul.x, g->x_ul.y, g->x_ul.z);

	g = g->next;
     }
}

