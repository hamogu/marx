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

static char *Pgm_Name = "sky2chip";
static double Focal_Length = 10065.5; /* 10079.0; */

static void usage (void)
{
   fprintf (stderr, "Usage: %s ACIS-S|ACIS-I|HRC-S|HRC-I  RA(arc-min) DEC(arc-min)\n", Pgm_Name);
   exit (1);
}


/*  A simple command line tool to convert sky positions to XY on chip.
 *
 *  RA and DEC are given in arcmin.
 *  It is not clear to me how Ra/Dec alone are sufficient to calculate
 *  X, Y. The pointig direction and roll angle are also important.
 *  I assume, this program putsi n some standard value.
 *  It was probably written to help with testing something.
 */
int main (int argc, char **argv)
{
   char *name;
   Marx_Detector_Type *det;
   double ra, dec;
   int i;
   JDMVector_Type mnc;

   if (argc != 4)
     usage ();

   name = argv[1];
   if ((1 != sscanf (argv[2], "%lf", &ra))
       || (1 != sscanf (argv[3], "%lf", &dec)))
     usage ();

   ra = PI/(180.0 * 60) * ra;
   dec = PI/(180.0 * 60) * dec;

   if (NULL == (det = marx_get_detector_info (name)))
     usage ();

   mnc.x = -cos (dec) * cos (ra);
   mnc.y = -cos (dec) * sin (ra);
   mnc.z = -sin (dec);

   for (i = det->first_facet_id; i <= det->last_facet_id; i++)
     {
	Marx_Chip_To_MNC_Type *chip_mnc;
	double x, y;

	fprintf (stdout, "%s Chip %d: ", name, i);

	if (NULL == (chip_mnc = marx_allocate_chip_to_mnc (name, i)))
	  exit (1);

	if (-1 == marx_init_chip_to_mnc (chip_mnc, Focal_Length, 0, 0, 0, 0))
	  exit (1);

	marx_mnc_to_chip (chip_mnc, &mnc, &x, &y);

	/* x and y are in pixels with (0,0) at corner of pixel.  SAO wants
	 * pixels with (0.5,0.5) at corner.  So, adjust
	 */
	x += 0.5;
	y += 0.5;
	fprintf (stdout, "\t% 8.2f\t%8.2f\n", x, y);

	marx_free_chip_to_mnc (chip_mnc);
     }

   return 0;
}
