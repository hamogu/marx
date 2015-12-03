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
/* Create an SAOSAC dpde file.
 */

#include "config.h"
#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <jdmath.h>

#include "marx.h"

/* This file is based on the following email: */

/* To: jiahong@head-cfa.harvard.edu, wise
 * Cc: edgar@head-cfa.harvard.edu, dj@head-cfa.harvard.edu
 * Subject: dpde format
 * Date: Mon, 01 Apr 1996 16:48:17 -0500
 * From: "Richard J. Edgar" <edgar@head-cfa.harvard.edu>
 * Status: RO
 *
 *
 * Mike & Jiahong,
 *
 * Here's the include file we use with all these "phot" photons,
 * which gives the C structures in question.  The one we want is
 * PhotonDPDE.  This file is in /proj/axaf/Simul/include/photon.h
 * on the CfA HEAD computer system.  The <fullray.h> file it
 * needs is in the same directory.  The energies are in keV, and
 * the coordinates are in mm, with our coordinate system;
 * z increases along the optical axis from the HRMA towards
 * the focus, with zero at a point of your choice (we typically
 * use 1000 mm in front of the front of the CAP); y is "up" in
 * the XRCF configuration, and x completes a right-handed coordinate
 * system.  x=y=0 on the optical axis of the telescope.
 * "wt" is the photon weight, from zero to one.
 */

/* This file consists of a header followed by structure defining a photon */

typedef struct /*{{{*/
{
   float32 exposure;
   int32 n;
   int32 ftype;			       /* = 12 */
}

/*}}}*/
DPDE_Header_Type;

typedef struct /*{{{*/
{
   float64 x, y, z;
   float64 dx, dy, dz;
   float64 wt;
   float64 energy;
}

/*}}}*/
DPDE_Type;

static void usage (void) /*{{{*/
{
   marx_error ("Usage: marx2dpde MARX-OUTPUT-DIR OUTPUT-FILE-NAME");
   marx_error ("   or: marx2dpde --dump DPDE-FILE-NAME");
   exit (1);
}

/*}}}*/

static FILE *open_input_file (char *dir, char *file, int32 *num_elements) /*{{{*/
{
   Marx_Dump_File_Type *dft;

   if (NULL == (file = marx_dircat (dir, file)))
     return NULL;

   if (NULL == (dft = marx_open_read_dump_file (file)))
     {
	free (file);
	return NULL;
     }

   /* Data file must be 4 byte float. */
   if ('E' == dft->type)
     {
	free (file);
	*num_elements = dft->num_rows;
	return dft->fp;
     }

   marx_error ("File %s does not contain floating point data.", file);
   free (file);
   fclose (dft->fp);
   return NULL;
}

/*}}}*/

static int dpde_create (char *dir, char *output_file) /*{{{*/
{
   FILE *fpout;
   FILE *fp_e, *fp_x, *fp_y, *fp_z, *fp_dx, *fp_dy, *fp_dz;
   int ret;
   DPDE_Header_Type header;
   float32 x, y, z, dx, dy, dz, en;
   float64 output_array[8];
   float64 x_offset = 0.0;
   int32 num_elements;

   fp_e = fp_x = fp_y = fp_z = fp_dx = fp_dy = fp_dz = fpout = NULL;

   ret = 1;

   if ((NULL == (fp_e = open_input_file (dir, "energy.dat", &num_elements)))
       || (NULL == (fp_x = open_input_file (dir, "xpos.dat", &num_elements)))
       || (NULL == (fp_y = open_input_file (dir, "ypos.dat", &num_elements)))
       || (NULL == (fp_z = open_input_file (dir, "zpos.dat", &num_elements)))
       || (NULL == (fp_dx = open_input_file (dir, "xcos.dat", &num_elements)))
       || (NULL == (fp_dy = open_input_file (dir, "ycos.dat", &num_elements)))
       || (NULL == (fp_dz = open_input_file (dir, "zcos.dat", &num_elements))))
     {
	goto return_error;
     }

   if (NULL == (fpout = fopen (output_file, "wb")))
     {
	fprintf (stderr, "Unable to open output file %s\n", output_file);
	goto return_error;
     }

   /* No we are set.  First write out the dpde header */
   header.exposure = 0.0;
   header.n = num_elements;
   header.ftype = 12;

   if ((1 != JDMwrite_float32 (&header.exposure, 1, fpout))
       || (1 != JDMwrite_int32 (&header.n, 1, fpout))
       || (1 != JDMwrite_int32 (&header.ftype, 1, fpout)))
     {
	marx_error ("Write error.");
	goto return_error;
     }

   while (num_elements > 0)
     {
	if ((1 != JDMread_float32 (&x, 1, fp_x))
	    || (1 != JDMread_float32 (&y, 1, fp_y))
	    || (1 != JDMread_float32 (&z, 1, fp_z))
	    || (1 != JDMread_float32 (&dx, 1, fp_dx))
	    || (1 != JDMread_float32 (&dy, 1, fp_dy))
	    || (1 != JDMread_float32 (&dz, 1, fp_dz))
	    || (1 != JDMread_float32 (&en, 1, fp_e)))
	  {
	     marx_error ("Read error.");
	     goto return_error;
	  }

	output_array[0] = -y;
	output_array[1] = z;
	output_array[2] = x_offset + x;
	output_array[3] = -dy;
	output_array[4] = dz;
	output_array[5] = dx;
	output_array[6] = 1.0;
	output_array[7] = en;

	if (8 != JDMwrite_float64 (output_array, 8, fpout))
	  {
	     marx_error ("Write error.");
	     goto return_error;
	  }
	num_elements--;
     }

   /* Success !! */
   ret = 0;
   /* drop */
   return_error:

   if (fp_e != NULL) fclose (fp_e);
   if (fp_x != NULL) fclose (fp_x);
   if (fp_y != NULL) fclose (fp_y);
   if (fp_z != NULL) fclose (fp_z);
   if (fp_dx != NULL) fclose (fp_dx);
   if (fp_dy != NULL) fclose (fp_dy);
   if (fp_dz != NULL) fclose (fp_dz);
   if (fpout != NULL) fclose (fpout);
   return ret;
}

/*}}}*/

static int dpde_dump (char *file) /*{{{*/
{
   FILE *fp;
   DPDE_Type dpde;
   float32 exposure;
   int32 n, ftype, ncounted;

   fp = fopen (file, "rb");

   if (fp == NULL)
     {
	marx_error ("Unable to open %s\n", file);
	return -1;
     }

   if ((1 != JDMread_float32 (&exposure, 1, fp))
       || (1 != JDMread_int32 (&n, 1, fp))
       || (1 != JDMread_int32 (&ftype, 1, fp)))
     {
	marx_error ("Error reading %s as DPDE file.", file);
	fclose (fp);
	return -1;
     }

   fprintf (stdout, "#Exposure: %f\n#number: %d\n#type: %d\n",
	    exposure, n, ftype);

   if (12 != ftype)
     {
	marx_error ("DPDE file mismatch.  Type %d expected.", ftype);
	fclose (fp);
	return -1;
     }

   ncounted = 0;

   while (8 == JDMread_float64 (&dpde.x, 8, fp))
     {
	fprintf (stdout, "%16e\t%16e\t%16e\t%16e\t%16e\t%16e\t%16e\t%16e\n",
		 dpde.x, dpde.y, dpde.z, dpde.dx, dpde.dy, dpde.dz,
		 dpde.wt, dpde.energy);
	ncounted++;
     }

   if ((n != 0)
       && (n != ncounted))
     {
	marx_error ("Warning: %d records read, %d expected.\n", ncounted, n);
     }

   fclose (fp);
   return 0;
}

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   if (argc != 3)
     {
	usage ();
     }

   if (!strcmp (argv[1], "--dump"))
     {
	return dpde_dump (argv[2]);
     }

   return dpde_create (argv[1], argv[2]);
}

/*}}}*/

