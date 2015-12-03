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
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <jdmath.h>
#include "marx.h"

static char *The_Marx_Dir;

static char *make_marx_filename (char *f) /*{{{*/
{
   static char file [1024];
   unsigned int len;

   strcpy (file, The_Marx_Dir);
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

static Marx_Dump_File_Type *open_marx_file (char *file)
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

static float *read_file (char *file, unsigned int *num)
{
   unsigned int this_num;
   float *f;

   Marx_Dump_File_Type *dft;

   dft = open_marx_file (file);
   if (dft == NULL)
     exit (1);

   if (dft->type != 'E')
     {
	marx_error ("File %s does not have correct type (E).", file);
	exit (1);
     }

   this_num = (unsigned int) dft->num_rows;
   if ((this_num == 0)
       || (*num && (this_num != *num)))
     {
	marx_error ("Marx file %s does not have the expected number of elements (%u).",
		    file, *num);
	exit (1);
     }
   *num = this_num;

   f = JDMfloat_vector (this_num);
   if (f == NULL)
     {
	marx_error ("Out of memory while reading %s.", file);
	exit (1);
     }

   if (this_num != JDMread_f_float32 (f, this_num, dft->fp))
     {
	marx_error ("Read error while reading %s.", file);
	exit (1);
     }

   marx_close_read_dump_file (dft);

   return f;
}

static int Grating_Used;
static char *Grating_Type;
static int Detector_Used;
static double Detector_X_Offset;

static void get_marx_par_values (char *file)
{
   Param_File_Type *pf;
   char value [PF_MAX_LINE_LEN];

   pf = pf_open_parameter_file (make_marx_filename (file), "r");
   if (pf == NULL)
     {
	fprintf (stderr, "***Unable to open %s.  Make sure the directory is correct.\n", file);
	exit (1);
     }

   if (-1 == pf_get_string (pf, "DetectorType", value, sizeof (value)))
     exit (1);

   Detector_Used = 1;

   if (0 == strcmp (value, "NONE"))
     {
	fprintf (stdout, "***The parameter file indicates that a detector was not used.\n");
	Detector_Used = 0;
     }

   Grating_Used = 1;
   if (-1 == pf_get_string (pf, "GratingType", value, sizeof (value)))
     exit (1);

   if (0 == strcmp (value, "NONE")) Grating_Used = 0;
   if (0 == strcmp (value, "LETG")) Grating_Used = -1;

   if (Detector_Used
       && (-1 == pf_get_double (pf, "DetOffsetX", &Detector_X_Offset)))
     exit (1);

   (void) pf_close_parameter_file (pf);
}

static signed char *get_marx_char_data_values (unsigned int *num,
					       char *file)
{
   signed char *orders;
   Marx_Dump_File_Type *dft;

   *num = 0;
   if (Grating_Used == 0)
     return NULL;

   dft = marx_open_read_dump_file (make_marx_filename (file));
   if (dft == NULL)
     exit (1);

   *num = (unsigned int) dft->num_rows;
   if (*num == 0)
     {
	fprintf (stderr, "***%s contains no values.\n", file);
	exit (1);
     }

   if (dft->type != 'A')
     {
	fprintf (stderr, "***%s has wrong type.\n", file);
	exit (1);
     }

   orders = (signed char *) malloc (*num);
   if (orders == NULL)
     {
	fprintf (stderr, "***Out of memory.\n");
	exit (1);
     }

   if (*num != fread ((char *) orders, 1, *num, dft->fp))
     {
	fprintf (stderr, "***Read error while reading %s\n", file);
	exit (1);
     }

   marx_close_read_dump_file (dft);

   return orders;
}

#include "argcargv.h"
#include "argcargv.c"

static void usage (void)
{
   fprintf (stderr, "Usage:\n\
bestfocus --dir MARX-DATA-DIR [--par parfile] [--order ORDER] [--grating TYPE]\n\
  where TYPE can must be one of HEG or MEG if the HETG was used.\n"
	    );
   exit (1);
}

static int The_Order;
static char *The_Marx_Par_File;

ArgcArgv_Type Arg_Table [] =
{
   {"--dir", ARGCARGV_STRING, (long) &The_Marx_Dir, NULL},
   {"-d", ARGCARGV_STRING, (long) &The_Marx_Dir, NULL},
   {"--par", ARGCARGV_STRING, (long) &The_Marx_Par_File, NULL},
   {"-p", ARGCARGV_STRING, (long) &The_Marx_Par_File, NULL},
   {"--order", ARGCARGV_INTEGER, (long) &The_Order, NULL},
   {"-o", ARGCARGV_INTEGER, (long) &The_Order, NULL},
   {"--grating", ARGCARGV_STRING, (long) &Grating_Type, NULL},
   {"-g", ARGCARGV_STRING, (long) &Grating_Type, NULL},
   {NULL}
};

int main (int argc, char **argv)
{
   float *xpos, *ypos, *zpos, *xcos, *ycos, *zcos;
   double qy_com, qz_com, q2_com;
   double y_com, z_com, x2_com;
   double xq_com;
   double xave;
   unsigned int num, i, actual_num;
   signed char *orders, *shells;
   int order, meg_orders = 0;
   double a, a_num, a_den, sigma;

   argc--; argv++;

   The_Marx_Par_File = "marx.par";
   The_Order = 0xFFFF;

   if ((-1 == argcargv (&argc, &argv, Arg_Table))
       || (argc != 0)
       || (The_Marx_Dir == NULL))
     usage ();

   order = The_Order;

   orders = NULL;
   shells = NULL;
   num = 0;

   get_marx_par_values (The_Marx_Par_File);

   if (Grating_Used)
     {
	if (Grating_Used != -1)
	  {
	     if (Grating_Type == NULL)
	       {
		  fprintf (stderr, "HETG Grating was used.  You must specifiy HEG or MEG.\n");
		  usage ();
	       }

	     if (0 == strcmp (Grating_Type, "MEG")) meg_orders = 1;
	     else if (0 == strcmp (Grating_Type, "HEG")) meg_orders = 0;
	     else
	       {
		  fprintf (stderr, "TYPE must be MEG or HEG.\n");
		  usage ();
	       }

	     shells = get_marx_char_data_values (&num, "mirror.dat");
	  }
	else shells = NULL;

	orders = get_marx_char_data_values (&num, "order.dat");
	if (order == 0xFFFF)
	  {
	     fprintf (stdout, "***No order specified.  Assuming 0.\n");
	     order = 0;
	  }
     }

   /* We need to open xpos, ypos, zpos, xcos, ycos, zcos */

   xcos = read_file ("xcos.dat", &num);
   ycos = read_file ("ycos.dat", &num);
   zcos = read_file ("zcos.dat", &num);
   xpos = read_file ("xpos.dat", &num);
   ypos = read_file ("ypos.dat", &num);
   zpos = read_file ("zpos.dat", &num);

   /* Compute COM values. */

   y_com = z_com = 0.0;
   qy_com = qz_com = 0.0;
   q2_com = 0.0;
   x2_com = 0.0;
   xq_com = 0.0;
   xave = 0.0;

   actual_num = 0;
   for (i = 0; i < num; i++)
     {
	double qy, qz, px, xi;
	double y_0, z_0;

	if (orders != NULL)
	  {
	     if (order != orders[i])
	       continue;
	     if (shells != NULL)
	       {
		  int s = shells [i];  /* values are 0, 1, 2, 3 */
		  if (meg_orders)
		    {
		       if (s > 1) continue;
		    }
		  else if (s < 2)
		    continue;
	       }
	  }

#if 0
	if (hypot (ypos[i], zpos[i]) > 0.1)
	  continue;
#endif
	if (actual_num != i)
	  {
	     /* Compress input vectors */
	     xcos [actual_num] = xcos [i];
	     ycos [actual_num] = ycos [i];
	     zcos [actual_num] = zcos [i];
	     xpos [actual_num] = xpos [i];
	     ypos [actual_num] = ypos [i];
	     zpos [actual_num] = zpos [i];
	  }

	px = xcos[actual_num];
	xi = xpos[actual_num];

	if (px == 0.0)
	  {
	     fprintf (stderr, "***Ray %u has 0 x-momentum component.\n", i);
	     exit (1);
	  }

	xave += xi;

	qy = ycos[actual_num] / px;
	qz = zcos[actual_num] / px;
	y_0 = ypos[actual_num] - xi * qy;
	z_0 = zpos[actual_num] - xi * qz;

	y_com += y_0;
	z_com += z_0;
	x2_com += (y_0 * y_0 + z_0 * z_0);

	qy_com += qy;
	qz_com += qz;
	q2_com += (1.0 + qy * qy + qz * qz);
	xq_com += (y_0 * qy + z_0 * qz);

	actual_num++;
     }

   if (actual_num == 0)
     {
	fprintf (stderr, "***There are no rays with order = %d.\n", order);
	exit (1);
     }
   fprintf (stdout, "\n%u rays read. %u of them used (%.2f%%).\n", num, actual_num, actual_num/(1.0*num));
   if (Grating_Used)
     fprintf (stdout, "Of those, %u have order = %d in specified grating.\n", actual_num, order);

   y_com = y_com / actual_num;
   z_com = z_com / actual_num;
   x2_com = x2_com / actual_num;
   qy_com = qy_com / actual_num;
   qz_com = qz_com / actual_num;
   q2_com = q2_com / actual_num;
   xq_com = xq_com / actual_num;
   xave = xave / actual_num;

   a_num = xq_com - (y_com * qy_com + z_com * qz_com);
   a_den = (1.0 + qy_com * qy_com + qz_com * qz_com - q2_com);
   a = a_num / a_den;

   fprintf (stdout, "\nBest focus position [mm]: (%f, %f, %f)\n",
	    a, (y_com + a * qy_com), (z_com + a * qz_com));

   sigma = (x2_com + 2.0 * a * xq_com + a * a * q2_com
	    - (y_com * y_com + z_com * z_com
	       + 2.0 * a * (y_com * qy_com + z_com * qz_com)
	       + a * a * (1.0 + qy_com * qy_com + qz_com * qz_com)));

   if (fabs (sigma) < 1.0e-10) sigma = 0.0;

   sigma = sqrt (sigma);

   fprintf (stdout, "At that position, the mean radius of the spot will be %f mm\n",
	    sigma);
   fprintf (stdout, "  (about %f ACIS Pixels)\n", sigma / 0.024);

   if (Detector_Used)
     {
	fprintf (stdout, "\nCurrent DetOffsetX value: %f\n", Detector_X_Offset);
	fprintf (stdout, "Suggested DetOffsetX value: %f (difference of %f)\n",
		 Detector_X_Offset + (a - xave), a - xave);
     }

#if 0
   fprintf (stdout, "\t(Ave X coordinate of input rays: %f mm)\n", xave);
#endif

   JDMfree_float_vector (xpos);
   JDMfree_float_vector (ypos);
   JDMfree_float_vector (zpos);
   JDMfree_float_vector (xcos);
   JDMfree_float_vector (ycos);
   JDMfree_float_vector (zcos);

   if (orders != NULL) free (orders);
   if (shells != NULL) free (shells);

   return 0;
}

