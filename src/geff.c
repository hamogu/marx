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
#include <stdio.h>
#include <limits.h>

#include <marx.h>

static void usage (char *pgm) /*{{{*/
{

   fprintf (stderr, "Usage: %s PFILE min-energy max-energy num-energies min-order max-order\n",
	    pgm);
   fprintf (stderr, "where: min-energy >= 0.03, max-energy < 10.0\n");
}

/*}}}*/

static Marx_Grating_Info_Type Geom;
static char *Optical_Constants;

static Param_Table_Type Grating_Parm_Table [] = /*{{{*/
{
   {"Gold",			PF_REAL_TYPE,	&Geom.gold},
   {"Nickel",			PF_REAL_TYPE,	&Geom.nickel},
   {"Chromium",			PF_REAL_TYPE,	&Geom.chromium},
   {"Polyimide",		PF_REAL_TYPE,	&Geom.polyimide},
   {"BarHeight",		PF_REAL_TYPE,	&Geom.bar_height},
   {"BarWidth",			PF_REAL_TYPE,	&Geom.bar_width},
   {"Period",			PF_REAL_TYPE,	&Geom.period},
   {"GratingOptConsts",		PF_FILE_TYPE,	&Optical_Constants},
   {NULL, 0, NULL}
};

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   Param_File_Type *pf;
   int order;
   char *pfile;
   double emin, emax, de;
   int min_order, max_order;
   int num_energies;

   if ((argc != 7)
       || (1 != sscanf (argv[2], "%lf", &emin))
       || (1 != sscanf (argv[3], "%lf", &emax))
       || (1 != sscanf (argv[4], "%d", &num_energies))
       || (1 != sscanf (argv[5], "%d", &min_order))
       || (1 != sscanf (argv[6], "%d", &max_order))
       || (emin > emax)
       || (min_order > max_order)
       || (emin < 0.03)
       || (num_energies <= 0)
       || (emax >= 10.0))
     {
	usage (argv[0]);
	return 1;
     }

   pfile = argv[1];

   if (NULL == (pf = pf_open_parameter_file (pfile, "r")))
     {
	fprintf (stderr, "Unable to open param file %s.\n", pfile);
	return 1;
     }

   if (-1 == pf_get_parameters (pf, Grating_Parm_Table))
     {
	fprintf (stderr, "Error obtain parameters from %s\n", pfile);
	pf_close_parameter_file (pf);
	return 1;
     }

   pf_close_parameter_file (pf);

   if (-1 == marx_create_grating_opt_const_tables (Optical_Constants))
     {
	fprintf (stderr, "Error getting optical constants from %s\n", Optical_Constants);
	return 1;
     }

   order = 0;
   de = (emax - emin) / num_energies;

   while (emin < emax)
     {
	fprintf (stdout, "%e", emin);
	for (order = min_order; order <= max_order; order++)
	  {
	     double eff;

	     eff = marx_compute_grating_efficiency (emin, order, &Geom);
	     if (emin < 0.0)
	       {
		  fprintf (stderr, "Error occured while computing efficiency\n");
		  return 1;
	       }

	     fprintf (stdout, "\t%e", eff);
	  }
	fputc ('\n', stdout);
	emin += de;
     }

   return 0;
}

/*}}}*/
