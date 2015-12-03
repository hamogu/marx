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
/* Dump a grating efficiency binary file. */
#include "marx.h"

int main (int argc, char **argv)
{
#if 1
   (void) argc;
   (void) argv;
   fprintf (stderr, "%s not implemented.\n", argv[0]);
   return 1;
#else
   char *file;
   Marx_Grating_Eff_Type *g;
   unsigned int num_energies, num;
   int order, max_order;

   if (argc != 2)
     {
	fprintf (stderr, "Usage: geffdump <FILE>\n");
	return 1;
     }
   file = argv[1];

   g = marx_read_grating_efficiencies (file);
   if (g == NULL) return 1;

   max_order = g->max_order;
   num_energies = g->num_energies;

   fprintf (stdout, "Max_Order: %d\n", max_order);
   fprintf (stdout, "Num_Energies: %d\n", num_energies);

   for (num = 0; num < num_energies; num++)
     {
	float *eff;

	fprintf (stdout, "%e", g->energies [num]);

	eff = g->efficiencies[num];

	for (order = -max_order; order <= max_order; order++)
	  fprintf (stdout, "\t%e", eff[order + max_order]);

	fputs ("\n", stdout);
     }

   /* marx_free_grating_efficiencies (g); */
   return 0;
#endif
}

