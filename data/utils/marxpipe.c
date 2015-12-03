/*
    This file is part of MARX

    Copyright (C) 2002 Massachusetts Institute of Technology

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
#include <math.h>
#include <time.h>

#include "marx.h"

static int read_pipe (FILE *fp, unsigned int *num_detectedp, unsigned int *num_inputp, double *total_timep)
{
   *num_detectedp = *num_inputp = 0;
   *total_timep = 0;
   while (1)
     {
	unsigned int num_sorted;
	unsigned int num_input;
	double tstart, duration;

	if (1 != fread (&num_sorted, sizeof (unsigned int), 1, fp))
	  return 0;

	if ((1 != fread (&num_input, sizeof (unsigned int), 1, fp))
	    || (1 != fread (&tstart, sizeof (double), 1, fp))
	    || (1 != fread (&duration, sizeof (double), 1, fp)))
	  {
	     fprintf (stderr, "Error reading from marx pipe\n");
	     return -1;
	  }

	*num_detectedp += num_sorted;
	*num_inputp += num_input;

	while (num_sorted != 0)
	  {
	     Marx_Photon_Attr_Type a;

	     if (1 != fread (&a, sizeof (Marx_Photon_Attr_Type), 1, fp))
	       {
		  fprintf (stderr, "Error reading from marx pipe\n");
		  return -1;
	       }
	     num_sorted--;
	  }
	*total_timep = tstart + duration;
     }
}

int main (int argc, char **argv)
{
   unsigned int num_detected, num_input;
   double total_time;

   if (-1 == read_pipe (stdin, &num_detected, &num_input, &total_time))
     return 1;

   fprintf (stdout, "%u\t%u\t%.16g\n", num_input, num_detected, total_time);
   return 0;
}

