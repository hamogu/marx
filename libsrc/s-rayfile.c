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

#include "marx.h"
#include "_marx.h"
#include "source.def"

static FILE *open_marx_rayfile (char *file, unsigned long *history) /*{{{*/
{
   FILE *fp;
   unsigned long magic;

   fp = fopen (file, "rb");
   if (fp == NULL)
     {
	marx_error ("Unable to open rayfile %s.", file);
	return NULL;
     }

   if ((1 != fread (&magic, sizeof (magic), 1, fp))
       || (magic != MARX_RAYFILE_MAGIC)
       || (1 != fread (history, sizeof (unsigned long), 1, fp)))
     {
	fclose (fp);
	marx_error ("Magic number mismatch in rayfile %s.", file);
	return NULL;
     }

   return fp;
}

/*}}}*/

int marx_dump_rayfile (char *file) /*{{{*/
{
   FILE *fp;
   Marx_Photon_Attr_Type at;
   unsigned long history;

   fp = open_marx_rayfile (file, &history);
   if (fp == NULL) return -1;
   /* This routine needs modified to use the history flag!!! */

   fprintf (stdout, "#Energy\t\tTime");
   fprintf (stdout, "\t\tX\t\tY\t\tZ");
   fprintf (stdout, "\t\tPX\t\tPY\t\tPZ");
   fprintf (stdout, "\tY-Pixel\tZ-Pixel");
   fprintf (stdout, "\tOrder\tMirrorShell\tCCD_Num\tHRC_Region");

   putc ('\n', stdout);

   while (1 == fread (&at, sizeof (Marx_Photon_Attr_Type), 1, fp))
     {
	fprintf (stdout, "%f\t%f\t", at.energy, at.arrival_time);
	fprintf (stdout, "%f\t%f\t%f\t", at.x.x, at.x.y, at.x.z);
	fprintf (stdout, "%f\t%f\t%f\t", at.p.x, at.p.y, at.p.z);
	fprintf (stdout, "%f\t%f\t", at.y_pixel, at.z_pixel);

	fprintf (stdout, "%d\t%u\t%d\t%d",
		 at.order, at.mirror_shell, at.ccd_num, at.detector_region);

	putc ('\n', stdout);
     }

   fclose (fp);
   return 0;
}

/*}}}*/

int marx_dump_to_rayfile (char *file, int new_file, /*{{{*/
			  Marx_Photon_Type *p,
			  double total_time)
{
   FILE *fp = NULL;
   unsigned long magic = MARX_RAYFILE_MAGIC;
   int ret = -1;
   unsigned int i, imax;
   Marx_Photon_Attr_Type *attr, *at;

   if (new_file)
     {
	fp = fopen (file, "wb");
	if (fp == NULL) goto return_error;
	if ((1 != fwrite (&magic, sizeof (unsigned long), 1, fp))
	    || (1 != fwrite (&p->history, sizeof (unsigned long), 1, fp)))
	  goto return_error;
     }
   else
     {
	fp = fopen (file, "ab");
	if (fp == NULL) goto return_error;
     }

   attr = p->attributes;
   imax = p->n_photons;

   for (i = 0; i < imax; i++)
     {
	at = attr + i;

	if (at->flags & BAD_PHOTON_MASK)
	  {
	     continue;
	  }

	at->arrival_time += total_time;
	if (1 != fwrite ((char *) at, sizeof (Marx_Photon_Attr_Type), 1, fp))
	  goto return_error;
     }

   ret = 0;
   /* Drop */

   return_error:

   if (fp != NULL) fclose (fp);

   return ret;
}

/*}}}*/

static int rayfile_open_source (Marx_Source_Type *st) /*{{{*/
{
   if (NULL == st->spectrum.s.ray.fp)
     return -1;

   return 0;
}

/*}}}*/

static int rayfile_close_source (Marx_Source_Type *st) /*{{{*/
{
   FILE *fp;

   fp = st->spectrum.s.ray.fp;
   if (fp == NULL)
     return 0;

   fclose (fp);
   st->spectrum.s.ray.fp = NULL;

   return 0;
}

/*}}}*/

static int rayfile_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
				   unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   Marx_Photon_Attr_Type *at;
   FILE *fp;
   double t, last_time;

   fp = st->spectrum.s.ray.fp;
   if (fp == NULL)
     return -1;

   at = pt->attributes;

   num = fread ((char *) at, sizeof (Marx_Photon_Attr_Type), num, fp);

   /* Now adjust arrival times. */
   last_time = st->spectrum.s.ray.last_time;
   t = 0.0;

   for (i = 0; i < num; i++)
     {
	t = at->arrival_time;
	at->arrival_time -= last_time;
	at++;
     }
   st->spectrum.s.ray.last_time = t;
   pt->history = st->spectrum.s.ray.history;

   *num_created = num;
   return 0;
}

/*}}}*/

int marx_select_rayfile_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
				char *name, unsigned int source_id)
{
   char buf [PF_MAX_LINE_LEN];
   int val;
   FILE *fp;

   (void) source_id; (void) name;
   st->open_source = rayfile_open_source;
   st->create_photons = rayfile_create_photons;
   st->close_source = rayfile_close_source;

   if ((-1 != pf_get_boolean (p, "DumpToRayFile", &val))
       && (val != 0))
     {
	marx_error ("You cannot specify a RAYFILE source AND dump to a ray file.");
	return -1;
     }

   if (-1 == pf_get_string (p, "RayFile", buf, sizeof (buf)))
     return -1;

   fp = open_marx_rayfile (buf, &st->spectrum.s.ray.history);
   if (fp == NULL)
     return -1;

   st->spectrum.s.ray.last_time = 0.0;
   st->spectrum.s.ray.fp = fp;
   st->spectrum.type = MARX_RAY_SPECTRUM;

   return 0;
}

/*}}}*/
