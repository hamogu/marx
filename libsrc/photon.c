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

#include "marx.h"
#include "_marx.h"

void marx_prune_photons (Marx_Photon_Type *pt) /*{{{*/
{
   double *energies;
   unsigned int *sorted_index;
   unsigned int i, num_sorted, new_num_sorted;
   Marx_Photon_Attr_Type *at, *attr;

   num_sorted = pt->num_sorted;
   attr = pt->attributes;
   energies = pt->sorted_energies;
   sorted_index = pt->sorted_index;
   new_num_sorted = 0;
   for (i = 0; i < num_sorted; i++)
     {
	int indx = sorted_index[i];
	at = attr + indx;
	if (at->flags & BAD_PHOTON_MASK) continue;

	energies[new_num_sorted] = at->energy;
	sorted_index[new_num_sorted] = indx;
	new_num_sorted++;
     }
   pt->num_sorted = new_num_sorted;
}

/*}}}*/

Marx_Photon_Type *marx_alloc_photon_type (unsigned int num) /*{{{*/
{
   Marx_Photon_Type *pt;

   if (NULL == (pt = (Marx_Photon_Type *) marx_malloc (sizeof (Marx_Photon_Type))))
     return NULL;

   memset ((char *) pt, 0, sizeof (Marx_Photon_Type));

   if (NULL == (pt->attributes = (Marx_Photon_Attr_Type *) marx_calloc (num, sizeof (Marx_Photon_Attr_Type))))
     goto error_return;

   if (NULL == (pt->sorted_energies = JDMdouble_vector (num)))
     goto error_return;

   pt->max_n_photons = num;
   return pt;

   error_return:
   marx_dealloc_photon_type (pt);
   return NULL;
}

/*}}}*/

int marx_dealloc_photon_type (Marx_Photon_Type *pt) /*{{{*/
{
   if (pt == NULL) return -1;

   marx_free ((char *) pt->attributes);

   if (pt->sorted_index != NULL)
     JDMfree_integer_vector ((int *) pt->sorted_index);

   if (pt->sorted_energies != NULL)
     JDMfree_double_vector (pt->sorted_energies);

   marx_free ((char *)pt);
   return 0;
}

/*}}}*/

