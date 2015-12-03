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

/* This value will be corrected by mirror routines upon initialization. */
double Marx_Mirror_Geometric_Area = 1.0;

static int Mirror = -1;

int marx_mirror_init (Param_File_Type *pf) /*{{{*/
{
   int type;
   char buf[PF_MAX_LINE_LEN];

   if (-1 == pf_get_string (pf, "MirrorType", buf, sizeof (buf)))
     return -1;

   if (!strcmp (buf, "HRMA"))
     type = MARX_MIRROR_HRMA;
   else if (!strcmp (buf, "EA-MIRROR"))
     type = MARX_MIRROR_EA;
   else if (!strcmp (buf, "FLATFIELD"))
     type = MARX_MIRROR_FFIELD;
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_MIRROR_IXO
   else if (!strcmp (buf, "IXO"))
     type = MARX_MIRROR_IXO;
# endif
#endif
   else
     {
	marx_error ("MirrorType %s not supported.", buf);
	return -1;
     }

   switch (type)
     {
      case MARX_MIRROR_EA:
	if (-1 == _marx_ea_mirror_init (pf))
	  type = -1;
	break;

      case MARX_MIRROR_HRMA:
	if (-1 == _marx_hrma_mirror_init (pf))
	  type = -1;
	break;

      case MARX_MIRROR_FFIELD:
	if (-1 == _marx_ff_mirror_init (pf))
	  type = -1;
	break;
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_MIRROR_IXO
      case MARX_MIRROR_IXO:
	if (-1 == _marx_ixo_mirror_init (pf))
	  type = -1;
	break;
# endif
#endif
     }

   if (type == -1)
     {
	marx_error ("Error initializing MirrorType %s.", buf);
     }

   Mirror = type;
   return type;
}

/*}}}*/

int marx_mirror_reflect (Marx_Photon_Type *p, int verbose) /*{{{*/
{
   switch (Mirror)
     {
      case MARX_MIRROR_HRMA:
	if (verbose) marx_message ("Reflecting from HRMA\n");
	return _marx_hrma_mirror_reflect (p);

      case MARX_MIRROR_EA:
	if (verbose) marx_message ("Reflecting from EA-MIRROR\n");
	return _marx_ea_mirror_reflect (p);

      case MARX_MIRROR_FFIELD:
	if (verbose) marx_message ("Flat Fielding Rays\n");
	return _marx_ff_mirror_reflect (p);
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_MIRROR_IXO
      case MARX_MIRROR_IXO:
	if (verbose) marx_message ("Reflecting from IXO\n");
	return _marx_ixo_mirror_reflect (p);
# endif
#endif
     }

   marx_error ("Mirror not initialized.");
   return -1;
}

/*}}}*/

int _marx_parse_shutter_string (char *str, unsigned int *bitmap, unsigned int *num_open) /*{{{*/
{
   unsigned int count;
   char ch;

   *bitmap = 0;			       /* all closed */
   *num_open = 0;

   count = 0;
   while ((ch = *str++) != 0)
     {
	if (ch == '0')
	  {
	     *bitmap |= (1 << count);
	     *num_open += 1;
	  }
	else if (ch != '1')
	  {
	     marx_error ("Ilegal character in shutter bitmap.");
	     return -1;
	  }
	count++;
     }
   if (count != 4)
     {
	marx_error ("Too many shutters specified in shutter bitmap.");
	return -1;
     }
   return 0;
}

/*}}}*/

