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

static int Grating = -1;

int marx_grating_init (Param_File_Type *pf) /*{{{*/
{
   int type;
   char buf[PF_MAX_LINE_LEN];

   if (-1 == pf_get_string (pf, "GratingType", buf, sizeof (buf)))
     return -1;

   if (!strcmp (buf, "LETG"))
     type = MARX_GRATING_LETG;
   else if (!strcmp (buf, "HETG"))
     type = MARX_GRATING_HETG;
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_GRATING_CATGS
   else if (!strcmp (buf, "CATGS"))
     type = MARX_GRATING_CATGS;
# endif
#endif
   else if (!strcmp (buf, "NONE"))
     type = 0;
   else
     {
	marx_error ("GratingType %s not supported.", buf);
	return -1;
     }

   switch (type)
     {
      case MARX_GRATING_HETG:
	if (-1 == _marx_hetg_init (pf))
	  type = -1;
	break;

      case MARX_GRATING_LETG:
	if (-1 == _marx_letg_init (pf))
	  type = -1;
	break;
#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_GRATING_CATGS
      case MARX_GRATING_CATGS:
	if (-1 == _marx_catgs_init (pf))
	  type = -1;
	break;
# endif
#endif
     }

   if (type == -1)
     {
	marx_error ("Error initializing GratingType %s.", buf);
     }

   Grating = type;
   return type;
}

/*}}}*/

int marx_grating_diffract (Marx_Photon_Type *pt, int verbose) /*{{{*/
{
   int status;

   switch (Grating)
     {
      case MARX_GRATING_LETG:
	if (verbose) marx_message ("Diffracting from LETG.\n");
	status = _marx_letg_diffract (pt);
	break;

      case MARX_GRATING_HETG:
	if (verbose) marx_message ("Diffracting from HETG.\n");
	status = _marx_hetg_diffract (pt);
	break;

#if MARX_HAS_IXO_SUPPORT
# ifdef MARX_GRATING_CATGS
      case MARX_GRATING_CATGS:
	if (verbose) marx_message ("Diffracting from CATGS.\n");
	status = _marx_catgs_diffract (pt);
	break;
# endif
#endif

      case 0:
	status = 0;
	break;

      default:
	marx_error ("Grating not initialized.");
	return -1;
     }

   return status;
}

/*}}}*/
