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
#include <string.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <stdarg.h>

#include <jdmath.h>

#include "marx.h"
#include "_marx.h"

int Marx_Error_Code;
int Marx_Verbose = 1;

void marx_error (char *fmt, ...)
{
   char buf[1024];
   va_list ap;

   va_start(ap, fmt);

   if (JDMath_Error)
     {
	vsprintf(buf, fmt, ap);
	JDMmsg_error (buf);
	JDMath_Error = 0;
     }
   else
     {
	vfprintf (stderr, fmt, ap);
	fputc ('\n', stderr);
     }

   va_end(ap);

#ifdef USE_SLANG
   if (SLang_Error == 0)
     SLang_Error = INTRINSIC_ERROR;
#endif

   if (Marx_Error_Code == 0)
     Marx_Error_Code = 1;
}

void marx_message (char *fmt, ...)
{
   va_list ap;

   if (Marx_Verbose == 0)
     return;

   va_start(ap, fmt);
   vfprintf (stdout, fmt, ap);
   va_end(ap);

   fflush (stdout);
}

char *marx_make_version_string (void) /*{{{*/
{
   static char buf[512];

   sprintf (buf, "%s#", MARX_VERSION);

#if defined(__DATE__) && defined(__TIME__)
   sprintf (buf + strlen (buf), " %s %s", __DATE__, __TIME__);
#endif

#ifdef MARX_HOSTNAME
   sprintf (buf + strlen (buf), " %s", MARX_HOSTNAME);
#endif

   return buf;
}
/*}}}*/
