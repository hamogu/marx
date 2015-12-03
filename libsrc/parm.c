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
#include <stdlib.h>
#endif

#include <string.h>

#include <pfile.h>

#include "marx.h"
#include "_marx.h"

#ifndef SLMALLOC
#define SLMALLOC malloc
#endif

int _marx_get_parameters (Param_File_Type *p, _Marx_Parameter_Table_Type *table) /*{{{*/
{
   char *name;
   int status;
   char buf[PF_MAX_LINE_LEN];

   while ((name = table->name) != NULL)
     {
	switch (table->type)
	  {
	   case PF_INTEGER_TYPE:
	     status = pf_get_integer (p, name, table->value);
	     break;

	   case PF_REAL_TYPE:
	     status = pf_get_double (p, name, table->value);
	     break;

	   case PF_BOOLEAN_TYPE:
	     status = pf_get_boolean (p, name, table->value);
	     break;

	   case PF_FILE_TYPE:
	     /* drop */
	   case PF_STRING_TYPE:
	     if (table->type == PF_STRING_TYPE)
	       status = pf_get_string (p, name, buf);
	     else status = pf_get_file (p, name, buf);

	     if (status != -1)
	       {
		  char *s = (char *) SLMALLOC (strlen (buf) + 1);
		  if (s == NULL)
		    {
		       marx_error ("Malloc error.");
		       return -1;
		    }
		  strcpy (s, buf);
		  *((char **)table->value) = s;
	       }
	     break;

	   default:
	     marx_error ("Unknown type");
	     return -1;
	  }

	if (status == -1)
	  {
	     marx_error ("Requested parameter %s not found in parameter file.",
			 name);
	     return -1;
	  }

	table++;
     }
   return 0;
}

/*}}}*/

int _marx_get_vector_parm (Param_File_Type *pf, char *parm, JDMVector_Type *v)
{
   char buf[256];
   JDMVector_Type vv;

   if (-1 == pf_get_string (pf, parm, buf, sizeof (buf)))
     {
	marx_error ("Error getting parameter value for %s", parm);
	return -1;
     }

   if (3 != sscanf (buf, "(%lf%*[ ,]%lf%*[ ,]%lf)", &vv.x, &vv.y, &vv.z))
     {
	marx_error ("Parameter %s does not have the proper format", parm);
	return -1;
     }

   *v = vv;
   return 0;
}
