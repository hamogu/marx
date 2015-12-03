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
#include "config.h"

#include "marx-feat.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"
#include "source.def"

#if MARX_HAS_DYNAMIC_LINKING
#include <dlfcn.h>

#ifndef isspace
#define isspace(x) (((x)==' ')||((x)=='\t')||((x)=='\n'))
#endif

static int (*User_Open_Source)(char **, int, double, double, double, double);
static void (*User_Close_Source)(void);
static int (*User_Generate_Ray)(double *, double *, double *,
				double *, double *);
static int (*User_Start_Iteration)(void);

static char *User_Handle;
static char User_Cmd_String [PF_MAX_LINE_LEN];
static char **User_Argv;
static int User_Argc;

static char *skip_white (char *p)
{
   while (isspace (*p))
     p++;

   return p;
}

static char *skip_non_white (char *p)
{
   while (*p && (0 == isspace (*p)))
     p++;

   return p;
}

static int cmdline_to_argv_argc (char *cmd, char ***argvp, int *argcp)
{
   char *p;
   int argc;
   char **argv;

   p = cmd = skip_white (cmd);

   argc = 0;
   while (*p != 0)
     {
	p = skip_non_white (p);
	p = skip_white (p);
	argc++;
     }

   if (NULL == (argv = (char **) marx_malloc (sizeof (char *) * (argc + 1))))
     return -1;

   argc = 0;
   p = cmd;
   while (*p != 0)
     {
	argv [argc] = p;
	argc++;

	p = skip_non_white (p);
	if (*p != 0) *p++ = 0;

	p = skip_white (p);
     }
   argv [argc] = NULL;

   *argcp = argc;
   *argvp = argv;

   return 0;
}

static int open_source (Marx_Source_Type *st) /*{{{*/
{

   if (User_Open_Source == NULL)
     return -1;

   if (-1 == cmdline_to_argv_argc (User_Cmd_String, &User_Argv, &User_Argc))
     return -1;

   return (*User_Open_Source)(User_Argv, User_Argc, Marx_Mirror_Geometric_Area,
			      st->p.x, st->p.y, st->p.z);
}

/*}}}*/

static int close_source (Marx_Source_Type *st) /*{{{*/
{
   (void) st;

   if (User_Close_Source != NULL)
     (*User_Close_Source)();

   if (User_Handle != NULL)
     dlclose (User_Handle);

   User_Handle = NULL;
   User_Open_Source = NULL;
   User_Close_Source = NULL;
   User_Generate_Ray = NULL;
   User_Start_Iteration = NULL;

   if (User_Argv != NULL)
     {
	marx_free ((char *) User_Argv);
	User_Argv = NULL;
     }
   User_Argc = 0;

   return 0;
}

/*}}}*/

static int create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
			   unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   Marx_Photon_Attr_Type *at;
   int (*efun) (Marx_Spectrum_Type *, double *);
   double t, energy, energy_drawn, last_time;

   if (User_Start_Iteration != NULL)
     {
	if (-1 == (*User_Start_Iteration) ())
	  {
	     *num_created = 0;
	     return 0;
	  }
     }

   at = pt->attributes;
   efun = st->spectrum.energy_function;

   t = -1.0;
   last_time = 0.0;

   for (i = 0; i < num; i++)
     {
	double p1, p2, p3;

	if (-1 == (*efun) (&st->spectrum, &energy_drawn))
	  return -1;
	
	energy = energy_drawn;

	if (-1 == (*User_Generate_Ray) (&t, &energy, &p1, &p2, &p3))
	  break;

	at->flags = 0;

	at->p.x = p1;
	at->p.y = p2;
	at->p.z = p3;

	if (t >= 0.0)
	  {
	     last_time += t;
	     at->arrival_time = last_time;
	  }

	/* This will not be called, unless the User_Generate_Ray resets the
           energy, but for backwards compatibility, we need it:
           Up to marx 5.4 the energy before user_generate_ray
           was undefined and the user source was expected to set it to
           -1, if it wanted marx to draw as usual.
	*/
	if (energy < 0.0)
	  {
	    energy = energy_drawn;
	  }

	at->energy = energy;
	at++;
     }

   *num_created = i;

   if (t >= 0.0)
     {
	pt->history = (MARX_ENERGY_OK
		       | MARX_TIME_OK
		       | MARX_X_VECTOR_OK
		       | MARX_P_VECTOR_OK);
     }

   return 0;
}

/*}}}*/

typedef void (*Fun_Ptr)(void);

static Fun_Ptr get_dlsym (char *name, int is_required)
{
   Fun_Ptr f;

   if (User_Handle == NULL)
     return NULL;

   f = (Fun_Ptr) dlsym (User_Handle, name);
   if ((f == NULL) && is_required)
     {
	char *err = (char *)dlerror ();

	if (err == NULL) err = "Unknown";
	marx_error ("Unable to get symbol '%s'\nReason: %s\n",
		    name, err);
	return NULL;
     }
   return f;
}

int marx_select_user_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			     char *name, unsigned int source_id)
{
   char file [PF_MAX_LINE_LEN];
   char *handle;

   (void) source_id;

   if (-1 == pf_get_file (p, "UserSourceFile", file, sizeof (file)))
     {
	marx_error ("Unable to find parameter 'UserSourceFile'");
	return -1;
     }

   User_Cmd_String [0] = 0;
   if (-1 == pf_get_string (p, "UserSourceArgs", User_Cmd_String, sizeof (User_Cmd_String)))
     return -1;

   marx_message ("Dynamically linking to file %s\n", file);

   handle = (char *) dlopen (file, RTLD_LAZY);
   if (handle == NULL)
     {
	char *err;

	err = (char *) dlerror ();
	if (err == NULL) err = "UNKNOWN";

	marx_error ("Error linking to %s\nReason: %s", file, err);
	return -1;
     }

   User_Handle = handle;
   if (NULL == (User_Open_Source = (int (*)(char **, int, double, double, double, double)) get_dlsym ("user_open_source", 1)))
     return -1;

   if (NULL == (User_Close_Source = (void (*)(void)) get_dlsym ("user_close_source", 1)))
     return -1;

   if (NULL == (User_Generate_Ray = (int (*)(double *, double *,
					     double *, double *, double *))
		get_dlsym ("user_create_ray", 1)))
     return -1;

   User_Start_Iteration = (int (*)(void)) get_dlsym ("user_start_iteration", 0);

   st->open_source = open_source;
   st->close_source = close_source;
   st->create_photons = create_photons;

   return _marx_get_simple_specrum_parms (p, st, name);
}

/*}}}*/

#else
int marx_select_user_source (Marx_Source_Type *st, Param_File_Type *p,
			     char *name, unsigned int source_id)
{
   (void) st;
   (void) p;
   (void) name;
   (void) source_id;

   marx_error ("This version of MARX does not support dynmamic linking");
   return -1;
}
#endif				       /* MARX_HAS_DYNAMIC_LINKING */
