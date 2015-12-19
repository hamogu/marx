/*
    This file is part of the MIT PFILE Parameter Library

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
/* Parse command line arguments */
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <ctype.h>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

#include "pfile.h"
#include "_pfile.h"

static int _pf_set_pf_value (Param_File_Type *p, Param_Type *pf,
			     char *value, unsigned int mode)
{
   mode |= _pf_get_effective_mode (p, pf);
   if (mode & PF_LEARN_MODE)
     {
	if (pf->value != NULL)
	  {
	     SLFREE (pf->value);
	     pf->value = NULL;
	  }

	if (NULL == (pf->value = _pf_create_string (value)))
	  return -1;

	pf->flags |= PF_PARAM_DIRTY;
	p->flags |= PFILE_DIRTY;
     }
   else
     {
	char *save_value = pf->value;
	pf->value = value;

	/* This is very ugly and ought to be cleaned up.  Note also
	 * that indirect values are not properly handled on the command line.
	 */
	_pf_free_current_value (pf);
	p->flags &= ~PF_INDIRECT_VALUE;

	if (*value == ')')
	  p->flags |= PF_INDIRECT_VALUE;

	if (-1 == _pf_get_value (p, pf, 0))
	  {
	     pf->value = save_value;
	     return -1;
	  }

	pf->value = save_value;
     }

   /* Look for the mode parameter and re-initialize the mode based on it. */
   if (!strcmp (pf->name, "mode"))
     {
	unsigned int force_bits = p->mode & (PF_FORCE_LEARN_MODE|PF_NEVER_QUERY_MODE);
	if (-1 == _pf_parse_mode (pf->value, &p->mode))
	  {
	     if (force_bits & PF_NEVER_QUERY_MODE)
	       {
		  pf_error ("Unable to parse the mode parameter");
		  return -1;
	       }
	     if (-1 == _pf_get_value (p, pf, PF_QUERY_MODE))
	       return -1;
	  }
	p->mode |= force_bits;
     }
   return 0;
}

static int set_value (Param_File_Type *p, char *name, char *value, unsigned int mode)
{
   Param_Type *pf;

   pf = _pf_locate_param_by_type (p, name, 0);
   if (pf == NULL)
     {
	pf_error ("Unable to find a parameter named %s.", name);
	return -1;
     }

   if (value == NULL)
     {
	return _pf_get_value (p, pf, PF_QUERY_MODE|mode);
     }
   else return _pf_set_pf_value (p, pf, value, mode);
}

int pf_set_value (Param_File_Type *p, char *name, char *value)
{
   return set_value (p, name, value, 0);
}

int pf_learn_value (Param_File_Type *p, char *name, char *value)
{
   return set_value (p, name, value, PF_LEARN_MODE);
}

Param_File_Type *pf_parse_cmd_line_no_exit (char *file, char *mode, int argc, char **argv)
{
   int i;
   Param_File_Type *p;
   Param_Type *pf;
   char *f;

   if (file == NULL)
     {
	/* derive from argv[0] */
	if (argc == 0) return NULL;
	file = _pf_rstrchr (argv[0], '/');
	if (file != NULL) file++;
	else file = argv[0];
     }

   /* argv[0] not needed anymore (only used if file is NULL) */
   if (argc)
     {
	argc--;
	argv++;
     }

   if ((argc == 1) && (0 == strcmp ("--help", argv[0])))
     {
	pf_usage (file, 0);
	return NULL;
     }

   if (argc)
     {
	/* If argv[0] is of the form @@file, use file as parameter file */
	f = argv[0];
	if ((f[0] == '@')
	    && (f[1] == '@'))
	  {
	     file = f + 2;
	     argv++;
	     argc--;
	  }
     }

   if (NULL == (p = pf_open_parameter_file (file, mode)))
     {
	pf_error ("Unable to open parameter file %s.", file);
	return NULL;
     }

   pf = p->pf;
   while ((pf != NULL) && (pf->type == PF_COMMENT_TYPE))
     pf = pf->next;

   /* The way this works is that as long as the options do not have = signs,
    * they serve to initialize the current pf.  I do not know how the
    * parameter interface handles string values with embedded equal signs.
    * It probably doesn't.  So, for simplicity I will not until I find out
    * otherwise.
    */
   for (i = 0; i < argc; i++)
     {
	char *arg = argv[i];
	char *value;
	char *str;
#if 0
	unsigned int nth = 0;
#endif
	/* The argument may be , separated.  break them out and handle
	 * them one by one.
	 *
	 * **** NOTE****
	 * The datamodel makes extensive use of commas in filename as part
	 * of dm filters.  Just skip comma handling below.
	 */
#if 0
	while (NULL != (value = _pf_extract_string_element (arg, ",", &nth)))
#else
	if (NULL != (value = _pf_create_string (arg)))
#endif
	  {
	     if (*value == 0)
	       {
		  /* ,, case */
		  if (pf != NULL) pf = pf->next;
		  while ((pf != NULL) && (pf->type == PF_COMMENT_TYPE))
		    pf = pf->next;

		  SLFREE (value);
		  continue;
	       }

	     if (pf != NULL)
	       {
		  if (NULL == _pf_strchr (value, '='))
		    {
		       unsigned int len = strlen (value);

		       if (len && (value[len - 1] == '!'))
			 {
			    value[len - 1] = 0;
			    if (-1 == pf_set_value (p, value, NULL))
			      {
				 pf_error ("Error setting %s's value.",
					   pf->name);
				 SLFREE (value);
				 (void) pf_close_parameter_file (p);
				 return NULL;
			      }
			 }
		       else if (-1 == _pf_set_pf_value (p, pf, value, 0))
			 {
			    pf_error ("Error setting %s's value to %s.",
				      pf->name, value);
			    SLFREE (value);
			    (void) pf_close_parameter_file (p);
			    return NULL;
			 }

		       SLFREE (value);
		       pf = pf->next;
		       while ((pf != NULL) && (pf->type == PF_COMMENT_TYPE))
			 pf = pf->next;
		       continue;
		    }
		  pf = NULL;
	       }

	     /* Now we have value=str pairs. */
	     str = _pf_strchr (value, '=');
	     if (str == NULL)
	       {
		  pf_error ("No more parameters: Expecting = sign for '%s'",
			    value);
		  SLFREE (value);
		  (void) pf_close_parameter_file (p);
		  return NULL;
	       }

	     *str++ = 0;	       /* value is malloced so knock out = */
	     if (-1 == pf_set_value (p, value, str))
	       {
		  pf_error ("Error setting value for %s=%s.", value,str);
		  SLFREE (value);
		  (void) pf_close_parameter_file (p);
		  return NULL;
	       }

	     SLFREE (value);
	  }
     }

   return p;
}

Param_File_Type *pf_parse_cmd_line (char *file, char *mode, int argc, char **argv)
{
   Param_File_Type *pf = pf_parse_cmd_line_no_exit (file, mode, argc, argv);
   if (pf == NULL)
     exit (1);
   return pf;
}

void pf_usage (char *file, int quit)
{
   fputs ("\
This program uses an IRAF-style parameter file interface.  It searches for\n\
the parameter file in the PFILES and UPARM directories.  If one is not\n\
found, it will look in the current directory.\n",
	  stderr);

   if (file != NULL)
     {
	fprintf (stderr, "\n\
The name of the parameter file that this program uses is %s.\n",
		 file);
     }

   fputs ("\n\
An alternative parameter file may be specified by prefixing the file name\n\
with \"@@\" and using the resulting expression as the first command line\n\
argument, e.g., program-name @@parameter-file-name.\n",
	  stderr);

   fputs ("\n\
Parameters may be set on the command line via the syntax:\n\
     PARAMETER-NAME=VALUE PARAMETER-NAME=VALUE ...\n\
The program will prompt for a parameter's value if VALUE is not specified,\n\
e.g., \"PARAMETER-NAME=\".  Note that there must be no whitespace surrounding\n\
the '=' sign.\n\
\n\
See your program's user manual for more information.\n",
	  stderr);

   fprintf (stderr, "pfile library version: %s\n\n",
	    PFILE_VERSION);

   if (quit) exit (1);
}

int pf_learn_double (Param_File_Type *pf, char *parm, double x)
{
   char buf[128];

   sprintf (buf, "%.16g", x);
   return pf_learn_value (pf, parm, buf);
}

int pf_learn_int (Param_File_Type *pf, char *parm, int x)
{
   char buf[128];

   sprintf (buf, "%d", x);
   return pf_learn_value (pf, parm, buf);
}

int pf_learn_string (Param_File_Type *pf, char *parm, char *s)
{
   return pf_learn_value (pf, parm, s);
}

