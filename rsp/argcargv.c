/* -*- mode: C; mode: fold; -*- */
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "argcargv.h"

static char *argcargv_strchr (char *name, char ch) /*{{{*/
{
   char ch1;
   while (1)
     {
	ch1 = *name;
	if (ch1 == ch) break;
	if (ch1 == 0) return NULL;
	name++;
     }
   return name;
}

/*}}}*/

static char *argcargv_make_string (char *str) /*{{{*/
{
   char *s;
   
   s = malloc (strlen (str) + 1);
   if (s != NULL) strcpy (s, str);
   return s;
}

/*}}}*/

int argcargv (int *argc_p, char ***argv_p, ArgcArgv_Type *table) /*{{{*/
{
   int argc;
   char **argv, *name, *eqs_value;
   ArgcArgv_Type *t;
   unsigned int name_len;
   char name_buf[256];
   int (*func)(int *, char ***);
   unsigned int or_value;
   
   argc = *argc_p;
   argv = *argv_p;
   
   while (argc > 0)
     {
	name = *argv++;
	argc--;
	
	/* *argc_p = argc; *argv_p = argv; */

	if (NULL != (eqs_value = argcargv_strchr (name, '=')))
	  {
	     name_len = (unsigned int) (eqs_value - name);
	     if (name_len >= sizeof (name_buf))
	       name_len = sizeof (name_buf) - 1;
	     
	     strncpy (name_buf, name, name_len);
	     name_buf[name_len] = 0;
	     name = name_buf;
	     eqs_value++;	       /* skip = */
	  }
	else name_len = strlen (name);

	t = table;
	while (t->name != NULL)
	  {
	     if (!strncmp (t->name, name, name_len))
	       {
		  /* Ok, we have a match upto name_len characters.  Is the
		   * match unique?  Let's assume that the table has been 
		   * ordered in such a way that this is not our concern.
		   */
		  
		  if (eqs_value == NULL)
		    {
		       or_value = 0;
		       eqs_value = *argv;
		    }
		  else or_value = ARGCARGV_EQS_FLAG;
		 
		  switch (t->type | or_value)
		    {
		     case ARGCARGV_STRING:
		       if (argc == 0)
			 {
			    fprintf (stderr, "The %s parameter takes an argument.\n", 
				     t->name);
			    return -1;
			 }
		       argv++; argc--;
		       /* drop */
		     case ARGCARGV_STRING | ARGCARGV_EQS_FLAG:
		       if (NULL == (eqs_value = argcargv_make_string (eqs_value)))
			 return -1;
		       *(char **) (t->addr) = eqs_value;
		       break;
		       
		     case ARGCARGV_INTEGER:
		       if (argc == 0)
			 {
			    fprintf (stderr, "The %s paramter takes an integer argument.\n",
				     t->name);
			    return -1;
			 }
		       argc--; argv++;
		       /* drop */
		     case ARGCARGV_INTEGER | ARGCARGV_EQS_FLAG:
		       if (1 != sscanf (eqs_value, "%d", (int *)t->addr))
			 {
			    fprintf (stderr, "Argument for %s must be an integer.\n", t->name);
			    return -1;
			 }
		       break;
		       
		     case ARGCARGV_FLOAT:
		       if (argc == 0)
			 {
			    fprintf (stderr, "The %s paramter takes a floating point argument.\n",
				     t->name);
			    return -1;
			 }
		       argc--; argv++;
		       /* drop */
		     case ARGCARGV_FLOAT | ARGCARGV_EQS_FLAG:
		       if (1 != sscanf (eqs_value, "%f", (float *)t->addr))
			 {
			    fprintf (stderr, "Argument for %s must be a float.\n", t->name);
			    return -1;
			 }
		       break;

		     case ARGCARGV_DOUBLE:
		       if (argc == 0)
			 {
			    fprintf (stderr, "The %s paramter takes a floating point argument.\n",
				     t->name);
			    return -1;
			 }
		       argc--; argv++;
		       /* drop */
		     case ARGCARGV_DOUBLE | ARGCARGV_EQS_FLAG:
		       if (1 != sscanf (eqs_value, "%lf", (double *)t->addr))
			 {
			    fprintf (stderr, "Argument for %s must be a float.\n", t->name);
			    return -1;
			 }
		       break;
		       
		     case ARGCARGV_BOOLEAN:
		       *(int *) t->addr = 1;
		       break;
		       
		     case ARGCARGV_TOGGLE_BOOLEAN:
		       *(int *) t->addr = ! *(int *) t->addr;
		       break;
		       
		     case ARGCARGV_FUNCTION:
		       func = (int (*)(int *, char ***)) t->addr;
		       if (func == NULL)
			 {
			    fprintf (stderr, "Function for %s is NULL.", t->name);
			    return -1;
			 }
		       argc--; argv++;
		       if (-1 == (*func)(&argc, &argv))
			 {
			    return -1;
			 }
		       break;
		       
		     default:
		       /* This should not happen */
		       break;
		    }
		  
		  break;	       /* break out of while loop */
	       }
	     
	     t++;
	  }
	
	if (t->name == NULL)
	  {
	     /* Not found */
	     argc++;
	     argv--;
	     break;
	  }
     }
   
   *argc_p = argc;
   *argv_p = argv;
   return 0;
}

/*}}}*/


#if 0

char *S1, *S2;
int I1, I2;

ArgcArgv_Type Table[] = /*{{{*/
{
   {"s1", ARGCARGV_STRING, (long) &S1, NULL},
   {"s2", ARGCARGV_STRING, (long) &S2, NULL},
   {"i1", ARGCARGV_INTEGER, (long) &I1, NULL},
   {"i2", ARGCARGV_INTEGER, (long) &I2, NULL},
   {NULL, 0, 0, NULL}
};

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   argc--; argv++;
   if (-1 == argcargv (&argc, &argv, Table))
     {
	fprintf (stderr, "Error.\n");
	return -1;
     }
   if (argc)
     {
	fprintf (stderr, "Unprocessed: %s\n", *argv);
     }
   
   fprintf (stdout, "i1 = %d\n", I1);
   fprintf (stdout, "i2 = %d\n", I2);
   fprintf (stdout, "s1 = %s\n", (S1 == NULL) ? "NULL" : S1);
   fprintf (stdout, "s2 = %s\n", (S2 == NULL) ? "NULL" : S2);
   return 0;
}

/*}}}*/

#endif
