#include <stdio.h>
#include <jdmath.h>
#include <slang.h>

#include "user.h"

static SLang_Name_Type *Open_Source;
static SLang_Name_Type *Create_Ray;

static SLang_Intrin_Fun_Type Intrinsics [] =
{
   SLANG_END_INTRIN_FUN_TABLE
};

static int init_slang (void)
{
   if ((-1 == SLang_init_all ())
       || (-1 == SLang_init_array_extra ())
       || (-1 == SLang_init_import ()) /* dynamic linking */
       || (-1 == SLadd_intrin_fun_table (Intrinsics, NULL)))
     {
	SLang_verror (0, "Unable to initialize S-Lang.\n");
	return -1;
     }
   return 0;
}

static int push_c_string_array (char **argv, int argc)
{
   SLang_Array_Type *at;
   char **strs;
   int i;

   if (NULL == (at = SLang_create_array (SLANG_STRING_TYPE, 1, NULL, &argc, 1)))
     return -1;
   
   strs = (char **) at->data;
   for (i = 0; i < argc; i++)
     {
	if (NULL == (strs[i] = SLang_create_slstring (argv[i])))
	  {
	     SLang_free_array (at);
	     return -1;
	  }
     }
   
   return SLang_push_array (at, 1);
}

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   char *file;
   int status;

   if (-1 == init_slang ())
     return -1;
   
   file = argv[0];
   if ((argc == 0) || (NULL == (file = argv[0])))
     {
	fprintf (stderr, "No filename specified for the slang source\n");
	return -1;
     }

   if (0 != SLang_load_file (file))
     {
	fprintf (stderr, "Encountered a problem loading %s\n", file);
	return -1;
     }

   if (NULL == (Open_Source = SLang_get_function ("user_open_source")))
     {
	fprintf (stderr, "%s failed to define user_open_source\n", file);
	return -1;
     }

   if (NULL == (Create_Ray = SLang_get_function ("user_create_ray")))
     {
	fprintf (stderr, "%s failed to define user_create_ray\n", file);
	return -1;
     }

   if ((-1 == SLang_start_arg_list ())
       || (-1 == push_c_string_array (argv, argc))
       || (-1 == SLang_push_double (area))
       || (-1 == SLang_push_double (cosx))
       || (-1 == SLang_push_double (cosy))
       || (-1 == SLang_push_double (cosz))
       || (-1 == SLang_end_arg_list ())
       || (-1 == SLexecute_function (Open_Source))
       || (-1 == SLang_pop_integer (&status)))
     {
	SLang_verror (0, "Error occured processing user_open_source in %s", file);
	return -1;
     }

   return status;
}

void user_close_source (void)
{
   (void) SLang_run_hooks ("user_close_source", 0);
}


typedef struct
{
   SLang_Array_Type *at;
   double *data;
   unsigned int next_i;
   unsigned int di;
   unsigned int num_elements;
}
Array_Type;

static Array_Type CosX_Array;
static Array_Type CosY_Array;
static Array_Type CosZ_Array;
static Array_Type Energy_Array;
static Array_Type dT_Array;
static unsigned int Num_Rays;

static int pop_array (Array_Type *a)
{
   SLang_Array_Type *at;

   if (a->at != NULL)
     SLang_free_array (a->at);

   if (-1 == SLang_pop_array_of_type (&at, SLANG_DOUBLE_TYPE))
     {
	SLang_verror (0, "Expecting an array");
	return -1;
     }
   a->at = at;
   a->data = (double *)at->data;
   a->num_elements = at->num_elements;
   a->next_i = 0;
   a->di = 1;
   if (Num_Rays < a->num_elements)
     Num_Rays = a->num_elements;

   return 0;
}

static double next_element (Array_Type *a)
{
   unsigned int i = a->next_i;
   double x = a->data[i];
   a->next_i = i + a->di;
   return x;
}


int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   if (Num_Rays == 0)
     {
	if (-1 == SLexecute_function (Create_Ray))
	  {
	     SLang_verror (0, "Encountered an error processing %s\n", "user_create_ray");
	     return -1;
	  }
	
	if (SLang_peek_at_stack () == SLANG_NULL_TYPE)
	  return -1;		       /* done */

	if ((-1 == pop_array (&CosZ_Array))
	    || (-1 == pop_array (&CosY_Array))
	    || (-1 == pop_array (&CosX_Array))
	    || (-1 == pop_array (&Energy_Array))
	    || (-1 == pop_array (&dT_Array)))
	  {
	     SLang_verror (0, "Encountered an error processing %s\n", "user_create_ray");
	     return -1;
	  }
	
	if (Num_Rays == 0)
	  return -1;
	
	if (CosX_Array.num_elements < Num_Rays)
	  CosX_Array.di = 0;
	if (CosY_Array.num_elements < Num_Rays)
	  CosY_Array.di = 0;
	if (CosZ_Array.num_elements < Num_Rays)
	  CosZ_Array.di = 0;
	if (dT_Array.num_elements < Num_Rays)
	  dT_Array.di = 0;
	if (Energy_Array.num_elements < Num_Rays)
	  Energy_Array.di = 0;
     }

   *cosx = next_element (&CosX_Array);
   *cosy = next_element (&CosY_Array);
   *cosz = next_element (&CosZ_Array);
   *delta_t = next_element (&dT_Array);
   *energy = next_element (&Energy_Array);
   
   Num_Rays--;
   return 0;
}

int main (int argc, char **argv)
{
   (void) argc;
   (void) argv;
   return 1;
}
