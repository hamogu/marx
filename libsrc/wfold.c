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
#include <string.h>
#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#include <jdmath.h>

/* Since the data files store the data as 4 byte floats, they will also be
 * stored as 4 byte floats in memory.  However, delta_p may be so small
 * that any calculation involving it must be performed in double precision.
 * For that reason, everything except the array values will be stored in
 * double precision.
 *
 * Also the probabilities in the array as well as p_min and p_max are
 * cumulative probabilities.
 */
typedef struct
{
   double e_alpha;		       /* energy * sin alpha */
   double p_min;		       /* probability of no scatter */
   double delta_p;		       /* probability increment of array */
   double p_max;		       /* max probability of array */
   double pow_law_norm;		       /* power law tail normalization */
   double pow_law_expon;	       /* exponent for power law tail */
   unsigned int num_theta_values;
   float *theta_values;
}
Fold_Array_Type;

struct _Marx_WFold_Table_Type
{
   unsigned int num_arrays;
   Fold_Array_Type **fold_arrays;      /* NULL terminated, num_arrays of non-NULL */
   double *e_alphas;
};

#define _MARX_FOLD_C_
#include "marx.h"

static void free_a_fold_array (Fold_Array_Type *fat)
{
   if (fat == NULL) return;

   if (fat->theta_values != NULL)
     JDMfree_float_vector (fat->theta_values);

   memset ((char *) fat, 0, sizeof (Fold_Array_Type));
}

static int sanity_check_array (Fold_Array_Type *fat)
{
   double p_min, p_max, delta_p;

   p_min = fat->p_min;
   p_max = fat->p_max;
   delta_p = fat->delta_p;

   if ((delta_p <= 0.0) || (delta_p > 1.0))
     {
	fprintf (stderr, "delta_p is out of range.\n");
	return -1;
     }

   if (fat->e_alpha < 0.0)
     {
	fprintf (stderr, "e_alpha is less than 0.\n");
	return -1;
     }

   if ((p_max > 1.0) || (p_max < 0.0))
     {
	fprintf (stderr, "p_max is out of range.\n");
	return -1;
     }

   if ((p_min > 1.0) || (p_min < 0.0))
     {
	fprintf (stderr, "p_min is out of range.\n");
	return -1;
     }

   return 0;
}

/* Returns 1 if array read in, 0 if at end of file, -1 on error. */

static int read_a_fold_array (Fold_Array_Type *ft, FILE *fp)
{
   int32 n_p_dist;
   unsigned int num_theta_values;
   float *theta_values;

   memset ((char *) ft, 0, sizeof (Fold_Array_Type));

   if (1 != JDMread_int32 (&n_p_dist, 1, fp))
     return 0;

   num_theta_values = (unsigned int) n_p_dist;
   if (num_theta_values == 0)
     return 0;

   if ((1 != JDMread_d_float32 (&ft->e_alpha, 1, fp))
	|| (1 != JDMread_d_float32 (&ft->p_min, 1, fp))
	|| (1 != JDMread_d_float32 (&ft->delta_p, 1, fp))
	|| (1 != JDMread_d_float32 (&ft->p_max, 1, fp))
	|| (1 != JDMread_d_float32 (&ft->pow_law_norm, 1, fp))
	|| (1 != JDMread_d_float32 (&ft->pow_law_expon, 1, fp)))
     {
	fprintf (stderr, "Read error.\n");
	return -1;
     }

   /* To get around possible truncation error due to reading in float values
    * into doubles, adjust delta_p.  This will guarantee with greater precision
    * that num_theta_values and delta_p are consistent.  My interpolation routines
    * that use these data depend upon this.
    */
   if (num_theta_values == 1)
     {
	ft->p_max = ft->p_min;
	ft->delta_p = 0.1;	       /* something greater than 0 */
     }
   else ft->delta_p = (ft->p_max - ft->p_min) / (double) (num_theta_values - 1);

   ft->num_theta_values = num_theta_values;

   if (-1 == sanity_check_array (ft))
     {
	fprintf (stderr, "Table has inconsistent values.\n");
	return -1;
     }

   if (NULL == (theta_values = JDMfloat_vector (num_theta_values)))
     {
	fprintf (stderr, "Out of memory\n");
	return -1;
     }

   if (num_theta_values != JDMread_f_float32 (theta_values, num_theta_values, fp))
     {
	fprintf (stderr, "Read error reading table.\n");
	JDMfree_float_vector (theta_values);
	return -1;
     }

   ft->theta_values = theta_values;

   return 1;
}

static void free_fold_arrays (Fold_Array_Type **fats)
{
   Fold_Array_Type **t;

   t = fats;
   if (t == NULL) return;

   while (*t != NULL)
     {
	free_a_fold_array (*t);
	t++;
     }

   free (fats);
}

static char *x_realloc (char *p, unsigned long size)
{
   if (p == NULL) return (char *) malloc (size);
   return (char *) realloc (p, size);
}

static Fold_Array_Type **read_fold_arrays (FILE *fp)
{
   Fold_Array_Type **arrays, *an_array;
   unsigned int num_arrays_allocated, num_arrays;
   int nread;
   char *err;

   num_arrays_allocated = 0;
   num_arrays = 0;
   an_array = NULL;
   arrays = NULL;
   err = "Not enough memory.";

   /* Note: Care is made to ensure that the arrays list remains NULL terminated
    * since the routines at end of loop assume it.
    */
   while (1)
     {
	if (num_arrays + 1 >= num_arrays_allocated)
	  {
	     Fold_Array_Type **t;

	     num_arrays_allocated += 128;
	     t = (Fold_Array_Type **) x_realloc ((char *)arrays, num_arrays_allocated * sizeof (Fold_Array_Type *));
	     if (t == NULL)
	       break;

	     arrays = t;
	     arrays [num_arrays] = NULL;
	  }

	an_array = (Fold_Array_Type *) malloc (sizeof (Fold_Array_Type));
	if (an_array == NULL)
	  break;

	nread = read_a_fold_array (an_array, fp);
	if (nread == 0)
	  {
	     free (an_array);
	     return arrays;
	  }

	if (nread == -1)
	  {
	     err = "Error reading file";
	     break;
	  }

	arrays [num_arrays++] = an_array;
	arrays [num_arrays] = NULL;
     }

   /* Get here only on error. */
   if (an_array != NULL) free (an_array);
   free_fold_arrays (arrays);
   fprintf (stderr, "read_fold_arrays: %s.\n", err);
   return NULL;
}

static Marx_WFold_Table_Type *read_fold_table (FILE *fp)
{
   unsigned int num, num_arrays;
   Fold_Array_Type **fold_arrays;
   Marx_WFold_Table_Type *table;
   double *e_alphas;

   if (NULL == (fold_arrays = read_fold_arrays (fp)))
     return NULL;

   num_arrays = 0;
   while (fold_arrays [num_arrays] != NULL) num_arrays++;

   if (NULL == (table = (Marx_WFold_Table_Type *) malloc (sizeof (Marx_WFold_Table_Type))))
     {
	free_fold_arrays (fold_arrays);
	fprintf (stderr, "Out of memory.\n");
	return NULL;
     }
   memset ((char *) table, 0, sizeof (Marx_WFold_Table_Type));

   table->fold_arrays = fold_arrays;
   table->num_arrays = num_arrays;

   if (NULL == (table->e_alphas = JDMdouble_vector (num_arrays)))
     {
	free_fold_arrays (fold_arrays);
	free (table);
	return NULL;
     }

   e_alphas = table->e_alphas;
   for (num = 0; num < num_arrays; num++)
     {
	e_alphas[num] = fold_arrays[num]->e_alpha;
     }
   return table;
}

void marx_free_wfold_table (Marx_WFold_Table_Type *table)
{
   if (table == NULL) return;
   if (table->fold_arrays != NULL) free_fold_arrays (table->fold_arrays);
   if (table->e_alphas != NULL) JDMfree_double_vector (table->e_alphas);
   free (table);
}

static double interpolate_theta (Fold_Array_Type *fat, double p)
{
   double delta_i;
   unsigned int i;
   float *t;

   if (p < fat->p_min)
     return 0.0;		       /* no scattering */

   if (p > fat->p_max)
     {
	return pow (fat->pow_law_norm * (1.0 - p), fat->pow_law_expon);
     }

   /* Interpolate in array. */
   t = fat->theta_values;

   delta_i = (p - fat->p_min) / fat->delta_p;
   i = (unsigned int) delta_i;

   /* Since fat->delta_p has been properly adjusted after it was read in,
    * i should be less than or equal num_theta_values and only equal if p = p_max.
    */
   if (i + 1 >= fat->num_theta_values)
     return (double) t[fat->num_theta_values - 1];

   delta_i -= (double) i;	       /* number between 0 and 1 */
   return (1.0 - delta_i) * t[i] + delta_i * t[i + 1];
}

double marx_wfold_table_interp (Marx_WFold_Table_Type *table, double energy, double sin_alpha, double r)
{
   double e_alpha, e_alpha_0, e_alpha_1;
   unsigned int i;
   double theta_0, theta_1, theta;

   if (table == NULL) return 0.0;

   if (table->num_arrays == 1)
     {
	return interpolate_theta (table->fold_arrays[0], r);
     }

   e_alpha = energy * sin_alpha;
   i = JDMbinary_search_d (e_alpha, table->e_alphas, table->num_arrays);

   if (i + 1 == table->num_arrays)
     {
	/* Back off one and extrapolate outside range */
	i--;
     }

   theta_0 = interpolate_theta (table->fold_arrays[i], r);
   theta_1 = interpolate_theta (table->fold_arrays[i + 1], r);
   e_alpha_0 = table->e_alphas[i];
   e_alpha_1 = table->e_alphas[i + 1];
   theta = theta_0 + (theta_1 - theta_0) * (e_alpha - e_alpha_0) / (e_alpha_1 - e_alpha_0);

   /* Arc seconds??  Assume radians */
   return theta /* * (PI / (180.0 * 3600.0))*/ ;
}

static int fold_array_sort_fun (Fold_Array_Type **a, Fold_Array_Type **b)
{
   if ((*a)->e_alpha > (*b)->e_alpha) return 1;
   if ((*a)->e_alpha < (*b)->e_alpha) return -1;
   return 0;
}

static void sort_fold_table (Marx_WFold_Table_Type *table)
{
   Fold_Array_Type **f0, **f1, **fmax;
   unsigned int new_num;
   double *e_alphas, e_alpha;
   /* This is a silly hack to make up for braindead compilers and the lack of
    * uniformity in prototypes for qsort.
    */
   void (*qsort_fun) (char *, unsigned int, unsigned int,
		      int (*)(Fold_Array_Type **,
			      Fold_Array_Type **));

   if ((table == NULL) || (table->num_arrays <= 1))
     return;

   qsort_fun = (void (*)(char *, unsigned int, unsigned int,
			 int (*)(Fold_Array_Type **,
				 Fold_Array_Type **))) qsort;

   qsort_fun ((char *) table->fold_arrays,
	      table->num_arrays, sizeof (Fold_Array_Type *),
	      fold_array_sort_fun);

   /* Now remove duplicates.  It it unlikely there are any. */
   f0 = f1 = table->fold_arrays;
   fmax = f1 + table->num_arrays;

   new_num = 0;
   e_alphas = table->e_alphas;

   while (f1 < fmax)
     {
	*f0 = *f1;
	e_alpha = (*f0)->e_alpha;
	e_alphas[new_num] = e_alpha;
	new_num++;
	f1++;
	f0++;
	while ((f1 < fmax) && (e_alpha == (*f1)->e_alpha))
	  {
	     free_a_fold_array (*f1);
	     f1++;
	  }
     }
   *f0 = NULL;			       /* NULL terminate */
   table->num_arrays = new_num;
}

static Marx_WFold_Table_Type *read_wfold_file (char *file)
{
   FILE *fp;
   Marx_WFold_Table_Type *table;

   if (NULL == (fp = fopen (file, "rb")))
     {
	fprintf (stderr, "marx_open_wfold_file: Unable to open %s\n", file);
	return NULL;
     }

   table = read_fold_table (fp);
   fclose (fp);
   return table;
}

Marx_WFold_Table_Type *marx_read_wfold_file (char *file)
{
   Marx_WFold_Table_Type *table;

   table = read_wfold_file (file);
   if (table != NULL) sort_fold_table (table);
   return table;
}

static void dump_fold_array (Fold_Array_Type *fat)
{
   double delta_p, p;
   float *theta, *theta_max;

   fprintf (stdout, "# e_alpha: %e\n", fat->e_alpha);
   fprintf (stdout, "# p_min: %e\n", fat->p_min);
   fprintf (stdout, "# delta_p: %e\n", fat->delta_p);
   fprintf (stdout, "# p_max: %e\n", fat->p_max);
   fprintf (stdout, "# pow_law_norm: %e\n", fat->pow_law_norm);
   fprintf (stdout, "# pow_law_expon: %e\n", fat->pow_law_expon);
   fprintf (stdout, "# num_theta_values: %u\n", fat->num_theta_values);

   theta = fat->theta_values;
   theta_max = theta + fat->num_theta_values;
   p = fat->p_min;
   delta_p = fat->delta_p;

   fprintf (stdout, "#\n#\tprobability\ttheta\n");
   while (theta < theta_max)
     {
	fprintf (stdout, "\t%12.8e\t%e\n", p, *theta);
	p += delta_p;
	theta++;
     }
   fputs ("\n\n", stdout);
}

int marx_wfold_dump_file (char *file)
{
   Marx_WFold_Table_Type *table;
   Fold_Array_Type **fatp;

   table = read_wfold_file (file);
   if (table == NULL)
     return -1;

   fatp = table->fold_arrays;	       /* should not be NULL */

   while (*fatp != NULL)
     {
	dump_fold_array (*fatp);
	fatp++;
     }

   marx_free_wfold_table (table);
   return 0;
}

