#include <stdio.h>
#include <math.h>
#include <jdmath.h>

static double funct (double x, void *v)
{
   double a = *(double *) v;
   
   return x * exp (-x) - a;
}

static int test_bisect (double a, double b, double c)
{
   double x, y;

   if (-1 == JDM_bisection (funct, a, b, (void *) &c, &x))
     {
	fprintf (stderr, "bisect failed on [%g,%g]: c = %g\n", a, b, c);
	return -1;
     }
   
   y = funct (x, (void *) &c);
   if (fabs (y) <= JDM_Mach_Eps)
     return 0;

   fprintf (stderr, "bisect [%g, %g] failed to find zero: y = %g, x = %g, a = %g\n",
	    a, b, y, x, c);
   return -1;
}
   
   
int main (int argc, char **argv)
{
   double max_y = exp (-1.0);
   unsigned int count = 10000;
   unsigned int i;
   double a;
   unsigned int num_errors = 0;

   JDMATH_INIT;

   for (i = 0; i < count; i++)
     {
	a = JDMrandom () * max_y;
	if (-1 == test_bisect (0, 1, a)) 
	  num_errors++;
	if (-1 == test_bisect (1, 10000, a)) 
	  num_errors++;
     }
   
   /* Check endpoints */
   if (-1 == test_bisect (0, 1, 0.0)) 
     num_errors++;
   if (-1 == test_bisect (1, 10000, 0.0)) 
     num_errors++;
   if (-1 == test_bisect (0, 1, max_y)) 
     num_errors++;
   if (-1 == test_bisect (1, 10000, max_y)) 
     num_errors++;

   if (num_errors == 0)
     fprintf (stderr, "%s passed.\n", argv[0]);
   else
     fprintf (stderr, "%s failed.\n", argv[0]);

   return 0;
}
