#include <stdio.h>
#include <jdmath.h>

#include "matrix.c"
/* #include "../ludecomp.c" */

static int my_reverse (double **a, double **b, unsigned int n)
{
   copy_matrix (b, a, n);
   if (-1 == JDM_ludecomp_inverse (b, n))
     {
	return -1;
     }
   
   return 1;
}

int main (int argc, char **argv)
{
   unsigned int n = 7;

   if (-1 == test_inverse (n, 1e-6, my_reverse))
     {
     }
   return 0;
}

