#include <math.h>

#define MIN(x,y) (((x) < (y)) ? (x) : (y))

static double random_matrix (unsigned int n, unsigned int i, unsigned int j)
{
   (void) n; (void) i; (void) j;
   return JDMrandom ();
}

static double hilbert (unsigned int n, unsigned int i, unsigned int j)
{
   (void) n;
   return 1.0 / (i+j+1);
}

static double ding_dong (unsigned int n, unsigned int i, unsigned int j)
{
   return 0.5 / ((double)n - (i + j) - 0.5);
}

static double moler (unsigned int n, unsigned int i, unsigned int j)
{
   (void) n;
   if (i == j) return i + 1;
   return MIN(i,j) - 1.0;
}

static double frank (unsigned int n, unsigned int i, unsigned int j)
{
   (void) n;
   return MIN(i,j) + 1.0;
}

static double bordered (unsigned int n, unsigned int i, unsigned int j)
{
   if (i == j) return 1.0;
   if (i + 1 == n) return pow (2.0, -1.0 * i);
   if (j + 1 == n) return pow (2.0, -1.0 * j);
   return 0.0;
}

static double diagonal (unsigned int n, unsigned int i, unsigned int j)
{
   (void) n;
   if (i == j) return i + 1;
   return 0;
}

static double wilkinson_wp (unsigned int n, unsigned int i, unsigned int j)
{
   if (i == j)
     return (n/2) + 1.0 - (double)MIN(i+1,n-i);
   if ((i + 1 == j) || (j + 1 == i))
     return 1.0;
   return 0.0;
}

static double wilkinson_wn (unsigned int n, unsigned int i, unsigned int j)
{
   if (i == j)
     return (double)(n/2) - (double)i;
   if ((i + 1 == j) || (j + 1 == i))
     return 1.0;
   return 0.0;
}

static double ones (unsigned int n, unsigned int i, unsigned int j)
{
   (void) n; (void) i; (void)j;
   return 1.0;
}



static void init_matrix (double **a, unsigned int n,
			 double (*f)(unsigned int, unsigned int, unsigned int))
{
   unsigned int i, j;
   
   for (i = 0; i < n; i++)
     for (j = 0; j < n; j++)
       a[i][j] = (*f)(n, i, j);
}

static void copy_matrix (double **a, double **b, unsigned int n)
{
   unsigned int i;

   for (i = 0; i < n; i++)
     memcpy ((char *) a[i], (char *) b[i], n * sizeof (double));
}

typedef struct 
{
   char *name;
   double (*f)(unsigned int, unsigned int, unsigned int);
   int is_singular;
}
Special_Matrix_Type;

static Special_Matrix_Type Special_Matrices [] = 
{
   {"random", random_matrix},
   {"hilbert", hilbert, 0},
   {"ding_dong", ding_dong, 0},
   {"moler", moler, 0},
   {"frank", frank, 0},
   {"bordered", bordered, 0},
   {"diagonal", diagonal, 0},
   {"wilkinson_wp", wilkinson_wp, 0},
   {"wilkinson_wn", wilkinson_wn, 0},
   {"ones", ones, 1},
   {NULL}
};

static int test_this_inverse (double **a, double **b, unsigned int n,
			      int (*fun) (double **, double **, unsigned int),
			      double *maxp)
{
   double max;
   unsigned int i, j;
   int status;

   status = (*fun) (a, b, n);
   if (status == -1)
     return -1;
   if (status == 0)
     return 0;
   
   max = 0.0;

   for (i = 0; i < n; i++)
     {
	unsigned int k;
	for (k = 0; k < n; k++)
	  {
	     double sum = 0;
	     for (j = 0; j < n; j++)
	       sum += a[i][j] * b[j][k];
	     
	     if (i == k)
	       sum -= 1.0;
		  
	     if (fabs (sum) > max) max = fabs(sum);
	  }
     }
   
   *maxp = max;
   return 1;
}


static int test_inverse (unsigned int n, double tol,
			 int (*fun) (double **, double **, unsigned int))
{
   Special_Matrix_Type *s;
   double **a, **b;
   
   if (NULL == (a = JDMdouble_matrix (n, n)))
     return -1;
   
   if (NULL == (b = JDMdouble_matrix (n, n)))
     {
	JDMfree_double_matrix (a, n);
	return -1;
     }

   s = Special_Matrices;
   while (s->f != NULL)
     {
	int status;
	double val;

	init_matrix (a, n, s->f);
	status = test_this_inverse (a, b, n, fun, &val);
	if (status == -1)
	  return -1;
	
	if (status == s->is_singular)
	  {
	     fprintf (stderr, "%s: failed singularity test\n", s->name);
	     s++;
	     continue;
	  }
	
	if (s->is_singular)
	  {
	     s++;
	     continue;
	  }
	
	if (val > tol)
	  fprintf (stderr, "%s failed: tolerance exceeded (%e)\n", s->name, val);
	
	s++;
     }
   
   return 0;
}
