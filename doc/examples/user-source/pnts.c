#include <stdio.h>
#include <stdlib.h>
#include <jdmath.h>
#include "user.h"

/* This user source implements many point sources via a file that 
 * specifies the source positions and energies.  The current implementation
 * assumes the format:
 *  RA  Dec weight energy
 * Here RA, Dec specifiy the source position, weight specifies the strength
 * of the source in relation to the others.
 */
typedef struct
{
   double cosx, cosy, cosz;
   double weight;
   double energy;
}
Point_Source_Type;

static unsigned int Num_Points;
static Point_Source_Type *Point_Sources;
static unsigned int Max_Num_Points;

static char *do_realloc (char *p, unsigned int len)
{
   if (p == NULL)
     p = malloc (len);
   else
     p = realloc (p, len);
   
   if (p == NULL)
     fprintf (stderr, "Not enough memory\n");
   
   return p;
}

static void free_sources (void)
{
   if (Point_Sources == NULL)
     return;
   
   free ((char *) Point_Sources);
   Point_Sources = NULL;
}

static int add_source (double ra, double dec, double weight, double energy)
{
   Point_Source_Type *p;
   double cosx, cosy, cosz;

   /* Convert to God's units from arc-min */
   ra = ra * (PI/(180.0 * 60.0));
   dec = dec * (PI/(180.0 * 60.0));
   
   if (Max_Num_Points == Num_Points)
     {
	Max_Num_Points += 32;
	p = (Point_Source_Type *)do_realloc ((char *)Point_Sources, Max_Num_Points * sizeof (Point_Source_Type));
	if (p == NULL)
	  {
	     free_sources ();
	     return -1;
	  }
	Point_Sources = p;
     }
   
   p = Point_Sources + Num_Points;
   /* Note the the minus sign is to generate a vector pointing from the
    * source to the origin
    */
   p->cosx = -cos (dec) * cos (ra);
   p->cosy = -cos (dec) * sin(ra);
   p->cosz = -sin (dec);
   
   p->weight = weight;
   p->energy = energy;
   Num_Points += 1;
   
   return 0;
}

static void normalize_sources (void)
{
   double total;
   unsigned int i;

   total = 0;
   for (i = 0; i < Num_Points; i++)
     {
	Point_Sources[i].weight += total;
	total = Point_Sources[i].weight;
     }

   for (i = 0; i < Num_Points; i++)
     Point_Sources[i].weight /= total;
   
   /* Make sure no round-off error affects the weight of the last point */
   Point_Sources[Num_Points - 1].weight = 1.0;
}

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   FILE *fp;
   char line[1024];
   char *file;
   unsigned int linenum;

   file = argv[0];
   if (file == NULL)
     {
	fprintf (stderr, "UserSource Model requires FILE as argument\n");
	return -1;
     }

   fp = fopen (file, "r");
   if (fp == NULL)
     {
	fprintf (stderr, "Unable to open %s\n", file);
	return -1;
     }

   linenum = 0;
   while (NULL != fgets (line, sizeof (line), fp))
     {
	double ra, dec, weight, energy;
	
	linenum++;
	if (4 != sscanf (line, "%lf %lf %lf %lf", &ra, &dec, &weight, &energy))
	  continue;
	
	if (weight <= 0.0)
	  {
	     fprintf (stderr, "weight on line %d of %s must be positive\n",
		      linenum, file);
	     free_sources ();
	     return -1;
	  }

	if (-1 == add_source (ra, dec, weight, energy))
	  {
	     fclose (fp);
	     return -1;
	  }
     }
   
   fclose (fp);
   if (Num_Points == 0)
     {
	fprintf (stderr, "%s contains no sources\n", file);
	return -1;
     }

   normalize_sources ();
   return 0;
}

void user_close_source (void)
{
   free_sources ();
}


int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   double r;
   Point_Source_Type *p;
   
   p = Point_Sources;
   
   r = JDMrandom ();
   while (r > p->weight)
     p++;
   
   *delta_t = -1.0;
   *energy = p->energy;
   *cosx = p->cosx;
   *cosy = p->cosy;
   *cosz = p->cosz;

   return 0;
}

int main (int a, char **b)
{
   (void) a;
   (void) b;
   return 1;
}
