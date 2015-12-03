#include <stdio.h>
#include <jdmath.h>
#include "user.h"

static double Source_CosX;
static double Source_CosY;
static double Source_CosZ;

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   Source_CosX = cosx;
   Source_CosY = cosy;
   Source_CosZ = cosz;

   return 0;
}

void user_close_source (void)
{
}

int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   *cosx = Source_CosX;
   *cosy = Source_CosY;
   *cosz = Source_CosZ;

   *delta_t = -1.0;
   return 0;
}

static double Theta;
int user_start_iteration (void)
{
   double t;
   
   if (Theta > 20.0)
     return -1;

   t = Theta * (PI/(60.0 * 180.0));

   Source_CosX = -cos (t);
   Source_CosY = sin (t);
   Source_CosZ = 0.0;
   
   Theta += 0.25;		       /* 1/4 an arc-minute */
   return 0;
}


int main (int a, char **b)
{
   (void) a;
   (void) b;
   return 1;
}
