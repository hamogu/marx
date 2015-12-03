#include <stdio.h>
#include <jdmath.h>
#include "user.h"

static double Source_CosX;
static double Source_CosY;
static double Source_CosZ;

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   return 0;
}

void user_close_source (void)
{
}

static double To_Radians = (PI / 180.0 / 3600.0);
#define ARC_SECONDS_PER_CELL 50

int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   static int last_i;
   double theta, phi;
   double cos_theta, sin_theta;

   if (last_i == 20) last_i = -20;
   
   theta = To_Radians * last_i * ARC_SECONDS_PER_CELL;
   phi = 2.0 * PI * JDMrandom ();

   sin_theta = sin(theta);

   *cosx = -cos (theta);
   *cosy = sin_theta * cos (phi);
   *cosz = sin_theta * sin (phi);

   *delta_t = -1.0;
   
   last_i++;  

   return 0;
}

int main (int a, char **b)
{
   (void) a;
   (void) b;
   return 1;
}
