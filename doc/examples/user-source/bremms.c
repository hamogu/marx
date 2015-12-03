#include <stdio.h>
#include <jdmath.h>
#include "user.h"

static double Source_CosX;
static double Source_CosY;
static double Source_CosZ;
static double Emin = 0.1;
static double Emax = 12.0;
static double Brems_kT = 2.0;
static double Brems_Amp = 1.0;
static double *Cum_Spect;
static double *Cum_Energies;
static unsigned int Num_Cum;


int user_open_source (char *cmd_args, double area,
		      double cosx, double cosy, double cosz)
{
   Source_CosX = cosx;
   Source_CosY = cosy;
   Source_CosZ = cosz;
   
   if (Cum_Energies != NULL)
     {
	free ((char *) Cum_Energies);
	
   return 0;
}

void user_close_source (void)
{
}

static define bremsstrahlung (double e, double amp, double kT)
{
   return amp * exp (-e/kT);
}


int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   *cosx = Source_CosX;
   *cosy = Source_CosY;
   *cosz = Source_CosZ;

   *delta_t = -1.0;
   *energy = -1.0;
   
   return 0;
}

int main (int a, char **b)
{
   (void) a;
   (void) b;
   return 1;
}
