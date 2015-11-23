#include <stdio.h>
#include <jdmath.h>

int main (int argc, char **argv)
{
   unsigned int num = 10000;
   double theta, phi, theta1, phi1;
   JDMVector_Type p, p1;
   unsigned int num_errors;
   
   (void) argc;

   num_errors = 0;
   while (num)
     {
	double error;

	phi = 2.0 * PI * JDMrandom ();
	theta = PI * JDMrandom ();
	
	p.x = sin (theta) * cos (phi);
	p.y = sin (theta) * sin (phi);
	p.z = cos (theta);
	
	JDMv_unit_vector_to_spherical (p, &theta1, &phi1);
	p1 = JDMv_spherical_to_vector (1.0, theta, phi);
	
	error = JDMv_distance (p, p1);
	if (error > 1e-6)
	  {
	     fprintf (stderr, "***Error: %e (theta=%e, phi=%e)\n", 
		      error, theta, phi);
	     num_errors++;
	  }

	num--;
     }
   
   if (num_errors == 0)
     fprintf (stderr, "%s passed.\n", argv[0]);
   else
     fprintf (stderr, "%s failed.\n", argv[0]);
     
   return 0;
}

	
	

	
   
   
