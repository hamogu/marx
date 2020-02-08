/* file: Draine_halo.c
 * author: Andrea Belfiore <andrea.belfiore@inaf.it>
 * date: 2018-10-23
 * github: https://github.com/andrea-belfiore/MARX-plugins
 * description:
 *   MARX plugin describing the brightness distribution expected from a dust
 *   scattering halo according to the model of Draine, 2003, ApJ 598, 1023.
 *   The model parameters are set, at the beginning of the script, to those
 *   of NGC5907 ULX-1 in Chandra observation 20830, under the assumptions
 *   described in Belfiore et al. (2019).
 * usage:
 *   this plugin can be compiled with a line like the following:
gcc -shared -fPIC -I${MARX_DIR}/include -L${MARX_DIR}/lib -ljdmath Draine_halo.c -o Draine_halo.so
 *   afterwards it can be called from MARX, with the following options:
marx SourceType=USER UserSourceFile=./Draine_halo.so UserSourceArgs="$d_10" ... (other options)
marx SourceType=USER UserSourceFile=./Draine_halo.so UserSourceArgs="$d_10 $E_1" ... (other options)
 *   the argument $d_10 is the distance between ULX1 and the dust layer
 *   in 10 kpc units, and $E_1 is the characteristic photon energy in keV,
 *   used to determine the differential scattering cross section.
 *   E_1 can be omitted, in which case, a default E_1=1 is taken.
 */
#include <stdio.h>
#include <math.h>
#include <jdmath.h>

#include "user.h"

/* This value applies to NGC 5907 (Tully et al, 2013 AJ 146, 86) */
static double D_10 = 1710.; // distance to the source in 10 kpc units
/* This value applies to Chandra obs. 20830 (Belfiore et al 2019) */
static double dt = 120.; // time since switch off in days
static double d_10, E_1; // input parameters
static double r_50, r_min; // halo shape parameters

/* Useful conversion factors */
static double days2pc = 0.3066014 / 365.25;
static double as2rad = M_PI / 180. / 3600.;
/* Source direction */
static double source_x, source_y;

/* Constructor routine, called when MARX is invoked */
int user_open_source (char **argv, int argc, double area,
    double cosx, double cosy, double cosz) {
  char *d_10_str, *E_1_str; // input parameter strings
  switch (argc) {
    case 1:
      d_10_str = argv[0];
      if (1 != sscanf (d_10_str, "%lf", &d_10))
        return Draine_halo_usage();
      E_1 = 1.;
      break;
    case 2:
      d_10_str = argv[0];
      E_1_str = argv[1];
      if (1 != sscanf (d_10_str, "%lf", &d_10))
        return Draine_halo_usage();
      if (1 != sscanf (E_1_str, "%lf", &E_1))
        return Draine_halo_usage();
      break;
    default:
      return Draine_halo_usage();
  }
  fprintf(stdout, "Generate a dust halo profile for NGC 5907 ULX-1 with d=%lf kpc at E=%lf keV\n", 10.*d_10, E_1);

  /* From Draine, 2003, ApJ 598, 1023 - eqq. 13 & 17
   * P(X<r) = A / (1 + A)  where A={r/[(1-x)*th_50]}**2
   * where (Draine's eq 9)  th_50 = 360" / E_1
   * it is convenient to define
   * r_50 = (1-x) * th_50  so that  A = (r/r_50)**2 
   */
  r_50 = (d_10 / D_10) *  360. / E_1;
  /* From Belfiore et al 2019 - eq. 13
   * B(r)=0 for r < r_min = sqrt(2c * dt * (1-x) / D)
   * note: the factor e-4 accounts for the fact that all distances
   * in this equation (d_10 and D_10) are in units of 10 kpc
   */
  r_min = sqrt(2.e-4 * days2pc * dt * d_10 / pow(D_10, 2.)) / as2rad;
  fprintf(stdout, "Dust halo shape parameters: r_50=%lf as and r_min=%lf as\n", r_50, r_min);

  /* store for later the source direction */
  source_y = asin(cosz);
  source_x = atan2(cosy, cosx);
  return 0;
}

/* Destructor routine, called before leaving */
void user_close_source (void) {
}

/* Main routine, called for each photon to be generated */
int user_create_ray (double *delta_t, double *energy,
    double *cosx, double *cosy, double *cosz) {
  double rand1; // temporary variable for the radial profile
  double r, phi; // circular coordinates TAN to the source
  double x, y, tmp; // CXO coordinates
  /* Return values are passed by reference for the 3 directional
   * cosines, time and energy. setting them to -1 means that the
   * default coming from other parameters passed to MARX are not
   * overridden. here we care only, as a first approximation,
   * for the spatial distribution of the source.
   */
  *delta_t = -1.0; // it means: do not override this parameter
  *energy = -1.0; // it means: do not override this parameter
  
  rand1 = JDMrandom();
  /* Compute the CDF of detecting a photon within a radius r adding
   * a hole with radius r_min to the Draine's profile (P(X<r)):
   * P = cdf(r) = {0 r<r_min; [A(r)-A(r_min)]/[1+A(r)] r>R_min}
   * that can be inverted as:
   * r = qf(P) = r_50 * sqrt((P + A(r_min))/(1 - P))
   * (mind that, as expected, qf(0) = r_min and qf(1) = infin )
   */
  r = r_50 * sqrt((rand1 + pow(r_min / r_50, 2.)) / (1 - rand1));

  /* Pick a uniformly random azimuth (central symmetry) */
  phi = 2. * M_PI * JDMrandom();

  /* All r's so far were expressed in arcsec */
  x = source_x + r * cos(phi) * as2rad;
  y = source_y + r * sin(phi) * as2rad;

  /* Switch from the TAN plane to the CXO directional cosines */
  tmp = cos(y);
  *cosz = sin(y);
  *cosx = tmp * cos(x);
  *cosy = tmp * sin(x);
  //fprintf(stdout, "raytrace an event with (r, phi)=(%lf as, %lf rad) to ray: cosx=%.16f cosy=%.16f cosz=%.16f\n", r, phi, *cosx, *cosy, *cosz);
  return 0;
}

/* Secondary routine, called from within user_create_ray() */
int Draine_halo_usage(void) {
   fprintf(stderr, "Draine (2003) Scattering Halo User Source Model Usage:\n");
   fprintf(stderr, "\targs: \"d_10[ E_1]\"\n");
   fprintf(stderr, "\td_10 distance between the source and the dust layer in 10 kpc units\n");
   fprintf(stderr, "\tE_1 characteristic photon energy for differential cross section in keV (default 1)\n");
   return -1; // If the user gets here, we take it as an error
}
