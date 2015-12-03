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
/* Note:  The code in this file is essentually a complete re-write of the
 * Lumb/Nousek/Townsley/Zacher CCD simulator.  The original simulator
 * was written in Fortran by Lumb and later modified by Nousek.  The Townsley
 * converted it to Fortran.  It was later converted to C by Bob Zacher.
 * Unfortunately, the original code had some bugs and those bugs were
 * propagated along with each translation.  Finally, the translations
 * themselves introduced new bugs.  The bugs I found were reported to Bob
 * Zacher at SAO.
 *
 * To produce this version, I went back to the original papers that the
 * previous version was based on and rewrote the simulator from scratch
 * by taking the ideas in the original papers as well as the ideas of the
 * Lumb/Nousek/Townsley/Zacher implementation.
 *
 * --John E. Davis
 * davis@space.mit.edu
 */
#include "config.h"

#include <stdio.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>
#include <jdmath.h>

#include <marx.h>
#include <pfile.h>

#include "acissim.h"
#include "henke.h"

#define USE_HOPKINSON_CHARGE_FRAC 1

/* Physical constants */

static double Fano_factor = 0.115;
/* The fano factor is defined as follows:
 * On the average, an X-ray of energy E will give rize to <N> e+e- pairs.
 * Then w = E / <N> is the average energy per pair.  In general, the number
 * <N> will fluctuate with some sigma dN.  The fano factor F is defined
 * by dN^2 = F<N>.   This means that in this monti-carlo simulation
 * we will assume that a photon of energy E will give rise to N e+e- pairs
 * where N = <N> + g dN
 * Here g is a gaussian distributed random number of unit variance.
 * Thus,
 * @       N = w<N>/w + g sqrt(Fw<N>/w)
 * @         = E/w + g sqrt (FE/w)
 * is the number of pairs liberated by the photon.
 * Note that F is unitless.
 */

static double Energy_Per_Pair = 3.68e-3;   /* KeV/e- */
/* Average energy per e+e- pair.  This is w in the above description.
 */

static double Fluor_Absorbtion_Energy = 1.839;
/* Absorbtion energy of fluorescent photon.  This is the binding energy of the
 * K shell electron.  It takes this amount of energy to remove the electron.
 */

static double Fluor_Emission_Prob = 0.043;
/* Fluorescence probability of silicon */

static double Fluor_Emission_Energy = 1.7396;
/* Energy of the fluoresced photon */

static double Fluor_Absorbtion_Coeff = 0.0816;
/* 1/absorbtion length of Si Ka */

static double Pixel_Size = 24.0;       /* microns */
static double Dead_Thickness = 0.6;
static double Field_Free_Thickness = 370.0;
static double Depletion_Thickness = 65.0;
static double Substrate_Thickness = 85.0;

static double Field_Free_Boundary;     /* depth of freefield/substrate interfacc */
static double Depletion_Boundary;      /* depth of depletion/freefield interface */
static double Substrate_Boundary;      /* depth of back of substrate */

static double Total_Thickness;

static double Diffusion_Constant = 35.0 * 1e8;   /* microns^2/sec */
static double Diffusion_Length = 75.0; /* microns */
static double Substrate_Diffusion_Length = 10.0;   /* microns */

static double Dopant_Concentration = 1e13 * 1e-12;   /* per um^3 */

static unsigned int Island_Size = 3;

static int apply_filter_qe (double);

#define JANESICK_FACTOR 1.0

/* Background.
 *
 * An xray photon will give up its energy via the photoelectric
 * effect a kick out an electron.  This electron will scatter giving
 * up its energy to electron-hole pairs.  This will happen at some
 * depth in the ccd determined by an energy-dependent absorbtion
 * length calculated below. This gives rise to an initial charge
 * cloud of pairs whose size is energy dependent. The size of the
 * cloud is given in Hopkinson 87 by
 * @  R = 0.5 * 0.0176 E^1.75
 * where E is the energy of the liberated electron.  The electron energy
 * will differ from the XRay energy-- the binding energy of K shell electrons
 * in silicon is 1.78 KeV.
 * The factor of 0.5 is necessary to convert diameter to radius.
 *
 * It is also possible for the initial photon to give rise to a fluorescence
 * photon.  This photon is also included in the simulation.
 *
 * The charge cloud then begins to diffuse until it has reached the channel
 * where it will be ``detected''.  The actual rate of diffusion will depend
 * upon the properties of the layer that it is travelling through.  It is
 * possible for the charge cloud to spread across several pixels before
 * detection.
 */

/* This routine is specific to silicon only.  It comes from Townsley's
 * code which was derived from a Fortran program by Ralph Kraft.
 */
static double absorbtion_coeff (double energy)
{
   double x;

   if ((energy < 0.03) || (energy > 14.0))
     fprintf(stderr, "**Warning: absorbtion_coeff: energy (%e) out of range.\n",
	     energy);

   x = log10 (energy);
   x += 3.0;			       /* convert to ev */

   if (energy < 0.1006)
     {
	x = -1.4716630 + x * (7.158886 + x * (-2.258985));
     }
   else if (energy < 0.16)
     {
	x = 7.146932 + x * (-1.752903 + x * (0.3575078));
     }
   else if (energy < 1.84)
     {
	x = -421.2560 + x * (800.2535 + x *
			     (-596.4018 + x *
			      (221.1723 + x * (-40.89891 + x * 3.014328))));
     }
   else				       /* < 14.0 */
     {
	x = -80.59615 + x * (91.47612 + x * (-35.69339 + x *
					     (6.007395 + x * (-0.3795995))));
     }

   x = pow(10.0, x);	/* absorption coeff in units of cm^2/gram. */

   /* Now we need to change units from cm^2/gm to 1/um.  We need the density
    * of silicon (2.33 gm/(cm^3)) to do this.  Also 1e-4 factor is conversion
    * from cm to microns.
    */

   return x * 2.33 * 1.0e-4;
}

/* In this function, xoff is the offset from the middle cell (ncells/2).
 * If the photon lands in the middle of the cell, xoff would be 0.5.  That is,
 * the x here is in pixel units whose size is pixel_size microns.
 *
 * ncells should be odd.  This is not checked.
 *
 * This routine is called alot.  Profiling indicates the program spends
 * about 10% of its time here.  For that reason the erf(x) call is
 * avoided if abs(x) is large.
 */
static void integrate_erf (double sigma, double *prob, unsigned int ncells,
			   double xoff, double pixel_size)
{
   double xmin, xmax;
   unsigned int i;
   double last_val, val;

   if (sigma == 0.0)
     {
	for (i = 0; i < ncells; i++) prob[i] = 0;
	prob [ncells/2] = 1.0;
	return;
     }

   pixel_size = pixel_size / sigma;

   xmin = -((ncells / 2) + xoff) * pixel_size;
   xmax = xmin + ncells * pixel_size;

   if (xmin < -4.0) last_val = 0.0;
   else if (xmin > 4.0) last_val = 1.0;
   else last_val = 0.5 * (1 + JDMerf (xmin));

   for (i = 0; i < ncells; i++)
     {
	val = xmin + (i + 1) * pixel_size;

	if (val > 4.0) val = 1.0;
	else if (val < -4.0) val = 0.0;
	else if (val == 0.0) val = 0.5;
	else val = 0.5 * (1 + JDMerf (val));

	prob [i] = val - last_val;
	last_val = val;
     }
}

static int compute_phas (double energy, double sigma, double charge_fraction,
			 AcisSim_Pixel_Island_Type *island)
{
   double yprobs [MAX_ISLAND_SIZE];
   double xprobs [MAX_ISLAND_SIZE];
   unsigned int i, j;
   double *phas;
   unsigned int island_size;

   /* Assume that sigma is the 1-sigma radius.  That is, it is the radius of
    * the charge cloud that hold 63% of the charge.  This is consistent with
    * a 2-d gaussian probability distribution of the form
    * @ Q(R) = \frac{1}{\pi \sigma^2} \int_0^R \dr 2\pi r exp (-r^2/\sigma^2)
    * Then, Q(\infty) = 1 and Q(\sigma) = 1 - 1/e
    *
    * In terms of cartesian coordinates:
    *
    * @ Q(X,Y) = 1/4 (1 + \erf(X/\sigma)) (1 + \erf(Y/\sigma))
    * where Q(X,Y) represents the charge fraction for
    * @ x < X and y < Y.
    */

   /* Here, the number of electrons is expressed as an energy via the relation
    * Energy = Energy_Per_Pair * Num_Pairs.
    * Thus, was simply multiply the energy by the charge-fraction to get the
    * fractional number of pairs.
    */
   energy = charge_fraction * energy;

   integrate_erf (sigma, xprobs, Island_Size, island->x - (int) island->x,
		  Pixel_Size);
   integrate_erf (sigma, yprobs, Island_Size, island->y - (int) island->y,
		  Pixel_Size);

   phas = island->phas;
   island_size = island->island_size = Island_Size;

   for (i = 0; i < island_size; i++)
     {
	double yprob = yprobs [i];

	/* Note: the memset near the top level of these routines
	 * ensures that island is zeroed.
	 */
	if (yprob != 0.0) for (j = 0; j < island_size; j++)
	  {
	     double energy_times_prob;
	     double pha;
	     double prob = yprob * xprobs [j];

	     /* The stddev of number of electrons in this cell is
	      * given by sqrt (n_electrons * prob * (1 - prob)) and there
	      * are Energy_Per_Pair KeV per electron.  Thus dE is
	      * Energy_Per_Pair * sqrt (n_electrons * prob * (1-prob))
	      * == sqrt (n_electrons * Energy_Per_Pair^2 prob * (1-prob))
	      * == sqrt (energy * prob * Energy_Per_Pair * (1 - prob))
	      */

	     if (prob == 0.0)
	       continue;

	     energy_times_prob = energy * prob;
	     pha = energy_times_prob
	       + JDMgaussian_random () * sqrt (energy_times_prob * Energy_Per_Pair * (1 - prob));

	     if (pha > 0.0)
	       phas [j] = pha;
	  }

	phas += island_size;
     }

   island->radius = sigma / Pixel_Size;
   island->charge_fraction = charge_fraction;
   return 0;
}

/* Here we assume that a photon depth comes from an exponential distribution
 */
static double compute_depth (double abs_coeff)
{
   double r;

   do
     {
	r = JDMrandom ();
     }
   while (r == 0.0);
   return log (r) / -abs_coeff;
}

#if !USE_HOPKINSON_CHARGE_FRAC

/* This function returns the square of the 1-sigma radius */
static double substrate_spread (double z)
{
   /* This routine solves equation 19 of Hopkinson via Newton's method.
    * The idea is to find the radius that encloses contains 60% of the charge
    * within the cloud.  This means that the equation should be integrated
    * between r = 0 and r = R to find the charge inside a circle of radius
    * R.  In the end, an equation is obtained of the form
    * @ Q(R) = Q_0 e^(-z/L) - Q_0 (z/L) e^(-u_R)/u_R
    * where u_R^2 = (1/L^2)(R^2 + z^2)
    * The total charge Q(inf) = Q_0 e^(-z/L)
    * Thus, the ratio is
    * @ Q(R)/Q(inf) = 1 - Q_0 (z/L) e^(-u_R + z/L)/u_R
    * or
    * @ e^(-u_R)/u_R = (1 - f)(L/z) e^(-z/L)
    * where f is the ratio Q(R)/Q(inf).
    * This equation will be solved via Newton's method.  It should converge
    * rather rapidly.
    */
   double u_R, new_u_R;
   double alpha;
   double f = 0.6;
   double eps = 1.0e-6, diff;
   unsigned int max_its;

   if ((z <= 0.0) || (Substrate_Diffusion_Length <= 0.0))
     return 0.0;

   u_R = z / Substrate_Diffusion_Length;   /* seed */
   alpha = (1 - f) * (Substrate_Diffusion_Length / z) * exp (-u_R);

   max_its = 20;
   do
     {
	double ex = exp (-u_R);
	new_u_R = (1.0 + u_R * ex) / (alpha + ex);
	diff = fabs (new_u_R - u_R);
	u_R = new_u_R;
	max_its--;
	if (max_its == 0)
	  {
	     if (diff > eps)
	       fprintf (stderr, "substrate_spread: failed to converge.\n");
	     break;
	  }
     }
   while (diff > eps);

   u_R = Substrate_Diffusion_Length * u_R;
   return u_R * u_R - z * z;
}
#endif				       /* NOT USE_HOPKINSON_CHARGE_FRAC */

static double initial_cloud_size (double energy)
{
   return (0.5 * 0.0176) * pow (energy, 1.75);
}

/* For an epitaxial fron-illuminated ccd, we have layers:
 *   dead layer (where gate structure is located)
 *   depletion region
 *   field-free region
 *   substrate
 */
static int epi_spread (double energy, double depth, double *radius, double *cfrac)
{
   double mobility = 1500.0 * 1.0e8;  /* convert from cm^2/V/s to um/V/sec */
   double echarge = 1.602e-19;	       /* coul */
   double epsilon = 1.044e-16;	       /* permit of silicon in F/um */
   double dz_cutoff;
   double dz;
   double r_depletion, r_ff;
   double r_squared;
   double charge_fraction;

   if (depth < Dead_Thickness)
     return -1;

   r_squared = initial_cloud_size (energy);
   r_squared = r_squared * r_squared;

   /* We know that it at least made it through the depletion region and
    * the cloud must spread there from the other layers.
    */

   /* In the depletion region, use equation 7 from Hopkinson.  If the photon
    * is too close to the boundary, the equation will not hold so use
    * a cutoff suggested by Hopkinson.  The cutoff also enters as a
    * factor in eq 7 of Hopkinson.
    */

   dz_cutoff = (Diffusion_Constant * epsilon) / (mobility * echarge * Dopant_Concentration);
   dz_cutoff = sqrt (dz_cutoff);

   dz = Depletion_Boundary - depth;
   if (dz < dz_cutoff) dz = dz_cutoff;

   r_depletion = dz_cutoff * sqrt (2.0 * log (Depletion_Boundary / dz));
   r_squared += r_depletion * r_depletion;

   charge_fraction = 1.0;

   /* Now check contribution from the other layers. */
   if (depth >= Depletion_Boundary)
     {
	double ratio;

	/* In the field free layer, use equation 7 from Janesick. */
	dz = Field_Free_Boundary - depth;
	if (dz < 0.0) dz = 0;
	ratio = dz / Field_Free_Thickness;

	/* Equation 7 follows below.  Note the factor of 0.5 that I added
	 * to account for the 2 sigma diameter --> 1 sigma radius.
	 *
	 * Another note: Bob Zacher reports that Janesick's 2-sigma diameter
	 * is not 2 * (2 * sigma) radius.  My own monti-carlo simulations
	 * indicate that it is a 1-sigma diameter.
	 */
	r_ff = JANESICK_FACTOR * Field_Free_Thickness * sqrt (1.0 - ratio * ratio);
	r_squared += r_ff * r_ff;

	/* Did the charge cloud start from the substrate?  If so, add its
	 * effect.  If not, be sure to compute the charge fraction lost
	 * here due to recombination
	 */
	if (depth < Field_Free_Boundary)
	  {
#if USE_HOPKINSON_CHARGE_FRAC
	     /* Event not in substrate.  Use equation 8 from Hopkinson.
	      * In it, take the reflection coeff R = 1 and the transmission
	      * coeff T = 0 since we expect the charge cloud to get reflected
	      * from the substrate boundary.
	      */
	     charge_fraction = cosh (dz / Diffusion_Length) / cosh (Field_Free_Thickness / Diffusion_Length);
#endif
	  }
	else return -1;
#if 0
	  {
	     /* We are in the substrate.  Equation 8 from Hopkinson holds
	      * again except now take R = 0 and T = 1.
	      */
	     dz = Substrate_Boundary - depth;
#if USE_HOPKINSON_CHARGE_FRAC
	     charge_fraction = sinh (dz / Substrate_Diffusion_Length) / sinh(Substrate_Thickness / Substrate_Diffusion_Length);
#endif
	     r_squared += substrate_spread (Substrate_Thickness - dz);
	  }
#endif
     }

   *radius = sqrt (r_squared);
   *cfrac = charge_fraction;

   return 0;
}

static int (*Spread_Function) (double, double, double *, double *);

static int charge_cloud_spread (double energy, double depth,
				AcisSim_Pixel_Island_Type *island)
{
   double radius, charge_fraction;

   /* Call the appropriate function to spread the cloud. */
   if (-1 == (*Spread_Function)(energy, depth, &radius, &charge_fraction))
     return -1;

   /* Convert the energy to an equivalent number of e+e- pairs */
   energy += JDMgaussian_random () * sqrt (Fano_factor * Energy_Per_Pair * energy);

   compute_phas (energy, radius, charge_fraction, island);

   return 0;
}

static int handle_fluorescent_photon (AcisSim_Ray_Type *ray, double energy, double depth,
				      AcisSim_Pixel_Island_Type *island)
{
   double theta, phi;
   double r;
   /* The fluorescent photon is going to travel some distance before it
    * stops.  First, compute its direction.
    */

   theta = PI * JDMrandom ();
   phi = 2.0 * PI * JDMrandom ();

   do
     {
	r = JDMrandom ();
     }
   while (r == 0.0);
   r = log (r) / -Fluor_Absorbtion_Coeff;

   depth += r * cos (theta);
   if ((depth < 0.0) || (depth >= Total_Thickness))
     return -1;

   r = r * sin (theta);

   r = r / Pixel_Size;

   island->x = ray->xpixel + r * cos (phi);
   island->y = ray->ypixel + r * sin (phi);

   return charge_cloud_spread (energy, depth, island);
}

static int handle_regular_event (AcisSim_Ray_Type *ray, double energy, double depth,
				 AcisSim_Pixel_Island_Type *island)
{
   island->x = ray->xpixel;
   island->y = ray->ypixel;

   return charge_cloud_spread (energy, depth, island);
}

static int process_ray_1 (AcisSim_Ray_Type *ray, AcisSim_Pixel_Event_Type *event,
			  double abs_coeff)
{
   double depth;
   double energy = ray->energy;

   depth = compute_depth (abs_coeff);

   if (depth >= Total_Thickness)
     return -1;

   /* Routines that this call assumes that this structure has been zeroed. */
   memset ((char *) event, 0, sizeof (AcisSim_Pixel_Event_Type));

   /* Does this photon give rise to a fluorescent photon? */
   if (0 != (event->did_fluoresc = ((energy > Fluor_Absorbtion_Energy)
				    && (Fluor_Emission_Prob >= JDMrandom ()))))
     {
	/* If it fluoresced, then we also need to consider what happened to
	 * the leftover energy.
	 */
	energy -= Fluor_Absorbtion_Energy;
	if (0 == handle_fluorescent_photon (ray, Fluor_Emission_Energy,
					    depth, &event->fluoresc_island))
	  event->flags |= FLUOR_EVENT_OK;
     }

   if (0 == handle_regular_event (ray, energy, depth, &event->regular_island))
     event->flags |= REGULAR_EVENT_OK;

   if (event->flags & (REGULAR_EVENT_OK | FLUOR_EVENT_OK))
     return 0;

   return -1;
}

int acissim_process_ray (AcisSim_Ray_Type *ray,
			 AcisSim_Pixel_Event_Type *event)
{
   double energy;
   double abs_coeff;

   energy = ray->energy;

   if (-1 == apply_filter_qe (energy))
     return -1;

   abs_coeff = absorbtion_coeff (energy);

   return process_ray_1 (ray, event, abs_coeff);
}

static char *Henke_Dir;

static double Lexan_Thickness;
static double Lexan_Density;
static double Aluminum_Thickness;
static double Aluminum_Density;
static double SiO2_Density;
static double SiO2_Thickness;

static float *Filter_Trans_Coeffs;
static float *Filter_Energies;
static unsigned int Num_Filter_Energies;

static int apply_filter_qe (double energy)
{
   float t;

   if (Filter_Trans_Coeffs == NULL)
     return 0;

   t = JDMlog_interpolate_f (energy, Filter_Energies, Filter_Trans_Coeffs, Num_Filter_Energies);

   if (JDMrandom () < t)
     return 0;

   return -1;
}

static int setup_filter (void)
{
   Henke_Type *al, *lexan, *sio2;
   float *sio2_betas;
   float *lexan_betas;
   float *al_betas;
   float emin, emax;
   float energy;
   unsigned int i;
   unsigned int num_sio2, num_lexan, num_al;

   sio2 = lexan = al = NULL;

   if ((Lexan_Thickness > 0.0) && (Lexan_Density > 0.0))
     {
	if (NULL == (lexan = henke_read_henke_table ("lexan")))
	  {
	     fprintf (stderr, "Unable to read henke tables for lexan.\n");
	     return -1;
	  }

	if (-1 == henke_beta_delta (lexan, Lexan_Density, &lexan_betas, NULL))
	  {
	     fprintf (stderr, "Error computing BETA for lexan.\n");
	     goto return_error;
	  }
     }

   if ((Aluminum_Thickness > 0.0) && (Aluminum_Density > 0.0))
     {
	if (NULL == (al = henke_read_henke_table ("Al")))
	  {
	     fprintf (stderr, "Unable to read henke tables for Al.\n");
	     goto return_error;
	  }

	if (-1 == henke_beta_delta (al, Aluminum_Density, &al_betas, NULL))
	  {
	     fprintf (stderr, "Error computing BETA for lexan.\n");
	     goto return_error;
	  }
     }

   if ((SiO2_Thickness > 0.0) && (SiO2_Density > 0.0))
     {
	if (NULL == (sio2 = henke_read_henke_table ("SiO2")))
	  {
	     fprintf (stderr, "Unable to read henke tables for SiO2.\n");
	     goto return_error;
	  }

	if (-1 == henke_beta_delta (sio2, SiO2_Density, &sio2_betas, NULL))
	  {
	     fprintf (stderr, "Error computing BETA for SiO2.\n");
	     goto return_error;
	  }

     }

   if ((lexan == NULL) && (al == NULL) && (sio2 == NULL))
     return 0;

   Num_Filter_Energies = 1024;
   if ((NULL == (Filter_Energies = JDMfloat_vector (Num_Filter_Energies)))
       || (NULL == (Filter_Trans_Coeffs = JDMfloat_vector (Num_Filter_Energies))))
     goto return_error;

   emin = 0.03;			       /* see absorbtion_coeff */
   emax = 14.0;

   if (sio2 != NULL)
     {
	num_sio2 = sio2->num_elements;
	emin = sio2->energy[0];
	emax = sio2->energy[num_sio2 - 1];
     }
   else num_sio2 = 0;

   if (al != NULL)
     {
	num_al = al->num_elements;
	if (emin < al->energy[0]) emin = al->energy[0];
	if (emax > al->energy[num_al - 1]) emax = al->energy[num_al - 1];
     }
   else num_al = 0;

   if (lexan != NULL)
     {
	num_lexan = lexan->num_elements;
	if (emin < lexan->energy[0])
	  emin = lexan->energy [0];
	if (emax > lexan->energy[num_lexan - 1])
	  emax = lexan->energy [num_lexan - 1];
     }
   else num_lexan = 0;

   if (emin <= 0.03) emin = 0.0301;       /* see absorbtion_coeff */
   if (emax > 14.0) emax = 14.0;

   (void) JDMlog_grid_f (Filter_Energies, Num_Filter_Energies, emin, emax);

   for (i = 0; i < Num_Filter_Energies; i++)
     {
	double a;
	float beta;
	double k;
	double factor = 2.0 / 1.9732858e-4; /* inv of 1/2 hbar c (KeV-microns) */

	energy = Filter_Energies [i];
	k = energy * factor;

	a = 0.0;

	if (al_betas != NULL)
	  {
	     beta = JDMlog_interpolate_f (energy, al->energy, al_betas, num_al);
	     a += beta * Aluminum_Thickness;
	  }

	if (lexan_betas != NULL)
	  {
	     beta = JDMlog_interpolate_f (energy, lexan->energy, lexan_betas, num_lexan);
	     a += beta * Lexan_Thickness;
	  }

	if (sio2_betas != NULL)
	  {
	     beta = JDMlog_interpolate_f (energy, sio2->energy, sio2_betas, num_sio2);
	     a += beta * SiO2_Thickness;
	  }

	Filter_Trans_Coeffs [i] = exp (-(k * a));
     }

   if (lexan != NULL) henke_free_henke_table (lexan);
   if (al != NULL) henke_free_henke_table (al);
   if (al_betas != NULL) free ((char *) al_betas);
   if (sio2 != NULL) henke_free_henke_table (sio2);
   if (sio2_betas != NULL) free ((char *) sio2_betas);
   if (lexan_betas != NULL) free ((char *) lexan_betas);
   return 0;

   return_error:

   if (lexan != NULL) henke_free_henke_table (lexan);
   if (al != NULL) henke_free_henke_table (al);

   if (al_betas != NULL) free ((char *) al_betas);
   if (lexan_betas != NULL) free ((char *) lexan_betas);

   if (sio2 != NULL) henke_free_henke_table (sio2);
   if (sio2_betas != NULL) free ((char *) sio2_betas);

   if (Filter_Energies != NULL)
     JDMfree_float_vector (Filter_Energies);
   if (Filter_Trans_Coeffs != NULL)
     JDMfree_float_vector (Filter_Trans_Coeffs);

   return -1;
}

static char *Illumination_Type;

static Param_Table_Type Acis_Parm_Table [] =
{
     {"Illum",			PF_STRING_TYPE,	&Illumination_Type},
     {"IslandSize",		PF_INTEGER_TYPE,&Island_Size},
     {"PixelSize",		PF_REAL_TYPE,	&Pixel_Size},
     {"DeadSize",		PF_REAL_TYPE,	&Dead_Thickness},
     {"DepletionSize",		PF_REAL_TYPE,	&Depletion_Thickness},
     {"FieldFieldSize",		PF_REAL_TYPE,	&Field_Free_Thickness},
     {"SubstrateSize",		PF_REAL_TYPE,	&Substrate_Thickness},
     {"DiffusionLength",	PF_REAL_TYPE,	&Diffusion_Length},
     {"SubstDiffusionLength",	PF_REAL_TYPE,	&Substrate_Diffusion_Length},
     {"DiffusionConst",		PF_REAL_TYPE,	&Diffusion_Constant},
     {"ImpurityConcent",	PF_REAL_TYPE,	&Dopant_Concentration},

     {"SiO2Thickness",		PF_REAL_TYPE,	&SiO2_Thickness},
     {"SiO2Density",		PF_REAL_TYPE,	&SiO2_Density},

     {"LexanThickness",		PF_REAL_TYPE,	&Lexan_Thickness},
     {"LexanDensity",		PF_REAL_TYPE,	&Lexan_Density},
     {"AluminumThickness",	PF_REAL_TYPE,	&Aluminum_Thickness},
     {"AluminumDensity",	PF_REAL_TYPE,	&Aluminum_Density},
     {"HenkeDir",		PF_FILE_TYPE,	&Henke_Dir},

     {"FanoFactor",		PF_REAL_TYPE,	&Fano_factor},
     {"PairEnergy",		PF_REAL_TYPE,	&Energy_Per_Pair},
     {"FluorAbsorbtionEnergy",	PF_REAL_TYPE,	&Fluor_Absorbtion_Energy},
     {"FluorEmissionEnergy",	PF_REAL_TYPE,	&Fluor_Emission_Energy},
     {"FluorProbability",	PF_REAL_TYPE,	&Fluor_Emission_Prob},
     {NULL, 0, NULL}
};

int acissim_init (Param_File_Type *p)
{
   char *dir;
   char dirbuf [1024];

   if (-1 == pf_get_parameters (p, Acis_Parm_Table))
     {
	pf_error ("acissim_init: error getting parameters.");
	return -1;
     }
   dir = Henke_Dir;

   if (*dir == '$')
     {
	char *env_end;
	char *env = dir + 1;

	env_end = strchr (env, '/');
	if (env_end != NULL)
	  {
	     unsigned int env_len = env_end - env;
	     strncpy (dirbuf, env, env_len);
	     dirbuf [env_len] = 0;
	     env = dirbuf;
	  }

	fprintf (stdout, "Setting HENKE data directory from %s environment variable\n", env);

	dir = getenv (env);
	if (dir == NULL)
	  {
	     fprintf (stderr, "Environment variable %s needs set\n", env);
	     return -1;
	  }

	if (env_end != NULL)
	  {
	     sprintf (dirbuf, "%s%s", dir, env_end);
	     dir = dirbuf;
	  }
     }

   fprintf (stdout, "Setting Henke Directory to %s.\n", dir);
   if (-1 == henke_set_data_dir (dir))
     {
	fprintf (stderr, "Unable to set Henke Directory.\n");
	return -1;
     }

   if (-1 == setup_filter ())
     return -1;

   Diffusion_Constant = Diffusion_Constant * 1.0e8;   /* from cm^2 to um^2 */
   Dopant_Concentration = Dopant_Concentration * 1.0e-12;   /* to 1/um^3 */

   Depletion_Boundary = Dead_Thickness + Depletion_Thickness;
   Field_Free_Boundary = Depletion_Boundary + Field_Free_Thickness;
   Substrate_Boundary = Field_Free_Boundary + Substrate_Thickness;

   Total_Thickness = Substrate_Boundary;

   Spread_Function = epi_spread;

   Fluor_Absorbtion_Coeff = absorbtion_coeff (Fluor_Emission_Energy);
   return 0;
}

