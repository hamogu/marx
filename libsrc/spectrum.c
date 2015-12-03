/* -*- mode: C; mode: fold; -*- */
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
#include "config.h"
#include "marx-feat.h"

/*{{{ #include files */

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#include <string.h>

#include "marx.h"
#include "_marx.h"

/*}}}*/

/*{{{ compute_cumulated_vector */

/* Returns cumulated spectrum.  If parameter 'total' is non-NULL, the
 * spectrum will be normalized and *total will be set to normalization
 * constant (integrated vector).
 */
static int compute_cumulated_vector (double *cum,
				     double *x, double *y,
				     unsigned int n,
				     double *total)
{
   unsigned int i;
   double tot, this_x, prev_x;

   /* This performs a simple naive integration.  Something better may be
    * needed.
    */
   prev_x = x[0];
   tot = 0.0;
   for (i = 0; i < n; i++)
     {
	double diff;

	this_x = x[i];
	diff = (this_x - prev_x);
	if (diff < 0.0)
	  {
	     marx_error ("Your file spectrum appears corrupt.  Perhaps you are using a wavelength grid");
	     return -1;
	  }

	tot += y[i] * (this_x - prev_x);
	cum[i] = tot;
	prev_x = this_x;
     }
   if (tot == 0.0)
     {
       marx_error("Your file spectrum seems to have zero flux in it. MARX needs a positive integrated flux for a spectrum.");
       return -1;
     }

   if (total != NULL)
     {
	if (tot != 0.0) for (i = 0; i < n; i++)
	  {
	     cum[i] = cum[i] / tot;
	  }
	*total = tot;
     }
   return 0;
}

/*}}}*/
/*{{{ compute_cum_spectrum */

/* If totalp != NULL, spectrum will be normalized and normalization constant
 * will be returned via *totalp
 */
static double *compute_cum_spectrum (double *energies, double *spect,
				     unsigned int n, double *totalp)
{
   double *cum;

   if (NULL == (cum = JDMdouble_vector (n))) return NULL;
   if (-1 == compute_cumulated_vector (cum, energies, spect, n, totalp))
     {
	JDMfree_double_vector (cum);
	return NULL;
     }
   return cum;
}

/*}}}*/

/*{{{ FLAT spectrum functions */

static int read_flat_spectrum_parameters (Param_File_Type *p,
					  double *min_energy,
					  double *max_energy,
					  double *flux)
{
   if ((-1 == pf_get_double (p, "MinEnergy", min_energy))
       || (-1 == pf_get_double (p, "MaxEnergy", max_energy))
       || (-1 == pf_get_double (p, "SourceFlux", flux)))
     return -1;

   if (*flux <= 0.0)
     {
	marx_error ("SourceFlux for FLAT spectrum must be positive, not %e",
		    *flux);
	return -1;
     }

   return 0;
}

static int flat_energy (Marx_Spectrum_Type *p, double *e)
{
   double emin, de;

   emin = p->s.flat.emin;
   de = p->s.flat.emax - emin;

   *e = emin + de * JDMrandom ();
   return 0;
}

int _marx_init_flat_spectrum (Param_File_Type *pf, Marx_Spectrum_Type *st)
{
   double de;
   double min_energy, max_energy, flux;

   memset ((char *) st, 0, sizeof (Marx_Spectrum_Type));

   if (-1 == read_flat_spectrum_parameters (pf, &min_energy,
					    &max_energy, &flux))
     return -1;

   de = max_energy - min_energy;

   if ((de < 0.0) || (max_energy <= 0.0))
     return -1;

   st->type = MARX_FLAT_SPECTRUM;
   st->s.flat.emin = min_energy;
   st->s.flat.emax = max_energy;
   st->s.flat.flux = flux;
#if 0
   if (de == 0.0)
     st->total_flux = flux;
   else st->total_flux = flux * de;
#else
   st->total_flux = flux;
#endif
   st->energy_function = flat_energy;
   st->close_spectrum = NULL;
   return 0;
}

/*}}}*/
/*{{{ FILE spectrum functions */

static int file_energy (Marx_Spectrum_Type *p, double *e)
{
   marx_get_random_event (p->s.file.energies, p->s.file.cum_flux, p->s.file.num,
			  e, 1);
   return 0;
}

static int read_file_spectrum_parameters (Param_File_Type *p, char *file, unsigned int filelen, double *flux)
{
   if ((-1 == pf_get_file (p, "SpectrumFile", file, filelen))
       || (-1 == pf_get_double (p, "SourceFlux", flux)))
     return -1;

   return 0;
}

static void close_file_spectrum (Marx_Spectrum_Type *sp)
{
   if (sp->s.file.energies != NULL)
     {
	JDMfree_double_vector (sp->s.file.energies);
	sp->s.file.energies = NULL;
     }
   if (sp->s.file.cum_flux != NULL)
     {
	JDMfree_double_vector (sp->s.file.cum_flux);
	sp->s.file.cum_flux = NULL;
     }
}

int _marx_init_file_spectrum (Param_File_Type *p, Marx_Spectrum_Type *sp)
{
   int cindex[2];
   double *data[2];
   double *cum_flux, total_flux, flux;
   unsigned int nread;
   char file [PF_MAX_LINE_LEN];

   memset ((char *) sp, 0, sizeof (Marx_Spectrum_Type));

   if (-1 == read_file_spectrum_parameters (p, file, sizeof (file), &flux))
     return -1;

   cindex[0] = 1;		       /* energy (KeV) */
   cindex[1] = 1;		       /* flux (Photons/Kev/cm^2) */

   marx_message ("Reading spectrum from:\n\t%s\n", file);

   if ((-1 == JDMread_column_ddata (file, data, cindex, 2, &nread))
       || (nread == 0))
     {
	marx_error ("Error reading spectrum from %s", file);
	return -1;
     }

   /* Note: this returns normalized flux */

   cum_flux = compute_cum_spectrum (data[0], data[1], nread, &total_flux);

   JDMfree_double_vector (data[1]);
   if (NULL == cum_flux)
     {
	JDMfree_double_vector (data[0]);
	return -1;
     }

   if (flux <= 0.0)
     flux = total_flux;

   marx_message ("Your spectrum has a total flux of %e photons/cm^2/sec\n", flux);
   if ((flux <= 0.0) || (0 == JDMfinite (flux)))
     {
	marx_error ("The flux must be positive and finite.  Check your spectrum file");
	JDMfree_double_vector (data[0]);
	return -1;
     }

   sp->s.file.energies = data[0];
   sp->s.file.num = nread;
   sp->s.file.cum_flux = cum_flux;

   sp->total_flux = flux;

   sp->type = MARX_FILE_SPECTRUM;
   sp->energy_function = file_energy;
   sp->close_spectrum = close_file_spectrum;

   return 0;
}

/*}}}*/

int _marx_get_simple_specrum_parms (Param_File_Type *p, Marx_Source_Type *st, char *name) /*{{{*/
{
   char buf[PF_MAX_LINE_LEN];

   if (-1 == pf_get_string (p, "SpectrumType", buf, sizeof (buf)))
     return -1;

   if (!strcmp (buf, "FLAT"))
     {
	if (-1 == _marx_init_flat_spectrum (p, &st->spectrum))
	  return -1;
     }
   else if (!strcmp (buf, "FILE"))
     {
	if (-1 == _marx_init_file_spectrum (p, &st->spectrum))
	  return -1;
     }
   else
     {
	marx_error ("SpectrumType: \"%s\" not supported by \"%s\" source",
		    buf, name);
	return -1;
     }
   return 0;
}

/*}}}*/

#if 0
static FILE *open_marx_rayfile (char *file) /*{{{*/
{
   FILE *fp;
   unsigned long magic;

   if (NULL == (fp = fopen (file, "rb")))
     return NULL;

   if ((1 != fread (&magic, sizeof (magic), 1, fp))
       || (magic != RAYFILE_MAGIC_LONG))
     {
	fclose (fp);
	marx_error ("Magic number mismatch: %s", file);
	return NULL;
     }

   return fp;
}

/*}}}*/

int marx_spectrum_from_ray_file (Spectrum_Type *s, char *file) /*{{{*/
{
   s->type = MARX_SPECTRUM_RAY;

   if (file == NULL) return -1;

   if (NULL == (s->fp = open_marx_rayfile (file)))
     return -1;

   return 0;
}

/*}}}*/

#endif
