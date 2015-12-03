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
/* interpreter interface */

#include "config.h"
#include "marx-feat.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <string.h>

#include <time.h>

#include <jdmath.h>
#include <slang.h>

#include "marx.h"

#ifdef USE_SLWDG
# include "slwdg.h"
#endif

#ifdef USE_GNUPLOT
# include "gnuplot.h"
#endif

#ifdef SIM_LIB_DIR
static char Sim_Lib_Dir[256] = SIM_LIB_DIR;
#else
static char Sim_Lib_Dir[256];
#endif

double Sim_Total_Time;

static void sim_print (char *txt)
{
   fputs (txt, stderr);
}

typedef struct
{
   char name[16];
   Sim_Name_Type *tbl;
} Variable_Table_Type;

#define MAX_TABLES 16
Variable_Table_Type Var_Tables[MAX_TABLES];

int register_variable_table (char *name, Sim_Name_Type *tbl)
{
   unsigned int i;
   for (i = 0; i < MAX_TABLES; i++)
     {
	if ((*Var_Tables[i].name == 0)
	    || !strcmp (Var_Tables[i].name, name))
	  {
	     Var_Tables[i].tbl = tbl;
	     strncpy (Var_Tables[i].name, name, 15);
	     Var_Tables[i].name[15] = 0;
	     return 0;
	  }
     }
   marx_error ("Max variable tables exceeded.");
   return -1;
}

static Sim_Name_Type *do_module_variable (void)
{
   char *name;
   char *tbl;
   int free_name = 0, free_tbl = 0;
   unsigned int i;
   Sim_Name_Type *t;

   if (SLang_pop_string (&name, &free_name)
       || SLang_pop_string (&tbl, &free_tbl))
     {
	t = NULL;
	goto cleanup;
     }

   for (i = 0; i < MAX_TABLES; i++)
     {
	if (*Var_Tables[i].name == 0) break;
	if (strcmp (Var_Tables[i].name, tbl)) continue;
	t = Var_Tables[i].tbl;

	while (t->name != NULL)
	  {
	     if (!strcmp (t->name, name))
	       {
		  goto cleanup;
	       }
	     t++;
	  }
	marx_error ("Unknown name %s for module %s.", name, tbl);
	break;
     }
   t = NULL;
   marx_error ("Unknown module: %s.", tbl);

   cleanup:
   if (free_tbl) SLFREE (tbl);
   if (free_name) SLFREE (name);

   return t;
}

static void set_module_variable (void)
{
   double val;
   int unused;
   Sim_Name_Type *t;

   if (SLang_pop_float (&val, &unused, &unused)) return;
   t = do_module_variable ();
   if (t != NULL) *t->ptr = val;
}

static double get_module_variable (void)
{
   Sim_Name_Type *t;

   t = do_module_variable ();
   if (t != NULL) return *t->ptr;
   return 0.0;
}

static SLang_Name_Type *Spectrum_Function;

static double user_spectrum_function (double en)
{
   double flux;
   int dum1, dum2;

   SLang_push_float (en);
   if (SLang_Error) goto error;

   SLexecute_function (Spectrum_Function);

   if (SLang_pop_float (&flux, &dum1, &dum2)) goto error;

   return flux;

   error:
   marx_error ("Error in user defined spectrum.");
   return 0.0;
}

static void user_spectrum (void)
{
   double emin, emax;
   int npts;
   char *function_name;
   int unused, do_free = 0;

   if (SLang_pop_string (&function_name, &do_free)
       || SLang_pop_integer (&npts)
       || SLang_pop_float (&emax, &unused, &unused)
       || SLang_pop_float (&emin, &unused, &unused))
     {
	goto clean_up;
     }

   if (npts <= 0)
     {
	marx_error ("Number of points must be grater than zero.");
	goto clean_up;
     }

   Spectrum_Function = SLang_get_function (function_name);
   if (Spectrum_Function == NULL)
     {
	marx_error ("Function %s is undefined.", function_name);
	goto clean_up;
     }

   Marx_Spectrum.n_energies = (unsigned int) npts;
   Marx_Spectrum.min_energy = emin;
   Marx_Spectrum.max_energy = emax;

   if (-1 == generic_spectrum (&Marx_Spectrum, user_spectrum_function))
     {
	marx_error ("Unable to initialize user defined spectrum (%s)", function_name);
	goto clean_up;
     }

   clean_up:
   if (do_free) SLFREE (function_name);
}

static void datafile_spectrum (char *file)
{
   if (-1 == spectrum_from_file (&Marx_Spectrum, file))
     {
	marx_error ("Unable to initialize spectrum from file %s", file);
     }
}

static void do_acis_init_detector (char *file)
{
#if 0
   if (-1 == acis_init_detector (file))
     {
	marx_error ("Unable to get detector efficiencies from %s.", file);
     }
#endif
}

static void do_hrc_init_detector (char *file1, char *file2)
{
#if 0
   if (-1 == hrc_init_detector (file1, file2))
     {
	marx_error ("Unable to get detector efficiencies from %s and %s",
		   file1, file2);
     }
#endif
}

static void do_flat_spectrum_init (void)
{
   double emin, emax, flux;
   int unused;

   if (SLang_pop_float (&flux, &unused, &unused)
       || SLang_pop_float (&emax, &unused, &unused)
       || SLang_pop_float (&emin, &unused, &unused))
     return;

   Marx_Spectrum.min_energy = emin;
   Marx_Spectrum.max_energy = emax;
   if (-1 == flat_spectrum (&Marx_Spectrum, flux))
     {
	marx_error ("Unable to initialize flat spectrum.");
     }
}

static void do_mirror_init (char *file)
{
#if 0
   if (-1 == mirror_init (file))
     {
	marx_error ("Unable to get mirror areas from %s.", file);
     }
#endif
}

static void do_grating_init (void)
{
#if 0
   if (-1 == grating_init ())
     {
	marx_error ("Unable to initialize grating.");
     }
#endif
}

static void do_mirror_blur_init (void)
{
#if 0
   if (-1 == init_mirror_blur ())
     {
	marx_error ("Unable to initial mirror blur.");
     }
#endif
}

static void do_collect_photons (int *num_photons)
{
   if (*num_photons <= 0)
     {
	marx_error ("collect_photons: number of photons must be positive.");
	return;
     }

   if (-1 == collect_photons ((unsigned int) *num_photons, &Marx_Photon_List, &Marx_Spectrum))
     {
	marx_error ("Error collecting photons.");
     }
   Sim_Total_Time += Marx_Photon_List.total_time;
}

static void append_string_to_file (char *file, char *str)
{
   FILE *fp;
   if (NULL == (fp = fopen (file, "a")))
     {
	marx_error ("Unable to open %s", file);
	return;
     }
   (void) fputs (str, fp);
   fclose (fp);
}

static void do_regrid_curves (char *file)
{
   unsigned int i, imax, col;
   double energy;
   Spectrum_Type *st = &Marx_Spectrum;
   FILE *fp;

   if ((st->cum_flux == NULL)
       || (st->energies == NULL))
     {
	marx_error ("regrid_curves: an energy grid has not been defined.");
	return;
     }
   imax = st->n_energies;

   if (NULL == (fp = fopen (file, "w")))
     {
	marx_error ("regrid_curves: unable to open %s", file);
	return;
     }

   col = 1;
   fprintf (fp, "#Column %d: energy in KeV", col);
   col = dump_mirror_curves_open (col, fp);
   col = dump_acis_curves_open (col, fp);
   col = dump_grating_curves_open (col, fp);
   fputc ('\n', fp);

   for (i = 0; i < imax; i++)
     {
	energy = st->energies[i];
	fprintf (fp, "%e", energy);

	if (
	    (-1 == dump_mirror_curves (energy, fp))
	    || (-1 == dump_acis_curves (energy, fp))
	    || (-1 == dump_grating_curves (energy, fp))
	    )
	  {
	     marx_error ("regrid_curves: unknown regrid error.");
	     break;
	  }

	if (EOF == fputc ('\n', fp))
	  {
	     marx_error ("regrid_curves: write error.");
	     break;
	  }
     }

   if (EOF == fclose (fp))
     marx_error ("regrid_curves: Error closing file.");
}

static char *get_time (void)
{
   char *the_time;
   time_t myclock;

   myclock = time((time_t *) 0);
   the_time = (char *) ctime(&myclock);
   /* returns the form Sun Sep 16 01:03:52 1985\n\0 */
   the_time[24] = '\0';
   return(the_time);
}

static void dump_to_log_file (char *file)
{
   FILE *fp;
   if (NULL == (fp = fopen (file, "a")))
     {
	fprintf (stderr, "Unable to create logfile %s", file);
	return;
     }
   fprintf (fp, "\n****%s: LOG FILE ENTRY CREATED by user %s\n",
	    get_time (),
	    (getenv ("USER") == NULL ? "unknown" : getenv("USER")));

   (void) dump_mirror_variables (fp);
   (void) dump_grating_variables (fp);
   (void) dump_acis_variables (fp);
#if 1
   (void) dump_hrc_variables (fp);
   /* (void) dump_disperse_variables (fp); */
#endif

   (void) fclose (fp);
}

static void do_mirror_reflect (void)
{
   if (-1 == mirror_reflect (&Marx_Photon_List))
     {
	marx_error ("Error during reflection.");
     }
}

static void do_mirror_perfect_reflect (void)
{
   if (-1 == mirror_perfect_reflect (&Marx_Photon_List))
     {
	marx_error ("Error during perfect reflection.");
     }
}

static void do_diffract (void)
{
   if (-1 == diffract (&Marx_Photon_List))
     {
	marx_error ("Unable to diffract photons.");
     }
}

static void do_mirror_blur (void)
{
   if (-1 == mirror_blur (&Marx_Photon_List))
     {
	marx_error ("Unable to blur photons.");
     }
}

static void do_disperse_photons (void)
{
#if 0
   if (-1 == disperse_photons (&Marx_Photon_List))
     {
	marx_error ("Unable to disperse photons.");
     }
#endif
}

static void do_acis_detect (void)
{
   if (-1 == acis_ccd_detect (&Marx_Photon_List))
     {
	marx_error ("Unable to detect photons.");
     }
}

static void do_hrc_detect (void)
{
#if 1
   if (-1 == hrc_detect (&Marx_Photon_List))
     {
	marx_error ("Unable to detect photons.");
     }
#endif
}

static int do_photon_number (void)
{
   prune_photons (&Marx_Photon_List);
   return Marx_Photon_List.num_sorted;
}

static void do_quit (void)
{
#ifdef USE_SLWDG
   /* Just in case the user forgot to do this */
   SLwdg_end ();
#endif
   exit (0);
}

static int write_formatted (FILE *fp, char *fmt)
{
   unsigned int i, imax;
   Photon_Attr_Type *at, *attr;
   double start_time;
   char output_buf[4096], *f, *p;
   unsigned int count;

   start_time = Sim_Total_Time - Marx_Photon_List.total_time;

   imax = Marx_Photon_List.n_photons;

   attr = Marx_Photon_List.attributes;
   count = 0;
   if (attr != NULL) for (i = 0; i < imax; i++)
     {
	char ch;

	at = attr + i;
	if (at->flags & BAD_PHOTON_MASK) continue;

	if (fmt == NULL)
	  {
	     fprintf (fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%u\t%u\t%u\t%u\t%d\t%f\n",
		      at->energy, at->arrival_time + start_time,
		      at->x.x, at->x.y, at->x.z,
		      at->p.x, at->p.y, at->p.z,
		      at->ccd_num, at->y_pixel, at->z_pixel,
		      at->mirror_shell, at->order, at->pulse_height
		      );
	     continue;
	  }

	f = fmt;
	count++;
	p = output_buf;

	while ((ch = *f++) != 0)
	  {
	     if (ch != '%')
	       {
		  *p++ = ch;
		  continue;
	       }

	     ch = *f++;

	     switch (ch)
	       {
		case 'e':
		  sprintf (p, "%f", at->energy);
		  break;
		case 't':
		  sprintf (p, "%f", at->arrival_time + start_time);
		  break;
		case 'x':
		  sprintf (p, "%f", at->x.x);
		  break;
		case 'y':
		  sprintf (p, "%f", at->x.y);
		  break;
		case 'z':
		  sprintf (p, "%f", at->x.z);
		  break;
		case 'a':
		  sprintf (p, "%f", at->p.x);
		  break;
		case 'b':
		  sprintf (p, "%f", at->p.y);
		  break;
		case 'g':
		  sprintf (p, "%f", at->p.z);
		  break;
		case 'c':
		  sprintf (p, "%u", at->ccd_num);
		  break;
		case 'Y':
		  sprintf (p, "%u", at->y_pixel);
		  break;
		case 'Z':
		  sprintf (p, "%u", at->z_pixel);
		  break;
		case 'p':
		  sprintf (p, "%f", at->pulse_height);
		  break;
		case 'm':
		  sprintf (p, "%u", at->mirror_shell);
		  break;
		case 'o':
		  sprintf (p, "%d", at->order);
		  break;
		case '%':
		  *p++ = '%';
		  *p = 0;
		  break;
		case 'd':
		  sprintf (p, "%u", count);
		  break;
		default:
		  marx_error ("Unknown format specifier: %c.", ch);
		  return -1;
	       }

	     p += strlen (p);
	  }
	*p = '\n';
	*(p + 1) = 0;
	fputs (output_buf, fp);
     }
   return 0;
}

static void write_formatted_photon_list (char *fmt, char *file, int *append)
{
   FILE *fp;
   char *mode;

   if (*append == 0) mode = "wt";
   else mode = "at";

   fp = fopen (file, mode);
   if (fp == NULL)
     {
	marx_error ("Unable to open output file %s.", file);
	return;
     }
   write_formatted (fp, fmt);
   fclose (fp);
}

void marx_write_photon_list (char *file, int *append)
{
   write_formatted_photon_list (NULL, file, append);
}

static void interp_marx_error (char *s)
{
   marx_error (s);
}

static char Simulator_Version[40];
static void make_simulator_version_string (void)
{
   int a, b, c;
   a = MARX_VERSION / 10000;
   b = MARX_VERSION % 10000;
   c = b % 100;
   b = b / 100;

   sprintf (Simulator_Version, "%d.%02d.%02d", a, b, c);
}

static char *make_date (int what)
{
   static char date[20];
   struct tm *tm_struct;
   time_t tloc;

   time (&tloc);
   tm_struct = gmtime (&tloc);
   if (what)
     {
	sprintf (date, "%02d/%02d/%02d",
		 tm_struct->tm_mday, 1 + tm_struct->tm_mon, tm_struct->tm_year);
     }
   else sprintf (date, "%02d:%02d:%02d",
		 tm_struct->tm_hour,
		 tm_struct->tm_min, tm_struct->tm_sec);
   return date;
}

static char *ut_date (void)
{
   return make_date (1);
}
static char *ut_time (void)
{
   return make_date (0);
}

#ifdef USE_GNUPLOT
/* I am going to ignore all errors in the plotting routines. */
static void open_gnuplot (int *idp)
{
   unsigned int id = (unsigned int) *idp;
   if (-1 == gnuplot_open (id, NULL))
     return;
}

static void close_gnuplot (int *idp)
{
   unsigned int id = (unsigned int) *idp;
   if (-1 == gnuplot_close (id))
     return;
}

static void cmd_gnuplot (int *idp, char *cmd)
{
   unsigned int id = (unsigned int) *idp;
   if (-1 == gnuplot_cmd (id, cmd))
     return;
}

static void plot_data (int *idp, char *fmt)
{
   char fmtbuf[10];
   unsigned int id = (unsigned int) *idp;
   FILE *fp;

   fp = gnuplot_get_fp (id);
   if (fp == NULL) return;

   if (-1 == gnuplot_cmd (id, "plot '-'"))
     return;

   if ((fmt[0] == 0)
       || (fmt[1] == 0))
     fmt = "xy";

   sprintf (fmtbuf, "%%%c %%%c", fmt[0], fmt[1]);

   write_formatted (fp, fmtbuf);
   gnuplot_cmd (id, "e");
}

#endif 				       /* USE GNUPLOT */

static SLang_Name_Type Sim_Intrinsics[] =
{
   MAKE_INTRINSIC(".error", interp_marx_error, VOID_TYPE, 1),
   MAKE_INTRINSIC(".print", sim_print, VOID_TYPE, 1),
   MAKE_INTRINSIC(".collect_photons", do_collect_photons, VOID_TYPE, 1),
   MAKE_INTRINSIC(".mirror_reflect", do_mirror_reflect, VOID_TYPE, 0),
   MAKE_INTRINSIC(".mirror_perfect_reflect", do_mirror_perfect_reflect, VOID_TYPE, 0),
   MAKE_INTRINSIC(".diffract", do_diffract, VOID_TYPE, 0),
   MAKE_INTRINSIC(".mirror_blur", do_mirror_blur, VOID_TYPE, 0),
   MAKE_INTRINSIC(".disperse_photons", do_disperse_photons, VOID_TYPE, 0),
   MAKE_INTRINSIC(".acis_detect", do_acis_detect, VOID_TYPE, 0),
   MAKE_INTRINSIC(".hrc_detect", do_hrc_detect, VOID_TYPE, 0),
   MAKE_INTRINSIC(".photon_number", do_photon_number, INT_TYPE, 0),
   MAKE_INTRINSIC(".quit", do_quit, VOID_TYPE, 0),
   MAKE_INTRINSIC(".acis_detector_init", do_acis_init_detector, VOID_TYPE, 1),
   MAKE_INTRINSIC(".hrc_detector_init", do_hrc_init_detector, VOID_TYPE, 2),
   MAKE_INTRINSIC(".flat_spectrum_init", do_flat_spectrum_init, VOID_TYPE, 0),
   MAKE_INTRINSIC(".mirror_init", do_mirror_init, VOID_TYPE, 1),
   MAKE_INTRINSIC(".grating_init", do_grating_init, VOID_TYPE, 0),
   MAKE_INTRINSIC(".mirror_blur_init", do_mirror_blur_init, VOID_TYPE, 0),
   MAKE_INTRINSIC(".write_photon_list", marx_write_photon_list, VOID_TYPE, 2),
   MAKE_INTRINSIC(".write_formatted_photon_list", write_formatted_photon_list, VOID_TYPE, 3),
   MAKE_INTRINSIC(".set_data_directory", set_data_directory, VOID_TYPE, 1),
   MAKE_INTRINSIC(".set_module_variable", set_module_variable, VOID_TYPE, 0),
   MAKE_INTRINSIC(".get_module_variable", get_module_variable, FLOAT_TYPE, 0),
   MAKE_INTRINSIC(".user_spectrum", user_spectrum, VOID_TYPE, 0),
   MAKE_INTRINSIC(".data_file_spectrum", datafile_spectrum, VOID_TYPE, 1),
   MAKE_INTRINSIC(".ut_date", ut_date, STRING_TYPE, 0),
   MAKE_INTRINSIC(".ut_time", ut_time, STRING_TYPE, 0),
   MAKE_INTRINSIC(".regrid_curves", do_regrid_curves, VOID_TYPE, 1),
   MAKE_INTRINSIC(".dump_to_log_file", dump_to_log_file, VOID_TYPE, 1),
   MAKE_INTRINSIC(".append_string_to_file", append_string_to_file, VOID_TYPE, 2),
#ifdef USE_GNUPLOT
   MAKE_INTRINSIC (".gnuplot_open", open_gnuplot, VOID_TYPE, 1),
   MAKE_INTRINSIC (".gnuplot_close", close_gnuplot, VOID_TYPE, 1),
   MAKE_INTRINSIC (".gnuplot_cmd", cmd_gnuplot, VOID_TYPE, 2),
   MAKE_INTRINSIC (".gnuplot_plot_data", plot_data, VOID_TYPE, 2),
#endif
   MAKE_VARIABLE(".USE_LETG", &Sim_Use_LETG, INT_TYPE, 0),
   MAKE_VARIABLE(".TOTAL_TIME", &Sim_Total_Time, FLOAT_TYPE, 0),
   MAKE_VARIABLE(".X_HRMA", &HRMA_Focal_Length, FLOAT_TYPE, 0),
   MAKE_VARIABLE(".X_GRATING", &Rowland_Torus_Diameter, FLOAT_TYPE, 0),
   MAKE_VARIABLE(".X_DETECTOR", &Focus_Detector_Offset, FLOAT_TYPE, 0),
   MAKE_VARIABLE(".GRATING_VIGNETTING_FACTOR", &Grating_Vignetting_Factor, FLOAT_TYPE, 0),
   MAKE_VARIABLE(".MIRROR_VIGNETTING_FACTOR", &Mirror_Vignetting_Factor, FLOAT_TYPE, 0),
   MAKE_VARIABLE(".SIMULATOR_VERSION", Simulator_Version, STRING_TYPE, 1),
   MAKE_VARIABLE(".SIM_LIB_DIR", Sim_Lib_Dir, STRING_TYPE, 1),

   SLANG_END_TABLE
};

static void run_startup_file (void)
{
   char file[256], *sld;
   unsigned int len = sizeof (Sim_Lib_Dir) - 3;

   sld = getenv ("SIM_LIB_DIR");

   if (sld != NULL)
     {
	strncpy (Sim_Lib_Dir, sld, len);
	Sim_Lib_Dir[len] = 0;
     }

   len = strlen (Sim_Lib_Dir);

   if (len)
     {
	struct stat st;

	if (Sim_Lib_Dir[len - 1] != '/')
	  {
	     Sim_Lib_Dir[len] = '/';
	     Sim_Lib_Dir[len + 1] = 0;
	  }

	sprintf (file, "%sstartup.sl", Sim_Lib_Dir);
	if (0 == stat (file, &st))
	  {
	     SLang_load_file (file);
	     if (SLang_Error) exit (-1);
	     return;
	  }
	fprintf (stderr, "Unable to read startup file %s.\n", file);
     }
   else fprintf (stderr, "\
Unable to read startup file because the simulator library directory has not\n\
been specified.\n"
		 );

   fprintf (stderr, "\
You may need to set the environment variable SIM_LIB_DIR to point to the\n\
directory where the simulator startup file 'startup.sl' can be found.\n"
	    );
   exit (-1);
}

int marx_interp_main (char *file)
{
   /* Initialize the library.  This is always needed. */

   if (!init_SLang()		       /* basic interpreter functions */
       || !init_SLmath() 	       /* sin, cos, etc... */
#ifdef unix
       || !init_SLunix()	       /* unix system calls */
#endif
       || !init_SLfiles()	       /* file i/o */
#ifdef USE_SLWDG
       || !init_SLwdg ()
#endif
       /* Now add intrinsics for this application */
       || !SLang_add_table(Sim_Intrinsics, "Sim")
       || (-1 == sim_init_fits ())
#ifdef USE_GNUPLOT
       || (0 == SLdefine_for_ifdef ("USE_GNUPLOT"))
#endif
       )
     {
	fprintf(stderr, "Unable to initialize Interpreter\n");
	exit(-1);
     }

   if ((-1 == init_mirror_module ())
       || (-1 == init_grating_module ())
       || (-1 == init_acis_module ()))
     {
	fprintf (stderr, "Unable to initialize module symbols.\n");
	exit (-1);
     }

   make_simulator_version_string ();

#ifdef USE_SLWDG
   SLwdg_initialize (0);
#endif
   /* Turn on debugging */
   SLang_Traceback = 1;

   run_startup_file ();

   /* Now load an initialization file and exit */
   SLang_load_file (file);

#ifdef USE_SLWDG
   SLwdg_end ();
#endif
   return (SLang_Error);
}
