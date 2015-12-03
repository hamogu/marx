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

/*{{{ Include Files */

#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

/* These are needed for mkdir */
#include <sys/types.h>
#include <sys/stat.h>

#include <errno.h>
#include <string.h>

#include <math.h>
#include <jdfits.h>

#ifdef __linux__
# ifdef __GLIBC__
#  include <fpu_control.h>
# else
#  include <i386/fpu_control.h>
# endif
# ifndef _FPU_SETCW
#  define _FPU_SETCW __setfpucw
# endif
#endif

#include "marx.h"

/*}}}*/

/*{{{ Static Variables */

#define PRINT_STATS_ARRAY 0

static char *Parameter_File;
static int Mirror_Module;
static int Grating_Module;
static int Detector_Module;

/* These are initialized from the marx.par file */
static double Exposure_Time = -0.0;
static char *Data_To_Write;
static char *Data_Directory;
/* See the note below where Output_Dir is initialized. */
static char *Output_Dir;
static char *Source_Name;
static double DetOffset_X;
static double DetOffset_Y;
static double DetOffset_Z;

static int Num_Rays_Per_Iteration;
static int Num_Rays_Marx_Par;

static int Random_Seed = -1;
static int Dump_To_Rayfile;
static char *Rayfile_Name;

static double MJDref = 50814.0;
static double MJDref_Years = 1998.0;
static time_t TStart_Unix;
static double TStart_MJDref;

/*}}}*/

/*{{{ Forward Function Declarations */

static Param_File_Type *setup_parms (int argc, char **argv);
static unsigned long compute_write_mask (void);
/* static int cp_par_file (char *, char *); */
static void version (void);
static int module_init (Param_File_Type *);
static int write_fits_info (double);

static int Is_Output_Piped;
static int write_photons_to_pipe (FILE *, unsigned int, Marx_Photon_Type *, double, double);
static FILE *open_pipe (char *);
static int close_pipe (FILE *);

/*}}}*/

static void message_copyright (FILE *fp)
{
   if (Marx_Verbose == 0)
     return;

   fprintf (fp, "MARX version %s, Copyright (C) 2002-2018 Massachusetts Institute of Technology\n\n",
	    MARX_VERSION_STRING);
}

static void usage (void) /*{{{*/
{
   message_copyright (stderr);

   pf_usage ("marx.par", 0);

   fprintf (stderr, "\
marx usage forms:\n\
   marx [parameters]\n\
   marx --dump file...\n\
   marx --raydump file\n\
   marx --version\n\
   marx --help\n\
");
#ifdef MARX_PFILE_DIR
   fprintf (stderr, "\nmarx parameter files may be found in:\n %s/\n", MARX_PFILE_DIR);
#endif

   exit (1);
}

/*}}}*/

static int Needs_Cp_Par_File = 0;
static int open_output_dir (Param_File_Type *pf)
{
   (void) pf;

   Needs_Cp_Par_File = 1;
   return 0;
}

static int write_to_output_dir (Marx_Photon_Type *pt, unsigned int num_collected,
				int open_mode, unsigned long write_mask, double tstart, double duration)
{
   (void) num_collected;
   (void) duration;

   marx_message ("\nWriting output to directory '%s' ...",
		 Output_Dir);

   return marx_write_photons (pt,
			      write_mask,
			      Output_Dir,
			      open_mode,
			      tstart);
}

static int close_output_dir (double total_time)
{
   return write_fits_info (total_time);
}

static int open_rayfile (Param_File_Type *p)
{
   return open_output_dir (p);
}

static int write_to_rayfile (Marx_Photon_Type *pt, unsigned int num_collected, int open_mode, unsigned long write_mask, double tstart, double duration)
{
   (void) num_collected;
   (void) write_mask;
   (void) duration;

   marx_message ("\nDumping to RayFile %s...", Rayfile_Name);

   if (-1 == marx_dump_to_rayfile (Rayfile_Name, open_mode, pt, tstart))
     {
	marx_error ("Error writing to rayfile.");
	return -1;
     }
   return 0;
}

static int close_rayfile (double total_time)
{
   return close_output_dir (total_time);
}

static FILE *Pipe_Fp;
static int open_output_pipe (Param_File_Type *p)
{
   (void) p;
   if (NULL == (Pipe_Fp = open_pipe (Output_Dir)))
     return -1;

   return 0;
}

static int write_to_pipe (Marx_Photon_Type *pt, unsigned int num_collected, int open_mode, unsigned long write_mask, double tstart, double duration)
{
   (void) open_mode;
   (void) write_mask;

   marx_message ("\nWriting %u photons to pipe...", num_collected);
   if (-1 == write_photons_to_pipe (Pipe_Fp, num_collected, pt, tstart, duration))
     {
	marx_error ("Writing to pipe failed");
	return -1;
     }
   return 0;
}

static int close_output_pipe (double total_time)
{
   int status = close_pipe (Pipe_Fp);
   (void) total_time;
   Pipe_Fp = NULL;
   return status;
}

#if PRINT_STATS_ARRAY
static unsigned long Stats[4];
#endif

static int process_photons (Marx_Photon_Type *pt)
{
#if PRINT_STATS_ARRAY
   marx_prune_photons (pt); Stats[0] += pt->num_sorted;
#endif
   if (-1 == marx_mirror_reflect (pt, 1))
     {
	marx_error ("Error during reflection from mirror");
	return -1;
     }

#if PRINT_STATS_ARRAY
   marx_prune_photons (pt); Stats[1] += pt->num_sorted;
#endif
   if (-1 == marx_grating_diffract (pt, 1))
     {
	marx_error ("Unable to diffract photons.");
	return -1;
     }

#if PRINT_STATS_ARRAY
   marx_prune_photons (pt); Stats[2] += pt->num_sorted;
#endif
   if (-1 == marx_detect (pt, 1))
     {
	marx_error ("Error during detection.");
	return -1;
     }

#if PRINT_STATS_ARRAY
   marx_prune_photons (pt); Stats[3] += pt->num_sorted;
#endif
   return 0;
}

static time_t marx_timegm (struct tm *tms)
{
#ifdef HAVE_TIMEGM
   return timegm (tms);
#else
   /* This is a hack.  It should be fixed. */
   return mktime (tms);
#endif
}

static int setup_tstart (Param_File_Type *pf)
{
   double secs_per_year = 365.25 * 24.0 * 3600.0;
   double tstart, yrs_tstart;
   time_t time_t_mjdref;
   struct tm tm;
#if 0
   char *asol_file = NULL;
   char buf[2048];
   if (-1 == pf_get_string (pf, "SourceType", buf, sizeof(buf)))
     return -1;
   if (0 == jdfits_strcasecmp (buf, "SAOSAC"))
     {
	if (-1 == pf_get_string (pf, "DitherModel", buf, sizeof(buf)))
	  return -1;
	if (0 == jdfits_strcasecmp (buf, "FILE"))
	  {
	     if (-1 == pf_get_string (pf, "DitherFile", buf, sizeof(buf)))
	       return -1;
	     asol_file = buf;
	  }
     }
#endif

   if (-1 == pf_get_double (pf, "TStart", &tstart))
     return -1;

   if (tstart < 1999)
     {
	marx_error ("A tstart value less than 1999 is not permitted");
	return -1;
     }

   if (tstart < 2100)
     {
	/* asssume that tstart is in years */
	yrs_tstart = tstart;
	tstart = (yrs_tstart - MJDref_Years)*secs_per_year;
     }
   else
     {
	/* assume tstart is secs from mjdref */
	yrs_tstart = MJDref_Years + tstart/secs_per_year;
     }
#if 0
   if (asol_file != NULL)
     {
	JDFits_Type *ft;
	double t0, t1;

	marx_message ("For dithered SAOSAC rays, the TSTART value will be taken from %s\n", asol_file);
	if (NULL == (ft = jdfits_open_binary_table (asol_file, "ASPSOL")))
	  return -1;

	if ((-1 == jdfits_read_keyword_dbl (ft, "TSTART", &t0))
	    || (-1 == jdfits_read_keyword_dbl (ft, "TSTOP", &t1)))
	  {
	     (void) jdfits_close_file (ft);
	     return -1;
	  }
	(void) jdfits_close_file (ft);
	if ((tstart < t0) || (tstart >= t1))
	  {
	     marx_error ("\
*** ERROR: For dithered SAOSAC rays, the marx.par TStart value must be in\n\
           the range %.17g <= TStart < %.17g.\n\
           The value in the marx.par file lies outside this range.\n\
",
			 t0, t1);
	  }
	return -1;
     }
#endif

   memset ((char *) &tm, 0, sizeof (tm));
   tm.tm_sec = 0;
   tm.tm_min = 0;
   tm.tm_hour = 0;
   tm.tm_mday = 1;
   tm.tm_mon = 0;
   tm.tm_year = 98;
   tm.tm_yday = 0;
   tm.tm_isdst = 0;
   time_t_mjdref = marx_timegm (&tm);

   TStart_Unix = (time_t) (time_t_mjdref + tstart);
   TStart_MJDref = TStart_Unix - time_t_mjdref;

   if (-1 == marx_set_time (yrs_tstart, tstart))
     return -1;

   return 0;
}

int main (int argc, char **argv) /*{{{*/
{
   Param_File_Type *p;
   Marx_Source_Type *st;
   Marx_Photon_Type *pt;
   double total_time;
   unsigned long total_num_collected, total_num_detected;
   unsigned long num_to_collect, num_to_detect;
   unsigned long write_mask;
   int open_mode;
   int (*open_function)(Param_File_Type *);
   int (*write_function)(Marx_Photon_Type *, unsigned int, int, unsigned long, double, double);
   int (*close_function)(double);
   int status = -1;

   JDMATH_INIT;

#ifdef __linux__
# if 0
     {
	unsigned int cw = 0x1372;
	_FPU_SETCW (cw);
     }
# endif
#endif

   if ((argc > 1) && (argv[1][0] == '-'))
     {
	if (!strcmp (argv[1], "--dump"))
	  {
	     return marx_dump (argc, argv);
	  }
#if 1
	if ((argc == 3) && !strcmp (argv[1], "--raydump"))
	  {
	     return marx_dump_rayfile (argv[2]);
	  }
#endif
	if (!strcmp (argv[1], "--version"))
	  version ();

	usage ();
     }

   p = setup_parms (argc, argv);

   message_copyright (stdout);

   if (-1 == setup_tstart (p))
     return 1;

   write_mask = compute_write_mask ();

   if (-1 == module_init (p))
     {
	return 1;
     }

   if ((NULL == (st = marx_create_source (p)))
       || (-1 == marx_open_source (st)))
     {
	marx_error ("marx: Error creating source.");
	return 1;
     }

   pt = marx_alloc_photon_type (Num_Rays_Per_Iteration);
   if (pt == NULL)
     {
	marx_error ("marx: Error allocating photon array.");
	return 1;
     }

   if (Dump_To_Rayfile)
     {
	open_function = open_rayfile;
	write_function = write_to_rayfile;
	close_function = close_rayfile;
     }
   else if (Is_Output_Piped)
     {
	open_function = open_output_pipe;
	write_function = write_to_pipe;
	close_function = close_output_pipe;
     }
   else
     {
	open_function = open_output_dir;
	write_function = write_to_output_dir;
	close_function = close_output_dir;
     }

   if (-1 == (*open_function) (p))
     {
	pf_close_parameter_file (p);
	marx_dealloc_photon_type (pt);
	return 1;
     }

   if (Needs_Cp_Par_File)
     {
	char *marx_par;

	if (NULL == (marx_par = marx_dircat (Output_Dir, "marx.par")))
	  {
	     pf_close_parameter_file (p);
	     marx_dealloc_photon_type (pt);
	     return 1;
	  }
	if (-1 == pf_set_output_filename (p, marx_par))
	  {
	     marx_free (marx_par);
	     pf_close_parameter_file (p);
	     marx_dealloc_photon_type (pt);
	     return 1;
	  }
	marx_free (marx_par);
     }
   pf_close_parameter_file (p);

#if 0
   if (Needs_Cp_Par_File
       && (-1 == cp_par_file (Parameter_File, Output_Dir)))
     {
	marx_error ("Unable to copy paramater file to %s", Output_Dir);
	marx_dealloc_photon_type (pt);
	return 1;
     }
#endif
   marx_message ("System initialized.\n\n");

   /*
    * Main loop
    */

   total_num_collected = 0;
   total_num_detected = 0;
   open_mode = 1;
   total_time = 0.0;
   num_to_collect = 0;
   num_to_detect = 0;

   if (Exposure_Time > 0.0)
     marx_message ("Starting simulation.  Exposure Time set to %e seconds\n", Exposure_Time);
   else if (Num_Rays_Marx_Par < 0)
     {
	num_to_detect = (unsigned long) (-Num_Rays_Marx_Par);
	if ((unsigned long)Num_Rays_Per_Iteration > num_to_detect)
	  Num_Rays_Per_Iteration = num_to_detect;
	marx_message ("Starting simulation.  NumRays to detect = %lu, dNumRays = %u\n",
		      num_to_detect, Num_Rays_Per_Iteration);
     }
   else
     {
	num_to_collect = (unsigned long) Num_Rays_Marx_Par;
	if ((unsigned long)Num_Rays_Per_Iteration > num_to_collect)
	  Num_Rays_Per_Iteration = num_to_collect;
	marx_message ("Starting simulation.  NumRays to collect = %lu, dNumRays = %u\n",
		      num_to_collect, Num_Rays_Per_Iteration);
     }

#if PRINT_STATS_ARRAY
   Stats[0] = 0;
   Stats[1] = 0;
   Stats[2] = 0;
   Stats[3] = 0;
#endif
   while (1)
     {
	unsigned int num_collected;
	double *exposure_time_ptr = NULL;
	double exposure_time_left;

	if (Exposure_Time > 0.0)
	  {
	     exposure_time_left = (Exposure_Time - total_time);
	     if (exposure_time_left <= 0.0)
	       break;

	     exposure_time_ptr = &exposure_time_left;
	  }
	else if (Num_Rays_Marx_Par > 0)
	  {
	     if (total_num_collected >= num_to_collect)
	       break;
	  }
	else if (total_num_detected >= num_to_detect)
	  break;

	marx_message ("Collecting %d photons...\n", Num_Rays_Per_Iteration);

	if (-1 == marx_create_photons (st, pt, Num_Rays_Per_Iteration,
				       &num_collected, exposure_time_ptr))
	  {
	     marx_error ("Error collecting photons.");
	     return 1;
	  }

	marx_message ("\t%d collected.\n", num_collected);
	if (num_collected == 0)
	  break;

	total_num_collected += num_collected;

	if (-1 == process_photons (pt))
	  goto return_error;

	if (-1 == (*write_function) (pt, num_collected, open_mode, write_mask, total_time, pt->total_time))
	  goto return_error;

	open_mode = 0;
	marx_message ("\n");

	marx_prune_photons (pt);
	total_num_detected += pt->num_sorted;
	total_time += pt->total_time;

 	marx_message ("Total photons: %lu, Total Photons detected: %lu, (efficiency: %f)\n",
		      total_num_collected, total_num_detected, (double)total_num_detected / total_num_collected);

 	marx_message ("  (efficiency this iteration  %f)  Total time: %f\n",
		      (double)pt->num_sorted / num_collected,
		      total_time);

	marx_message ("\n");

	if ((unsigned long)Num_Rays_Per_Iteration != num_collected)
	  break;
     }

   /* drop */
   status = 0;

   return_error:

   if ((-1 == marx_dealloc_photon_type (pt))
       && (status == 0))
     status = -1;

   if ((-1 == (*close_function)(total_time))
       && (status == 0))
     status = -1;

   if ((-1 == marx_close_source (st))
       && (status == 0))
     status = -1;

   marx_message("main(): status = %d\n",status);

   if (status == -1)
     return 1;

#if PRINT_STATS_ARRAY
   fprintf (stdout, "Stats[0] = %ld\n", Stats[0]);
   fprintf (stdout, "Stats[1] = %ld\n", Stats[1]);
   fprintf (stdout, "Stats[2] = %ld\n", Stats[2]);
   fprintf (stdout, "Stats[3] = %ld\n", Stats[3]);
#endif
   return 0;
}  // end of main()

/*}}}*/

static int module_init (Param_File_Type *p) /*{{{*/
{
   Mirror_Module = marx_mirror_init (p);
   if (Mirror_Module == -1)
     return -1;

   Grating_Module = marx_grating_init (p);
   if (Grating_Module == -1)
     return -1;

   Detector_Module = marx_detector_init (p);
   if (Detector_Module == -1)
     return -1;

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ Version/Feature information */

static void print_feature (char *name, int has)
{
   fprintf (stderr, "%20s : %s\n", name, (has ? "yes" : "no"));
}

static void version (void) /*{{{*/
{
   message_copyright (stderr);
#if 0
   fprintf (stderr, "MARX version: %s\n", MARX_VERSION_STRING);
#endif
   fprintf (stderr, "JDMATH library version: %d.%02d\n",
	    (JDMATH_VERSION / 100), (JDMATH_VERSION % 100));
   fprintf (stderr, "PFILE library version: %d.%02d\n",
	    (PFILE_VERSION / 100), (PFILE_VERSION % 100));
   fprintf (stderr, "JDFITS library version: %d.%02d\n",
	    (JDFITS_VERSION / 100), (JDFITS_VERSION % 100));

   fprintf (stderr, "\n");

   fprintf (stderr, "Supported Sources: %s\n", Marx_Supported_Sources);
   fprintf (stderr, "Supported Detectors: %s\n", Marx_Supported_Detectors);

   fprintf (stderr, "\nOther Features:\n");

   print_feature ("HRMA Pitch/Yaw", MARX_HAS_HRMA_PITCH_YAW);
   print_feature ("Wfold Scattering", MARX_HAS_WFOLD);
   print_feature ("Drake Flat", MARX_HAS_DRAKE_FLAT);
   print_feature ("Dynamic Linking", MARX_HAS_DYNAMIC_LINKING);
   print_feature ("ACIS Streak", MARX_HAS_ACIS_STREAK);
   print_feature ("ACIS FEF Support", MARX_HAS_ACIS_FEF);
   print_feature ("Dither Support", MARX_HAS_DITHER);
   exit (0);
}

/*}}}*/

/*}}}*/
#if 0
static int cp_par_file (char *file, char *dir) /*{{{*/
{
   char *ofile = NULL;
   FILE *fpin = NULL, *fpout = NULL;
   char buf[1024];
   unsigned int readlen;
   char *name;
   struct stat in_stat, out_stat;

#if 1
   /* use marx.par no matter what original file is */
   name = "marx.par";
#else
   name = file + strlen (file);
   while ((name >= file) && (*name != '/'))
     name--;
   name++;
#endif

   if (NULL == (ofile = marx_dircat (dir, name)))
     {
	return -1;
     }

   if ((0 == strcmp (file, ofile))
       || ((0 == stat (file, &in_stat))
	   && (0 == stat (ofile, &out_stat))
	   && (in_stat.st_dev == out_stat.st_dev)
	   && (in_stat.st_ino == out_stat.st_ino)))
     {
	marx_error ("Cannot copy a file onto itself");
	marx_error ("Don't shoot yourself in the foot, choose a different OutputDir name");
	goto return_error;
     }

   if (NULL == (fpin = fopen (file, "rb")))
     {
	pf_error ("Unable to open file %s for reading.", file);
	goto return_error;
     }

   if (NULL == (fpout = fopen (ofile, "wb")))
     {
	marx_error ("Unable to open %s for writing.", ofile);
	goto return_error;
     }

   if (EOF == fprintf (fpout, "#@#MARX VERSION: %s\n", marx_make_version_string ()))
     goto write_error;

   do
     {
	readlen = fread (buf, 1, sizeof(buf), fpin);
	if (readlen)
	  {
	     if (readlen != fwrite (buf, 1, readlen, fpout))
	       goto write_error;
	  }
     }
   while (readlen == sizeof (buf));

   if (EOF == fclose (fpout))
     {
	fpout = NULL;
	goto write_error;
     }
   marx_free (ofile);
   fclose (fpin);
   return 0;

   /* Get here only if error occurs. */
   write_error:

   marx_error ("Write to %s failed.", ofile);

   return_error:
   if (ofile != NULL) marx_free (ofile);
   if (fpout != NULL) fclose (fpout);
   if (fpin != NULL) fclose (fpin);
   return -1;
}

/*}}}*/
#endif
static Param_Table_Type Control_Parm_Table [] = /*{{{*/
{
     {"NumRays",	PF_INTEGER_TYPE, 	&Num_Rays_Marx_Par},
     {"dNumRays",	PF_INTEGER_TYPE, 	&Num_Rays_Per_Iteration},
     {"OutputDir",	PF_FILE_TYPE,	 	&Output_Dir},
     {"OutputVectors",	PF_STRING_TYPE,	 	&Data_To_Write},
     {"RandomSeed",	PF_INTEGER_TYPE,	&Random_Seed},
     {"DataDirectory",	PF_STRING_TYPE,		&Data_Directory},
     {"DumpToRayFile",	PF_BOOLEAN_TYPE,	&Dump_To_Rayfile},
     {"RayFile",	PF_STRING_TYPE,		&Rayfile_Name},
     {"SourceType",	PF_STRING_TYPE,		&Source_Name},
     {"ExposureTime",	PF_REAL_TYPE,		&Exposure_Time},
     {"Verbose",	PF_BOOLEAN_TYPE,	&Marx_Verbose},
     {"FocalLength",	PF_REAL_TYPE,		&Marx_Focal_Length},
     {"DetOffsetX",	PF_DOUBLE_TYPE,		&DetOffset_X},
     {"DetOffsetY",	PF_DOUBLE_TYPE,		&DetOffset_Y},
     {"DetOffsetZ",	PF_DOUBLE_TYPE,		&DetOffset_Z},

     {NULL, 0, NULL}
};

/*}}}*/

static Param_File_Type *setup_parms (int argc, char **argv) /*{{{*/
{
   Param_File_Type *p;

   if (NULL == (p = marx_pf_parse_cmd_line ("marx.par", "rwL", argc, argv)))
     {
	fprintf (stderr, "marx: Error opening parameter file.\n");
	return NULL;
     }

   if (NULL == (Parameter_File = pf_get_output_filename (p)))
     {
	fprintf (stderr, "marx: Cannot get output parameter file name.\n");
	pf_close_parameter_file (p);
	return NULL;
     }

   if (-1 == pf_get_parameters (p, Control_Parm_Table))
     {
	pf_error ("marx: error getting parameters.");
	exit (1);
     }

   if (Exposure_Time <= 0.0)
     {
	if (Num_Rays_Per_Iteration <= 0)
	  {
	     (void) fprintf (stderr, "dNumRays must be > 0\n");
	     exit (1);
	  }
	if (Num_Rays_Marx_Par == 0)
	  {
	     (void) fprintf (stderr, "NumRays must be non-zero\n");
	     exit (1);
	  }
     }

   if (Random_Seed == -1)
     JDMsrandom((unsigned long) time (NULL));
   else
     JDMsrandom ((unsigned long) Random_Seed);

   if (-1 == marx_set_data_directory (Data_Directory))
     exit (1);

   if ((Output_Dir != NULL)
       && (Output_Dir[0] == '$')
       && (Output_Dir[1] != 0))

     {
	char *od = getenv (Output_Dir + 1);
	if (od == NULL)
	  {
	     fprintf (stderr, "Environment variable %s does not exist.  Check your OutputDir parameter.\n",
		      Output_Dir + 1);
	     exit (1);
	  }
	marx_free (Output_Dir);

	/* Note: if Output_Dir if freed elsewhere, then we will need to
	 * malloc this.
	 */
	Output_Dir = od;
     }

   if ((Output_Dir == NULL)
       || (*Output_Dir == 0))
     {
	fprintf (stderr, "OutputDir not specified.\n");
	exit (1);
     }

   if (Output_Dir[0] == '|')
     {
	Output_Dir++;
	Is_Output_Piped = 1;
	return p;
     }

   if ((-1 == mkdir (Output_Dir, 0777))
       && (errno != EEXIST))
     {
	fprintf (stderr, "Error creating directory %s\n", Output_Dir);
	exit (1);
     }
   return p;
}

/*}}}*/

static unsigned long compute_write_mask (void) /*{{{*/
{
   char *v, ch;
   unsigned long mask;

   v = Data_To_Write;
   mask = 0;

   mask |= MARX_DET_UV_PIXEL_OK;
   while ((ch = *v++) != 0)
     {
	switch (ch)
	  {
	   case '#':
	     mask |= MARX_TAG_OK;
	     break;

	   case 'E':
	     mask |= MARX_ENERGY_OK;
	     break;
	   case 'T':
	     mask |= MARX_TIME_OK;
	     break;
	   case 'X':
	   case 'Y':
	   case 'Z':
	     mask |= MARX_X_VECTOR_OK;
	     break;
	   case '1':
	   case '2':
	   case '3':
	     mask |= MARX_P_VECTOR_OK;
	     break;
	   case 'P':
	     mask |= MARX_PULSEHEIGHT_OK;
	     break;
	   case 'x':
	   case 'y':
	     mask |= MARX_DET_PIXEL_OK;
	     break;
	   case 'M':
	     mask |= MARX_MIRROR_SHELL_OK;
	     break;
	   case 'O':
	     mask |= MARX_ORDER_OK;
	     break;
	   case 'a':
	     mask |= MARX_ORDER1_OK;
	     break;
	   case 'b':
	     mask |= MARX_ORDER2_OK;
	     break;
	   case 'c':
	     mask |= MARX_ORDER3_OK;
	     break;
	   case 'd':
	     mask |= MARX_ORDER4_OK;
	     break;
	   case 'D':
	     mask |= MARX_DET_NUM_OK;
	     break;
	   case 'S':
	     mask |= MARX_SKY_DITHER_OK|MARX_DET_DITHER_OK;
	     break;
	   case 'r':
	     mask |= MARX_DET_REGION_OK;
	     break;
	   case 'B':
	     mask |= MARX_PI_OK;
	     break;
	   // AML October 22, 2018:
	   case 'G':
	     mask |= MARX_DET_GRADE_OK;
	     break;

	   default:
	     marx_error ("Output vector '%c' not understood.", ch);
	     exit (1);
	  }
     }

   return mask;
}

/*}}}*/

static double normalize_angle (double theta)
{
   theta = fmod (theta, 360.0);
   if (theta < 0.0)
     theta += 360.0;
   return theta;
}

/** Write a set of observation parameters to file obs.par
 *  These parameters will later be used to populate the fits
 *  header keywords when marx2fits is called.
 *  Some keywords have constant values and are just selected based on the 
 *  detector and grating used. Those are loaded from mode dependent parameter
 *  files that are part of marx (e.g. acis_s_obs_par). Others are taken from
 *  global parameters that are set elsewhere (e.g. the nominal pointing).
 */
static int write_obs_par (double total_time)
{
   char *obs_file;
   char *file;
   char *obsdir_file;
   Param_File_Type *pf;
   double ra_nom, dec_nom, roll_nom;
   double ra_pnt, dec_pnt, roll_pnt;
   double dy, dz, dtheta;
   char *detector, *grating;
   char *observer;
   struct tm *tms;
   time_t tstop;
   char date_obs[32];
   char date_end[32];
   double sim_x, sim_y, sim_z;
   Marx_Detector_Type *d;
   char *datamode;
   int with_grating;

   observer = getenv ("USER");
   if (observer == NULL)
     observer = getenv ("LOGNAME");
   if (observer == NULL)
     observer = "DR AXAF";

   switch (Grating_Module)
     {
      case MARX_GRATING_LETG:
	grating = "LETG";
	with_grating = 1;
	break;

      case MARX_GRATING_HETG:
	grating = "HETG";
	with_grating = 1;
	break;

      default:
	grating = "NONE";
	with_grating = 0;
	break;
     }

   switch (Detector_Module)
     {
      case MARX_DETECTOR_HRC_I:
	obsdir_file = "obs/hrc_i_obs.par";
	detector = "HRC-I";
	datamode = "IMAGING";
	if (with_grating) datamode = "SPECTROSCOPIC";
	d = marx_get_detector_info (detector);
	break;

      case MARX_DETECTOR_HRC_S:
	obsdir_file = "obs/hrc_s_obs.par";
	detector = "HRC-S";
	datamode = "IMAGING";
	if (with_grating) datamode = "SPECTROSCOPIC";
	d = marx_get_detector_info (detector);
	break;

      case MARX_DETECTOR_ACIS_I:
	obsdir_file = "obs/acis_i_obs.par";
	detector = "ACIS-I";
	datamode = "GRADED";
	d = marx_get_detector_info (detector);
	break;

      case MARX_DETECTOR_ACIS_S:
	obsdir_file = "obs/acis_s_obs.par";
	detector = "ACIS-S";
	datamode = "GRADED";
	d = marx_get_detector_info (detector);
	break;

      default:
	obsdir_file = "obs/generic_obs.par";
	detector = "NONE";
	datamode = "IMAGING";
	d = NULL;
	break;
     }

   if (d != NULL)
     {
	sim_x = DetOffset_X + d->aimpoint_offset.x;
	sim_y = DetOffset_Y + d->aimpoint_offset.y;
	sim_z = DetOffset_Z + d->aimpoint_offset.z;
     }
   else
     {
	sim_x = DetOffset_X;
	sim_y = DetOffset_Y;
	sim_z = DetOffset_Z;
     }

   if (NULL == (obs_file = marx_dircat (Output_Dir, "obs.par")))
     return -1;

   if (NULL == (file = marx_make_data_file_name (obsdir_file)))
     {
	marx_free (obs_file);
	return -1;
     }

   if (NULL == (pf = pf_open_parameter_file (file, "r")))
     {
	marx_free (obs_file);
	marx_free (file);
	return -1;
     }
   marx_free (file);

   if (-1 == pf_set_output_filename (pf, obs_file))
     {
	marx_free (obs_file);
	pf_close_parameter_file (pf);
	return -1;
     }
   marx_free (obs_file);

   (void) marx_get_nominal_pointing (&ra_nom, &dec_nom, &roll_nom);
   /* Convert to evil degrees */
   ra_nom *= 180.0/PI;
   ra_nom = normalize_angle (ra_nom);
   dec_nom *= 180.0/PI;
   roll_nom *= 180.0/PI;
   roll_nom = normalize_angle (roll_nom);
   (void) marx_get_pointing (&ra_pnt, &dec_pnt, &roll_pnt);
   /* Convert to evil degrees */
   ra_pnt *= 180.0/PI;
   ra_pnt = normalize_angle (ra_pnt);
   dec_pnt *= 180.0/PI;
   roll_pnt *= 180.0/PI;
   roll_pnt = normalize_angle (roll_pnt);

   tms = gmtime (&TStart_Unix);
   sprintf (date_obs, "%4d-%02d-%02dT%02d:%02d:%02d",
	    1900 + tms->tm_year, tms->tm_mon + 1, tms->tm_mday,
	    tms->tm_hour, tms->tm_min, tms->tm_sec);
   tstop = TStart_Unix + (unsigned long) total_time;

   tms = gmtime (&tstop);
   sprintf (date_end, "%4d-%02d-%02dT%02d:%02d:%02d",
	    1900 + tms->tm_year, tms->tm_mon + 1, tms->tm_mday,
	    tms->tm_hour, tms->tm_min, tms->tm_sec);

   dy=0.;
   dz=0.;
   dtheta=0.;

   // marx_message ("write_obs_par(): A\n");
   if (-1 == marx_average_dither(total_time, &dy, &dz, &dtheta))
      return -1;
   // marx_message ("write_obs_par(): B\n");

   if ((-1 == pf_learn_double (pf, "ra_nom", ra_nom))
       || (-1 == pf_learn_double (pf, "roll_nom", roll_nom))
       || (-1 == pf_learn_double (pf, "dec_nom", dec_nom))
       || (-1 == pf_learn_double (pf, "ra_pnt", ra_pnt))
       || (-1 == pf_learn_double (pf, "roll_pnt", roll_pnt))
       || (-1 == pf_learn_double (pf, "dec_pnt", dec_pnt))
       || (-1 == pf_learn_double (pf, "foc_len", Marx_Focal_Length))
       || (-1 == pf_learn_string (pf, "datamode", datamode))
       || (-1 == pf_learn_string (pf, "detnam", detector))
       || (-1 == pf_learn_string (pf, "grating", grating))
       || (-1 == pf_learn_string (pf, "object", Source_Name))
       || (-1 == pf_learn_string (pf, "title", Source_Name))
       || (-1 == pf_learn_string (pf, "observer", observer))
       || (-1 == pf_learn_string (pf, "date-obs", date_obs))
       || (-1 == pf_learn_string (pf, "date-end", date_end))
       || (-1 == pf_learn_double (pf, "mjdref", MJDref))
       || (-1 == pf_learn_double (pf, "tstart", TStart_MJDref))
       || (-1 == pf_learn_double (pf, "tstop", TStart_MJDref + total_time))

       || (-1 == pf_learn_double (pf, "sim_x", sim_x))
       || (-1 == pf_learn_double (pf, "sim_y", sim_y))
       || (-1 == pf_learn_double (pf, "sim_z", sim_z))
       || (-1 == pf_learn_double (pf, "defocus", DetOffset_X))

       || (-1 == pf_learn_double (pf, "DY_AVG", dy))
       || (-1 == pf_learn_double (pf, "DZ_AVG", dz))
       || (-1 == pf_learn_double (pf, "DTH_AVG", dtheta))

       /* Marx specific additions */
       || (-1 == pf_learn_double (pf, "MJD-OBS", TStart_MJDref/86400.0 + MJDref))
       || (-1 == pf_learn_double (pf, "LIVETIME", total_time))
       || (-1 == pf_learn_double (pf, "ONTIME", total_time))
       || (-1 == pf_learn_double (pf, "EXPOSURE", total_time))
       || (-1 == pf_learn_double (pf, "TELAPSE", total_time)))
       {
	  // marx_message ("write_obs_par(): C\n");
	  pf_close_parameter_file (pf);
	  return -1;
       }

   if (-1 == pf_close_parameter_file (pf))
     return -1;

   return 0;
}

/** call write_obs_par, which write obs.par that stores info later needed for fits headers
 */
static int write_fits_info (double total_time)
{
  int iret;

  iret = write_obs_par (total_time);
  // marx_message ("write_fits_info(): return val. of write_obs_par() = %d\n",iret);

  if (iret == -1)
     return -1;

  return 0;
}

/* Pipe routines */

static FILE *open_pipe (char *cmd)
{
   FILE *fp;

   /* fp = popen (cmd, "wb"); */
   fp = popen (cmd, "w");
   if (fp == NULL)
     {
	marx_error ("Error opening pipe to %s", cmd);
	return NULL;
     }
   return fp;
}

static int close_pipe (FILE *fp)
{
   if (EOF == pclose (fp))
     {
	marx_error ("Error closing pipe");
	return -1;
     }
   return 0;
}

static int write_photons_to_pipe (FILE *fp, unsigned int num_input,
				  Marx_Photon_Type *pt,
				  double tstart, double duration)
{
   Marx_Photon_Attr_Type *at;
   unsigned int num_sorted;

   marx_prune_photons (pt);
   num_sorted = pt->num_sorted;
   /*
    * Do not try this optimization since num_input and total_time may
    * be used even when num_sorted is 0
    * if (num_sorted == 0)
    *   return 0;
    */
   if (0 == fwrite ((char *) &num_sorted, sizeof (unsigned int), 1, fp))
     goto write_error;

   if (0 == fwrite ((char *) &num_input, sizeof (unsigned int), 1, fp))
     goto write_error;

   if (0 == fwrite ((char *) &tstart, sizeof (double), 1, fp))
     goto write_error;

   if (0 == fwrite ((char *) &duration, sizeof (double), 1, fp))
     goto write_error;

   at = pt->attributes;

   while (num_sorted)
     {
	if (0 == (at->flags & BAD_PHOTON_MASK))
	  {
	     if (1 != fwrite ((char *) at, sizeof (Marx_Photon_Attr_Type), 1, fp))
	       goto write_error;
	     num_sorted--;
	  }

	at++;
     }

   if (EOF == fflush (fp))
     goto write_error;

   return 0;

   write_error:
   marx_error ("Error writing to pipe");
   return -1;
}
