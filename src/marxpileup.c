/* -*- mode: C; mode: fold -*- */
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

#include <stdio.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <string.h>

/* These are needed for mkdir */
#include <sys/types.h>
#include <sys/stat.h>

#include <errno.h>
#include <jdmath.h>

#include <marx.h>
#include <pfile.h>

#include "acissim.h"

#define NUM_X_PIXELS	1024
#define NUM_Y_PIXELS	1024

#define USE_RMF_CODE	0

static char *Output_Marx_Dir;
static char *Output_Pileup_Dir;
static char *Parameter_File;

static double Alpha;
static double Frame_Time;
static int Verbose;
static int Simulation_Used_Dither = 1;
static unsigned int Total_Num_Input_Rows;

static char *make_marx_filename (char *dir, char *file)
{
   char *dirfile;

   dirfile = marx_dircat (dir, file);

   if (dirfile == NULL)
     {
	fprintf (stderr, "Error make filename for %s.\n", file);
	return NULL;
     }

   return dirfile;
}

static Marx_Dump_File_Type *open_marx_input_file (char *file, int type)
{
   Marx_Dump_File_Type *d;

   if (NULL == (file = make_marx_filename (Output_Marx_Dir, file)))
     return NULL;

   if (Verbose > 1)
     fprintf (stdout, "Opening %s for read\n", file);

   d = marx_open_read_dump_file (file);
   if (d == NULL)
     {
	fprintf (stderr, "Error opening %s\n", file);
	free (file);
	return NULL;
     }

   if ((int) d->type != type)
     {
	fprintf (stderr, "%s is not correct type.  Expecting type '%c'\n",
		 file, type);
	marx_close_read_dump_file (d);
	d = NULL;
     }

   free (file);

   return d;
}

static int open_output_file (Marx_Dump_File_Type *dft,
			     char *filename,
			     char *colname,
			     char type,
			     unsigned int num_cols)
{
   memset ((char *) dft, 0, sizeof (Marx_Dump_File_Type));
   strncpy (dft->colname, colname, sizeof (dft->colname));
   dft->colname[sizeof(dft->colname) - 1] = 0;
   dft->type = type;
   dft->num_cols = (int32) num_cols;

   filename = make_marx_filename (Output_Pileup_Dir, filename);
   if (filename == NULL)
     return -1;

   if (Verbose > 1)
     fprintf (stdout, "Opening %s for write\n", filename);

   if (NULL == marx_create_write_dump_file (filename, dft))
     {
	free ((char *) filename);
	return -1;
     }
   free ((char *) filename);
   return 0;
}

typedef struct
{
   char *name;
   int type;
   Marx_Dump_File_Type *dft;
}
Input_File_Type;

static Input_File_Type Input_Files [] =
{
     {"b_energy.dat", 'E', NULL},
#define INPUT_BENERGY_DAT		0
     {"time.dat", 'E', NULL},
#define INPUT_TIME_DAT		1
     {"xpixel.dat", 'E', NULL},
#define INPUT_XPIXEL_DAT	2
     {"ypixel.dat", 'E', NULL},
#define INPUT_YPIXEL_DAT	3
     {"detector.dat", 'A', NULL},
#define INPUT_DETECTOR_DAT	4
     {"energy.dat", 'E', NULL},
#define INPUT_ENERGY_DAT	5
     {NULL, 0, NULL}
};

static FILE *Detector_Dat_Fp;
static FILE *BEnergy_Dat_Fp;
static FILE *Energy_Dat_Fp;
static FILE *Time_Dat_Fp;
static FILE *XPixel_Dat_Fp;
static FILE *YPixel_Dat_Fp;
static FILE *Sky_RA_Dat_Fp;
static FILE *Sky_Dec_Dat_Fp;
static FILE *Sky_Roll_Dat_Fp;
static FILE *Det_DY_Dat_Fp;
static FILE *Det_DZ_Dat_Fp;
static FILE *Det_Theta_Dat_Fp;

static Input_File_Type Dither_Input_files [] =
{
#define INPUT_SKY_RA_DAT	0
#define INPUT_SKY_DEC_DAT	1
#define INPUT_SKY_ROLL_DAT	2
#define INPUT_DET_DY_DAT	3
#define INPUT_DET_DZ_DAT	4
#define INPUT_DET_THETA_DAT	5
     {"sky_ra.dat", 'E', NULL},
     {"sky_dec.dat", 'E', NULL},
     {"sky_roll.dat", 'E', NULL},
     {"det_dy.dat", 'E', NULL},
     {"det_dz.dat", 'E', NULL},
     {"det_theta.dat", 'E', NULL},
     {NULL, 0, NULL}
};

typedef struct
{
   char *name;			       /* sans .dat */
   char *colname;
   int type;
   int requires_dither;
   Marx_Dump_File_Type dft;
}
Output_File_Type;

static Output_File_Type Output_Files [] =
{
#define XPIXEL_DAT	0
     {"xpixel", "CHIPX", 'E', 0},
#define YPIXEL_DAT	1
     {"ypixel", "CHIPY", 'E', 0},
#define TIME_DAT	2
     {"time", "TIME", 'E', 0},
#define FRAME_DAT	3
     {"frame", "FRAME", 'J', 0},
#define BENERGY_DAT	4
     {"b_energy", "B_ENERGY", 'E', 0},
#define DETECTOR_DAT	5
     {"detector", "CCDID", 'A', 0},
#define SKY_RA_DAT	6
     {"sky_ra", "RA", 'E', 1},
#define SKY_DEC_DAT	7
     {"sky_dec", "DEC", 'E', 1},
#define SKY_ROLL_DAT	8
     {"sky_roll", "ROLL", 'E', 1},
#define DET_DY_DAT	9
     {"det_dy", "DET_DY", 'E', 1},
#define DET_DZ_DAT	10
     {"det_dz", "DET_DZ", 'E', 1},
#define DET_THETA_DAT	11
     {"det_theta", "DET_THETA", 'E', 1},
#define NPHOTONS_DAT	12
     {"nphotons", "NUM_PHOTONS", 'I', 0},
#define PHA_DAT		13
     {"pha", "PHA", 'I', 0},

     {NULL}
};

static FILE *XPixel_Fp;
static FILE *YPixel_Fp;
static FILE *Time_Fp;
static FILE *Frame_Fp;
static FILE *BEnergy_Fp;
static FILE *Detector_Fp;
static FILE *Sky_RA_Fp;
static FILE *Sky_Dec_Fp;
static FILE *Sky_Roll_Fp;
static FILE *Det_DY_Fp;
static FILE *Det_DZ_Fp;
static FILE *Det_Theta_Fp;
static FILE *NPhotons_Fp;
static FILE *Pha_Fp;

static unsigned int Num_Detected;
static unsigned int Num_Input;

static int close_marx_output_files (void)
{
   Output_File_Type *oft;
   int ret = 0;
   unsigned int nrows;

   nrows = Num_Detected;
   oft = Output_Files;
   while (oft->name != NULL)
     {
	FILE *fp;

	fp = oft->dft.fp;
	if (fp != NULL)
	  {
	     if (-1 == marx_close_write_dump_file (fp, nrows))
	       ret = -1;
	  }

	oft++;
     }

   Frame_Fp = Time_Fp = YPixel_Fp = XPixel_Fp = NULL;
   BEnergy_Fp = NULL;

   return ret;
}

static int open_marx_output_files (void)
{
   Output_File_Type *oft;
   char file [1024];

   if (Output_Pileup_Dir == NULL)
     {
	if (NULL == (Output_Pileup_Dir = make_marx_filename (Output_Marx_Dir, "pileup")))
	  return -1;
     }

   if ((-1 == mkdir (Output_Pileup_Dir, 0777))
       && (errno != EEXIST))
     {
	fprintf (stderr, "Error creating directory %s\n", Output_Pileup_Dir);
	return -1;
     }

   oft = Output_Files;
   while (oft->name != NULL)
     {
	int cols = 0;

	if (oft->requires_dither && (Simulation_Used_Dither == 0))
	  {
	     oft++;
	     continue;
	  }

	sprintf (file, "%s.dat", oft->name);
	if (-1 == open_output_file (&oft->dft, file, oft->colname,
				    oft->type, cols))
	  {
	     close_marx_output_files ();
	     return -1;
	  }
	oft++;
     }

   Time_Fp = Output_Files [TIME_DAT].dft.fp;
   XPixel_Fp = Output_Files [XPIXEL_DAT].dft.fp;
   YPixel_Fp = Output_Files [YPIXEL_DAT].dft.fp;
   BEnergy_Fp = Output_Files [BENERGY_DAT].dft.fp;
   Frame_Fp = Output_Files [FRAME_DAT].dft.fp;
   Detector_Fp = Output_Files[DETECTOR_DAT].dft.fp;
   NPhotons_Fp = Output_Files[NPHOTONS_DAT].dft.fp;
   Pha_Fp = Output_Files[PHA_DAT].dft.fp;

   if (Simulation_Used_Dither)
     {
	Sky_RA_Fp = Output_Files[SKY_RA_DAT].dft.fp;
	Sky_Dec_Fp = Output_Files[SKY_DEC_DAT].dft.fp;
	Sky_Roll_Fp = Output_Files[SKY_ROLL_DAT].dft.fp;
	Det_DY_Fp = Output_Files[DET_DY_DAT].dft.fp;
	Det_DZ_Fp = Output_Files[DET_DZ_DAT].dft.fp;
	Det_Theta_Fp = Output_Files[DET_THETA_DAT].dft.fp;
     }

   return 0;
}

static void close_some_input_files (Input_File_Type *ift)
{
   while (ift->name != NULL)
     {
	if (ift->dft != NULL)
	  {
	     marx_close_read_dump_file (ift->dft);
	     ift->dft = NULL;
	  }
	ift++;
     }
}

static void close_marx_input_files (void)
{
   close_some_input_files (Input_Files);
   close_some_input_files (Dither_Input_files);

   Detector_Dat_Fp = NULL;
   BEnergy_Dat_Fp = NULL;
   Energy_Dat_Fp = NULL;
   Time_Dat_Fp = NULL;
   XPixel_Dat_Fp = NULL;
   YPixel_Dat_Fp = NULL;
   NPhotons_Fp = NULL;
   Pha_Fp = NULL;

   Sky_RA_Dat_Fp = NULL;
   Sky_Dec_Dat_Fp = NULL;
   Sky_Roll_Dat_Fp = NULL;
   Det_DY_Dat_Fp = NULL;
   Det_DZ_Dat_Fp = NULL;
   Det_Theta_Dat_Fp = NULL;
}

static int open_some_input_files (Input_File_Type *ift, int *num_rowsp)
{
   int num_rows;

   num_rows = *num_rowsp;

   while (ift->name != NULL)
     {
	if (NULL == (ift->dft = open_marx_input_file (ift->name, ift->type)))
	  {
	     close_marx_input_files ();
	     return -1;
	  }

	if (num_rows == -1)
	  num_rows = ift->dft->num_rows;
	if (num_rows != ift->dft->num_rows)
	  {
	     fprintf (stderr, "Input data files are not of same length.\n");
	     close_marx_input_files ();
	     return -1;
	  }
	ift++;
     }

   *num_rowsp = num_rows;
   return 0;
}

static int open_marx_input_files (void)
{
   int num_rows = -1;

   if (-1 == open_some_input_files (Input_Files, &num_rows))
     return -1;

   Detector_Dat_Fp = Input_Files[INPUT_DETECTOR_DAT].dft->fp;
   BEnergy_Dat_Fp = Input_Files[INPUT_BENERGY_DAT].dft->fp;
   Energy_Dat_Fp = Input_Files[INPUT_ENERGY_DAT].dft->fp;
   Time_Dat_Fp = Input_Files[INPUT_TIME_DAT].dft->fp;
   XPixel_Dat_Fp = Input_Files[INPUT_XPIXEL_DAT].dft->fp;
   YPixel_Dat_Fp = Input_Files[INPUT_YPIXEL_DAT].dft->fp;

   if (Simulation_Used_Dither)
     {
	if (-1 == open_some_input_files (Dither_Input_files, &num_rows))
	  return -1;
	Sky_RA_Dat_Fp = Dither_Input_files[INPUT_SKY_RA_DAT].dft->fp;
	Sky_Dec_Dat_Fp = Dither_Input_files[INPUT_SKY_DEC_DAT].dft->fp;
	Sky_Roll_Dat_Fp = Dither_Input_files[INPUT_SKY_ROLL_DAT].dft->fp;
	Det_DY_Dat_Fp = Dither_Input_files[INPUT_DET_DY_DAT].dft->fp;
	Det_DZ_Dat_Fp = Dither_Input_files[INPUT_DET_DZ_DAT].dft->fp;
	Det_Theta_Dat_Fp = Dither_Input_files[INPUT_DET_THETA_DAT].dft->fp;
     }

   if (Verbose > 1)
     fprintf (stdout, "%d input events available (across all CCDs)\n", num_rows);

   Total_Num_Input_Rows = num_rows;
   return 0;
}

#define MAX_NUM_CCDS	10

typedef struct
{
   float benergy;
   unsigned int num_photons;
   float island_benergy;
   unsigned int island_num_photons;
   float x, y;
   /* The xpixel values are energy weighted as follows:
    *    X_k = (e0*x0 + e1*x1 + ... + ek*xk)/E_k
    *    E_k = e0+e1+...e_k
    *    X_{k+1} = (E_k*X_k + e_{k+1}*x_{k+1})/(E_k+e_{k+1})
    */
}
CCD_Pixel_Type;

static CCD_Pixel_Type *CCD_Arrays[MAX_NUM_CCDS];

typedef struct Input_Event_Type
{
   float benergy;
   float energy;
   float t;
   float x;
   float y;
   int ccdid;
   float sky_ra, sky_dec, sky_roll;
   float det_dy, det_dz, det_theta;
   CCD_Pixel_Type *ccd_pixels;
   int can_be_center;		       /* if it can be the center of a 3x3 island */
   struct Input_Event_Type *next;
}
Input_Event_Type;

static int allocate_ccd_array (unsigned int ccdid)
{
   unsigned int size;
   CCD_Pixel_Type *p;

   if (ccdid >= MAX_NUM_CCDS)
     {
	fprintf (stderr, "Invalid CCD_ID (%u)\n", ccdid);
	return -1;
     }

   size = NUM_X_PIXELS * NUM_Y_PIXELS * sizeof (CCD_Pixel_Type);

   if (NULL == (p = (CCD_Pixel_Type *) marx_malloc (size)))
     return -1;

   memset ((char *) p, 0, size);
   CCD_Arrays[ccdid] = p;
   return 0;
}

static void free_ccd_pixels (CCD_Pixel_Type *p)
{
   if (p != NULL)
     marx_free ((char *) p);
}

static void deallocate_buffers (void)
{
   unsigned int i;

   for (i = 0; i < MAX_NUM_CCDS; i++)
     {
	free_ccd_pixels (CCD_Arrays[i]);
	CCD_Arrays[i] = NULL;
     }
}

static CCD_Pixel_Type *get_pixel_list (unsigned int ccdid)
{
   if (ccdid >= MAX_NUM_CCDS)
     {
	fprintf (stderr, "Invalid CCD_ID (%u)\n", ccdid);
	return NULL;
     }

   if ((CCD_Arrays[ccdid] == NULL)
       && (-1 == allocate_ccd_array (ccdid)))
     return NULL;

   return CCD_Arrays[ccdid];
}

static Input_Event_Type *allocate_input_event (int ccd, float x, float y, float t, float benergy,
					       float energy,
					       float ra, float dec, float roll,
					       float dy, float dz, float theta)
{
   Input_Event_Type *evt;
   CCD_Pixel_Type *ccd_pixels;

   if (NULL == (ccd_pixels = get_pixel_list (ccd)))
     return NULL;

   evt = (Input_Event_Type *) marx_malloc (sizeof (Input_Event_Type));
   if (evt == NULL)
     return NULL;

   evt->ccd_pixels = ccd_pixels;
   evt->ccdid = ccd;
   evt->x = x;
   evt->y = y;
   evt->t = t;
   evt->benergy = benergy;
   evt->energy = energy;
   evt->sky_ra = ra;
   evt->sky_dec = dec;
   evt->sky_roll = roll;
   evt->det_dy = dy;
   evt->det_dz = dz;
   evt->det_theta = theta;

   evt->can_be_center = ((x >= 1) && (x < NUM_X_PIXELS - 1)
			 && (y >= 1) && (y < NUM_Y_PIXELS - 1));

   if ((evt->can_be_center == 0)
       && ((x >= NUM_Y_PIXELS) || (y >= NUM_Y_PIXELS)))
     {
	marx_free ((char *) evt);
	fprintf (stderr, "Corrupt file?  X or Y pixel coord is out of range.\n");
	return NULL;
     }

   return evt;
}

static int
read_input_event (Input_Event_Type **evtp, unsigned int *frame)
{
   float32 x, y, t, ra, dec, roll, dy, dz, theta, benergy, energy;
   char ccd;

   if (1 != fread (&ccd, 1, 1, Detector_Dat_Fp))
     return 0;

   if ((1 != JDMread_float32 (&x, 1, XPixel_Dat_Fp))
       || (1 != JDMread_float32 (&y, 1, YPixel_Dat_Fp))
       || (1 != JDMread_float32 (&t, 1, Time_Dat_Fp))
       || (1 != JDMread_float32 (&benergy, 1, BEnergy_Dat_Fp))
       || (1 != JDMread_float32 (&energy, 1, Energy_Dat_Fp)))
     {
	fprintf (stderr, "Read Error.\n");
	return -1;
     }

   if (Simulation_Used_Dither)
     {
	if ((1 != JDMread_float32 (&ra, 1, Sky_RA_Dat_Fp))
	    || (1 != JDMread_float32 (&dec, 1, Sky_Dec_Dat_Fp))
	    || (1 != JDMread_float32 (&roll, 1, Sky_Roll_Dat_Fp))
	    || (1 != JDMread_float32 (&dy, 1, Det_DY_Dat_Fp))
	    || (1 != JDMread_float32 (&dz, 1, Det_DZ_Dat_Fp))
	    || (1 != JDMread_float32 (&theta, 1, Det_Theta_Dat_Fp)))
	  {
	     fprintf (stderr, "Read Error.\n");
	     return -1;
	  }
     }
   else
     {
	ra = dec = roll = dy = dz = theta = 0.0;
     }

   if (NULL == (*evtp = allocate_input_event (ccd, x, y, t,
					      benergy, energy,
					      ra, dec, roll,
					      dy, dz, theta
					     )))
	  return -1;

   *frame = (unsigned int) (t / Frame_Time);
   Num_Input++;
   return 1;
}

static int write_event (char ccd, float32 xpix, float32 ypix, float32 energy,
			int32 frame, float32 ra, float32 dec, float32 roll,
			float32 dy, float32 dz, float32 theta,
			int16 nphotons)
{
   float32 t = (float32) (frame * Frame_Time);
   float32 benergy;
   short pha;
   int16 pha16;

   benergy = (float32) energy;
   if (-1 == marx_map_energy_to_acis_pha (ccd, xpix, ypix, benergy, &pha))
     return -1;
   pha16 = (int16) pha;

   if ((1 != fwrite (&ccd, 1, 1, Detector_Fp))
       || (1 != JDMwrite_float32 (&xpix, 1, XPixel_Fp))
       || (1 != JDMwrite_float32 (&ypix, 1, YPixel_Fp))
       || (1 != JDMwrite_int32 (&frame, 1, Frame_Fp))
       || (1 != JDMwrite_float32 (&t, 1, Time_Fp))
       || (1 != JDMwrite_int16 (&nphotons, 1, NPhotons_Fp))
       || (1 != JDMwrite_int16 (&pha16, 1, Pha_Fp))
       || (1 != JDMwrite_float32 (&benergy, 1, BEnergy_Fp)))
     {
	fprintf (stderr, "Write Error.\n");
	return -1;
     }

   if (Simulation_Used_Dither)
     {
	if ((1 != JDMwrite_float32 (&ra, 1, Sky_RA_Fp))
	    || (1 != JDMwrite_float32 (&dec, 1, Sky_Dec_Fp))
	    || (1 != JDMwrite_float32 (&roll, 1, Sky_Roll_Fp))
	    || (1 != JDMwrite_float32 (&dy, 1, Det_DY_Fp))
	    || (1 != JDMwrite_float32 (&dz, 1, Det_DZ_Fp))
	    || (1 != JDMwrite_float32 (&theta, 1, Det_Theta_Fp)))
	  {
	     fprintf (stderr, "Write Error\n");
	     return -1;
	  }
     }

   Num_Detected++;
   return 0;
}

static int will_grade_migrate (unsigned int num)
{
   double prob;

   prob = pow (Alpha, num - 1);
   return (JDMrandom () >= prob);
}

static int
event_detect (Input_Event_Type *event_list, unsigned int frame)
{
   Input_Event_Type *evt;

   evt = event_list;
   while (evt != NULL)
     {
	CCD_Pixel_Type *xy0, *xy1, *xy2;
	unsigned int island_num_photons;
	double island_benergy, benergy;
	double x, y;

	if (0 == evt->can_be_center)
	  {
	     evt = evt->next;
	     continue;
	  }
	x = evt->x;
	y = evt->y;
	xy0 = evt->ccd_pixels + ((unsigned int)y-1) * NUM_X_PIXELS + ((unsigned int)x-1);
	xy1 = xy0 + NUM_X_PIXELS;
	xy2 = xy1 + NUM_X_PIXELS;

	benergy = xy1[1].benergy;
	island_benergy = xy1[1].island_benergy;
	island_num_photons = xy1[1].island_num_photons;
#if 1
	if (((xy0[0].benergy >= benergy) || (xy0[1].benergy >= benergy) || (xy0[2].benergy >= benergy))
	    || (xy1[0].benergy > benergy) ||                                             (xy1[2].benergy >= benergy)
	    || (xy2[0].benergy > benergy) || (xy2[1].benergy > benergy) || (xy2[2].benergy > benergy))
	  {
	     evt = evt->next;
	     continue;
	  }
#endif
#if 1
	if (((xy0[0].island_benergy >= island_benergy) || (xy0[1].island_benergy >= island_benergy) || (xy0[2].island_benergy >= island_benergy))
	    || (xy1[0].island_benergy > island_benergy) ||                                             (xy1[2].island_benergy >= island_benergy)
	    || (xy2[0].island_benergy > island_benergy) || (xy2[1].island_benergy > island_benergy) || (xy2[2].island_benergy > island_benergy))
	  {
	     evt = evt->next;
	     continue;
	  }
#endif

	if (island_num_photons >= 2)
	  {
	     unsigned int i;
	     if (will_grade_migrate (island_num_photons))
	       {
		  evt = evt->next;
		  continue;
	       }
	     /* Compute a weighted average */
	     x = 0; y = 0;
	     for (i = 0; i < 3; i++)
	       {
		  x += xy0[i].x*xy0[i].benergy; y += xy0[i].y*xy0[i].benergy;
		  x += xy1[i].x*xy1[i].benergy; y += xy1[i].y*xy1[i].benergy;
		  x += xy2[i].x*xy2[i].benergy; y += xy2[i].y*xy2[i].benergy;
	       }
	     x /= xy1[1].island_benergy;
	     y /= xy1[1].island_benergy;
	  }

	if (-1 == write_event (evt->ccdid, x, y, island_benergy, frame,
			       evt->sky_ra, evt->sky_dec, evt->sky_roll,
			       evt->det_dy, evt->det_dz, evt->det_theta,
			       island_num_photons))
	  return -1;

	evt = evt->next;
     }

   return 0;
}

static int store_event (Input_Event_Type *evt, int *is_duplicatep)
{
   CCD_Pixel_Type *xy0, *xy1, *xy2;
   unsigned int x, y;

   x = (unsigned int) evt->x;
   y = (unsigned int) evt->y;

   if (evt->can_be_center)
     {
	double e, cent_e;
	double cent, corn, side;
	CCD_Pixel_Type *cp;

#if USE_RMF_CODE
	e = evt->energy;
#else
	e = evt->benergy;
#endif
	/* This ought to be replaced by real energy-dependent grade
	 * splitting patterns.  However, such data do not appear to be
	 * available.
	 */
	cent = 1.0;
	corn = (1 - cent) / 8.0;
	side = (1 - cent) / 8.0;

	xy0 = evt->ccd_pixels + (y-1) * NUM_X_PIXELS + (x-1);
	xy1 = xy0 + NUM_X_PIXELS;
	xy2 = xy1 + NUM_X_PIXELS;

	cp = xy1+1;
	*is_duplicatep = (cp->num_photons != 0);
	cp->num_photons += 1;

	cent_e = cent * e;
	cp->x = cp->x*cp->benergy + cent_e*evt->x;
	cp->y = cp->y*cp->benergy + cent_e*evt->y;
	cp->benergy += cent_e;
	cp->x /= cp->benergy;
	cp->y /= cp->benergy;

	xy0[0].benergy += corn*e; xy0[1].benergy += side*e; xy0[2].benergy += corn*e;
	xy1[0].benergy += side*e;                           xy1[2].benergy += side*e;
	xy2[0].benergy += corn*e; xy2[1].benergy += side*e; xy2[2].benergy += corn*e;

	return 0;
     }
   *is_duplicatep = 0;
#if 0
   xy->benergy += evt->benergy;
   xy->num_photons += 1;
#endif
   return 0;
}

#define CORNER_FRACTION (9.0/9.0)
#define SIDE_FRACTION	(9.0/9.0)
#define CENTER_FRACTION	(1.0)

static int collect_charge (Input_Event_Type *event_list)
{
   Input_Event_Type *evt;

   evt = event_list;
   while (evt != NULL)
     {
	CCD_Pixel_Type *xy0, *xy1, *xy2;
	unsigned int x, y;

	x = (unsigned int) evt->x;
	y = (unsigned int) evt->y;
	if (evt->can_be_center)
	  {
	     xy0 = evt->ccd_pixels + (y-1) * NUM_X_PIXELS + (x-1);
	     xy1 = xy0 + NUM_X_PIXELS;
	     xy2 = xy1 + NUM_X_PIXELS;

	     xy1[1].island_benergy =
	       CORNER_FRACTION*xy0[0].benergy + SIDE_FRACTION*xy0[1].benergy + CORNER_FRACTION*xy0[2].benergy
	       + SIDE_FRACTION*xy1[0].benergy + CENTER_FRACTION*xy1[1].benergy + SIDE_FRACTION*xy1[2].benergy
	       + CORNER_FRACTION*xy2[0].benergy + SIDE_FRACTION*xy2[1].benergy + CORNER_FRACTION*xy2[2].benergy;
	     xy1[1].island_num_photons =
	       xy0[0].num_photons + xy0[1].num_photons + xy0[2].num_photons
	       + xy1[0].num_photons + xy1[1].num_photons + xy1[2].num_photons
	       + xy2[0].num_photons + xy2[1].num_photons + xy2[2].num_photons;
	  }

	evt = evt->next;
     }
   return 0;
}

static void free_event_list (Input_Event_Type *evt)
{
   CCD_Pixel_Type *pmin, *pmax;

   while (evt != NULL)
     {
	CCD_Pixel_Type *xy;
	Input_Event_Type *evt1;

	pmin = evt->ccd_pixels;
	pmax = pmin + NUM_X_PIXELS*NUM_Y_PIXELS;

	xy = pmin + (unsigned int)(evt->y-1) * NUM_X_PIXELS + (unsigned int)(evt->x-1);
	if (evt->can_be_center)
	  {
	     memset ((char *) xy, 0, 3 * sizeof (CCD_Pixel_Type));
	     xy += NUM_X_PIXELS;
	     memset ((char *) xy, 0, 3 * sizeof (CCD_Pixel_Type));
	     xy += NUM_X_PIXELS;
	     memset ((char *) xy, 0, 3 * sizeof (CCD_Pixel_Type));
	  }
	else
	  {
	     unsigned int i, j;

	     for (j = 0; j < 3; j++)
	       {
		  for (i = 0; i < 3; i++)
		    {
		       CCD_Pixel_Type *xyi = xy + i;
		       if ((xyi >= pmin) && (xyi < pmax))
			 memset ((char *) xyi, 0, sizeof (CCD_Pixel_Type));
		    }
		  xy += NUM_X_PIXELS;
	       }
	  }

	evt1 = evt->next;
	marx_free ((char *) evt);
	evt = evt1;
     }
}

static int process_frame (Input_Event_Type *event_list, unsigned int frame)
{
   Input_Event_Type *evt, *levt;

   levt = NULL;
   evt = event_list;

   while (evt != NULL)
     {
	int is_dup;
	if (-1 == store_event (evt, &is_dup))
	  return -1;

	if (is_dup)
	  {
	     levt->next = evt->next;
	     marx_free ((char *)evt);
	     evt = levt;
	  }
	else
	  levt = evt;

	evt = evt->next;
     }

   if (-1 == collect_charge (event_list))
     return -1;

   if (-1 == event_detect (event_list, frame))
     return -1;

   return 0;
}

static double Frame_Transfer_Time;

static Param_Table_Type Pileup_Parm_Table [] =
{
     {"MarxOutputDir",		PF_FILE_TYPE,	&Output_Marx_Dir},
     {"Verbose",		PF_INTEGER_TYPE,&Verbose},
     {"FrameTime",		PF_REAL_TYPE,	&Frame_Time},
     {"FrameTransferTime",	PF_REAL_TYPE,	&Frame_Transfer_Time},
     {"Alpha",			PF_REAL_TYPE,	&Alpha},
     {NULL, 0, NULL}
};

static int cp_file (char *file, char *dir, char *name) /*{{{*/
{
   char *ofile = NULL;
   FILE *fpin = NULL, *fpout = NULL;
   char buf[1024];
   unsigned int readlen;
   struct stat in_stat, out_stat;

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
	marx_error ("Cannot copy a file onto itself (%s --> %s)", file, ofile);
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

static int copy_files (void)
{
   char *marx_dir_file;

   if (NULL == (marx_dir_file = make_marx_filename (Output_Marx_Dir, "marx.par")))
     return -1;

   if (-1 == cp_file (marx_dir_file, Output_Pileup_Dir, "marx.par"))
     {
	free ((char *) marx_dir_file);
	return -1;
     }
   free ((char *) marx_dir_file);

   if (NULL == (marx_dir_file = make_marx_filename (Output_Marx_Dir, "obs.par")))
     return -1;

   if (-1 == cp_file (marx_dir_file, Output_Pileup_Dir, "obs.par"))
     {
	free ((char *) marx_dir_file);
	return -1;
     }
   free ((char *) marx_dir_file);

   /* Now copy the parameter file for this program to the output dir. */
   if (-1 == cp_file (Parameter_File, Output_Pileup_Dir, "marxpileup.par"))
     {
	free ((char *) marx_dir_file);
	return -1;
     }

   return 0;
}

static int initialize (int argc, char **argv)
{
   Param_File_Type *p;
   char detector[64];
   char *file;
   int status;

   p = marx_pf_parse_cmd_line ("marxpileup.par", NULL, argc, argv);

   if (p == NULL)
     {
	fprintf (stderr, "pileup: Error opening parameter file.\n");
	return -1;
     }

   if (-1 == pf_get_parameters (p, Pileup_Parm_Table))
     {
	pf_error ("pileup: error getting parameters.");
	pf_close_parameter_file (p);
	return -1;
     }

   /* Get a copy of the parameter file output name so we can copy it
    * to the output directory later on.
    */
   if (NULL == (Parameter_File = pf_get_output_filename (p)))
     {
	fprintf (stderr, "pileup: Cannot get output parameter file name.\n");
	pf_close_parameter_file (p);
	return -1;
     }

   pf_close_parameter_file (p);

   if ((Alpha <= 0) || (Alpha > 1.0))
     {
	fprintf (stderr, "Alpha is out of range.  It is required to be in the range (0,1]\n");
	return -1;
     }

   if (Frame_Transfer_Time > 0.0)
     Frame_Time += Frame_Transfer_Time;

   if (NULL == (file = make_marx_filename (Output_Marx_Dir, "marx.par")))
     return -1;

   if (NULL == (p = pf_open_parameter_file (file, "r")))
     {
	marx_error ("*** Unable to open %s\n", file);
	free (file);
	return -1;
     }

   if (-1 == pf_get_string (p, "DetectorType", detector, sizeof(detector)))
     {
	free (file);
	(void) pf_close_parameter_file (p);
	return -1;
     }

   free (file);

   status = 0;
   if (0 == strcmp (detector, "ACIS-S"))
     {
	status = marx_init_acis_s_rmf (p);
     }
   else if (0 == strcmp (detector, "ACIS-I"))
     {
	status = marx_init_acis_i_rmf (p);
     }
   else
     {
	marx_error ("DetectorType=%s is not supported\n", detector);
	status = -1;
     }

   pf_close_parameter_file (p);
   return status;
}

int main (int argc, char **argv)
{
   Input_Event_Type *event_list;
   unsigned int frame;
   unsigned int last_frame;
   int last_percent;

   if (-1 == initialize (argc, argv))
     return 1;

   if (-1 == open_marx_input_files ())
     return 1;

   if (-1 == open_marx_output_files ())
     {
	close_marx_input_files ();
	return 2;
     }

   last_frame = 0;
   event_list = NULL;
   last_percent = -1;
   while (1)
     {
	Input_Event_Type *evt;
	int status;

        status = read_input_event (&evt, &frame);
	if (status == -1)
	  goto return_error;

	if ((status == 1)
	    && (frame == last_frame))
	  {
	     evt->next = event_list;
	     event_list = evt;
	     continue;
	  }

	if ((event_list != NULL)
	    || (status == 0))
	  {
	     if (-1 == process_frame (event_list, last_frame))
	     goto return_error;
	     free_event_list (event_list);

	     if (Verbose > 1)
	       {
		  int percent = (Num_Input * 100)/Total_Num_Input_Rows;
		  if (last_percent != percent)
		    {
		       fprintf (stdout, "Percent Completed: % 2d   \r", percent);
		       fflush (stdout);
		    }
		  last_percent = percent;
	       }
	     if (status == 0)
	       break;
	  }

	event_list = evt;
	evt->next = NULL;
	last_frame = frame;
     }

   close_marx_input_files ();
   close_marx_output_files ();
   deallocate_buffers ();

   if (-1 == copy_files ())
     return 1;

   if (Verbose)
     {
	fprintf (stdout, "Total Number Input: %u\n",
		 Num_Input);
	fprintf (stdout, "Total Number Detected: %u\n",
		 Num_Detected);
	if (Num_Input != 0)
	  fprintf (stdout, "Efficiency: %e\n", (double) Num_Detected/(double) Num_Input);
     }
   return 0;

   return_error:

   close_marx_input_files ();
   if (event_list != NULL)
     free_event_list (event_list);
   deallocate_buffers ();

   return 1;
}

