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
#include <string.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>

#include "marx.h"
#include "_marx.h"

#include "_jdfits.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

/*}}}*/

typedef struct /*{{{*/
{
   int32 num_rows;
   int32 num_cols;
   FILE *fp;
}
/*}}}*/
Write_File_Type;

#define NBITS_LONG (8 * sizeof (long))
static Write_File_Type File_Pointers [NBITS_LONG];

static unsigned char MagicChars [4] = /*{{{*/
{
   0x83, 0x13, 0x89, 0x8D
};

/*}}}*/
#define HEADER_SIZE 32
/* Num element offset MUST follow 4 magic bytes + 16 byte column name. */
#define NUM_ROWS_OFFSET 20
#define NUM_COLS_OFFSET 24
#define RESERVED_OFFSET (NUM_COLS_OFFSET + 4)
#define HEADER_RESERVED (HEADER_SIZE - RESERVED_OFFSET)

int marx_close_write_dump_file (FILE *fp, unsigned long n)
{
#define MCWDF_ERRFMT "marx_close_write_dump_file: %s (errno = %d)"
   int ret = 0;
   int32 i32;

   if (fp == NULL)
     return -1;

   i32 = (int32) n;

   errno = 0;
#ifndef SEEK_SET
# define SEEK_SET	0
#endif
   if (-1 == FSEEK (fp, NUM_ROWS_OFFSET, SEEK_SET))
     {
	marx_error (MCWDF_ERRFMT, "seek error", errno);
	ret = -1;
     }
   else if (1 != JDMwrite_int32 (&i32, 1, fp))
     {
	marx_error (MCWDF_ERRFMT, "write error", errno);
	ret = -1;
     }

   if (EOF == fclose (fp))
     {
	marx_error (MCWDF_ERRFMT, "write error", errno);
	ret = -1;
     }

   if (ret == -1)
     {
#ifdef ENOSPC
	marx_error ("No space left on device.");
#endif
     }

   return ret;
}

static int close_write_files (int do_num) /*{{{*/
{
   unsigned int i;
   FILE *fp;
   int ret = 0;

   for (i = 0; i < NBITS_LONG; i++)
     {
	if (NULL != (fp = File_Pointers[i].fp))
	  {
	     File_Pointers[i].fp = NULL;

	     if (do_num)
	       {
		  unsigned long num_rows = File_Pointers[i].num_rows;

		  if (-1 == marx_close_write_dump_file (fp, num_rows))
		    ret = -1;
	       }
	     else if (EOF == fclose (fp))
	       ret = -1;

	     /* File_Pointers[i].num_rows = 0; */
	  }
     }
   return ret;
}

/*}}}*/

FILE *marx_create_write_dump_file (char *filename, Marx_Dump_File_Type *dft)
{
   char header[HEADER_SIZE];
   char *h, *hmax, *p;
   FILE *fp;
   unsigned int len;

   dft->fp = NULL;

   if (NULL == (fp = fopen (filename, "w+b")))
     {
	marx_error ("Unable to open %s for writing.", filename);
#ifdef EMFILE
	if (errno == EMFILE)
	  marx_error (" (Too many open files.)");
#endif
	return NULL;
     }

   memcpy (header, (char *) MagicChars, 4);
   header[4] = dft->type;

   h = header + 5;
   hmax = header + HEADER_SIZE;

   p = dft->colname;
   while (h < hmax)
     {
	if (*p != 0) *h++ = *p++;
	else *h++ = 0;
     }

   if ((NUM_ROWS_OFFSET != fwrite (header, 1, NUM_ROWS_OFFSET, fp))
       || (1 != JDMwrite_int32 (&dft->num_rows, 1, fp))
       || (1 != JDMwrite_int32 (&dft->num_cols, 1, fp)))
     {
	fclose (fp);
	marx_error ("Write to %s failed.", filename);
	return NULL;
     }

   /* Now write remaining bytes. */
   h = header + RESERVED_OFFSET;
   len = HEADER_RESERVED;

   if (len != fwrite (h, 1, len, fp))
     {
	marx_error ("Write to %s failed.", filename);
	fclose (fp);
	return NULL;
     }

   dft->fp = fp;
   return fp;
}

typedef struct _Outfile_Info_Type
{
   unsigned long mask;
   char *filename;
   char *colname;
   char data_type;
   int (*write_func)(struct _Outfile_Info_Type *, FILE *fp, Marx_Photon_Attr_Type *at, double);
}
Outfile_Info_Type;

#define MAKE_WRITE_FLOAT_FUNC(name, expr) \
   static int name (Outfile_Info_Type *info, FILE *fp, Marx_Photon_Attr_Type *at, double total_time) \
   { \
      float f = expr; \
      (void) info; (void) total_time; \
      if (1 != JDMwrite_float32 (&f, 1, fp)) \
	return -1; \
      return 0; \
   }

#define MAKE_WRITE_INT32_FUNC(name, expr) \
   static int name (Outfile_Info_Type *info, FILE *fp, Marx_Photon_Attr_Type *at, double total_time) \
   { \
      int32 i = expr; \
      (void) info; (void) total_time; \
      if (1 != JDMwrite_int32 (&i, 1, fp)) \
	return -1; \
      return 0; \
   }

#define MAKE_WRITE_INT16_FUNC(name, expr) \
   static int name (Outfile_Info_Type *info, FILE *fp, Marx_Photon_Attr_Type *at, double total_time) \
   { \
      int16 i = expr; \
      (void) info; (void) total_time; \
      if (1 != JDMwrite_int16 (&i, 1, fp)) \
	return -1; \
      return 0; \
   }

#define MAKE_WRITE_INT8_FUNC(name, expr) \
   static int name (Outfile_Info_Type *info, FILE *fp, Marx_Photon_Attr_Type *at, double total_time) \
   { \
      SIGNED_CHAR i = expr; \
      (void) info; (void) total_time; \
      if (1 != fwrite(&i, 1, 1, fp)) \
	return -1; \
      return 0; \
   }

MAKE_WRITE_FLOAT_FUNC(write_pi, at->pi)
MAKE_WRITE_FLOAT_FUNC(write_energy, at->energy)
MAKE_WRITE_FLOAT_FUNC(write_time, at->arrival_time + total_time)
MAKE_WRITE_FLOAT_FUNC(write_x_x, at->x.x)
MAKE_WRITE_FLOAT_FUNC(write_x_y, at->x.y)
MAKE_WRITE_FLOAT_FUNC(write_x_z, at->x.z)
MAKE_WRITE_FLOAT_FUNC(write_p_x, at->p.x)
MAKE_WRITE_FLOAT_FUNC(write_p_y, at->p.y)
MAKE_WRITE_FLOAT_FUNC(write_p_z, at->p.z)
MAKE_WRITE_FLOAT_FUNC(write_ypixel, at->y_pixel)
MAKE_WRITE_FLOAT_FUNC(write_zpixel, at->z_pixel)
MAKE_WRITE_FLOAT_FUNC(write_upixel, at->u_pixel)
MAKE_WRITE_FLOAT_FUNC(write_vpixel, at->v_pixel)
#if MARX_HAS_DITHER
MAKE_WRITE_FLOAT_FUNC(write_sky_ra, at->dither_state.ra)
MAKE_WRITE_FLOAT_FUNC(write_sky_dec, at->dither_state.dec)
MAKE_WRITE_FLOAT_FUNC(write_sky_roll, at->dither_state.roll)
MAKE_WRITE_FLOAT_FUNC(write_det_dy, at->dither_state.dy)
MAKE_WRITE_FLOAT_FUNC(write_det_dz, at->dither_state.dz)
MAKE_WRITE_FLOAT_FUNC(write_det_theta, at->dither_state.dtheta)
#endif
MAKE_WRITE_INT16_FUNC(write_pha, at->pulse_height)
MAKE_WRITE_INT8_FUNC(write_det_num, at->ccd_num)
MAKE_WRITE_INT16_FUNC(write_mirror_shell, at->mirror_shell)
MAKE_WRITE_INT8_FUNC(write_det_region, at->detector_region)
MAKE_WRITE_INT8_FUNC(write_order, at->order)
MAKE_WRITE_INT8_FUNC(write_order1, at->support_orders[0])
MAKE_WRITE_INT8_FUNC(write_order2, at->support_orders[1])
MAKE_WRITE_INT8_FUNC(write_order3, at->support_orders[2])
MAKE_WRITE_INT8_FUNC(write_order4, at->support_orders[3])
// AML October 22, 2018
MAKE_WRITE_INT16_FUNC(write_grade, at->grade)

MAKE_WRITE_INT32_FUNC(write_photon_tag, at->tag)

Outfile_Info_Type Outfile_Info_Table[] =
{
   {MARX_PI_OK, "b_energy.dat", "B_ENERGY", 'E', &write_pi},
   {MARX_ENERGY_OK, "energy.dat", "ENERGY", 'E', &write_energy},
   {MARX_TIME_OK, "time.dat", "TIME", 'E', &write_time},
   {MARX_TAG_OK, "tag.dat", "TAG", 'J', &write_photon_tag},
   {MARX_X_VECTOR_OK, "xpos.dat", "XPOS", 'E', &write_x_x},
   {MARX_X_VECTOR_OK, "ypos.dat", "YPOS", 'E', &write_x_y},
   {MARX_X_VECTOR_OK, "zpos.dat", "ZPOS", 'E', &write_x_z},
   {MARX_P_VECTOR_OK, "xcos.dat", "COSX", 'E', &write_p_x},
   {MARX_P_VECTOR_OK, "ycos.dat", "COSY", 'E', &write_p_y},
   {MARX_P_VECTOR_OK, "zcos.dat", "COSZ", 'E', &write_p_z},
   {MARX_PULSEHEIGHT_OK, "pha.dat", "PHA", 'I', &write_pha},
   {MARX_DET_NUM_OK, "detector.dat", "CCDID", 'A', &write_det_num},
   {MARX_DET_PIXEL_OK, "xpixel.dat", "CHIPX", 'E', &write_ypixel},
   {MARX_DET_PIXEL_OK, "ypixel.dat", "CHIPY", 'E', &write_zpixel},
   {MARX_DET_UV_PIXEL_OK, "hrc_u.dat", "U", 'E', &write_upixel},
   {MARX_DET_UV_PIXEL_OK, "hrc_v.dat", "V", 'E', &write_vpixel},
   {MARX_MIRROR_SHELL_OK, "mirror.dat", "MIRROR", 'I', &write_mirror_shell},
   {MARX_DET_REGION_OK, "hrcregion.dat", "HRCREGION", 'A', &write_det_region},
   {MARX_ORDER_OK, "order.dat", "ORDER", 'A', &write_order},
   {MARX_ORDER1_OK, "ofine.dat", "FINE", 'A', &write_order1},
   {MARX_ORDER2_OK, "ocoarse1.dat", "COARSE1", 'A', &write_order2},
   {MARX_ORDER3_OK, "ocoarse2.dat", "COARSE2", 'A', &write_order3},
   {MARX_ORDER4_OK, "ocoarse3.dat", "COARSE3", 'A', &write_order4},
   {MARX_DET_GRADE_OK, "grade.dat", "GRADE", 'I', &write_grade}, // AML October 22, 2018
#if MARX_HAS_DITHER
   {MARX_SKY_DITHER_OK, "sky_ra.dat", "RA", 'E', &write_sky_ra},
   {MARX_SKY_DITHER_OK, "sky_dec.dat", "DEC", 'E', &write_sky_dec},
   {MARX_SKY_DITHER_OK, "sky_roll.dat", "ROLL", 'E', &write_sky_roll},
   {MARX_DET_DITHER_OK, "det_dy.dat", "DET_DY", 'E', &write_det_dy},
   {MARX_DET_DITHER_OK, "det_dz.dat", "DET_DZ", 'E', &write_det_dz},
   {MARX_DET_DITHER_OK, "det_theta.dat", "DET_THETA", 'E', &write_det_theta},
#endif
};
#define MAX_NUM_FILES (sizeof(Outfile_Info_Table)/sizeof(Outfile_Info_Type))

static int open_files (char *dir, unsigned long write_mask, int new_file) /*{{{*/
{
   char *filename = NULL;
   char *open_mode;
   unsigned int nbits;

   if (MAX_NUM_FILES > NBITS_LONG)
     {
	marx_error ("%s", "marx internal error: The number of table entries must be less than NBITS_LONG");
	return -1;
     }

   open_mode = "r+b";
   /* Note: "ab" will not work because of the system 5 behavior of writes
    * to a file opened this way: all writes will go to end of file
    * regardless of seek!!!
    */

   for (nbits = 0; nbits < MAX_NUM_FILES; nbits++)
     {
	FILE *fp;
	Outfile_Info_Type *info;

	info = Outfile_Info_Table + nbits;
	if ((info->mask & write_mask) == 0)
	  continue;

	if (NULL == (filename = marx_dircat (dir, info->filename)))
	  goto return_error;

	if (new_file)
	  {
	     Marx_Dump_File_Type dft;

	     memset ((char *) &dft, 0, sizeof (Marx_Dump_File_Type));
	     dft.type = info->data_type;
	     strncpy (dft.colname, info->colname, sizeof(dft.colname));
	     dft.colname[sizeof(dft.colname)-1] = 0;
	     fp = marx_create_write_dump_file (filename, &dft);
	  }
	else fp = fopen (filename, open_mode);

	if (fp == NULL)
	  {
	     marx_error ("Unable to open %s.", filename);
#ifdef EMFILE
	     if (errno == EMFILE)
	       marx_error (" (Too many open files.)");
#endif
	     goto return_error;
	  }

	File_Pointers[nbits].fp = fp;

	if (new_file == 0)
	  {
	     /* See comment above regarding r+t mode. */
	     if (-1 == FSEEK (fp, 0, SEEK_END))
	       {
		  marx_error ("seek error.");
		  goto return_error;
	       }
	  }
	SLFREE (filename);
	filename = NULL;
     }

   return 0;

   return_error:
   if (filename != NULL) SLFREE (filename);
   (void) close_write_files (0);
   return -1;
}

/*}}}*/

int marx_write_photons (Marx_Photon_Type *p, unsigned long write_mask, /*{{{*/
			char *dir, int open_mode, double total_time)
{
   unsigned int i, imax, nbits;
   Marx_Photon_Attr_Type *attr, *at;

#if MARX_HAS_DITHER
   if (_Marx_Dither_Mode != _MARX_DITHER_MODE_NONE)
     write_mask &= (p->history | MARX_SKY_DITHER_OK | MARX_DET_DITHER_OK);
   else
#endif
     write_mask &= p->history;

   if (-1 == open_files (dir, write_mask, open_mode))
     return -1;

   imax = p->n_photons;
   attr = p->attributes;

   for (i = 0; i < imax; i++)
     {
#if MARX_HAS_DITHER
	double ra = 0, dec = 0;
#endif
	at = attr + i;
	if (at->flags & BAD_PHOTON_MASK)
	  {
#define OUTPUT_ABSORBED	0
#if OUTPUT_ABSORBED
	     if (0 == (at->flags & OUTPUT_ABSORBED))
#endif
	       continue;
	  }
#if OUTPUT_ABSORBED
	else continue;
#endif

#if MARX_HAS_DITHER
	if (write_mask & MARX_SKY_DITHER_OK)
	  _marx_ray_to_sky_ra_dec (at, &ra, &dec);
#endif

	for (nbits = 0; nbits < MAX_NUM_FILES; nbits++)
	  {
	     FILE *fp;
	     Outfile_Info_Type *info;
	     unsigned long mask;

	     info = Outfile_Info_Table + nbits;
	     mask = info->mask;
	     if ((mask & write_mask) == 0)
	       continue;

	     fp = File_Pointers[nbits].fp;
	     if (fp == NULL)
	       continue;

	     if (-1 == info->write_func (info, fp, at, total_time))
	       {
		  marx_error ("write error.");
		  (void) close_write_files (0);
		  return -1;
	       }
	     File_Pointers[nbits].num_rows += 1;
	  }
     }

   if (-1 == close_write_files (1))
     {
	marx_error ("write error.");
	return -1;
     }
   return 0;
}

/*}}}*/

static void dump_usage (void) /*{{{*/
{
   fprintf (stderr, "Usage: marx --dump filename [filename...]\n");
}

/*}}}*/

static Marx_Dump_File_Type *Data_File_Type_Root;

static int chain_read_data_file_type (Marx_Dump_File_Type *d) /*{{{*/
{
   Marx_Dump_File_Type *dlast;

   if (Data_File_Type_Root == NULL)
     {
	Data_File_Type_Root = d;
	return 0;
     }

   dlast = Data_File_Type_Root;
   while (dlast->next != NULL) dlast = dlast->next;
   dlast->next = d;
   return 0;
}

/*}}}*/

static int close_read_data_files (void) /*{{{*/
{
   Marx_Dump_File_Type *d, *dnext;
   int ret = 0;

   d = Data_File_Type_Root;

   while (d != NULL)
     {
	dnext = d->next;
	marx_close_read_dump_file (d);
	d = dnext;
     }

   Data_File_Type_Root = NULL;
   return ret;
}

/*}}}*/

static int dump_data_files (void) /*{{{*/
{
   int not_done;
   Marx_Dump_File_Type *dft;
   int32 num_rows;
   unsigned int num_cols;

   dft = Data_File_Type_Root;

   fputc ('#', stdout);
   while (dft != NULL)
     {
	if (0 == (num_cols = dft->num_cols))
	  num_cols = 1;

	while (num_cols--)
	  fprintf (stdout, "%16s", dft->colname);

	dft = dft->next;
     }
   fputc ('\n', stdout);

   num_rows = 0;
   while (1)
     {
	not_done = 1;

	dft = Data_File_Type_Root;

	while (dft != NULL)
	  {
	     FILE *fp = dft->fp;
	     float64 dval;
	     float32 fval;
	     int32 ival;
	     int16 sval;
	     SIGNED_CHAR cval;

	     num_cols = dft->num_cols;
	     if (num_cols == 0)
	       num_cols = 1;

	     switch (dft->type)
	       {
		default:
		  marx_error ("Unknown data type: '%c'", dft->type);
		  (void) close_read_data_files ();
		  return -1;

		case 'A':
		  while (num_cols && not_done)
		    {
		       if (1 == fread (&cval, 1, 1, fp))
			 fprintf (stdout, "%16d", cval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;

		case 'I':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_int16 (&sval, 1, fp))
			 fprintf (stdout, "%16d", sval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;

		case 'J':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_int32 (&ival, 1, fp))
			 fprintf (stdout, "%16ld", (long) ival);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;

		case 'E':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_float32 (&fval, 1, fp))
			 fprintf (stdout, "%16e", fval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;

		case 'D':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_float64 (&dval, 1, fp))
			 fprintf (stdout, "%16e", dval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;
	       }

	     if (not_done == 0)
	       {
		  if (num_rows != dft->num_rows)
		    {
		       marx_error ("Read error. (%d/%d rows)",
				  num_rows, dft->num_rows);
		       close_read_data_files ();
		       return -1;
		    }
		  break;
	       }
	     dft = dft->next;
	  }
	num_rows++;
	putc ('\n', stdout);
	if (not_done == 0)
	  break;
     }

   close_read_data_files ();
   return 0;
}

/*}}}*/

Marx_Dump_File_Type *marx_open_read_dump_file (char *file) /*{{{*/
{
   FILE *fp;
   char header[HEADER_SIZE];
   char reserved[HEADER_RESERVED];
   Marx_Dump_File_Type *dft;

   dft = (Marx_Dump_File_Type *) SLMALLOC (sizeof (Marx_Dump_File_Type));
   if (dft == NULL)
     {
	marx_error ("Memory allocation failure.");
	return NULL;
     }

   memset ((char *) dft, 0, sizeof (Marx_Dump_File_Type));

   if (NULL == (fp = fopen (file, "rb")))
     {
	marx_error ("Unable to open %s.", file);
#ifdef EMFILE
	if (errno == EMFILE)
	  marx_error (" (Too many open files.)");
#endif
	SLFREE (dft);
	return NULL;
     }

   if ((NUM_ROWS_OFFSET != fread (header, 1, NUM_ROWS_OFFSET, fp))
       || (1 != JDMread_int32 (&dft->num_rows, 1, fp))
       || (1 != JDMread_int32 (&dft->num_cols, 1, fp))
       || (HEADER_RESERVED != fread(reserved, 1, HEADER_RESERVED, fp)))
     {
	marx_error ("Error reading header information for %s.", file);
	fclose (fp);
	SLFREE (dft);
	return NULL;
     }

   if (memcmp (header, (char *)MagicChars, 4))
     {
	marx_error ("File %s does not appear to be a marx output file.  Bad Magic Number.",
		   file);
	fclose (fp);
	SLFREE (dft);
	return NULL;
     }

   dft->type = *(header + 4);
   strncpy (dft->colname, header + 5, 15);
   dft->colname[15] = 0;
   dft->fp = fp;

   return dft;
}

/*}}}*/

/* Some programs may pass NULL to this.  No error please. */
void marx_close_read_dump_file (Marx_Dump_File_Type *dft) /*{{{*/
{
   if (dft == NULL)
     return;

   if (dft->fp != NULL)
     (void) fclose (dft->fp);

   SLFREE (dft);
}

/*}}}*/

int marx_dump (int argc, char **argv) /*{{{*/
{
   char *file;
   int i;
   Marx_Dump_File_Type *dft;

   if (argc < 3)
     {
	dump_usage ();
	return -1;
     }

   for (i = 2; i < argc; i++)
     {
	file = argv[i];

	dft = marx_open_read_dump_file (file);

	if ((dft == NULL)
	    || (-1 == chain_read_data_file_type (dft)))

	  {
	     marx_error ("Unable to open %s.", file);
#ifdef EMFILE
	     if (errno == EMFILE)
	       marx_error (" (Too many open files.)");
#endif
	     close_read_data_files ();
	     if (dft != NULL) fclose (dft->fp);
	     return -1;
	  }
     }

   return dump_data_files ();
}

/*}}}*/
