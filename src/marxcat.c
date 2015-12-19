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
#include <string.h>
#if HAVE_STDLIB_H
# include <stdlib.h>
#endif
/* These are needed for mkdir */
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <errno.h>

#include <marx.h>

#define NUM_DATA_FILE_NAMES	27
static char *Data_Files [NUM_DATA_FILE_NAMES] =
{
   "time.dat",			       /* required */
   "energy.dat",
   "xpos.dat",
   "ypos.dat",
   "zpos.dat",
   "xcos.dat",
   "ycos.dat",
   "zcos.dat",
   "pha.dat",
   "detector.dat",
   "xpixel.dat",
   "ypixel.dat",
   "mirror.dat",
   "order.dat",
   "ofine.dat",
   "ocoarse1.dat",
   "ocoarse2.dat",
   "ocoarse3.dat",
   "sky_ra.dat",
   "sky_dec.dat",
   "sky_roll.dat",
   "det_dy.dat",
   "det_dz.dat",
   "det_theta.dat",
   "hrc_u.dat",
   "hrc_v.dat",
   "b_energy.dat"
};

typedef struct
{
   char *filename;
   Marx_Dump_File_Type *mdft;
   FILE *fp;			       /* pointer to fp in mdft struct */
   unsigned int data_size;
}
Data_File_Type;

typedef struct
{
   char *dir;
   Data_File_Type dft_list [NUM_DATA_FILE_NAMES];
   float32 current_time;
   unsigned int flags;
   FILE *fp_time;			       /* pointer to time.dat fp */
   unsigned long num_read;
}
Data_Dir_Type;

static void free_data_file (Data_File_Type *dft)
{
   if (dft == NULL)
     return;

   if (dft->mdft != NULL)
     marx_close_read_dump_file (dft->mdft);

   marx_free (dft->filename);

   dft->mdft = NULL;
   dft->filename = NULL;
}

static void close_read_data_files (Data_File_Type *dft_list)
{
   unsigned int i;

   if (dft_list == NULL) return;

   for (i = 0; i < NUM_DATA_FILE_NAMES; i++)
     {
	free_data_file (dft_list + i);
     }
}

static void free_data_dir_type (Data_Dir_Type *ddt)
{
   if (ddt == NULL) return;
   marx_free (ddt->dir);
   close_read_data_files (ddt->dft_list);
   marx_free ((char *) ddt);
}

static Data_Dir_Type *allocate_data_dir_type (char *dir)
{
   Data_Dir_Type *ddt;
   unsigned int len;

   ddt = (Data_Dir_Type *) marx_malloc (sizeof (Data_Dir_Type));
   if (ddt == NULL)
     return NULL;
   memset ((char *) ddt, 0, sizeof(Data_Dir_Type));

   len = strlen (dir);
   if (NULL == (ddt->dir = (char *)marx_malloc (len + 1)))
     {
	free_data_dir_type (ddt);
	return NULL;
     }
   strcpy (ddt->dir, dir);

   return ddt;
}

static Data_Dir_Type *open_read_data_dir (char *dir)
{
   Data_Dir_Type *ddt;
   Data_File_Type *dft_list;
   unsigned int i;

   ddt = allocate_data_dir_type (dir);
   if (ddt == NULL)
     return NULL;

   dft_list = ddt->dft_list;

   marx_message ("Examining files in %s\n", dir);

   for (i = 0; i < NUM_DATA_FILE_NAMES; i++)
     {
	char *filename;
	unsigned int data_size;

	filename = marx_dircat (dir, Data_Files [i]);
	if (filename == NULL)
	  {
	     free_data_dir_type (ddt);
	     return NULL;
	  }

	if (0 == marx_file_exists (filename))
	  {
	     dft_list[i].mdft = NULL;
	     dft_list[i].filename = NULL;
	     marx_free (filename);
	     continue;
	  }

	dft_list[i].filename = filename;
	if (NULL == (dft_list[i].mdft = marx_open_read_dump_file (filename)))
	  {
	     free_data_dir_type (ddt);
	     return NULL;
	  }

	dft_list[i].fp = dft_list[i].mdft->fp;

	data_size = 0;
	switch (dft_list[i].mdft->type)
	  {
	   case 'A':
	     data_size = sizeof (char);
	     break;
	   case 'I':
	     data_size = sizeof (int16);
	     break;
	   case 'J':
	     data_size = sizeof (int32);
	     break;
	   case 'E':
	     data_size = sizeof (float32);
	     break;
	   case 'D':
	     data_size = sizeof (float64);
	     break;
	   default:
	     fprintf (stderr, "%s has type '%c'; it is not supported.\n",
		      filename, dft_list[i].mdft->type);
	     free_data_dir_type (ddt);
	     return NULL;
	  }

	dft_list[i].data_size = data_size;
     }

   return ddt;
}

static void free_directory_list (Data_Dir_Type **dir_list, unsigned int ndirs)
{
   unsigned int i;
   Data_Dir_Type *ddt;

   if (dir_list == NULL) return;

   for (i = 0; i < ndirs; i++)
     {
	if (NULL == (ddt = dir_list [i]))
	  continue;

	free_data_dir_type (ddt);
	dir_list [i] = NULL;
     }

   marx_free ((char *) dir_list);
}

static Data_Dir_Type **open_read_directories (char **dirs, unsigned int ndirs,
					      struct stat *new_st)
{
   unsigned int i;
   Data_Dir_Type **dir_list;

   /* Make sure the directories already exist. */
   for (i = 0; i < ndirs; i++)
     {
	struct stat st;
	char *dir = dirs[i];

	if (2 != marx_file_exists (dir))
	  {
	     marx_error ("%s is not a directory.\n", dir);
	     return NULL;
	  }
	if (-1 == stat (dir, &st))
	  {
	     marx_error ("Unable to stat %s", dir);
	     return NULL;
	  }

	if ((st.st_dev == new_st->st_dev)
	    && (st.st_ino == new_st->st_ino))
	  {
	     marx_error ("You cannot marxcat %s directory onto itself", dir);
	     return NULL;
	  }
     }

   if (NULL == (dir_list = (Data_Dir_Type **) marx_calloc (ndirs, sizeof (Data_Dir_Type *))))
     return NULL;

   for (i = 0; i < ndirs; i++)
     {
	if (NULL == (dir_list[i] = open_read_data_dir (dirs[i])))
	  {
	     free_directory_list (dir_list, ndirs);
	     return NULL;
	  }
     }
   return dir_list;
}

static int check_dirlist_consistency (Data_Dir_Type **dir_list, unsigned int ndirs)
{
   unsigned int i, j;

   for (i = 0; i < NUM_DATA_FILE_NAMES; i++)
     {
	unsigned int num_present;
	Data_Dir_Type *ddt;
	char *file;

	num_present = 0;

	for (j = 0; j < ndirs; j++)
	  {
	     ddt = dir_list[j];
	     if (ddt->dft_list[i].mdft != NULL) num_present++;
	  }

	if (num_present == ndirs)
	  continue;

	file = Data_Files [i];
	if (i == 0)
	  {
	     fprintf (stderr, "%s must exist in all directories.\n", file);
	     return -1;
	  }

	if (num_present == 0)
	  continue;

	fprintf (stderr, "%s does not exist in all directories.  Skipping it.\n",
		 file);

	for (j = 0; j < ndirs; j++)
	  {
	     if (NULL == (ddt = dir_list [j]))
	       continue;

	     free_data_file (ddt->dft_list + i);
	  }
     }

   /* Now loop through the file list and compress those. */

   return 0;
}

static int close_output_files (Marx_Dump_File_Type *mdft, char **files)
{
   int ret;
   unsigned int i;

   ret = 0;
   for (i = 0; i < NUM_DATA_FILE_NAMES; i++)
     {
	marx_free (files [i]);
	if (mdft[i].fp != NULL)
	  {
	     int status;

	     status = marx_close_write_dump_file (mdft[i].fp, mdft[i].num_rows);
	     if (status == -1)
	       ret = status;
	  }
	mdft[i].fp = NULL;
	files[i] = NULL;
     }

   return ret;
}

static int check_read_error (Data_Dir_Type *ddt)
{
   if (ddt->num_read != (unsigned long) ddt->dft_list[0].mdft->num_rows)
     {
	fprintf (stderr, "Read error processing %s in %s\n",
		 Data_Files[0], ddt->dir);
	return -1;
     }

   return 0;
}

static unsigned int compress_dir_list (Data_Dir_Type **dirs, unsigned int ndirs)
{
   Data_Dir_Type **dir0, **dir1, **dirmax;

   dir0 = dir1 = dirs;
   dirmax = dirs + ndirs;

   while (dir1 < dirmax)
     {
	if (*dir1 == NULL)
	  {
	     ndirs--;
	     dir1++;
	     continue;
	  }
	if (dir1 != dir0)
	  {
	     *dir0 = *dir1;
	     *dir1 = NULL;
	  }

	dir0++;
	dir1++;
     }

   return ndirs;
}

static int write_row_of_data (Data_Dir_Type *dir,
			      Marx_Dump_File_Type *new_mdfts, char **filenames,
			      unsigned int num_files)
{
   Data_File_Type *dft;
   unsigned int i;

   if (1 != JDMwrite_float32 (&dir->current_time, 1, new_mdfts [0].fp))
     {
	fprintf (stderr, "Write error: %s\n", filenames [0]);
	return -1;
     }
   new_mdfts[0].num_rows += 1;

   for (i = 1; i < num_files; i++)
     {
	unsigned char buf[32];
	unsigned int data_size;

	dft = dir->dft_list + i;

	data_size = dft->data_size;
	if (1 != fread ((char *) buf, data_size, 1, dft->fp))
	  {
	     fprintf (stderr, "Read Error: %s\n", dft->filename);
	     return -1;
	  }

	if (1 != fwrite ((char *) buf, data_size, 1, new_mdfts [i].fp))
	  {
	     fprintf (stderr, "Write error: %s\n", filenames [i]);
	     return -1;
	  }
	new_mdfts[i].num_rows += 1;
     }
   return 0;
}

static int do_marxcat (char *odir, Data_Dir_Type **dirs, unsigned int ndirs)
{
   Marx_Dump_File_Type new_mdfts [NUM_DATA_FILE_NAMES];
   char *filenames [NUM_DATA_FILE_NAMES];
   unsigned int i, j;
   unsigned int num_files;

   memset ((char *) new_mdfts, 0, sizeof (new_mdfts));
   memset ((char *) filenames, 0, sizeof (filenames));

   num_files = 0;
   for (i = 0; i < NUM_DATA_FILE_NAMES; i++)
     {
	char *file;
	Marx_Dump_File_Type *mdftp;

	if (NULL == (mdftp = dirs[0]->dft_list[i].mdft))
	  continue;

	new_mdfts[num_files].type = mdftp->type;
	strcpy (new_mdfts[num_files].colname, mdftp->colname);

	file = marx_dircat (odir, Data_Files[i]);
	if (file == NULL)
	  {
	     close_output_files (new_mdfts, filenames);
	     return -1;
	  }

	marx_message ("Creating %s\n", file);
	if (NULL == marx_create_write_dump_file (file, new_mdfts + num_files))
	  {
	     close_output_files (new_mdfts, filenames);
	     return -1;
	  }

	filenames [num_files] = file;
	if (num_files != i)
	  {
	     /* Squeeze the dat_file_type for each directory also. */
	     for (j = 0; j < ndirs; j++)
	       {
		  dirs[j]->dft_list [num_files] = dirs[j]->dft_list [i];
		  memset ((char *) (dirs [j]->dft_list + i),
			  0, sizeof (Data_File_Type));
	       }
	  }
	num_files++;
     }

   /* Now start the real work */

   j = 0;
   while (j < ndirs)
     {
	dirs[j]->fp_time = dirs[j]->dft_list[0].mdft->fp;

	if (1 == JDMread_float32 (&(dirs[j]->current_time), 1, dirs[j]->fp_time))
	  {
	     dirs[j]->num_read = 1;
	     j++;
	     continue;
	  }

	if (-1 == check_read_error (dirs[j]))
	  {
	     close_output_files (new_mdfts, filenames);
	     return -1;
	  }
	free_data_dir_type (dirs[j]);
	dirs [j] = NULL;
	ndirs = compress_dir_list (dirs, ndirs);
     }

   while (ndirs)
     {
	unsigned int min_j = 0;
	double tmin = dirs [0]->current_time;

	for (j = 1; j < ndirs; j++)
	  {
	     if (dirs [j]->current_time < tmin)
	       {
		  min_j = j;
		  tmin = dirs [j]->current_time;
	       }
	  }

	if (-1 == write_row_of_data (dirs[min_j], new_mdfts, filenames, num_files))
	  {
	     close_output_files (new_mdfts, filenames);
	     return -1;
	  }

	if (1 == JDMread_float32 (&(dirs[min_j]->current_time), 1, dirs[min_j]->fp_time))
	  {
	     dirs[min_j]->num_read++;
	     continue;
	  }

	if (-1 == check_read_error (dirs[min_j]))
	  {
	     close_output_files (new_mdfts, filenames);
	     return -1;
	  }
	free_data_dir_type (dirs[min_j]);
	dirs [min_j] = NULL;
	ndirs = compress_dir_list (dirs, ndirs);
     }

   return close_output_files (new_mdfts, filenames);
}

static int cp_files (FILE *fpin, FILE *fpout)
{
   char buf[1024];
   unsigned int readlen;
   do
     {
	readlen = fread (buf, 1, sizeof(buf), fpin);
	if (readlen)
	  {
	     if (readlen != fwrite (buf, 1, readlen, fpout))
	       return -1;
	  }
     }
   while (readlen == sizeof (buf));

   return 0;
}

static int cp_to_new_file (char *name,
			   char *odir,
			   Data_Dir_Type **dirlist, unsigned int ndirs)
{
   char *ofile, *ifile;
   FILE *fpin, *fpout;
   unsigned int i;

   if (NULL == (ofile = marx_dircat (odir, name)))
     return -1;

   if (NULL == (ifile = marx_dircat (dirlist[0]->dir, name)))
     {
	marx_free (ofile);
	return -1;
     }

   if (NULL == (fpin = fopen (ifile, "rb")))
     {
	marx_error ("Unable to open file %s for reading.", ifile);
	goto return_error;
     }

   if (NULL == (fpout = fopen (ofile, "wb")))
     {
	fclose (fpin);
	marx_error ("Unable to open %s for writing.", ofile);
	goto return_error;
     }

   if (EOF == fputs ("\
#@#The simulation represented by this directory is a composite simulation.\n\
#@#The following directories were used:\n",
		       fpout))
     goto return_write_error;

   for (i = 0; i < ndirs; i++)
     {
	if (EOF == fprintf (fpout, "#@# %s\n", dirlist[i]->dir))
	  goto return_write_error;
     }

   if (EOF == fputs ("\
#@#The rest of this file is a copy of the file from the\n\
#@#first directory in the above list.\n",
		     fpout))
     goto return_write_error;

   if ((-1 == cp_files (fpin, fpout))
       || (EOF == fclose (fpout)))
     goto return_write_error;

   (void) fclose (fpin);
   marx_free (ofile);
   marx_free (ifile);
   return 0;

   /* Get here only if error occurs. */
   return_write_error:
   marx_error ("Write to %s failed.", ofile);
   (void) fclose (fpout);
   (void) fclose (fpin);

   /* drop */
   return_error:
   marx_free (ofile);
   marx_free (ifile);
   return -1;
}

static char *Pgm_Name;
static void usage (void)
{
   fprintf (stderr, "marxcat version: %s\n", MARX_VERSION);
   fprintf (stderr, "Usage: %s [--help] DIR1 [DIR2 [DIR3 ...]] NEWDIR\n", Pgm_Name);
   exit (1);
}

int main (int argc, char **argv)
{
   Data_Dir_Type **dir_list;
   unsigned int ndirs;
   char *output_dir;
   int status;
   struct stat st;

   Pgm_Name = argv[0];

   if (argc < 3) usage ();

   if (!strcmp (argv[1], "--help"))
     usage ();

   ndirs = (unsigned int)argc - 2;
   output_dir = argv [argc - 1];
   if ((-1 == mkdir (output_dir, 0777))
       && (errno != EEXIST))
     {
	marx_error ("Error creating directory %s", output_dir);
	return 1;
     }
   if (-1 == stat (output_dir, &st))
     {
	marx_error ("stat failed on directory %s", output_dir);
	return 1;
     }

   if (NULL == (dir_list = open_read_directories (argv + 1, ndirs, &st)))
     return 1;

   if (-1 == check_dirlist_consistency (dir_list, ndirs))
     {
	free_directory_list (dir_list, ndirs);
	return 1;
     }

   if (-1 == cp_to_new_file ("marx.par", output_dir, dir_list, ndirs))
     {
	free_directory_list (dir_list, ndirs);
	return 1;
     }

   if (-1 == cp_to_new_file ("obs.par", output_dir, dir_list, ndirs))
     {
	free_directory_list (dir_list, ndirs);
	return 1;
     }

   fprintf (stderr, "Merging directories...\n");
   status = do_marxcat (output_dir, dir_list, ndirs);
   fputs ("\n", stderr);

   free_directory_list (dir_list, ndirs);

   if (status == -1)
     return 1;
   return 0;
}
