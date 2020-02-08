/* -*- mode: C; mode: fold; -*- */
#include <stdio.h>

#include <stdlib.h>
#include <string.h>

#include "cfits.h"
#include <errno.h>

#include <marx.h>
#include <jdmath.h>
#include <pfile.h>

#include "argcargv.h"

#define MAX_FILENAME_LEN 1024

static void integrate_response_vector (CFits_Response_Vector_Type *v) /*{{{*/
{
   CFits_Response_Element_Type *elem, *elem_max;
   double sum;
   
   elem = v->elem;
   if (elem == NULL) return;
   
   elem_max = elem + v->num_grps;
   sum = 0.0;
   while (elem < elem_max)
     {
	double *rsp, *rsp_max;
	
	rsp = elem->response;
	rsp_max = rsp + elem->num_channels;
	
	while (rsp < rsp_max)
	  {
	     sum += *rsp;
	     *rsp = sum;
	     
	     rsp++;
	  }
	elem++;
     }
}

/*}}}*/

static int compute_response_from_vector (CFits_Response_Vector_Type *v) /*{{{*/
{
   double r;
   CFits_Response_Element_Type *elem, *elem_max;
   
   r = JDMrandom ();
   
   elem = v->elem;
   if (elem == NULL) return -1;
   
   elem_max = elem + v->num_grps;
   while (elem < elem_max)
     {
	double *rsp, *rsp_max;
	
	rsp = elem->response;
	rsp_max = rsp + elem->num_channels;
	
	while (rsp < rsp_max)
	  {
	     if (r <= *rsp)
	       {
		  return elem->first_channel 
		    +  (unsigned int)(rsp - elem->response);
	       }
	     rsp++;
	  }
	
	elem++;
     }
   
   return -1;
}
   
/*}}}*/

static int compute_response (CFits_Response_File_Type *rsp, /*{{{*/
			     float *energies, unsigned int *sort_index,
			     char *filter, int *channels, unsigned int num_energies)
{
   float emin, emax, en;
   long row, num_rows;
   unsigned int i;
   CFits_Type *ft;
   
   num_rows = rsp->num_rows;
   ft = rsp->ft;
   
   
   /* Find out which energy bin to use. */
   row = 0;
   emax = -1.0e30;

   i = 0;
   while (i < num_energies)
     {
	CFits_Response_Vector_Type rv;
	unsigned int j;
	
	j = sort_index[i];
	i++;

	if ((filter != NULL) && (filter[j] == 0))
	  continue;

	en = energies[j];
	
	while ((row < num_rows) && (en >= emax))
	  {
	     row++;
	     if ((-1 == cfits_read_column_floats (ft, rsp->energ_lo_col, row, 1, &emin, 1))
		 || (-1 == cfits_read_column_floats (ft, rsp->energ_hi_col, row, 1, &emax, 1)))
	       return -1;
	  }

	if ((en >= emax) || (en < emin))
	  {
	     channels[j] = -1;
	     continue;
	  }

	if (-1 == cfits_read_response_vector (rsp, row, &rv))
	  return -1;
	
	integrate_response_vector (&rv);
	
	while (1)
	  {
	     if ((filter == NULL) || (filter[j] != 0))
	       channels [j] = compute_response_from_vector (&rv);

	     if (i == num_energies) break;

	     j = sort_index [i];
	     en = energies[j];
	     
	     if (en >= emax)
	       break;
	     
	     i++;
	  }
	
	cfits_free_response_vector (&rv);
     }

   return 0;	
}

/*}}}*/

#define MAX_ENERGIES 8192
static float Energies[MAX_ENERGIES];
static int Channels[MAX_ENERGIES];
static char Filter [MAX_ENERGIES];
static char *Rsp_File, *Marx_Dir;
static int CCD_Id;
static float Chip_X_Min = 0;
static float Chip_X_Max = 1024;
static float Chip_Y_Min = 0;
static float Chip_Y_Max = 1024;
static int Force_Flag;

static ArgcArgv_Type ArgcArgv_Table [] = /*{{{*/
{
   {"--marx", ARGCARGV_STRING, (long) &Marx_Dir, "Name of marx output dir"},
   {"--rmf", ARGCARGV_STRING, (long) &Rsp_File, "Name of RMF file"},
   {"--chip", ARGCARGV_INTEGER, (long) &CCD_Id, "CCD chip number"},
   {"--xmin", ARGCARGV_FLOAT, (long) &Chip_X_Min, "min X pixel"},
   {"--ymin", ARGCARGV_FLOAT, (long) &Chip_Y_Min, "min Y pixel"},
   {"--xmax", ARGCARGV_FLOAT, (long) &Chip_X_Max, "max X pixel"},
   {"--ymax", ARGCARGV_FLOAT, (long) &Chip_Y_Max, "max Y pixel"},
   {"--force", ARGCARGV_BOOLEAN, (long) &Force_Flag, NULL},
   {NULL, 0, (long) NULL, NULL}
};

/*}}}*/

   
static void usage (char *pgm) /*{{{*/
{
   fprintf (stderr, "Usage:\n%s [optional args] --rmf <rmf-file> --marx <marx-dir>\n", pgm);
   fprintf (stderr, "Optional arguments include:\n");
   fprintf (stderr, "  --chip <ccdid>\n");
   fprintf (stderr, "  --xmin <min x pixel>\n");
   fprintf (stderr, "  --ymin <min y pixel>\n");
   fprintf (stderr, "  --xmax <max x pixel>\n");
   fprintf (stderr, "  --ymax <max y pixel>\n");
   fprintf (stderr, "  --force\n");

   fprintf (stderr, "\n\
The --force option will allow the use of fits files that have the HDUCLAS3\n\
keyword set to FULL.  Keep in mind that such files already have the effective\n\
folded in, or are not up to spec.\n");
	    
   fprintf (stderr, "\nA cartoon of the chip numbers for the ACIS-I and ACIS-S follows:\n\n");
   fprintf (stderr, "  ACIS-I:   01\n");
   fprintf (stderr, "            23\n");
   fprintf (stderr, "\n");
   fprintf (stderr, "  ACIS-S: 456789\n\n");
#if 1
   fprintf (stderr, "\
***Note: This program modifies the file called pha.dat and copies the original\n\
         version to pha.dat.bak\n");
#endif
}

/*}}}*/


static int check_simulation_parameters (void)
{
   char filename [MAX_FILENAME_LEN];
   char buf [PF_MAX_LINE_LEN];
   Marx_Detector_Type *det;
   int ideal;
   Param_File_Type *pf;

   sprintf (filename, "%s/%s", Marx_Dir, "marx.par");
   if (NULL == (pf = pf_open_parameter_file (filename, "r")))
     {
	fprintf (stderr, "Unable to open your MARX parameter file in %s\n",
		 Marx_Dir);
	return -1;
     }
   
   if ((-1 == pf_get_string (pf, "DetectorType", buf, sizeof (buf)))
       || (-1 == pf_get_boolean (pf, "DetIdeal", &ideal)))
     {
	pf_close_parameter_file (pf);
	return -1;
     }
   
   pf_close_parameter_file (pf);

   if ((0 != strncmp (buf, "ACIS", 4))
       || (NULL == (det = marx_get_detector_info (buf))))
     {
	fprintf (stderr, "DetectorType %s not supported by this program\n",
		 buf);
	return -1;
     }

   if ((CCD_Id != -1)
       && ((CCD_Id < det->first_facet_id)
	   || (CCD_Id > det->last_facet_id)))
     {
	fprintf (stderr, "Your `--chip' parameter is out of range for the detector\n");
	return -1;
     }

   if (ideal == 0)
     fprintf (stderr, "Warning: this simulation was not run with DetIdeal set to yes.\n");

   return 0;
}

static int 
read_energies_etc (unsigned int *num_energies, 
		   FILE *en_fp, FILE *pha_fp, FILE *det_fp, FILE *chipx_fp, FILE *chipy_fp)
{
   unsigned int i;
   unsigned char det;
   unsigned int num_read;

   *num_energies = num_read = JDMread_f_float32 (Energies, MAX_ENERGIES, en_fp);

   if (num_read == 0)
     return 0;

   memset (Filter, 1, num_read);       /* nothing gets filtered */

   if (pha_fp == NULL)
     return 0;			       /* only needed if we filter */
	
   if (num_read != JDMread_i_int16 (Channels, num_read, pha_fp))
     {
	fprintf (stderr, "*** Read of pha.dat failed\n");
	return -1;
     }

   if (det_fp != NULL) 
     {
	for (i = 0; i < num_read; i++)
	  {
	     if (1 != fread ((char *) &det, 1, 1, det_fp))
	       {
		  fprintf (stderr, "*** Read of detector.dat failed\n");
		  return -1;
	       }
	     if ((int) det != CCD_Id)
	       Filter[i] = 0;
	  }
     }
   if (chipx_fp != NULL) 
     {
	for (i = 0; i < num_read; i++)
	  {
	     float x;
	     if (1 != JDMread_f_float32 (&x, 1, chipx_fp))
	       {
		  fprintf (stderr, "*** Read of xpixel.dat failed\n");
		  return -1;
	       }
	     if ((x < Chip_X_Min) || (x > Chip_X_Max))
	       Filter[i] = 0;
	  }
     }

   if (chipy_fp != NULL) 
     {
	for (i = 0; i < num_read; i++)
	  {
	     float y;
	     if (1 != JDMread_f_float32 (&y, 1, chipy_fp))
	       {
		  fprintf (stderr, "*** Read of ypixel.dat failed\n");
		  return -1;
	       }
	     if ((y < Chip_Y_Min) || (y > Chip_Y_Max))
	       Filter[i] = 0;
	  }
     }

   return 0;
}

static Marx_Dump_File_Type *open_marx_file (char *file, int type, FILE **fp)
{
   char filename [MAX_FILENAME_LEN];
   Marx_Dump_File_Type *d;
   
   *fp = NULL;

   sprintf (filename, "%s/%s", Marx_Dir, file);
   if (NULL == (d = marx_open_read_dump_file (filename)))
     {
	fprintf (stderr, "*** Unable to open %s\n", filename);
	return NULL;
     }

   if (type != d->type)
     {
	fprintf (stderr, "%s has incorrect type, expecting '%c'\n", filename, type);
	marx_close_read_dump_file (d);
	return NULL;
     }

   *fp = d->fp;
   return d;
}

static int rename_marx_file (char *f, char *t)
{
   char from[MAX_FILENAME_LEN];
   char to[MAX_FILENAME_LEN];
   
   sprintf (from, "%s/%s", Marx_Dir, f);
   sprintf (to, "%s/%s", Marx_Dir, t);
   
   if (0 == rename (from, to))
     return 0;
   
   if (errno == EEXIST)
     {
	if ((0 == remove (to))
	    && (0 == rename (from, to)))
	  return 0;
     }

   fprintf (stderr, "Unable to rename %s to %s (errno=%d)\n", from, to, errno);
   return -1;
}

     
int main (int argc, char **argv) /*{{{*/
{
   Marx_Dump_File_Type *energy_df;
   Marx_Dump_File_Type *detector_df;
   Marx_Dump_File_Type *chipx_df;
   Marx_Dump_File_Type *chipy_df;
   Marx_Dump_File_Type *pha_df;
   Marx_Dump_File_Type out_cdf;
   CFits_Response_File_Type *rsp;
   unsigned int num_energies;
   char *pgm;
   int use_stdout = 0;
   FILE *fpout, *en_fp, *det_fp, *chipx_fp, *chipy_fp, *pha_fp;
   unsigned int *sort_index;
   unsigned int num_output;
   char filename[MAX_FILENAME_LEN];

   energy_df = detector_df = chipx_df = chipy_df = pha_df = NULL;
   fpout = en_fp = det_fp = chipx_fp = chipy_fp = pha_fp = NULL;
   sort_index = NULL;

   /* pgm = argv[0]; */
   pgm = "marxrsp";

   argv++; argc--;
   CCD_Id = -1;

   if ((-1 == argcargv (&argc, &argv, ArgcArgv_Table))
       || (Rsp_File == NULL)
       || (Marx_Dir == NULL))
     {
	usage (pgm);
	return 1;
     }

   if (argc)
     {
	fprintf (stderr, "%s: command line option %s not supported.\n",
		 pgm, argv[0]);
	usage (pgm);
	return 1;
     }

   if (Force_Flag == 0)
     {
	if (-1 == check_simulation_parameters ())
	  return 1;
     }

   if (NULL == (energy_df = open_marx_file ("energy.dat", 'E', &en_fp)))
     goto return_error;

   if (CCD_Id != -1)
     {
	if (NULL == (detector_df = open_marx_file ("detector.dat", 'A', &det_fp)))
	  goto return_error;
     }

   if ((Chip_X_Min > 0.0) || (Chip_X_Max < 1024.0))
     {
	if (NULL == (chipx_df = open_marx_file ("xpixel.dat", 'E', &chipx_fp)))
	  goto return_error;
     }

   if ((Chip_Y_Min > 0.0) || (Chip_Y_Max < 1024.0))
     {
	if (NULL == (chipy_df = open_marx_file ("ypixel.dat", 'E', &chipy_fp)))
	  goto return_error;
     }
   
   if ((detector_df != NULL) || (chipy_df != NULL) || (chipx_df != NULL))
     {
	if (NULL == (pha_df = open_marx_file ("pha.dat", 'I', &pha_fp)))
	  goto return_error;
     }

   
   if (NULL == (rsp = cfits_open_rmf_response_file (Rsp_File, 
						    (Force_Flag ? 2 : 1))))
     goto return_error;

   /* Now the output filename */
   sprintf (filename, "%s/%s", Marx_Dir, "rmf_pha.dat");
   memset ((char *) &out_cdf, 0, sizeof (Marx_Dump_File_Type));
   out_cdf.type = 'I';
   out_cdf.num_rows = energy_df->num_rows;
   strcpy (out_cdf.colname, "PHA");

   if (NULL == (fpout = marx_create_write_dump_file (filename, &out_cdf)))
     {
	fprintf (stderr, "*** Unable to create %s\n", filename);
	goto return_error;
     }

   num_output = 0;
   while (1)
     {
	unsigned int i;

	if (-1 == read_energies_etc (&num_energies, en_fp, pha_fp, det_fp, chipx_fp, chipy_fp))
	  goto return_error;

	if (num_energies == 0)
	  break;
	
	/* Sort energies to optimize response lookup. */
	sort_index = JDMsort_floats (Energies, num_energies);
	if (sort_index == NULL)
	  goto return_error;
	 
	if (-1 == compute_response (rsp, Energies, sort_index, Filter, Channels, num_energies))
	  goto return_error;

	for (i = 0; i < num_energies; i++)
	  {
	     if (use_stdout)
	       {
		  unsigned int j = sort_index[i];
		  fprintf (stdout, "%e\t%4d\n", Energies[j], Channels[j]);
	       }
	     else 
	       {
		  int16 i16 = Channels[i];
		  
		  if (1 != JDMwrite_int16 (&i16, 1, fpout))
		    {
		       fprintf (stderr, "write failed.\n");
		       goto return_error;
		    }
	       }
	  }
	num_output += num_energies;
	JDMfree_integer_vector ((int *) sort_index);
	sort_index = NULL;
     }

   if (-1 == marx_close_write_dump_file (fpout, num_output))
     {
	fpout = NULL;
	goto return_error;
     }
   cfits_close_response_file (rsp);
   marx_close_read_dump_file (detector_df);   /* NULL ok */
   marx_close_read_dump_file (energy_df);   /* NULL ok */
   marx_close_read_dump_file (chipx_df);   /* NULL ok */
   marx_close_read_dump_file (chipy_df);   /* NULL ok */
   marx_close_read_dump_file (pha_df); /* NULL ok */
   
   if (-1 == rename_marx_file ("pha.dat", "pha.dat.BAK"))
     return 1;
   if (-1 == rename_marx_file ("rmf_pha.dat", "pha.dat"))
     return 1;
   
   return 0;
   
   
   return_error:
   
   if (sort_index != NULL) JDMfree_integer_vector ((int *)sort_index);
   marx_close_read_dump_file (detector_df);   /* NULL ok */
   marx_close_read_dump_file (energy_df);   /* NULL ok */
   marx_close_read_dump_file (chipx_df);   /* NULL ok */
   marx_close_read_dump_file (chipy_df);   /* NULL ok */
   marx_close_read_dump_file (pha_df); /* NULL ok */
   if (rsp != NULL) cfits_close_response_file (rsp);
   if (fpout != NULL) fclose (fpout);
   return 1;
}

/*}}}*/
