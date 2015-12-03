#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <jdmath.h>
#include <marx.h>
#include <fitsio.h>

#include "user.h"

static fitsfile *Events_Ptr;
static int X_Col_Num, Y_Col_Num, Pi_Col_Num, Expno_Col_Num;
static double X_CrVal, X_CrPix, X_CDelt;
static double Y_CrVal, Y_CrPix, Y_CDelt;
static double TStart;
static double Frame_Time = 3.240;
static double Nominal_Ra, Nominal_Dec, Nominal_Roll;
static JDMVector_Type Nominal_Pointing;
static JDMVector_Type Rotation_Vector;
static double Rotation_Sin_Theta, Rotation_Cos_Theta;
static JDMVector_Type Roll_Vector;
static double Roll_Sin_Theta, Roll_Cos_Theta;
static JDMVector_Type Marx_Pointing;

/* Use a simple CCD model based upon Fano statistics */
static double Read_Noise = 10.0;
static double Energy_Gain = 0.00368;
static double Fano_Factor = 0.115;

typedef struct
{
   int pi_chan;		       /* Used if no line list given */

   float *cum_strengths;	       /* malloced */
   float *energies;		       /* pointers to line list-- not malloced */
   unsigned int num_energies;
}
Spectral_Bin_Type;

static float *Energy_Line_List;
static float *Line_Strengths;
static unsigned int Num_Lines;
static int Have_Line_List = 0;

#define MAX_NUM_CHANNELS	2048
static double KeV_Per_Bin = 14.6e-3;

static Spectral_Bin_Type Spectral_Bin_Table[MAX_NUM_CHANNELS];

static void do_free (char *s)
{
   if (s != NULL) free (s);
}

static char *do_malloc (unsigned int n)
{
   char *s = malloc (n);
   if (s == NULL)
     fprintf (stderr, "Out of memory\n");
   return s;
}

static double channel_to_kev (int chan)
{
   return KeV_Per_Bin * (chan - 1 + JDMrandom ());
}

static double compute_ccd_sigma (double en)
{
   return Energy_Gain * sqrt (Read_Noise * Read_Noise
			      + Fano_Factor * en / Energy_Gain);
}

static int init_spectral_bin (int chan, Spectral_Bin_Type *sb)
{
   double en, emin, emax, sigma;
   unsigned int min_line, max_line;
   unsigned int num, i;
   double sum, last;
   float *energies, *cum_strengths, *strengths;

   if (Have_Line_List == 0)
     {
	sb->cum_strengths = NULL;
	sb->energies = NULL;
	sb->num_energies = 0;
	sb->pi_chan = chan;
	return 0;
     }

   en = channel_to_kev (chan);
   sigma = compute_ccd_sigma (en);
   emin = en - 3 * sigma;
   emax = en + 3 * sigma;

   if (emin < 0.0) emin = 0.0;
   min_line = JDMbinary_search_f (emin, Energy_Line_List, Num_Lines);
   max_line = JDMbinary_search_f (emax, Energy_Line_List, Num_Lines);

   if (min_line == Num_Lines-1)
     {
	sb->cum_strengths = NULL;
	sb->energies = NULL;
	sb->num_energies = 0;
	return 0;
     }

   num = max_line - min_line + 1;
   if (NULL == (cum_strengths = (float *) do_malloc (num * sizeof (float))))
     return -1;

   sb->cum_strengths = cum_strengths;
   sb->energies = energies = Energy_Line_List + min_line;
   sb->num_energies = num;
   strengths = Line_Strengths + min_line;

   if (num == 1)
     {
	cum_strengths[0] = 1;
	return 0;
     }

   sigma *= sqrt(2);
   sum = 0.0;
   last = JDMerf ((energies[0] - en)/sigma);
   for (i = 1; i < num; i++)
     {
	double next = JDMerf ((energies[i] - en)/sigma);
	sum += (next - last) * strengths[i-1];
	cum_strengths[i-1] = sum;
	last = next;
     }

   if (sum == 0.0)
     {
	for (i = 1; i < num; i++)
	  {
	     sum += 1.0;
	     cum_strengths[i-1] = sum;
	  }
     }

   num--;
   for (i = 0; i < num; i++)
     cum_strengths[i] /= sum;

   return 0;
}

static void free_line_list (void)
{
   do_free ((char *) Energy_Line_List);
   do_free ((char *) Line_Strengths);
   Energy_Line_List = NULL;
   Line_Strengths = NULL;
   Num_Lines = 0;
}

static void free_spectral_bins (void)
{
   Spectral_Bin_Type *b, *bmax;

   b = Spectral_Bin_Table;
   bmax = b + MAX_NUM_CHANNELS;

   while (b < bmax)
     {
	do_free ((char *) b->cum_strengths);
	b->cum_strengths = NULL;
	b++;
     }
}

static int read_line_list (char *file)
{
   float *energies, *strengths;
   unsigned int num;

   if (0 == strcmp (file, "NONE"))
     {
	Have_Line_List = 0;
	Num_Lines = 0;
	return 0;
     }
   Have_Line_List = 1;

   if (-1 == JDMread_float_xy (file, &energies, &strengths, 1, 2, &num))
     {
	fprintf (stderr, "Unable to read spectral file %s\n", file);
	return -1;
     }

   Num_Lines = num;
   Line_Strengths = strengths;
   Energy_Line_List = energies;

   if (Num_Lines == 0)
     {
	fprintf (stderr, "Spectral file %s is invalid", file);
	free_line_list ();
	return -1;
     }

   return 0;
}

static int init_spectral_bins (char *file)
{
   int i;

   if (-1 == read_line_list (file))
     return -1;

   for (i = 1; i < MAX_NUM_CHANNELS; i++)
     {
	if (-1 == init_spectral_bin (i, Spectral_Bin_Table + i))
	  {
	     free_line_list ();
	     free_spectral_bins ();
	  }
     }

   return 0;

}

static void dump_fits_error (int status)
{
   if (status)
     fits_report_error (stderr, status);
}

static int read_key_double (fitsfile *f, char *key, double *v)
{
   int status = 0;

   if (0 == fits_read_key (f, TDOUBLE, key, v, NULL, &status))
     return 0;

   fprintf (stderr, "Unable to read %s keyword value:\n", key);
   dump_fits_error (status);
   return -1;
}

static int read_wcs_info (fitsfile *f, int col,
			  double *crval, double *crpix, double *cdelt)
{
   char key[32];

   sprintf (key, "TCRVL%d", col);
   if (-1 == read_key_double (f, key, crval))
     return -1;

   sprintf (key, "TCRPX%d", col);
   if (-1 == read_key_double (f, key, crpix))
     return -1;

   sprintf (key, "TCDLT%d", col);
   if (-1 == read_key_double (f, key, cdelt))
     return -1;

   return 0;
}

static int get_col_num (fitsfile *f, char *name, int *col)
{
   int status = 0;

   if (0 == fits_get_colnum (f, CASEINSEN, name, col, &status))
     return 0;

   dump_fits_error (status);
   return -1;
}

static void setup_rotations (double new_ra, double new_dec)
{
   JDMVector_Type p1, p2, n;
   double cos_theta, sin_theta;

   p1 = JDMv_spherical_to_vector (1.0, PI/2 - new_dec, new_ra);
   cos_theta = JDMv_dot_prod (p1, Marx_Pointing);
   if (cos_theta > 1.0) cos_theta = 1.0;
   else if (cos_theta < -1.0) cos_theta = -1.0;

   n = JDMv_cross_prod (p1, Marx_Pointing);
   sin_theta = JDMv_length (n);
   if (sin_theta > 1.0)
     sin_theta = 1.0;

   Rotation_Sin_Theta = sin_theta;
   Rotation_Cos_Theta = cos_theta;

   Rotation_Vector = JDMv_unit_vector (n);

   /* See what this vector does to the tangent plane vectors.  In particular,
    * look what it does to the vector in the x tangent plane direction here.
    */
   p1.x = -sin(new_ra);
   p1.y = cos(new_ra);
   p1.z = 0;

   p1 = JDMv_rotate_unit_vector1 (p1, Rotation_Vector, Rotation_Cos_Theta, Rotation_Sin_Theta);
   /* This should be perpindicular to Marx_Pointing.  So look at the y and z
    * components.  Back in the marx system, this should be rolled to appear
    * to lie along the +y axis
    */
   p2.x = 0; p2.y = 1; p2.z = 0;

   cos_theta = JDMv_dot_prod (p1, p2);
   if (cos_theta > 1.0) cos_theta = 1.0;
   else if (cos_theta < -1.0) cos_theta = -1.0;

   n = JDMv_cross_prod (p1, p2);
   sin_theta = JDMv_length (n);
   if (sin_theta > 1.0)
     sin_theta = 1.0;

   Roll_Sin_Theta = sin_theta;
   Roll_Cos_Theta = cos_theta;
   Roll_Vector = JDMv_unit_vector (n);
}

static void close_file (fitsfile *f)
{
   int status = 0;

   if (f != NULL)
     (void) fits_close_file (f, &status);
}

static int usage (void)
{
   fprintf (stderr, "EventList User Source Model Usage:\n");
   fprintf (stderr, " FitsEventsFile SpectralFile [RA_NOM DEC_NOM]\n");
   return -1;
}

/* This function is called with (cosx,cosy,cosz) denoting the unit vector
 * in marx spacecraft coordinates that corresponds to the desired
 * location of the optical axis.  This will place a source at this
 * position on axis (Duhh!).  To simulate an event file as-is, the user
 * should use the nominal ra and dec as given in the event file for the
 * source position.
 */
int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   char *file, *extname;
   fitsfile *f;
   int status = 0;
   double want_ra, want_dec;
   char *spectral_file, *ra_str, *dec_str;

   file = spectral_file = ra_str = dec_str = NULL;

   switch (argc)
     {
      default:
	return usage ();

      case 4:
	ra_str = argv[2];
	dec_str = argv[3];
	if (1 != sscanf (ra_str, "%lf", &want_ra))
	  return usage ();

	if (1 != sscanf (dec_str, "%lf", &want_dec))
	  return usage ();

	want_ra *= PI/180;
	want_dec *= PI/180;

	/* drop */

      case 2:
	file = argv[0];
	spectral_file = argv[1];
	break;
     }

   extname = "EVENTS";

   (void) fits_open_file (&f, file, READONLY, &status);
   if (status)
     {
	fprintf (stderr, "Unable to open fits file %s\n", file);
	dump_fits_error (status);
	return -1;
     }

   if (0 != fits_movnam_hdu (f, BINARY_TBL, extname, 1, &status))
     {
	dump_fits_error (status);
	close_file (f);
	return -1;
     }

   if ((-1 == get_col_num (f, "X", &X_Col_Num))
       || (-1 == get_col_num (f, "Y", &Y_Col_Num))
       || (-1 == get_col_num (f, "EXPNO", &Expno_Col_Num))
       || (-1 == get_col_num (f, "PI", &Pi_Col_Num)))
     {
	close_file (f);
	return -1;
     }

   if ((-1 == read_wcs_info (f, X_Col_Num, &X_CrVal, &X_CrPix, &X_CDelt))
      || (-1 == read_wcs_info (f, Y_Col_Num, &Y_CrVal, &Y_CrPix, &Y_CDelt)))
     {
	close_file (f);
	return -1;
     }

   if ((-1 == read_key_double (f, "RA_NOM", &Nominal_Ra))
       || (-1 == read_key_double (f, "DEC_NOM", &Nominal_Dec))
       || (-1 == read_key_double (f, "ROLL_NOM", &Nominal_Roll))
       || (-1 == read_key_double (f, "TSTART", &TStart)))
     {
	close_file (f);
	return -1;
     }
   Nominal_Ra *= PI/180.0;
   Nominal_Dec *= PI/180.0;
   Nominal_Roll *= PI/180.0;

   Nominal_Pointing = JDMv_spherical_to_vector (1.0, PI/2 - Nominal_Dec, Nominal_Ra);

   Marx_Pointing.x = -cosx;
   Marx_Pointing.y = -cosy;
   Marx_Pointing.z = -cosz;

   if (ra_str == NULL)
     want_ra = Nominal_Ra;

   if (dec_str == NULL)
     want_dec = Nominal_Dec;

   setup_rotations (want_ra, want_dec);

   if (-1 == init_spectral_bins (spectral_file))
     {
	close_file (f);
	return -1;
     }

   Events_Ptr = f;
   return 0;
}

void user_close_source (void)
{
   close_file (Events_Ptr);
   free_line_list ();
   free_spectral_bins ();
   Events_Ptr = NULL;
}

static int read_row_value (fitsfile *f, int row, int col, double *val)
{
   int status = 0;
   int anynul = 0;
   if (0 == fits_read_col (f, TDOUBLE, col, row, 1, 1, NULL, val, &anynul, &status))
     return 0;
   /* dump_fits_error (status); */
   return -1;
}
static int read_row_value_short (fitsfile *f, int row, int col, short *val)
{
   int status = 0;
   int anynul = 0;
   if (0 == fits_read_col (f, TSHORT, col, row, 1, 1, NULL, val, &anynul, &status))
     return 0;
   /* dump_fits_error (status); */
   return -1;
}
static int read_row_value_long (fitsfile *f, int row, int col, long *val)
{
   int status = 0;
   int anynul = 0;
   if (0 == fits_read_col (f, TLONG, col, row, 1, 1, NULL, val, &anynul, &status))
     return 0;
   /* dump_fits_error (status); */
   return -1;
}

static void map_xy_to_unit_vector (double x, double y,
				   double *cosx, double *cosy, double *cosz)
{
   JDMVector_Type p;

   /* Upon input, (x,y) are in the tangent plane system.  Map them to offsets
    * upon the unit circle.
    *
    * Do this in 2 stages.
    */

   /* First, apply the WCS to get RA/Dec tangent plane coordinates. */
   x = X_CrVal + (x - X_CrPix) * X_CDelt;
   y = Y_CrVal + (y - Y_CrPix) * Y_CDelt;

   /* Sadly, these are in degrees so convert them to Nature's preferred unit */
   x *= PI/180.0;
   y *= PI/180.0;

   /* Now, subtract off the nominal pointing so that we are left with a an
    * offset from that pointing in the tangent plane.
    */
   x = (x-Nominal_Ra);/* *cos(Nominal_Dec); */
   y -= Nominal_Dec;

   /* In the tangent plane, positive RA runs to the left.  This corresponds
    * to the MARX +Y axis.
    */

   (void) marx_tan_plane_to_vector (x, y, &Nominal_Pointing, &p);

   /* Now rotate this vector to align it with the desired pointing */

   p = JDMv_rotate_unit_vector1 (p, Rotation_Vector,
				 Rotation_Cos_Theta, Rotation_Sin_Theta);

   p = JDMv_rotate_unit_vector1 (p, Roll_Vector,
				 Roll_Cos_Theta, Roll_Sin_Theta);

   /* Make ray point from sky to origin */
   *cosx = -p.x;
   *cosy = -p.y;
   *cosz = -p.z;
}

static int get_time (double *timeptr, int perform_reset)
{
   static int row = 1;
   static int num_events = 0;
   static long next_expno = 0;
   static double dt, t = 0;
   static long expno_offset = 0;

   if (num_events == 0)
     {
	long expno;

	if (perform_reset)
	  {
	     row = 1;
	     expno_offset = next_expno;
	     next_expno = 0;
	  }

	if (-1 == read_row_value_long (Events_Ptr, row, Expno_Col_Num, &expno))
	  return -1;

	num_events++;
	row++;

	while (0 == read_row_value_long (Events_Ptr, row, Expno_Col_Num, &next_expno))
	  {
	     row++;
	     if (next_expno != expno)
	       break;

	     num_events++;
	  }

	expno += expno_offset;
	t = (expno - 1) * Frame_Time;
	dt = Frame_Time / num_events;
     }

   num_events--;
   *timeptr = t + JDMrandom () * dt;
   t += dt;

   return 0;
}

static double map_pi_to_line (short pi)
{
   Spectral_Bin_Type *b;
   unsigned int i;

   if ((pi < 0) || (pi >= MAX_NUM_CHANNELS))
     return -1.0;

   b = Spectral_Bin_Table + pi;
   if (Have_Line_List == 0)
     return channel_to_kev (b->pi_chan);

   if (NULL == b->cum_strengths)
     {
	if (NULL == b->energies)
	  return -1.0;
	return b->energies[0];
     }
   i = JDMbinary_search_f (JDMrandom (), b->cum_strengths, b->num_energies);
   return b->energies[i];
}

static int read_data (double *tp, double *xp, double *yp, short *chanp)
{
   static int row = 0;
   static double last_t = -1;
   double t;

   if (-1 == get_time (&t, 0))
     {
	row = 0;
	if (-1 == get_time (&t, 1))
	  return -1;
     }

   row++;

   if (last_t < 0.0)
     last_t = t;
   *tp = (t - last_t);
   last_t = t;

   if ((-1 == read_row_value_short (Events_Ptr, row, Pi_Col_Num, chanp))
       || (-1 == read_row_value (Events_Ptr, row, X_Col_Num, xp))
       || (-1 == read_row_value (Events_Ptr, row, Y_Col_Num, yp)))
     return -1;

   return 0;
}

int user_create_ray (double *delta_t, double *energy,
		     double *cosx, double *cosy, double *cosz)
{
   double x, y;
   short chan;

   while (1)
     {
	if (-1 == read_data (delta_t, &x, &y, &chan))
	  return -1;
	if (chan >= 1024)
	  continue;

	if ((*energy = map_pi_to_line (chan)) < 0.0)
	  continue;

	map_xy_to_unit_vector (x, y, cosx, cosy, cosz);
	return 0;
     }
}

#ifdef TESTING
int main (int argc, char **argv)
{
   double dt, en, cx, cy, cz;
   double t = 0.0;
   static char *args[] =
     {
	"test.fits", "absorb.out", "350.86418", "58.81611", NULL
     };

   if (-1 == user_open_source (args, 4, 0, 1, 0, 0))
     return 1;

   while (0 == user_create_ray (&dt, &en, &cx, &cy, &cz))
     {
	t += dt;
	if (t < 7000)
	  continue;
	fprintf (stdout, "%e\t%e\t%e\t%e\t%e\n", t, en, cx, cy, cz);
	break;
     }

   user_close_source ();
   return 0;
}
#endif
