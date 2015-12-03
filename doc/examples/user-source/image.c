/* This is an example of a user defined source that reads a fits image.
 * Author: John Houck (houck@space.mit.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jdmath.h>
#include <fitsio.h>

#include "user.h"

static float *Image;		/* pointer to Image allocated array */
static unsigned int XSize;	/* dimensions of Image in pixels */
static unsigned int YSize;
static unsigned int Image_Size;	/* Total number of pixels */
static float Rad_Per_Pixel;	/* angular scale of each axis in rad/pixel */

int user_open_source (char **argv, int argc, double area,
		      double cosx, double cosy, double cosz)
{
   double rad_per_deg = PI/180.0;
   unsigned int i;
   int xsize, ysize;
   float sum;
   float cdelt1,cdelt2;           /* angular scale of each axis in deg/pixel */
   float cd1_1,cd2_2,cd1_2,cd2_1;
   fitsfile *fp;                  /* pointer to FITS file */
   int status = 0;                /* FITSIO status flag */
   int iomode = 0;                /* FITSIO I/O flag 0=readonly, 1=readwrite */
   long firstelem = 1;            /* index of first image pixel */
   char comment[FLEN_COMMENT];	  /* comment string from FITS header */
   char err_text[FLEN_STATUS];    /* FITSIO error message string */
   int anynul;                    /* various nulls that FITSIO wants */
   float nullval = 0.0;
   float *pnullval = &nullval;

   (void) area; (void) cosx; (void) cosy; (void) cosz;
   
   if (argc != 1)
     {
	fprintf (stderr, "User Source requires a filename\n");
	return -1;
     }

   fprintf (stdout, "Initializing user source.\n");
   fprintf (stdout, "Reading FITS image %s\n", argv[0]);
   fflush (stdout);
   
   fits_open_file (&fp, argv[0], iomode, &status);
   if (status != 0) 
     {
	fits_get_errstatus (status, err_text);
	fprintf(stderr, "Fitsio Error: %s\n", err_text);
	return -1;
     }

   fits_read_key (fp, TINT, "NAXIS1", &xsize, comment, &status);
   if (status != 0) 
     {
	fits_get_errstatus (status, err_text);
	fprintf(stderr, "Fitsio Error: (NAXIS1) %s\n", err_text);
	return -1;
     }
	    
   fits_read_key (fp, TINT, "NAXIS2", &ysize, comment, &status);
   if (status != 0) 
     {
	fits_get_errstatus (status,err_text);
	fprintf(stderr, "Fitsio Error: (NAXIS2) %s\n", err_text);
	return -1;
     }

   /* Get angular scale of image and convert to radians per pixel */

   fits_read_key (fp, TFLOAT, "CDELT1", &cdelt1, comment, &status);
   cdelt1 = fabs(cdelt1);
   if (status == 0)
     {
       fits_read_key (fp, TFLOAT, "CDELT2", &cdelt2, comment, &status);
       if (status!=0) 
	 {
	   fprintf(stderr, "Error: Could not find scale keywords in header\n");
	   fits_get_errstatus(status,err_text);
	   fprintf(stderr, "Fitsio Error: (CDELT2) %s\n",err_text);
	   return -1;
	 }

	cdelt2 = fabs(cdelt2);

	fprintf(stdout, "Image is %4.2f x %4.2f arcmin\n", 
		cdelt1*xsize*60.0, cdelt2*ysize*60.0);
	
	Rad_Per_Pixel = 0.5*(cdelt1 + cdelt2) * rad_per_deg;

     } 
   else 
     {
	status = 0;
	fits_read_key(fp, TFLOAT, "CD1_1", &cd1_1, comment, &status);       
	cd1_1 = fabs(cd1_1);
	fits_read_key(fp, TFLOAT, "CD2_2", &cd2_2, comment, &status);       
	cd2_2 = fabs(cd2_2);
	fits_read_key(fp, TFLOAT, "CD1_2", &cd1_2, comment, &status);       
	cd1_2 = fabs(cd1_2);
	fits_read_key(fp, TFLOAT, "CD2_1", &cd2_1, comment, &status);       
	cd2_1 = fabs(cd2_1);
	if (status != 0)
	  {
	     fprintf(stderr, "Error: could not find scale keywords in header\n");
	     fits_get_errstatus(status, err_text);
	     fprintf(stderr, "Fitsio Error: %s\n",err_text);
	     return -1;
	  }

	Rad_Per_Pixel = 0.5 * rad_per_deg * (sqrt(cd1_1*cd1_1 + cd1_2*cd1_2) 
					     + sqrt(cd2_1*cd2_1 + cd2_2*cd2_2));

	fprintf(stdout,"Image is %4.2f x %4.2f arcmin\n", 
		xsize * 60.0*Rad_Per_Pixel/rad_per_deg, 
		ysize * 60.0*Rad_Per_Pixel/rad_per_deg);
     }


   fprintf(stdout, "Allocating space for %d x %d image\n", xsize, ysize);
   
   XSize = (unsigned int) xsize;
   YSize = (unsigned int) ysize;

   Image_Size = XSize * YSize;

   Image = (float *) malloc (Image_Size * sizeof(float));
   if (Image==NULL) 
     {
	fprintf(stderr, "Error: malloc failed trying to allocate %u element float array\n",
		Image_Size);
       return -1;
     }

   fprintf (stdout, "Reading image data...\n");
   fflush (stdout);

   fits_read_img(fp, TFLOAT, firstelem, Image_Size, pnullval, Image,
		 &anynul, &status);

   if (status != 0) 
     {
	fits_get_errstatus (status, err_text);
	fprintf(stderr, "Fits Error: %s\n",err_text);
	free ((char *) Image);
	Image = NULL;
	return -1;
     }

   fits_close_file(fp, &status);

   if (status != 0)
     {
	fits_get_errstatus (status,err_text);
	fprintf (stderr, "Fitsio Error: %s\n",err_text);
	free ((char *) Image);
	Image = NULL;
	return -1;
     }

   /* compute a cumulative distribution */
   sum = 0;
   for (i = 0; i < Image_Size; i++)
     {
        sum += Image [i];
        Image [i] = sum;
     }
   
   /* Normalize it. */
   for (i = 0; i < Image_Size; i++)
     {
        Image [i] /= sum;
     }
   
   return 0;
}

void user_close_source (void)
{
   if (Image == NULL) return;
   free (Image);
   Image = NULL;
}

int user_create_ray (double *delta_t, double *energy,
                     double *cosx, double *cosy, double *cosz)
{
   double r;
   unsigned int n;
   double x, y, tmp;

   *delta_t = -1.0;
   
   r = JDMrandom ();
   n = JDMbinary_search_f (r, Image, Image_Size);
   
   y = (double) (n / XSize);
   x = (double) (n % XSize);
   
   y -= 0.5 * YSize;
   x -= 0.5 * XSize;

   y *= Rad_Per_Pixel;
   x *= Rad_Per_Pixel;
   
   tmp = cos (y);

   *cosz = sin (y);
   *cosx = -tmp * cos (x);
   *cosy = tmp * sin (x);

   return 0;
}


#ifdef TESTING

int main (int argc, char **argv)
{
   double delta_t, energy, cosx, cosy, cosz;
   
   if (argc != 2)
     {
	fprintf (stderr, "Usage: %s <FITSIMAGE>\n", argv[0]);
	return 1;
     }
   
   if (-1 == user_open_source (argv[1], 1.e3, 0.0, 0.0, 0.0))
     return 1;

   if (-1 == user_create_ray (&delta_t, &energy, &cosx, &cosy, &cosz))
     return 1;

   user_close_source ();
   return 0;
}
#endif
