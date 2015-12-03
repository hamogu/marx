#include <stdio.h>
#include <math.h>

#include "user.h"

/* This is a math library that will be used to
 * generate the random numbers and binary 
 * search.
 */
#include <jdmath.h>

static float Image[512 * 512];

/* 1/2 angular width of region of image */
static double Field_Width = 5.0/60.0*3.1415926/180.0;

int user_open_source (char *cmd_args, double area,
		      double cosx, double cosy, double cosz)
{
   int ch;
   float sum;
   unsigned int x;
   FILE *fp;
   
   fprintf (stdout, "Reading 512x512 grayscale image %s\n", 
            cmd_args);

   fp = fopen (cmd_args, "rb");
   if (fp == NULL)
     {
        fprintf (stderr, "Unable to open %s.\n", cmd_args);
        return -1;
     }
   
   /* read the array and compute a cumulative 
    * distribution.
    */
   for (x = 0; x < (512 * 512); x++)
     {
        unsigned char uch;

        if (EOF == (ch = getc (fp)))
          {
             fprintf (stderr, "Error reading %s.\n", 
                      cmd_args);
             fclose (fp);
             return -1;
          }
             
        uch = (unsigned char) ch;
        
        sum += (float) uch;
        Image [x] = sum;
     }
   
   fclose (fp);

   /* Normalize it. */
   for (x = 0; x < (512 * 512); x++)
     {
        Image [x] = Image [x] / sum;
     }
   return 0;
}

void user_close_source (void)
{
}

int user_create_ray (double *delta_t, double *energy,
                     double *cosx, double *cosy, double *cosz)
{
   double r;
   unsigned int n;
   double x, y, tmp;

   *delta_t = -1.0;
   
   r = JDMrandom ();
   n = JDMbinary_search_f (r, Image, (512 * 512));
   
   y = (double) (n / 512);
   x = (double) (n % 512);
   
   y -= 256.0;
   x -= 256.0;
   
   x = (x / 256) * Field_Width;
   y = (y / 256) * Field_Width;
   
   tmp = cos(y);

   *cosz = sin (y);
   *cosx = tmp * cos (x);
   *cosy = tmp * sin (x);

   return 0;
}


