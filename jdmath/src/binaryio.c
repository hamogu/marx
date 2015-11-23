/* -*- mode: C; mode: fold; -*- */
/*
 Copyright (c) 2002,2013 John E. Davis

 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the Free
 Software Foundation; either version 2 of the License, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc., 675
 Mass Ave, Cambridge, MA 02139, USA.
*/
#include "config.h"

#include <stdio.h>
#include <string.h>

#include "jdmath.h"
#include "jdsys.h"

#ifndef WORDS_BIGENDIAN
static void byte_swap64 (unsigned char *ss, unsigned int n) /*{{{*/
{
   unsigned char *p, *pmax, ch;

   if (n <= 0) return;
   p = (unsigned char *) ss;
   pmax = p + 8 * n;
   while (p < pmax)
     {
	ch = *p;
	*p = *(p + 7);
	*(p + 7) = ch;

	ch = *(p + 6);
	*(p + 6) = *(p + 1);
	*(p + 1) = ch;

	ch = *(p + 5);
	*(p + 5) = *(p + 2);
	*(p + 2) = ch;

	ch = *(p + 4);
	*(p + 4) = *(p + 3);
	*(p + 3) = ch;

	p += 8;
     }
}

/*}}}*/
static void byte_swap32 (unsigned char *ss, unsigned int n) /*{{{*/
{
   unsigned char *p, *pmax, ch;

   p = (unsigned char *) ss;
   pmax = p + 4 * n;
   while (p < pmax)
     {
	ch = *p;
	*p = *(p + 3);
	*(p + 3) = ch;

	ch = *(p + 1);
	*(p + 1) = *(p + 2);
	*(p + 2) = ch;
	p += 4;
     }
}

/*}}}*/
static void byte_swap16 (unsigned char *p, unsigned int nread) /*{{{*/
{
   unsigned char *pmax, ch;

   pmax = p + 2 * nread;
   while (p < pmax)
     {
	ch = *p;
	*p = *(p + 1);
	*(p + 1) = ch;
	p += 2;
     }
}

/*}}}*/
#endif

/*{{{ File Read Functions */

unsigned int JDMread_int32 (int32 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nread = fread (ss, 4, n, fp);
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, nread);
#endif
   return nread;
}

/*}}}*/

unsigned int JDMread_int16 (int16 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nread = fread (ss, 2, n, fp);
#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) ss, nread);
#endif
   return nread;
}

/*}}}*/

unsigned int JDMread_float64 (float64 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nread = fread (ss, 8, n, fp);

#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *)ss, nread);
#endif
   return nread;
}

/*}}}*/

unsigned int JDMread_float32 (float32 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nread = fread (ss, 4, n, fp);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *)ss, nread);
#endif
   return nread;
}

/*}}}*/

unsigned int JDMread_f_float32 (float *f, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_FLOAT == 4
   return JDMread_float32 ((float32 *) f, n, fp);
#else
   /* Do it the hard way. */
   float32 f32;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	if (1 != JDMread_float32 (&f32, 1, fp))
	  break;

	f[i] = (float) f32;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMread_d_float64 (double *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_DOUBLE == 8
   return JDMread_float64 ((float64 *) d, n, fp);
#else
   float64 f64;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	if (1 != JDMread_float64 (&f64, 1, fp))
	  break;

	d[i] = (double) f64;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMread_d_float32 (double *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_DOUBLE == 4
   /* Pity the person who has this system */
   return JDMread_float32 ((float32 *) d, n, fp);
#else
   /* Do it the hard way. */
   float32 f32;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	if (1 != JDMread_float32 (&f32, 1, fp))
	  break;

	d[i] = (double) f32;
     }

   return i;
#endif
}

/*}}}*/

unsigned int JDMread_f_float64 (float *f, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_FLOAT == 8
   return JDMread_float64 ((float64 *) f, n, fp);
#else
   float64 f64;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	if (1 != JDMread_float64 (&f64, 1, fp))
	  break;

	f[i] = (float) f64;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMread_l_int32 (long *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_LONG == 4
   return JDMread_int32 ((int32 *) d, n, fp);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 i32;
	if (1 != JDMread_int32 (&i32, 1, fp))
	  break;

	d[i] = (long) i32;
     }

   return i;
#endif
}

/*}}}*/

unsigned int JDMread_l_int16 (long *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_LONG == 2
   return JDMread_int16 ((int16 *) d, n, fp);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 x;
	if (1 != JDMread_int16 (&x, 1, fp))
	  break;

	d[i] = (long) x;
     }

   return i;
#endif
}

/*}}}*/

unsigned int JDMread_i_int16 (int *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_INT == 2
   /* Pity the person who has this system */
   return JDMread_int16 ((int16 *) d, n, fp);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 i16;
	if (1 != JDMread_int16 (&i16, 1, fp))
	  break;

	d[i] = (int) i16;
     }

   return i;
#endif
}

/*}}}*/

unsigned int JDMread_i_int32 (int *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_INT == 4
   /* Pity the person who has this system */
   return JDMread_int32 ((int32 *) d, n, fp);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 x;
	if (1 != JDMread_int32 (&x, 1, fp))
	  break;

	d[i] = (int) x;
     }

   return i;
#endif
}

/*}}}*/

unsigned int JDMread_s_int16 (short *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_SHORT == 2
   return JDMread_int16 ((int16 *) d, n, fp);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 x;
	if (1 != JDMread_int16 (&x, 1, fp))
	  break;

	d[i] = (short) x;
     }

   return i;
#endif
}

/*}}}*/

unsigned int JDMread_s_int32 (short *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_SHORT == 4
   return JDMread_int32 ((int32 *) d, n, fp);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 x;
	if (1 != JDMread_int32 (&x, 1, fp))
	  break;

	d[i] = (short) x;
     }

   return i;
#endif
}

/*}}}*/

/*}}}*/

/*{{{ File Write Functions */

unsigned int JDMwrite_float32 (float32 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nwrote;

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, n);
#endif
   nwrote = fwrite (ss, 4, n, fp);
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, n);
#endif
   return nwrote;
}

/*}}}*/

unsigned int JDMwrite_int32 (int32 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nwrote;

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, n);
#endif
   nwrote = fwrite (ss, 4, n, fp);
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, n);
#endif
   return nwrote;
}

/*}}}*/

unsigned int JDMwrite_int16 (int16 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nwrote;

#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) ss, n);
#endif
   nwrote = fwrite (ss, 2, n, fp);
#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) ss, n);
#endif
   return nwrote;
}

/*}}}*/

unsigned int JDMwrite_float64 (float64 *ss, unsigned int n, FILE *fp) /*{{{*/
{
   unsigned int nwrote;

#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *) ss, n);
#endif
   nwrote = fwrite (ss, 8, n, fp);
#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *) ss, n);
#endif
   return nwrote;
}

/*}}}*/

unsigned int JDMwrite_f_float32 (float *f, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_FLOAT == 4
   return JDMwrite_float32 ((float32 *) f, n, fp);
#else
   float32 f32;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	f32 = (float32) f[i];
	if (1 != JDMwrite_float32 (&f32, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_d_float32 (double *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_DOUBLE == 4
   return JDMwrite_float32 ((float32 *) d, n, fp);
#else
   float32 f32;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	f32 = (float32) d[i];
	if (1 != JDMwrite_float32 (&f32, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_f_float64 (float *f, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_FLOAT == 8
   return JDMwrite_float64 ((float64 *) f, n, fp);
#else
   float64 f64;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	f64 = (float64) f[i];
	if (1 != JDMwrite_float64 (&f64, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_d_float64 (double *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_DOUBLE == 8
   return JDMwrite_float64 ((float64 *) d, n, fp);
#else
   float64 f64;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	f64 = (float64) d[i];
	if (1 != JDMwrite_float64 (&f64, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_s_int32 (short *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_SHORT == 4
   return JDMwrite_int32 ((int32 *) d, n, fp);
#else
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 i32 = (int32) d[i];
	if (1 != JDMwrite_int32 (&i32, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_i_int32 (int *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_INT == 4
   return JDMwrite_int32 ((int32 *) d, n, fp);
#else
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 i32 = (int32) d[i];
	if (1 != JDMwrite_int32 (&i32, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_l_int32 (long *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_LONG == 4
   return JDMwrite_int32 ((int32 *) d, n, fp);
#else
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 i32 = (int32) d[i];
	if (1 != JDMwrite_int32 (&i32, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_i_int16 (int *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_INT == 2
   return JDMwrite_int16 ((int16 *) d, n, fp);
#else
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 i16 = (int16) d[i];
	if (1 != JDMwrite_int16 (&i16, 1, fp))
	  break;
     }
   return i;
#endif
}

/*}}}*/

unsigned int JDMwrite_s_int16 (short *d, unsigned int n, FILE *fp) /*{{{*/
{
#if SIZEOF_SHORT == 2
   return JDMwrite_int16 ((int16 *) d, n, fp);
#else
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 i16 = (int16) d[i];
	if (1 != JDMwrite_int16 (&i16, 1, fp))
	  break;
     }
   return i;
#endif
}
/*}}}*/

/*}}}*/

/*{{{ Reading from a string Functions */

unsigned char *JDMstr_read_int32 (int32 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 4 * n;
   if (s != (unsigned char *) ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_read_int16 (int16 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 2 * n;

   if (s != (unsigned char *)ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_read_float64 (float64 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 8 * n;

   if (s != (unsigned char *)ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *)ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_read_float32 (float32 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 4 * n;
   if (s != (unsigned char *)ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *)ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_read_f_float32 (float *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_FLOAT == 4
   return JDMstr_read_float32 ((float32 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	float32 x;
	s = JDMstr_read_float32 (&x, 1, s);
	f[i] = (float) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_d_float32 (double *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_DOUBLE == 4
   return JDMstr_read_float32 ((float32 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	float32 x;
	s = JDMstr_read_float32 (&x, 1, s);
	f[i] = (double) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_f_float64 (float *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_FLOAT == 8
   return JDMstr_read_float64 ((float64 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	float64 x;
	s = JDMstr_read_float64 (&x, 1, s);
	f[i] = (float) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_d_float64 (double *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_DOUBLE == 8
   return JDMstr_read_float64 ((float64 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	float64 x;
	s = JDMstr_read_float64 (&x, 1, s);
	f[i] = (double) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_s_int16 (short *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_SHORT == 2
   return JDMstr_read_int16 ((int16 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 x;
	s = JDMstr_read_int16 (&x, 1, s);
	f[i] = (short) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_s_int32 (short *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_SHORT == 4
   return JDMstr_read_int16 ((int32 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 x;
	s = JDMstr_read_int32 (&x, 1, s);
	f[i] = (short) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_i_int16 (int *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_INT == 2
   return JDMstr_read_int16 ((int16 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 x;
	s = JDMstr_read_int16 (&x, 1, s);
	f[i] = (int) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_i_int32 (int *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_INT == 4
   return JDMstr_read_int32 ((int32 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 x;
	s = JDMstr_read_int32 (&x, 1, s);
	f[i] = (int) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_l_int16 (long *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_LONG == 2
   return JDMstr_read_int16 ((int16 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int16 x;
	s = JDMstr_read_int16 (&x, 1, s);
	f[i] = (long) x;
     }
   return s;
#endif
}

/*}}}*/

unsigned char *JDMstr_read_l_int32 (long *f, unsigned int n, unsigned char *s) /*{{{*/
{
#if SIZEOF_LONG == 4
   return JDMstr_read_int32 ((int32 *) f, n, s);
#else
   /* Do it the hard way. */
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	int32 x;
	s = JDMstr_read_int32 (&x, 1, s);
	f[i] = (long) x;
     }
   return s;
#endif
}

/*}}}*/

/*}}}*/

/*{{{ Writing to a string Functions */

unsigned char *JDMstr_write_int32 (int32 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 4 * n;

   memcpy ((char *) s, (char *) ss, len);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) s, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_write_int16 (int16 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 2 * n;

   memcpy ((char *) s, (char *) ss, len);

#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) s, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_write_float32 (float32 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 4 * n;

   memcpy ((char *) s, (char *) ss, len);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) s, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *JDMstr_write_float64 (float64 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 8 * n;

   memcpy ((char *) s, (char *) ss, len);

#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *) s, n);
#endif
   return s + len;
}

/*}}}*/

/*}}}*/

