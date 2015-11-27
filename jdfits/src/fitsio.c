/*
    Copyright (C) 2002 MIT Center For Space Research

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
/* The routines here are for reading data sections of a fits file */
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "jdfits.h"
#include "_jdfits.h"

int jdfits_read_open_data (JDFits_Type *ft)
{
   unsigned long size, dsize;

   if ((ft->mode & JDFITS_READ_MODE) == 0)
     {
	jdfits_error ("jdfits_read_open_data: Fits file not open for read.");
	return -1;
     }

   if (ft->header == NULL)
     {
	jdfits_error ("jdfits_read_open_data: Fits header is NULL!");
	return -1;
     }

   size = ft->header->size;

   dsize = size % JDFITS_RECORD_SIZE;
   if (dsize) dsize = JDFITS_RECORD_SIZE - dsize;

   ft->bytes_left_to_read = size;
   ft->bytes_padded = dsize;

   return 0;
}

static int read_pad (JDFits_Type *ft)
{
   unsigned char tmpbuf [JDFITS_RECORD_SIZE];
   unsigned int nread;
   unsigned int n = ft->bytes_padded;
#if 0
   unsigned char *p, *pmax;
#endif
   ft->bytes_padded = 0;
   nread = fread (tmpbuf, 1, n, ft->fp);
   if (n != nread)
     {
	jdfits_error ("read_pad: Error reading padded record.  Is this record padded?");
	return -1;
     }
#if 0
   p = tmpbuf;
   pmax = p + n;
   while (p < pmax)
     {
	if (*p)
	  {
	     jdfits_warning ("Data fill area is not padded with null characters.");
	     break;
	  }
	p++;
     }
#endif
   return 0;
}

#ifndef WORDS_BIGENDIAN
static void byte_swap32 (unsigned char *ss, unsigned int n)
{
   unsigned char *p, *pmax, ch;

   if (n <= 0) return;
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

static void byte_swap64 (unsigned char *ss, unsigned int n)
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

static void byte_swap16 (unsigned char *p, unsigned int nread)
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
#endif

unsigned int jdfits_read_bytes (JDFits_Type *ft, unsigned char *b, unsigned int n)
{
   unsigned int nread;

   if (ft->bytes_left_to_read < n)
     {
	n = ft->bytes_left_to_read;
     }

   if (n != (nread = fread (b, 1, n, ft->fp)))
     {
	jdfits_error ("jdfits_read_bytes: read error.");
     }

   ft->bytes_left_to_read -= n;

   if ((ft->bytes_left_to_read == 0) && (ft->bytes_padded))
     {
	(void) read_pad (ft);
     }

   return nread;
}

unsigned int jdfits_read_int16 (JDFits_Type *ft, int16 *ss, unsigned int n)
{
   unsigned int nread;

   if (ft->bytes_left_to_read < 2 * n)
     {
	n = ft->bytes_left_to_read / 2;
     }

   if (n != (nread = fread (ss, 2, n, ft->fp)))
     {
	jdfits_error ("jdfits_read_int16: read error.");
     }

#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) ss, nread);
#endif

   ft->bytes_left_to_read -= n * 2;

   if ((ft->bytes_left_to_read == 0) && (ft->bytes_padded))
     {
	(void) read_pad (ft);
     }

   return nread;
}

unsigned int jdfits_read_int32 (JDFits_Type *ft, int32 *ss, unsigned int n)
{
   unsigned int nread;

   if (ft->bytes_left_to_read < 4 * n)
     {
	n = ft->bytes_left_to_read / 4;
     }

   if (n != (nread = fread (ss, 4, n, ft->fp)))
     {
	jdfits_error ("jdfits_read_int32: read error.");
     }

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, nread);
#endif

   ft->bytes_left_to_read -= n * 4;

   if ((ft->bytes_left_to_read == 0) && (ft->bytes_padded))
     {
	(void) read_pad (ft);
     }

   return nread;
}

unsigned int jdfits_read_float64 (JDFits_Type *ft, float64 *ss, unsigned int n)
{
   unsigned int nread;

   if (ft->bytes_left_to_read < 8 * n)
     {
	n = ft->bytes_left_to_read / 8;
     }

   if (n != (nread = fread (ss, 8, n, ft->fp)))
     {
	jdfits_error ("jdfits_read_float64: read error.");
     }

#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *) ss, nread);
#endif

   ft->bytes_left_to_read -= n * 8;

   if ((ft->bytes_left_to_read == 0) && (ft->bytes_padded))
     {
	(void) read_pad (ft);
     }

   return nread;
}

unsigned int jdfits_read_float32 (JDFits_Type *ft, float32 *ss, unsigned int n)
{
   unsigned int nread;

   if (ft->bytes_left_to_read < 4 * n)
     {
	n = ft->bytes_left_to_read / 4;
     }

   if (n != (nread = fread (ss, 4, n, ft->fp)))
     {
	jdfits_error ("jdfits_read_float32: read error.");
     }

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, nread);
#endif

   ft->bytes_left_to_read -= n * 4;

   if ((ft->bytes_left_to_read == 0) && (ft->bytes_padded))
     {
	(void) read_pad (ft);
     }

   return nread;
}

int jdfits_read_close_data (JDFits_Type *ft)
{
   off_t n = ft->bytes_padded + ft->bytes_left_to_read;
   ft->bytes_left_to_read = 0;
   ft->bytes_padded = 0;
#ifndef SEEK_CUR
# define SEEK_CUR 1
#endif
   if (-1 == FSEEK (ft->fp, n, SEEK_CUR))
     {
	jdfits_error ("jdfits_read_close_data: fseek error.");
	return -1;
     }
   return 0;
}

int jdfits_write_float32 (JDFits_Type *ft, float32 *x, unsigned int n)
{
   int ret;
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) x, n);
#endif
   ret = jdfits_write (ft, (unsigned char *) x, 4 * n);
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) x, n);
#endif
   return ret;
}

int jdfits_write_float64 (JDFits_Type *ft, float64 *x, unsigned int n)
{
   int ret;
#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *) x, n);
#endif
   ret = jdfits_write (ft, (unsigned char *) x, 8 * n);
#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *) x, n);
#endif
   return ret;
}

int jdfits_write_int32 (JDFits_Type *ft, int32 *x, unsigned int n)
{
   int ret;
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) x, n);
#endif
   ret = jdfits_write (ft, (unsigned char *) x, 4 * n);
#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) x, n);
#endif
   return ret;
}

int jdfits_write_int16 (JDFits_Type *ft, int16 *x, unsigned int n)
{
   int ret;
#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) x, n);
#endif
   ret = jdfits_write (ft, (unsigned char *) x, 2 * n);
#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) x, n);
#endif
   return ret;
}

static int jdfits_write_pad (JDFits_Type *ft, unsigned int len, char fill_char)
{
   if (-1 == jdfits_check_mode (ft, JDFITS_WRITE_MODE))
     return -1;

   while (len--)
     {
	if (EOF == putc (fill_char, ft->fp))
	  {
	     jdfits_error ("Error writing to file.");
	     return -1;
	  }
     }
   return 0;
}

int jdfits_flush_output (JDFits_Type *ft, char fill_char)
{
   if (-1 == jdfits_check_mode (ft, JDFITS_WRITE_MODE))
     return -1;

   if (ft->write_buffer_len == 0) return 0;

   if (ft->write_buffer_len != fwrite (ft->write_buffer, 1, ft->write_buffer_len, ft->fp))
     {
	jdfits_error ("Error writing to file.");
	return -1;
     }

   if (-1 == jdfits_write_pad (ft, JDFITS_RECORD_SIZE - ft->write_buffer_len, fill_char))
     {
	return -1;
     }
   ft->write_buffer_len = 0;
   return 0;
}

int jdfits_write (JDFits_Type *ft, unsigned char *buf, unsigned int len)
{
   unsigned int space;
   unsigned int n;
   unsigned char *fbuf;

   if (-1 == jdfits_check_mode (ft, JDFITS_WRITE_MODE)) return -1;
   fbuf = ft->write_buffer;

   while (len)
     {
	n = ft->write_buffer_len;
	space = JDFITS_RECORD_SIZE - n;
	if (space == 0)
	  {
	     if (JDFITS_RECORD_SIZE != fwrite (fbuf, 1, JDFITS_RECORD_SIZE, ft->fp))
	       {
		  jdfits_error ("Error writing to file.");
		  return -1;
	       }
	     n = ft->write_buffer_len = 0;
	     space = JDFITS_RECORD_SIZE;
	  }

	if (len <= space)
	  {
	     space = len;
	  }

	memcpy ((char *) fbuf + n, (char *) buf, space);
	ft->write_buffer_len += space;
	len -= space;
	buf += space;
     }
   return 0;
}

/* String read functions */
unsigned char *jdfits_str_read_int32 (int32 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 4 * n;
   if (s != (unsigned char *) ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *) ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *jdfits_str_read_int16 (int16 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 2 * n;

   if (s != (unsigned char *)ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap16 ((unsigned char *) ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *jdfits_str_read_float64 (float64 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 8 * n;

   if (s != (unsigned char *)ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap64 ((unsigned char *)ss, n);
#endif
   return s + len;
}

/*}}}*/

unsigned char *jdfits_str_read_float32 (float32 *ss, unsigned int n, unsigned char *s) /*{{{*/
{
   unsigned int len = 4 * n;
   if (s != (unsigned char *)ss) memcpy ((char *) ss, (char *) s, len);

#ifndef WORDS_BIGENDIAN
   byte_swap32 ((unsigned char *)ss, n);
#endif
   return s + len;
}

/*}}}*/

#define STR_READ_X_Y_FUN(from,to) \
   unsigned char *jdfits_read_ ## from ## _ ## to (to *y, unsigned int n, unsigned char *s) \
       { unsigned int i; \
	  for (i = 0; i < n; i++) \
	    { from x; s = jdfits_str_read_ ## from(&x,1,s); y[i]=(to)x; } \
	  return s; \
       }

#define STR_READ_TO_Y_FUNS(to) \
   STR_READ_X_Y_FUN(int16,to) STR_READ_X_Y_FUN(int32,to) \
       STR_READ_X_Y_FUN(float32,to) STR_READ_X_Y_FUN(float64,to)

STR_READ_TO_Y_FUNS(short)
STR_READ_TO_Y_FUNS(int)
STR_READ_TO_Y_FUNS(long)
STR_READ_TO_Y_FUNS(float)
STR_READ_TO_Y_FUNS(double)
