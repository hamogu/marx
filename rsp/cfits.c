/* -*- mode: C; mode: fold; -*- */
#include <stdio.h>
#include <string.h>

#include <stdlib.h>

#include <stdarg.h>

#include "cfits.h"

#define CHECK_FT(ft)   if (((ft) == NULL) || ((ft)->status)) return -1

void cfits_error (char *fmt, ...) /*{{{*/
{
   va_list ap;

   fprintf (stderr, "ERROR: ");
   va_start(ap, fmt);
   (void) vfprintf(stderr, fmt, ap);
   va_end(ap);

   putc('\n', stderr);
}

/*}}}*/

void cfits_clear_error (CFits_Type *ft) /*{{{*/
{
   ft->status = 0;
   ffcmsg ();
}

/*}}}*/

static void error (CFits_Type *ft, char *fun) /*{{{*/
{
   char err [FLEN_ERRMSG];
   char sterr [FLEN_STATUS];

   ffgerr (ft->status, sterr);
   fprintf (stderr, "ERROR in %s: %s\nReasons follow:\n", fun, sterr);
   while (0 != ffgmsg (err))
     {
	fprintf (stderr, "\t%s\n", err);
     }
}

/*}}}*/

char *cfits_malloc (unsigned int len) /*{{{*/
{
   char *p = (char *)malloc (len);

   if (p == NULL)
     {
	cfits_error ("Out of memory.");
     }
   return p;
}

/*}}}*/

CFits_Type *cfits_open_file (char *filename, int mode) /*{{{*/
{
   CFits_Type *ft;

   if (filename == NULL)
     {
	cfits_error ("cfits_open_file: NULL passed.");
	return NULL;
     }

   ft = (CFits_Type *) malloc (sizeof (CFits_Type));
   if (ft == NULL)
     {
	cfits_error ("cfits_open_file: malloc error");
	return NULL;
     }

   ft->status = 0;
   ft->fio = NULL;

   if (ffopen (&ft->fio, filename, mode, &ft->status))
     {
	error (ft, "cfits_open_file");
	free (ft);
	return NULL;
     }

   return ft;
}

/*}}}*/

CFits_Type *cfits_create_file (char *filename, int null_table)
{
   CFits_Type *ft;

   if (filename == NULL)
     {
	cfits_error ("cfits_create_file: NULL passed.");
	return NULL;
     }

   ft = (CFits_Type *) malloc (sizeof (CFits_Type));

   if (ft == NULL)
     {
	cfits_error ("cfits_create_file: malloc error");
	return NULL;
     }
   memset ((char *) ft, 0, sizeof (CFits_Type));

   if (ffinit (&ft->fio, filename, &ft->status))
     goto return_error;

   if (null_table)
     {
	if (ffcrim (ft->fio, -32, 0, NULL, &ft->status))
	  {
	     cfits_error ("cfits_create_file: unable to create NULL extension.");
	     goto return_error;
	  }
     }

   return ft;

   return_error:

   if (ft != NULL)
     {
	error (ft, "cfits_create_file");
	free (ft);
     }
   return NULL;
}

int cfits_get_header_string (CFits_Type *ft, char *name, char *value, char *comment) /*{{{*/
{
   char cbuf [FLEN_COMMENT];

   CHECK_FT(ft);

   if (comment == NULL) comment = cbuf;

   if (0 != ffgkys (ft->fio, name, value, comment, &ft->status))
     {
	if (ft->status != KEY_NO_EXIST)
	  error (ft, "cfits_get_header_string");
	return -1;
     }
   return 0;
}

/*}}}*/

int cfits_get_header_integer (CFits_Type *ft, char *name, long *value, char *comment) /*{{{*/
{
   char cbuf [FLEN_COMMENT];

   CHECK_FT(ft);

   if (comment == NULL) comment = cbuf;

   if (0 != ffgkyj (ft->fio, name, value, comment, &ft->status))
     {
	if (ft->status != KEY_NO_EXIST)
	  error (ft, "cfits_get_header_integer");
	return -1;
     }
   return 0;
}

/*}}}*/

int cfits_close_file (CFits_Type *ft) /*{{{*/
{
   if (ft == NULL)
     return -1;

   if (ffclos (ft->fio, &ft->status))
     {
	error (ft, "cfits_close_file");
	return -1;
     }
   free (ft);
   return 0;
}

/*}}}*/

int cfits_next_hdu (CFits_Type *ft) /*{{{*/
{
   int type;

   CHECK_FT(ft);

   ffmrhd (ft->fio, 1, &type, &ft->status);
   if (ft->status)
     {
	if (ft->status != END_OF_FILE)
	  error (ft, "cfits_next_hdu");
	else cfits_clear_error (ft);

	return -1;
     }
   return type;
}

/*}}}*/

int cfits_prev_hdu (CFits_Type *ft) /*{{{*/
{
   int type;

   CHECK_FT(ft);

   ffmrhd (ft->fio, -1, &type, &ft->status);
   if (ft->status)
     {
	error (ft, "cfits_prev_hdu");
	return -1;
     }
   return type;
}

/*}}}*/

int cfits_locate_vextension (CFits_Type *ft, int argc, char **ext_names, /*{{{*/
			     int (*fun) (CFits_Type *))
{
   int type;
   char buf[FLEN_VALUE];

   CHECK_FT(ft);

   /* Move to fits extension */
   ffmahd (ft->fio, 2, &type, &ft->status);
   if (ft->status)
     {
	error (ft, "cfits_locate_extension");
	return -1;
     }

   do
     {
	if (0 == cfits_get_header_string (ft, "EXTNAME", buf, NULL))
	  {
	     int i;

	     for (i = 0; i < argc; i++)
	       {
		  if (!strcmp (buf, ext_names[i]))
		    {
		       if ((fun == NULL) || (-1 != (*fun) (ft)))
			 return type;
		    }
	       }
	  }

	if (ft->status == KEY_NO_EXIST)
	  cfits_clear_error (ft);
     }
   while (-1 != (type = cfits_next_hdu (ft)));

   return -1;
}

/*}}}*/

int cfits_locate_extension (CFits_Type *ft, char *extname, int (*fun) (CFits_Type *)) /*{{{*/
{
   int type = cfits_locate_vextension (ft, 1, &extname, fun);

   if (type == -1)
     cfits_error ("cfits_locate_extension: Unable to find extension: %s", extname);

   return type;
}

/*}}}*/

int cfits_get_column_numbers (CFits_Type *ft, unsigned int num, char **names, int *cols, /*{{{*/
			      int error_flag)
{
   unsigned int i;
   CHECK_FT(ft);

   for (i = 0; i < num; i++)
     {
	if (ffgcno (ft->fio, 1, names[i], cols + i, &ft->status))
	  {
	     if (error_flag)
	       {
		  error (ft, "cfits_get_column_numbers");
		  return -1;
	       }
	     else cfits_clear_error (ft);
	     cols[i] = -1;
	  }
     }
   return 0;
}

/*}}}*/

long cfits_get_num_rows (CFits_Type *ft) /*{{{*/
{
   long nrows;
   CHECK_FT(ft);

   if (-1 == cfits_get_header_integer (ft, "NAXIS2", &nrows, NULL))
     return -1;

   return nrows;
}

/*}}}*/

int cfits_read_column_floats (CFits_Type *ft, int col, long row, long ofs, /*{{{*/
			      float *data, int nrows)
{
   int anynul;

   CHECK_FT(ft);

   if (nrows == 0) return 0;

   if (ffgcve (ft->fio, col, row, ofs, nrows, 0, data, &anynul, &ft->status))
     {
	error (ft, "cfits_read_column_floats");
	return -1;
     }
   return 0;
}

/*}}}*/

int cfits_read_column_doubles (CFits_Type *ft, int col, long row, long ofs, /*{{{*/
			       double *data, int nrows)
{
   int anynul;

   CHECK_FT(ft);

   if (nrows == 0) return 0;

   if (ffgcvd (ft->fio, col, row, ofs, nrows, 0, data, &anynul, &ft->status))
     {
	error (ft, "cfits_read_column_floats");
	return -1;
     }
   return 0;
}

/*}}}*/

int cfits_read_column_shorts (CFits_Type *ft, int col, int row, long ofs, /*{{{*/
			      short *data, int nrows)
{
   int anynul;

   CHECK_FT(ft);

   (void) ofs;

   if (nrows == 0) return 0;

   if (ffgcvi (ft->fio, col, row, 1, nrows, 0, data, &anynul, &ft->status))
     {
	error (ft, "cfits_read_column_shorts");
	return -1;
     }
   return 0;
}

/*}}}*/

int cfits_write_header_string (CFits_Type *ft, char *kw, char *val, char *com)
{
   CHECK_FT(ft);

   if (com == NULL) com = "";

   if (ffpkys (ft->fio, kw, val, com, &ft->status))
     {
	cfits_error ("cfits_write_header_string");
	return -1;
     }
   return 0;
}

int cfits_write_header_long (CFits_Type *ft, char *kw, long val, char *com)
{
   CHECK_FT(ft);

   if (com == NULL) com = "";

   if (ffpkyj (ft->fio, kw, val, com, &ft->status))
     {
	cfits_error ("cfits_write_header_string");
	return -1;
     }
   return 0;
}

int cfits_init_btable_extension (CFits_Type *ft, char *ext_name)
{
   CHECK_FT(ft);

   if (ffcrhd (ft->fio, &ft->status))
     {
	cfits_error ("cfits_init_btable_extension");
	return -1;
     }
#if 1
   if (-1 == cfits_write_header_string (ft, "XTENSION", "BINTABLE", NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "BITPIX", 8, NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "NAXIS", 2, NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "NAXIS1", 0, NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "NAXIS2", 0, NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "PCOUNT", 0, NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "GCOUNT", 1, NULL))
     return -1;

   if (-1 == cfits_write_header_long (ft, "TFIELDS", 0, NULL))
     return -1;

#endif

   if (-1 == cfits_write_header_string (ft, "EXTNAME", ext_name, NULL))
     return -1;

   return 0;
}
