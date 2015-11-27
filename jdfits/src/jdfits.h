/*
    Copyright (C) 2002-2010 MIT Center For Space Research

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
#ifndef JDFITS_H_INCLUDED
#define JDFITS_H_INCLUDED
#include <stdio.h>
#include <time.h>
#include <sys/types.h>		       /* off_t */

#define JDFITS_VERSION 166	       /* 1.66 */

#ifndef HAS_BASIC_TYPEDEFS_DEFINED
# if defined(__alpha__) || defined(__ALPHA__) || defined(__alpha)
#  define INT16_BASIC_TYPE	short
#  define INT32_BASIC_TYPE	int
#  define FLOAT32_BASIC_TYPE	float
#  define FLOAT64_BASIC_TYPE	double
# else
#  include <limits.h>
/* These provide defaults.  They break on a 64 bit system. */
#  if defined(INT_MAX) && defined(LONG_MAX)
#   if (LONG_MAX != 2147483647) && (INT_MAX == 2147483647)
#    define INT32_BASIC_TYPE	int
#   endif
#  endif
#  ifndef INT32_BASIC_TYPE
#   define INT32_BASIC_TYPE	long
#  endif
#  define INT16_BASIC_TYPE	short
#  define FLOAT32_BASIC_TYPE	float
#  define FLOAT64_BASIC_TYPE	double
# endif
# define HAS_BASIC_TYPEDEFS_DEFINED
#endif

#ifdef SIZEOF_INT
# if SIZEOF_INT == 2
#  undef INT16_BASIC_TYPE
#  define INT16_BASIC_TYPE int
# else
#  if SIZEOF_INT == 4
#   undef INT32_BASIC_TYPE
#   define INT32_BASIC_TYPE int
#  endif
# endif
#endif

#ifdef SIZEOF_SHORT
# if SIZEOF_SHORT == 2
#  undef INT16_BASIC_TYPE
#  define INT16_BASIC_TYPE short
# else
#  if SIZEOF_SHORT == 4
#   undef INT32_BASIC_TYPE
#   define INT32_BASIC_TYPE short
#  endif
# endif
#endif

#ifdef SIZEOF_FLOAT
# if SIZEOF_FLOAT == 4
#  undef FLOAT32_BASIC_TYPE
#  define FLOAT32_BASIC_TYPE float
# else
#  if SIZEOF_FLOAT == 8
#   undef FLOAT64_BASIC_TYPE
#   define FLOAT64_BASIC_TYPE float
#  endif
# endif
#endif

#ifdef SIZEOF_DOUBLE
# if SIZEOF_DOUBLE == 4
#  undef FLOAT32_BASIC_TYPE
#  define FLOAT32_BASIC_TYPE double
# else
#  if SIZEOF_DOUBLE == 8
#   undef FLOAT64_BASIC_TYPE
#   define FLOAT64_BASIC_TYPE double
#  endif
# endif
#endif

#ifndef INT16_TYPEDEFED
# define INT16_TYPEDEFED
  typedef INT16_BASIC_TYPE int16;
  typedef unsigned INT16_BASIC_TYPE uint16;
#endif

#ifndef INT32_TYPEDEFED
# define INT32_TYPEDEFED
  typedef INT32_BASIC_TYPE int32;
  typedef unsigned INT32_BASIC_TYPE uint32;
#endif

#ifndef FLOAT32_TYPEDEFED
# define FLOAT32_TYPEDEFED
  typedef FLOAT32_BASIC_TYPE float32;
#endif

#ifndef FLOAT64_TYPEDEFED
# define FLOAT64_TYPEDEFED
  typedef FLOAT64_BASIC_TYPE float64;
#endif

#define JDFITS_FMT_16 "% 8hd"
#define JDFITS_FMT_32 "% 11ld"

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

extern void jdfits_error (char *fmt, ...);
extern void jdfits_warning (char *fmt, ...);
extern int JDFits_Error_Num;

/* Important note: ALL char * types in this structure point to NON-MALLOCED
 * data in JDFits_Header_Type.
 */
typedef struct
{
   char *name;
   unsigned int type;
#define JDFITS_INT32_TYPE		0x0001U
#define JDFITS_INT16_TYPE		0x0002U
#define JDFITS_INT_TYPE			0x0004U
#define JDFITS_INT_MASK			0x0007U
#define JDFITS_FLOAT32_TYPE		0x0008U
#define JDFITS_FLOAT64_TYPE		0x0010U
#define JDFITS_FLOAT_MASK		0x0018U
#define JDFITS_NUMBER_MASK		0x00FFU

#define JDFITS_COMMENT_TYPE		0x0100U
#define JDFITS_BOOL_TYPE		0x0200U
#define JDFITS_STRING_TYPE		0x0400U
#define JDFITS_BYTE_TYPE		0x0800U
#define JDFITS_BIT_TYPE			0x1000U

   /* If this bit is set, the value is a 4 byte integer that represents
    * a pointer into the heap.
    */
#define JDFITS_POINTER_FLAG		0x4000U

#define JDFITS_ALL_TYPES		0x7FFFU
#define JDFITS_UNPARSED_TYPE		0x8000U
   union
     {
	float64 dval;		       /* 64 bit floats */
	char *sval;		       /* pointer to a string */
	float32 fval;		       /* 32 bit floats */
	int32 lval;		       /* 32 bit */
	int16 hval;		       /* 16 bit */
	int ival;
     }
   v;
   /* This is a pointer to the comment field for the record.
    * For an unparsed card, it is a pointer to the END of the card such that
    * comment - v.sval == length of the unparsed NON-null terminated string.
    * This is NOT null terminated.
    */
   char *comment;
   int comment_len;
}
JDFits_Keyword_Type;

typedef struct
{
   unsigned int type;
   unsigned int nelements;	       /* this is the repeat count * naxis2 */
   union
     {
	float64 *dval;		       /* 64 bit floats */
	char *sval;		       /* pointer to a string */
	float32 *fval;		       /* 32 bit floats */
	int32 *lval;		       /* 32 bit */
	int16 *hval;		       /* 16 bit */
	int *ival;
     }
   v;
}
JDFits_Data_Array_Type;

typedef struct
{
   unsigned int type;		       /* field type from tform */
   unsigned int repeat;		       /* repeat count from tform */
   unsigned int size;		       /* bit size of field (does not include repeat) */
   char *tunit;			       /* field units */
   int has_scaling;		       /* if zero, tscal, tcrvl... not present */
   double tscal;		       /* scale factor */
   double tzero;		       /* offset.  These are applied before display */
   int tnull;
   char tdisp[20];		       /* f90 format strings converted to C */
   char *tdim;			       /*  */
   char *ttype;			       /* label for field */
   int theap;
   double tcrpx;		       /* These are more scaling factors. */
   double tcrvl;		       /* y = TCRVL + TCDLT * (x - TCRVL) */
   double tcdlt;
   JDFits_Data_Array_Type data_array;
}
JDFits_Bintable_Field_Type;

typedef struct
{
   int tfields;			       /* number of fields */
   JDFits_Bintable_Field_Type *finfo;    /*  */
   char *extname;		       /* name of the extension */
   int extver;			       /* extension version number */
   int naxis1, naxis2;
}
JDFits_Bintable_Type;

typedef struct JDFits_Header_Type
{
   unsigned int num_keywords;
   JDFits_Keyword_Type *keys;	       /* these contain char * which
					* point to header_data
					*/
   unsigned char *header_data_buf;	       /* This is malloced. */
   off_t size;		       /* size of the data */

   /* Although this information can be obtained from the keys structure, it
    * is convenient to put it here as well.
    */
   int bitpix, naxis;
   int gcount, pcount;
   JDFits_Keyword_Type *kw_naxis1;       /* pointer to the first axis */

   unsigned int type;		       /* type of header */
#define JDFITS_EXTENSION_HEADER 1	       /* unknown extension type */
#define JDFITS_SIMPLE_HEADER 2	       /* basic fits header */
#define JDFITS_TABLE_HEADER 3	       /* ascii table */
#define JDFITS_BINTABLE_HEADER 4	       /* binary table */
   char *name;			       /* name of header if extension type */
   union
     {
	JDFits_Bintable_Type *bintable;
     }
   ext;
   void (*free_routine)(struct JDFits_Header_Type *);
}
JDFits_Header_Type;

typedef struct
{
   FILE *fp;
   unsigned int mode;
#define JDFITS_READ_MODE  1
#define JDFITS_WRITE_MODE 2
   JDFits_Header_Type *header;

   /* private members. */
   off_t bytes_left_to_read;
   off_t bytes_padded;

   unsigned char *write_buffer;	       /* JDFITS_RECORD_SIZE long */
   unsigned int write_buffer_len;
}
JDFits_Type;

extern unsigned int JDFits_Message_Type;
#define JDFITS_ERRORS   1
#define JDFITS_WARNINGS 2

extern JDFits_Type *jdfits_open_file (char *, int);
extern int jdfits_close_file (JDFits_Type *);
extern int jdfits_skip_to_next_header (JDFits_Type *);

extern JDFits_Keyword_Type *jdfits_parse_keyword (JDFits_Header_Type *, char *, unsigned int);
extern int jdfits_parse_key (JDFits_Keyword_Type *, unsigned int);
int jdfits_keyword_exists (JDFits_Type *ft, char *key);
/* This function returns a malloced copy of the string */
extern int jdfits_read_keyword_string (JDFits_Type *ft, char *key, char **value);
extern int jdfits_read_keyword_dbl (JDFits_Type *ft, char *key, double *);
extern int jdfits_read_keyword_int (JDFits_Type *ft, char *key, int *);

extern int jdfits_bintable_parse_headers (JDFits_Type *);
extern int jdfits_bintable_dump_data (JDFits_Type *, int, FILE *, char *);
extern int jdfits_bintable_read_data (JDFits_Type *);

extern unsigned int jdfits_read_bytes (JDFits_Type *, unsigned char *, unsigned int);
extern unsigned int jdfits_read_int32 (JDFits_Type *, int32 *, unsigned int);
extern unsigned int jdfits_read_int16 (JDFits_Type *, int16 *, unsigned int);
extern unsigned int jdfits_read_float32 (JDFits_Type *, float32 *, unsigned int);
extern unsigned int jdfits_read_float64 (JDFits_Type *, float64 *, unsigned int);

extern unsigned char *jdfits_str_read_int32 (int32 *, unsigned int, unsigned char *);
extern unsigned char *jdfits_str_read_int16 (int16 *, unsigned int, unsigned char *);
extern unsigned char *jdfits_str_read_float64 (float64 *, unsigned int, unsigned char *);
extern unsigned char *jdfits_str_read_float32 (float32 *, unsigned int, unsigned char *);

extern unsigned char *jdfits_read_int16_short (short *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int32_short (short *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float32_short (short *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float64_short (short *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int16_int (int *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int32_int (int *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float32_int (int *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float64_int (int *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int16_long (long *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int32_long (long *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float32_long (long *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float64_long (long *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int16_float (float *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int32_float (float *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float32_float (float *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float64_float (float *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int16_double (double *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_int32_double (double *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float32_double (double *, unsigned int, unsigned char *);
extern unsigned char *jdfits_read_float64_double (double *, unsigned int, unsigned char *);

extern int jdfits_parse_header (JDFits_Type *);
extern int jdfits_read_header (JDFits_Type *);
extern int jdfits_read_open_data (JDFits_Type *);
extern int jdfits_read_close_data (JDFits_Type *);

extern char *jdfits_malloc (unsigned int);
extern void jdfits_free (char *);
extern char *jdfits_make_string (char *);

extern int jdfits_flush_output (JDFits_Type *, char);
extern int jdfits_check_mode (JDFits_Type *, unsigned int);
extern int jdfits_strcasecmp (char *, char *);

extern int jdfits_write_header_logical (JDFits_Type *, char *, int, char *);
extern int jdfits_write_header_integer (JDFits_Type *, char *, int, char *);
extern int jdfits_write_header_string (JDFits_Type *, char *, char *, char *);
extern int jdfits_write_header_double (JDFits_Type *, char *, double, char *);
extern int jdfits_write_header_comment (JDFits_Type *, char *, char *);

extern int jdfits_write_float32 (JDFits_Type *, float32 *, unsigned int);
extern int jdfits_write_float64 (JDFits_Type *, float64 *, unsigned int);
extern int jdfits_write_int32 (JDFits_Type *, int32 *, unsigned int);
extern int jdfits_write_int16 (JDFits_Type *, int16 *, unsigned int);
extern int jdfits_write (JDFits_Type *, unsigned char *, unsigned int);

extern int jdfits_end_header (JDFits_Type *ft);
extern int jdfits_end_data (JDFits_Type *ft);

extern int jdfits_ffmt_to_cfmt (char *, char *);

typedef struct
{
   JDFits_Keyword_Type *kws_start;
   JDFits_Keyword_Type *kws_end;
   JDFits_Keyword_Type *kw_next;
}
JDFits_Read_Keyword_Type;

extern JDFits_Keyword_Type *jdfits_find_keyword (JDFits_Type *, char *);
extern JDFits_Read_Keyword_Type *jdfits_open_keywords (JDFits_Type *);
extern JDFits_Keyword_Type *jdfits_read_keyword (JDFits_Read_Keyword_Type *);
extern void jdfits_rewind_keywords (JDFits_Read_Keyword_Type *);
extern void jdfits_close_keywords (JDFits_Read_Keyword_Type *);
extern int jdfits_write_keyword (JDFits_Type *, JDFits_Keyword_Type *);
extern int jdfits_extract_comment (JDFits_Keyword_Type *, char **);
extern int jdfits_extract_logical (JDFits_Keyword_Type *, int *);
extern int jdfits_extract_string (JDFits_Keyword_Type *, char **);
extern int jdfits_extract_double (JDFits_Keyword_Type *, double *);
extern int jdfits_extract_integer (JDFits_Keyword_Type *, int *);

extern int jdfits_copy_data (JDFits_Type *, JDFits_Type *);
extern int jdfits_copy_header (JDFits_Type *, JDFits_Type *);

   typedef struct
     {
	char *ttype, *ttype_comment;
	char *tform, *tform_comment;
	char *tunit, *tunit_comment;

	/* If the min_max_type is non-zero, the TLMIN/TLMAX fields will be
	 * written out.  Valid non-zero values are: 'A', 'J', 'I', 'E', 'D'.
	 * These should be consistent with tform.
	 */
	char min_max_type;
	union
	  {
	     int32 j_val;
	     float64 d_val;
	  }
	min_value;
	char *min_comment;
	union
	  {
	     int32 j_val;
	     float64 d_val;
	  }
	max_value;
	char *max_comment;

	/* WCS --- if ctype is non-NULL */
	char *ctype;
	double cdelt, crpix, crval;
     }
   JDFits_BTable_Keyword_Type;

   extern int jdfits_create_btable_extension (JDFits_Type *,
					      JDFits_BTable_Keyword_Type *,
					      int,/* naxis2 */
					      int,/* pcount */
					      int,/* gcount */
					      char *);   /* extname */

extern int jdfits_init_null_primary_header (JDFits_Type *);

extern int jdfits_add_headers_from_file (JDFits_Type *, char *);
extern int jdfits_add_comments_from_file (JDFits_Type *, char *, char *,
					  char *, int);

   extern double jdfits_time_t_to_mjd (time_t);

   /* Routines in fitsuser.c */
   typedef struct
     {
	unsigned int type;
#define JDFITS_USER_KW_INT 1
#define JDFITS_USER_KW_STR 2
#define JDFITS_USER_KW_BOOL 3
#define JDFITS_USER_KW_FLOAT 4
#define JDFITS_USER_KW_COMMENT 5
	char *keyword;
	union
	  {
	     double d;
	     char *str;
	     int i;
	  }
	value;
	char *comment;
     }
   JDFits_User_KW_Type;

   typedef struct _JDFits_User_KW_Table_Type JDFits_User_KW_Table_Type;

   /* This function reads a keyword table and returns a pointer to it.  If the
    * second parameter is not NULL, the new keyword table will be appended to
    * the table defined by second parameter and it will be returned.
    */
   extern JDFits_User_KW_Table_Type *jdfits_read_user_kw_table (char *,
								JDFits_User_KW_Table_Type *);
   extern void jdfits_free_user_kw_table (JDFits_User_KW_Table_Type *);
extern JDFits_User_KW_Type *jdfits_find_user_kw (JDFits_User_KW_Table_Type *, char *);
extern int jdfits_write_user_ky (JDFits_Type *, JDFits_User_KW_Type *);

/* These are very crude interface designed to facilitate reading of
 * binary tables.  Yuk!
 */
extern JDFits_Type *jdfits_open_binary_table (char *file, char *extnam);
extern int jdfits_bintable_column_exists (JDFits_Type *f, char *column);

extern JDFits_Type *jdfits_find_binary_table (char *file,
					      int (*fun)(void *, JDFits_Type *),
					      void *clientdata);

typedef struct
{
   unsigned char *row_data;		       /* buffer used to store RAW row bytes */
   unsigned int row_data_len;	       /* same as NAXIS1 parameter */
   unsigned int num_rows;
   unsigned int num_rows_to_read;

   /* The following are constructed from what the user wants to read */
   unsigned int *data_types;
   unsigned int *data_offsets;
   unsigned int num_data_types;

   JDFits_Type *ft;
}
JDFits_BTable_Read_Type;

extern JDFits_BTable_Read_Type *jdfits_simple_aopen_btable (char *,
							    char *,
							    unsigned int,
							    char **);

extern JDFits_BTable_Read_Type *jdfits_simple_open_btable (char *,
							   char *,
							   unsigned int, ...);

extern int jdfits_simple_d_read_btable (JDFits_BTable_Read_Type *, double *);
extern void jdfits_simple_close_btable (JDFits_BTable_Read_Type *);

typedef struct
{
   unsigned int repeat;
   int data_type;		       /* 'a', 's', 'i', 'l', 'f', or 'd' */
   union
     {
	char *a;
	short *s;
	int *i;
	long *l;
	float *f;
	double *d;
     }
   data;
   unsigned int data_offset;
   unsigned char * (*read_fun)(void *, unsigned int, unsigned char *);
}
JDFits_Col_Data_Type;

typedef struct
{
   unsigned int num_columns;
   unsigned int num_rows;
   JDFits_Col_Data_Type *col_data;

   /* private */
   unsigned char *row_bytes;
   off_t num_bytes;
   off_t num_rows_to_read;
}
JDFits_Row_Type;

extern JDFits_Row_Type *jdfits_bintable_aopen_rows (JDFits_Type *f,
						    unsigned int ncols,
						    char **column_names);
extern JDFits_Row_Type *
  jdfits_bintable_open_rows (JDFits_Type *ft, unsigned int ncols, ...);

extern void jdfits_bintable_close_rows (JDFits_Row_Type *rt);
extern int jdfits_read_next_row (JDFits_Type *f, JDFits_Row_Type *r);

#define JDFITS_RECORD_SIZE 2880U
#define JDFITS_CARD_SIZE 80U

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /* JDFITS_H_INCLUDED */
