#include <fitsio.h>

typedef struct
{
   fitsfile *fio;
   int status;
}
CFits_Type;

/* These functions return TYPE of hdu or -1 upon failure.  Valid types are:
 * HDU_IMAGE, HDU_ATABLE, HDU_BTABLE
 */
extern int cfits_next_hdu (CFits_Type *);
extern int cfits_prev_hdu (CFits_Type *);


extern void cfits_error (char *, ...);
extern void cfits_clear_error (CFits_Type *);
extern char *cfits_malloc (unsigned int);

extern CFits_Type *cfits_open_file (char *, int);
extern CFits_Type *cfits_create_file (char *, int);
extern int cfits_get_header_string (CFits_Type *, char *, char *, char *);
extern int cfits_get_header_integer (CFits_Type *, char *, long *, char *);
extern int cfits_close_file (CFits_Type *);

/* Returns type of extension. */
extern int cfits_locate_extension (CFits_Type *, char *, int (*) (CFits_Type *));
extern int cfits_locate_vextension (CFits_Type *ft, int, char **, int (*)(CFits_Type *));

extern int cfits_get_column_numbers (CFits_Type *, unsigned int, char **, int *, int);
extern long cfits_get_num_rows (CFits_Type *);

extern int cfits_read_column_doubles (CFits_Type *, int, long, long, double *, int);
extern int cfits_read_column_floats (CFits_Type *, int, long, long, float *, int);
extern int cfits_read_column_shorts (CFits_Type *, int, int, long, short *, int);

/* response file stuff */

typedef struct 
{
   CFits_Type *ft;
   long num_rows;
   int energ_lo_col;
   int energ_hi_col;
   int n_grp_col;
   int f_chan_col;
   int n_chan_col;
   int matrix_col;
   
   /* If the column values are constant, they may be moved to keyword values
    * to save space.  In that case, the _col fields of this structure will be 
    * set to -1 and the actual value is specified below.
    */
   long n_grp_val;
   long f_chan_val;
   long n_chan_val;
}
CFits_Response_File_Type;

extern CFits_Response_File_Type *cfits_open_rmf_response_file (char *, int);
extern void cfits_close_response_file (CFits_Response_File_Type *);


typedef struct
{
   unsigned int first_channel;
   unsigned int num_channels;
   double *response;
}
CFits_Response_Element_Type;

typedef struct
{
   unsigned int num_grps;
   CFits_Response_Element_Type *elem;
}
CFits_Response_Vector_Type;

extern void cfits_free_response_vector (CFits_Response_Vector_Type *);
extern int cfits_read_response_vector (CFits_Response_File_Type *, int,
				       CFits_Response_Vector_Type *);

typedef struct 
{
   char *ttype;
   char *tform;
   char *tunit;
}
CFits_Binary_Column_Type;

typedef struct
{
   char *extname;
   long num_rows;
   long num_columns;
   CFits_Binary_Column_Type *column_info;
}
CFits_Binary_Table_Type;

extern int cfits_write_header_string (CFits_Type *, char *, char *, char *);
extern int cfits_write_header_long (CFits_Type *, char *, long, char *);
extern int cfits_init_btable_extension (CFits_Type *, char *);



