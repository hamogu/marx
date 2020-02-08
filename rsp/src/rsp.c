/* -*- mode: C; mode: fold; -*- */
#include <stdio.h>

#include <stdlib.h>
#include <string.h>

#include "cfits.h"
#include <jdmath.h>


#define MAX_GRPS 100

   
static int check_rmf_extension (CFits_Type *ft) /*{{{*/
{
   char value [80];

   if ((-1 == cfits_get_header_string (ft, "HDUCLAS2", value, NULL))
       || strcmp (value, "RSP_MATRIX"))
     {
	cfits_clear_error (ft);
	return -1;
     }

   if (-1 == cfits_get_header_string (ft, "HDUCLAS3", value, NULL))
     {
	cfits_clear_error (ft);
	return -1;
     }
   
   if (!strcmp (value, "REDIST"))
     return 1;
   
   if (!strcmp (value, "DETECTOR"))
     return 2;
   
   if (!strcmp (value, "FULL"))
     {
	return 3;
     }
   
   
   return -1;
}

/*}}}*/

/* This routine is complicated by the fact that CAL/GEN/92-002 allows 
 * columns to be specified as a keyword value.
 */
CFits_Response_File_Type *cfits_open_rmf_response_file (char *file, int strict) /*{{{*/
{
   static char *names[6] =
     {
	"ENERG_LO", "ENERG_HI", "N_GRP", "F_CHAN", "N_CHAN", "MATRIX"
     };
   static char *ext_list [2] = 
     {
	"SPECRESP MATRIX", "MATRIX"
     };
   CFits_Type *ft = NULL;
   int columns[6];
   long allow_as_keyword [6];
   CFits_Response_File_Type *rmf = NULL;
   unsigned int i;
   
   if (NULL == (rmf = (CFits_Response_File_Type *) cfits_malloc (sizeof (CFits_Response_File_Type))))
     {
	return NULL;
     }
   memset ((char *) rmf, 0, sizeof (CFits_Response_File_Type));
   
   if (NULL == (ft = cfits_open_file (file, READONLY)))
     goto return_error;
   
   if (-1 == cfits_locate_vextension (ft, 2, ext_list, 
				      (strict ? check_rmf_extension : NULL)))
     {
	if (strict)
	  cfits_error ("This RMF file does not appear to be OGIP approved.");
	else
	  cfits_error ("cfits_open_rmf_response_file: Unable to find RMF extension.");
	goto return_error;
     }
     
   if (strict != 2)
     {
	if (3 == check_rmf_extension (ft))
	  {
	     cfits_error ("*** Warning: RMF file already has effective area folded in.");
	     if (strict)
	       {
		  cfits_error ("Use the --force command line parameter to override.");
		  goto return_error;
	       }
	  }
     }
   

   rmf->num_rows = cfits_get_num_rows (ft);
   
   if (-1 == cfits_get_column_numbers (ft, 6, names, columns, 0))
     goto return_error;
   
   /* This is ugly !!!!!!!! */
   allow_as_keyword [0] = 0;
   allow_as_keyword [1] = 0;
   allow_as_keyword [5] = 0;
   allow_as_keyword [2] = 1;
   allow_as_keyword [3] = 1;
   allow_as_keyword [4] = 1;
   
   for (i = 0; i < 6; i++)
     {
	if (columns[i] != -1) continue;
	if ((allow_as_keyword [i] == 0)
	    || (-1 == cfits_get_header_integer (ft, names[i], &allow_as_keyword[i], NULL)))
	  {
	     cfits_error ("cfits_open_rmf_file: %s column not found.",
			  names[i]);
	     goto return_error;
	  }
     }

   rmf->ft = ft;
   
   rmf->energ_lo_col = columns[0];
   rmf->energ_hi_col = columns[1];
   rmf->matrix_col = columns[5];
   
   if (-1 == (rmf->n_grp_col = columns[2]))
     rmf->n_grp_val = allow_as_keyword[2];
   
   if (-1 == (rmf->f_chan_col = columns[3]))
     rmf->f_chan_val = allow_as_keyword[3];
   
   if (-1 == (rmf->n_chan_col = columns[4]))
     rmf->n_chan_val = allow_as_keyword[4];

   return rmf;
   
   return_error:
   
   if (rmf != NULL) free (rmf);
   if (ft != NULL) cfits_close_file (ft);
   return NULL;
}

/*}}}*/

void cfits_close_response_file (CFits_Response_File_Type *rsp) /*{{{*/
{
   if (rsp == NULL) return;
   if (rsp->ft != NULL)
     (void) cfits_close_file (rsp->ft);
   free (rsp);
}

/*}}}*/


void cfits_free_response_vector (CFits_Response_Vector_Type *v) /*{{{*/
{
   CFits_Response_Element_Type *elem, *elem_max;
   
   if ((v == NULL) || ((elem = v->elem) == NULL))
     return;
   
   elem_max = elem + v->num_grps;
   
   while (elem < elem_max)
     {
	if (elem->response != NULL) JDMfree_double_vector (elem->response);
	elem++;
     }
   
   free (v->elem);
   v->elem = NULL;
   v->num_grps = 0;
}

/*}}}*/
   
int cfits_read_response_vector (CFits_Response_File_Type *rft, int row, /*{{{*/
				CFits_Response_Vector_Type *v)
{
   short ngrps;
   unsigned int i;
   short nchan[MAX_GRPS], fchan[MAX_GRPS];
   CFits_Response_Element_Type *elem, *elem_max;
   CFits_Type *ft;
   int matrix_col;
   long offset;
   
   
   v->num_grps = 0;
   v->elem = NULL;
   
   ft = rft->ft;

   if (rft->n_grp_col == -1)
     ngrps = (short) rft->n_grp_val;
   else if (-1 == cfits_read_column_shorts (ft, rft->n_grp_col, row, 1, &ngrps, 1))
     return -1;
   
   if (ngrps == 0)
     return 0;

   elem = (CFits_Response_Element_Type *) cfits_malloc (ngrps * sizeof (CFits_Response_Element_Type));
   if (elem == NULL)
     return -1;
   
   memset ((char *) elem, ngrps * sizeof(CFits_Response_Element_Type), 0);
   
   v->num_grps = ngrps;
   v->elem = elem;
   
   /* Hopefully the cfitsio library will cache things. In general things will 
    * be in the heap and whether it has a cache or not may be irrelevant.
    */
	
   if (-1 == rft->f_chan_col)
     {
	unsigned int imax = ngrps;
	short f_chan_val = rft->f_chan_val;
	
	for (i = 0; i < imax; i++) fchan [i] = f_chan_val;
     }
   else
     if (-1 == cfits_read_column_shorts (ft, rft->f_chan_col, row, 1, fchan, ngrps))
       goto return_error;

   if (-1 == rft->n_chan_col)
     {
	unsigned int imax = ngrps;
	short n_chan_val = rft->n_chan_val;
	
	for (i = 0; i < imax; i++) nchan [i] = n_chan_val;
     }
   else
     if (-1 == cfits_read_column_shorts (ft, rft->n_chan_col, row, 1, nchan, ngrps))
       goto return_error;
   
   elem = v->elem;
   elem_max = elem + ngrps;
   i = 0;
   
   matrix_col = rft->matrix_col;
   offset = 1;
   
   while (elem < elem_max)
     {
	elem->first_channel = fchan[i];
	elem->num_channels = nchan[i];
	
	elem->response = JDMdouble_vector (elem->num_channels);
	if ((elem->response == NULL)
	    || (-1 == cfits_read_column_doubles (ft, matrix_col, row, offset,
						elem->response, elem->num_channels)))
	  goto return_error;
	
	offset += elem->num_channels;
	i++;
	elem++;
     }
   
   return 0;
   
   return_error:
   cfits_free_response_vector (v);
   return -1;
}

/*}}}*/
