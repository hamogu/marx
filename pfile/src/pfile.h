/*
    This file is part of the MIT PFILE Parameter Library

    Copyright (C) 2002 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research
    under contract SV1-61010 from the Smithsonian Institution.

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
#ifndef _PF_H_LOADED_
#define _PF_H_LOADED_

#define PFILE_VERSION "2.4.1"

#define PF_MAX_LINE_LEN 1024

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif
typedef struct _Param_Type Param_Type;
typedef struct _Param_File_Type Param_File_Type;

typedef struct
{
   char *name;
   unsigned int type;
#define PF_INTEGER_TYPE		1
#define PF_BOOLEAN_TYPE		2
#define PF_REAL_TYPE		3      /* also DOUBLE type */
#define PF_DOUBLE_TYPE		4      /* also REAL type */
#define PF_STRING_TYPE		5
#define PF_COMMENT_TYPE		6
#define PF_FILE_TYPE		7
#define PF_INT_TYPE		PF_INTEGER_TYPE
#define PF_UINT_TYPE		8
#define PF_STRING0_TYPE		20     /* non-empty, "" not accepted */
   void *value;
}
Param_Table_Type;

extern int pf_get_parameters (Param_File_Type *, Param_Table_Type *);
extern int pf_set_output_filename (Param_File_Type *, char *);

/* Note: These functions return a malloced pointer. */
extern char *pf_get_output_filename (Param_File_Type *);
extern char *pf_get_input_filename (Param_File_Type *);

extern Param_File_Type *pf_open_parameter_file (char *, char *);
extern int pf_close_parameter_file (Param_File_Type *);
extern void pf_error (char *, ...);
extern int PF_Errno;
#define PF_MALLOC_ERROR			1
#define PF_NUMBER_FORMAT_BAD		2
#define PF_CORRUPT_FIELD		3
#define PF_NOT_IMPLEMENTED		4
#define PF_RANGE_ERROR			5
#define PF_UNKNOWN_PARAMETER		7
#define PF_ACCESS_ERROR			8
#define PF_FILE_NOT_FOUND		9
#define PF_FILE_OPEN_ERROR		10
#define PF_BAD_ARGUMENT			11
#define PF_IO_ERROR			12
#define PF_INDIRECT_ERROR		13
#define PF_UNKNOWN_ERROR		14

extern int pf_get_boolean (Param_File_Type *, char *, int *);
extern int pf_get_integer (Param_File_Type *, char *, int *);
extern int pf_get_double (Param_File_Type *, char *, double *);
extern int pf_get_string (Param_File_Type *, char *, char *, unsigned int);
extern int pf_get_file (Param_File_Type *, char *, char *, unsigned int);
extern int pf_get_uint (Param_File_Type *, char *, unsigned int *);

extern int pf_get_type (Param_File_Type *, char *);

extern void pf_usage (char *, int);
extern int pf_set_value (Param_File_Type *, char *, char *);

extern int pf_learn_value (Param_File_Type *, char *, char *);
extern int pf_learn_string (Param_File_Type *, char *, char *);
extern int pf_learn_double (Param_File_Type *, char *, double);
extern int pf_learn_int (Param_File_Type *, char *, int);

extern int pf_parameter_exists (Param_File_Type *, char *);

extern int pf_add_to_search_path (char *path, int append_flag);

extern Param_File_Type *pf_parse_cmd_line (char *, char *, int, char **);
extern Param_File_Type *pf_parse_cmd_line_no_exit (char *, char *, int, char **);
#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* _PF_H_LOADED_ */

