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
extern void *_JDMmalloc (unsigned long, char *);
extern void _JDMfree (void *);
extern void *_JDMrealloc (void *, unsigned int);

extern void _JDMswap_dvector (double *, double *, unsigned int);
extern double _JDM_innerprod (double *, double *, unsigned int);
extern void _JDM_add_to_vector (double *, double *, double, unsigned int);
extern double _JDM_innerprod_col (double *, double **, unsigned int, unsigned int);
extern double *_JDM_equilibrate (double **, unsigned int, double *);
extern double _JDM_nvector_max_abs (double *, unsigned int);
extern double _JDM_nvector_max (double *, unsigned int);
extern double _JDM_nvector_sum (double *, unsigned int);

extern int _JDM_init_machine_constants (void (*)(double *, double *));

#define JDMATH_SING_MATRIX_ERROR	9

