/*
    This file is part of MARX

    Copyright (C) 1999 Massachusetts Institute of Technology

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
typedef struct
{
   char name[80];
   char *formula;		       /* if molecule */
   float *energy;
   float *f1;
   float *f2;
   unsigned int num_elements;
   float density;
   float atomic_mass;
}
Henke_Type;

extern void henke_free_henke_table (Henke_Type *);

extern int henke_beta_delta (Henke_Type *, float,
			     float **, float **);

extern Henke_Type *henke_read_henke_table (char *);

extern int henke_set_data_dir (char *);

