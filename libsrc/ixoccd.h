#ifndef _MARX_IXOCCD_H_INCLUDED
#define _MARX_IXOCCD_H_INCLUDED
/*
    This file is part of MARX

    Copyright (C) 2011-2013 Massachusetts Institute of Technology

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

typedef struct _IXO_CCD_QE_Type IXO_CCD_QE_Type;
static double Rowland_R, Rowland_Theta;
static int CatGS_Init_Called;

static void _marx_catgs_init_variables (void);

#define MARX_DET_FACET_PRIVATE_DATA		\
   IXO_CCD_QE_Type *qeinfo; \
   double read_noise; \
   double energy_gain; \
   double fano_factor;

#endif				       /* _MARX_IXOCCD_H_INCLUDED */
