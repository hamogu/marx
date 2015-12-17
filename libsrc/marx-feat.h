/*
    This file is part of MARX

    Copyright (C) 2002-2018 Massachusetts Institute of Technology

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
/* Set this to 1 to enable the use of wfold scattering tables */
#define MARX_HAS_WFOLD			1

/* Set this to 1 to allow HRMA pitch and YAW to be specified. */
#define MARX_HAS_HRMA_PITCH_YAW		1

/* Set this to 1 for HRC Drake Flat support */
#define MARX_HAS_DRAKE_FLAT		1

#define MARX_HAS_DITHER			1

#define MARX_HAS_ACIS_STREAK		1

/* One or the other of these may be defined, but not both. */
#define MARX_HAS_ACIS_GAIN_MAP		0
#define MARX_HAS_ACIS_FEF		1
#define MARX_HRMA_HAS_STRUTS		1

#define MARX_HAS_IXO_SUPPORT		1
