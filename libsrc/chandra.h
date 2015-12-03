#ifndef _MARX_CHANDRA_H_INCLUDED
#define _MARX_CHANDRA_H_INCLUDED
/* Chandra-specific structures */
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


/* stt_lsi_offset: origin of LSI in STT system.
 *  This comes from table 18 of JMcD's coordinate memo.
 *
 * stf_stt_offset: origin of STT in STF at nominal aimpoint for
 * the detector.  This value comes from table 19 of JMcD's coord
 * memo.
 */
/* The following quantities have units of mm */
#define MARX_DETECTOR_TYPE_PRIVATE_DATA \
   JDMVector_Type stt_lsi_offset; \
   JDMVector_Type stf_stt_offset;

#endif				       /* _MARX_CHANDRA_H_INCLUDED */
