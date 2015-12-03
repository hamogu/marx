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
   double xpixel, ypixel, zpixel;		       /* pixel units */
   double px, py, pz;		       /* direction cosines */
   double energy;		       /* KeV */
   double time;			       /* secs */
}
AcisSim_Ray_Type;

typedef struct
{
#define MAX_ISLAND_SIZE 7
   unsigned int island_size;
   double phas [MAX_ISLAND_SIZE * MAX_ISLAND_SIZE];   /* KeV */

   double x, y;			       /* x, y pixel location of center
					* pixel.  For example, if x is 12.5
					* then the center pixel is number 12
					* and the event happened at its center.
					*/
   double radius;		       /* num pixels */
   double charge_fraction;
}
AcisSim_Pixel_Island_Type;

typedef struct
{
   AcisSim_Pixel_Island_Type regular_island;
   AcisSim_Pixel_Island_Type fluoresc_island;
   unsigned int flags;
#define REGULAR_EVENT_OK 1
#define FLUOR_EVENT_OK 2

   int did_fluoresc;		       /* Non-zero if it fluoresced.  Note,
					* this does not mean the fluorescent
					* event is good.
					*/
}
AcisSim_Pixel_Event_Type;

extern int acissim_process_ray (AcisSim_Ray_Type *ray, AcisSim_Pixel_Event_Type *);
extern int acissim_init (Param_File_Type *);
