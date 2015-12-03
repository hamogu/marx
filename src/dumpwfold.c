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
#include <stdio.h>

#include "marx.h"

int main (int argc, char **argv)
{
   char *file;
   /* Marx_WFold_Table_Type *table; */

   if (argc != 2)
     {
	fprintf (stderr, "Usage: %s FILENAME\n", argv[0]);
	return 1;
     }

   file = argv[1];
#if 0
   if (NULL == (table = marx_read_wfold_file (file)))
     {
	return 1;
     }

   marx_free_wfold_table (table);
#endif
   return marx_wfold_dump_file (file);
}

