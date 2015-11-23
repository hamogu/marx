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
#include "config.h"

#include <stdio.h>


#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "jdmath.h"

static void error (char *s)
{
   fprintf (stderr, "\n\n****ERROR: %s\n", s);
   fputs ("\
\n\
An error has been encountered while trying to determine information about\n\
your computer system.  You may have to edit the file jdmath.h to port the\n\
library to your system.\n\n", 
	  stderr);
	  
   exit (-1);
}


int main (int argc, char **argv)
{
   union
     {
	uint32 n;
	float32 f;
     }
   u;
   unsigned short s;
   
   (void) argc; (void) argv;
   
   fprintf (stdout, "#ifndef INT16_DEFINED\n#define INT16_DEFINED\n");
   
   if (sizeof (int16) != 2) 
     error ("Sizeof int16 != 2");
   if (sizeof (int32) != 4)
     error ("Sizeof int32 != 4");
   if (sizeof (float32) != 4)
     error ("Sizeof float32 != 4");
   if (sizeof (float64) != 8)
     error ("Sizeof float64 != 8");

   /* Now determine floating point format and byte swapping */
   
   s = 0x1234;
   if (*(char *) &s == 0x34)
     {
	fprintf (stdout, "#define NEEDS_BYTE_SWAP\n");
     }
   
   /* Now determine the floating point style */
   
   u.f = 1.2345678f;
   if (u.n == 0x3F9E0651)
     {
	fprintf (stdout, "#define USES_IEEE_FLOAT\n");
     }
   else if (u.n == 0x0651409E)
     {
	fprintf (stdout, "#define USES_VAX_FLOAT\n");
     }
   else 
     {
        error ("Unknown floating point format.  Not supported.");
     }
   
   fprintf (stdout, "#endif /* INT16_DEFINED */\n");
   return 0;
}

