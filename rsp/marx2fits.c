#include <stdio.h>

#include <stdlib.h>
#include <string.h>

#include <cfits.h>

#include <marx.h>
#include <jdmath.h>

#include "argcargv.h"

int main (int argc, char **argv)
{
   char *file;
   CFits_Type *ft;

   (void) argc;
   file = argv[1];
   
   ft = cfits_create_file (file, 1);
   if (ft == NULL)
     return 1;
   
   (void) cfits_init_btable_extension (ft, "SPEC RESPONSE");
   
   cfits_close_file (ft);
   return 0;
}

   
   
