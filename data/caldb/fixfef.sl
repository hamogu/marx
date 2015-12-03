#!/usr/bin/env slsh
require ("fits");

define slsh_main ()
{
   if (__argc != 3)
     {
	() = fprintf (stderr, "Usage: %s <infile> <outfile>\n", __argv[0]);
	exit (1);
     }
   variable infile = __argv[1];
   variable outfile = __argv[2];

   % slang's sort is stable.
   variable tbl = fits_read_table (infile);
   struct_filter (tbl, array_sort(tbl.energy));
   struct_filter (tbl, array_sort(tbl.regnum));

   variable keys = struct
     {
	function = fits_read_key (infile, "FUNCTION"),
     };

   fits_write_binary_table (outfile, "FUNCTION", tbl, keys,
			    sprintf ("This file differs from %s in that it is sorted differently",
				     infile));
}
