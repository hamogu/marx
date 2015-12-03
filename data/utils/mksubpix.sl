#!/usr/bin/env slsh
require ("fits");

private variable Grade_Map = 
  [
   0,   1,   2,   5,   1,   1,   5,   7,   3,   5,   6,   6,   3,   5,   7,   7,
   4,   4,   6,   7,   5,   5,   6,   7,  -1,   7,   7,   7,   7,   7,   7,   7,
   1,   1,   2,   5,   1,   1,   5,   7,   5,   7,   7,   7,   5,   7,   7,   7,
   4,   4,   6,   7,   5,   5,   6,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   2,   2,  -1,   7,   2,   2,   7,   7,   6,   7,   7,   7,   6,   7,   7,   7,
   6,   6,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   5,   5,   7,   7,   5,   5,   7,   7,   6,   7,   7,  -1,   6,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   1,   1,   2,   5,   1,   1,   5,   7,   3,   5,   6,   6,   3,   5,   7,   7,
   5,   5,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   1,   1,   2,   5,   1,   1,   5,   7,   5,   7,   7,   7,   5,   7,   7,   7,
   5,   5,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   5,   5,   7,   7,   5,   5,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   6,   6,   7,   7,   7,   7,  -1,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,  -1,
  ];

private define read_file (file, ccd)
{
   variable t = fits_read_table (sprintf ("%s+%d", file, ccd+1));
   variable grade = Grade_Map[t.fltgrade];
   struct_filter (t, where(63 mod (grade+2)); dim=0);

   % Impose a fixed order so the units can be specified.
   return struct
     {
	fltgrade = t.fltgrade,
	npoints = t.npoints,
	energy = t.energy * 0.001f,
	chipx_offset = t.chipx_offset,
	chipy_offset = t.chipy_offset,
     };
}


define slsh_main ()
{
   if (__argc != 3)
     {
	() = fprintf (stderr, "Usage: %s <in-caldb-file> <out-marx-file>\n", __argv[0]);
	exit (1);
     }
   variable infile = __argv[1];
   variable outfile = __argv[2];

   variable fp = fits_open_file (outfile, "c");
   fits_create_image_hdu (fp, NULL, Char_Type, Int_Type[0]);
   fits_write_history (fp, sprintf ("This file was created using the marx script `%s %s %s`",
				    __argv[0], __argv[1], __argv[2]));

   variable keywords = struct
     {
	mission = "AXAF",
	instrum = "ACIS",
	marxvers = 5.0,
	detnam,
	tunit3 = "keV",
     };

   variable data;

   data = read_file (infile, 0);
   keywords.detnam = "ACIS-01234689";
   fits_write_binary_table (fp, "MARX_ACIS_SUBPIX_FI"$, data, keywords);

   data = read_file (infile, 5);
   keywords.detnam = "ACIS-57";
   fits_write_binary_table (fp, "MARX_ACIS_SUBPIX_BI"$, data, keywords);

   fits_close_file (fp);
}

