#!/usr/bin/env slsh

require ("readascii");
require ("fits");

define slsh_main ()
{
   if (__argc != 3)
     {
	() = fprintf (stderr, "Usage: %s <infile.ascii> <outfile.fits>\n", __argv[0]);
	exit (1);
     }
   variable infile = __argv[1];
   variable outfile = __argv[2];

   variable fp = fopen (infile, "r");
   variable dhistory, hist = struct {history={}};
   while (-1 != fgets (&dhistory, fp))
     {
	if (dhistory[0] != '#')
	  break;
	list_append (hist.history, strtrim_end(dhistory[[1:]]));
	%history += strtrim_end (dhistory[[1:]]);
     }
   () = fclose (fp);

   variable d = struct
     {
	energy,
	% effm2,
	effm1, eff0, eff1, eff2, eff3, eff4, eff5,
	eff6, eff7, eff8, eff9, eff10, eff11, eff12,
	eff13
     };

   () = readascii (infile,
		   &d.energy,
		   %&d.effm2,
		   &d.effm1, &d.eff0, &d.eff1,
		   &d.eff2, &d.eff3, &d.eff4, &d.eff5,
		   &d.eff6, &d.eff7, &d.eff8, &d.eff9,
		   &d.eff10, &d.eff11, &d.eff12
		   , &d.eff13
		   ; type="%f");

   %if (d.effm2 == NULL) d.effm2 = 0*d.eff0;

   % The values need to be sorted on energy
   struct_filter (d, [length(d.energy)-1:0:-1]);

   variable keys = struct
     {
	hduname="IXO_GREFF", tunit1="keV",
     };

   fits_write_binary_table (outfile, "IXO_GREFF", d, keys, hist);
}
