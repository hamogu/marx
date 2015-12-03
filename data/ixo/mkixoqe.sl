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
	dhistory = strtrim_end (dhistory);
	dhistory = strtrans (dhistory, "\t", " ");
	list_append (hist.history, strtrim_end(dhistory[[1:]]));
	%history += strtrim_end (dhistory[[1:]]);
     }
   () = fclose (fp);

   variable d = struct
     {
	energy, qe, filter,
     };
   
   () = readascii (infile,
		   &d.energy, &d.qe, &d.filter
		   ; type="%f", cols=[1,2,3]);
   d.energy *= 0.001f;

   variable keys = struct 
     {
	hduname="IXO_CCD_QE",
	tunit1 = "keV",
     };

   fits_write_binary_table (outfile, "IXO_CCD_QE", d, keys, hist);
}
