#!/usr/bin/env slsh
require ("fits");
require ("histogram");

%require ("imageplot");

% On 20 Feb 2010, I asked Alexey about the uncertaintly on the fxy
% image values.  Here is his response:
% 
%    "I'd say, a 10% relative uncertainty is a good estimate, probably
%    slightly on the concervative side. Or, 10% relative and ~ +- 0.05
%    absolute, whichever is greater. (The absolute uncertainty is in
%    units of = tau for 0.67 keV at the given date and location).
private variable Rel_Tolerance = 0.05;
private variable Abs_Tolerance = 0.025;

define block_image (img)
{
   variable xgrid = [0:1023];
   variable ygrid = [0:1023];
   variable b = 1;
   variable type = _typeof(img);
   variable new_img = img;
   variable img1 = @img;
   variable new_img1 = @img1;

   while (1)
     {
	variable tmp_b = b*2;
	variable b2 = tmp_b*tmp_b;
	variable tmp_xgrid = [xgrid[0] : xgrid[-1] : tmp_b];
	variable tmp_ygrid = [ygrid[0] : ygrid[-1] : tmp_b];
	variable n = length (tmp_ygrid);
	variable tmp_img = type[n,n];
	variable i,j;
	_for i (0,n-1,1)
	  {
	     variable ii = [tmp_b*i:tmp_b*(i+1)-1];
	     _for j (0, n-1, 1)
	       {
		  variable jj = [tmp_b*j:tmp_b*(j+1)-1];
		  tmp_img[i, j] = sum(img[ii,jj])/b2;
		  img1[ii,jj] = tmp_img[i,j];
	       }
	  }

	i = where (fneqs (img1, img, Rel_Tolerance, Abs_Tolerance));
	if (length(i))
	  {
	     variable diff = abs(img1-img);
	     j = i[wherefirstmax (diff[i])];
	     i = j/1024;
	     j = j mod 1024;
	     vmessage ("diff=%g@[%d,%d] (%g vs %g) rdiff=%g", diff[i,j],i,j,
		       img[i,j], img1[i,j], diff[i,j]/img[i,j]);
	     break;
	  }

	b = tmp_b;
	new_img = tmp_img;
	new_img1 = @img1;
     }
   return new_img, b, new_img1;
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
	detnam,
	ccd_id,
	marxvers = 5.2,
	fxyblk,
     };

   _for (0, 9, 1)
     {
	variable ccd = ();
	variable vfile = sprintf ("%s[AXAF_CONTAM%d]", infile, ccd+1);
	variable t = fits_read_table (vfile,
				      ["component", "n_energy", "energy", "mu",
				       "n_time", "time", "tau0", "tau1", "fxy"]);
	if (array_shape(t.fxy)[0] != 1)
	  {
	     () = fprintf (stderr, "\
*** ERROR: This script requires a single component fxy image\n\
           Found %d components for CCD_ID=%d\n",
			   array_shape(t.fxy)[0], ccd);
	     exit (1);
	  }

	keywords.detnam = "ACIS-$ccd"$;
	keywords.ccd_id = ccd;
	variable b, img, ubimg;
	(img, b, ubimg) = block_image (t.fxy[0,*,*]);
	vmessage ("b=%d", b);
	keywords.fxyblk = b;

#ifexists imageplot
	variable zmin = _min(min(img), min(ubimg));
	variable zmax = _max(min(img), max(ubimg));
	imageplot ("contam_img_${ccd}_old.png"$, t.fxy[0,*,*], "ds9sls";
		   zmin=zmin, zmax=zmax);
	imageplot ("contam_img_${ccd}_new.png"$, ubimg, "ds9sls";
		   zmin=zmin, zmax=zmax);
#endif

	reshape (img, [1, array_shape (img)]);
	t.fxy = img;
	fits_write_binary_table (fp, "ACIS${ccd}_CONTAM"$, t, keywords);
     }

   fits_close_file (fp);
}
