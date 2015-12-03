#!/usr/bin/env isis-script
require ("fits");
require ("linearfit");

if (__argc != 5)
{
   () = fprintf (stderr, "Usage: %s FEF CCDID CHIPX CHIPY\n", __argv[0]);
   exit (1);
}

private variable File = __argv[1];
private variable Chip_Id = integer (__argv[2]);
private variable Chip_X = atof (__argv[3]);
private variable Chip_Y = atof (__argv[4]);

private define compute_fef (t, p, kth)
{
   variable i, imax = 10;

   variable y = Double_Type[length (kth),length(p)];

   for (i = 1; i <= imax; i++)
     {
	variable sigma = get_struct_field (t, sprintf ("g%d_fwhm", i))[kth]
	  * 0.4246609001440095;
	variable pos = get_struct_field (t, sprintf ("g%d_pos", i))[kth];
	variable amp = get_struct_field (t, sprintf ("g%d_ampl", i))[kth];
	variable k0=0;
	vmessage ("[%d %d]: fwhm=%S, pos=%S, ampl=%S (regnum=%d)", i, k0, sigma[k0], pos[k0], amp[k0],
		  t.regnum[kth[0]]);

	variable j;

	_for (0, length (kth)-1,1)
	  {
	     j = ();
	     if (sigma[j] != 0)
	       {
		  variable x = (p - pos[j])/sigma[j];
		  %y[j,*] += amp[j] * exp (-0.5*x*x)/(sqrt(2*PI) * sigma[j]);
		  % The FEF uses unnormalized gaussians
		  y[j,*] += amp[j] * exp (-0.5*x*x);
	       }
	  }
     }

   variable e = t.energy[kth];
   variable c = t.channel[kth];
   e = interpol (p, c, e);
   return (e,y);
   %return (p,y);
}

private define read_fef (file, ccdid, x, y)
{
   variable t = fits_read_table (file);
   variable i = where ((t.chipx_lo <= x)
		       and (t.chipx_hi >= x)
		       and (t.chipy_lo <= y)
		       and (t.chipy_hi >= y)
		       and (t.ccd_id == ccdid));

   variable p = [-4096:4096];
   variable e;
   (e,y) = compute_fef (t, p, i);
   variable id = plot_open ("out.ps/CPS");

   ylabel ("PHA");
   xlabel ("ENERGY");
   xrange (0.0, 10.0);
   connect_points (1);
   %variable m, b;
   pointstyle(3);
   %ylog;
   variable m, b;
   (m, b) = linear_fit (t.energy[i], t.channel[i]);
   plot (t.energy[i], t.channel[i] - (m*t.energy[i]+b));
   %plot (t.energy[i], 1000.0*t.energy[i]/t.channel[i]);
   plot_pause ();
   ylin;

   connect_points(1);
   variable j;
   for (j = 0; j < length(i); j++)
     {
	variable en = t.energy[i[j]];
	variable yy = y[j,*];
	variable de = shift (e,1)-e;
	variable mean_en = sum(yy * e *de)/sum(yy*de);
	variable mean_chan = sum (p*yy)/sum(yy);
	variable k = wherefirstmax(yy);
	title (sprintf ("En=%g, Mean En=%g, Peak En=%g, Mean Chan=%g=%g keV",
			en, mean_en, e[k],
			mean_chan,
			interpol (mean_chan, t.channel[i], t.energy[i])[0]));

	xrange (_max(e[k]*0.95,0), e[k]*1.1);
	%ylog(); yrange (max(yy)/1e6,);
	variable ipos = where(p>0);
	yy = cumsum(yy[ipos]); yy/=yy[-1];
	plot (e[ipos],yy);
	plot_pause;
     }
   plot_close (id);
}

read_fef (File, Chip_Id, Chip_X, Chip_Y);
