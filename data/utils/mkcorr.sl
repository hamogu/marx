require ("fits");
require ("gslinterp");
require ("histogram");
require ("xfig");

private variable Flux = 1.0;	       %  simulation SourceFlux value
private variable Prefix = "scaled";

define read_file (shell)
{
   variable file = "${Prefix}_hrma_ea_$shell.dat"$;
   variable e, n, t;
   () = readascii (file, &e, &n, &t; cols=[1,4,5]);
   %(e,n,t) = readcol (file, 1, 4, 5);
   t *= Flux;
   return e, n/t, sqrt(n)/t;
}

define read_arf (shell)
{
   variable file = sprintf ("hrma%d.arf", shell);
   variable t = fits_read_table (file);
   return (0.5*(t.energ_lo+t.energ_hi), t.specresp);
}

define compute_ratio (shell)
{
   variable e,a,da;
   variable e1, a1;

   (e,a,da) = read_file (shell);
   (e1, a1) = read_arf (shell);

   a1 = interp_linear (e, e1, a1);
   variable i = where ((a != 0) and (a1 != 0));
   return e[i], a1[i]/a[i], da[i]/a1[i];
}

define prune (x, y, ysigma, tol)
{
   variable min_slope, max_slope;
   variable len = length (x);
   variable p = Int_Type[len];
   variable n = 0, i;
   variable infinity = 1e38;

   p[n] = 0;
   n++;
   variable x0 = x[0];
   variable y0 = y[0];

   min_slope = infinity;
   max_slope = -infinity;

   i = 1;
   while (i < len)
     {
	variable y1 = y[i];
	variable x1 = x[i];
	variable dx = x1 - x0;
	variable dy = y1 - y0;
	if (x1 == x0)
	  {
	     i++;
	     continue;
	  }
	variable m = dy/dx;
	if ((m < min_slope) or (m > max_slope))
	  {
	     i--;
	     p[n] = i; n++;
	     x0 = x[i];
	     y0 = y[i];
	     min_slope = -infinity;
	     max_slope = infinity;
	     i++;
	     continue;
	  }

	variable dm = abs ((y1*tol)/dx);
	if (ysigma != NULL)
	  dm += 0.5*abs (ysigma[i]/dx);

	if (m + dm < max_slope)
	  max_slope = m + dm;
	if (m - dm > min_slope)
	  min_slope = m - dm;

	i++;
     }
   if (n == len)
     n--;
   else
     p[n] = len-1;
   return p[[1:n]];
}

define smooth (r)
{
   variable r0, r1, r2;
   variable n = length (r);
   variable new_r = Double_Type [n];
   r0 = r[0];
   new_r[0] = r0;
   variable i;
   _for (1, n-2, 1)
     {
	i = ();
	r1 = r[i];
	r2 = r[i+1];
	new_r[i] = (r0 + r1 + r2)/3.0;
	r0 = r1;
     }
   new_r[-1] = r[-1];
   return new_r;
}

define smooth (r)
{
   variable r0, r1, r2;
   variable n = length (r);
   variable new_r = Double_Type [n];
   r0 = r[0];
   new_r[0] = r0;
   variable i;
   _for (1, n-2, 1)
     {
	i = ();
	r1 = r[i];
	r2 = r[i+1];
	if (0)
	  new_r[i] = [r0, r1, r2][array_sort([r0,r1,r2])[1]];
	else
	  new_r[i] = (r0 + r1 + r2)/3.0;
	r0 = r1;
     }
   new_r[-1] = r[-1];
   return new_r;
}

define process_shell (shell, nsmooth)
{
   variable e, r, dr;
   variable i;
   (e,r,dr) = compute_ratio (shell);
   variable e1, e2, r1, r2;

   variable cutoff;
   variable cutoffs = [0.0, 6.8, 0.0, 6.8, 8.0, 0.0, 10.0];

   if (nsmooth < 0)
     return e, r;

   loop (nsmooth)
     r = smooth (r);

   cutoff = cutoffs[shell];
   i = where (e <= cutoff);
   e1 = e[i]; r1 = r[i];
   i = prune (e1, r1, NULL, 0.005/2);
   e1 = e1[i]; r1 = r1[i];

   i = where (e > cutoff);
   e2 = [min(e[i]):max(e[i]):0.05];
   % The 0.003 factor must match the bin width in mk_hrma_ea.sl
   r2 = hist1d_rebin (e2, e[i], r[i])*0.003/0.05;
   i = prune (e2, r2, NULL, 0.01/2);
   e2 = e2[i]; r2 = r2[i];

   vmessage ("reduction of %d/%d [%d,%d]",
	     length (e1)+length(e2), length(e), length (e1), length (e2));
   return ([e1, e2, [e2[-1]+1]], [r1, r2, [r2[-1]]]);
}

define plot_shell (shell, nsmooth)
{
   variable w = xfig_plot_new (15,7);
   w.title ("HRMA shell $shell"$);
   w.xlabel ("Energy [keV]");
   w.ylabel ("Correction Factor");

   variable e1, r1, e, r, dr;
   (e, r, dr) = compute_ratio (shell);
   (e1,r1) = process_shell (shell, nsmooth);
   w.plot (e,r;color="blue");
   w.plot (e1,r1;color="red");

   variable w1 = xfig_plot_new (15,3);
   r1 = interp_linear (e, e1, r1);
   w1.plot (e, (r - r1)/dr; color="blue");
   w1.xlabel ("Energy [keV]");
   w1.ylabel ("Residual");
   return xfig_multiplot (w,w1);
}

define write_shell (shell, nsmooth)
{
   variable fp = fopen ("${Prefix}corr_$shell.tbl"$, "w");
   () = array_map (Int_Type, &fprintf, fp, "%g %g\n", process_shell (shell, nsmooth));
   () = fclose (fp);
}

define slsh_main ()
{
   variable wlist = {};
   variable nsmooth = 2;
   foreach ([1,3,4,6])
     {
	variable shell = ();
	write_shell (shell, nsmooth);
	list_insert (wlist, plot_shell (shell, nsmooth));
     }

   xfig_new_vbox_compound
     (
	 xfig_new_hbox_compound (wlist[0], wlist[1], 1),
	 xfig_new_hbox_compound (wlist[2], wlist[3], 1),
	 2
     ).render ("${Prefix}corr.pdf"$);
}


