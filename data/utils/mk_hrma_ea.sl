require ("fits");

private variable Marx_Cmd = "../../src/objs/marx";
private variable Marx_Args= "SourceRA=0 SourceDec=0"
  + " SpectrumType=FLAT  RA_Nom=0 Dec_Nom=0"
  + " Roll_Nom=0 DitherModel=NONE";
Marx_Args += " ExposureTime=0";
Marx_Args += " GratingType=NONE Verbose=no";
Marx_Args += " DetectorType=NONE";
Marx_Args += " SourceFlux=1.0";
% If HRMA_Use_Scale_Factors is yes, then prefix the output filename
% with "scaled".

private variable Use_Scale_Factors = 0;
private variable Prefix = Use_Scale_Factors ? "scaled" : "";
Marx_Args += " HRMA_Use_Scale_Factors=" + (Use_Scale_Factors ? "yes" : "no");

putenv ("MARX_DATA_DIR=..");

private define run_marx_sim (numrays, en_lo, en_hi, file, shell)
{
   variable shutters;

   if (shell == 1346)
     shutters = "Shutters1=0000 Shutters3=0000 Shutters4=0000 Shutters6=0000";
   else
     {
	shutters = "Shutters1=1111 Shutters3=1111 Shutters4=1111 Shutters6=1111";
	shutters += sprintf (" Shutters%d=0000", shell);
     }

   shutters += sprintf (" NumRays=%d", numrays);

   variable cmd = sprintf ("%s %s %s OutputDir=\"|./marxpipe\" MinEnergy=%g MaxEnergy=%g >> %s",
			   Marx_Cmd, Marx_Args, shutters, en_lo, en_hi, file);

   vmessage ("Running %s", cmd);
   if (0 != system (cmd))
     {
	verror ("%s failed", cmd);
     }
}

define slsh_main ()
{
   variable de = 0.003;		       %  if changed, update mkcorr.sl
   variable en_lo = [0.03:12.0:de];
   variable en_hi = en_lo + de;

   () = system ("cp ../../par/marx.par .");

   foreach ([1, 3, 4, 6])
     {
	variable shell = ();
	variable i;
	variable Out_Dat = "${Prefix}_hrma_ea_${shell}.dat"$;
	() = remove (Out_Dat);
	variable numrays = -100000;
	for (i = 0; i < length (en_lo); i++)
	  {
	     variable fp = fopen (Out_Dat, "a");
	     if (en_lo[i] > 6.0)
	       numrays = 5000000;
	     () = fprintf (fp, "%g\t%g\t", en_lo[i], en_hi[i]);
	     () = fclose (fp);
	     run_marx_sim (numrays, en_lo[i], en_hi[i], Out_Dat, shell);
	  }
     }
}

