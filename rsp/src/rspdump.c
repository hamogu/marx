/* -*- mode: C; mode: fold; -*- */
#include <stdio.h>

#include <stdlib.h>
#include <string.h>

#include "cfits.h"

#include "argcargv.h"

static int dump_response_file (char *file, /*{{{*/
			       int summary,
			       double energy,
			       int min_channel,
			       int max_channel,
			       int dump_stats)
{
   long num_rows, i;
   CFits_Response_File_Type *rsp;
   CFits_Type *ft;

   rsp = cfits_open_rmf_response_file (file, 0);
   if (rsp == NULL) return -1;

   num_rows = rsp->num_rows;
   ft = rsp->ft;

   fprintf (stdout, "#NUM ROWS: %ld\n", num_rows);

   if (dump_stats)
     {
	fprintf (stdout, "#%19s %10s %10s %10s\n",
		 "Energy", "Norm", "Peak-Chan", "Mean");
     }

   if (summary)
     {
	(void) fprintf (stdout, "#%19s\t%s\t%s\t%s\t%s\t%20s\t%s\n",
		       "ENERG_LO", "NGRPS", "MINCHAN", "MAXCHAN",
			"XXX", "INTEG_RESP", "MEANCHAN"
		       );
     }

   for (i = 1; i <= num_rows; i++)
     {
	CFits_Response_Vector_Type rv;
	CFits_Response_Element_Type *elem, *elem_max;
	double emin, emax, sum, sum1;
	unsigned int min_chan, max_chan;
	double max_chan_value;
	unsigned int max_val_chan;

	if (-1 == cfits_read_column_doubles (ft, rsp->energ_lo_col, i, 1, &emin, 1))
	  break;

	if (-1 == cfits_read_column_doubles (ft, rsp->energ_hi_col, i, 1, &emax, 1))
	  break;

	if (energy >= 0.0)
	  {
	     if (emax <= energy)
	       continue;
	     if (emin > energy)
	       return 0;
	  }

	if (-1 == cfits_read_response_vector (rsp, i, &rv))
	  {
	     cfits_close_response_file (rsp);
	     return -1;
	  }

	if (rv.num_grps == 0)
	  {
	     if (dump_stats)
	       fprintf (stdout, "%20.14e\t0\t-1\t0\n", 0.5*(emin+emax));
	     else if (summary)
	       fprintf (stdout, "%20.14e\t0\t-1\t-1\t-1\t0\t-1\n", emin);
	     else
	       fprintf (stdout, "\n#ENERG_LO: %20.14e\tENERG_HI: %20.14e\tN_GRP: %d\n",
			emin, emax, rv.num_grps);
	     continue;
	  }

	elem = rv.elem;
	elem_max = elem + rv.num_grps;

	min_chan = elem->first_channel;
	max_chan = (elem_max - 1)->first_channel + (elem_max - 1)->num_channels - 1;

	if (min_channel >= 0)
	  {
	     if (max_chan < (unsigned int) min_channel)
	       continue;

	     if ((unsigned int) min_channel > min_chan)
	       min_chan = min_channel;
	  }

	if (max_channel >= 0)
	  {
	     if (min_chan > (unsigned int) max_channel)
	       return 0;
	     if (max_chan > (unsigned int) max_channel)
	       max_chan = max_channel;
	  }

	if (dump_stats)
	  fprintf (stdout, "%20.14e", 0.5*(emin+emax));
	else if (summary)
	  fprintf (stdout, "%20.14e\t%d\t%4u\t%4u", emin, rv.num_grps,
		   min_chan, max_chan);
	else
	  fprintf (stdout, "\n#ENERG_LO: %20.14e\tENERG_HI: %20.14e\tN_GRP: %d\n",
		   emin, emax, rv.num_grps);

	sum = 0.0;
	sum1 = 0.0;
	max_val_chan = 0;
	max_chan_value = -1.0;

	while (elem < elem_max)
	  {
	     double *val, *val_max;
	     unsigned int chan_num;

	     chan_num = elem->first_channel;

	     if (chan_num + (elem->num_channels - 1) < min_chan)
	       {
		  elem++;
		  continue;
	       }

	     if (chan_num > max_chan)
	       break;

	     val = elem->response;
	     val_max = elem->response + elem->num_channels;

	     if (summary || dump_stats)
	       {
		  while (val < val_max)
		    {
		       if (chan_num > max_chan)
			 break;

		       if (chan_num >= min_chan)
			 {
			    if (*val > max_chan_value)
			      {
				 max_chan_value = *val;
				 max_val_chan = chan_num;
			      }
			    sum += *val;
			    sum1 += *val * chan_num;
			 }

		       val++;
		       chan_num++;
		    }
	       }
	     else
	       {
		  while (val < val_max)
		    {
		       if (chan_num > max_chan)
			 break;

		       if (chan_num >= min_chan)
			 fprintf (stdout, "\t%04u\t%20.14e\n", chan_num, *val);

		       val++;
		       chan_num++;
		    }
	       }

	     elem++;
	  }

	if (dump_stats)
	  {
	     fprintf (stdout, " %10g %10u %10.4e\n",
		      sum, max_val_chan, ((sum == 0)?0.0:sum1/sum));
	  }
	else if (summary) fprintf (stdout, "\t%u\t%20.14e\t%g\n", max_val_chan, sum,
			      ((sum == 0)?0.0:sum1/sum));
	else
	  fputc ('\n', stdout);
     }

   cfits_close_response_file (rsp);

   return 0;
}

/*}}}*/

static double Opt_Energy = -1.0;
static int Opt_Summary = 0;
static int Opt_Min_Chan = -1;
static int Opt_Max_Chan = -1;
static int Opt_Stats = 0;

static ArgcArgv_Type ArgcArgv_Table [] = /*{{{*/
{
   {"--integrate", ARGCARGV_BOOLEAN, (long) &Opt_Summary, "Integrate across channels"},
   {"--stats", ARGCARGV_BOOLEAN, (long) &Opt_Stats, "Compute moments of the distribution"},
   {"--energy", ARGCARGV_DOUBLE, (long) &Opt_Energy, "Specifies energy of channel to dump"},
   {"--min-chan", ARGCARGV_INTEGER, (long) &Opt_Min_Chan, "Specifies minimum channel"},
   {"--max-chan", ARGCARGV_INTEGER, (long) &Opt_Max_Chan, "Specifies maximum channel"},
   {NULL, 0, (long) NULL, NULL}
};

/*}}}*/

static void usage (char *pgm) /*{{{*/
{
   ArgcArgv_Type *a;

   fprintf (stderr, "Usage:\n%s [--integrate] [--energy E] [--min-chan f] [--max-chan l] rsp-file\n",
	    pgm);

   a = ArgcArgv_Table;
   while (a->name != NULL)
     {
	if (a->help != NULL)
	  {
	     fprintf (stderr, "\t%s\t\t%s\n", a->name, a->help);
	  }
	a++;
     }
   exit (1);
}

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   char *rsp_file;
   char *pgm;

   pgm = argv[0];
   argv++; argc--;

   if ((-1 == argcargv (&argc, &argv, ArgcArgv_Table))
       || (argc != 1))
     usage (pgm);

   rsp_file = argv[argc - 1];

   if (-1 == dump_response_file (rsp_file, Opt_Summary, Opt_Energy, Opt_Min_Chan, Opt_Max_Chan, Opt_Stats))
     return 1;

   return 0;
}

/*}}}*/
