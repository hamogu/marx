/*
    This file is part of MARX

    Copyright (C) 2002-2018 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research
    under contract SV1-61010 from the Smithsonian Institution.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
/* This file uses MARX's internal model to output a PCAD file. */
#include "config.h"
#include "marx-feat.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <sys/types.h>
#include <time.h>

#include <marx.h>
#include <jdfits.h>

#include "argcargv.h"
#include "argcargv.c"

static char *Marx_Dir;
static char MarxAsp_Pgm[80];
static char *Program_Name = "marxasp";

static int Simulation_Detector_Type = 0;
static int Simulation_Grating_Type = 0;

static char *DetectorType;
static char *GratingType;
static char *Dither_Model;
static char *SourceType;

static char Todays_Date [64];	       /*  */
static double Sim_X;
static double Sim_Y;
static double Sim_Z;
static double Focal_Length;
static double Time_Start;
static double Time_Stop;
static double Delta_Time = 0.256;      /* seconds */
static double Nominal_RA;
static double Nominal_Dec;
static double Nominal_Roll;
static double Nominal_Roll_In_Radians;
static double ACIS_Exposure_Time;
static double ACIS_Frame_Transfer_Time;
static double TimeDel;
static double DT_Cor;
static double Exposure;

static double DetOffset_X;
static double DetOffset_Y;
static double DetOffset_Z;

static double Ra_Amp;
static double Dec_Amp;
static double Roll_Amp;
static double Ra_Period;
static double Dec_Period;
static double Roll_Period;
static double Ra_Phase;
static double Dec_Phase;
static double Roll_Phase;
#if 0
static double RA_Sigma;
static double Dec_Sigma;
static double Roll_Sigma;
#endif

/* The following unit vectors form an orthonormal basis */
static JDMVector_Type Nominal_Pointing;
static JDMVector_Type RA_Hat;
static JDMVector_Type Dec_Hat;

static int Simulation_Detector_Type;   /* bitmapped */
#define DETECTOR_NONE	0x00
#define DETECTOR_ACIS_S	0x01
#define DETECTOR_ACIS_I 0x02
#define DETECTOR_ACIS	(DETECTOR_ACIS_I|DETECTOR_ACIS_S)
#define DETECTOR_HRC_S	0x04
#define DETECTOR_HRC_I	0x08
#define DETECTOR_HRC	(DETECTOR_HRC_I|DETECTOR_HRC_S)

static char *Instrum_Name;
static char *HDU_Class;
static char *HDU_Class1;
static char *HDU_Class2;
static char *Content_Hdr;
static char *HDU_Name_Hdr;
static Param_File_Type *Obs_Par_Parms;

typedef struct
{
   unsigned int location;	       /* a bitmapped quantity */
#define FULL_COMPONENT 1
#define SHORT_COMPONENT        2

   char *keyword;
   int type;
#define H_PINT	1
#define H_PFLT	2
#define H_STR	3
#define H_PSTR	4
#define H_COM	5
#define H_ENV	6
/* #define H_INT	7 */
#define H_SMARX	8		       /* read from parameter file */
#define H_DMARX	8		       /* read from parameter file */
#define H_FILE	10   		       /* read from a header info file */
#define H_LOG	11		       /* logical */
   void *value;
   char *comment;
}
Fits_Header_Table_Type;

static int Int_Value_0 = 0;
static int Int_Value_1 = 1;
static int Int_Value_1024 = 1024;

static Fits_Header_Table_Type CC_NULL_Component [] =
{
     {3,"COMMENT",	H_COM,	NULL,	"\n------- Configuration Control Component -------\n\n"},
     {3,"ORIGIN",	H_STR,	"ASC",		NULL},
     {3,"CREATOR",	H_STR,	MarxAsp_Pgm, "Program creating this file"},
     {3,"HDUNAME",	H_PSTR,	&HDU_Name_Hdr,		NULL},
     {3,"HDUDOC",	H_STR,	"ASC-FITS-1.4: McDowell, ROTS: ASC FIts FIle Designers Guide",		NULL},
     {3,"HDUVERS",	H_STR,	"1.0.0",		NULL},

     {3,"HDUCLASS",	H_PSTR,	&HDU_Class,		NULL},
     {3,"HDUCLAS1",	H_PSTR,	&HDU_Class1,		NULL},
     {3,"HDUCLAS2",	H_PSTR,	&HDU_Class2,		NULL},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type CC_Component [] =
{
     {3,"COMMENT",	H_COM,	NULL,	"\n------- Configuration Control Component -------\n\n"},
     {3,"ORIGIN",	H_STR,	"ASC",		NULL},
     {3,"CREATOR",	H_STR,	MarxAsp_Pgm, "Program creating this file"},
     {1,"REVISION",	H_FILE,	NULL,		"Processing revision"},
     {3,"CONTENT",	H_PSTR,	&Content_Hdr,	NULL},

     {3,"HDUNAME",	H_PSTR,	&HDU_Name_Hdr,		NULL},
     {1,"HDUSPEC",	H_STR,	"ASPSOL ICD V2.4",		NULL},
     {3,"HDUDOC",	H_STR,	"ASC-FITS-1.4",		NULL},
     {3,"HDUVERS",	H_STR,	"1.0.0",		NULL},

     {3,"HDUCLASS",	H_PSTR,	&HDU_Class,		NULL},
     {3,"HDUCLAS1",	H_PSTR,	&HDU_Class1,		NULL},
     {3,"HDUCLAS2",	H_PSTR,	&HDU_Class2,		NULL},

     {1,"LONGSTRN",	H_STR, "OGIP 1.0",	"Unofficial Convention"},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Timing_Component [] =
{
     {3,"COMMENT",	H_COM,	NULL,	"\n-------- Timing Component -------\n\n"},
     {3,"DATE",		H_STR,	Todays_Date,	"Date and time of file creation (UTC)"},
     {3,"DATE-OBS",	H_FILE,	NULL,		"TT, with clock correction if CLOCKAPP"},
     {3,"DATE-END",	H_FILE,	NULL,		"TT, with clock correction if CLOCKAPP"},
     {3,"TIMESYS",	H_STR,	"TT",		"AXAF time will be Terrestrial Time"},
     {3,"MJDREF",		H_FILE,	NULL,		"MJD of clock start"},

     {3,"TIMEZERO",	H_FILE,	NULL,	"Clock Correction"},
     {3,"TIMEUNIT",	H_STR,	"s",	"seconds"},
     {1,"BTIMNULL",	H_FILE,	NULL,	"Basic Time offset (s)"},
     {1,"BTIMRATE",	H_FILE,	NULL,	"Basic Time clock rate (s/VCDUcount)"},
     {1,"BTIMDRFT",	H_FILE,	NULL,	"Basic Time clock drift (s/VCDUcount^2)"},
     {1,"BTIMCORR",	H_FILE,	NULL,	"Correction applied to Basic Time rate (s)"},

     {1,"TIMEREF",	H_STR,	"LOCAL",	"Time is local for data"},
     {1,"TASSIGN",	H_STR,	"SATELLITE",	"Source of time assignment"},
     {3,"CLOCKAPP",	H_LOG,	&Int_Value_1,	"Clock correction applied"},

     {1,"TIERRELA",	H_FILE,	NULL,	"Short term clock stability"},
     {1,"TIERABSO",	H_FILE,	NULL,	"Absolute precision of clock correction"},

     {1,"TIMVERSN",	H_STR,	"ASC-FITS-1.1", "AXAF Fits design document"},

     {3,"TSTART",	H_FILE,	NULL,		"As in the TIME column: raw space craft clock"},
     {3,"TSTOP",	H_FILE,	NULL,		"  add TIMEZERO and MJDREF for absolute TT"},

     {1,"TIMEPIXR",	H_PINT,	&Int_Value_0,	"Time stamp refers to start of bin"},
     {3,"TIMEDEL", 	H_PFLT,	&Delta_Time,	"Time resolution of data in seconds"},
     {3,"DTASPSOL", 	H_PFLT,	&Delta_Time,	"Time resolution of data in seconds"},
     {3,"ASPTYPE", 	H_STR,	"kalman",	"Simulated by marxasp"},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_Timing_Component [] =
{
     {1,"STARTMJF", H_PINT, &Int_Value_0,	"Major frame count at start"},
     {1,"STARTMNF", H_PINT, &Int_Value_0,	"Minor frame count at start"},
     {1,"STARTOBT", H_PINT, &Int_Value_0,	"Onboard MET close to STARTMJF and STARTMNF"},
     {1,"STOPMJF", H_PINT, &Int_Value_0,	"Major frame count at stop"},
     {1,"STOPMNF", H_PINT, &Int_Value_0,	"Minor frame count at stop"},
     {0,NULL, 0, NULL, NULL}
};

static char *Canonical_Detnam;

static Fits_Header_Table_Type Obs_Info_Component [] =
{
     {3,"COMMENT",	H_COM, NULL, "\n------- Observation Information -------\n\n"},
     {1,"OBSERVER",	H_ENV,	"USER",		"Observer or PI"},
     {1,"TITLE",	H_FILE,	NULL,	"Title of Observation"},
     {3,"OBS_ID",	H_STR,	"0",		"Observation ID (*)"},
     {3,"MISSION",	H_STR,	"AXAF",	"Advanced X-ray Astrophysics Facility"},
     {3,"TELESCOP",	H_STR,	"CHANDRA",	"Telescope used"},
     {3,"INSTRUME",	H_PSTR,	&Instrum_Name,	NULL},
     {1,"DETNAM",	H_PSTR,	&Canonical_Detnam,	NULL},
     {1,"GRATING",	H_PSTR,	&GratingType,	"Grating type used"},
     {1,"OBS_MODE",	H_STR,	"POINTING",	"Observation mode"},

     {1,"SIM_X",	H_PFLT,	&Sim_X,		"SIM offset, mm"},
     {1,"SIM_Y",	H_PFLT,	&Sim_Y,		"SIM offset, mm"},
     {1,"SIM_Z",	H_PFLT,	&Sim_Z,		"SIM offset, mm"},
     {1,"FOC_LEN",	H_PFLT,	&Focal_Length,	NULL},
     {1,"ONTIME",	H_FILE,	NULL,		"Ontime in seconds"},
     {1,"LIVETIME",	H_FILE,	NULL,		"Livetime in seconds"},
     {1,"EXPOSURE",	H_PFLT,	&Exposure,		NULL},
     {1,"DTCOR",	H_PFLT,	&DT_Cor,		"Dead time correction factor"},
     {1,"OBJECT",	H_FILE, NULL,	"Source Name"},
     {1,"RA_NOM",	H_PFLT,	&Nominal_RA,		"Nominal RA in degrees"},
     {1,"DEC_NOM",	H_PFLT,	&Nominal_Dec,		"Nominal Dec in degrees"},
     {1,"ROLL_NOM",	H_PFLT,	&Nominal_Roll,		"Nominal Roll in degrees"},
     {1,"EQUINOX",	H_FILE,	NULL,		"Equinox"},
     {1,"RADECSYS",	H_STR,	"ICRS",		"WCS system"},
     {1,"DATACLAS",	H_STR,	"SIMULATED",	"File contains simulated data produced by MARX"},

     {0,NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_Obs_Info_Component [] =
{
   {1,"DATAMODE",	H_FILE,	NULL,	"telemetry mode"},
   {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Acis_Faint_Header_Keywords [] =
{
     {3,"COMMENT",	H_COM,	NULL, "\nAXAF FITS Event File: ACIS Level 1\n\n"},
     {3,"READMODE",	H_STR,	"TIMED",	"CCD exposure mode"},

     {3,"FIRSTROW",	H_PINT,	&Int_Value_1,	"Index of first row of CCD readout"},
     {3,"NROWS",	H_PINT,	&Int_Value_1024,	"Number of rows in readout"},
     {3, "EXPTIME",	H_PFLT,	(void *)&ACIS_Exposure_Time,	"Commanded exposure time in secs"},

     {3,"COMMENT",	H_COM,	NULL, "\nApplied event correction/flagging reference files\n\n"},

     {3,"BIASFIL0",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL1",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL2",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL3",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL4",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL5",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL6",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL7",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL8",	H_STR,	"/dev/null",	NULL},
     {3,"BIASFIL9",	H_STR,	"/dev/null",	NULL},
     {3,"BPIXFILE",	H_STR,	"/dev/null",	NULL},

     {3,"COMMENT",	H_COM,	NULL, "\nApplied event calibration/transform reference files/systems\n\n"},

#if 0
     {3,"ACSYSCHP",	H_STR,	"AXAF-ACIS-1.0","TBR"},
     {3,"ACSYSDFP",	H_STR,	"ASC-FP-STF-1.0","TBR"},
     {3,"ACSYSSKY",	H_STR,	"ASC-SKY-STF-1.0","TBR"},
#endif
     {3,"GAINFILE",	H_STR,	"/dev/null",	NULL},
     {3,"GRD_FILE",	H_STR,	"/dev/null",	NULL},
     {3,"GRD_SCHM",	H_STR,	"ACIS",	NULL},

     {0, NULL, 0, NULL, NULL}
};

static Fits_Header_Table_Type Hrc_Header_Keywords [] =
{
     {0,NULL, 0, NULL, NULL}
};

#if 0
static Fits_Header_Table_Type GoodTime_Header_Keywords [] =
{
     {2,"EXTVER",	H_PINT,	&Int_Value_1,	"Convention for HDUCLASS"},
#if 0
     {2,"HDUCLASS",	H_STR,	"OGIP",	"Convention for HDUCLASS"},
     {2,"HDUCLAS1",	H_STR,	"GTI",	"File contains good time intervals"},
     {2,"HDUCLAS2",	H_STR,	"STANDARD","File contains standard good time intervals"},
#endif
#if 0
     {2,"COMMENT",	H_COM,	NULL,		"\nData model support keywords"},
     {2,"CFIELDS",	H_PINT,	&Int_Value_1,"Number of ASC Table Columns"},
     {2,"CNAM1",	H_STR,	"TIME",	"Time"},

     {2,"CNC1",		H_PINT,	&Int_Value_2,"Number of FITS cols for ASC col 1"},
     {2,"CETYP1",	H_STR,	"R",	"Data is an interval"},
     {2,"CITYP1",	H_STR,	"[)",	"Interval is semi-open"},
     {2,"TDISP1",	H_STR,	"F20.6","Display format for FITS col 1"},
     {2,"TDISP2",	H_STR,	"F20.6","Display format for FITS col 2"},
#endif
#if 0
     {2,"MTYPE1",	H_STR,	"TIME", "Data Model Keyword"},
     {2,"MFORM1",	H_STR,	"START,STOP", "Data Model Keyword"},
     {2,"METYP1",	H_STR,	"R", "Data Model Keyword"},
#endif
     {0,NULL, 0, NULL, NULL}
};
#endif

static char *make_marx_filename (char *f)
{
   static char file [1024];
   unsigned int len;

   strcpy (file, Marx_Dir);
   if (0 != (len = strlen (file)))
     {
	if (file[len - 1] != '/')
	  {
	     file[len++] = '/';
	     file[len] = 0;
	  }
     }
   strcpy (file + len, f);
   return file;
}

static Param_Table_Type Parm_Table [] =
{
     {"SourceType",	PF_STRING_TYPE,		&SourceType},
     {"GratingType",	PF_STRING_TYPE,		&GratingType},
     {"DetectorType",	PF_STRING_TYPE,		&DetectorType},
     {"DetOffsetX",	PF_DOUBLE_TYPE,		&DetOffset_X},
     {"DetOffsetY",	PF_DOUBLE_TYPE,		&DetOffset_Y},
     {"DetOffsetZ",	PF_DOUBLE_TYPE,		&DetOffset_Z},
     {"DitherModel",	PF_STRING_TYPE,		&Dither_Model},
     /* {"MirrorType",	PF_STRING_TYPE,		&Mirror_Type}, */
     {"FocalLength",	PF_DOUBLE_TYPE,		&Focal_Length},

     /* {"DitherModel",		PF_STRING_TYPE, &Dither_Model}, */

     {"DitherAmp_RA",		PF_REAL_TYPE,	&Ra_Amp},
     {"DitherAmp_Dec",		PF_REAL_TYPE,	&Dec_Amp},
     {"DitherAmp_Roll",		PF_REAL_TYPE,	&Roll_Amp},

     {"DitherPeriod_RA",	PF_REAL_TYPE,	&Ra_Period},
     {"DitherPeriod_Dec",	PF_REAL_TYPE,	&Dec_Period},
     {"DitherPeriod_Roll",	PF_REAL_TYPE,	&Roll_Period},

     {"DitherPhase_RA",		PF_REAL_TYPE,	&Ra_Phase},
     {"DitherPhase_Dec",	PF_REAL_TYPE,	&Dec_Phase},
     {"DitherPhase_Roll",	PF_REAL_TYPE,	&Roll_Phase},
#if 0
     {"Roll_Nom",		PF_REAL_TYPE,	&Nominal_Roll},
#endif

     {NULL, 0, NULL}
};

static int setup_dither (void)
{
   double ra, dec;

   Ra_Amp *= PI/(180.0 * 3600);
   Dec_Amp *= PI/(180.0 * 3600);
   Roll_Amp *= PI/(180.0 * 3600);

   Nominal_Roll_In_Radians = Nominal_Roll * PI/180.0;
   ra = Nominal_RA * PI/180.0;
   dec = Nominal_Dec * PI/180.0;

   Nominal_Pointing = JDMv_spherical_to_vector (1.0, PI/2.0 - dec, ra);

   RA_Hat.x = -sin(ra);
   RA_Hat.y = cos (ra);
   RA_Hat.z = 0;

   Dec_Hat.x = -sin(dec) * cos (ra);
   Dec_Hat.y = -sin(dec) * sin (ra);
   Dec_Hat.z = cos (dec);

   return 0;
}

static int get_marx_pfile_info (Param_File_Type *pf)
{
   Marx_Detector_Type *dt;
   int used_dither = 0;

   if (-1 == pf_get_parameters (pf, Parm_Table))
     return -1;

   if (0 == strcmp (Dither_Model, "INTERNAL"))
     used_dither = 1;
   else
     {
	if (0 != strcmp (Dither_Model, "NONE"))
	  {
	     marx_error ("\
*** This simulation did not use the INTERNAL dither model.  Re-run the\n\
    simulation with DitherModel=INTERNAL.\n");
	     return -1;
	  }
     }

   if (used_dither && (0 == strcmp (SourceType, "SAOSAC")))
     used_dither = 0;

   if (used_dither == 0)
     {
	Ra_Amp = Dec_Amp = Roll_Amp = 0.0;
	/* RA_Sigma = Dec_Sigma = Roll_Sigma = 0.0; */
     }

   if (-1 == marx_set_data_directory ("$MARX_DATA_DIR"))
     return -1;

   Canonical_Detnam = DetectorType;

   if (0 == strcmp (DetectorType, "ACIS-S"))
     {
	Simulation_Detector_Type = DETECTOR_ACIS_S;
	Canonical_Detnam = "ACIS-456789";
     }
   else if (0 == strcmp (DetectorType, "ACIS-I"))
     {
	Canonical_Detnam = "ACIS-0123";
	Simulation_Detector_Type = DETECTOR_ACIS_I;
     }
   else if (0 == strcmp (DetectorType, "HRC-S"))
     Simulation_Detector_Type = DETECTOR_HRC_S;
   else if (0 == strcmp (DetectorType, "HRC-I"))
     Simulation_Detector_Type = DETECTOR_HRC_I;

  if (0 == strcmp (GratingType, "HETG"))
     Simulation_Grating_Type = 1;
   else if (0 == strcmp (GratingType, "LETG"))
     Simulation_Grating_Type = 2;
   else Simulation_Grating_Type = 0;

   if (NULL == (dt = marx_get_detector_info (DetectorType)))
     {
	pf_close_parameter_file (pf);
	marx_error ("*** DetectorType %s not supported.\n", DetectorType);
	return -1;
     }

   /* FIXME!!!  This may not be currently what is used in marx2fits!!! */
   Sim_X = DetOffset_X + dt->aimpoint_offset.x;
   Sim_Y = DetOffset_Y + dt->aimpoint_offset.y;
   Sim_Z = DetOffset_Z + dt->aimpoint_offset.z;

   Instrum_Name = "PCAD";

   if (Simulation_Detector_Type & DETECTOR_ACIS)
     {
	if ((-1 == pf_get_double (pf, "ACIS_Exposure_Time", &ACIS_Exposure_Time))
	    || (-1 == pf_get_double (pf, "ACIS_Frame_Transfer_Time", &ACIS_Frame_Transfer_Time)))
	  {
	     pf_close_parameter_file (pf);
	     return -1;
	  }
	if (ACIS_Exposure_Time > 0.0)
	  {
	     TimeDel = ACIS_Exposure_Time;
	     if (ACIS_Frame_Transfer_Time > 0.0)
	       TimeDel += ACIS_Frame_Transfer_Time;

	     DT_Cor = ACIS_Exposure_Time/TimeDel;
	  }
	else
	  {
	     TimeDel = 0.0;
	     DT_Cor = 1.0;
	  }
     }

   return 0;
}

static Param_Table_Type ObsPar_Parm_Table [] =
{
     {"RA_Nom",			PF_REAL_TYPE,	&Nominal_RA},
     {"Dec_Nom",		PF_REAL_TYPE,	&Nominal_Dec},
     {"Roll_Nom",		PF_REAL_TYPE,	&Nominal_Roll},
     {"TSTART",			PF_REAL_TYPE,	&Time_Start},
     {"TSTOP",			PF_REAL_TYPE,	&Time_Stop},
     {"EXPOSURE",		PF_REAL_TYPE,	&Exposure},

     {NULL, 0, NULL}
};

static Param_File_Type *read_obspar_file (void)
{
   char *file;
   Param_File_Type *pf;

   file = make_marx_filename ("obs.par");
   pf = pf_open_parameter_file (file, "r");
   if (pf == NULL)
     {
	fprintf (stderr, "Unable to open obs.par file %s\n", file);
	return NULL;
     }

   if (-1 == pf_get_parameters (pf, ObsPar_Parm_Table))
     return NULL;

   Exposure *= DT_Cor;

   return pf;
}

static int get_simulation_info (void)
{
   Param_File_Type *pf;
   char *file;

   file = make_marx_filename ("marx.par");
   pf = pf_open_parameter_file (file, "r");
   if (pf == NULL)
     {
	marx_error ("*** Unable to open %s\n", file);
	return -1;
     }

   if (-1 == get_marx_pfile_info (pf))
     {
	pf_close_parameter_file (pf);
	return -1;
     }

   pf_close_parameter_file (pf);

   if (NULL == (Obs_Par_Parms = read_obspar_file ()))
     return -1;

   if (-1 == setup_dither ())
     return -1;

   return 0;
}

static int get_date (char *str)
{
   time_t tloc;
   struct tm *tms;

   time (&tloc);
   tms = gmtime (&tloc);

   sprintf (str, "%4d-%02d-%02dT%02d:%02d:%02d",
	    1900 + tms->tm_year, tms->tm_mon + 1, tms->tm_mday,
	    tms->tm_hour, tms->tm_min, tms->tm_sec);

   return 0;
}


/** Copy a single value h from a parameter file p to a fits file f.
 */
static int write_parfile_value (JDFits_Type *f, Param_File_Type *p,
				Fits_Header_Table_Type *h)
{
   char buf[1024];
   char *name;
   int ret, type;
   int i;
   double d;

   name = h->keyword;

   if ((p == NULL)
       || (-1 == (type = pf_get_type (p, name))))
     {
	fprintf (stderr, "**Warning: %s not found in obs.par file.\n", name);
	return 0;
     }

   switch (type)
     {
      default:
	fprintf (stderr, "**Warning: obs.par parameter %s has unsupport type.\n", name);
	return 0;

      case PF_FILE_TYPE:
      case PF_STRING_TYPE:
	if (-1 == pf_get_string (p, name, buf, sizeof (buf)))
	  return -1;
	ret = jdfits_write_header_string (f, name, buf, h->comment);
	break;

      case PF_BOOLEAN_TYPE:
	if (-1 == pf_get_boolean (p, name, &i))
	  return -1;
	ret = jdfits_write_header_logical (f, name, i, h->comment);
	break;

      case PF_INTEGER_TYPE:
	if (-1 == pf_get_integer (p, name, &i))
	  return -1;
	ret = jdfits_write_header_integer (f, name, i, h->comment);
	break;

      case PF_REAL_TYPE:
      case PF_DOUBLE_TYPE:
	if (-1 == pf_get_double (p, name, &d))
	  return -1;
	ret = jdfits_write_header_double (f, name, d, h->comment);
	break;
     }

   return ret;
}

/** Write some keywords in h to fitsfile ft 
 *  The mask integer selects which keywords are written.
 *  It is combined with the location as the keywords as defined in h; only
 *  keywords where location & mask are true will be written.
 */
static int write_extra_headers (JDFits_Type *ft,
				Fits_Header_Table_Type *h,
				unsigned int mask) /*{{{*/
{
   get_date (Todays_Date);

   while (h->keyword != NULL)
     {
	int ret;
	char *str;
	double d;
	int i;

	if ((h->location & mask) == 0)
	  {
	     h++;
	     continue;
	  }

	switch (h->type)
	  {
	   case H_FILE:
	     ret = write_parfile_value (ft, Obs_Par_Parms, h);
	     break;

	   case H_ENV:
	     str = getenv ((char *) h->value);
	     if (str == NULL)
	       {
		  fprintf (stderr, "Warning: Environment variable %s not set.\n", (char *) h->value);
		  str = "Unknown";
	       }
	     ret = jdfits_write_header_string (ft, h->keyword, str, h->comment);
	     break;

	   case H_STR:
	     str = (char *) h->value;
	     if (str == NULL)
	       str = "Unknown";

	     ret = jdfits_write_header_string (ft, h->keyword, str, h->comment);
	     break;

	   case H_PSTR:
	     if (h->value == NULL) str = "Unknown";
	     else str = *(char **)h->value;
	     ret = jdfits_write_header_string (ft, h->keyword, str, h->comment);
	     break;

	   case H_PFLT:
	     if (h->value == NULL) d = 0.0;
	     else d = *(double *) h->value;
	     ret = jdfits_write_header_double (ft, h->keyword, d, h->comment);
	     break;
#ifdef H_INT
	   case H_INT:
	     if (h->value == NULL) i = 0;
	     else i = (int) h->value;

	     ret = jdfits_write_header_integer (ft, h->keyword, i, h->comment);
	     break;
#endif

	   case H_LOG:
	     if (h->value == NULL) i = 0;
	     else i = *(int *) h->value;

	     ret = jdfits_write_header_logical (ft, h->keyword, i, h->comment);
	     break;

	   case H_PINT:
	     if (h->value == NULL) i = 0;
	     else i = *(int *) h->value;

	     ret = jdfits_write_header_integer (ft, h->keyword, i, h->comment);
	     break;

	   case H_COM:
	     ret = jdfits_write_header_comment (ft, h->keyword, h->comment);
	     break;

	   default:
	     fprintf (stderr, "write_extra_headers: %s: type %d not supported.\n",
		      h->keyword, h->type);
	     ret = 0;
	  }

	if (ret == -1)
	  return -1;

	if (h->value == NULL)
	  {
	     if ((h->type != H_COM) && (h->type != H_FILE)
#ifdef H_INT
		 && (h->type != H_INT)
#endif
		)
	       fprintf (stderr, "Warning: Fits keyword %s has no value\n",
			h->keyword);
	  }

	h++;
     }
   return 0;
}

/*}}}*/

#if 0
static int add_goodtime_extension (JDFits_Type *ft)
{
   JDFits_BTable_Keyword_Type keywords[3];
   float64 tmin, tmax;

   memset ((char *) keywords, 0, sizeof (keywords));

   keywords[0].tform = "1D";
   keywords[0].ttype = "START";
   keywords[0].tunit = "s";

   keywords[1].tform = "1D";
   keywords[1].ttype = "STOP";
   keywords[1].tunit = "s";

   if (-1 == jdfits_create_btable_extension (ft, keywords,
					     1, 0, 1,
					     "GTI"))
     return -1;

   HDU_Name_Hdr = "GTI";
   HDU_Class = "OGIP";
   HDU_Class1 = "GTI";
   HDU_Class2 = "ALL";
   Content_Hdr = "GTI";

   if (-1 == write_extra_headers (ft, CC_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, Timing_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, Obs_Info_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, GoodTime_Header_Keywords, 0xFF))
     return -1;

   if (-1 == jdfits_end_header (ft))
     return -1;

   tmin = Time_Start;
   tmax = Time_Stop;

   (void) jdfits_write_float64 (ft, &tmin, 1);
   (void) jdfits_write_float64 (ft, &tmax, 1);

   return jdfits_end_data (ft);
}
#endif

static void compute_dither (double t, double *rap, double *decp, double *rollp,
			    double *dyp, double *dzp, double *dthetap)
{
   double ra, dec, roll;
   double cos_dec;
   JDMVector_Type p;

   t = (2.0*PI) * t;

   ra = Ra_Amp * sin (t/Ra_Period + Ra_Phase);
   dec = Dec_Amp * sin (t/Dec_Period + Dec_Phase);
   roll = Roll_Amp * sin (t/Roll_Period + Roll_Phase);
#if 0
   ra += RA_Sigma * JDMgaussian_random ();
   dec += Dec_Sigma * JDMgaussian_random ();
   roll += Roll_Sigma * JDMgaussian_random ();
#endif

   /* The above are actually offsets.  We need to convert them to absolute values
    * in RA and Dec with the roll properly taken care of.  This is achieved
    * in several steps:
    *
    * 1.  Convert ra/dec offsets to absolute pointing.
    * 2.  Roll resulting vector about nominal pointing
    * 3.  Convert result to ra/dec.
    */

   /* Step 1. */
   cos_dec = cos (dec);
   p = JDMv_ax1_bx2_cx3 (cos (ra) * cos_dec, Nominal_Pointing,
			 cos_dec * sin (ra), RA_Hat,
			 sin (dec), Dec_Hat);

   /* Step 2 */
   roll += Nominal_Roll_In_Radians;
   p = JDMv_rotate_unit_vector (p, Nominal_Pointing, roll);

   /* Step 3.  Note that the JDMv_unit_vector_to_spherical returns values in
    * a traditional spherical coordinate system.  This differs from the ra/dec
    * system in the way dec is defined.  So it needs tweeked. */
   JDMv_unit_vector_to_spherical (p, &dec, &ra);
   dec = PI/2 - dec;

   /* Convert to evil degrees */
   ra *= 180.0/PI;
   dec *= 180.0/PI;
   roll *= 180.0/PI;

   if (ra < 0) ra += 360.0;
   if (roll < 0) roll += 360.0;

   /* Make sure dec is somewhere in range -90, 90 */
   if (dec > 180) dec -= 360;
   else if (dec < -180) dec += 360;
   /* Now dec should be somewhere in range -180 --> 180 */
   if (dec >= 0)
     {
	if (dec > 90) dec = 180 - dec;
     }
   else if (dec < -90) dec = -180 - dec;

   *rap = ra;
   *decp = dec;
   *rollp = roll;

   *dyp = 0;
   *dzp = 0;
   *dthetap = 0;
}

static void compute_quaternion (double ra, double dec, double roll, float64 q[4])
{
   double cos_ra, cos_dec, cos_roll, sin_ra, sin_dec, sin_roll;
   double q0, q1, q2, q3, len;

   /* Note the use of PI/360.  This is because we need cos(ra/2), etc.... */
   ra *= PI/360; dec *= PI/360; roll = (180.0-roll)*PI/360.0;

   cos_ra = cos(ra); cos_dec = cos(dec); cos_roll = cos(roll);
   sin_ra = sin(ra); sin_dec = sin(dec); sin_roll = sin(roll);

   q0 = cos_ra*cos_dec*cos_roll + sin_ra*sin_dec*sin_roll;
   q1 = sin_ra*cos_dec*cos_roll - cos_ra*sin_dec*sin_roll;
   q2 = cos_ra*sin_dec*cos_roll + sin_ra*cos_dec*sin_roll;
   q3 = cos_ra*cos_dec*sin_roll - sin_ra*sin_dec*cos_roll;
   len = sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
   q[0] = q0/len; q[1] = q1/len; q[2] = q2/len; q[3] = q3/len;
}

static int write_marxasp (JDFits_Type *ft)
{
   unsigned int num, i;
   JDFits_BTable_Keyword_Type columns[8+1];   /* last is NULL */

   memset ((char *) columns, 0, sizeof (columns));

   columns[0].ttype = "time";
   columns[0].tform = "1D";
   columns[0].tunit = "s";

   columns[1].ttype = "ra";
   columns[1].tform = "1D";
   columns[1].tunit = "degrees";

   columns[2].ttype = "dec";
   columns[2].tform = "1D";
   columns[2].tunit = "degrees";

   columns[3].ttype = "roll";
   columns[3].tform = "1D";
   columns[3].tunit = "degrees";

   columns[4].ttype = "dy";
   columns[4].tform = "1E";
   columns[4].tunit = "mm";

   columns[5].ttype = "dz";
   columns[5].tform = "1E";
   columns[5].tunit = "mm";

   columns[6].ttype = "dtheta";
   columns[6].tform = "1E";
   columns[6].tunit = "degrees";

   columns[7].ttype = "q_att";
   columns[7].tform = "4D";

   num = (Time_Stop - Time_Start + 1.0) / Delta_Time;

   if (-1 == jdfits_create_btable_extension (ft,
					     columns,
					     num,
					     0, 1,
					     "ASPSOL"))
     return -1;

   HDU_Name_Hdr = "ASPSOL";
   HDU_Class = "ASC";
   HDU_Class1 = "TEMPORALDATA";
   HDU_Class2 = "ASPSOL";
   Content_Hdr = "ASPSOL";

   if (-1 == write_extra_headers (ft, CC_Component, 3))
     return -1;

   if (-1 == write_extra_headers (ft, Timing_Component, 3))
     return -1;

   if (Simulation_Detector_Type & DETECTOR_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Timing_Component, 3))
	  return -1;
     }

   if (-1 == write_extra_headers (ft, Obs_Info_Component, 3))
     return -1;

   if (Simulation_Detector_Type & DETECTOR_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Obs_Info_Component, 3))
	  return -1;
     }

   if (Simulation_Detector_Type & DETECTOR_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Faint_Header_Keywords, 3))
	  return -1;
     }
   else
     {
	if (-1 == write_extra_headers (ft, Hrc_Header_Keywords, 3))
	  return -1;
     }

   if (-1 == jdfits_end_header (ft))
     return -1;

   for (i = 0; i < num; i++)
     {
	double ra, dec, roll, dy, dz, dtheta;
	double t;
	float64 v64[4], q64[4];
	float32 v32[3];
	t = i * Delta_Time;

	compute_dither (t, &ra, &dec, &roll, &dy, &dz, &dtheta);
	v64 [0] = (float64) t + Time_Start;
	v64 [1] = (float64) ra;
	v64 [2] = (float64) dec;
	v64 [3] = (float64) roll;
	v32 [0] = (float32) dy;
	v32 [1] = (float32) dz;
	v32 [2] = (float32) dtheta;

	compute_quaternion (ra, dec, roll, q64);

	if ((-1 == jdfits_write_float64 (ft, v64, 4))
	    || (-1 == jdfits_write_float32 (ft, v32, 3))
	    || (-1 == jdfits_write_float64 (ft, q64, 4)))
	  {
	     fprintf (stderr, "Error writing row to fits file\n");
	     return -1;
	  }
     }

   if (-1 == jdfits_end_data (ft))
     return -1;

   return 0;
}

static int init_fits_file (JDFits_Type *ft)
{
   if (-1 == jdfits_init_null_primary_header (ft))
     return -1;

   HDU_Name_Hdr = "NULL";
   HDU_Class = "ASC";
   HDU_Class1 = "";
   HDU_Class2 = "";
   Content_Hdr = "";

   if (-1 == write_extra_headers (ft, CC_NULL_Component, 1))
     return -1;
   if (-1 == write_extra_headers (ft, Timing_Component, 2))
     return -1;

   if (Simulation_Detector_Type & DETECTOR_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Timing_Component, 2))
	  return -1;
     }

   if (-1 == write_extra_headers (ft, Obs_Info_Component, 2))
     return -1;

   if (Simulation_Detector_Type & DETECTOR_ACIS)
     {
	if (-1 == write_extra_headers (ft, Acis_Obs_Info_Component, 2))
	  return -1;
     }

   if (-1 == jdfits_end_header (ft))
     return -1;

   return 0;
}

static char *Output_File;

static Param_Table_Type Marxasp_Parm_Table [] =
{
     {"TimeDel",	PF_DOUBLE_TYPE,	&Delta_Time},
#if 0
     {"RA_Sigma",	PF_DOUBLE_TYPE,	&RA_Sigma},
     {"Dec_Sigma",	PF_DOUBLE_TYPE,	&Dec_Sigma},
     {"Roll_Sigma",	PF_DOUBLE_TYPE,	&Roll_Sigma},
#endif
     {"OutputFile",	PF_STRING_TYPE,	&Output_File},
     {"MarxDir",	PF_STRING_TYPE,	&Marx_Dir},

     {NULL,		0, NULL}
};

static int marxasp_init (Param_File_Type *p)
{
   if (-1 == pf_get_parameters (p, Marxasp_Parm_Table))
     {
	pf_error ("%s: error getting parameters.", Program_Name);
	return -1;
     }
#if 0
   RA_Sigma *= PI/(180.0*3600);
   Dec_Sigma *= PI/(180.0*3600);
   Roll_Sigma *= PI/(180.0*3600);
#endif
   if (Delta_Time <= 0)
     Delta_Time = 0.256;

   if (2 != marx_file_exists (Marx_Dir))
     {
	fprintf (stderr, "MarxDir directory '%s' does not exist\n", Marx_Dir);
	return -1;
     }

   if (-1 == get_simulation_info ())
     return -1;

   return 0;
}

int main (int argc, char **argv)
{
   JDFits_Type *f;
   Param_File_Type *p;

   sprintf (MarxAsp_Pgm, "%s v%s", Program_Name, MARX_VERSION);

   p = marx_pf_parse_cmd_line ("marxasp.par", NULL, argc, argv);

   if (p == NULL)
     {
	fprintf (stderr, "%s: Error opening parameter file.\n", Program_Name);
	return 1;
     }

   if (-1 == marxasp_init (p))
     {
	pf_close_parameter_file (p);
	return -1;
     }
   pf_close_parameter_file (p);

   if (NULL == (f = jdfits_open_file (Output_File, JDFITS_WRITE_MODE)))
     {
	marx_error ("*** Unable to open output file %s\n", Output_File);
	return 1;
     }

   if (-1 == init_fits_file (f))
     {
	(void) jdfits_close_file (f);
	return 1;
     }

   if (-1 == write_marxasp (f))
     {
	(void) jdfits_close_file (f);
	return 1;
     }
#if 0
   if (-1 == add_goodtime_extension (f))
     {
	(void) jdfits_close_file (f);
	return 1;
     }
#endif
   if (-1 == jdfits_close_file (f))
     {
	return 1;
     }

   return 0;
}

