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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <jdmath.h>

#include <ctype.h>

#include "henke.h"

typedef struct
{
   char name [29];		       /* lowercase */
   char symbol [3];		       /* lowercase */
   float atomic_mass;		       /* gm/mole */
   float density;		       /* gm/cm^3 */
   char *comment;		       /* most likely NULL */
}
Atom_Type;

#define MAX_ATOMIC_NUMBER 103
/* These are arranged in atomic order and may be indexed by atomic number.
 * The first element is a dummy.
 */
static Atom_Type Atoms [MAX_ATOMIC_NUMBER + 1] =
{
     {"",		"",	0,		0,		NULL},
     {"hydrogen",	"H",	1.00794,	0.00008988,	"gas @ STP"},
     {"helium",		"He",	4.00260,	0.0001785,	"gas @ STP"},
     {"lithium",	"Li",	6.941,		0.534,	NULL},
     {"beryllium",	"Be",	9.01218,	1.848,	NULL},
     {"boron",		"B",	10.81,		2.34,	"crystalline"},
     {"carbon",		"C",	12.011,		2.1,	"graphite"},
     {"nitrogen",	"N",	14.0067,	0.0012506,	"gas @ STP"},
     {"oxygen",		"O",	15.9994,	0.001492,	"gas @ STP"},
     {"fluorine",	"F",	18.9984,	0.001696,	"gas @ STP"},
     {"neon",		"Ne",	20.1179,	0.00089990,	"gas @ STP"},
     {"sodium",		"Na",	22.9898,	0.971,	NULL},
     {"magnesium",	"Mg",	24.305,		1.738,	NULL},
     {"aluminum",	"Al",	26.9815,	2.6989,	NULL},
     {"silicon",	"Si",	28.0855,	2.33,	NULL},
     {"phosphorous",	"P",	30.9738,	1.82,	"white form"},
     {"sulphur",	"S",	32.06,		2.07,	"rhombic form"},
     {"chlorine",	"Cl",	35.453,		0.003214,	"gas @ STP"},
     {"argon",		"Ar",	39.948,		0.0017837,	"gas @ STP"},
     {"potassium",	"K",	39.0983,	0.862,	NULL},
     {"calcium",	"Ca",	40.08,		1.55,	NULL},
     {"scandium",	"Sc",	44.9559,	2.989,	NULL},
     {"titanium",	"Ti",	47.88,		4.54,	NULL},
     {"vanadium",	"V",	50.9415,	6.11,	NULL},
     {"chromium",	"Cr",	51.996,		7.19,	NULL},
     {"manganese",	"Mn",	54.9380,	7.32,	NULL},
     {"iron",		"Fe",	55.847,		7.874,	NULL},
     {"cobalt",		"Co",	51.9332,	8.9,	NULL},
     {"nickel",		"Ni",	58.69,		8.902,	NULL},
     {"copper",		"Cu",	63.546,		8.96,	NULL},
     {"zinc",		"Zn",	65.39,		7.133,	NULL},
     {"gallium",	"Ga",	69.72,		6.095,	"liquid"},
     {"germanium",	"Ge",	72.59,		5.323,	NULL},
     {"arsenic",	"As",	74.9216,	5.73,	"grey form"},
     {"selenium",	"Se",	78.96,		4.79,	"grey form"},
     {"bromine",	"Br",	79.904,		3.12,	"liquid"},
     {"krypton",	"Kr",	83.80,		0.003733,	"gas @ STP"},
     {"rubidium",	"Rb",	85.4678,	1.532,	NULL},
     {"strontium",	"Sr",	87.62,		2.54,	NULL},
     {"yttrium",	"Y",	88.9059,	4.469,	NULL},
     {"zirconium",	"Zr",	91.224,		6.506,	NULL},
     {"niobium",	"Nb",	92.9064,	8.57,	NULL},
     {"molybdenum",	"Mo",	95.94,		10.22,	NULL},
     {"technetium",	"Tc",	98.0,		11.50,	NULL},
     {"ruthenium",	"Ru",	101.07,		12.41,	NULL},
     {"rhodium",	"Rh",	102.906,	12.41,	NULL},
     {"palladium",	"Pd",	106.42,		12.02,	NULL},
     {"silver",		"Ag",	107.868,	10.50,	NULL},
     {"cadmium",	"Cd",	112.41,		8.65,	NULL},
     {"indium",		"In",	114.82,		7.31,	NULL},
     {"tin",		"Sn",	118.71,		7.31,	"white form"},
     {"antimony",	"Sb",	121.75,		6.691,	NULL},
     {"tellurium",	"Te",	127.60,		6.24,	NULL},
     {"iodine",		"I",	126.905,	4.93,	NULL},
     {"xenon",		"Xe",	131.29,		0.005887,	"gas @ STP"},
     {"cesium",		"Cs",	132.905,	1.873,	NULL},
     {"barium",		"Ba",	137.33,		3.5,	NULL},
     {"lanthanum",	"La",	138.906,	6.145,	NULL},
     {"cerium",		"Ce",	140.12,		6.657,	NULL},
     {"praseodymium",	"Pr",	140.908,	6.773,	"alpha form"},
     {"neodymium",	"Nd",	144.24,		6.9,	NULL},
     {"promethium",	"Pm",	145.0,		7.22,	NULL},
     {"samarium",	"Sm",	150.36,		7.520,	"alpha form"},
     {"europium",	"Eu",	151.96,		5.243,	NULL},
     {"gadolinium",	"Gd",	157.25,		7.9004,	NULL},
     {"terbium",	"Tb",	158.925,	8.229,	NULL},
     {"dysprosium",	"Dy",	162.50,		8.550,	NULL},
     {"holmium",	"Ho",	164.930,	8.795,	NULL},
     {"erbium",		"Er",	167.26,		9.066,	NULL},
     {"thulium",	"Tm",	168.934,	9.321,	NULL},
     {"ytterbium",	"Yb",	173.04,		6.965,	"alpha form"},
     {"lutetium",	"Lu",	174.967,	9.840,	NULL},
     {"hafnium",	"Hf",	178.49,		13.31,	NULL},
     {"tantalum",	"Ta",	180.9479,	16.654,	NULL},
     {"tungsten",	"W",	183.85,		19.3,	NULL},
     {"rhenium",	"Re",	186.207,	21.02,	NULL},
     {"osmium",		"Os",	190.2,		22.57,	NULL},
     {"iridium",	"Ir",	192.22,		22.42,	NULL},
     {"platinum",	"Pt",	195.08,		21.45,	NULL},
     {"gold",		"Au",	196.967,	19.3,	NULL},
     {"mercury",	"Hg",	200.59,		13.546,	NULL},
     {"thallium",	"Tl",	204.383,	11.85,	NULL},
     {"lead",		"Pb",	207.2,		11.35,	NULL},
     {"bismuth",	"Bi",	208.980,	9.747,	NULL},
     {"polonium",	"Po",	209.0,		9.32,	NULL},
     {"astatine",	"At",	210.0,		0.00001,	"unknown"},
     {"radon",		"Rn",	222,		0.00973,	"gas @ STP"},
     {"francium",	"Fr",	223.0,		0.00001,	"unknown"},
     {"radium",		"Ra",	226.025,	5.0,	NULL},
     {"actinium",	"Ac",	227.028,	10.07,	NULL},
     {"thorium",	"Th",	232.038,	11.72,	NULL},
     {"protactinium",	"Pa",	231.0359,	15.37,	NULL},
     {"uranium",	"U",	238.029,	18.95,	NULL},
     {"neptunium",	"Np",	237.048,	20.25,	NULL},
     {"plutonium",	"Pu",	244.0,		19.84,	NULL},
     {"americium",	"Am",	243.0,		13.67,	NULL},
     {"curium",		"Cm",	247.0,		13.51,	NULL},
     {"berkelium",	"Bk",	247.0,		14.0,	NULL},
     {"californium",	"Cf",	251.0,		0.00001,	"unknown"},
     {"einsteinium",	"Es",	252.0,		0.00001,	"unknown"},
     {"fermium",	"Fm",	257.0,		0.00001,	"unknown"},
     {"mendelevium",	"Md",	258.0,		0.00001,	"unknown"},
     {"nobelium",	"No",	259.0,		0.00001,	"unknown"},
     {"lawrencium",	"Lr",	260.0,		0.00001,	"unknown"}
};

typedef struct Molecule_Table_Type
{
   Henke_Type *h;
   char *formula;
   char *name;
   float density;
   struct Molecule_Table_Type *next;
}
Molecule_Table_Type;

Molecule_Table_Type *Molecule_Table;

static int case_insensitive_strcmp (char *a, char *b)
{
   return strcmp (a, b);
}

static char *do_malloc (unsigned int len, int do_memset)
{
   char *s;

   s = (char *)malloc (len);
   if (s == NULL)
     fprintf (stderr, "Out of memory.\n");
   if (do_memset)
     memset (s, 0, len);
   return s;
}

static char *make_string (char *s)
{
   char *s1;

   s1 = do_malloc (strlen (s) + 1, 0);
   if (s1 != NULL)
     strcpy (s1, s);
   return s1;
}

/* Ok for outside routines to pass NULL. */
static void free_henke_type (Henke_Type *h)
{
   if (h == NULL)
     return;

   if (h->energy != NULL)
     free ((char *) h->energy);
   if (h->f1 != NULL)
     free ((char *) h->f1);
   if (h->f2 != NULL)
     free ((char *) h->f2);
   free ((char *) h);
}

typedef struct Cached_Henke_Table_Type
{
   unsigned int ref_count;
   Henke_Type *table;
   struct Cached_Henke_Table_Type *next;
}
Cached_Henke_Table_Type;

Cached_Henke_Table_Type *Cached_Henke_Tables;

static void free_henke_cache_table (Cached_Henke_Table_Type *h)
{
   if (h != NULL)
     free ((char *) h);
}

void henke_free_henke_table (Henke_Type *h)
{
   Cached_Henke_Table_Type *ch, *last;

   last = NULL;
   ch = Cached_Henke_Tables;

   while (ch != NULL)
     {
	if (ch->table == h)
	  {
	     ch->ref_count -= 1;
	     if (ch->ref_count != 0)
	       return;

	     if (last != NULL)
	       last->next = ch->next;
	     else
	       Cached_Henke_Tables = ch->next;

	     free_henke_cache_table (ch);
	     break;
	  }
	last = ch;
	ch = ch->next;
     }

   free_henke_type (h);
}

static Henke_Type *read_henke_file (char *file)
{
   Henke_Type *h;
   unsigned int max_num_elements, num_elements, i;
   char buf[256];
   FILE *fp;

   if (NULL == (h = (Henke_Type *) do_malloc (sizeof (Henke_Type), 1)))
     return NULL;

   /* Cheat.  Read the file to determine the number of elements. */
   fp = fopen (file, "r");
   if (fp == NULL)
     {
	fprintf (stderr, "Unable to open %s\n", file);
	goto return_error;
     }
#if 1
   fprintf (stdout, "Opening %s\n", file);
#endif

   max_num_elements = 0;
   /* This is a lazy approach */
   while (NULL != fgets (buf, sizeof (buf), fp))
     {
	max_num_elements++;
     }

   if (max_num_elements)
     {
	if ((NULL == (h->energy = (float *) do_malloc (max_num_elements * sizeof (float), 0)))
	    || (NULL == (h->f1 = (float *) do_malloc (max_num_elements * sizeof (float), 0)))
	    || (NULL == (h->f2 = (float *) do_malloc (max_num_elements * sizeof (float), 0))))
	  goto return_error;
     }

   rewind (fp);

   i = 0;
   num_elements = 0;
   while (num_elements < max_num_elements)
     {
	float en;
	float f1, f2;

	if (NULL == fgets (buf, sizeof (buf), fp))
	  break;

	if (3 != sscanf (buf, "%f %f %f", &en, &f1, &f2))
	  continue;

	if (f1 < -9998.0)
	  continue;

	h->energy [num_elements] = en * 1.0e-3;   /* convert to KeV */
	h->f1 [num_elements] = f1;
	h->f2 [num_elements] = f2;

	num_elements++;
     }

   if (num_elements == 0)
     {
	fprintf (stderr, "File %s contains no data.\n", file);
	goto return_error;
     }

   h->num_elements = num_elements;
   fclose (fp);
   return h;

   return_error:
   if (h != NULL) free_henke_type (h);
   if  (fp != NULL) fclose (fp);
   return NULL;
}

#define CLASSICAL_ELECTRON_RADIUS (2.8179380e-9)    /* microns */
#define AVOGADROS_NUMBER (6.022045e23)    /* per mole */
#define HBAR_C (1.9732858e-4)	       /* KeV-microns */

int henke_beta_delta (Henke_Type *h, float density,
		      float **beta_ptr, float **delta_ptr)
{
   unsigned int i;
   float *f1, *f2, *beta, *delta;
   double k;
   float *energy;
   unsigned int num_elements;

   if ((h == NULL)
       || (h->num_elements == 0)
       || (h->f1 == NULL)
       || (h->f2 == NULL)
       || (h->energy == NULL))
     return -1;

   f1 = h->f1;
   f2 = h->f2;
   num_elements = h->num_elements;
   energy = h->energy;

   if (density <= 0.0)
     {
	density = h->density;
	if (density <= 0.0)
	  {
	     fprintf (stderr, "Density of %s not specified.  Assuming 1.0 gm/cm^3.\n",
		      h->name);
	     density = 1.0;
	  }
     }

   if (beta_ptr != NULL)
     {
	if (NULL == (beta = (float *)do_malloc (num_elements * sizeof (float), 0)))
	  return -1;
     }
   else beta = NULL;

   if (delta_ptr != NULL)
     {
	if (NULL == (delta = (float *)do_malloc (num_elements * sizeof (float), 0)))
	  {
	     if (beta != NULL) free ((char *) beta);
	     return -1;
	  }
     }
   else delta = NULL;

   k = 2.0 * PI * CLASSICAL_ELECTRON_RADIUS * (HBAR_C * HBAR_C * AVOGADROS_NUMBER);
   k = k * (density * 1.0e-12) / h->atomic_mass;

   for (i = 0; i < num_elements; i++)
     {
	double en;
	double factor;

	en = energy[i];
	factor = k / (en * en);
	if (beta != NULL) beta [i] = factor * f2[i];
	if (delta != NULL) delta [i] = factor * f1[i];
     }

   if (beta != NULL) *beta_ptr = beta;
   if (delta != NULL) *delta_ptr = delta;

   return 0;
}

static char *Henke_Dir;

int henke_set_data_dir (char *dir)
{
   char *dirbuf;

   if (dir == NULL)
     {
	dir = getenv ("HENKE_DATA_DIR");
	if (dir == NULL)
	  {
	     fprintf (stderr, "HENKE_DATA_DIR not specified. Assuming \".\"\n");
	     dir = ".";
	  }
     }

   if (NULL ==  (dirbuf = make_string (dir)))
     return -1;

   if (Henke_Dir != NULL)
     free (Henke_Dir);

   Henke_Dir = dirbuf;
   return 0;
}

static char *henke_dircat (char *dir, char *name)
{
   static char file [1024];
   unsigned int dirlen;

   if (dir == NULL)
     {
	dir = Henke_Dir;
	if ((dir == NULL)
	    && (-1 == henke_set_data_dir (NULL)))
	  return NULL;
	dir = Henke_Dir;
     }

   if (name != NULL)
     {
	if (name [0] == '/')
	  dir = "";
     }
   else name = "";

   dirlen = strlen (dir);
   if (dirlen + strlen (name) + 2 >= sizeof (file))
     {
	fprintf (stderr, "Filename too long.\n");
	return NULL;
     }

   strcpy (file, dir);
   if (dirlen && (file [dirlen - 1] != '/'))
     {
	file [dirlen] = '/';
	dirlen++;
     }
   strcpy (file + dirlen, name);
   return file;
}

static char *make_henke_data_filename (char *symbol)
{
   char file [1024];

   if (strlen (symbol) + 5 > sizeof (file))
     {
	fprintf (stderr, "Filename too long.\n");
	return NULL;
     }

   sprintf (file, "%s.nff", symbol);
   file[0] = tolower (file[0]);
   return henke_dircat (NULL, file);
}

static Henke_Type *read_henke_atom_table (unsigned int z)
{
   Atom_Type *a;
   char *file;
   Henke_Type *h;

   if (z > MAX_ATOMIC_NUMBER)
     {
	fprintf (stderr, "Element with Z = %u is not supported.\n", z);
	return NULL;
     }

   a = Atoms + z;

   file = make_henke_data_filename (a->symbol);
   if (file == NULL)
     return NULL;

   h = read_henke_file (file);
   if (h != NULL)
     {
	h->density = a->density;
	h->atomic_mass = a->atomic_mass;
	h->formula = a->symbol;
     }
   return h;
}

static Cached_Henke_Table_Type *allocate_henke_cache_table (void)
{
   Cached_Henke_Table_Type *h;

   h = (Cached_Henke_Table_Type *) do_malloc (sizeof (Cached_Henke_Table_Type),
					      1);
   return h;
}

static void free_molcule_table_type (Molecule_Table_Type *m)
{
   if (m == NULL)
     return;
   if (m->name != NULL)
     free (m->name);
   if (m->formula != NULL)
     free (m->formula);
   if (m->h != NULL)
     henke_free_henke_table (m->h);
   free ((char *) m);
}

static int add_molecule (char *name, float density, char *formula)
{
   Molecule_Table_Type *m;

   m = (Molecule_Table_Type *) do_malloc (sizeof (Molecule_Table_Type), 1);
   if (m == NULL)
     return -1;

   if ((NULL == (m->name = make_string (name)))
       || (NULL == (m->formula = make_string (formula))))
     {
	free_molcule_table_type (m);
	return -1;
     }

   m->density = density;
   m->h = NULL;
   m->next = Molecule_Table;
   Molecule_Table = m;
   return 0;
}

static char *skip_non_whitespace (char *s)
{
   char ch;

   while (((ch = *s) != 0)
	  && (ch != ' ')
	  && (ch != '\t')
	  && (ch != '\r')
	  && (ch != '\f')
	  && (ch != '\n'))
     s++;

   return s;
}

static char *skip_whitespace (char *s)
{
   char ch;

   while (((ch = *s) != 0)
	  && ((ch == ' ')
	      || (ch == '\t')
	      || (ch == '\r')
	      || (ch == '\f')
	      || (ch == '\n')))
     s++;

   return s;
}

static void trim_comment (char *str, char *comment_chars)
{
   char ch;
   char *s;

   s = str;

   while ((ch = *s) != 0)
     {
	if (NULL != strchr (comment_chars, ch))
	  break;
	s++;
     }
   *s = 0;

   while (s > str)
     {
	ch = *s;

	if ((ch == 0) || (ch == '\n') || (ch == ' ') || (ch == '\t'))
	  *s-- = 0;
	else break;
     }
}

static char *parse_density (float *density, char *s)
{
   s = skip_whitespace (s);
   if (1 != sscanf (s, "%f", density))
     return NULL;

   return skip_non_whitespace (s);
}

static int read_henke_molecule_definitions (char *file)
{
   char line [512];
   FILE *fp;
   unsigned int line_num;
   char *err;

   file = henke_dircat (NULL, file);
   if (file == NULL)
     return -1;

   fp = fopen (file, "r");
   if (fp == NULL)
     {
	fprintf (stderr, "Unable to open %s\n", file);
	return -1;
     }

   line_num = 0;
   err = NULL;

   while (NULL != fgets (line, sizeof (line), fp))
     {
	char *s;
	char *name;
	float density;
	char *formula;

	line_num++;

	trim_comment (line, "#;!%");
	s = skip_whitespace (line);
	if (*s == 0)
	  continue;

	name = s;
	s = skip_non_whitespace (s);
	if (*s != 0)
	  {
	     *s++ = 0;
	     s = skip_whitespace (s);
	  }

	if (NULL == (s = parse_density (&density, s)))
	  {
	     err = "Error parsing density";
	     goto parse_error;
	  }

	s = skip_whitespace (s);
	if (*s == 0)
	  {
	     err = "Expecting molecular formula";
	     goto parse_error;
	  }

	formula = s;

	if (-1 == add_molecule (name, density, formula))
	  {
	     goto parse_error;
	  }
     }

   fclose (fp);
   return 0;

   parse_error:
   if (err == NULL) err = "Unknown Error";
   fprintf (stderr, "*** %s: %u: %s\n",
	    file, line_num, err);
   return -1;
}

static int float_cmp (float *a, float *b)
{
   if (*a > *b)
     return 1;
   if (*a == *b)
     return 0;
   return -1;
}

static Henke_Type *compute_molecular_henke_type (Henke_Type **atoms, double *weights, unsigned int num_atoms)
{
   float *energies, *ftmp;
   unsigned int num_energies;
   unsigned int i, j;
   float last_energy;
   void (*qsort_fun) (float *, unsigned int, unsigned int, int (*)(float *, float *));
   Henke_Type *h;

   /* We must make sure all edges are taken into account.  One way to do this
    * is to use all energies defined by the individal atom henke tables
    * and remove duplicates.  This is the method I will use.  First count
    * the number of unique energies.  This can be done by sorting.
    */

   num_energies = 0;
   for (i = 0; i < num_atoms; i++)
     num_energies += atoms[i]->num_elements;

   energies = (float *) do_malloc (num_energies * sizeof (float), 0);
   if (energies == NULL)
     {
	return NULL;
     }

   num_energies = 0;
   for (i = 0; i < num_atoms; i++)
     {
	memcpy ((char *) (energies + num_energies),
		(char *) (atoms[i]->energy),
		sizeof (float) * atoms[i]->num_elements);

	num_energies += atoms[i]->num_elements;
     }

   qsort_fun = (void (*)(float *, unsigned int, unsigned int,
			 int (*)(float *, float *))) qsort;
   (*qsort_fun) (energies, num_energies, sizeof (float), float_cmp);

   j = 0;
   last_energy = -1.0;
   for (i = 0; i < num_energies; i++)
     {
	if (last_energy == energies[i])
	  continue;

	last_energy = energies[i];
	energies[j] = last_energy;
	j++;
     }
   num_energies = j;

   if (NULL == (energies = (float *) realloc ((char *) energies,
					      num_energies * sizeof (float))))
     {
	fprintf (stderr, "Not enough memory.\n");
	return NULL;
     }

   if (NULL == (h = (Henke_Type *) do_malloc (sizeof(Henke_Type), 1)))
     {
	free ((char *) energies);
	return NULL;
     }

   h->energy = energies;
   h->num_elements = num_energies;
   if ((NULL == (h->f1 = (float *) do_malloc (sizeof (float) * num_energies, 1)))
       || (NULL == (h->f2 = (float *) do_malloc (sizeof (float) * num_energies, 1)))
       || (NULL == (ftmp = (float *) do_malloc (sizeof (float) * num_energies, 1))))
     {
	free_henke_type (h);
	return NULL;
     }

   for (i = 0; i < num_atoms; i++)
     {
	double z = weights [i];
	Henke_Type *a = atoms[i];

	(void) JDMlog_interpolate_fvector (energies, ftmp, num_energies,
					   a->energy, a->f1, a->num_elements);
	for (j = 0; j < num_energies; j++)
	  h->f1[j] += ftmp [j] * z;

	(void) JDMlog_interpolate_fvector (energies, ftmp, num_energies,
					   a->energy, a->f2, a->num_elements);

	for (j = 0; j < num_energies; j++)
	  h->f2[j] += ftmp [j] * z;

	h->atomic_mass += z * a->atomic_mass;
     }

   free ((char *) ftmp);

   return h;
}

static char *skip_number (char *s)
{
   char ch;
   s = skip_whitespace (s);

   while ((ch = *s++) != 0)
     {
	if (ch == '.')
	  continue;
	if (isdigit (ch))
	  continue;

	break;
     }

   return s - 1;
}

/* This function parses constructs such as
 * CH4 and Si(OH)4SI
 */
static int parse_partial_formula (char *formula, double *z, double factor)
{
   char ch;
   char element [3];
   double new_factor;
   unsigned int i;

   while (1)
     {
	formula = skip_whitespace (formula);

	ch = *formula++;

	if (ch == 0)
	  break;

	if (ch == ')')
	  break;

	if (ch == '(')
	  {
	     char *f = formula;
	     int level = 1;

	     while ((ch = *f) != 0)
	       {
		  f++;
		  if (ch == '(')
		    level++;
		  else if (ch == ')')
		    {
		       level--;
		       if (level == 0)
			 break;
		    }
	       }
	     if (level)
	       {
		  fprintf (stderr, "Molecule parse error: unbalanced parenthesis.\n");
		  return -1;
	       }

	     if (1 == sscanf (f, "%lf", &new_factor))
	       {
		  f = skip_number (f);
	       }
	     else new_factor = 1.0;

	     if (-1 == parse_partial_formula (formula, z, factor * new_factor))
	       return -1;

	     formula = f;
	     continue;
	  }

	if (0 == isupper (ch))
	  {
	     fprintf (stderr, "Molecule parse error: expecting uppercase letter.\n");
	     return -1;
	  }

	element [0] = ch;
	element [1] = 0;

	ch = *formula;

	if (islower (ch))
	  {
	     element[1] = ch;
	     element[2] = 0;
	     formula++;
	  }

	if (1 == sscanf (formula, "%lf", &new_factor))
	  formula = skip_number (formula);
	else
	  new_factor = 1;

	i = 0;
	while (1)
	  {
	     Atom_Type *a;

	     if (i >= MAX_ATOMIC_NUMBER)
	       {
		  fprintf (stderr, "Parse Molecule: %s is not an atom.\n", element);
		  return -1;
	       }

	     a = (Atoms + 1) + i;
	     if (0 == strcmp (a->symbol, element))
	       {
		  z[i] += factor * new_factor;
		  break;
	       }

	     i++;
	  }
     }

   return 0;
}

static Henke_Type *parse_molecular_formula (char *formula)
{
   char *f, ch;
   int paren_level;
   unsigned int i, num_atom_types;
   double z[MAX_ATOMIC_NUMBER];
   Henke_Type *henke_types [MAX_ATOMIC_NUMBER];
   Henke_Type *h;
   double atomic_mass;

   memset ((char *) z, 0, sizeof (z));
   memset ((char *) henke_types, 0, sizeof (henke_types));

   /* Sanity check.  Check for balanced parenthesis */
   paren_level = 0;
   f = formula;
   while ((ch = *f) != 0)
     {
	if (ch == '(') paren_level++;
	else if (ch == ')')
	  {
	     paren_level--;
	     if (paren_level < 0)
	       break;
	  }
	f++;
     }

   if (paren_level)
     {
	fprintf (stderr, "Unbalanced parenthesis: %s\n", formula);
	return NULL;
     }

   if (-1 == parse_partial_formula (formula, z, 1.0))
     {
	return NULL;
     }

   num_atom_types = 0;

   h = NULL;

   atomic_mass = 0.0;
   for (i = 0; i < MAX_ATOMIC_NUMBER; i++)
     {
	if (z[i] == 0.0)
	  continue;

	if (NULL == (henke_types[i] = read_henke_atom_table (i + 1)))
	  goto free_and_return;

	atomic_mass += z[i] * Atoms[i + 1].atomic_mass;
	num_atom_types++;
     }

   if (num_atom_types == 0)
     {
	fprintf (stderr, "Formula %s contain no atoms!", formula);
	goto free_and_return;
     }

   /* Compress the list.  We nolonger need the atomic number */
   num_atom_types = 0;
   for (i = 0; i < MAX_ATOMIC_NUMBER; i++)
     {
	if (henke_types[i] == NULL)
	  continue;

	if (i != num_atom_types)
	  {
	     henke_types[num_atom_types] = henke_types [i]; henke_types [i] = NULL;
	     z[num_atom_types] = z[i]; z[i] = 0.0;
	  }
	num_atom_types++;
     }

   h = compute_molecular_henke_type (henke_types, z, num_atom_types);

   /* drop */

   free_and_return:

   for (i = 0; i < MAX_ATOMIC_NUMBER; i++)
     {
	if (henke_types [i] != NULL)
	  henke_free_henke_table (henke_types [i]);
     }
   return h;
}

static Henke_Type *parse_molecule (Molecule_Table_Type *m)
{
   fprintf (stdout, "Molecule: %s (%s) %f\n", m->name, m->formula, m->density);

   if (m->h != NULL)
     return m->h;

   if (NULL == (m->h = parse_molecular_formula (m->formula)))
     return NULL;

   m->h->density = m->density;
   m->h->formula = m->formula;

   return m->h;
}

static Henke_Type *read_henke_molecule_table (char *name)
{
   Molecule_Table_Type *m;
   Henke_Type *h;

   if (Molecule_Table == NULL)
     {
	if (-1 == read_henke_molecule_definitions ("molecules.def"))
	  return NULL;
     }

   m = Molecule_Table;
   while (m != NULL)
     {
	if (0 == strcmp (m->name, name))
	  {
	     return parse_molecule (m);
	  }

	m = m->next;
     }

   h = parse_molecular_formula (name);

   if (h == NULL)
     fprintf (stderr, "Molecule %s is unknown.\n", name);
   else
     h->formula = "";

   return h;
}

Henke_Type *henke_read_henke_table (char *atom)
{
   Cached_Henke_Table_Type *h;
   unsigned int z;
   Henke_Type *table;

   h = Cached_Henke_Tables;
   while (h != NULL)
     {
	if (0 == strcmp (atom, h->table->name))
	  {
	     h->ref_count += 1;
	     return h->table;
	  }
	h = h->next;
     }

   for (z = 1; z <= MAX_ATOMIC_NUMBER; z++)
     {
	Atom_Type *a = Atoms + z;
	if ((0 == case_insensitive_strcmp (a->name, atom))
	     || (0 == strcmp (a->symbol, atom)))
	  break;
     }

   if (z > MAX_ATOMIC_NUMBER)
     table = read_henke_molecule_table (atom);
   else
     table = read_henke_atom_table (z);

   if (table == NULL)
     return NULL;

   if (NULL == (h = allocate_henke_cache_table ()))
     {
	free_henke_type (table);
	return NULL;
     }

   strncpy (table->name, atom, sizeof (table->name));
   table->name [sizeof (table->name) - 1] = 0;

   h->ref_count = 1;
   h->table = table;
   h->next = Cached_Henke_Tables;
   Cached_Henke_Tables = h;

   return table;
}

#if 0
int main (int argc, char **argv)
{
   char *name;
   Henke_Type *h;
   unsigned int i;
   float *beta, *delta;

   if (argc == 1)
     {
	fprintf (stderr, "Usage: %s element\n", argv [0]);
	return 1;
     }

   name = argv[1];

#if 0
   fprintf (stderr, "#\
Name: %s\n#\
Symbol: %s\n#\
Atomic Number: %u\n#\
Atomic Mass: %g\n#\
Density: %g\n#\
Comment: %s\n\
",
	    a->name,
	    a->symbol,
	    z,
	    a->atomic_mass,
	    a->density,
	    ((a->comment == NULL) ? "" : a->comment));
#endif

   h = henke_read_henke_table (name);
   if (h == NULL)
     return 1;

   if (-1 == henke_beta_delta (h, 0.0, NULL, 0, &beta, &delta))
     {
	free_henke_type (h);
	return 1;
     }

   fprintf (stdout, "#Object: %s [%s]\n#Density: %e gm/cm^3\n",
	    h->name, h->formula, h->density);

   for (i = 0; i < h->num_elements; i++)
     {
#if 1
	fprintf (stdout, "\t%e\t%e\t%e\n",
		 h->energy [i], beta[i], delta [i]);
#else
	fprintf (stdout, "\t%e\t%e\t%e\n",
		 h->energy [i], h->f1[i], h->f2[i]);
#endif
     }
   free ((char *) beta);
   free ((char *) delta);

   free_henke_type (h);
   return 0;
}
#endif
