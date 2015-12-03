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
#ifndef JD_ARGC_ARGV_H
#define JD_ARGC_ARGV_H

typedef struct
{
   char *name;
   unsigned int type;
#define ARGCARGV_STRING		1      /* -name string_value or -name=str*/
#define ARGCARGV_INTEGER	2      /* -name integer_value */
#define ARGCARGV_BOOLEAN	3      /* -name */
#define ARGCARGV_TOGGLE_BOOLEAN	4
#define ARGCARGV_FLOAT		5
#define ARGCARGV_DOUBLE		6
#define ARGCARGV_FUNCTION	10
#define ARGCARGV_EQS_FLAG	0x100  /* name=value */

   long addr;			       /* address of object to set or function */
   char *help;
}
ArgcArgv_Type;

/* The first argument argc will contain the number of unparsed arguments
 * upon return.  The second argument, argv, will be left pointing at the
 * unparsed argument.  If called from main using argc, argv from main,
 * be sure to bump argc, argv so that argv[0] is not passed.
 *
 * Returns 0 upon success, -1 upon failure.
 */
extern int argcargv (int * /* argc */, char *** /* argv */ , ArgcArgv_Type *);

#endif 				       /* JD_ARGC_ARGV_H */
