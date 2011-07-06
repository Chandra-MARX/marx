/*
    This file is part of the MIT PFILE Parameter Library

    Copyright (C) 2002 Massachusetts Institute of Technology

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
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "pfile.h"
#include "_pfile.h"

static void usage (char *pgm)
{
   fprintf (stderr, "Usage: %s [--verbose] pfile\n", pgm);
   exit (1);
}

int main (int argc, char **argv)
{
   Param_File_Type *p;
   int verbose = 0;
   char *file = NULL, *pgm;

   argc--;
   pgm = *argv++;

   while (argc)
     {
	char *arg;

	argc--;
	arg = *argv++;

	if (!strcmp (arg, "--verbose")) verbose = 1;
	else if (argc == 0)
	  {
	     file = arg;
	  }
	else usage (pgm);
     }

   if (file == NULL) usage (pgm);

   p = pf_open_parameter_file (file, (verbose ? "rWV" : "r"));

   if (p == NULL)
     {
	fprintf (stderr, "Error opening %s as a parameter file.\n", file);
	return 1;
     }

   fprintf (stderr, "Input: %s\nOutput: %s\n",
	    p->input_filename,
	    ((p->output_filename == NULL) ? p->input_filename : p->output_filename));

   pf_close_parameter_file (p);

   return 0;
}
