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

#include "pfile.h"
#include "_pfile.h"

int main (int argc, char **argv)
{
   Param_Type *pf;
   Param_File_Type *p;
   
   if (argc == 1)
     {
	fprintf (stderr, "Usage: %s pfile [p-assign...]\n", argv[0]);
	return -1;
     }
   
   argc--;
   argv++;
   p = pf_parse_cmd_line (argv[0], "rwL", argc, argv);
   if (p == NULL)
     return 1;

   if (argc == 1)
     {
	pf = p->pf;
	while (pf != NULL)
	  {
	     if (-1 == _pf_get_value (p, pf, PF_QUERY_MODE | PF_LEARN_MODE))
	       {
		  return -1;
	       }
	     pf = pf->next;
	  }
     }
   pf_close_parameter_file (p);
	
   return 0;
}
