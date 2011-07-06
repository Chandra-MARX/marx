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
/* Parse command line arguments */
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include "pfile.h"
#include "_pfile.h"

int pf_get_parameters (Param_File_Type *p, Param_Table_Type *table)
{
   char *name;
   int status;
   char buf[PF_MAX_LINE_LEN];

   while ((name = table->name) != NULL)
     {
	switch (table->type)
	  {
	   case PF_INTEGER_TYPE:
	     status = pf_get_integer (p, name, (int *)table->value);
	     break;

	   case PF_UINT_TYPE:
	     status = pf_get_uint (p, name, (unsigned int *)table->value);
	     break;

	   case PF_REAL_TYPE:
	   case PF_DOUBLE_TYPE:
	     status = pf_get_double (p, name, (double *)table->value);
	     break;

	   case PF_BOOLEAN_TYPE:
	     status = pf_get_boolean (p, name, (int *)table->value);
	     break;

	   case PF_FILE_TYPE:
	   case PF_STRING_TYPE:
	   case PF_STRING0_TYPE:
	     if (table->type == PF_FILE_TYPE)
	       status = pf_get_file (p, name, buf, sizeof (buf));
	     else
	       status = pf_get_string (p, name, buf, sizeof (buf));

	     if (status != -1)
	       {
		  char *s;

		  if ((table->type != PF_STRING_TYPE)
		      && (*buf == 0))
		    {
		       pf_error ("pf_get_parameters: Parameter %s's value must be non-empty", name);
		       return -1;
		    }

		  s = _pf_create_string (buf);

		  if (s == NULL)
		    {
		       return -1;
		    }
		  *((char **)table->value) = s;
	       }
	     break;

	   default:
	     pf_error ("pf_get_parameters: unspported type for %s.", name);
	     return -1;
	  }

	if (status == -1)
	  {
	    pf_error ("Requested parameter %s not found in parameter file.",
		      name);
	     return -1;
	  }

	table++;
     }
   return 0;
}
