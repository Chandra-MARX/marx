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
#include <stdlib.h>
#endif

#include "pfile.h"
#include "_pfile.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

static int output (Param_Type *);

int main (int argc, char **argv)
{
   Param_Type *pf;
   Param_File_Type *p;
   
   if (argc == 1)
     {
	fprintf (stderr, "Usage: %s pfiles...\n", argv[0]);
	return -1;
     }
   
   argc--;
   
   while (argc--)
     {
	argv++;
	if (NULL == (p = pf_open_parameter_file (*argv, "r")))
	  {
	     fprintf (stderr, "plist: Unable to open %s.\n", *argv);
	     continue;
	  }
	
	fprintf (stdout, "\nParameters for %s\n\n", p->input_filename);
	
	/* Make two passes: once for non-hidden then do hidden. */
	pf = p->pf;
	while (pf != NULL)
	  {
	     if ((pf->mode & PF_HIDDEN_MODE) == 0)
	       output (pf);
	     pf = pf->next;
	  }

	pf = p->pf;
	while (pf != NULL)
	  {
	     if (pf->mode & PF_HIDDEN_MODE)
	       output (pf);
	     pf = pf->next;
	  }
	
	pf_close_parameter_file (p);
     }
   return 0;
}


static int output (Param_Type *pf)
{
   char buf[PF_MAX_LINE_LEN];
   char valbuf[PF_MAX_LINE_LEN];
   char *indirect_value = NULL;
   
   if (pf->type == PF_COMMENT_TYPE)
     {
	fprintf (stdout, "%s\n", pf->name);
	return 0;
     }
   
   if (pf->flags & PF_INDIRECT_VALUE)
     {
	if (-1 == _pf_get_indirect_value (pf->pfile, pf->value, &indirect_value))
	  {
	     pf_error ("Error encountered getting indirect value %s.", 
		       pf->value);
	     indirect_value = NULL;
	  }
     }

   *buf = 0;
   if (pf->mode & PF_HIDDEN_MODE)
     {
	buf[0] = '(';
	buf[1] = 0;
     }
   strcat (buf, pf->name);

   *valbuf = 0;
   
   if (pf->value != NULL) 
     strcpy (valbuf, pf->value);
   
   if (pf->flags & PF_INDIRECT_VALUE)
     {
	strcat (valbuf, " -> ");
	strcat (valbuf, 
		((indirect_value == NULL) ? "INDEF" : indirect_value));
     }
   
   if (pf->mode & PF_HIDDEN_MODE)
     strcat (valbuf, ")");
	
   fprintf (stdout, "%14s = %-16s %s\n",
	    buf, valbuf, ((pf->prompt == NULL) ? "" : pf->prompt));
   
   if (indirect_value != NULL) SLFREE (indirect_value);
   return 0;
}

	    
   
   
