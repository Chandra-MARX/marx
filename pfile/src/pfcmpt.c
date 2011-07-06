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
/* Compatability routines for interacting with programs using the SAO
 * interface.
 */

#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "pfile.h"
#include "_pfile.h"
#include "parameter.h"

paramfile paramopen (char *file, char **argv, int argc, char *mode)
{
   return (paramfile) pf_parse_cmd_line (file, mode, argc, argv);
}

void paramclose (paramfile pfile)
{
   if (-1 == pf_close_parameter_file ((Param_File_Type *) pfile))
     {
     }
}

short pgets (paramfile pfile, char *pname)
{
   int i;
   if (-1 == pf_get_integer ((Param_File_Type *) pfile, pname, &i))
     {
	fprintf (stderr, "pgets: Error getting value for %s\n", pname);
	exit (1);
     }
   return (short) i;
}

int pgeti (paramfile pfile, char *pname)
{
   int i;
   if (-1 == pf_get_integer ((Param_File_Type *) pfile, pname, &i))
     {
	fprintf (stderr, "pgeti: Error getting value for %s\n", pname);
	exit (1);
     }
   return i;
}

int pgetb (paramfile pfile, char *pname)
{
   int i;
   if (-1 == pf_get_boolean ((Param_File_Type *) pfile, pname, &i))
     {
	fprintf (stderr, "pgetb: Error getting value for %s\n", pname);
	exit (1);
     }
   return i;
}

long pgetl (paramfile pfile, char *pname)
{
   int i;
   if (-1 == pf_get_integer ((Param_File_Type *) pfile, pname, &i))
     {
	fprintf (stderr, "pgetl: Error getting value for %s\n", pname);
	exit (1);
     }
   return (long) i;
}

float pgetf (paramfile pfile, char *pname)
{
   double d;
   if (-1 == pf_get_double ((Param_File_Type *) pfile, pname, &d))
     {
	fprintf (stderr, "pgetf: Error getting value for %s\n", pname);
	exit (1);
     }
   return (float) d;
}

double pgetd (paramfile pfile, char *pname)
{
   double d;
   if (-1 == pf_get_double ((Param_File_Type *) pfile, pname, &d))
     {
	fprintf (stderr, "pgetd: Error getting value for %s\n", pname);
	exit (1);
     }
   return d;
}

char *pgetstr (paramfile pfile, char *pname, char *string, int length)
{
   char buf[PF_MAX_LINE_LEN];
   unsigned int len;

   if (-1 == pf_get_string ((Param_File_Type *)pfile, pname, buf, sizeof (buf)))
     {
	fprintf (stderr, "pgetstr: Error getting value for %s\n", pname);
	exit (1);
     }
   len = (unsigned int) length;
   if (len >= PF_MAX_LINE_LEN) len = PF_MAX_LINE_LEN;
   strncpy (string, buf, len);
   if (len) len--;
   string[len] = 0;
   return string;
}

int paccess (paramfile pfile, char *pname)
{
   return pf_parameter_exists (pfile, pname);
}

char *paramgetpath (paramfile pfile)
{
   return pf_get_input_filename (pfile);
}

static void copy_if_non_null (char *to, char *from)
{
   if (to == NULL)
     return;
   if (from == NULL)
     from = "";

   strncpy (to, from, PF_MAX_LINE_LEN);
   if (strlen (from) >= PF_MAX_LINE_LEN)
     to[PF_MAX_LINE_LEN-1] = 0;
}

/* yuk */
int ParamInfo (paramfile pfile, char *name, char *mode, char *type, char *value,
	       char *min, char *max, char *prompt)
{
   Param_Type *pf;
   char buf[PF_MAX_LINE_LEN];
   char *b;

   if (pfile == NULL)
     return 0;
   if (name == NULL)
     return 0;

   if (NULL == (pf = _pf_locate_param_by_type ((Param_File_Type *)pfile, name, 0)))
     return 0;

   /* type */
   b = buf;
   if (pf->type & PF_LIST_TYPE)
     *b++ = '*';
   switch (pf->type & 0xFF)
     {
      case PF_BOOLEAN_TYPE:
	*b++ = 'b';
	break;

      case PF_INTEGER_TYPE:
	*b++ = 'i';
	break;

      case PF_REAL_TYPE:
	*b++ = 'r';
	break;

      case PF_DOUBLE_TYPE:
	*b++ = 'd';
	break;

      case PF_STRING_TYPE:
	*b++ = 's';
	break;

      case PF_FILE_TYPE:
	*b++ = 'f';
	if (pf->type & PF_FILE_EXISTS) *b++ = 'e';
	if (pf->type & PF_FILE_NEXISTS) *b++ = 'n';
	if (pf->type & PF_FILE_READABLE) *b++ = 'r';
	if (pf->type & PF_FILE_WRITABLE) *b++ = 'w';
	break;

      default:
	break;
     }
   *b = 0;
   copy_if_non_null (type, buf);

   /* mode */
   b = buf;
   if (pf->mode & PF_AUTO_MODE) *b++ = 'a';
   if (pf->mode & PF_QUERY_MODE) *b++ = 'q';
   if (pf->mode & PF_HIDDEN_MODE) *b++ = 'h';
   if (pf->mode & PF_LEARN_MODE) *b++ = 'l';
   *b = 0;
   copy_if_non_null (mode, buf);

   /* value */
   if ((pf->flags & PF_CURRENT_VALUE) == 0)
     b = pf->value;
   else
     {
	b = buf;
	switch (pf->type & 0xFF)
	  {
	   case PF_INTEGER_TYPE:
	     sprintf (buf, "%d", pf->current_value.ival);
	     break;

	   case PF_UINT_TYPE:
	     sprintf (buf, "%u", pf->current_value.uval);
	     break;

	   case PF_REAL_TYPE:
	   case PF_DOUBLE_TYPE:
	     sprintf (buf, "%.16g", pf->current_value.dval);
	     break;

	   case PF_FILE_TYPE:
	   case PF_STRING_TYPE:
	     b = pf->current_value.sval;
	     break;

	   case PF_BOOLEAN_TYPE:
	     strcpy (buf, (pf->current_value.ival ? "yes" : "no"));
	     break;

	   default:
	     *b = 0;
	  }
     }
   copy_if_non_null (value, buf);
   copy_if_non_null (min, pf->min);
   copy_if_non_null (max, pf->max);
   copy_if_non_null (prompt, pf->prompt);

   return 1;
}
