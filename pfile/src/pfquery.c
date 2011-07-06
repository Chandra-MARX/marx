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

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include "pfile.h"
#include "_pfile.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

static char *trim_string (char *s, int ltw)
{
   char *s1;
   unsigned char ch;

   if (ltw) s = _pf_skip_whitespace (s);
   s1 = s + strlen (s);

   while (s1 > s)
     {
	s1--;
	ch = (unsigned char) *s1;
	if ((ch == '\n') || (ch == ' ') || (ch == '\t'))
	  *s1 = 0;
	else
	  break;
     }
   return s;
}

/* This function only performs syntax checking.  It does not perform
 * range testing.
 */
int _pf_query_current_value (Param_File_Type *p, Param_Type *pf)
{
   char prompt_buf[2 * PF_MAX_LINE_LEN];
   char input_buf[PF_MAX_LINE_LEN];
   int len;
   char *prompt;

   (void) p;

   if (0 == isatty(fileno (stdin)))
     {
	pf_error ("stdin is NOT a tty.");
	return -1;
     }

   if (NULL == (prompt = pf->prompt))
     prompt = pf->name;

   sprintf (prompt_buf, "%s ", prompt);
   len = strlen (prompt_buf);

   if ((pf->min != NULL)
       || (pf->max != NULL))
     {
	char *min = pf->min, *max = pf->max;
	if (min != NULL) min = _pf_escape_string (min);
	if (max != NULL) max = _pf_escape_string (max);

	if (NULL != _pf_strchr (min, '|'))
	  sprintf (prompt_buf + len, "(%s)", min);
	else sprintf (prompt_buf + len,
		      "(%s:%s) ",
		      ((min != NULL) ? min : "-"),
		      ((max != NULL) ? max : "-"));

	if (min != NULL) SLFREE (min);
	if (max != NULL) SLFREE (max);

	len = strlen (prompt_buf);
     }

   if (pf->flags & PF_CURRENT_VALUE)
     {
	char *buf = prompt_buf + len;
	char *tmp;

	switch (pf->type & 0xFF)
	  {
	   case PF_UINT_TYPE:
	     sprintf (buf, "[%u]:", pf->current_value.uval);
	     break;

	   case PF_INTEGER_TYPE:
	     sprintf (buf, "[%d]:", pf->current_value.ival);
	     break;

	   case PF_REAL_TYPE:
	   case PF_DOUBLE_TYPE:
	     sprintf (buf, "[%.16g]:", pf->current_value.dval);
	     break;

	   case PF_FILE_TYPE:
	   case PF_STRING_TYPE:
	     tmp = _pf_escape_string (pf->current_value.sval);
	     sprintf (buf, "[%s]:",
		      ((tmp == NULL) ? "" : tmp));
	     SLFREE (tmp);
	     break;

	   case PF_BOOLEAN_TYPE:
	     sprintf (buf, "[%s]:", (pf->current_value.ival ? "yes" : "no"));
	     break;

	   default:
	     pf_error ("query_value: Type %c is not supported",
		       pf->type & 0xFF);
	     PF_Errno = PF_NOT_IMPLEMENTED;
	     return -1;
	  }
     }
   else strcpy (prompt_buf + len, ":");

   while (1)
     {
	int ival;
	double dval;
	char *str;

	fputs (prompt_buf, stdout);
	fflush (stdout);
	if (NULL == fgets (input_buf, sizeof (input_buf), stdin))
	  {
	     PF_Errno = PF_IO_ERROR;
	     return -1;
	  }

	if (*input_buf == '\n')
	  {
	     if ((pf->flags & PF_CURRENT_VALUE) == 0) continue;
	     return 0;
	  }

	pf->flags |= PF_PARAM_DIRTY;

	switch (pf->type & 0xFF)
	  {
	   case PF_INTEGER_TYPE:
	     str = trim_string (input_buf, 1);
	     if (0 == _pf_parse_single_number (str, &ival, NULL))
	       {
		  pf->current_value.ival = ival;
		  pf->flags |= PF_CURRENT_VALUE;
		  return 0;
	       }
	     break;

	   case PF_UINT_TYPE:
	     str = trim_string (input_buf, 1);
	     if (0 == _pf_parse_single_number (str, &ival, NULL))
	       {
		  pf->current_value.uval = (unsigned int) ival;
		  pf->flags |= PF_CURRENT_VALUE;
		  return 0;
	       }
	     break;

	   case PF_REAL_TYPE:
	   case PF_DOUBLE_TYPE:
	     str = trim_string (input_buf, 1);
	     if (0 == _pf_parse_single_number (str, NULL, &dval))
	       {
		  pf->current_value.dval = dval;
		  pf->flags |= PF_CURRENT_VALUE;
		  return 0;
	       }
	     break;

	   case PF_STRING_TYPE:
	   case PF_FILE_TYPE:
	     str = trim_string (input_buf, 0);
	     if (pf->current_value.sval != NULL)
	       SLFREE (pf->current_value.sval);

	     if (NULL == (pf->current_value.sval = _pf_unescape_string (str)))
	       {
		  pf->flags &= ~PF_CURRENT_VALUE;
		  return -1;
	       }

	     pf->flags |= PF_CURRENT_VALUE;
	     return 0;

	   case PF_BOOLEAN_TYPE:
	     str = trim_string (input_buf, 1);
	     if (0 == _pf_parse_boolean (str, &ival))
	       {
		  pf->current_value.ival = ival;
		  pf->flags |= PF_CURRENT_VALUE;
		  return 0;
	       }
	     break;

	   default:
	     pf_error ("type not implemented.");
	     PF_Errno = PF_NOT_IMPLEMENTED;
	     return -1;
	  }
     }
}
