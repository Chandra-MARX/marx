/*
    Copyright (C) 2002 MIT Center For Space Research

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

#include <ctype.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "jdfits.h"

/* The routines here are for converting fortran format strings to C ones. */

int jdfits_ffmt_to_cfmt (char *ffmt, char *cfmt)
{
   char ch, *p, fmt_type;
   char wstr[16], pstr[16], *s;
   char *smax;

   /* skip whitespace */
   while ((*ffmt == ' ') || (*ffmt == '\t')) ffmt++;

   if (*ffmt == 0)
     {
	jdfits_error ("jdfits_ffmt_to_cfmt: Format is not specified.");
	return -1;
     }
    
    p = ffmt;
    
    /* Skip by repeat count */
    while (isdigit (*p)) p++;
    if (ffmt != p) jdfits_warning ("jdfits_ffmt_to_cfmt: Repeat count ignored. (%s)", ffmt);
    
    fmt_type = *p++ | 0x20;

    *wstr = *pstr = 0;
    /* get the width */
    s = wstr;
   smax = wstr + (sizeof(wstr)-1);
    while ((ch = *p), isdigit (ch) && (s < smax))
      {
	  *s++ = ch;
	  p++;
      }
    *s = 0;

    /* Now precision specifier */
    if (ch == '.')
      {
	  p++;
	  s = pstr;
	 smax = pstr + (sizeof(pstr)-1);
	  while ((ch = *p++), isdigit (ch) && (s < smax))
	    {
		*s++ = ch;
	    }
	  *s = 0;
      }

    *cfmt++ = '%';
    switch (fmt_type)
      {
       case 'd':
	  /* drop */
       case 'e':
	  fmt_type = 'E';
	  /* drop */
       case 'g':
       case 'f':
       case 'i':
	  if (fmt_type == 'i') fmt_type = 'd';
	  
	  s = wstr;
	  while ((ch = *s++) != 0) *cfmt++ = ch;
	  s = pstr;
	  if (*s)
	    {
		*cfmt++ = '.';
		while ((ch = *s++) != 0) *cfmt++ = ch;
	    }
	  *cfmt++ = fmt_type;
	  break;
	  
       case 'a':
	  s = wstr;
	  while ((ch = *s++) != 0) *cfmt++ = ch;
	  *cfmt++ = 's';
	  break;
	  
       case 'z':		       /* hex */
       case 'o':		       /* octal */
       default:
	 jdfits_error ("jdfits_ffmt_to_cfmt: Format type '%c' not implemented.", fmt_type);
	 return -1;
      }
    *cfmt = 0;
    return 0;
}
