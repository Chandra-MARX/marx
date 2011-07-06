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

#include <string.h>
#include <memory.h>
#include <ctype.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <stdarg.h>
#include <ctype.h>

#include "jdfits.h"
#include "_jdfits.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLCALLOC calloc
# define SLREALLOC realloc
# define SLFREE free
#endif

#ifndef tolower
# define tolower(x) ((x)|0x20)
#endif

unsigned int Fits_Message_Type = JDFITS_ERRORS | JDFITS_WARNINGS;
int Fits_Error_Num;

void jdfits_error (char *fmt, ...)
{
   va_list ap;

   Fits_Error_Num = 1;

   if ((Fits_Message_Type & JDFITS_ERRORS) == 0) return;

   fprintf (stderr, "ERROR: ");
   va_start(ap, fmt);
   (void) vfprintf(stderr, fmt, ap);
   va_end(ap);

   putc('\n', stderr);
}

void jdfits_warning (char *fmt, ...)
{
   va_list ap;

   if ((Fits_Message_Type & JDFITS_WARNINGS) == 0) return;
   fprintf (stderr, "WARNING: ");
   va_start(ap, fmt);
   (void) vfprintf(stderr, fmt, ap);
   va_end(ap);

   putc('\n', stderr);
}

void jdfits_free (char *s)
{
   if (s != NULL) SLFREE (s);
}

char *jdfits_malloc (unsigned int size)
{
   char *buf;
   buf = (char *) SLMALLOC(size);
   if (buf == NULL)
     jdfits_error ("Memory allocation failure.  Requested size: %d bytes.", size);
   return buf;
}

int jdfits_check_mode (JDFits_Type *ft, unsigned int mode)
{
   if (ft->mode != mode)
     {
	jdfits_error ("File not open for requested operation.");
	return -1;
     }

   return 0;
}

int jdfits_strcasecmp (char *a, char *b)
{
   while (1)
     {
	char cha, chb;

	cha = *a;
	chb = *b;
	if (tolower(cha) != tolower(chb))
	  return (int)cha - (int)chb;

	if (cha == 0)
	  return 0;

	a++;
	b++;
     }
}

char *jdfits_make_string (char *s)
{
   char *s1 = jdfits_malloc (strlen (s) + 1);
   if (s1 != NULL)
     strcpy (s1, s);

   return s1;
}

char *_jdfits_skip_whitespace (char *b)
{
   char ch;
   while (1)
     {
	ch = *b;
	if ((ch == ' ') || (ch == '\t') || (ch == '\n')
	    || (ch == '\f'))
	  {
	     b++;
	     continue;
	  }
       return b;
     }
}

