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

#include <ctype.h>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include "pfile.h"
#include "_pfile.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

char *_pf_skip_whitespace (char *s)
{
   unsigned char ch;
   
   while (((ch = (unsigned char) *s) != 0)
	  && (ch <= ' '))
     s++;
   return s;
}

static char *malloc_internal (unsigned int len)
{
   char *p;
   
   p = (char *) SLMALLOC (len);
   if (p == NULL)
     {
	PF_Errno = PF_MALLOC_ERROR;
	pf_error ("Memory allocation failure.");
     }
   return p;
}

   
char *_pf_malloc (unsigned int len)
{
   char *p;
   
   p = (char *) malloc_internal (len);
   if (p != NULL) memset (p, 0, len);
   return p;
}

   

char *_pf_create_string (char *s)
{
   char *p;
   unsigned int len;
   
   if (s == NULL) return NULL;
   
   len = strlen (s);
   
   if (NULL != (p = malloc_internal (len + 1))) strcpy (p, s);
   return p;
}

char *_pf_create_nstring (char *s, unsigned int len)
{
   char *p;
   
   if (s == NULL) return NULL;
   
   if (NULL != (p = malloc_internal (len + 1)))
     {
	strncpy (p, s, len);
	p[len] = 0;
     }
   return p;
}

/* Note: Some routines will pass NULL as str. */
char *_pf_strchr (char *str, int ch)
{
   int ch1;
   if (str == NULL) return NULL;
   while (1)
     {
	ch1 = *str;
	if (ch1 == ch) return str;
	if (ch1 == 0) return NULL;
	str++;
     }
}

char *_pf_rstrchr (char *str, int ch)
{
   char *s;
     
   s = str + strlen (str);
   while (1)
     {
	if (*s == ch) return s;
	if (str == s) return NULL;
	s--;
     }
}

int _pf_parse_single_number (char *s, int *ival, double *dval)
{	
   char *save_s;
   
   if (s == NULL) 
     return -1;
   
   s = _pf_skip_whitespace (s);
   save_s = s;
   
   if ((*s == '-') || (*s == '+')) s++;
   
   if ((dval != NULL) && (*s == '.')) s++;
   
   if (0 == (isdigit ((unsigned char)*s)))
     {
	pf_error ("Error parsing field as %s.", ((dval != NULL) ? "real" : "int"));
	PF_Errno = PF_NUMBER_FORMAT_BAD;
	return -1;
     }
   
   while (1)
     {
	if (0 == isdigit ((unsigned char)*s)) break;
	s++;
     }
   
   if (NULL == dval)
     {
	s = _pf_skip_whitespace (s);
	
	if (*s != 0)
	  {
	     pf_error ("Error parsing integer field.");
	     PF_Errno = PF_NUMBER_FORMAT_BAD;
	     return -1;
	  }
	
	*ival = atoi (save_s);
	return 0;
     }
   
   if (*s == '.')
     {
	s++;
	while (1)
	  {
	     if (0 == isdigit ((unsigned char)*s)) break;
	     s++;
	  }
     }
   
   if ((*s == 'e') || (*s == 'E'))
     {
	s++;
	if ((*s == '-') || (*s == '+')) s++;
	while (1)
	  {
	     if (0 == isdigit ((unsigned char)*s)) break;
	     s++;
	  }
     }
   
   s = _pf_skip_whitespace (s);
   if (*s != 0)
     {
	pf_error ("Error parsing field as real.");
	PF_Errno = PF_NUMBER_FORMAT_BAD;
	return -1;
     }
   
   *dval = atof (save_s);
   return 0;
}


int _pf_parse_boolean (char *str, int *ival)
{
   unsigned char *b, ch;
   char buf[4];
   
   *buf = 0;
   
   if (strlen (str) < 4)
     {
	strcpy (buf, str);
     }
   
   b = (unsigned char *) buf;
   while ((ch = *b) != 0) *b++ = (ch | 0x20);   /* lowercase it */
   
   if (!strcmp (buf, "yes"))
     {
	*ival = 1;
     }
   else if (!strcmp (buf, "no"))
     {
	*ival = 0;
     }
   else
     {
	pf_error ("Boolean value requires 'yes' or 'no'.");
	PF_Errno = PF_NUMBER_FORMAT_BAD;
	return -1;
     }
   return 0;
}


char *_pf_unescape_string (char *str)
{
   char *new_string, *s;
   char ch;
   
   if (NULL == (new_string = _pf_create_string (str)))
     return NULL;
   
   s = str = new_string;
   while ((ch = *str++) != 0)
     {
	if (ch == '\\')
	  {
	     ch = *str++;
	     switch (ch)
	       {
		case 0:
		  pf_error ("String has extra \\ character.");
		  SLFREE (new_string);
		  return NULL;
		  
		case 't': ch = '\t'; break;
		case 'n': ch = '\n'; break;
		case 'r': ch = '\r'; break;
		case 'e': ch = 27; break;
		default:
		  break;
	       }
	  }
	*s++ = ch;
     }
   *s = 0;
   return new_string;
}

char *_pf_escape_string (char *str)
{
   char *new_string, *s;
   char ch;
   unsigned int len;
   
   /* Find out how long the result should be. */
   len = strlen (str);
   s = str;
   while ((ch = *s) != 0)
     {
	switch (ch)
	  {
	   case '\t':
	   case '\n':
	   case '\r':
	   case 27:		       /* escape */
	   case '\\':
	   case '\'':
	   case '"':
	     len++;
	  }
	s++;
     }
   
   if (NULL == (new_string = malloc_internal (len + 1)))
     return NULL;
   
   s = new_string;

   while ((ch = *str++) != 0)
     {
	switch (ch)
	  {
	   case '\t':
	     *s++ = '\\'; ch = 't'; break;
	   case '\n':
	     *s++ = '\\'; ch = 'n'; break;
	   case '\r':
	     *s++ = '\\'; ch = 'r'; break;
	   case 27:
	     *s++ = '\\'; ch = 'e'; break;
	   case '\\':
	   case '"':
	   case '\'':
	     *s++ = '\\'; break;
	  }
	*s++ = ch;
     }
   *s = 0;
   return new_string;
}

   
char *_pf_strcat (char *a, char *b)
{
   char *c;
   unsigned int lena = strlen (a);
   
   c = malloc_internal (lena + strlen (b) + 1);
   if (c != NULL)
     {
	strcpy (c, a);
	strcpy (c + lena, b);
     }
   return c;
}

/* Search for characters from list in string str.  If found, return a pointer
 * to the first occurrence.  If not found, return a pointer to the end of the 
 * string (I find this more useful than returning NULL).
 */
char *_pf_strbrk (char *str, char *list)
{
   char ch, ch1, *p;
   
   while ((ch = *str) != 0)
     {
	p = list;
	while ((ch1 = *p) != 0)
	  {
	     if (ch == ch1) return str;
	     p++;
	  }
	str++;
     }
   return str;
}

char *_pf_extract_string_element (char *list, char *delim, unsigned int *pos)
{
   char *last_list;
   unsigned int i;
   
   if ((list == NULL) || (pos == NULL))
     {
	PF_Errno = PF_BAD_ARGUMENT;
	return NULL;
     }

   if (*list == 0)
     return NULL;

   i = *pos;
   while (i > 0)
     {
	/* Note the behavior of _pf_strbrk is exploited here. */
	list = _pf_strbrk (list, delim);
	if (*list == 0) return NULL;
	list++;
	i--;
     }
   
   last_list = _pf_strbrk (list, delim);
   
   *pos += 1;
   
   return _pf_create_nstring (list, (unsigned int) (last_list - list));
}

int _pf_strcasecmp (char *a, char *b)
{
   while (1)
     {
	char cha, chb;
	
	cha = *a;
	chb = *b;
	if (toupper(cha) != toupper(chb))
	  return (int)cha - (int)chb;
	
	if (cha == 0)
	  return 0;
	
	a++;
	b++;
     }
}
