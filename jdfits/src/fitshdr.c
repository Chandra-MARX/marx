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

#include "jdfits.h"
#include "_jdfits.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLCALLOC calloc
# define SLREALLOC realloc
# define SLFREE free
#endif

JDFits_Keyword_Type *jdfits_find_keyword (JDFits_Type *ft, char *name)
{
   return jdfits_parse_keyword (ft->header, name, JDFITS_ALL_TYPES);
}

JDFits_Read_Keyword_Type *jdfits_open_keywords (JDFits_Type *ft)
{
   JDFits_Read_Keyword_Type *r;
   
   if ((ft == NULL) 
       || (ft->header == NULL))
     {
	jdfits_error ("jdfits_open_keywords: header is NULL");
	return NULL;
     }
   
   if (NULL == (r = (JDFits_Read_Keyword_Type *) jdfits_malloc (sizeof (JDFits_Read_Keyword_Type))))
     return NULL;
   memset ((char *) r, 0, sizeof (JDFits_Read_Keyword_Type));
   
   r->kws_start =  r->kw_next = ft->header->keys;
   r->kws_end = r->kws_start + ft->header->num_keywords;
   return r;
}

JDFits_Keyword_Type *jdfits_read_keyword (JDFits_Read_Keyword_Type *r)
{
   JDFits_Keyword_Type *k, *next;
   
   k = r->kw_next;
   if (k == NULL) return NULL;
   next = k + 1;
   if (next >= r->kws_end) next = NULL;
   r->kw_next = next;
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   return k;
}

void jdfits_rewind_keywords (JDFits_Read_Keyword_Type *r)
{
   if (r == NULL) return;
   r->kw_next = r->kws_start;
}

void jdfits_close_keywords (JDFits_Read_Keyword_Type *r)
{
   if (r != NULL) SLFREE (r);
}


int jdfits_write_keyword (JDFits_Type *ft, JDFits_Keyword_Type *k)
{
   char comment_buf[JDFITS_CARD_SIZE + 1];
   char *comment, *name;
   
   if (k == NULL) return -1;
   
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   
   if ((k->comment_len > 0) && (k->comment != NULL))
     {
	unsigned int len;
	len = k->comment_len;
	if (len > JDFITS_CARD_SIZE) len = JDFITS_CARD_SIZE;
	
	strncpy (comment_buf, (char *) k->comment, len);
	comment_buf[len] = 0;
	comment = comment_buf;
     }
   else comment = NULL;
   
   name = k->name;
   
   switch (k->type & JDFITS_ALL_TYPES)
     {
      case JDFITS_FLOAT32_TYPE:
	return jdfits_write_header_double (ft, name, k->v.fval, comment);
	
      case JDFITS_FLOAT64_TYPE:
	return jdfits_write_header_double (ft, name, k->v.dval, comment);
	
      case JDFITS_INT32_TYPE:
	return jdfits_write_header_integer (ft, name, (int) k->v.lval, comment);
	
      case JDFITS_INT_TYPE:
	return jdfits_write_header_integer (ft, name, k->v.ival, comment);
	
      case JDFITS_INT16_TYPE:
	return jdfits_write_header_integer (ft, name, (int) k->v.hval, comment);
	
      case JDFITS_BOOL_TYPE:
	return jdfits_write_header_logical (ft, name, k->v.ival, comment);
	
      case JDFITS_STRING_TYPE:
	return jdfits_write_header_string (ft, name, (char *)k->v.sval, comment);
	
      case JDFITS_COMMENT_TYPE:
	return jdfits_write_header_comment (ft, name, comment);
	  
      default:
	break;
     }

   return -1;
}

   
int jdfits_extract_comment (JDFits_Keyword_Type *k, char **s)
{   
   static char comment[JDFITS_CARD_SIZE + 1];
   int len;
   
   if (k == NULL) return -1;
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   
   len = k->comment_len;
   if ((k->comment == NULL) || (len <= 0) || ((unsigned int) len > JDFITS_CARD_SIZE))
     return -1;
   
   strncpy (comment, (char *)k->comment, (unsigned int) len);
   comment[len] = 0;
   *s = comment;
   return -1;
}


static int check_type (char *name, 
		       JDFits_Keyword_Type *k, unsigned int desired_type)
{
   if (0 == (k->type & desired_type))
     {
	jdfits_error ("%s: type mismatch.", name);
	return -1;
     }
   return 0;
}


int jdfits_extract_string (JDFits_Keyword_Type *k, char **s)
{   
   if (k == NULL) return -1;
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   if (-1 == check_type ("jdfits_extract_string", k, JDFITS_STRING_TYPE))
     return -1;
   *s = (char *) k->v.sval;
   return 0;
}

int jdfits_extract_integer (JDFits_Keyword_Type *k, int *i)
{   
   if (k == NULL) return -1;
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   if (-1 == check_type ("jdfits_extract_integer", k, JDFITS_INT_MASK))
     return -1;

   if (k->type & JDFITS_INT32_TYPE) *i = (int) k->v.lval;
   else if (k->type & JDFITS_INT16_TYPE) *i = (int) k->v.hval;
   else *i = k->v.ival;
   
   return 0;
}

int jdfits_extract_double (JDFits_Keyword_Type *k, double *d)
{   
   if (k == NULL) return -1;
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   if (-1 == check_type ("jdfits_extract_double", k, JDFITS_NUMBER_MASK))
     return -1;

   switch (k->type)
     {
      case JDFITS_FLOAT64_TYPE: *d = k->v.dval; return 0;
      case JDFITS_FLOAT32_TYPE: *d = (double) k->v.fval; return 0;
      case JDFITS_INT32_TYPE: *d = (double) k->v.lval; return 0;
      case JDFITS_INT16_TYPE: *d = (double) k->v.hval; return 0;
      case JDFITS_INT_TYPE: *d = (double) k->v.ival; return 0;
     }
   
   jdfits_error ("jdfits_extract_double: type %d is not supported", k->type);
   return -1;
}

int jdfits_extract_logical (JDFits_Keyword_Type *k, int *i)
{
   if (k == NULL) return -1;
   (void) jdfits_parse_key (k, JDFITS_ALL_TYPES);
   if (-1 == check_type ("jdfits_extract_logical", k, JDFITS_BOOL_TYPE))
     return -1;
   
   *i = k->v.ival;
   return 0;
}

int jdfits_keyword_exists (JDFits_Type *ft, char *key)
{
   return (NULL != _jdfits_find_keyword (ft->header, key));
}

int jdfits_read_keyword_string (JDFits_Type *ft, char *key, char **value)
{   
   JDFits_Keyword_Type *k;
   char *s;

   if (NULL == (k = _jdfits_find_keyword (ft->header, key)))
     return -1;
   
   if (-1 == jdfits_extract_string (k, &s))
     return -1;
   
   if (NULL == (*value = jdfits_make_string (s)))
     return -1;
   
   return 0;
}

int jdfits_read_keyword_dbl (JDFits_Type *ft, char *key, double *d)
{   
   JDFits_Keyword_Type *k;

   if (NULL == (k = _jdfits_find_keyword (ft->header, key)))
     return -1;
   
   return jdfits_extract_double (k, d);
}

int jdfits_read_keyword_int (JDFits_Type *ft, char *key, int *d)
{   
   JDFits_Keyword_Type *k;

   if (NULL == (k = _jdfits_find_keyword (ft->header, key)))
     return -1;
   
   return jdfits_extract_integer (k, d);
}

   
   
