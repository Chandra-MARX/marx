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
#include <ctype.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "jdfits.h"
#include "_jdfits.h"

struct _JDFits_User_KW_Table_Type
{
   JDFits_User_KW_Type kw;
   struct _JDFits_User_KW_Table_Type *next;
};


static FILE *open_file_for_read (char *file)
{
   FILE *fp;
   
   fp = fopen (file, "r");

   if (fp == NULL)
     jdfits_error ("Unable to open file %s.", file);
   
   return fp;
}

int jdfits_add_comments_from_file (JDFits_Type *ft, char *file,
				   char *comment_name,
				   char *com_prefix, int not_flag)
{
   FILE *fp;
   char buf [4096];
   unsigned int len;
   
   if ((ft == NULL) || (file == NULL))
     return -1;
   
   if (com_prefix == NULL) com_prefix = "";
   len = strlen (com_prefix);

   if (comment_name == NULL) comment_name = "";
   
   if ((len == 0) && not_flag)
     return 0;
   
   if (NULL == (fp = open_file_for_read (file)))
     return -1;

   while (NULL != fgets (buf, sizeof (buf), fp))
     {
	int matches;
	
	matches = ((len == 0)
		   || ((*buf == *com_prefix) 
		       && (0 == strncmp (buf, com_prefix, len))));
	if (not_flag) matches = !matches;
	
	if (matches == 0) continue;
	
	if (-1 == jdfits_write_header_comment (ft, comment_name, buf + len))
	  {
	     fclose (fp);
	     return -1;
	  }
     }
   
   fclose (fp);
   return 0;
}

static char *skip_word (char *b)
{
   char ch;
   while (1)
     {
	ch = *b;
	if ((ch == ' ') || (ch == '\t') || (ch == '\n')
	    || (ch == '\f') || (ch == '=') || (ch == 0))
	  return b;
	b++;
     }
}


void jdfits_free_user_kw_table (JDFits_User_KW_Table_Type *t)
{
   JDFits_User_KW_Table_Type *next;
   
   while (t != NULL)
     {
	if (t->kw.comment != NULL)
	  free (t->kw.comment);

	if ((t->kw.type == JDFITS_STRING_TYPE)
	    && (t->kw.value.str != NULL))
	  free (t->kw.value.str);

	next = t->next;
	
	memset ((char *) t, 0, sizeof (JDFits_User_KW_Table_Type));
	free ((char *) t);

	t = next;
     }
}

static char *dup_str (char *str)
{
   char *s;
   if (str == NULL)
     return NULL;
   
   s = (char *) jdfits_malloc (1 + strlen (str));
   if (s != NULL)
     strcpy (s, str);
   
   return s;
}

static int set_keyword_comment (JDFits_User_KW_Table_Type *ht, char *keyword, char *comment)
{
   if (comment != NULL)
     {
	if (*comment == 0) comment = NULL;
	if ((comment != NULL)
	    && (NULL == (comment = dup_str (comment))))
	  return -1;
     }
   
   keyword = dup_str (keyword);
   if (keyword == NULL)
     {
	if (comment != NULL) free (comment);
	return -1;
     }
     
   ht->kw.keyword = keyword;
   ht->kw.comment = comment;
   ht->kw.type = JDFITS_USER_KW_COMMENT;

   return 0;
}

   
static int set_keyword_string (JDFits_User_KW_Table_Type *ht, char *keyword, 
				char *str, char *comment)
{
   if (-1 == set_keyword_comment (ht, keyword, comment))
     return -1;
   
   if (NULL == (ht->kw.value.str = dup_str (str)))
     return -1;
   
   ht->kw.type = JDFITS_USER_KW_STR;
   return 0;
}

static int set_keyword_integer (JDFits_User_KW_Table_Type *ht, char *keyword, 
				int i, char *comment)
{
   if (-1 == set_keyword_comment (ht, keyword, comment))
     return -1;
   
   ht->kw.value.i = i;
   ht->kw.type = JDFITS_USER_KW_INT;
   return 0;
}

static int set_keyword_float (JDFits_User_KW_Table_Type *ht, char *keyword, 
			      double d, char *comment)
{
   if (-1 == set_keyword_comment (ht, keyword, comment))
     return -1;
   
   ht->kw.value.d = d;
   ht->kw.type = JDFITS_USER_KW_FLOAT;
   return 0;
}

static int set_keyword_logical (JDFits_User_KW_Table_Type *ht, char *keyword, 
				int i, char *comment)
{
   if (-1 == set_keyword_comment (ht, keyword, comment))
     return -1;
   
   ht->kw.value.i = i;
   ht->kw.type = JDFITS_USER_KW_BOOL;
   return 0;
}

JDFits_User_KW_Table_Type *jdfits_read_user_kw_table (char *file,
							   JDFits_User_KW_Table_Type *old_root)
{
   FILE *fp;
   char buf[512];
   unsigned int linenum;
   JDFits_User_KW_Table_Type *ht, *ht_root;
   
   fp = open_file_for_read (file);
   if (fp == NULL)
     return NULL;
   
   ht = ht_root = NULL;
   
   linenum = 0;
   while (NULL != fgets (buf, sizeof (buf), fp))
     {
	char *b;
	char ch;
	char *keyword;
	char *str;
	char *comment;
	int i;
	JDFits_User_KW_Table_Type *new_ht;
	
	linenum++;
	
	b = _jdfits_skip_whitespace (buf);
	ch = *b;
	
	if ((ch == 0) || (ch == '#') || (ch == ';') || (ch == '@')
	    || (ch == '%'))
	  continue;
	
	new_ht = (JDFits_User_KW_Table_Type *) jdfits_malloc (sizeof (JDFits_User_KW_Table_Type));
	if (new_ht == NULL)
	  goto return_error;

	memset ((char *) new_ht, 0, sizeof (JDFits_User_KW_Table_Type));
	
	if (ht_root == NULL) 
	  ht_root = new_ht;
	else 
	  ht->next = new_ht;
	ht = new_ht;
       
	keyword = b;
	b = skip_word (b);
	if ((*b != 0) && (*b != '='))
	  {
	     *b++ = 0;
	     b = _jdfits_skip_whitespace (b);
	  }
	
	if (*b != '=') 
	  {
	     /* Assume comment. */
	     if (-1 == set_keyword_comment (ht, keyword, b))
	       goto return_error;

	     continue;
	  }

	*b = 0;
	b = _jdfits_skip_whitespace (b + 1);

	ch = *b;
	if ((ch == '\'') || (ch == '"'))
	  {
	     b++;
	     str = b;
	     if (*b) b++;
	     
	     comment = strchr (b, ch);
	     if (comment == NULL)
	       {
		  comment = b + strlen (b);
		  if (*(comment - 1) == '\n')
		    {
		       comment--;
		       *comment = 0;
		    }
	       }
	     else *comment++ = 0;
	     
	     comment = _jdfits_skip_whitespace (comment);
	     
	     if (-1 == set_keyword_string (ht, keyword, str, comment))
	       goto return_error;
	     
	     continue;
	  }
	
	
	/* Now we have either boolean, integer, or float */
	if ((ch == 'T') || (ch == 'F'))
	  {
	     comment = _jdfits_skip_whitespace (b + 1);
	  
	     if (-1 == set_keyword_logical (ht, keyword, (ch == 'T'), comment))
	       goto return_error;

	     continue;
	  }
	
	str = b;
	
	if ((ch == '+') || (ch == '-'))
	  b++;
	while (isdigit (*b)) b++;
	if ((*b == '.') || (*b == 'e') || (*b == 'E'))
	  {
	     double d;
	     if (1 != sscanf (str, "%lf", &d))
	       {
		  jdfits_error ("Error parsing value as float.");
		  goto return_error;
	       }
	     comment = _jdfits_skip_whitespace (skip_word (b));
	     
	     if (-1 == set_keyword_float (ht, keyword, d, comment))
	       goto return_error;
	     
	     continue;
	  }
	
	if (1 != sscanf (str, "%d", &i))
	  {
	     jdfits_error ("Error parsing value as integer.");
	     goto return_error;
	  }
	
	comment = _jdfits_skip_whitespace (skip_word (b));
	
	if (-1 == set_keyword_integer (ht, keyword, i, comment))
	  goto return_error;
	
     }
   
   fclose (fp);
   
   if (old_root != NULL)
     {
	ht = old_root;
	while (ht->next != NULL) ht = ht->next;
	ht->next = ht_root;
	ht_root = old_root;
     }
   return ht_root;
   
   return_error:
   
   jdfits_error ("This error occured on line %u of %s.", linenum, file);
   fclose (fp);
   jdfits_free_user_kw_table (ht_root);

   return NULL;
}

JDFits_User_KW_Type *jdfits_find_user_kw (JDFits_User_KW_Table_Type *t, char *kw)
{
   if (kw == NULL) return NULL;
   
   while (t != NULL)
     {
	if (0 == strcmp (kw, t->kw.keyword))
	  return &t->kw;
	
	t = t->next;
     }
   return NULL;
}

int jdfits_write_user_ky (JDFits_Type *ft, JDFits_User_KW_Type *kw)
{
   if (kw == NULL) return -1;
   
   switch (kw->type)
     {
      case JDFITS_USER_KW_INT:
	return jdfits_write_header_integer (ft, kw->keyword, kw->value.i, kw->comment);
      case JDFITS_USER_KW_STR:
	return jdfits_write_header_string (ft, kw->keyword, kw->value.str, kw->comment);
      case JDFITS_USER_KW_BOOL:
	return jdfits_write_header_logical (ft, kw->keyword, kw->value.i, kw->comment);
      case JDFITS_USER_KW_FLOAT:
	return jdfits_write_header_double (ft, kw->keyword, kw->value.d, kw->comment);
      case JDFITS_USER_KW_COMMENT:
	return jdfits_write_header_comment (ft, kw->keyword, kw->comment);
      default:
	jdfits_error ("jdfits_write_user_ky: type %d unknown.", kw->type);
     }
   return -1;
}

	
       
	  



		  
	     
	
	
       
	
