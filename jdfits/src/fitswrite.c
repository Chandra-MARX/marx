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
/* Basic routines for writing fits files. */
#include "config.h"

#include <stdio.h>
#include <string.h>


#include <memory.h>
#include <ctype.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include "jdfits.h"

static int check_ascii (char *name, int allow_newline)
{
   unsigned int ch;
   char *save_name;
   
   if (name == NULL) return 0;	       /* NULL values mapped to "" by calling routine */
   
   save_name = name;
   
   while ((ch = (unsigned int) *name++) != 0)
     {
	if ((ch >= 127) || (ch < 32))
	  {
	     if (allow_newline && (ch == '\n'))
	       continue;

	     jdfits_error ("%s\n\
  A Fits card must contain only printable ascii characters.",
			 save_name);
	     return -1;
	  }
     }
   return 0;
}


static int check_keyword (char *name)
{
   char *n, ch;
   
   if (name == NULL) return 0;
   if (check_ascii (name, 0) == -1) return -1;
   n = name;
   while ((ch = *n++) != 0)
     {
	if (((ch > 'Z') || (ch < 'A'))
	    && ((ch > '9') || (ch < '0'))
	    && (ch != '-') && (ch != '_'))
	  {
	     jdfits_error ("%s\n\
  A FITS keyword must contain either UPPERCASE characters or '-' or '_'",
			 name);
	     return -1;
	  }
     }
   
   if ((int)strlen (name) > 8)
     {
	jdfits_error ("Keyword (%s) must be less than 8 characters.\n", name);
	return -1;
     }
   return 0;
}

static unsigned int check_comment (char *comment, unsigned int len)
{
   unsigned int clen = strlen (comment);
   
   if (clen && (comment [clen - 1] == '\n'))
     clen--;

   if (clen > len) 
     {
	fprintf (stderr, "Comment too long-- truncated.\n");
	return len;
     }
   return clen;
}

static void jdfits_append_comment (char *buf, char *comment)
{
   unsigned int comment_len;
   unsigned int buf_len = strlen (buf);
   char *b, *bmax;
   
   b = buf + buf_len;
   
   if (comment != NULL)
     {
	int slash_ok;
	unsigned int space;
	
	slash_ok = 0;
	space = JDFITS_CARD_SIZE - buf_len;
	
	if (space > 2)
	  {
	     slash_ok = 1;
	     space -= 2;
	  }
	if (*comment == '/') comment++;
	
	comment_len = check_comment (comment, space);

	if (comment_len)
	  {
	     if (comment_len + 30 < JDFITS_CARD_SIZE)
	       {
		  bmax = buf + 30;
		  while (b < bmax) *b++ = ' ';
	       }
	     if (slash_ok)
	       {
		  *b++ = ' ';
		  *b++ = '/';
	       }
	     strncpy (b, comment, comment_len);
	  }
	b += comment_len;
     }
   bmax = buf + JDFITS_CARD_SIZE;
   while (b < bmax) *b++ = ' ';
}

int jdfits_write_header_string (JDFits_Type *ft, char *name, char *s, char *comment)
{
   char buf[JDFITS_CARD_SIZE + 1];
   char *b, *bmax, *bsave, ch;
   
   if ((-1 == check_keyword (name))
       || (-1 == check_ascii (s, 0))
       || (-1 == check_ascii (comment, 1)))
     return -1;
   
   /* Strings must begin with a ' character in column 11 and have a closing one
    * before the end of the card.  Actually the closing one must be at column
    * 20 or beyond.
    */
   sprintf (buf, "%-8s= '", name);
   bsave = b = buf + strlen(buf);
   bmax = buf + JDFITS_CARD_SIZE;
   
   while (b < bmax)
     {
	ch = *s++;
	if (ch == 0)
	  {
	     /* Now pad it with blanks */
	     bmax = bsave + 8;
	     while (b < bmax) *b++ = ' ';
	     *b = '\'';
	     *(b + 1) = 0;
	     jdfits_append_comment (buf, comment);
	     return jdfits_write (ft, (unsigned char *) buf, JDFITS_CARD_SIZE);
	  }
	
	if (ch == '\'')
	  {
	     *b++ = ch;
	     if (b == bmax) break;
	  }
	*b++ = ch;
     }
   jdfits_error ("String is too long for keyword %s", name);
   return -1;
}

int jdfits_write_header_logical (JDFits_Type *ft, char *name, int val, char *comment)
{
   char buf[JDFITS_CARD_SIZE + 1];
   
   if ((-1 == check_keyword (name))
       || (-1 == check_ascii(comment, 1)))
     return -1;

   if (val) val = 'T'; else val = 'F';
   
   sprintf (buf, "%-8s=%21c", name, val);
   jdfits_append_comment (buf, comment);
   
   return jdfits_write (ft, (unsigned char *) buf, JDFITS_CARD_SIZE);
}

int jdfits_write_header_integer (JDFits_Type *ft, char *name, int val, char *comment)
{
   char buf[JDFITS_CARD_SIZE + 1];
   
   if ((-1 == check_keyword (name))
       || (-1 == check_ascii(comment, 1)))
     return -1;
   
   sprintf (buf, "%-8s= %20d", name, val);
   jdfits_append_comment (buf, comment);
   
   return jdfits_write (ft, (unsigned char *) buf, JDFITS_CARD_SIZE);
}

static int format_double (double d, int prec, char *buf)
{
   char *expon, *decim;

   if (sprintf (buf, "%.*G", prec, d) <= 0)
     return -1;

   /* Unfortunately, I have no idea what form has been used.  So, let's probe. */
   expon = strchr (buf, 'E');
   decim = strchr (buf, '.');
	
   if (expon == NULL)
     {
	if (decim == NULL)
	  strcat (buf, ".0");
     }
   else if (decim == NULL)
     {
	/* Looks like 1E-9, force a decimal point */
	if (sprintf (buf, "%.1E", d) <= 0.0)
	  return -1;
     }
   return 0;
}

int jdfits_write_header_double (JDFits_Type *ft, char *name, double val, char *comment)
{
   char buf[JDFITS_CARD_SIZE + 1];
   char numbuf[64];
   int prec;

   if ((-1 == check_keyword (name))
       || (-1 == check_ascii(comment, 1)))
     return -1;
   
   /* According to the FITS NOST document, the resulting value should be right
    * justified in columns 11-30.  It must contain a decimal point and E must
    * be used if with in exponential form.
    */
   prec = 16;
   do
     {
	if (-1 == format_double (val, prec, numbuf))
	  {
	     strcpy (numbuf, "Floating point exception");
	     jdfits_warning ("jdfits_write_header_double: unable to convert value\n");
	     prec = 0;
	  }
	prec--;
	sprintf (buf, "%-8s= %20s", name, numbuf);
     }
   while ((prec > 5) && (strlen (buf) > 30));

   jdfits_append_comment (buf, comment);
   
   return jdfits_write (ft, (unsigned char *) buf, JDFITS_CARD_SIZE);
}


int jdfits_end_header (JDFits_Type *ft)
{
   char buf [JDFITS_CARD_SIZE + 1];
   strcpy (buf, "END");
   jdfits_append_comment (buf, NULL);
   if (-1 == jdfits_write (ft, (unsigned char *) buf, JDFITS_CARD_SIZE)) return -1;
   return jdfits_flush_output (ft, ' ');
}

int jdfits_write_header_comment (JDFits_Type *ft, char *name, char *comment)
{
   char buf[JDFITS_CARD_SIZE + 1];
   char *bmin, *bmax, *b, ch;

   if (name == NULL) name = "";
   
   /* Check for all blanks which is ok for comments */
   b = name;
   while (((ch = *b++) != 0) && (ch == ' '))
     ;
   
   if (ch == 0) 
     {
	if (strlen (name) > 8)
	  {
	     jdfits_error ("A keyword must be less than 9 characters long.");
	     return -1;
	  }
     }
   else if (-1 == check_keyword (name))
     return -1;
   
   if (comment == NULL) comment = "";
   
   if (-1 == check_ascii (comment, 1)) return -1;
   
   sprintf (buf, "%-8s ", name);
   bmin = buf + 9;
   bmax = buf + JDFITS_CARD_SIZE;
   
   
   while (1)
     {
	b = bmin;
	while (((ch = *comment) != 0) && (b < bmax))
	  {
	     comment++;
	     if (ch == '\n') break;
	     *b++ = ch;
	  }
	
	if ((b == bmax) && (ch != 0))
	  {
	     char *bmin_2;
	     unsigned int diff;
	     
	     diff = bmax - bmin;
	     bmin_2 = bmin + (diff / 2);
	     /* Try to break it at a space. */
	     b--;
	     while ((b > bmin_2) && (*b != ' ')) b--;
	     if (b > bmin_2)
	       {
		  comment -= (bmax - b);
		  comment++; /* +1 to avoid the space */
		  ch = ' ';
	       }
	     else
	       {
		  comment--;
		  b = bmax;
		  *(b - 1) = '\\';
	       }
	  }
	
	while (b < bmax) *b++ = ' ';
	
	*b = 0;
	
	if (-1 == jdfits_write (ft, (unsigned char *) buf, JDFITS_CARD_SIZE))
	  return -1;
	
	if ((ch == '\n') && (*comment == 0)) break;
	if (ch == 0) break;
     }
   return 0;
}

int jdfits_end_data (JDFits_Type *ft)
{   
   return jdfits_flush_output (ft, 0);
}



