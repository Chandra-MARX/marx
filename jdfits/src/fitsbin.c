/* -*- mode: C; mode: fold; -*- */
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
/* This file implements the binary table inteface to fits as described
 */
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>


#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <ctype.h>

#include "jdfits.h"
#include "_jdfits.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLCALLOC calloc
# define SLREALLOC realloc
# define SLFREE free
#endif

char *_jdfits_allocate_bytes_of_type (int type, unsigned int nelements)
{
   char *s;
   unsigned int len;

   switch (type)
     {
      case JDFITS_INT16_TYPE:
	len = sizeof (int16);
	break;
	
      case JDFITS_INT32_TYPE:
	len = sizeof (int32);
	break;
	
      case JDFITS_FLOAT32_TYPE:
	len = sizeof (float32);
	break;
	
      case JDFITS_FLOAT64_TYPE:
	len = sizeof (float64);
	break;
	
      case JDFITS_STRING_TYPE:
	len = 1;
	nelements++;		       /* allow for null termination */
	break;

      case JDFITS_BOOL_TYPE:
      case JDFITS_BIT_TYPE:
      case JDFITS_BYTE_TYPE:
      default:
	jdfits_error ("_jdfits_allocate_bytes_of_type: Field type '%d' is not implemented.", type);
	return NULL;
     }
	
   if (nelements == 0)
     nelements = 1;

   if (NULL == (s = (char *) SLCALLOC (nelements, len)))
     jdfits_error ("_jdfits_allocate_bytes_of_type: calloc error.");
   
   return s;
}

   
static JDFits_Keyword_Type *bt_parse_index_kw (JDFits_Header_Type *h, char *kwd, /*{{{*/
					       int i, unsigned int type)
{
   char buf[9];
   sprintf (buf, "%s%d", kwd, i);
   return jdfits_parse_keyword (h, buf, type);
}

/*}}}*/

static void free_bintable (JDFits_Header_Type *fh) /*{{{*/
{
   JDFits_Bintable_Type *bt = fh->ext.bintable;
   if (bt == NULL) return;
   if (bt->finfo != NULL) SLFREE (bt->finfo);
   SLFREE (bt);
   fh->ext.bintable = NULL;
}

/*}}}*/

int jdfits_bintable_parse_headers (JDFits_Type *ft) /*{{{*/
{
   JDFits_Header_Type *h = NULL;
   JDFits_Keyword_Type *kwt = NULL;
   JDFits_Bintable_Type *bt = NULL;
   JDFits_Bintable_Field_Type *bft = NULL;
   int tfields, i;
   
   /* Check to make sure this is a bintable */
   
   if ((ft == NULL) || (NULL == (h = ft->header)))
     {
	jdfits_error ("jdfits_bintable_parse_headers: header is NULL.");
	return -1;
     }
   
   if (h->type == JDFITS_EXTENSION_HEADER)
     {
	if ((h->name != NULL) && !strcmp ((char *) h->name, "BINTABLE"))
	  h->type = JDFITS_BINTABLE_HEADER;
     }
   
   if (h->type != JDFITS_BINTABLE_HEADER)
     {
	jdfits_error ("jdfits_bintable_parse_headers: Not a BINTABLE header.");
	return -1;
     }

   /* Check to make sure the header meets the other basic requirements. */
   if ((h->bitpix != 8) 
       || (h->naxis != 2)
       || (h->gcount != 1))
     {
	jdfits_error ("jdfits_bintable_parse_headers: Header does not meet bintable requirements.");
	return -1;
     }

   if (h->ext.bintable != NULL)
     return 0;			       /* already parsed */

   /* Now look for the tfields header. */
   if ((NULL == (kwt = jdfits_parse_keyword (h, "TFIELDS", JDFITS_INT_TYPE)))
       || ((tfields = kwt->v.ival) < 0) || (tfields > 999))
     {
	jdfits_error ("jdfits_bintable_parse_headers: Bad or missing TFIELDS entry.");
	goto free_error;
     }
   
   if (NULL == (bt = (JDFits_Bintable_Type *) SLMALLOC (sizeof (JDFits_Bintable_Type))))
     {
	jdfits_error ("jdfits_bintable_parse_headers: Malloc Error.");
	goto free_error;
     }
   
   memset ((char *) bt, 0, sizeof (JDFits_Bintable_Type));
   
   if (0 != tfields)
     {
	bft = (JDFits_Bintable_Field_Type *) SLCALLOC (tfields, sizeof (JDFits_Bintable_Field_Type));
	if (bft == NULL) 
	  {
	     jdfits_error ("jdfits_bintable_parse_headers: Calloc Error.");
	     goto free_error;
	  }
     }

   bt->finfo = bft;
   bt->tfields = tfields;
   
   if (NULL != (kwt = jdfits_parse_keyword (h, "EXTNAME", JDFITS_STRING_TYPE)))
     bt->extname = (char *) kwt->v.sval;
     
   for (i = 1; i <= tfields; i++)
     {
	char buf[9], *b;
	unsigned char ch;
	unsigned int repeat, type, size;
	char *fmt, *default_format;
	
	JDFits_Bintable_Field_Type *bfti = bft + (i - 1);
	
	
	sprintf (buf, "TFORM%d", i);
	
	kwt = jdfits_parse_keyword (h, buf, JDFITS_STRING_TYPE);
	if (kwt == NULL)
	  {
	     jdfits_error ("jdfits_bintable_parse_headers: Bad or missing %s entry.", buf);
	     goto free_error;
	  }
	
	b = (char *) kwt->v.sval;
	ch = *b;
	
	/* Some files produced by ftools use spaces like:  "  1J  " */
	while (ch == ' ')
	  {
	     b++;
	     ch = *b;
	  }

	if (isdigit (ch) == 0)
	  repeat = 1;
	else
	  {
	     repeat = 0;
	     while (((ch = (unsigned char) *b) != 0) && isdigit (ch))
	       {
		  repeat = 10 * repeat + (ch - '0');
		  b++;
	       }
	  }

	bfti->repeat = repeat;
	
	default_format = "% 16.9E";
	
	type = 0;
	pointer_type_label:
	
	switch (ch)
	  {
	   case 'L':
	     type |= JDFITS_BOOL_TYPE;
	     size = 8;
	     fmt = "%d";
	     break;
	   case 'X':
	     switch (repeat)
	       {
		  /* FIXME: I need to use hex formats, but that requires
		   * proper parsing of the format string in the dump function.
		   */
		case 8:
		  type |= JDFITS_BYTE_TYPE;
		  bfti->repeat = 1;
		  size = 8;
		  fmt = "%d";
		  break;
		case 16:
		  type |= JDFITS_INT16_TYPE;
		  bfti->repeat = 1;
		  size = 16;
		  fmt = JDFITS_FMT_16;
		  break;
		case 32:
		  type |= JDFITS_INT32_TYPE;
		  bfti->repeat = 1;
		  size = 32;
		  fmt = JDFITS_FMT_32;
		  break;

		default:
		  type |= JDFITS_BIT_TYPE; 
		  size = 1;
		  fmt = "%X";
		  break;
	       }
	     break;
	   case 'B':
	     type |= JDFITS_BYTE_TYPE;
	     size = 8;
	     fmt = "%c";
	     break;
	   case 'I':
	     type |= JDFITS_INT16_TYPE;
	     size = 16;
	     fmt = JDFITS_FMT_16;
	     break;
	     
	   case 'P':		       /* variable length array */
	     if (repeat > 1)
	       {
		  jdfits_error ("jdfits_bintable_parse_headers: The repeat count must be 0 or 1 for P type.");
		  goto free_error;
	       }
	     if (repeat) repeat = 2;
	     bfti->repeat = repeat;
	     b++;
	     ch = *b;
	     type = JDFITS_POINTER_FLAG;
	     goto pointer_type_label;
	       
	   case 'J':
	     type |= JDFITS_INT32_TYPE;
	     size = 32;
	     fmt = JDFITS_FMT_32;
	     break;
	   case 'A':
	     type |= JDFITS_STRING_TYPE;
	     size = 8;
	     fmt = "%s";
	     break;
	   case 'E':
	     type |= JDFITS_FLOAT32_TYPE;
	     fmt = "% 15.7E";
	     size = 32;
	     break;
	   case 'D':
	     type |= JDFITS_FLOAT64_TYPE;
	     size = 64;
	     fmt = "% 16.9E";
	     break;

	   case 0:
	     ch = '0';
	     /* drop */
	   case 'C':		       /* complex */
	   case 'M':		       /* double precision complex */
	   default:
	     jdfits_error ("jdfits_bintable_parse_headers: Field type '%c' is not implemented. (%s)", 
			 ch, buf);
	     goto free_error;
	  }
	bfti->size = size;
	bfti->type = type;
	
	
	/* Now pick up the remaining indexed keywords. */
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TTYPE", i, JDFITS_STRING_TYPE)))
	  bfti->ttype = NULL;
	else bfti->ttype = kwt->v.sval;
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TUNIT", i, JDFITS_STRING_TYPE)))
	  bfti->tunit = NULL;
	else bfti->tunit = kwt->v.sval;
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TDIM", i, JDFITS_STRING_TYPE)))
	  bfti->tdim = NULL;
	else bfti->tdim = kwt->v.sval;
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TNULL", i, JDFITS_INT_TYPE)))
	  bfti->tnull = 0;
	else bfti->tnull = kwt->v.ival;
	
	bfti->has_scaling = 0;
	if (NULL == (kwt = bt_parse_index_kw (h, "TSCAL", i, JDFITS_FLOAT64_TYPE)))
	  {
	     bfti->tscal = 1.0;
	  }
	else 
	  {
	     bfti->tscal = kwt->v.dval;
	     bfti->has_scaling = 1;
	  }
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TZERO", i, JDFITS_FLOAT64_TYPE)))
	  bfti->tzero = 0.0;
	else
	  {
	     bfti->tzero = kwt->v.dval;
	     bfti->has_scaling = 1;
	  }
	 
	if (NULL == (kwt = bt_parse_index_kw (h, "TCRPX", i, JDFITS_FLOAT64_TYPE)))
	  bfti->tcrpx = 0.0;
	else 
	  {
	     bfti->tcrpx = kwt->v.dval;
	     bfti->has_scaling = 1;
	  }
	
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TCRVL", i, JDFITS_FLOAT64_TYPE)))
	  bfti->tcrvl = 0.0;
	else 
	  {
	     bfti->tcrvl = kwt->v.dval;
	     bfti->has_scaling = 1;
	  }
	
	if (NULL == (kwt = bt_parse_index_kw (h, "TCDLT", i, JDFITS_FLOAT64_TYPE)))
	  bfti->tcdlt = 1.0;
	else 
	  {
	     bfti->tcdlt = kwt->v.dval;
	     bfti->has_scaling = 1;
	  }
#if 0
	/* We do not yet handled the unsigned int hack.  So avoid float format */
	if (bfti->has_scaling) fmt = default_format;
#endif		
	if ((NULL == (kwt = bt_parse_index_kw (h, "TDISP", i, JDFITS_STRING_TYPE)))
	    || (-1 == jdfits_ffmt_to_cfmt (kwt->v.sval, bfti->tdisp)))
	  strcpy ((char *) bfti->tdisp, fmt);
     }
   
   h->ext.bintable = bt;
   h->free_routine = free_bintable;

   /* Now set up the remaining fields.  Note that the positions of 
    * naxis1 and naxis2 are fixed.
    */
   bt->naxis1 = h->keys[3].v.ival;
   bt->naxis2 = h->keys[4].v.ival;
   
   return 0;
   
   /* Only get here if something goes wrong.  Here all malloced items in this
    * function are freed.
    */
   free_error:
   h->ext.bintable = NULL;
   if (bt != NULL) SLFREE (bt);
   if (bft != NULL) SLFREE (bft);
   return -1;
}

/*}}}*/

static unsigned char *check_buf_size (unsigned char **bufp, unsigned int *max_lenp, unsigned int len)
{
   unsigned char *buf;

   if (len == 0)
     len = 1;

   if (len <= *max_lenp)
     return *bufp;
   
   if (*bufp == NULL)
     buf = (unsigned char *) malloc (len);
   else
     buf = (unsigned char *) realloc ((char *) *bufp, len);
   
   if (buf == NULL)
     {
	jdfits_error ("Out of memory: %u bytes requested", len);
	return NULL;
     }
   *bufp = buf;
   *max_lenp = len;
   return buf;
}

   

int jdfits_bintable_dump_data (JDFits_Type *ft, int scale, FILE *fp,
			       char *colname) /*{{{*/
{
   JDFits_Bintable_Field_Type *btf;
   JDFits_Bintable_Type *bt;
   int rows, cols, i;
   unsigned char *buf;
   unsigned int buf_size;
   
   buf = NULL;
   buf_size = 0;

   if ((ft->header->type != JDFITS_BINTABLE_HEADER)
       || ((bt = ft->header->ext.bintable) == NULL))
     {
	jdfits_error ("jdfits_bintable_dump_data: Not a binary table.");
	return -1;
     }
   
   if (-1 == jdfits_read_open_data (ft)) return -1;
   
   btf = bt->finfo;
   rows = bt->naxis2;
   cols = bt->tfields;
   
   /* dump out the table heading */
   putc ('#', fp);
   
   for (i = 0; i < cols; i++)
     {
	char fmt[20], *f, *f1, ch;
	unsigned int repeat = btf[i].repeat;
	
	/* look at the format to guess the field width */
	if (btf[i].type & JDFITS_POINTER_FLAG)
	  f = JDFITS_FMT_32;
	else
	  f = btf[i].tdisp;

	f1 = fmt;
	
	if (*f == '%') f++;
	/* Take the extra space associated with, e.g., "% 12f" into account. */
	if (*f == ' ')
	  f++;

	*f1++ = '%';
	
	while (((ch = *f++) != 0) && isdigit (ch)) *f1++ = ch;
	*f1 = 's';
	*(f1 + 1) = 0;
	
	if ((f = btf[i].ttype) == NULL) f = "??";

	if ((colname != NULL)
	    && jdfits_strcasecmp (colname, f))
	  continue;

	/* Do not print large arrays */
	if ((colname == NULL) && (repeat > 128) && (cols > 1)
	    && (btf[i].type != JDFITS_STRING_TYPE))
	  repeat = 16;

	if (btf[i].type == JDFITS_STRING_TYPE)
	  {
	     sprintf (fmt, "%%-%ds", repeat);
	     repeat = 1;
	  }

	while (repeat--)
	  {
	     if (fprintf (fp, (char *) fmt, f) < 0)
	       return -1;

	     if (EOF == putc (' ', fp))
	       return -1;
	     f = "";
	  }
     }
   putc ('\n', fp);
   
   /* dump out the units */
   putc ('#', fp);
   
   for (i = 0; i < cols; i++)
     {
	char fmt[20], *f, *f1, ch;
	unsigned int repeat = btf[i].repeat;

	if ((colname != NULL)
	    && ((btf[i].ttype == NULL)
		|| jdfits_strcasecmp (colname, btf[i].ttype)))
	  continue;

	if ((colname == NULL) && (repeat > 128) && (cols > 1)
	    && (btf[i].type != JDFITS_STRING_TYPE))
	  repeat = 16;

	/* look at the format to guess the field width */
	if (btf[i].type & JDFITS_POINTER_FLAG)
	  f = JDFITS_FMT_32;
	else
	  f = btf[i].tdisp;

	f1 = fmt;
	if (*f == '%') f++;

	/* Take the extra space associated with, e.g., "% 12f" into account */
	if (*f == ' ')
	  f++;

	*f1++ = '%';
	
	while (((ch = *f++) != 0) && isdigit (ch)) *f1++ = ch;
	*f1 = 's';
	*(f1 + 1) = 0;
	
	if (btf[i].type == JDFITS_STRING_TYPE)
	  {
	     sprintf (fmt, "%%-%ds", repeat);
	     repeat = 1;
	  }

	if ((f = btf[i].tunit) == NULL) f = "";
	
	while (repeat--)
	  {
	     if (0 > fprintf (fp, (char *) fmt, f))
	       return -1;
	     if (EOF == putc (' ', fp))
	       return -1;
	  }
     }
   putc ('\n', fp);
   
	
   while (rows-- > 0)
     {
	putc (' ', fp);	       /* to match # in column names */
	for (i = 0; i < cols; i++)
	  {
	     char *fmt;
	     double tscal, tzero;
	     unsigned int repeat;
	     int32 *i32, *i32max;
	     int16 *i16, *i16max;
	     float32 *f32, *f32max;
	     float64 *f64, *f64max;
	     unsigned char *b, *bmax;
	     int has_scaling;
	     unsigned int type;
	     int do_output = 1;
	     unsigned int num_to_read;

	     if ((colname != NULL)
		 && ((btf[i].ttype == NULL)
		     || jdfits_strcasecmp (colname, btf[i].ttype)))
	       do_output = 0;

	     repeat = btf[i].repeat;
	     fmt = btf[i].tdisp;
	     num_to_read = repeat;

	     if ((colname == NULL) && (repeat > 128) && (cols > 1)
		 && (btf[i].type != JDFITS_STRING_TYPE))
	       repeat = 16;
	     
#if 0
	     has_scaling = btf[i].has_scaling;
#else				       /* avoid unsiged int hack */
	     has_scaling = 0;
#endif
	     tscal = btf[i].tscal;
	     tzero = btf[i].tzero;
	     if ((has_scaling) && (scale > 0))
	       {
		  tzero = btf[i].tcrvl + btf[i].tcdlt * (tzero - btf[i].tcrpx);
		  tscal = btf[i].tcdlt * tscal;
	       }

	     type = btf[i].type;
	     
	     if (type & JDFITS_POINTER_FLAG)
	       {
		  type = JDFITS_INT32_TYPE;
	       }

	     
	     switch (type)
	       {
		case JDFITS_INT32_TYPE:
		  if (NULL == (i32 = (int32 *) check_buf_size (&buf, &buf_size, num_to_read * 4)))
		    goto return_error;
		  if (num_to_read != jdfits_read_int32 (ft, i32, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;

		  i32max = i32 + repeat;
		  if (has_scaling && (scale >= 0)) while (i32 < i32max)
		    {
		       if ((0 > fprintf (fp, (char *) fmt, tzero + tscal * *i32))
			   || (EOF == putc (' ', fp)))
			 goto return_error;

		       i32++;
		    }
		  else while (i32 < i32max)
		    {
		       if ((0 > fprintf (fp, JDFITS_FMT_32, (long) *i32))
			   || (EOF == putc (' ', fp)))
			 goto return_error;

		       i32++;
		    }
		  break;
		  
		case JDFITS_INT16_TYPE:
		  if (NULL == (i16 = (int16 *) check_buf_size (&buf, &buf_size, num_to_read * 2)))
		    goto return_error;
		  if (num_to_read != jdfits_read_int16 (ft, i16, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;
		  i16max = i16 + repeat;
		  if ((tzero == 32768) && (tscal == 1))
		    {		       /* unsigned hack */
		       fmt = "%hu";
		       has_scaling = 1;
		       scale = 1;
		    }

		  if (has_scaling && (scale >= 0)) while (i16 < i16max)
		    {
		       if ((0 > fprintf (fp, (char *) fmt, (short) (tzero + tscal * *i16)))
			   || (EOF == putc (' ', fp)))
			 goto return_error;

		       i16++;
		    }
		  else while (i16 < i16max)
		    {
		       if ((0 > fprintf (fp, JDFITS_FMT_16, *i16))
			   || (EOF == putc (' ', fp)))
			 goto return_error;

		       i16++;
		    }
		  break;
		  
		case JDFITS_FLOAT64_TYPE:
		  if (NULL == (f64 = (float64 *) check_buf_size (&buf, &buf_size, num_to_read * 8)))
		    goto return_error;
		  if (num_to_read != jdfits_read_float64 (ft, f64, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;
		  f64max = f64 + repeat;
		  while (f64 < f64max)
		    {
		       double the_number = *f64;
		       int eof;

		       if (isnan (the_number))
			 eof = fprintf (fp, "%16s", "NaN");
		       else
			 eof = fprintf (fp, (char *) fmt, tzero + tscal * the_number);

		       if ((eof < 0) || (EOF == putc (' ', fp)))
			 goto return_error;

		       f64++;
		    }
		  break;
		  
		case JDFITS_FLOAT32_TYPE:
		  if (NULL == (f32 = (float32 *) check_buf_size (&buf, &buf_size, num_to_read * 4)))
		    goto return_error;
		  if (num_to_read != jdfits_read_float32 (ft, f32, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;
		  f32max = f32 + repeat;
		  while (f32 < f32max)
		    {
		       float the_number = *f32;
		       int eof;

		       if (isnan (the_number))
			 eof = fprintf (fp, "%15s", "NaN");
		       else
			 eof = fprintf (fp, (char *) fmt, tzero + tscal * the_number);

		       if ((eof < 0) || (EOF == putc (' ', fp)))
			 goto return_error;

		       f32++;
		    }
		  break;
		  
		case JDFITS_BYTE_TYPE:
		  if (NULL == (b = (unsigned char *) check_buf_size (&buf, &buf_size, num_to_read * 1)))
		    goto return_error;
		  if (num_to_read != jdfits_read_bytes (ft, b, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;
		  bmax = b + repeat;
		  while (b < bmax)
		    {
		       if ((0 > fprintf (fp, "%4d", (int) (tzero + tscal * (float) *b)))
			   || (EOF == putc (' ', fp)))
			 goto return_error;

		       b++;
		    }
		  break;
		  
		case JDFITS_STRING_TYPE:
		  if (NULL == (b = (unsigned char *) check_buf_size (&buf, &buf_size, num_to_read * 1)))
		    goto return_error;
		  if (num_to_read != jdfits_read_bytes (ft, b, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;
		  bmax = b + repeat;
		  while (b < bmax)
		    {
		       if (*b == 0) 
			 {
			    do
			      {
				 if (EOF == putc (' ', fp))
				   goto return_error;

				 b++;
			      }
			    while (b < bmax);
			    break;
			 }
		       if (EOF == putc (*b, fp))
			 goto return_error;
		       b++;
		    }
		  if (EOF == putc (' ', fp))
		    goto return_error;
		  break;
		  
		  
		case JDFITS_BOOL_TYPE:
		  if (NULL == (b = (unsigned char *) check_buf_size (&buf, &buf_size, num_to_read * 1)))
		    goto return_error;
		  if (num_to_read != jdfits_read_bytes (ft, b, num_to_read))
		    {
		       jdfits_error ("Read error.");
		       goto return_error;
		    }
		  if (do_output == 0)
		    break;
		  bmax = b + repeat;
		  while ((b < bmax) && (*b == ' '))
		    b++;
		  if (b < bmax)
		    {
		       if (EOF == putc (*b, fp))
			 goto return_error;
		    }
		  else if (EOF == putc (*b, fp))
		    goto return_error;

		  if (EOF == putc (' ', fp))
		    goto return_error;
		  break;

		case JDFITS_BIT_TYPE:
		default:
		  jdfits_error ("read Not implemented (0x%X)", btf[i].type);
		  goto return_error;
	       }
	  }
	putc ('\n', fp);
     }
   jdfits_read_close_data (ft);
   jdfits_free ((char *) buf);
   return 0;
   
   return_error:
   jdfits_read_close_data (ft);
   jdfits_free ((char *) buf);

   return -1;
}

/*}}}*/
     
int jdfits_bintable_read_data (JDFits_Type *ft) /*{{{*/
{
   unsigned int tfields, i, j;
   unsigned int rows;
   JDFits_Bintable_Type *fbt;
   
   fbt = ft->header->ext.bintable;
   
   rows = fbt->naxis2;
   tfields = (unsigned int) fbt->tfields;
		    
   for (i = 0; i < tfields; i++)
     {
	unsigned int nelements;
	unsigned int type;
	JDFits_Bintable_Field_Type *finfo;
	JDFits_Data_Array_Type *fdat;
	
	finfo = &fbt->finfo[i];
	nelements = rows * finfo->repeat;
	type = finfo->type;
	fdat = &(finfo->data_array);
	
	fdat->nelements = nelements;
	fdat->type = type;

	if (NULL == (fdat->v.sval = _jdfits_allocate_bytes_of_type (type, nelements)))
	  goto close_error;
     }
   
   /* Now data arrays are malloced, so read in the data.  There is no scaling
    * performed.
    */

   if (-1 == jdfits_read_open_data (ft)) goto close_error;

   for (i = 0; i < rows; i++)
     {
	for (j = 0; j < tfields; j++)
	  {
	     JDFits_Bintable_Field_Type *finfo;
	     JDFits_Data_Array_Type *fdat;
	     unsigned int repeat; 
	     unsigned int offset; 
	     
	     finfo = &fbt->finfo[j];
	     fdat = &(finfo->data_array);
	     repeat = finfo->repeat;
	     offset = i * repeat;
	     
	     switch (fdat->type)
	       {
		case JDFITS_INT16_TYPE:
		  if (repeat != jdfits_read_int16 (ft, fdat->v.hval + offset, repeat))
		    {
		       jdfits_error ("Read error.");
		       goto close_error;
		    }
		  break;
		  
		case JDFITS_INT32_TYPE:
		  if (repeat != jdfits_read_int32 (ft, fdat->v.lval + offset, repeat))
		    {
		       jdfits_error ("Read error.");
		       goto close_error;
		    }
		  break;
	     
		case JDFITS_FLOAT32_TYPE:
		  if (repeat != jdfits_read_float32 (ft, fdat->v.fval + offset, repeat))
		    {
		       jdfits_error ("Read error.");
		       goto close_error;
		    }
		  break;
	     
		case JDFITS_FLOAT64_TYPE:
		  if (repeat != jdfits_read_float64 (ft, fdat->v.dval + offset, repeat))
		    {
		       jdfits_error ("Read error.");
		       goto close_error;
		    }
		  break;
		  
		case JDFITS_STRING_TYPE:
		case JDFITS_BOOL_TYPE:
		case JDFITS_BIT_TYPE:
		case JDFITS_BYTE_TYPE:
		default:
		  jdfits_error ("jdfits_bintable_read_data: Field type '%d' is not implemented.", fdat->type);
		  goto close_error;
	       }
	  }
     }
   jdfits_read_close_data (ft);
   return 0;
   
   close_error:
   jdfits_read_close_data (ft);
   for (i = 0; i < tfields; i++)
     {
	if (fbt->finfo[i].data_array.v.ival != NULL)
	  {
	     SLFREE (fbt->finfo[i].data_array.v.ival);
	     fbt->finfo[i].data_array.v.ival = NULL;
	  }
     }
   return -1;
}

/*}}}*/


