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

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLCALLOC calloc
# define SLREALLOC realloc
# define SLFREE free
#endif

int jdfits_init_null_primary_header (JDFits_Type *ft)
{
   if ((-1 == jdfits_write_header_logical (ft, "SIMPLE", 1, "FITS STANDARD"))
       || (-1 == jdfits_write_header_integer (ft, "BITPIX", 8, "Binary Data"))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS", 0, "No DATA image array present"))
       || (-1 == jdfits_write_header_logical (ft, "EXTEND", 1, "There MAY be standard extensions")))
     return -1;

  return 0;
}

int jdfits_create_btable_extension (JDFits_Type *ft,
				    JDFits_BTable_Keyword_Type *bkw,
				    int naxis2, int pcount, int gcount,
				    char *extname)
{
   int tfields, size, naxis1;
   unsigned char *p;
   unsigned char ch;
   int repeat;
   JDFits_BTable_Keyword_Type *k;

   naxis1 = tfields = 0;
   k = bkw;

   while (k->tform != NULL)
     {
	tfields++;

	repeat = 0;
	p = (unsigned char *) k->tform;

	ch = *p;
	if (0 == isdigit (ch))
	  {
	     repeat = 1;
	  }
	else
	  {
	     repeat = 0;
	     while (((ch = *p) != 0) && isdigit (ch))
	       {
		  p++;
		  repeat = 10 * repeat + (ch - '0');
	       }
	  }

	switch (ch)
	  {
	   case 'A':
	   case 'B':
	   case 'L':
	     size = 1;
	     break;

	   case 'I':
	     size = 2;
	     break;

	   case 'E':
	   case 'J':
	     size = 4;
	     break;

	   case 'D':
	     size = 8;
	     break;

	   case 'X':
	     switch (repeat)
	       {
		case 8:
		  size = 1;
		  break;
		case 16:
		  size = 2;
		  break;
		case 32:
		  size = 4;
		  break;
		default:
		  jdfits_error ("jdfits_create_btable_extension: type %dX not supported", repeat);
		  return -1;
	       }
	     repeat = 1;
	     break;

	   default:
	     jdfits_error ("Type %c not implemented.", ch);
	     return -1;
	  }
	naxis1 += repeat * size;
	k++;
     }

   if ((-1 == jdfits_write_header_string (ft, "XTENSION", "BINTABLE", "FITS BINARY TABLE"))
       || (-1 == jdfits_write_header_integer (ft, "BITPIX", 8, "Binary data"))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS", 2, "Table is a matrix"))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS1", naxis1, "Width of table in bytes"))
       || (-1 == jdfits_write_header_integer (ft, "NAXIS2", naxis2, "Number of entries in table"))
       || (-1 == jdfits_write_header_integer (ft, "PCOUNT", pcount, "Random parameter count"))
       || (-1 == jdfits_write_header_integer (ft, "GCOUNT", gcount, "Group count"))
       || (-1 == jdfits_write_header_integer (ft, "TFIELDS", tfields, "Number of fields in each row")))
     {
	return -1;
     }

   k = bkw;

   tfields = 0;
   while (k->tform != NULL)
     {
	char tform[12];
	char ttype[12];
	char tunit[12];
	char tlmin[12];
	char tlmax[12];
	char *comment, *comment1;

	tfields++;

	sprintf (tform, "TFORM%d", tfields);
	sprintf (ttype, "TTYPE%d", tfields);
	sprintf (tunit, "TUNIT%d", tfields);
	sprintf (tlmin, "TLMIN%d", tfields);
	sprintf (tlmax, "TLMAX%d", tfields);

	if (NULL == (comment = k->tform_comment))
	  comment = "Data type for field";
	if (-1 == jdfits_write_header_string (ft, tform, k->tform, comment))
	  return -1;

	if (k->ttype != NULL)
	  {
	     if (NULL == (comment = k->ttype_comment))
	       comment = "Label for field";
	     if (-1 == jdfits_write_header_string (ft, ttype, k->ttype, comment))
	       return -1;
	  }

	if (k->tunit != NULL)
	  {
	     if (NULL == (comment = k->tunit_comment))
	       comment = "Physical units for field";

	     if (-1 == jdfits_write_header_string (ft, tunit, k->tunit, comment))
	       return -1;
	  }

	if (NULL == (comment = k->min_comment))
	  comment = "Min Value for field";
	if (NULL == (comment1 = k->min_comment))
	  comment1 = "Max Value for field";

	switch (k->min_max_type)
	  {
	   case 'A':
	   default:
	     break;

	   case 'I':
	   case 'J':
	     if (-1 == jdfits_write_header_integer (ft, tlmin, k->min_value.j_val, comment))
	       return -1;
	     if (-1 == jdfits_write_header_integer (ft, tlmax, k->max_value.j_val, comment1))
	       return -1;
	     break;

	   case 'E':
	   case 'D':
	     if (-1 == jdfits_write_header_double (ft, tlmin, k->min_value.d_val, comment))
	       return -1;
	     if (-1 == jdfits_write_header_double (ft, tlmax, k->max_value.d_val, comment1))
	       return -1;
	     break;

	  }

	if (k->ctype != NULL)
	  {
	     char buf[12];

	     sprintf (buf, "TCTYP%d", tfields);
	     if (-1 == jdfits_write_header_string (ft, buf, k->ctype, NULL))
	       return -1;

	     sprintf (buf, "TCRVL%d", tfields);
	     if (-1 == jdfits_write_header_double (ft, buf, k->crval, NULL))
	       return -1;

	     sprintf (buf, "TCDLT%d", tfields);
	     if (-1 == jdfits_write_header_double (ft, buf, k->cdelt, NULL))
	       return -1;

	     sprintf (buf, "TCRPX%d", tfields);
	     if (-1 == jdfits_write_header_double (ft, buf, k->crpix, NULL))
	       return -1;
	  }
	k++;
     }

   if (-1 == jdfits_write_header_string (ft, "EXTNAME", extname, "Table name"))
     return -1;

   return 0;
}

