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
/* The routines in this file are used to create IDL data files from fits files.
 * exact interpretation of the binary data depends upon what is in the header.
 * The binary data is arranged in such a way that IDL can read into IDL arrays
 * as efficiently as possible.
 */

#include "config.h"

#include <stdio.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <string.h>

#include <jdfits.h>


/* Header looks like:
 * %%%fits-to-idl%%% NFIELDS NROWS TYPE1 TYPE2 ... TYPEN
 */
static int jdfits_idl_write_header (JDFits_Type *ft)
{
   int i;
   char *s = "%%%fits-to-idl%%% ";
   char buf[256];
   Fits_Bintable_Field_Type *finfo;
   Fits_Bintable_Type *bt;
   
   bt = ft->header->ext.bintable;
   
   sprintf (buf, "%s %d %d ", s, bt->tfields, bt->naxis2);
   s = buf + strlen (buf);
   
   for (i = 0; i < bt->tfields; i++)
     {
	finfo = &bt->finfo[i];
	sprintf (s, "%d ", finfo->type);
	s += strlen (s);
     }
   
   fputs (buf, stdout);
   i = strlen (buf);
   while (i++ < 256) putc (0, stdout);
   
   return 0;
}

static int jdfits_idl_dump_table (JDFits_Type *ft)
{
   return -1;
}

static int jdfits_to_idl (char *file, char *header)
{
   JDFits_Type *ft;
   int tfields, i;
   char *name;
   unsigned int rows;
   Fits_Bintable_Type *fbt;
   
   ft = jdfits_open_file (file, JDFITS_READ_MODE);
   if (ft == NULL)
     {
	jdfits_error ("jdfits_to_idl: Unable to open %s.", file);
	return -1;
     }
   
   while (1)
     {
	if (ft->header->type == JDFITS_BINTABLE_HEADER)
	  {
	     if (-1 == jdfits_bintable_parse_headers (ft))
	       {
		  jdfits_close_file (ft);
		  return -1;
	       }
	     
	     if (header == NULL) break;/* ok */
	     
	     name = ft->header->ext.bintable->extname;
	     if ((name != NULL) && !strcmp (header, name)) break;
	  }
	
	if (jdfits_skip_to_next_header (ft) == -1) 
	  {
	     jdfits_close_file (ft);
	     return -1;
	  }
	if (jdfits_read_header (ft) == -1) 
	  {
	     jdfits_error ("jdfits_to_idl: appropriate fits header not found.");
	     jdfits_close_file (ft);
	     return -1;
	  }
     }

   /* We only get here if we have read and parsed the appropriate header. */
   
   if (-1 == jdfits_bintable_read_data (ft))
     {
	jdfits_close_file (ft);
	return -1;
     }

   fbt = ft->header->ext.bintable;
   
   rows = fbt->naxis2;
   tfields = fbt->tfields;

   jdfits_idl_write_header (ft);
   
   for (i = 0; i < tfields; i++)
     {
	unsigned int type;
	unsigned int nelements;
	Fits_Bintable_Field_Type *finfo;
	Fits_Data_Array_Type *fdat;
	
	finfo = &fbt->finfo[i];
	fdat = &(finfo->data_array);
	
	type = fdat->type;
	nelements = fdat->nelements;
	
	switch (type)
	  {
	   case JDFITS_INT16_TYPE:
#if 0
	       {
		  /* byte swap */
		  register unsigned char *p, *pmax, ch;
		  p = (unsigned char *) fdat->v.hval;
		  pmax = p + 2 * nelements;
		  while (p < pmax)
		    {
		       ch = *p;
		       *p = *(p + 1);
		       *(p + 1) = ch;
		       p += 2;
		    }
	       }
#endif
	     if (nelements != fwrite (fdat->v.hval, 2, nelements, stdout))
	       {
		  jdfits_error ("jdfits_to_idl: write error.");
		  jdfits_close_file (ft);
		  return -1;
	       }
	     break;
	     
	   case JDFITS_INT32_TYPE:
	     if (nelements != fwrite (fdat->v.lval, 4, nelements, stdout))
	       {
		  jdfits_error ("jdfits_to_idl: write error.");
		  jdfits_close_file (ft);
		  return -1;
	       }
	     break;
	     
	   case JDFITS_FLOAT32_TYPE:
	     if (nelements != fwrite (fdat->v.fval, 4, nelements, stdout))
	       {
		  jdfits_error ("jdfits_to_idl: write error.");
		  jdfits_close_file (ft);
		  return -1;
	       }
	     break;
	     
	   case JDFITS_FLOAT64_TYPE:
	     if (nelements != fwrite (fdat->v.dval, 8, nelements, stdout))
	       {
		  jdfits_error ("jdfits_to_idl: write error.");
		  jdfits_close_file (ft);
		  return -1;
	       }
	     break;
	     
	   case JDFITS_STRING_TYPE:
	   case JDFITS_BOOL_TYPE:
	   case JDFITS_BIT_TYPE:
	   case JDFITS_BYTE_TYPE:
	   default:
	     jdfits_error ("jdfits_to_idl: Field type '%d' is not implemented.", type);
	     jdfits_close_file (ft);
	     return -1;
	  }
     }
   return 0;
}




	
	
