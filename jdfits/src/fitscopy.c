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


int jdfits_copy_header (JDFits_Type *ftin, JDFits_Type *ftout)
{
   JDFits_Read_Keyword_Type *r;
   JDFits_Keyword_Type *k;

   if (NULL == (r = jdfits_open_keywords (ftin)))
     {
	jdfits_error ("Unable to open keyword section.");
	return -1;
     }
   
   while (NULL != (k = jdfits_read_keyword (r)))
     {
	if (-1 == jdfits_write_keyword (ftout, k))
	  {
	     jdfits_error ("Error writing keyword: %s.", k->name);
	     jdfits_close_keywords (r);
	     return -1;
	  }
     }
   jdfits_close_keywords (r);
   
   return 0;
}

int jdfits_copy_data (JDFits_Type *ftin, JDFits_Type *ftout)
{
   unsigned int nread;
   unsigned char buf[JDFITS_RECORD_SIZE];

   if (((ftout->mode & JDFITS_WRITE_MODE) == 0)
       || ((ftin->mode & JDFITS_READ_MODE) == 0))
     {
	jdfits_error ("jdfits_copy_data: file is open in incorrect mode.");
	return -1;
     }
   
   if (-1 == jdfits_read_open_data (ftin))
     {
	jdfits_error ("jdfits_copy_data: Unable to open data section.");
	return -1;
     }
   
   while (0 != (nread = jdfits_read_bytes (ftin, buf, sizeof(buf))))
     {
	if (-1 == jdfits_write (ftout, buf, nread))
	  {
	     jdfits_error ("jdfits_copy_data: write error.");
	     return -1;
	  }
     }
   return jdfits_end_data (ftout);
}
