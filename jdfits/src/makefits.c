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

#include <math.h>
#include "jdfits.h"

static int add_btable_extension (JDFits_Type *);
static int add_additional_headers (JDFits_Type *);

int main (int argc, char **argv)
{
   JDFits_Type *ft;
   
   ft = jdfits_open_file ("test.fits", JDFITS_WRITE_MODE);
   if (ft == NULL)
     {
	jdfits_error ("Unable to create test.fits.");
	return -1;
     }
   
   /* Write out the primary header */
   if (-1 == jdfits_init_null_primary_header (ft))
     return 1;
   
   if (-1 == add_additional_headers (ft))
     return 1;
       
   if (-1 == jdfits_end_header (ft))
     return 1;

   if (-1 == add_btable_extension (ft))
     return 1;
   
   jdfits_close_file (ft);
   return 0;
}


static int add_additional_headers (JDFits_Type *ft)
{
   if ((-1 == jdfits_write_header_string (ft, "CONTENT", "SIMULATION", NULL))
       || (-1 == jdfits_write_header_string (ft, "ORIGIN", "MIT", "Institution Name"))
       || (-1 == jdfits_write_header_string (ft, "OBSERVER", "DAVIS", NULL))
       || (-1 == jdfits_write_header_string (ft, "TELESCOP", "S-AXAF", NULL))
       || (-1 == jdfits_write_header_string (ft, "PROGNAME", "csim", NULL))
       || (-1 == jdfits_write_header_string (ft, "OBJECT", "pointsource", NULL))
       || (-1 == jdfits_write_header_comment (ft, "", "This is my comment")))
     return -1;
   
   return 0;
}

     

JDFits_BTable_Keyword_Type Binary_Table_Keywords [] = 
{
   {"X",	NULL,	"D",	NULL,	"None", NULL},
   {"SIN(X)",	NULL,	"D",	NULL,	"None",	NULL},
   {"Sign of X",NULL,	"J",	NULL,	"None",	NULL},
   {NULL, NULL, NULL, NULL, NULL}
};

static int add_btable_extension (JDFits_Type *ft)
{
   double x, y;
   int i, imax;
   
   imax = 100;

   if (-1 == jdfits_create_btable_extension (ft, 
					     Binary_Table_Keywords,
					     imax,
					     0, 
					     1,
					     "RAYTRACE"))
     return -1;
					     
		     
   if (-1 == (add_additional_headers (ft)))
     return -1;
   
   if (-1 == jdfits_end_header (ft))
     return -1;
   
   /* And the data */
   for (i = 0; i < imax; i++)
     {
	int32 s;
	x = (double)i / 10.0;
	y = sin (x);
	if (y >= 0.0) s = 1; else s = -1;
	
	if ((-1 == jdfits_write_float64 (ft, &x, 1))
	    || (-1 == jdfits_write_float64 (ft, &y, 1))
	    || (-1 == jdfits_write_int32 (ft, &s, 1)))
	  break;
     }
   
   if (-1 == jdfits_end_data (ft))
     return -1;

   return 0;
}
