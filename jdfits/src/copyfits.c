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

int main (int argc, char **argv)
{
   JDFits_Type *ftin, *ftout;
   char *infile, *outfile;
   
   if (argc != 3)
     {
	fprintf (stderr, "Usage: copyfits input-file output-file\n");
	return -1;
     }
   
   infile = argv[1];
   outfile = argv[2];

   ftin = jdfits_open_file (infile, JDFITS_READ_MODE);
   if (ftin == NULL)
     {
	jdfits_error ("Unable to open %s.", infile);
	return -1;
     }
   
   ftout = jdfits_open_file (outfile, JDFITS_WRITE_MODE);
   if (ftout == NULL)
     {
	jdfits_error ("Unable to create %s.", outfile);
	return -1;
     }
   
   /* This is a do loop because the act of opening a fits file for read 
    * causes the header to be read.   We keep looping while there are 
    * extension headers.
    */
   do
     {
	if (-1 == jdfits_copy_header (ftin, ftout))
	  {
	     jdfits_error ("Error copying header.");
	     return -1;
	  }
	
	/* Add more keywords here....*/
	
	/* Now close the header. */
	if (-1 == jdfits_end_header (ftout))
	  {
	     jdfits_error ("Error writing header.");
	     return -1;
	  }
	
	/* copy the data section */
	if (-1 == jdfits_copy_data (ftin, ftout))
	  {
	     jdfits_error ("Error copying data.");
	     return -1;
	  }
     }
   while (0 == jdfits_read_header (ftin));
   
   jdfits_close_file (ftin);
   jdfits_close_file (ftout);
   return 0;
}

   
