/*
 Copyright (c) 2002,2013 John E. Davis

 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the Free
 Software Foundation; either version 2 of the License, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc., 675
 Mass Ave, Cambridge, MA 02139, USA.
*/
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"

#ifndef SLMALLOC
#define SLMALLOC malloc
#define SLFREE free
#define SLCALLOC calloc
#define SLREALLOC realloc
#endif

int JDMread_column_fdata (char *file, float **data, int *cindex, unsigned int n, unsigned int *num_read)
{
   char buf[512];
   unsigned int space = 1024;
   unsigned int nread = 0;
   int do_close;
   unsigned int i;
#define MAX_COLUMNS 20
   float dum[MAX_COLUMNS];
   FILE *fp;

   *num_read = 0;
   if (n > MAX_COLUMNS)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("JDMread_column_fdata");
	return -1;
     }

   if (file == NULL)
     {
	fp = stdin;
	do_close = 0;
	file = "<stdin>";
     }
   else
     {
	if (NULL == (fp = fopen (file, "r")))
	  {
	     JDMath_Error = JDMATH_FILE_OPEN_ERROR;
	     JDMmsg_error2 ("JDMread_column_fdata", file);
	     return -1;
	  }
	do_close = 1;
     }

   /* initialize everything to NULL */
   for (i = 0; i < n; i++) data[i] = NULL;

   for (i = 0; i < n; i++)
     {
	float *ptr;
	if (cindex[i] == 0) continue;
	ptr = (float *) SLMALLOC (space * sizeof (float));
	if (ptr == NULL) goto return_error;
	data[i] = ptr;
     }

   while (NULL != fgets (buf, sizeof (buf) - 1, fp))
     {
	if ((int)n > sscanf (buf, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
			dum+0, dum+1, dum+2, dum+3, dum+4,
			dum+5, dum+6, dum+7, dum+8, dum+9,
			dum+10, dum+11, dum+12, dum+13, dum+14,
			dum+15, dum+16, dum+17, dum+18, dum+19
			))
	  continue;

	if (nread == space)
	  {
	     space += 1024;
	     for (i = 0; i < n; i++)
	       {
		  float *ptr;
		  if (cindex[i] == 0) continue;
		  ptr = (float *) SLREALLOC (data[i], space * sizeof (float));
		  if (ptr == NULL) goto return_error;
		  data[i] = ptr;
	       }
	  }

	for (i = 0; i < n; i++)
	  {
	     if (cindex[i] == 0) continue;
	     data[i][nread] = dum[i];
	  }

	nread++;
     }

   if (nread == 0)
     {
	JDMath_Error = JDMATH_CORRUPT_FILE_ERROR;
	JDMmsg_error2 ("JDMread_column_fdata", file);
	goto read_error;
     }

   if (file != NULL) fclose (fp);
   *num_read = nread;
   return 0;

   return_error:
   JDMath_Error = JDMATH_MALLOC_ERROR;

   read_error:

   if (do_close) fclose (fp);

   JDMmsg_error2 ("JDMread_column_fdata", file);

   for (i = 0; i < n; i++)
     {
	if (data[i] != NULL) SLFREE (data[i]);
     }
   return -1;
}

int JDMread_float_xy (char *file, float **x, float **y, int colx, int coly, unsigned int *num_read)
{
   int n;
   int cindex [MAX_COLUMNS];
   float *data [MAX_COLUMNS];

   *num_read = 0;
   if (x != NULL) *x = NULL;
   if (y != NULL) *y = NULL; else coly = 2;

   colx--;
   coly--;
   if ((colx < 0) || (coly < 0) || (x == NULL)
       || (colx >= MAX_COLUMNS)
       || (coly >= MAX_COLUMNS))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("JDMread_float_xy");
	return -1;
     }

   memset ((char *) cindex, 0, sizeof (cindex));

   cindex [colx] = 1;
   n = colx;

   if (y != NULL)
     {
	cindex [coly] = 1;
	if (coly > colx) n = coly;
     }

   if (-1 == JDMread_column_fdata (file, data, cindex, n + 1, num_read))
     return -1;

   *x = data[colx];
   if (y != NULL) *y = data[coly];

   return 0;
}
