/*
 Copyright (c) 2002 John E. Davis

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
#ifndef JDMHISTOGRAM_FUNCTION
/* file included for float/double */

#include "config.h"

#include <stdio.h>
#include <math.h>


#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>

#include "jdmath.h"

#define JDMHISTOGRAM_FUNCTION JDMhistogram_d
#define BINARY_SEARCH JDMbinary_search_d
#define FLOAT_TYPE double
#include "hist.c"

#define JDMHISTOGRAM_FUNCTION JDMhistogram_f
#define BINARY_SEARCH JDMbinary_search_f
#define FLOAT_TYPE float
#include "hist.c"
#else
#if 0
static unsigned int binary_search_d (double x, double *xp, unsigned int n)
{
   unsigned int n0, n1, n2;
   
   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
	n2 = (n0 + n1) / 2;
	if (xp[n2] >= x) 
	  {
	     if (xp[n2] == x) return n2;
	     n1 = n2;
	  }
	else n0 = n2;
     }
   return n0;
}
#endif
/*
 * histogram routines
 * 
 * A 1-d histogram is specified by a set of N grid points X_k and some set
 * of values y_i to be grouped into the histogram.  The bin size of the
 * nth bin is given by x_{i+1} - x_i, except for the last bin, which is
 * assumed to be of infinite width.
 */

/* If reverse_indices is NON-NULL, it is assumed to point to an array of 
 * size npts.
 */
int JDMHISTOGRAM_FUNCTION (FLOAT_TYPE *pts, unsigned int npts,
			   FLOAT_TYPE *bin_edges, unsigned int nbins,
			   unsigned int *histogram,
			   int *reverse_indices)
{
   unsigned int i;
   FLOAT_TYPE xlo;

   for (i = 0; i < nbins; i++)
     histogram[i] = 0;

   if (reverse_indices != NULL)
     for (i = 0; i < npts; i++)
       reverse_indices[i] = -1;

   if (nbins == 0)
     return 0;

   xlo = bin_edges[0];

   for (i = 0; i < npts; i++)
     {
	FLOAT_TYPE val = pts[i];
	unsigned int j;

	if (val < xlo)
	  continue;
	
	j = BINARY_SEARCH (val, bin_edges, nbins);
	histogram[j] += 1;
	
	if (reverse_indices != NULL)
	  reverse_indices[i] = (int) j;
     }
   return 0;
}

#undef FLOAT_TYPE
#undef JDMHISTOGRAM_FUNCTION
#undef BINARY_SEARCH
#endif /* JDMHISTOGRAM_FUNCTION */
