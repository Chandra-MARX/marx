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
#include "config.h"

#include <stdio.h>
#include <math.h>


#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"

#define MAX_POINTS 25

int JDMpoly_interp (double *xdat, double *ydat, unsigned int npts,
		    double x, double *yp, double *dyp)
{
   double c[MAX_POINTS], d[MAX_POINTS];
   unsigned int i, m;
   unsigned int closest_i;
   double diff, y, dy, dx;
   
   if (npts > MAX_POINTS)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error2 ("JDMpoly_interp", "polynomial degree not supported.");
	return -1;
     }
   
   /* Find the point that is closest to x.  This point will be used below to
    * compute the error *dyp.
    */
   closest_i = 0;
   diff = fabs (x - xdat[0]);

   for (i = 0; i < npts; i++)
     {
	d[i] = c[i] = ydat[i];
	dx = fabs (x - xdat[i]);
	if (dx < diff)
	  {
	     diff = dx;
	     closest_i = i;
	  }
     }
   
   dy = y = ydat [closest_i];
   
   /* Since the first column was filled-in above, we only need to fill in 
    * npts - 1 more.
    */
   
   for (m = 1; m < npts; m++)
     {
	unsigned int imax = npts - m;
	
	for (i = 0; i < imax; i++)
	  {
	     double den, num;

	     diff = xdat[i] - x;
	     dx = xdat[i + m] - x;
	     
	     den = diff - dx;
	     if (den == 0.0)
	       {
		  JDMath_Error = JDMATH_DIVIDE_ZERO_ERROR;
		  JDMmsg_error ("JDMpoly_interp");
		  return -1;
	       }
	     
	     num = (c[i + 1] - d[i]) / den;
	     
	     c[i] = diff * num;
	     d[i] = dx * num;
	  }
	
	/* Now correct y value by this iteration's results. */
	if (closest_i != 0)
	  dy = d[--closest_i];
	else
	  dy = c[0];
	
	  y += dy;
     }
   
   *dyp = dy;
   *yp = y;
   return 0;
}

#if 0

int main (int argc, char **argv)
{
#define NPTS 10
   
   double xa[NPTS], ya[NPTS];
   double y, dy, x;
   unsigned int i;
   
   for (i = 0; i < NPTS; i++)
     {
	xa[i] = i + 0.5 * JDMRANDOM;
	ya[i] = JDMRANDOM;
	fprintf (stdout, "%f\t%f\n", xa[i], ya[i]);
     }
   
   for (x = 0; x < NPTS; x += 0.1)
     {
	(void) JDMpoly_interp (xa, ya, NPTS, x, &y, &dy);
	fprintf (stdout, "%f\t%f\t%f\n", x, y, dy);
     }
   return 0;
}
#endif
