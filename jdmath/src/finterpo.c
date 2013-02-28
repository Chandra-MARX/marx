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
#include <math.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>

#include "jdmath.h"

unsigned int JDMbinary_search_f (float x, float *xp, unsigned int n)
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
	else
	  {
	     n0 = n2;
	  }
     }
   return n0;
}

float JDMinterpolate_f (float x, float *xp, float *yp, unsigned int n)
{
   unsigned int n0, n1;
   double x0, x1;

   n0 = JDMbinary_search_f (x, xp, n);

   x0 = xp[n0];
   n1 = n0 + 1;

   if (x == x0)
     return yp[n0];
   if (n1 == n)
     {
	if (n0 == 0)
	  return yp[n0];
	n1 = n0 - 1;
     }

   x1 = xp[n1];
   if (x1 == x0) return yp[n0];

   return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * (x - x0);
}

int JDMinterpolate_fvector (float *xp, float *yp, unsigned int n,
			float *oldxp, float *oldyp, unsigned int oldn)
{
   double y_0, y_1, x_0, x_1, x;
   float *xpmax, *oldxpmax;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("JDMinterpolate_fvector");
	return -1;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);	       /* last value */

   oldxp++;
   oldyp++;

   while (xp < xpmax)
     {
	x = *xp++;

	/* Move along the old axis until x is between two values. */
	while ((oldxp < oldxpmax)
	       && (x > *oldxp))
	  {
	     oldxp++; oldyp++;
	  }

	y_0 = *(oldyp - 1);
	y_1 = *oldyp;
	x_0 = *(oldxp - 1);
	x_1 = *oldxp;

	/* linear interpolation -- more generally it may be better to do:
	 * *yp++ = (*interp_fun) (x, x0, x1, y0, y1);
	 */

	/* We have to form the test because the only thing that is assumed
	 * is that the x values are ordered.  They may not be unique, */
	if (x_1 == x_0) *yp++ = y_0;
	else *yp++ = y_0 + (y_1 - y_0) * (x - x_0) / (x_1 - x_0);
     }
   return 0;
}

int JDMinterpolate_n_fvector (float *xp, float **yp, unsigned int n,
			      float *oldxp, float **oldyp, unsigned int oldn,
			      unsigned int n_yp)
{
   double y_0, y_1, x_0, x_1, dx_10, x;
   float *xpmax, *oldxpmax;
   unsigned int i;
   unsigned int count, yp_count;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("JDMinterpolate_n_fvector");
	return -1;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);	       /* last value */

   oldxp++;
   count = 1;
   yp_count = 0;

   while (xp < xpmax)
     {
	x = *xp++;

	/* Move along the old axis until x is between two values. */
	while ((oldxp < oldxpmax)
	       && (x > *oldxp))
	  {
	     oldxp++;
	     count++;
	  }

	x_0 = *(oldxp - 1);
	x_1 = *oldxp;
	dx_10 = (x_1 - x_0);
	if (dx_10 == 0.0)
	  {
	     for (i = 0; i < n_yp; i++)
	       yp[i][yp_count] = oldyp[i][count - 1];
	  }
	else for (i = 0; i < n_yp; i++)
	  {
	     y_0 = oldyp[i][count - 1];
	     y_1 = oldyp[i][count];

	     /* linear interpolation -- more generally it may be better to do:
	      * *yp++ = (*interp_fun) (x, x0, x1, y0, y1);
	      */

	     /* We have to form the test because the only thing that is assumed
	      * is that the x values are ordered.  They may not be unique, */
	     yp[i][yp_count] = y_0 + (y_1 - y_0) * (x - x_0) / dx_10;
	  }
	yp_count++;
     }
   return 0;
}

/* Log versions of the above */

float JDMlog_interpolate_f (float x, float *xp, float *yp, unsigned int n)
{
   unsigned int n0, n1;
   double x0, x1;

   n0 = JDMbinary_search_f (x, xp, n);

   x0 = xp[n0];
   n1 = n0 + 1;

   if ((x == x0) || (n1 == n)) return yp[n0];

   x1 = xp[n1];
   if (x1 == x0) return yp[n0];

   return yp[n0] + (yp[n1] - yp[n0]) * (log(x/x0) / log (x1/x0));
}

int JDMlog_interpolate_fvector (float *xp, float *yp, unsigned int n,
			float *oldxp, float *oldyp, unsigned int oldn)
{
   double x, y_0, y_1, x_0, x_1;
   float *xpmax, *oldxpmax;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("JDMinterpolate_fvector");
	return -1;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);	       /* last value */

   oldxp++;
   oldyp++;

   while (xp < xpmax)
     {
	x = *xp++;

	/* Move along the old axis until x is between two values. */
	while ((oldxp < oldxpmax)
	       && (x > *oldxp))
	  {
	     oldxp++; oldyp++;
	  }

	y_0 = *(oldyp - 1);
	y_1 = *oldyp;
	x_0 = *(oldxp - 1);
	x_1 = *oldxp;

	/* linear interpolation -- more generally it may be better to do:
	 * *yp++ = (*interp_fun) (x, x0, x1, y0, y1);
	 */

	/* We have to form the test because the only thing that is assumed
	 * is that the x values are ordered.  They may not be unique, */
	if (x_1 == x_0) *yp++ = y_0;
	else *yp++ = y_0 + (y_1 - y_0) * (log(x/x_0) / log(x_1/x_0));
     }
   return 0;
}

int JDMlog_interpolate_n_fvector (float *xp, float **yp, unsigned int n,
				  float *oldxp, float **oldyp, unsigned int oldn,
				  unsigned int n_yp)
{
   double y_0, y_1, x_0, x_1, x, dx_10;
   float *xpmax, *oldxpmax;
   unsigned int i;
   unsigned int count, yp_count;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("JDMinterpolate_n_fvector");
	return -1;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);	       /* last value */

   oldxp++;
   count = 1;
   yp_count = 0;

   while (xp < xpmax)
     {
	x = *xp++;

	/* Move along the old axis until x is between two values. */
	while ((oldxp < oldxpmax)
	       && (x > *oldxp))
	  {
	     oldxp++;
	     count++;
	  }

	x_0 = *(oldxp - 1);
	x_1 = *oldxp;
	if (x_0 == x_1)
	  {
	     for (i = 0; i < n_yp; i++)
	       yp[i][yp_count] = oldyp[i][count - 1];
	  }
	else
	  {
	     dx_10 = log (x_1/x_0);

	     for (i = 0; i < n_yp; i++)
	       {
		  y_0 = oldyp[i][count - 1];
		  y_1 = oldyp[i][count];

		  yp[i][yp_count] = y_0 + (y_1 - y_0) * log(x/x_0) / dx_10;
	       }
	  }
	yp_count++;
     }
   return 0;
}
