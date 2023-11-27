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

unsigned int JDMbinary_search_d (double x, double *xp, unsigned int n)
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
   if (x >= xp[n0]) return n1;
   return n0;
}

double JDMinterpolate_d (double x, double *xp, double *yp, unsigned int n)
{
	unsigned int n0, n1;
	double x0, x1;

	if (n == 1)
		return yp[0];
	n1 = JDMbinary_search_f(x, xp, n);
	n0 = n1 - 1;

	if (x == xp[n1])
		return yp[n1];
	if (n1 == n)
	{
		n1--;
		n0--;
	}
	if (n1 == 0)
	{
		n0 = 1;
	}

	x0 = xp[n0];
	x1 = xp[n1];
	if (x1 == x0)
		return yp[n1];

	return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * (x - x0);
}

int JDMinterpolate_dvector (double *xp, double *yp, unsigned int n,
			    double *oldxp, double *oldyp, unsigned int oldn)
{
   double y_0, y_1;
   double x_0, x_1, *xpmax, *oldxpmax;
   double x;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   /* I think that this will win if log(oldn) < oldn/n */
   if (n * 10 < oldn)
     {
	unsigned int i;

	for (i = 0; i < n; i++)
	  yp[i] = JDMinterpolate_d (xp[i], oldxp, oldyp, oldn);
	return 0;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);	       /* last value */

   oldxp++;
   oldyp++;

   /* Find out where to begin */
   x = *xp;
   while ((oldxp + 5 < oldxpmax)
	  && (x > *(oldxp + 5)))
     {
	oldxp += 5;
	oldyp += 5;
     }

   while (xp < xpmax)
     {
	double *oldxp_save = oldxp;

	x = *xp++;

	/* Move along the old axis until x is between two values. */
	while ((oldxp < oldxpmax)
	       && (x > *oldxp))
	  {
	     oldxp++;
	  }

	oldyp += (oldxp - oldxp_save);

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

/* This interpolation routine is like JDMinterpolate_dvector except that it
 * simultaneously handles n_yp vectors of y values (given by oldyp).
 */
int JDMinterpolate_n_dvector (double *xp, double **yp, unsigned int n,
			      double *oldxp, double **oldyp, unsigned int oldn,
			      unsigned int n_yp)
{
   double y_0, y_1, x_0, x_1, *xpmax, *oldxpmax, dx_10;
   double x;
   unsigned int i;
   unsigned int count, yp_count;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
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
	else
	  {
	     double factor = (x - x_0) / dx_10;

	     for (i = 0; i < n_yp; i++)
	       {
		  register double *ypptr;

		  ypptr = oldyp[i] + count;

		  y_1 = *ypptr;
		  y_0 = *(ypptr - 1);
		  /* y_0 = oldyp[i][count - 1];
		  y_1 = oldyp[i][count]; */

		  /* linear interpolation -- more generally it may be better to do:
		   * *yp++ = (*interp_fun) (x, x0, x1, y0, y1);
		   */

		  /* We have to form the test because the only thing that is assumed
		   * is that the x values are ordered.  They may not be unique, */
		  yp[i][yp_count] = y_0 + (y_1 - y_0) * factor;
	       }
	  }
	yp_count++;
     }
   return 0;
}

int JDMinterpolate_dfvector (double *xp, double *yp, unsigned int n,
			     float *oldxp, float *oldyp, unsigned int oldn)
{
   double y_0, y_1;
   double x_0, x_1, *xpmax;
   float *oldxpmax;
   double x;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);	       /* last value */

   oldxp++;
   oldyp++;

   /* Find out where to begin */
   x = *xp;
   while ((oldxp + 5 < oldxpmax)
	  && (x > *(oldxp + 5)))
     {
	oldxp += 5;
	oldyp += 5;
     }

   while (xp < xpmax)
     {
	float *oldxp_save = oldxp;

	x = *xp++;

	/* Move along the old axis until x is between two values. */
	while ((oldxp < oldxpmax)
	       && (x > *oldxp))
	  {
	     oldxp++;
	  }

	oldyp += (oldxp - oldxp_save);

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

int JDMinterpolate_n_dfvector (double *xp, double **yp, unsigned int n,
			       float *oldxp, float **oldyp, unsigned int oldn,
			       unsigned int n_yp)
{
   double y_0, y_1, x_0, x_1, dx_10;
   double *xpmax;
   float *oldxpmax;
   double x;
   unsigned int i;
   unsigned int count, yp_count;

   if (oldn < 2)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
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
	else
	  {
	     double factor = (x - x_0) / dx_10;

	     for (i = 0; i < n_yp; i++)
	       {
		  register float *ypptr;

		  ypptr = oldyp[i] + count;

		  y_1 = *ypptr;
		  y_0 = *(ypptr - 1);
		  /* y_0 = oldyp[i][count - 1];
		  y_1 = oldyp[i][count]; */

		  /* linear interpolation -- more generally it may be better to do:
		   * *yp++ = (*interp_fun) (x, x0, x1, y0, y1);
		   */

		  /* We have to form the test because the only thing that is assumed
		   * is that the x values are ordered.  They may not be unique, */
		  yp[i][yp_count] = y_0 + (y_1 - y_0) * factor;
	       }
	  }
	yp_count++;
     }
   return 0;
}
