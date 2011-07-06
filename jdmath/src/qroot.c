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
#include <stdlib.h>
#endif
#include <string.h>

#include "jdmath.h"

/* Returns -1 is a == 0,
 * returns +1 if roots are real
 * returns 0 if roots are complex.  Then roots are:
 *    x1 = rplus + i rminus;
 *    x2 = rplus - i rminus;
 */
int JDMquadratic_root (double a, double b, double c,
		       double *rplus, double *rminus)
{
   double bsqr, ac4, factor, neg_b_over_2a;

   if (a == 0.0)
     {
	if (JDMath_Error == 0)
	  JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   bsqr = b * b;
   ac4 = a * c * 4;
   neg_b_over_2a = -b / (2.0 * a);

   if (bsqr > ac4)
     {
	factor = 1.0 + sqrt (1.0 - ac4 / bsqr);
	*rplus = -2.0 * c / (b * factor);
	*rminus = neg_b_over_2a * factor;
	return 1;
     }

   if (bsqr == ac4)
     {
	*rplus = *rminus = neg_b_over_2a;
	return 1;
     }

   /* complex */
   *rplus = neg_b_over_2a;
   *rminus = sqrt (c/a) * sqrt(1.0 - bsqr/ac4);
   return 0;
}

#if 0

static int eval(double a, double b, double c)
{
   double x, y;
   int ret;

   ret = JDMquadratic_root (a, b, c, &x, &y);
   if (ret == -1)
     {
	fprintf (stderr, "returned -1.  Bad!\n");
	return -1;
     }
   if (ret == 0)
     {
	fprintf (stderr, "%f %f %f ==> COMPLEX: (%f plus/minus i * %f)\n",
		 a, b, c, x, y);
	return 0;
     }

   fprintf (stderr, "%f %f %f ==> REAL: %f, %f\n",
	    a, b, c, x, y);
   return 0;
}

int main ()
{
   if (
       eval (1, 2, 1)
       || eval (1, -2, 1)
       || eval (1, 0, 1)
       || eval (1, 1, -20)
       || eval (1, 1, -12)
       )
     return 1;

   return 0;
}

#endif
