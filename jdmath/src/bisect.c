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
#include "_jdmath.h"

int JDM_bisection (double (*func)(double, void *), double a, double b, void *cd, double *xp)
{
   unsigned int count = 1;
   unsigned int bisect_count = 5;
   double fa, fb;
   double x;

   if (a > b)
     {
	double tmp = a; a = b; b = tmp;
     }

   fa = (*func)(a, cd);
   fb = (*func)(b, cd);

   if (fa * fb > 0)
     {
	/* "bisection: the interval may not bracket a root: f(a)*f(b)>0" */
	*xp = a;
	return -1;
     }

   x = a;
   while (b > a)
     {
	double fx;

	if (fb == 0.0)
	  {
	     x = b;
	     break;
	  }
	if (fa == 0.0)
	  {
	     x = a;
	     break;
	  }

	if (count % bisect_count)
	  {
	     x = (a*fb - b*fa) / (fb - fa);
	     if ((x <= a) || (x >= b))
	       x = 0.5 * (a + b);
	  }
	else
	  x = 0.5 * (a + b);

	if ((x <= a) || (x >= b))
	  break;

	fx = (*func)(x, cd);
	count++;
	
	if (fx*fa < 0.0)
	  {
	     fb = fx;
	     b = x;
	  }
	else
	  {
	     fa = fx;
	     a = x;
	  }
     }
   
   *xp = x;
   return 0;
}
