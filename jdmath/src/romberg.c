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

void JDMint_trapezoid (double (*f)(double, JDM_User_Type *),
		       JDM_User_Type *user, double a, double b, 
		       unsigned int *itnum, double *value)
{
   
   double ba;
   double delta, x, y;
   unsigned int i, npts;
   
   npts = *itnum;
   ba = b - a;
   if (npts == 0)
     {
	*itnum = 1;
	*value = 0.5 * ba * ((*f)(a, user) + (*f)(b, user));
	return;
     }
   
   delta = ba / (double) npts;
   x = a + 0.5 * delta;
   y = 0.0;
   
   for (i = 0; i < npts; i++)
     {
	y += (*f)(x, user);
	x += delta;
     }
   
   *value = 0.5 * (*value + delta * y);
   *itnum = 2 * npts;
}


#define EPSILON 1.0e-12
#define MAX_STEPS 20
#define MAX_DEGREE_IN_HSQR 5

int JDMint_romberg (double (*f)(double, JDM_User_Type *), 
		JDM_User_Type *user, double a, double b, 
		double *result)
{
   double value [MAX_STEPS], v;
   double harray[MAX_STEPS], h;
   unsigned int i;
   unsigned int context;
   
   if (a == b) 
     {
	*result = 0.0;
	return 0;
     }

   context = 0;
   h = 1.0;

   for (i = 0; i < MAX_STEPS; i++)
     {
	/* Note: JDMtrapezoid doubles the number of intervals with each call.
	 */
	JDMint_trapezoid (f, user, a, b, &context, &v);
	harray[i] = h;
	value[i] = v;
	
	if (i >= MAX_DEGREE_IN_HSQR)
	  {
	     double y, dy;
	     
	     /* Perform a polynomial extrapolation to the h = 0 limit.
	      */
	     if (-1 == JDMpoly_interp (harray + (i - MAX_DEGREE_IN_HSQR),
				       value + (i - MAX_DEGREE_IN_HSQR),
				       MAX_DEGREE_IN_HSQR,
				       0.0,   /* seeking result at h = 0 */
				       &y, &dy))
	       {
		  return -1;
	       }
	     
	     if (fabs (dy) < EPSILON * fabs (y))
	       {
		  *result = y;
		  return 0;
	       }
	  }
	
	
	/* error polynomial will be a function of h^2 for the trapezoidal 
	 * algorithm.
	 */
	h = 0.25 * h;
     }
   
   JDMmsg_error2 ("JDMromberg", "Too many steps.");
   return -1;
}

#if 0
static int count;
double my_fun (double x, JDM_User_Type *u)
{
   count++;
   /* fprintf (stdout, "x = %f\n", x); */
   return 1+exp (x);
}


int main (int argc, char **argv)
{
   unsigned int i;
   double y;
   int context;
   
   context = 0;
   for (i = 0; i < 10; i++)
     {
	JDMint_trapezoid (my_fun, NULL, 0, 1, &context, &y);
	
	fprintf (stderr, "%20.16f\n", y - exp(1.0));
     }
   fprintf (stdout, "count = %d\n", count);
   count = 0;
   y = 0;
   
   JDMint_romberg (my_fun, NULL, 0, 1, &y);
   fprintf (stderr, "%e\n", y - exp(1.0));
   fprintf (stdout, "count = %d\n", count);
   
   return 0;
}
#endif
