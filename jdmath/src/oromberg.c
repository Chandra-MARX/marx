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

/* First time through, set *itnum to 0.  On subsequent iterations, do NOT
 * modify *itnum nor *value.
 */
void JDMint_midpoint (double (*f)(double, JDM_User_Type *),
		      JDM_User_Type *user, double a, double b, 
		      unsigned int *itnum, double *value)
{
   double ba;
   double delta1, delta2, x, y;
   unsigned int i, npts;	       /* actually half as many points */
   
   npts = *itnum;
   ba = b - a;
   if (npts == 0)
     {
	*itnum = 1;
	*value = ba * (*f)(0.5 * (a + b), user);
	return;
     }
   
   delta1 = ba / (3.0 * (double) npts);
   delta2 = 2.0 * delta1;
   
   x = a + 0.5 * delta1;
   y = 0.0;
   
   for (i = 0; i < npts; i++)
     {
	y += (*f)(x, user);
	x += delta2;
	y += (*f)(x, user);
	x += delta1;
     }
   
   *value = (*value/3.0 + delta1 * y);
   *itnum = 3 * npts;
}

/* This routine is like the JDMint_midpoint except it makes a change of
 * variable x->1/x.
 */
void JDMint_midpoint_inf (double (*f)(double, JDM_User_Type *),
			  JDM_User_Type *user, double a1, double b1,
			  unsigned int *itnum, double *value)
{
   double ba;
   double delta1, delta2, x, y;
   unsigned int i, npts;	       /* actually half as many points */
   double a, b, t;
   
   /* change of variables. */
   b = 1.0 / a1;
   a = 1.0 / b1;
   
   npts = *itnum;
   ba = b - a;
   if (npts == 0)
     {
	*itnum = 1;
	t = 2.0 / (a + b);
	*value = ba * t * t * (*f)(t, user);
	return;
     }
   
   delta1 = ba / (3.0 * (double) npts);
   delta2 = 2.0 * delta1;
   
   x = a + 0.5 * delta1;
   y = 0.0;
   
   for (i = 0; i < npts; i++)
     {
	t = 1.0 / x;
	y += t * t * (*f)(t, user);
	x += delta2;
	t = 1.0 / x;
	y += t * t * (*f)(t, user);
	x += delta1;
     }
   
   *value = (*value/3.0 + delta1 * y);
   *itnum = 3 * npts;
}


#define EPSILON 1.0e-6
#define MAX_STEPS 14
#define MAX_DEGREE_IN_HSQR 5

static int oromberg (double (*f)(double, JDM_User_Type *),
		     JDM_User_Type *user, double a, double b,
		     void (*intf)(double (*)(double, JDM_User_Type *),
				  JDM_User_Type *, double, double,
				  unsigned int *, double *),
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
	/* Note: JDMint_midpoint TRIPLES the number of intervals with each call.
	 */
	(*intf)(f, user, a, b, &context, &v);
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
	 * algorithm.  The open routines triple the number of intervals.
	 */
	h = h / 9.0;
     }
   
   JDMmsg_error2 ("JDMrombergo", "Too many steps.");
   return -1;
}

int JDMint_oromberg (double (*f)(double, JDM_User_Type *), 
		     JDM_User_Type *user, double a, double b,
		     int infinite_flag,
		     double *result)
{
   if (infinite_flag)
     {
	if ((a == 0.0) || (b == 0.0)
	    || ((a > 0.0) && (b < 0.0))
	    || ((a < 0.0) && (b > 0.0)))
	  {
	     JDMath_Error = JDMATH_INVALID_PARAMETER;
	     JDMmsg_error2 ("JDMint_oromberg", 
			    "endpoint must have same sign.");
	     return -1;
	  }
	
	return oromberg (f, user, a, b, JDMint_midpoint_inf, result);
     }

   return oromberg (f, user, a, b, JDMint_midpoint, result);
}


#if 0
double my_fun (double x, JDM_User_Type *u)
{
   return 1.0 + exp (x);
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
	
	fprintf (stderr, "%20.16f", y - exp(1.0));
     }
   
   y = 0;
   
   JDMint_romberg (my_fun, NULL, 0, 1, &y);
   fprintf (stderr, "%20.16f\n", y - exp(1.0));
   
   return 0;
}
#endif
