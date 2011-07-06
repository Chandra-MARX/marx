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

JDMComplex_Type JDMc_add (JDMComplex_Type z1, JDMComplex_Type z2)
{
   JDMComplex_Type z;
   z.r = z1.r + z2.r;
   z.i = z1.i + z2.i;
   return z;
}

JDMComplex_Type JDMc_mul (JDMComplex_Type z1, JDMComplex_Type z2)
{
   JDMComplex_Type z;
   z.r = z1.r * z2.r - z1.i * z2.i;
   z.i = z1.r * z2.i + z1.i * z2.r;
   return z;
}

JDMComplex_Type JDMc_sub (JDMComplex_Type z1, JDMComplex_Type z2)
{
   JDMComplex_Type z;
   z.r = z1.r - z2.r;
   z.i = z1.i - z2.i;
   return z;
}

JDMComplex_Type JDMc_div (JDMComplex_Type z1, JDMComplex_Type z2)
{
   JDMComplex_Type z;
   double r1, r2, i1, i2;
   double ratio, denom;

   r1 = z1.r; i1 = z1.i;
   r2 = z2.r; i2 = z2.i;

   /* DO it this way to avoid overflow in the denomenator */
   if (fabs(r2) > fabs(i2))
     {
	ratio = i2 / r2;
	denom = r2 + i2 * ratio;
	z.r = (r1 + ratio * i1) / denom;
	z.i = (i1 - r1 * ratio) / denom;
     }
   else
     {
	ratio = r2 / i2;
	denom = r2 * ratio + i2;
	z.r = (r1 * ratio + i1) / denom;
	z.i = (i1 * ratio - r1) / denom;
     }
   return z;
}

JDMComplex_Type JDMc_exp (JDMComplex_Type z1)
{
   JDMComplex_Type z;
   double amp = exp (z1.r);
   z.r = amp * cos (z1.i);
   z.i = amp * sin (z1.i);
   return z;
}

/* exp (i*a*z) */

JDMComplex_Type JDMc_iexp (double a, JDMComplex_Type z1)
{
   JDMComplex_Type z;
   double theta = a * z1.r;
   double amp = exp (-a * z1.i);
   z.r = amp * cos (theta);
   z.i = amp * sin (theta);
   return z;
}

JDMComplex_Type JDMc_smul (double a, JDMComplex_Type z1)
{
   JDMComplex_Type z;
   z.r = z1.r * a;
   z.i = z1.i * a;
   return z;
}

/* Multiple by imaginary scalar */
JDMComplex_Type JDMc_imul (double a, JDMComplex_Type z1)
{
   JDMComplex_Type z;
   z.r = -z1.i * a;
   z.i = z1.r * a;
   return z;
}

JDMComplex_Type JDMc_sadd (double a, JDMComplex_Type z1)
{
   JDMComplex_Type z;
   z.r = z1.r + a;
   z.i = z1.i;
   return z;
}

JDMComplex_Type JDMComplex (double a, double b)
{
   JDMComplex_Type z;
   z.r = a;
   z.i = b;
   return z;
}

JDMComplex_Type JDMc_a_bz (double a, double b, JDMComplex_Type z1)
{
   JDMComplex_Type z;
   z.r = a + b * z1.r;
   z.i = b * z1.i;
   return z;
}

JDMComplex_Type JDMc_az1_bz2 (double a, JDMComplex_Type z1,
			      double b, JDMComplex_Type z2)
{
   JDMComplex_Type z;
   z.r = a * z1.r + b * z2.r;
   z.i = a * z1.i + b * z2.i;
   return z;
}

double JDMc_abs (JDMComplex_Type z)
{
   double r, i, fr, fi;
   double ratio;

   r = z.r;
   i = z.i;

   fr = fabs(r);
   fi = fabs(i);

   if (fr > fi)
     {
	ratio = i / r;
	return fr * sqrt (1.0 + ratio * ratio);
     }

   if (fi == 0.0) return 0.0;

   ratio = r / i;
   return fi * sqrt (1.0 + ratio * ratio);
}

void JDMc_inc (JDMComplex_Type *z, JDMComplex_Type dz)
{
   z->r += dz.r;
   z->i += dz.i;
}

int JDMc_eqs (JDMComplex_Type a, JDMComplex_Type b)
{
   return ((a.r == b.r) && (a.i == b.i));
}

int JDMc_neqs (JDMComplex_Type a, JDMComplex_Type b)
{
   return ((a.r != b.r) || (a.i != b.i));
}

double JDMc_real (JDMComplex_Type a)
{
   return a.r;
}

double JDMc_imag (JDMComplex_Type a)
{
   return a.i;
}

JDMComplex_Type JDMc_sqrt (JDMComplex_Type a)
{
   double r;

   r = JDMc_abs (a);

   if (r == 0.0) return a;

   if (a.r >= 0.0)
     {
	a.r = sqrt (0.5 * (r + a.r));
	a.i = 0.5 * a.i / a.r;
     }
   else
     {
	r = sqrt (0.5 * (r - a.r));
	a.r = 0.5 * a.i / r;
	a.i = r;

	if (a.r < 0.0)
	  {
	     a.r = -a.r;
	     a.i = -a.i;
	  }
     }
   return a;
}

#if 0

static int test_complex_sqrt (void)
{
   unsigned int n = 100000;
   double x, y;
   JDMComplex_Type z, z1;
   double r = 3.0, theta;
   double epsilon = 1.0e-10;

   while (n)
     {
	double diff;
	n--;

	theta = 2.0 * PI * JDMrandom ();
	x = r * cos (theta);
	y = r * sin (theta);

	z = JDMc_sqrt (JDMComplex (x, y));

	z1.r = (z.r * z.r - z.i * z.i);
	z1.i = (2 * z.r * z.i);

        z1.r -= x;
	z1.i -= y;
	diff = JDMc_abs (z1);

	if ((z.r < 0.0)
	    || ((z.r == 0.0) && (z.i < 0.0))
	    || (diff > epsilon))
	  {
	     fprintf (stderr, "***Error: JDMc_sqrt(%e,%e) ==> (%e,%e)\n",
		      x,y,z.r,z.i);
	     fprintf (stderr, "   Difference is: %e\n", diff);
	     diff = JDMc_abs (z1);
	  }
     }
   return 0;
}

int main ()
{
   if (-1 == test_complex_sqrt ())
     return 1;

   return 0;
}

#endif

