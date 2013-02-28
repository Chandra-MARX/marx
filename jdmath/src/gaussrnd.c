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
# include <stdlib.h>
#endif
#include <string.h>

#include "jdmath.h"

double JDMgaussian_random (void)
{
   static int one_available = 0;
   double g1, g, s;
   static double g2;

   if (one_available)
     {
	one_available = 0;
	return g2;
     }

   one_available = 1;

   do
     {
	g1 = 2.0 * JDMrandom () - 1.0;
	g2 = 2.0 * JDMrandom () - 1.0;
	g = g1 * g1 + g2 * g2;
     }
   while ((g >= 1.0) || (g == 0.0));

   s = sqrt (-2.0 * log (g) / g);
   g2 = g2 * s;
   return g1 * s;
}

double JDMexpn_random (void)
{
   double r;

   do
     {
	r = JDMrandom ();
     }
   while (r == 0.0);
   return -log(r);
}

#if 0

int main ()
{
   unsigned int i, imax = 1000000;
   double x, x2, xi;

   x = x2 = 0.0;
   for (i = 0; i < imax; i++)
     {
	xi = JDMgaussian_random ();

	x = (i * x + xi) / (i + 1);
	x2 = (i * x2 + xi * xi) / (i + 1);

	fprintf (stdout, "%f\n", xi);
     }

   fprintf (stderr, "<x> = %f, <x^2> = %f, <x^2> - <x>^2 = %f\n",
	    x, x2, x2 - x * x);
   return 0;
}

#endif
