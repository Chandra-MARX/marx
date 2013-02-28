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

int JDMlog_grid (double *x, unsigned int num, double a, double b)
{
   double ratio, factor;
   unsigned int i;

   if ((a <= 0.0) || (a > b) || (num < 2))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   ratio = b / a;
   factor = 1.0 / (num - 1);

   for (i = 0; i < num; i++)
     {
	x[i] = a * pow (ratio, i * factor);
     }

   return 0;
}

int JDMlog_grid_f (float *x, unsigned int num, float a, float b)
{
   double ratio, factor;
   unsigned int i;

   if ((a <= 0.0) || (a > b) || (num < 2))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   ratio = b / a;
   factor = 1.0 / (num - 1);

   for (i = 0; i < num; i++)
     {
	x[i] = (float) (a * pow (ratio, i * factor));
     }

   return 0;
}
