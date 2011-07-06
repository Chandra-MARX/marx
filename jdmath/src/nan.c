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
#include "jdmath.h"

unsigned int JDMcheck_nan (double *x, unsigned int n)
{
#ifdef HAVE_ISNAN
   unsigned int i;
   for (i = 0; i < n; i++)
     {
	if (isnan (x[i]))
	  return i + 1;
     }
#endif
   return 0;
}

unsigned int JDMcheck_inf (double *x, unsigned int n, int *sign)
{
#ifdef HAVE_ISINF
   unsigned int i;
   for (i = 0; i < n; i++)
     {
	if (0 == isinf (x[i]))
	  continue;
	*sign = isinf (x[i]);
	return i + 1;
     }
#endif
   return 0;
}

unsigned int JDMcheck_finite (double *x, unsigned int n)
{
#ifdef HAVE_FINITE
   unsigned int i;
   for (i = 0; i < n; i++)
     {
	if (0 == finite (x[i]))
	  return i + 1;
     }
#else
# if defined(HAVE_ISINF) && defined(HAVE_ISNAN)
   unsigned int i;
   for (i = 0; i < n; i++)
     {
	double xi = x[i];
	if (isinf (xi) || isnan (xi))
	  return i + 1;
     }
# endif
#endif
   return 0;
}

int JDMisnan (double x)
{
#ifdef HAVE_ISNAN
   return isnan (x);
#else
   return 0;
#endif
}

int JDMisinf (double x)
{
#ifdef HAVE_ISINF
   return isinf (x);
#else
   return 0;
#endif
}

int JDMfinite (double x)
{
#ifdef HAVE_FINITE
   return finite (x);
#else
# if defined(HAVE_ISINF) && defined(HAVE_ISNAN)
   return ((0 == isinf(x)) && (0 == isnan(x)));
# else
   return 1;
# endif
#endif
}
