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

double JDMerf (double x)
{
#ifdef HAVE_ERF
   return erf(x);
#else
   int sign;
   double t;

   if (x < 0)
     {
	sign = -1;
	x = -x;
     }
   else sign = 1;

   /* This form uses equation 7.126 of A&S.  It is accurate to
    * 1.5e-7.
    */

   t = 1 / (1 + 0.3275911 * x);

   if (x < 1.0e17)
     {
	t = 1.0 - (t * (0.254829592 + t *
			(-0.284496736 + t *
			 (1.421413741 + t *
			  (-1.453152027 + t *
			   (1.061405429)))))) * exp (-x * x);
     }
   else t = 1.0;

   t = t * sign;

   return t;
#endif
}
