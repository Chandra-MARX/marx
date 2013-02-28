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
#include "_jdmath.h"

struct _JDM_Bilinear_Interp_Type
{
   double inv_dx, inv_dy;	       /* Inverse length of sizes */
   double x,y;			       /* lower left corner */
   unsigned int num_values;
   double *f0;
   double *f1;
   double *f2;
   double *f3;
};

JDM_Bilinear_Interp_Type *JDM_bilinear_open (double x, double y,
					     double dx, double dy,
					     unsigned int num,
					     double *f_ll, double *f_lr, double *f_ul, double *f_ur)
{
   unsigned int i;
   JDM_Bilinear_Interp_Type *b;

   if ((dx == 0.0) || (dy == 0.0) || (num == 0))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return NULL;
     }

   b = (JDM_Bilinear_Interp_Type *)_JDMmalloc (sizeof (JDM_Bilinear_Interp_Type), "JDM_bilinear_open");
   if (b == NULL)
     return b;

   memset ((char *) b, 0, sizeof (JDM_Bilinear_Interp_Type));

   if ((NULL == (b->f0 = JDMdouble_vector (num)))
       || (NULL == (b->f1 = JDMdouble_vector (num)))
       || (NULL == (b->f2 = JDMdouble_vector (num)))
       || (NULL == (b->f3 = JDMdouble_vector (num))))
     {
	JDM_bilinear_close (b);
	return NULL;
     }

   b->x = x;
   b->y = y;
   b->inv_dx = 1.0 / dx;
   b->inv_dy = 1.0 / dy;

   for (i = 0; i < num; i++)
     {
	double lr, ll, ul, ur;

	lr = f_lr[i];
	ll = f_ll[i];
	ul = f_ul[i];
	ur = f_ur[i];

	lr -= ll;
	ur -= ll;
	ul -= ll;

	b->f0[i] = ll;
	b->f1[i] = lr;
	b->f2[i] = ul;
	b->f3[i] = (ur - lr - ul);
     }

   b->num_values = num;

   return b;
}

int JDM_bilinear_interp (JDM_Bilinear_Interp_Type *b, double x, double y,
			 double *result)
{
   unsigned int i, num_values;

   if (b == NULL)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }
   x = b->inv_dx * (x - b->x);
   y = b->inv_dy * (y - b->y);

   num_values = b->num_values;

   for (i = 0; i < num_values; i++)
     result[i] = b->f0[i] + x * (b->f1[i] + y * b->f3[i]) + y * b->f2[i];

   return 0;
}

void JDM_bilinear_close (JDM_Bilinear_Interp_Type *b)
{
   if (b == NULL)
     return;

   _JDMfree ((char *) b->f0);
   _JDMfree ((char *) b->f1);
   _JDMfree ((char *) b->f2);
   _JDMfree ((char *) b->f3);
   _JDMfree ((char *) b);
}
