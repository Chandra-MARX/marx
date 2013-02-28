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
#include <string.h>

#include "jdmath.h"
#include "_jdmath.h"

void _JDMswap_dvector (double *a, double *b, unsigned int n)
{
   double tmp;
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	tmp = a[i];
	a[i] = b[i];
	b[i] = tmp;
     }
}

/* Compute y=y+a*x */
void _JDM_add_to_vector (double *y, double *x, double a, unsigned int n)
{
   unsigned int n8;

   n8 = n / 8;
   n = n % 8;

   while (n8--)
     {
	y[0] += a * x[0];
	y[1] += a * x[1];
	y[2] += a * x[2];
	y[3] += a * x[3];
	y[4] += a * x[4];
	y[5] += a * x[5];
	y[6] += a * x[6];
	y[7] += a * x[7];

	y += 8;
	x += 8;
     }

   while (n--)
     {
	y[0] += a * x[0];
	y++;
	x++;
     }
}

double _JDM_innerprod (double *a, double *b, unsigned int n)
{
   unsigned int n8;
   double sum;

   n8 = n / 8;
   n = n % 8;
   sum = 0.0;
   while (n8--)
     {
	sum += a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
	  + a[4]*b[4] + a[5]*b[5] + a[6]*b[6] + a[7]*b[7];
	a += 8;
	b += 8;
     }

   while (n--)
     {
	sum += a[0] * b[0];
	a++;
	b++;
     }
   return sum;
}

double _JDM_innerprod_col (double *a, double **b, unsigned int col, unsigned int n)
{
   unsigned int n8;
   double sum;

   n8 = n / 8;
   n = n % 8;
   sum = 0.0;
   while (n8--)
     {
	sum += a[0]*b[0][col] + a[1]*b[1][col] + a[2]*b[2][col]
	  + a[3]*b[3][col] + a[4]*b[4][col] + a[5]*b[5][col]
	  + a[6]*b[6][col] + a[7]*b[7][col];
	a += 8;
	b += 8;
     }

   while (n--)
     {
	sum += a[0] * b[0][col];
	a++;
	b++;
     }

   return sum;
}

double _JDM_nvector_sum (double *a, unsigned int n)
{
   unsigned int i;
   double sum;

   sum = 0.0;
   for (i = 0; i < n; i++) sum += a[i];
   return sum;
}

double _JDM_nvector_max (double *a, unsigned int n)
{
   unsigned int i;
   double max;

   if (n == 0)
     return 0;

   max = a[0];
   for (i = 1; i < n; i++)
     {
	if (a[i] > max)
	  max = a[i];
     }
   return max;
}

double _JDM_nvector_max_abs (double *a, unsigned int n)
{
   unsigned int i;
   double max;

   if (n == 0)
     return -1.0;

   max = fabs (a[0]);
   for (i = 1; i < n; i++)
     {
	if (fabs(a[i]) > max)
	  max = fabs (a[i]);
     }
   return max;
}

double *_JDM_equilibrate (double **a, unsigned int n, double *eq)
{
   unsigned int i;
   int is_malloced;

   is_malloced = 0;
   if (eq == NULL)
     {
	eq = JDMdouble_vector (n + 1);
	if (eq == NULL)
	  return eq;
	is_malloced = 1;
     }

   for (i = 0; i < n; i++)
     {
	double max = _JDM_nvector_max_abs (a[i], n);
	if (max == 0)
	  {
	     JDMath_Error = JDMATH_SING_MATRIX_ERROR;
	     if (is_malloced) JDMfree_double_vector (eq);
	     return NULL;
	  }
	eq[i] = 1.0 / max;
     }

   return eq;
}

