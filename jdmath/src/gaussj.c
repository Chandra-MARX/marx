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
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"

/* Important note:  Occasionally it is necessary to swap rows.  This may
 * usually be accomplished by swapping pointers, however this will produce
 * incorrect results if the row pointers is imposed on a linear array.
 * So, _JDMswap_dvector is called.
 */

int JDMgauss_jordan_n (double **a, double **b, unsigned int n, unsigned int m)
{
   int *piv, *row, *col;
   unsigned int i;

   if ((n == 0) || (m == 0))
     return 0;

   if ((a == NULL) || (b == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   if (NULL == (piv = (int *)_JDMmalloc (3 * n * sizeof (int), "JDMgauss_jordan")))
     return -1;

   row = piv + n;
   col = row + n;

   memset ((char *) piv, 0, n * sizeof (int));

   for (i = 0; i < n; i++)
     {
	unsigned int j;
	double largest = 0;
	double *a_piv, *b_piv;
	double inv;
	unsigned int piv_row = 0, piv_col = 0;

	for (j = 0; j < n; j++)
	  {
	     unsigned int k;

	     if (piv[j] == 1)
	       continue;

	     for (k = 0; k < n; k++)
	       {
		  if (piv[k] != 0)
		    {
		       if (piv[k] == 1)
			 continue;

		       goto return_error;
		    }

		  if (fabs (a[j][k]) < largest)
		    continue;

		  largest = fabs (a[j][k]);
		  piv_row = j;
		  piv_col = k;
	       }
	  }

	piv[piv_col] += 1;
	if (piv_row != piv_col)
	  _JDMswap_dvector (a[piv_row], a[piv_col], n);

	a_piv = a[piv_col];
	b_piv = b[piv_col];

	inv = a_piv[piv_col];
	if (inv == 0)
	  goto return_error;

	inv = 1/inv;
	a_piv[piv_col] = 1;

	row[i] = piv_row;
	col[i] = piv_col;

	for (j = 0; j < n; j++) a_piv[j] = a_piv[j] * inv;
	for (j = 0; j < m; j++) b_piv[j] = b_piv[j] * inv;

	for (j = 0; j < n; j++)
	  {
	     unsigned int k;
	     double val;
	     double *aj;

	     if (j == piv_col) continue;

	     aj = a[j];
	     val = aj[piv_col];
	     aj[piv_col] = 0;

	     for (k = 0; k < n; k++) aj[k] -= a_piv[k] * val;
	     aj = b[j]; for (k = 0; k < m; k++) aj[k] -= b_piv[k] * val;
	  }
     }

   i = n;
   while (i)
     {
	unsigned int k;

	i--;
	if (row[i] == col[i]) continue;
	for (k = 0; k < n; k++)
	  {
	     double tmp = a[k][row[i]];
	     a[k][row[i]] = a[k][col[i]];
	     a[k][col[i]] = tmp;
	  }
     }

   _JDMfree ((char *)piv);
   return 0;

   return_error:
   JDMmsg_error ("JDMgauss_jordan: singular matrix");
   _JDMfree ((char *)piv);
   return -1;
}

int JDMgauss_jordan (double **a, double *b, unsigned int n)
{
   int *piv, *row, *col;
   unsigned int i;

   if (n == 0)
     return 0;

   if ((a == NULL) || (b == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   if (NULL == (piv = (int *)_JDMmalloc (3 * n * sizeof (int), "JDMgauss_jordan")))
     return -1;

   row = piv + n;
   col = row + n;

   memset ((char *) piv, 0, n * sizeof (int));

   for (i = 0; i < n; i++)
     {
	unsigned int j;
	double largest = 0;
	double *a_piv;
	double inv;
	unsigned int piv_row = 0, piv_col = 0;

	for (j = 0; j < n; j++)
	  {
	     unsigned int k;

	     if (piv[j] == 1)
	       continue;

	     for (k = 0; k < n; k++)
	       {
		  if (piv[k] != 0)
		    {
		       if (piv[k] == 1)
			 continue;

		       goto return_error;
		    }

		  if (fabs (a[j][k]) < largest)
		    continue;

		  largest = fabs (a[j][k]);
		  piv_row = j;
		  piv_col = k;
	       }
	  }

	piv[piv_col] += 1;
	if (piv_row != piv_col)
	  {
	     double t;
	     t = b[piv_row]; b[piv_row] = b[piv_col]; b[piv_col] = t;
	     _JDMswap_dvector (a[piv_row], a[piv_col], n);
	  }

	a_piv = a[piv_col];

	inv = a_piv[piv_col];
	if (inv == 0)
	  goto return_error;

	inv = 1/inv;
	a_piv[piv_col] = 1;

	row[i] = piv_row;
	col[i] = piv_col;

	for (j = 0; j < n; j++) a_piv[j] = a_piv[j] * inv;
	b[piv_col] = b[piv_col] * inv;

	for (j = 0; j < n; j++)
	  {
	     unsigned int k;
	     double val;
	     double *aj;

	     if (j == piv_col) continue;

	     aj = a[j];
	     val = aj[piv_col];
	     aj[piv_col] = 0;

	     for (k = 0; k < n; k++) aj[k] -= a_piv[k] * val;
	     b[j] -= b[piv_col] * val;
	  }
     }

   i = n;
   while (i)
     {
	unsigned int k;

	i--;
	if (row[i] == col[i]) continue;
	for (k = 0; k < n; k++)
	  {
	     double tmp = a[k][row[i]];
	     a[k][row[i]] = a[k][col[i]];
	     a[k][col[i]] = tmp;
	  }
     }

   _JDMfree ((char *)piv);
   return 0;

   return_error:
   JDMmsg_error ("JDMgauss_jordan: singular matrix");
   _JDMfree ((char *)piv);
   return -1;
}

int JDMgauss_jordan_inverse (double **a, unsigned int n)
{
   int *piv, *row, *col;
   unsigned int i;

   if (n == 0)
     return 0;

   if (a == NULL)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   if (NULL == (piv = (int *)_JDMmalloc (3 * n * sizeof (int), "JDMgauss_jordan")))
     return -1;

   row = piv + n;
   col = row + n;

   memset ((char *) piv, 0, n * sizeof (int));

   for (i = 0; i < n; i++)
     {
	unsigned int j;
	double largest = 0;
	double *a_piv;
	double inv;
	unsigned int piv_row = 0, piv_col = 0;

	for (j = 0; j < n; j++)
	  {
	     unsigned int k;

	     if (piv[j] == 1)
	       continue;

	     for (k = 0; k < n; k++)
	       {
		  if (piv[k] != 0)
		    {
		       if (piv[k] == 1)
			 continue;

		       goto return_error;
		    }

		  if (fabs (a[j][k]) < largest)
		    continue;

		  largest = fabs (a[j][k]);
		  piv_row = j;
		  piv_col = k;
	       }
	  }

	piv[piv_col] += 1;
	if (piv_row != piv_col)
	  _JDMswap_dvector (a[piv_row], a[piv_col], n);

	a_piv = a[piv_col];

	inv = a_piv[piv_col];
	if (inv == 0)
	  goto return_error;

	inv = 1/inv;
	a_piv[piv_col] = 1;

	row[i] = piv_row;
	col[i] = piv_col;

	for (j = 0; j < n; j++) a_piv[j] = a_piv[j] * inv;

	for (j = 0; j < n; j++)
	  {
	     unsigned int k;
	     double val;
	     double *aj;

	     if (j == piv_col) continue;

	     aj = a[j];
	     val = aj[piv_col];
	     aj[piv_col] = 0;

	     for (k = 0; k < n; k++) aj[k] -= a_piv[k] * val;
	  }
     }

   i = n;
   while (i)
     {
	unsigned int k;

	i--;
	if (row[i] == col[i]) continue;
	for (k = 0; k < n; k++)
	  {
	     double tmp = a[k][row[i]];
	     a[k][row[i]] = a[k][col[i]];
	     a[k][col[i]] = tmp;
	  }
     }

   _JDMfree ((char *)piv);
   return 0;

   return_error:
   JDMmsg_error ("JDMgauss_jordan: singular matrix");
   _JDMfree ((char *)piv);
   return -1;
}

#ifdef TEST_GAUSSJ

int main (int argc, char **argv)
{
   unsigned int n;
   unsigned int i, j;
   double **a, *b, **a1, *b1;

   n = 50;

   a = JDMdouble_matrix (n, n);
   a1 = JDMdouble_matrix (n, n);
   b = JDMdouble_vector (n);
   b1 = JDMdouble_vector (n);

   for (i = 0; i < n; i++)
     {
	for (j = 0; j < n; j++)
	  a[i][j] = a1[i][j] = JDMrandom ();

	b[i] = b1[i] = JDMrandom ();
     }

   if (0 == JDMgauss_jordan (a, b, n))
     {
	for (i = 0; i < n; i++)
	  {
	     double sum = 0;
	     for (j = 0; j < n; j++)
	       sum += a[i][j] * b1[j];

	     fprintf (stdout, "%e\n", sum - b[i]);
	  }
     }

   return 0;
}
#endif
