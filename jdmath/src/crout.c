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
/* Crout LU factorization */
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"

/*
 * Crout's algorithm with partial pivoting.
 *
 * This algorithm performs a Lower-Upper factorization of a matrix A.  The
 * method is actually quite simple.  Write:
 *
 *     L_ij = L_ij Q(i>=j)
 *     U_ij = U_ij Q(i<=j)
 *
 * where Q(x) is 1,0 if x is TRUE,FALSE.  Then,
 *
 *     A_ij = \sum_k L_ik U_kj Q(i>=k) Q(k<=j)
 * or
 *
 *     A_ij = \sum_{k<=min(i,j)} L_ik U_kj
 *
 * or
 *
 *     A_ij = L_ii U_ij + \sum_{k<i} L_ik U_kj     (i <= j)
 *          = L_ij U_jj + \sum_{k<j} L_ik U_kj     (j <= i)
 *
 * We are free to choose L_ii = 1 for all i (since there are n^2+n unknown
 * matrix elements.  Then for i <= j, we obtain:
 *
 *     U_ij = A_ij - \sum{k<i} L_ik U_kj    (i <= j)
 *
 * and for j > i:
 *
 *     L_ij = (1/U_jj) (A_ij - \sum{k<j} L_ik U_kj)  (j >= i)
 *
 * These equations allow a straightforward recursive solution for U_ij and L_ij
 * since the sums on the RHS contain L_ik and U_kj that are known.  In fact,
 * by not storing the diagonal elements of L, which have been set to 1, we can
 * do the factorization in place since each A_ij in the original matrix is
 * explicitly one time.
 *
 * The only problem is that because of the presence of (1/U_jj), this algorithm
 * is unstable without pivoting.  Pivoting can be easily accomplished by noting
 * that the computation of U_ij involves the element A_ij, the values of
 * L_ij from the same row, and the values of U_kj from earlier rows.  The key
 * point is that L_ij involves values of L from the same row.  Thus, we may
 * perform simple row swapping.
 *
 * Note that the algorithm equilibrates the matrix in order to find the
 * appropriate pivot.
 */

#define TOLERANCE 1e-23

int JDM_lu_decomp (double **a, unsigned int n, unsigned int *pivot, double epsilon, int *parityp)
{
   unsigned int j;
   int det = 1;
   double *eq;

   if ((n == 0) || (a == NULL) || (pivot == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   if (NULL == (eq = _JDM_equilibrate (a, n, NULL)))
     return -1;

   for (j = 0; j < n; j++)
     {
	unsigned int i, ipiv;
	double big;

	/* Do i < j case */
	for (i = 0; i < j; i++)
	  {
	     double *ai = a[i];
	     ai[j] -= _JDM_innerprod_col (ai, a, j, i);
	  }

	/* Now i >= j.
	 * We also pivot in this step.  Since we do not know what value
	 * goes on the diagonal, compute all possible diagonal values, then
	 * pivot, then divide.
	 */

	big = 0;
	ipiv = j;

	for (i = j; i < n; i++)
	  {
	     double *ai = a[i];
	     ai[j] -= _JDM_innerprod_col (ai, a, j, j);

	     if (eq[j] * fabs (ai[j]) > big)
	       {
		  big = eq[j] * fabs (ai[j]);
		  ipiv = i;
	       }
	  }
	if (ipiv != j)
	  {
	     double eq_j;
	     _JDMswap_dvector (a[ipiv], a[j], n);
	     eq_j = eq[j]; eq[j] = eq[ipiv]; eq[ipiv] = eq_j;
	     det = -det;
	  }
	pivot[j] = ipiv;

	big = a[j][j];
	if (fabs (big) < epsilon)
	  {
	     JDMfree_double_vector (eq);
	     return -1;
	  }

	big = 1.0/big;
	for (i = j + 1; i < n; i++)
	  a[i][j] *= big;
     }

   if (parityp != NULL)
     *parityp = det;

   JDMfree_double_vector (eq);
   return 0;
}

static void backsubst (double **a, unsigned int n, double *b)
{
   unsigned int i;

   i = n;
   do
     {
	double *ai;
	unsigned int i1;

	i1 = i;
	i--;
	ai = a[i];
	b[i] = (b[i] - _JDM_innerprod (ai + i1, b + i1, n - i1))/ai[i];
     }
   while (i != 0);
}

static void forwsubst (double **a, unsigned int n, double *b, unsigned int *pivot)
{
   unsigned int i;
   /* Do forward substitution and unscramble b in the process */
   for (i = 0; i < n; i++)
     {
	double x_i;
	unsigned int j;

	j = pivot[i];
	x_i = b[j] - _JDM_innerprod (a[i], b, i);
	if (i != j)
	  b[j] = b[i];
	b[i] = x_i;
     }
}

int
JDM_lu_backsubst (double **a, unsigned int n, unsigned int *pivot, double *b)
{
   if ((a == NULL) || (b == NULL) || (pivot == NULL) || (n == 0))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   forwsubst (a, n, b, pivot);
   backsubst (a, n, b);

   return 0;
}

static void forwsubst_unit_vector (double **a, unsigned int n,
				   double *b, unsigned int *pivot,
				   unsigned int col)
{
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	if (i > col)
	  {
	     double *bcol = b + col;
	     while (i < n)
	       {
		  b[i] = -_JDM_innerprod (a[i] + col, bcol, i - col);
		  i++;
	       }
	     break;
	  }

	if (col == i)
	  col = pivot[i];
	else if (col == pivot[i])
	  col = i;

	b[i] = (i == col) ? 1.0 : 0.0;
     }
}

int JDM_ludecomp_inverse (double **a, unsigned int n)
{
   unsigned int *pivot;
   double *col;
   double **inv_a;
   unsigned int i, j;
   int ret;

   if (n == 0)
     return 0;

   if (a == NULL)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   pivot = NULL;
   col = NULL;
   inv_a = NULL;
   ret = -1;

   if ((NULL == (pivot = (unsigned int *)JDMinteger_vector (n)))
       || (NULL == (col = JDMdouble_vector (n)))
       || (NULL == (inv_a = JDMdouble_matrix (n, n))))
     goto return_error;

   if (-1 == JDM_lu_decomp (a, n, pivot, TOLERANCE, NULL))
     goto return_error;

   for (j = 0; j < n; j++)
     {
	forwsubst_unit_vector (a, n, col, pivot, j);
	backsubst (a, n, col);

	for (i = 0; i < n; i++)
	  inv_a [i][j] = col[i];
     }

   for (i = 0; i < n; i++)
     {
	double *ai, *inv_ai;

	ai = a[i];
	inv_ai = inv_a[i];
	for (j = 0; j < n; j++)
	  ai[j] = inv_ai[j];
     }

   ret = 0;
   /* drop */

   return_error:
   JDMfree_integer_vector ((int *) pivot);
   JDMfree_double_vector (col);
   JDMfree_double_matrix (inv_a, n);

   return ret;
}

#ifdef TEST_LUDECOMP
#include <stdio.h>
#include "jdmath.h"

int main (int argc, char **argv)
{
   unsigned int n;
   unsigned int i, j;
   double **a, **a1;
   unsigned int count;

   n = 100;

   a = JDMdouble_matrix (n, n);
   a1 = JDMdouble_matrix (n, n);

   count = 5;

   while (count)
     {
	double max;

	count--;

	for (i = 0; i < n; i++)
	  {
	     for (j = 0; j < n; j++)
	       a[i][j] = a1[i][j] = JDMrandom ();
	  }

	if (-1 == JDM_ludecomp_inverse (a, n))
	  continue;

	max = 0;

	for (i = 0; i < n; i++)
	  {
	     unsigned int k;
	     for (k = 0; k < n; k++)
	       {
		  double sum = 0;
		  for (j = 0; j < n; j++)
		    sum += a[i][j] * a1[j][k];

		  if (i == k)
		    sum -= 1.0;

		  if (fabs (sum) > max) max = fabs(sum);
	       }
	  }
	fprintf (stdout, "%e\n", max);
     }

   return 0;
}

#endif

