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
/* Gaussian elimination */
#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"

/* The Gaussian elminiation algorithm involves solving the system of equations
 * @
 * @ A_ij x_j = b_i   (i,j = 1..N)
 * @
 * by converting this system into an equivalent one where the matrix A_ij
 * is upper triangular.  Then all one needs to do is back substitution.
 * Converting this system to a triangular system is rather straightforward.
 * One can easily show that the transformations may be done in place with:
 * @
 * @    A_ij -> A_ij - A_kj * A_ik/A_kk
 * @     b_i -> b_i - b_k * A_ik/A_kk
 * @
 * where k = 1..N and (i,j) = k+1 .. N.
 *
 * During this iteration A_kk could be small and give rise to errors.  For this
 * reason, the algorithm is more complicated and the order of the equations
 * gets changed such that A_kk will be the largest of A_jk for (j>=k).  This
 * is known as partial pivoting.
 *
 * Finally, we allow b_i to be a vector, e.g., b_i = (B_l)_i; l = 1...M
 *
 * OF course, this algorithm will be formaulated in terms of 0 offset arrays.
 *
 * Profiling indicates that there is a slight performance increase by
 * unrolling the loops.  To this end, small functions like combine_vectors
 * are used.
 */

#define TOLERANCE 1e-23

int JDMgauss_elimin (double **a, double *b, unsigned int n,
		     double epsilon, int *parityp)
{
   unsigned int k, nm1;
   int parity;

   if ((n == 0) || (a == NULL) || (b == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   nm1 = n - 1;
   parity = 1;
   for (k = 0; k < nm1; k++)
     {
	double f, bk;
	double *ak;
	unsigned int pivot;
	unsigned int k1;
	unsigned int i;

	/* Find the largest element in this column and swap rows */
	f = fabs (a[k][k]);
	pivot = k;
	k1 = k + 1;
	for (i = k1; i < n; i++)
	  {
	     if (fabs (a[i][k]) > f)
	       {
		  pivot = i;
		  f = fabs(a[i][k]);
	       }
	  }
	if (pivot != k)
	  {
	     double tmp;
	     /* Since we have already eliminated the variables up to the
	      * kth one, the first k elements are 0.  So, do not include
	      * them in the swap.
	      */
	     _JDMswap_dvector (a[pivot] + k, a[k] + k, n - k);
	     tmp = b[pivot]; b[pivot] = b[k]; b[k] = tmp;
	     parity = -parity;
	  }

	if (f <= epsilon)
	  {
	     JDMmsg_error ("JDMgauss_elimin: singular matrix");
	     return -1;
	  }

	ak = a[k];
	f = (1.0 / ak[k]);
	bk = b[k];

	/* Now make the elimination step to get the modified elements. */
	for (i = k1; i < n; i++)
	  {
	     double *ai;
	     double factor;
	     /* unsigned int j; */

	     ai = a[i];
	     factor = -ai[k] * f;
	     _JDM_add_to_vector (ai + k1, ak + k1, factor, n - k1);
	     /* for (j = k1; j < n; j++) ai[j] -= ak[j] * factor; */

	     b[i] -= bk * factor;
	  }
     }

   if (fabs(a[nm1][nm1]) <= epsilon)
     {
	JDMmsg_error ("JDMgauss_elimin: singular matrix");
	return -1;
     }

   if (parityp != NULL)
     *parityp = parity;

   return 0;
}

int JDMgauss_elimin_n (double **a, unsigned int n,
		       double **b, unsigned int m,   /* b[n][m] */
		       double epsilon, int *parityp)
{
   unsigned int k, nm1;
   int parity;

   if ((n == 0) || (m == 0) || (a == NULL) || (b == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   nm1 = n - 1;
   parity = 1;
   for (k = 0; k < nm1; k++)
     {
	double f;
	double *ak, *bk;
	unsigned int pivot;
	unsigned int k1;
	unsigned int i;

	/* Find the largest element in this column and swap rows */
	f = fabs (a[k][k]);
	pivot = k;
	k1 = k + 1;
	for (i = k1; i < n; i++)
	  {
	     if (fabs (a[i][k]) > f)
	       {
		  pivot = i;
		  f = fabs(a[i][k]);
	       }
	  }
	if (pivot != k)
	  {
	     /* Since we have already eliminated the variables up to the
	      * kth one, the first k elements are 0.  So, do not include
	      * them in the swap.
	      */
	     _JDMswap_dvector (a[pivot] + k, a[k] + k, n - k);
	     _JDMswap_dvector (b[pivot], b[k], m);
	     parity = -parity;
	  }

	if (f <= epsilon)
	  {
	     JDMmsg_error ("JDMgauss_elimin_n: singular matrix");
	     return -1;
	  }

	ak = a[k];
	bk = b[k];
	f = (1.0 / ak[k]);

	/* Now make the elimination step to get the modified elements. */
	for (i = k1; i < n; i++)
	  {
	     double *ai;
	     double factor;

	     ai = a[i];
	     factor = -ai[k] * f;
	     /* for (j = k1; j < n; j++) ai[j] -= ak[j] * factor; */
	     _JDM_add_to_vector (ai + k1, ak + k1, factor, n - k1);

	     /* ai = b[i]; for (j = 0; j < m; j++) ai[j] -= bk[j] * factor; */
	     _JDM_add_to_vector (b[i], bk, factor, m);
	  }
     }

   if (fabs(a[nm1][nm1]) <= epsilon)
     {
	JDMmsg_error ("JDMgauss_elimin_n: singular matrix");
	return -1;
     }

   if (parityp != NULL)
     *parityp = parity;

   return 0;
}

int JDMback_subst (double **a, double *b, unsigned int n)
{
   unsigned int k;

   if ((n == 0) || (a == NULL) || (b == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   k = n;
   do
     {
	double *ak;
	unsigned int k1;

	k1 = k;
	k--;
	ak = a[k];
	b[k] = (b[k] - _JDM_innerprod (ak + k1, b + k1, n - k1))/ak[k];
     }
   while (k != 0);

   return 0;
}

int JDMback_subst_n (double **a, unsigned int n,
			    double **b, unsigned int m)
{
   unsigned int k;
   if ((n == 0) || (a == NULL) || (b == NULL))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   k = n;
   do
     {
	unsigned int l, k1;
	double *ak;
	double akk;

        k1 = k;
	k--;
	ak = a[k];
	akk = ak[k];

	for (l = 0; l < m; l++)
	  b[k][l] = (b[k][l] - _JDM_innerprod_col (ak+k1, b+k1, l, n-k1))/akk;
     }
   while (k != 0);

   return 0;
}

int JDMmatrix_inverse (double **a, unsigned int n)
{
   double **inv;
   unsigned int i;
   unsigned int nbytes;

   if (a == NULL)
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return -1;
     }

   inv = JDMdouble_matrix (n, n);
   if (inv == NULL)
     return -1;

   nbytes = n * sizeof (double);
   for (i = 0; i < n; i++)
     {
	memset ((char *) inv[i], 0, nbytes);
	inv[i][i] = 1;
     }

   if ((-1 == JDMgauss_elimin_n (a, n, inv, n, TOLERANCE, NULL))
       || (-1 == JDMback_subst_n (a, n, inv, n)))
     {
	JDMfree_double_matrix (inv, n);
	return -1;
     }

   for (i = 0; i < n; i++)
     memcpy ((char *)a[i], (char *)inv[i], nbytes);

   JDMfree_double_matrix (inv, n);
   return 0;
}

#ifdef TEST_GAUSSELIM

static int solve (double **a, double *b, unsigned int n, int method)
{
   double *bb;
   unsigned int i;

   if (method == 0)
     {
	if (-1 == JDMgauss_elimin (a, b, n, 1e-23, NULL))
	  return -1;

	return JDMback_subst (a, b, n);
     }

   if (-1 == JDMmatrix_inverse (a, n))
     return -1;

   bb = JDMdouble_vector (n);
   memcpy ((char *) bb, (char *)b, n * sizeof (double));

   for (i = 0; i < n; i++)
     {
	unsigned int j;
	double sum = 0;
	for (j = 0; j < n; j++)
	  sum += a[i][j] * bb[j];
	b[i] = sum;
     }
   JDMfree_double_vector (bb);
   return 0;
}

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

   if (0 == solve (a, b, n, 1))
     {
	for (i = 0; i < n; i++)
	  {
	     double sum = 0;
	     for (j = 0; j < n; j++)
	       sum += a1[i][j] * b[j];

	     fprintf (stdout, "%e\n", sum - b1[i]);
	  }
     }

   return 0;
}
#endif

