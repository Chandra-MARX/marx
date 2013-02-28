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
/* This routine implements the simplex algorithm for the minimization of
 * functions.  It is based on my understanding of the simplex concept as
 * described in Numerical Recipes; however, it is not derived from the code
 * there.
 */
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "jdmath.h"

static double Simplex_Reflect = 1.5;
static double Simplex_Contract = 0.7;
static double Simplex_Stretch = 2.0;
static int Max_Iterations = 2500;
static double Simplex_Precision = 1.0e-7;
static void (*Progress_Function) (int, double, double, unsigned int, unsigned int);
static int Simplex_Report_Freq = 50;
static int (*Validity_Function) (double *, unsigned int);

#define MAX_CHISQR 1e36

void JDMset_simplex_parameters (JDMSimplex_Type *p, int set)
{
   if (set == 0)
     {
	if (p == NULL) return;

	p->reflect = Simplex_Reflect;
	p->contract = Simplex_Contract;
	p->stretch = Simplex_Stretch;
	p->max_iterations = Max_Iterations;
	p->precision = Simplex_Precision;
	p->report_fun = Progress_Function;
	p->report_freq = Simplex_Report_Freq;
	p->check_validity_fun = Validity_Function;
	return;
     }

   if (p != NULL)
     {
	if (p->reflect > 1.0)
	  Simplex_Reflect = p->reflect;
	if ((p->contract > 0.0) && (p->contract < 1.0))
	  Simplex_Contract = p->contract;
	if (p->stretch > 1.0)
	  Simplex_Stretch = p->stretch;

	Max_Iterations = p->max_iterations;
	Simplex_Precision = p->precision;
	Progress_Function = p->report_fun;
	Simplex_Report_Freq = p->report_freq;
	Validity_Function = p->check_validity_fun;
     }
   else
     {
	Simplex_Reflect = 1.5;
	Simplex_Contract = 0.7;
	Simplex_Stretch = 2.0;
	Max_Iterations = 2500;
	Simplex_Precision = 1.0e-11;
	Progress_Function = NULL;
	Simplex_Report_Freq = 50;
	Validity_Function = NULL;
     }
}

static void find_highest_point (double *y_vector, unsigned int npts,
				unsigned int *ihi, double *hi,
				unsigned int *ilo, double *lo,
				unsigned int *inext_hi, double *next_hi)
{
   unsigned int i;
   double highest, lowest, next_highest;
   unsigned int i_highest, i_lowest, i_next_highest;
   double y;

   highest = y_vector[0];
   lowest  = y_vector[1];
   if (highest < lowest)
     {
	i_highest = 1;
	i_lowest = 0;
        next_highest = lowest;
	lowest = highest;
	highest = next_highest;
     }
   else
     {
	i_highest = 0;
	i_lowest = 1;
	next_highest = highest;
     }

   i_next_highest = i_highest;

   for (i = 1; i < npts; i++)
     {
	y = y_vector[i];
	if (y < lowest)
	  {
	     i_lowest = i;
	     lowest = y;
	  }
	else if (y > highest)
	  {
	     i_next_highest = i_highest;
	     next_highest = highest;
	     i_highest = i;
	     highest = y;
	  }
	else if (y >= next_highest)
	  {
	     i_next_highest = i;
	     next_highest = y;
	  }
     }
   *ihi = i_highest; *hi = highest;
   *ilo = i_lowest; *lo = lowest;
   *inext_hi = i_next_highest; *next_hi = next_highest;
}

/* We do not want to call this routine two much since it does not work well
 * with the CPU cache.  If possible, the COM is adjusted inline.  It ignores
 * the contribution from the point specified by 'hi'.
 */
static void compute_com (double **pnts,
			 unsigned int npts,
			 unsigned int ndim, double *com, unsigned int hi,
			 unsigned int lo,
			 double *y)
{
   unsigned int i, j;
   double sum;
   double dnpts;

   (void) lo; (void) y;

   dnpts = npts - 1.0;

   for (i = 0; i < ndim; i++)
     {
	sum = 0.0;
	for (j = 0; j < npts; j++)
	  {
	     /* if (j != hi) */
	     sum += pnts[j][i];
	  }
	com[i] = (sum - pnts[hi][i]) / dnpts;
     }
}

/* stretch nth point about the COM by lambda. */
static void stretch (double **pnts, unsigned int npts, unsigned int ndim,
		     unsigned int nth,
		     double lambda, double *com)
{
   double *xn, *xnmax, oldx, c;

   (void) npts;

   xn = pnts[nth];
   xnmax = xn + ndim;

   while (xn < xnmax)
     {
	oldx = *xn;
	c = *com++;
	*xn++ = c  + lambda * (oldx - c);
     }
}

/* contract all points toward the nth one by a factor of lambda.
 * It does NOT compute the COM.
 */
static void contract (double **pnts, unsigned int npts, unsigned int ndim,
		      unsigned int nth,
		      double lambda)
{
   double *xn, *xnmax, oldx;
   double *xnth_ptr, *xnth_save_ptr, xnth;
   unsigned int n;

   xnth_save_ptr = pnts[nth];

   for (n = 0; n < npts; n++)
     {
	if (n == nth) continue;

	xn = pnts[n];
	xnmax = xn + ndim;
	xnth_ptr = xnth_save_ptr;

	while (xn < xnmax)
	  {
	     oldx = *xn;
	     xnth = *xnth_ptr++;

	     *xn++ = xnth + lambda * (oldx - xnth);
	  }
     }
}

static int simplex (double **pnts, unsigned int npts, unsigned int ndim,
		    double (*f)(double *, unsigned int, JDM_User_Type *),
		    double *lowest, unsigned int *lowest_i, JDM_User_Type *user)
{
   unsigned int i_hi, i_lo, i_next_hi, i;
   double hi, lo, next_hi;
   int max_iter;
   double *com, *y_vector;
   double y;

   /* allocate space for the COM array */
   com = JDMdouble_vector (ndim);
   if (com == NULL)
     {
	return -1;
     }

   /* Now the array for function evaluations */
   y_vector = JDMdouble_vector (npts);
   if (y_vector == NULL)
     {
	JDMfree_double_vector (com);
	return -1;
     }

   for (i = 0; i < npts; i++) y_vector[i] = (*f)(pnts[i], ndim, user);

   max_iter = 0;
   while (max_iter++ < Max_Iterations)
     {
	double hi_lo_diff;

	find_highest_point (y_vector, npts,
			    &i_hi, &hi, &i_lo, &lo, &i_next_hi, &next_hi);

	if ((Progress_Function != NULL)
	    && ((max_iter % Simplex_Report_Freq) == 0))
	  {
	     (*Progress_Function) (max_iter, lo, hi, i_lo, i_hi);
	  }

	if (JDMUser_Break) break;

	hi_lo_diff = hi - lo;

	if (hi_lo_diff <= (fabs(hi) + fabs(lo)) * Simplex_Precision)
	  break;

	compute_com (pnts, npts, ndim, com, i_hi, i_lo, y_vector);

	/* Try to reflect highest about COM */
	stretch (pnts, npts, ndim, i_hi, -Simplex_Reflect, com);

	y = (*f)(pnts[i_hi], ndim, user);

	if (y <= lo)
	  {
	     hi = y_vector[i_hi] = y;
	     stretch (pnts, npts, ndim, i_hi, Simplex_Stretch, com);
	     y = (*f)(pnts[i_hi], ndim, user);
	     if (y >= hi)
	       {
		  /* extra stretch failed so put it back */
		  stretch (pnts, npts, ndim, i_hi, 1.0 / Simplex_Stretch, com);
	       }
	     else y_vector[i_hi] = y;
	     continue;
	  }
	else if (y >= next_hi)
	  {
	     stretch (pnts, npts, ndim, i_hi, Simplex_Contract, com);
	     y = (*f)(pnts[i_hi], ndim, user);

             if (y < next_hi)	       /* was: y < next_hi */
	       {
		  y_vector[i_hi] = y;
		  continue;
	       }

	     stretch (pnts, npts, ndim, i_hi, -1.0 / Simplex_Reflect, com);

	     y = (*f)(pnts[i_hi], ndim, user);

	     if (y < hi)
	       {
		  y_vector[i_hi] = y;
		  continue;
	       }

	     stretch (pnts, npts, ndim, i_hi, 1.0 / Simplex_Contract, com);
	     contract (pnts, npts, ndim, i_lo, Simplex_Contract);

#if 0
	     if (y < hi)
	       {
		  y_vector[i_hi] = y;
		  continue;
	       }

	     /* put it back and simply contract */
	     stretch (pnts, npts, ndim, i_hi,
		      -1.0 / (Simplex_Contract * Simplex_Reflect), com);

	     contract (pnts, npts, ndim, i_lo, Simplex_Contract, com);
#endif
	     /* Now adjust y_vector */
	     for (i = 0; i < npts; i++)
	       {
		  if (i != i_lo)
		    y_vector[i] = (*f)(pnts[i], ndim, user);
	       }
	  }
	else y_vector[i_hi] = y;
     }

   /* Re-evaluate in-case last operation screwed up lo */
   for (i = 0; i < npts; i++) y_vector[i] = (*f)(pnts[i], ndim, user);
   find_highest_point (y_vector, npts,
		       &i_hi, &hi, &i_lo, &lo, &i_next_hi, &next_hi);

   JDMfree_double_vector (com);
   JDMfree_double_vector (y_vector);
   *lowest = lo;
   *lowest_i = i_lo;
   return 0;
}

typedef struct
{
   JDMFit_Funct_Type f;
   double *xp, *yp, *dyp;
   unsigned int ndata_points;
   double *parm_list;
   unsigned int *vary_list;
   unsigned int num_parms;
   JDM_User_Type *luser;
}
Simplex_Chi_Square_Type;

double (*JDMSimplex_Chisqr)
(
 double *,			       /* x array */
 double *,			       /* y array */
 double *,			       /* dy array */
 unsigned int,			       /* num x points */
 double *,			       /* parm array */
 unsigned int,			       /* num parms */
 JDM_User_Type *,		       /* user defined */
 JDMFit_Funct_Type		       /* function to compute fit */
 );

static double compute_chisqr (double *parms, unsigned int nparms, JDM_User_Type *luser)
{
   Simplex_Chi_Square_Type *user;
   JDMFit_Funct_Type f;
   unsigned int i, imax, *v, *vmax;
   unsigned int num_parms;
   double *xp, *yp, *dyp, *p, *u, *parm_list;
   double ans, chisqr;

   (void) nparms;

   user = (Simplex_Chi_Square_Type *) luser;
   f = user->f;
   imax = user->ndata_points;
   xp = user->xp;
   yp = user->yp;
   dyp = user->dyp;
   luser = user->luser;
   parm_list = user->parm_list;
   num_parms = user->num_parms;

   /* Now, fixup the parm list to the order that the function expects. */

   v = user->vary_list;
   vmax = v + user->num_parms;
   p = parms;
   u = parm_list;
   while (v < vmax)
     {
	if (*v++) *u = *p++;
	u++;
     }

   if (Validity_Function != NULL)
     {
	if (0 == (*Validity_Function) (parm_list, num_parms))
	  return MAX_CHISQR;
     }

   if (JDMSimplex_Chisqr != NULL)
     return (*JDMSimplex_Chisqr) (xp, yp, dyp, imax, parm_list, num_parms, luser, f);

   chisqr = 0.0;
#if 1
   for (i = 0; i < imax; i++)
     {
	double sig = dyp[i];
	if (sig == 0) continue;
	ans = (*f) (i, xp[i], parm_list, num_parms, luser) - yp[i];
	ans = ans / sig;
	chisqr = chisqr + ans * ans;
     }
#else
   for (i = 0; i < imax; i++)
     {
	double sig = dyp[i];
	if (sig == 0) continue;
	ans = (*f) (i, xp[i], parm_list, num_parms, luser) - yp[i];
	ans = ans / sig;
	chisqr = chisqr + fabs (ans);
     }
#endif
#if 0
   if (chisqr > 1000)
     {
	fprintf (stderr, "Chisqr is huge: %f\n", chisqr);
     }
#endif
   return chisqr;
}

static void init_simplex_point (double *vector, double *parms,
				unsigned int *vary_list,
				unsigned int nparms, double scale,
				unsigned int vparm)
{
   unsigned int i;
   unsigned int n;
   unsigned int i_to_vary = 0;
   double p;

   n = 0;
   /* Make a copy of the parameter list but mark one we want to vary */
   for (i = 0; i < nparms; i++)
     {
	vector[i] = parms[i];
	if (vary_list[i])
	  {
	     if (vparm == n) i_to_vary = i;
	     n++;
	  }
     }

   p = parms[i_to_vary];
   do
     {
	vector[i_to_vary] = p + scale;
	scale = -0.9 * scale;
     }
   while ((scale != 0.0)
	  && (Validity_Function != NULL)
	  && (0 == (*Validity_Function)(vector, nparms)));

   n = 0;
   for (i = 0; i < nparms; i++)
     {
	if (vary_list[i])
	  {
	     if (i != n) vector[n] = vector[i];
	     n++;
	  }
     }
}

int JDMsimplex_best_fit (double *xp, double *yp, double *dyp, unsigned int npts,
			 double *parm_list, unsigned int n_parms,
			 unsigned int *vary_list, double scale,
			 double *chisqr, JDMFit_Funct_Type f,
			 JDM_User_Type *user)
{
   unsigned int vparms;			       /* actual number to vary */
   unsigned int *v, *vmax;
   Simplex_Chi_Square_Type scst;
   double **simplex_pnt_matrix;
   int return_value;
   unsigned int best_i, i;
   unsigned int matrix_dim0, matrix_dim1;

   v = vary_list;
   vmax = v + n_parms;
   vparms = 0;
   while (v < vmax) if (*v++) vparms++;

   scst.f = f;
   scst.xp = xp;
   scst.yp = yp;
   scst.dyp = dyp;
   scst.ndata_points = npts;
   scst.parm_list = parm_list;
   scst.vary_list = vary_list;
   scst.num_parms = n_parms;
   scst.luser = user;

   matrix_dim0 = vparms + 1;
   matrix_dim1 = n_parms;	       /* Use n_parms instead of vparms because
					* init_simplex_point will use the extra
					* space.
					*/
   simplex_pnt_matrix = JDMdouble_matrix (matrix_dim0, matrix_dim1);
   if (simplex_pnt_matrix == NULL) return -1;

   /* Now initialize the simplex.  The vary_list together with the scale
    * parameter determines the simplex.
    */
   v = vary_list;
   vmax = v + n_parms;
   vparms = 0;
   init_simplex_point (simplex_pnt_matrix[0],
		       parm_list, vary_list, n_parms, 0.0, 0);
   while (v < vmax)
     {
	if (*v++)
	  {
	     init_simplex_point (simplex_pnt_matrix[vparms + 1],
				 parm_list, vary_list, n_parms, scale, vparms);
	     vparms++;

	  }
     }

   return_value =  simplex (simplex_pnt_matrix, matrix_dim0, matrix_dim1,
			    compute_chisqr,
			    chisqr, &best_i,
			    (JDM_User_Type *) (void *)&scst);

   /* The next call is to simply make sure that the caller gets
    * the best value.
    */
   *chisqr = compute_chisqr (simplex_pnt_matrix[best_i],
			     vparms, (JDM_User_Type *) (void *)&scst);

   vparms = 0;
   for (i = 0; i < n_parms; i++)
     {
	if (vary_list[i])
	  {
	     parm_list[i] = simplex_pnt_matrix[best_i][vparms];
	     vparms++;
	  }
     }

   JDMfree_double_matrix (simplex_pnt_matrix, matrix_dim0);
   return return_value;
}

