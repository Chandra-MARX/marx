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
#include "jdmath.h"

/* This implementation is derived from
 * Spouge JL. Computation of the gamma, digamma, and trigamma
 *    functions. SIAM J Numerical Anal 1994;31:931-44.
 * as pointed out by Glynne Casteel <glynnec@ix.netcom.com> in the
 * sci.math.num-analysis article <6potip$l56@sjx-ixn4.ix.netcom.com>.
 */

#define NUM_COEFFS 10		       /* make this even */
static double Param = NUM_COEFFS + 0.840;
/* determined empirically using the driver in main at end of the file. */

static double Coeffs[NUM_COEFFS+1];

static int Coeffs_Initialized = 0;

static void init_coefficients (void)
{
   int i;
   double e = 2.718281828459045235360287471352;
   double f = 2.506628274631000502415765284809;   /* sqrt(2pi) */

   Coeffs[0] = exp (-Param)*f;
   Coeffs[1] = sqrt (Param-1)/e;
   for (i = 1; i < NUM_COEFFS; i++)
     {
	register double x = Param - i;

	Coeffs[i+1] = Coeffs[i]
	  * ((x - 1)*pow (1-1/x,i-0.5)/(i*e));
     }
   Coeffs_Initialized = 1;
}

double JDMlog_gamma (double x)
{
   register double sum;
   unsigned int i;

   if (Coeffs_Initialized == 0)
     init_coefficients ();

   x -= 1.0;

   sum = Coeffs[0];
   i = 1;
   while (i < NUM_COEFFS)
     {
	register double dsum;

	dsum = Coeffs[i]/(x+i);
	i++;
	sum += (dsum - Coeffs[i]/(x+i));
	i++;
     }

   return log(sum) + ((x+0.5)*log(x+Param) - x);
}

/* See A&S 6.5.4 for a definition.  Here, the series expansion given by
 * A&S 6.5.29 is used:
 *
 *    gamma_star(a,x) = exp(-x) \sum_n x^n/Gamma(a+n+1)
 *
 * Here, Gamma(a+1)=a*Gamma(a), Gamma(a+2)=(a+1)*a*Gamma(a), ....
 * Thus,
 *
 *    gamma_star(a,x) = exp(-x)/Gamma(a) \sum_n x^n c_n
 *
 * where c_n = 1/(a(a+1)...(a+n)).  For a > 0, c_n --> 0, hence based upon
 * this c_n, the sum will terminate.  However, x^n may overflow before this
 * happens.  However, as long as |x| < a, this should not happen.  So avoid
 * this function if |x| > a.
 */
static double JDMlog_gamma_star (double a, double x)
{
   double c_n, xn;
   double sum;
   int n;

   if (a == 0.0)
     return 0.0;

   if (a < 0.0)
     {
	/* FIXME!!! Handle a < 0 */
     }

   c_n = 1.0 / a;
   xn = 1.0;
   sum = c_n * xn;
   n = 0;
   while (n < 500)
     {
	n++;
	c_n /= (a + n);
	if (c_n == 0.0)
	  break;
	xn *= x;		       /* may overflow */
	sum += c_n * xn;
     }

   return log (sum) - x - JDMlog_gamma (a);
}

/* This incomplete gamma function has a continued fraction approximation
 * given by A&S 6.5.31.  This expansion may be used to compute the incomplete
 * Gamma function (A&S 6.5.1) via A&S 6.5.2,3:
 *    Gamma(a,x) = Gamma(a) - P(a,x) Gamma(a)
 * or
 *    P(a,x) = 1 - Gamma(a,x)/Gamma(a).
 *
 * The continued fraction expansion (A&S 6.5.31) is given by:
 *
 *    Gamma(a,x) = exp(-x)x^a (1/x+ (1-a)/1+ 1/x+ (2-a)/1+ 2/x+ ...)
 *
 * Use the recursion relations in theorem 2 of A&S 3.10.1 to evaluate this.
 * That is, let f_n be defined as
 *
 *    f_n = b_0 + a_1/b_1+ a_2/b_2+...a_n/b_n := A_n/B_n
 *
 * Then:
 *
 *    A_n = b_n A_{n-1} + a_n A_{n-2}
 *    B_n = b_n B_{n-1} + a_n B_{n-2}
 *
 * where A_{-1}=1, A_0=b_0, B_{-1}=0, B_0=1
 */
static double JDMlog_CapGamma (double a, double x)
{
   double a1, a0, b0, b1;
   double f;
   double renorm_factor;
   double eps = 1e-14;
   int n;

   if (x <= 0.0)
     return log(x);		       /* NaN/-Inf */

   /* Set up the recursion at the initial 1/x+ piece */
   a1 = 1;
   b1 = x;
   a0 = 0;
   b0 = 1;

   renorm_factor = 1.0/b1;
   f = a1 * renorm_factor;
   n = 1;

   if (renorm_factor != 0.0) while (n < 500)
     {
	double f1, aa;

	/* Note that the renormalization factor is 1/b1.
	 * So, replace renorm_factor*b1 combinations by 1 */

	aa = n - a;
	/* Now the (n-a)/(1+??) piece */
	a0 = (a1 + aa * a0) * renorm_factor;
	b0 = (b1 + aa * b0) * renorm_factor;

	/* Since a0,b0 now have the renorm_factor applied, omit it below */
	/* Handle {1,2,3...}/(x+??) piece */
	a1 = x * a0 + n * a1 * renorm_factor;
	/* b1 = x * b0 + n * b1 * renorm_factor;*/
	b1 = x * b0 + n;

	n++;

	if (b1 == 0.0)
	  continue;

	renorm_factor = 1.0 / b1;
	f1 = a1 * renorm_factor;

	if (fabs (f - f1) < eps * fabs(f1))
	  {
	     f = f1;
	     break;
	  }
	f = f1;
     }

   return a * log(x) - x + log (f);
}

/* See A&S 6.5.1 */
double JDMincomplete_gamma (double a, double x)
{
   if (a <= 0.0)
     return log(a);

   if (x > a)
     return 1.0 - exp (JDMlog_CapGamma (a,x) - JDMlog_gamma (a));

   return exp (a * log(x) + JDMlog_gamma_star (a, x));
}

#if 0
#include <stdlib.h>
int main (int argc, char **argv)
{
   double a = atof(argv[1]);
   double x = atof(argv[2]);

   /* fprintf (stdout, "%e\n", JDMincomplete_gamma (a, x)); */
   fprintf (stdout, "%e\n", JDMincomplete_gamma (a/2, x/2));
   return 0;
}
#endif

#if 0
static double doit (int x)
{
   double y, z;
   double f;
   int i;

   f = 1;
   for (i = 1; i <= x; i++)
     f = f * i;

   y = JDMlog_gamma (x+1);
   z = exp (y)/f-1;
   return z;
}

int main (int argc, char **argv)
{
   int i, j;
   double d;

   for (j = 0; j < 1000; j++)
     {
	double sum, max;
#if 0
	d = 0.001 * j;
	Param = NUM_COEFFS + d;
	init_coefficients ();
#endif
	sum = 0.0;
	max = 0.0;

	for (i = 1; i < 100; i += 1)
	  {
	     double z = fabs(doit (i));
	     if (z > max)
	       max = z;
	     sum += z;
	  }
	fprintf (stdout, "%g\t%g\t%g\n", d, sum*0.01, max);
	break;
     }

   return 0;
}

#endif
