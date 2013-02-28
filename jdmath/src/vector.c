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

JDMVector_Type JDMv_vector (double x, double y, double z)
{
   JDMVector_Type a;
   a.x = x;
   a.y = y;
   a.z = z;
   return a;
}

JDMVector_Type JDMv_cross_prod (JDMVector_Type a, JDMVector_Type b)
{
   JDMVector_Type c;
   c.z = a.x * b.y - a.y * b.x;
   c.x = a.y * b.z - a.z * b.y;
   c.y = a.z * b.x - a.x * b.z;

   return c;
}

JDMVector_Type JDMv_pcross_prod (JDMVector_Type *a, JDMVector_Type *b)
{
   JDMVector_Type c;
   c.z = a->x * b->y - a->y * b->x;
   c.x = a->y * b->z - a->z * b->y;
   c.y = a->z * b->x - a->x * b->z;

   return c;
}

JDMVector_Type JDMv_smul (double x, JDMVector_Type a)
{
   a.x *= x;
   a.y *= x;
   a.z *= x;
   return a;
}

double JDMv_dot_prod (JDMVector_Type a, JDMVector_Type b)
{
   return a.x * b.x + a.y * b.y + a.z * b.z;
}

double JDMv_pdot_prod (JDMVector_Type *a, JDMVector_Type *b)
{
   return a->x * b->x + a->y * b->y + a->z * b->z;
}

double JDMv_length (JDMVector_Type a)
{
   double x, y, z, tmp;

   x = fabs (a.x);
   y = fabs (a.y);
   z = fabs (a.z);

   if (z < x)
     {
	tmp = z; z = x; x = tmp;
     }
   if (z < y)
     {
	tmp = z; z = y; y = tmp;
     }
   if (z == 0.0) return 0.0;
   x = x / z; y = y / z;
   z = z * sqrt (1.0 + x * x + y * y);
   return z;
}

extern void JDMv_normalize (JDMVector_Type *a)
{
   double len;

   if ((len = JDMv_length (*a)) != 0.0)
     {
	a->x = a->x / len;
	a->y = a->y / len;
	a->z = a->z / len;
     }
}

JDMVector_Type JDMv_unit_vector (JDMVector_Type a)
{
   JDMVector_Type n;

   n = a;
   JDMv_normalize (&n);
   return n;
}

JDMVector_Type JDMv_ax1_bx2 (double a, JDMVector_Type x1,
			     double b, JDMVector_Type x2)
{
   JDMVector_Type c;

   c.x = a * x1.x + b * x2.x;
   c.y = a * x1.y + b * x2.y;
   c.z = a * x1.z + b * x2.z;

   return c;
}

JDMVector_Type JDMv_ax1_bx2_cx3 (double a, JDMVector_Type x1,
				 double b, JDMVector_Type x2,
				 double c, JDMVector_Type x3)
{
   JDMVector_Type d;

   d.x = a * x1.x + b * x2.x + c * x3.x;
   d.y = a * x1.y + b * x2.y + c * x3.y;
   d.z = a * x1.z + b * x2.z + c * x3.z;

   return d;
}

JDMVector_Type JDMv_pax1_bx2 (double a, JDMVector_Type *x1,
			      double b, JDMVector_Type *x2)
{
   JDMVector_Type c;

   c.x = a * x1->x + b * x2->x;
   c.y = a * x1->y + b * x2->y;
   c.z = a * x1->z + b * x2->z;

   return c;
}

JDMVector_Type JDMv_sum (JDMVector_Type x1, JDMVector_Type x2)
{
   JDMVector_Type c;

   c.x = x1.x + x2.x;
   c.y = x1.y + x2.y;
   c.z = x1.z + x2.z;

   return c;
}

JDMVector_Type JDMv_diff (JDMVector_Type x1, JDMVector_Type x2)
{
   JDMVector_Type c;

   c.x = x1.x - x2.x;
   c.y = x1.y - x2.y;
   c.z = x1.z - x2.z;

   return c;
}

double JDMv_distance (JDMVector_Type a, JDMVector_Type b)
{
   a.x -= b.x;
   a.y -= b.y;
   a.z -= b.z;

   return JDMv_length (a);
}

JDMVector_Type JDMv_rotate_vector1 (JDMVector_Type p, JDMVector_Type n,
					 double cos_theta, double sin_theta)
{
   double pn = JDMv_dot_prod (p, n);

   return JDMv_ax1_bx2_cx3 (cos_theta, p,
			    pn * (1.0 - cos_theta), n,
			    sin_theta, JDMv_pcross_prod(&n,&p));
}

JDMVector_Type JDMv_rotate_unit_vector1 (JDMVector_Type p, JDMVector_Type n,
					 double cos_theta, double sin_theta)
{
   JDMVector_Type u;

   u = JDMv_rotate_vector1 (p, n, cos_theta, sin_theta);
   JDMv_normalize (&u);
   return u;
}

JDMVector_Type JDMv_rotate_unit_vector (JDMVector_Type p, JDMVector_Type n,
					double theta)
{
   return JDMv_rotate_unit_vector1 (p, n, cos(theta), sin(theta));
}

void JDMv_unit_vector_to_spherical (JDMVector_Type p,
				    double *thetap, double *phip)
{
   double theta, phi;
   double sin_theta;

   if (fabs (p.z) >= 1.0)
     {
	if (p.z >= 1.0)
	  *thetap = 0;
	else
	  *thetap = PI;
	*phip = 0;
	return;
     }

   theta = acos (p.z);		       /* [0, PI] */
   sin_theta = sin (theta);

   if (fabs (p.x) <= fabs (p.y))
     {
	phi = acos (p.x / sin_theta);       /* [0, PI] */
	if (p.y < 0.0)
	  phi = -phi;			       /* [-PI to PI] */
     }
   else
     {
	phi = asin (p.y / sin_theta);   /* [-PI/2,PI/2] */
	if (p.x < 0)
	  {
	     if (phi >= 0) phi = PI - phi;
	     else phi = -PI - phi;
	  }
     }

   *thetap = theta;
   *phip = phi;
}

void JDMv_vector_to_spherical (JDMVector_Type p,
			       double *rp, double *thetap, double *phip)
{
   double r;

   r = JDMv_length (p);
   if (r == 0.0)
     {
	*thetap = *phip = *rp = 0;
	return;
     }

   p.x /= r;
   p.y /= r;
   p.z /= r;

   *rp = r;
   JDMv_unit_vector_to_spherical (p, thetap, phip);
}

JDMVector_Type JDMv_spherical_to_vector (double r, double theta, double phi)
{
   JDMVector_Type p;
   double s;

   s = sin(theta);
   p.x = r * s * cos (phi);
   p.y = r * s * sin (phi);
   p.z = r * cos (theta);
   return p;
}

/* Find the rotation axis and angle such that when a is rotated about that
 * axis by the angle, it rotates into b.  a and b are assumed to be
 * unit vectors.
 */
double JDMv_find_rotation_axis (JDMVector_Type a, JDMVector_Type b, JDMVector_Type *c)
{
   double theta, len;
   JDMVector_Type cc;

   theta = JDMv_dot_prod (a, b);
   if (fabs (theta) > 1.0)
     {
	if (theta > 1.0)
	  theta = 1.0;
	else
	  theta = -1.0;
     }
   theta = acos (theta);
   cc = JDMv_cross_prod (a, b);
   len = JDMv_length (cc);
   if (len > 0.0)
     {
	c->x = cc.x/len;
	c->y = cc.y/len;
	c->z = cc.z/len;
	return theta;
     }

   /* Singular case: a = const b. Just pick any vector orthogonal to a. */

   if ((a.x != 0.0) || (a.y != 0.0))
     {
	c->x = -a.y;
	c->y = -a.x;
	c->z = 0.0;
     }
   else
     {
	c->x = -a.z;
	c->y = 0.0;
	c->z = 0.0;
     }

   JDMv_normalize (c);
   return theta;
}

void
JDMv_spherical_to_triad (double theta, double phi,
			 JDMVector_Type *r_hat, JDMVector_Type *theta_hat,
			 JDMVector_Type *phi_hat)
{
   double ct, st, cp, sp;

   ct = cos (theta);
   st = sin (theta);
   cp = cos (phi);
   sp = sin (phi);

   r_hat->x = st * cp;
   r_hat->y = st * sp;
   r_hat->z = ct;

   phi_hat->x = -sp;
   phi_hat->y = cp;
   phi_hat->z = 0.0;

   theta_hat->x = ct * cp;
   theta_hat->y = ct * sp;
   theta_hat->z = -st;
}
