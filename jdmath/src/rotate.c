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
#include <math.h>


#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#include <string.h>

#include "jdmath.h"

/* These define ACTIVE rotations of the coordinate system.  */

void JDM3m_rot_x_matrix (JDM_3Matrix_Type m, double theta)
{
   double c, s;
   
   c = cos (theta);
   s = sin (theta);
   
   s = -s;

   m[0][0] = 1.0;  m[0][1] = 0.0;  m[0][2] = 0.0;
   m[1][0] = 0.0;  m[1][1] = c;    m[1][2] = s;
   m[2][0] = 0.0;  m[2][1] = -s;   m[2][2] = c;
}

void JDM3m_rot_y_matrix (JDM_3Matrix_Type m, double theta)
{
   double c, s;
   
   c = cos (theta);
   s = sin (theta);
   
   s = -s;

   m[0][0] = c;    m[0][1] = 0.0;  m[0][2] = -s;
   m[1][0] = 0.0;  m[1][1] = 1.0;  m[1][2] = 0.0;
   m[2][0] = s;    m[2][1] = 0.0;  m[2][2] = c;
}

void JDM3m_rot_z_matrix (JDM_3Matrix_Type m, double theta)
{
   double c, s;
   
   c = cos (theta);
   s = sin (theta);
   
   s = -s;

   m[0][0] = c;    m[0][1] = s;    m[0][2] = 0.0;
   m[1][0] = -s;   m[1][1] = c;    m[1][2] = 0.0;
   m[2][0] = 0.0;  m[2][1] = 0.0;  m[2][2] = 1.0;
}

/* This function creates a rotation matrix appropriate for a rotation
 * about an arbitrary direction.  It is easy to show that the matrix R_ij
 * is given by:
 *    R_ij = cos(t) d_ij + n_i n_j (1 - cos(t))  sin t e_ijk n_k
 * where d_ij is the Kronecker delta and e_ijk is the completely 
 * anti-symmetric tensor.  Also n_i are the components of a unit vector.
 */
void JDM3m_rot_matrix (JDM_3Matrix_Type m, JDMVector_Type n, double theta)
{
   double c, s, n_0, n_1, n_2, n_01, n_02, n_12, c1;
   
   c = cos (theta);
   s = sin (theta);
   
   c1 = 1 - c;

   n_0 = n.x;
   n_1 = n.y;
   n_2 = n.z;
   n_01 = n_0 * n_1 * c1;
   n_02 = n_0 * n_2 * c1;
   n_12 = n_1 * n_2 * c1;
   
   m[0][0] = n_0 * n_0 * c1 + c;
   m[0][1] = n_01 - s * n_2;
   m[0][2] = n_02 + s * n_1;

   m[1][0] = n_01 + s * n_2;
   m[1][1] = n_1 * n_1 * c1 + c;
   m[1][2] = n_12 - s * n_0;

   m[2][0] = n_02 - s * n_1;
   m[2][1] = n_12 + s * n_0;
   m[2][2] = n_2 * n_2 * c1 + c;
}


void JDM3m_mul (JDM_3Matrix_Type c, JDM_3Matrix_Type a, JDM_3Matrix_Type b)
{
   unsigned int i, j;
   
   for (i = 0; i < 3; i++)
     {
	for (j = 0; j < 3; j++)
	  c[i][j] = a[i][0] * b[0][j]  + a[i][1] * b[1][j] + a[i][2] * b[2][j];
     }
}

JDMVector_Type JDM3m_vector_mul (JDM_3Matrix_Type a, JDMVector_Type v)
{
   JDMVector_Type b;
   
   b.x = a[0][0] * v.x + a[0][1] * v.y + a[0][2] * v.z;
   b.y = a[1][0] * v.x + a[1][1] * v.y + a[1][2] * v.z;
   b.z = a[2][0] * v.x + a[2][1] * v.y + a[2][2] * v.z;
   
   return b;
}

   
