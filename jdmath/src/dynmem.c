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

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"


double *JDMdouble_vector (unsigned int n)
{
   return (double *) _JDMmalloc (n * sizeof (double), "JDMdouble_vector");
}

float *JDMfloat_vector (unsigned int n)
{
   return (float *) _JDMmalloc (n * sizeof (float), "JDMfloat_vector");
}

int *JDMinteger_vector (unsigned int n)
{
   return (int *) _JDMmalloc (n * sizeof (int), "JDMinteger_vector");
}

void JDMfree_integer_vector (int *p)
{
   _JDMfree ((char *) p);
}

void JDMfree_double_vector (double *v)
{
   _JDMfree ((char *) v);
}

void JDMfree_float_vector (float *v)
{
   _JDMfree ((char *) v);
}

void JDMfree_double_matrix (double **matrix, unsigned int n)
{  
   unsigned int i;
   if (matrix == NULL) return;
   
   for (i = 0; i < n; i++)
     {
	_JDMfree ((char *) matrix[i]);
     }
   _JDMfree ((char *) matrix);
}


double **JDMdouble_matrix (unsigned int n, unsigned int m)
{   
   double **matrix;
   unsigned int i;
   
   matrix = (double **) _JDMmalloc (n * sizeof (double **), "JDMdouble_matrix");
   if (matrix == NULL)
     return NULL;
   
   /* initialize everything to NULL */
   for (i = 0; i < n; i++) matrix[i] = NULL;

   for (i = 0; i < n; i++)
     {
	double *ptr;
	ptr = (double *) _JDMmalloc (m * sizeof (double), "JDMdouble_matrix");
	if (ptr == NULL) 
	  {
	     JDMfree_double_matrix (matrix, n);
	     return NULL;
	  }
	matrix[i] = ptr;
     }
   
   return matrix;
}

void JDMfree_float_matrix (float **matrix, unsigned int n)
{  
   unsigned int i;
   if (matrix == NULL) return;
   
   for (i = 0; i < n; i++)
     {
	_JDMfree ((char *) matrix[i]);
     }
   _JDMfree ((char *) matrix);
}

float **JDMfloat_matrix (unsigned int n, unsigned int m)
{   
   float **matrix;
   unsigned int i;
   
   matrix = (float **) _JDMmalloc (n * sizeof (float **), "JDMfloat_matrix");
   if (matrix == NULL)
     return NULL;
   
   /* initialize everything to NULL */
   for (i = 0; i < n; i++) matrix[i] = NULL;

   for (i = 0; i < n; i++)
     {
	float *ptr;
	ptr = (float *) _JDMmalloc (m * sizeof (float), "JDMfloat_matrix");
	if (ptr == NULL) 
	  {
	     JDMfree_float_matrix (matrix, n);
	     return NULL;
	  }
	matrix[i] = ptr;
     }
   
   return matrix;
}

   
   
   

