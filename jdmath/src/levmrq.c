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

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"

typedef struct
{
   unsigned int num_vary;
   unsigned int *parm_index;
   double *tmp_parms;
   double *delta_parms;
   double lambda;
   double last_chisqr;
}
JDMLevMarq_Type;

static void free_lm (JDMLevMarq_Type *lm)
{
   if (lm == NULL) return;
   _JDMfree ((char *) lm->parm_index);
   _JDMfree ((char *) lm->tmp_parms);
   _JDMfree ((char *) lm->delta_parms);
   _JDMfree ((char *) lm);
}

JDMLevMarq_Type *JDMleven_marq_open (unsigned int *vary_list,
				     unsigned int num_parms,
				     unsigned int num_to_vary)
{
   JDMLevMarq_Type *lm;
   unsigned int i, j;

   if ((num_parms == 0) || (num_to_vary == 0)
       || (num_to_vary > num_parms))
     {
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	return NULL;
     }

   if (NULL == (lm = (JDMLevMarq_Type *) _JDMmalloc (sizeof (JDMLevMarq_Type),
						     "JDMleven_marq_open")))
     return NULL;
   memset ((char *) lm, 0, sizeof (JDMLevMarq_Type));

   if ((NULL == (lm->parm_index = (unsigned int *) JDMinteger_vector (num_to_vary)))
       || (NULL == (lm->tmp_parms = JDMdouble_vector (num_to_vary)))
       || (NULL == (lm->delta_parms = JDMdouble_vector (num_to_vary))))
     {
	free_lm (lm);
	return NULL;
     }

   lm->num_vary = num_to_vary;

   j = 0;
   for (i = 0; i < num_parms; i++)
     {
	if (vary_list [i] == 0)
	  continue;

	lm->parm_index[j] = i;
	j++;
     }

   lambda = 0.001;
   return lm;
}

int JDMleven_marq (JDMLevMarq_Type *lm,
		   double *x, double *y, double *dy, unsigned int npts,
		   double *parm_list, double **covar, double **alpha,
