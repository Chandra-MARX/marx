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
#include <stdlib.h>
#endif
#include <string.h>

#include "jdmath.h"
#include "_jdmath.h"

double JDM_Mach_Eps = 2.22044604925031308084726e-16;

/* This is convoluted to fool the optimizer */
int _JDM_init_machine_constants (void (*silly)(double *, double *))
{
   JDM_Mach_Eps = 1.0;

   while (1)
     {
	double ans;

	(*silly) (&JDM_Mach_Eps, &ans);

	if (ans == 1.0)
	  break;

	JDM_Mach_Eps *= 0.5;
     }

   return 0;
}

