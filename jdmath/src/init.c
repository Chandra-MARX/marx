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

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"

static void half_add_one (double *x, double *y)
{
   *y = 1.0 + 0.5 * *x;
}

int JDMcheck_types (unsigned long i16, 
		    unsigned long i32, 
		    unsigned long f32, 
		    unsigned long f64)
{
   if ((2 != i16)
       || (4 != i32)
       || (4 != f32)
       || (8 != f64))
     {
	JDMmsg_error ("JDMcheck_types: configuration error.  Type sizes are not correct.");
	exit (1);
     }

   /* See machar.c for an explanation */
   if (-1 == _JDM_init_machine_constants (half_add_one))
     exit (1);

   return 0;
}
