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

int JDMath_Error;
int JDMUser_Break;

static char *map_math_error (void)
{
   char *err = NULL;
   
   switch (JDMath_Error)
     {
      case JDMATH_MALLOC_ERROR:
	err = "Memory Allocation Failure";
	break;
	
      case JDMATH_INVALID_PARAMETER:
	err = "Invalid Parameter";
	break;
	
      case JDMATH_FILE_OPEN_ERROR:
	err = "File Open Error";
	break;
	
      case JDMATH_FILE_READ_ERROR:
	err = "File Read Error";
	break;
	
      case JDMATH_FILE_WRITE_ERROR:
	err = "File Write Error";
	break;
	
      case JDMATH_FILE_CLOSE_ERROR:
	err = "File Close Error";
	break;

      case JDMATH_CORRUPT_FILE_ERROR:
	err = "File appears corrupt";
	break;
	
      case JDMATH_DIVIDE_ZERO_ERROR:
	err = "Division by zero.";
	break;
	
      default:
	if (JDMUser_Break) err = "User Break";
	else err = "Unknown Error";
     }
   return err;
}

void JDMmsg_error (char *s)
{
   char *err;

   err = map_math_error ();
   
   if (s == NULL) s = "JDMath Lib Error";
   
   fprintf (stderr, "%s: %s\n", s, err);
   fflush (stderr);
}

	  

void JDMmsg_error2 (char *s1, char *s2)
{
   char *err;
   
   if (s2 == NULL)
     {
	JDMmsg_error (s1);
	return;
     }
   
   err = map_math_error ();
   if (s1 == NULL) s1 = "JDMath Lib Error";
   
   fprintf (stderr, "%s: %s: %s\n", s1, s2, err);
   fflush (stderr);
}


#ifndef SLMALLOC
#define SLMALLOC malloc
#define SLFREE free
#define SLCALLOC calloc
#define SLREALLOC realloc
#endif

char *_JDMmalloc (unsigned long nbytes, char *err)
{
   char *p;
   
   p = SLMALLOC (nbytes);
   if (p == NULL)
     {
	JDMath_Error = JDMATH_MALLOC_ERROR;
	JDMmsg_error (err);
     }
   
   return p;
}

void _JDMfree (char *s)
{   
   if (s != NULL) SLFREE (s);
}

char *_JDMrealloc (char *s, unsigned int len)
{   
   if (s == NULL)
     return _JDMmalloc (len, "realloc");
   
   s = SLREALLOC (s, len);
   if (s == NULL)
     JDMath_Error = JDMATH_MALLOC_ERROR;
   
   return s;
}

