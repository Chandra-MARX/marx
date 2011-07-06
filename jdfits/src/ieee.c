/*
    Copyright (C) 2002 MIT Center For Space Research

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
/* Conversions to/from IEEE floating point formats. */

/* IEEE is fairly simple.  It looks like:
 *
 *   seeeeee emmmmmmm mmmmmmmm mmmmmmmm
 *
 * where 's' is the sign bit, e represents the exponent and m is the mantissa.
 * Then, the number is given by:
 *
 *    value = (-1)^s 2^(e - 127) * 1.mmmmmmmmm
 */

#include "config.h"

#include <stdio.h>

void vax_to_ieee (float *xp)
{
   float x = *xp;
   unsigned char *chp1, *chp2;

   chp1 = (unsigned char *) xp;
   chp2 = (unsigned char *) &x;

   chp1[1] = chp2[0];
   chp1[2] = chp2[3];
   chp1[3] = chp2[2];
   chp1[0] = chp2[1] + 1;

   fprintf (stderr, "%X --> %X\n", *(unsigned int *) &x, *(unsigned int *) xp);
}

void error (char *s)
{
   fprintf (stderr, s);
   /* exit (-1); */
}

int main (int argc, char **argv)
{
   float x;
   int use_long = 0;
   unsigned long n;
   unsigned char *xp, *xpmax, ch;
   unsigned short s;

   /* First of all, determine something about the data sizes */
   fprintf (stdout, "/* sizeof (short) = %d */\n", sizeof (short));
   fprintf (stdout, "/* sizeof (int) = %d */\n", sizeof (int));
   fprintf (stdout, "/* sizeof (long) = %d */\n", sizeof (long));
   fprintf (stdout, "/* sizeof (float) = %d */\n", sizeof (float));
   fprintf (stdout, "/* sizeof (double) = %d */\n", sizeof (double));

   if (sizeof (short) == 2)
     {
	fprintf (stdout, "typedef short int16;\n");
     }
   else if (sizeof (int) == 2)
     {
	fprintf (stdout, "typedef int int16;\n");
     }
   else error ("Machine cannot support 16 bit integers.");

   if (sizeof (int) == 4)
     {
	fprintf (stdout, "typedef int int32;\n");
     }
   else if (sizeof (long) == 4)
     {
	fprintf (stdout, "typedef long int32;\n");
	use_long = 1;
     }
   else error ("Machine cannot support 32 bit integers.");

   if (sizeof (float) == 4)
     {
	fprintf (stdout, "typedef float float32;\n");
     }
   else if (sizeof (double) == 4)
     {
	fprintf (stdout, "typedef double float32;\n");
     }
   else error ("Machine cannot support 32 bit floats.");

   if (sizeof (float) == 8)
     {
	fprintf (stdout, "typedef float float64;\n");
     }
   else if (sizeof (double) == 8)
     {
	fprintf (stdout, "typedef double float64;\n");
     }
   else error ("Machine cannot support 64 bit floats.");

   s = 0x1234;
   if (*(char *) &s == 0x34)
     {
	fprintf (stdout, "#define JDFITS_BYTE_SWAP\n");
     }

   /* Now determine the floating point style */

   x = 1.2345678;

   if (use_long)
     {
	n = *(unsigned long *) &x;
     }
   else n = (unsigned long) *(unsigned int *) &x;

   if (n == 0x3F9E0651)
     {
	/* IEEE --- Do nothing. */
     }
   else if (n == 0x0651409E)
     {
	fprintf (stdout, "#define JDFITS_VAX_FLOAT\n");
     }
   else
     {
        error ("Unknown floating point format.  Not supported.");
     }

   xp = (unsigned char *) &x;
   xpmax = xp + sizeof (float);

   x = 27.0004;
   while (xp < xpmax)
     {
	int i = 128;
	ch = *xp++;

	while (i)
	  {
	     if (ch & i) putc ('1', stdout);
	     else putc ('0', stdout);
	     i = i >> 1;
	  }
	putc (' ', stdout);
     }
   putc ('\n', stdout);

   return 0;
}

