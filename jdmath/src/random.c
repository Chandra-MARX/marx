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
/* Returns a random double precision number between 0.0 and 1.0 */

#include "config.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdmath.h"
#include "_jdmath.h"

/* The random number generator used here was derived from the article:
 *
 *      G. Marsaglia and A. Zaman, ``Some portable very-long-period
 *           random number generators'', Computers in Physics, Vol. 8, No. 1,
 *           Jan/Feb 1994
 *
 * However, I believe the algorithm (mzran13) published there contained a few
 * typos.  The version coded here passes all the DIEHARD battery of tests
 * as well as the 2d random walk and n-block tests by I. Vattulainen et al.,
 * Phys. Rev. Lett. 73, 2513 (1994).  I used 'nblocktest.f' and '2drwtest.f'
 * available from netlib.
 *
 * The main diffeerence between mzran13 and my version is that I have replaced
 * the 69069 generator with Marsaglia's MWC1616 generator.
 *
 * The resulting period is something like 6 * 10^17 * 4 *10^28
 * or 24 * 10^45 which should be contrasted with the 10^18 period of
 * Numerical Recipes $1000 RAN2 generator.  This generator is also twice
 * as fast as RAN2.
 *
 * Note: George Marsaglia is _the_ leading authority on the generation
 * of pseudorandom numbers.
 */

#define SEED_X	521288629U
#define SEED_Y	362436069U
#define SEED_Z	16163801U
#define SEED_C  1    /* (SEED_Y > SEED_Z) */
#define SEED_N	1131199209U

#define SEED_X1	1245364
#define SEED_Y1 46347

struct _JDMRandom_Type
{
   uint32 x, y, z, c, n;
   uint32 x1, y1;
};

static JDMRandom_Type This_Rand =
{
   SEED_X, SEED_Y, SEED_Z, SEED_C, SEED_N,
   SEED_X1, SEED_Y1
};

int JDMseed_random (JDMRandom_Type *rt, unsigned long seed)
{
   if (rt == NULL)
     rt = &This_Rand;

   rt->x = (uint32) SEED_X;
   rt->y = (uint32) SEED_Y;
   rt->z = (uint32) SEED_Z;

   rt->x += (uint32) seed;
   rt->y += (uint32) seed / 2;
   rt->z += 2 * (uint32) seed;

   rt->n = (uint32) SEED_N;
   rt->n += (uint32) (seed & 0xFFFF);

   rt->x1 = (uint32) SEED_X1 + 23343 * seed;
   rt->y1 = (uint32) SEED_Y1 - 12 * seed;

   rt->c = (rt->y > rt->z) ? 1 : 0;
   return 0;
}

uint32 JDMgenerate_uint32_random (JDMRandom_Type *rt)
{
   register uint32 s, x, y, c;

   if (rt == NULL)
     rt = &This_Rand;

   x = rt->x;
   y = rt->y;
   c = rt->c;

   if (y > x + c)
     {
	s = y - (x + c);
	rt->c = 0;
     }
   else
     {
	s = y - (x + c) - 18;
	rt->c = 1;
     }

   rt->x = y;
   rt->y = rt->z;
   rt->z = s;

   x = rt->x1;
   y = rt->y1;

   /* This part is the Marsaglia MWC1616 generator.
    * The 18000 and 30903 numbers are from Marsaglia's diehard program.
    *
    * FIXME: priority=low: Make these numbers customizable from a table
    *        that Marsaglia provides in Diehard.
    */
   rt->x1 = x = 18000U * (x & 0xFFFFU) + (x >> 16);
   rt->y1 = y = 30903U * (y & 0xFFFFU) + (y >> 16);

   return s + ((x << 16) + (y * 0xFFFFU));
}

double JDMgenerate_random (JDMRandom_Type *rt)
{
   return ((double)JDMgenerate_uint32_random (rt) * (1.0 / (double)(uint32)0xFFFFFFFFU));
}

uint32 JDMuint32_random (void)
{
   return JDMgenerate_uint32_random (&This_Rand);
}

double JDMrandom (void)
{
   return ((double)JDMgenerate_uint32_random (&This_Rand) * (1.0 / ((double)(uint32)0xFFFFFFFFU)));
}

/* Seed the random number generator */
int JDMsrandom (unsigned long s)
{
   return JDMseed_random (&This_Rand, s);
}

JDMRandom_Type *JDMcreate_random(void)
{
   JDMRandom_Type *rt;

   rt = (JDMRandom_Type *) _JDMmalloc (sizeof (JDMRandom_Type), NULL);
   if (rt == NULL)
     return NULL;

   JDMseed_random (rt, 0);
   return rt;
}

void JDMfree_random (JDMRandom_Type *rt)
{
   _JDMfree ((char *)rt);
}

/* This generator is not the best but it is fast.  For a good generator,
 * use JDMrandom.
 */
static uint32 Fast_Random;
void JDMseed_fast_random (unsigned long s)
{
   Fast_Random = (uint32) s;
}

uint32 JDMfast_uint32_random (void)
{
   return (Fast_Random = Fast_Random * 69069U + 1013904243U);
}

double JDMfast_random (void)
{
   Fast_Random = Fast_Random * 69069U + 1013904243U;
#if 1
   return ((double)Fast_Random/4294967296.0);
#else
   return ((double)Fast_Random/((double)(uint32)0xFFFFFFFFU));
#endif
}
