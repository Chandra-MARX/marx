/* -*- C -*- */

unsigned char *STR_READ_X_TO_Y_FUN (TO_TYPE *y, unsigned int n, unsigned char *s)
{
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	FROM_TYPE x;
	s = STR_READ_X_FUN (&x, 1, s);
	y[i] = (TO_TYPE) x;
     }
   return s;
}
#undef STR_READ_X_FUN
#undef STR_READ_X_TO_Y_FUN
#undef FROM_TYPE
#undef TO_TYPE
