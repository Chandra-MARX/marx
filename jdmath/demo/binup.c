#include <stdio.h>
#include <jdmath.h>

#include <string.h>
#include <stdlib.h>
#include <unistd.h>

static void usage (char *pgm)
{
   char *p;
   p = pgm + strlen (pgm);
   while (p > pgm)
     {
	if (*p == '/')
	  {
	     p++;
	     break;
	  }
	p--;
     }
   
   fprintf (stderr, "%s Usage:\n", p);
   fprintf (stderr, "\
%s [--help] [-z] [-i] [-c COL-NUM] [-min MIN] [-max MAX] [-numbins NUMBINS]\n\
%s [--help] [-z] [-i] [-c COL-NUM] [-min MIN] [-max MAX] [-dx DX]\n",
	    p, p);
   
   exit (1);
}

#define MAX_COLUMNS 10

int main (int argc, char **argv)
{
   double fmin, fmax, dx;
   unsigned int i;
   unsigned int num_bins;
   char *pgm;
   int column;
   double user_min, user_max;
   int has_min, has_max, has_dx;
   float *data[MAX_COLUMNS], *array;
   int cindex[MAX_COLUMNS], nread;
   unsigned int *bin;
   int zero_flag = 0;
   float *a0, *a1, *amax;
   
   pgm = argv[0];
   
   has_min = 0;
   has_max = 0;
   has_dx = 0;
   num_bins = 1024;
   column = 1;
   dx = 1.0;

   for (i = 0; i < MAX_COLUMNS; i++) cindex[i] = 0;

   i = 1;
   while (i < (unsigned int) argc)
     {
	char *arg;
   
	arg = argv[i++];
	if (!strcmp (arg, "-z")) zero_flag = 1;
	else if (i < (unsigned int) argc)
	  {
	     if (!strcmp (arg, "-c"))
	       {
		  if ((1 != sscanf (argv[i], "%d", &column))
		      || (column <= 0)
		      || (column > MAX_COLUMNS))
		    usage (pgm);
		  i++;
	       }
	     else if (!strcmp (arg, "-numbins"))
	       {
		  if (1 != sscanf (argv[i], "%u", &num_bins))
		    usage (pgm);
		  i++;
	       }
	     else if (!strcmp (arg, "-min"))
	       {
		  if (1 != sscanf (argv[i], "%lf", &user_min))
		    usage (pgm);
		  has_min = 1;
		  i++;
	       }
	     else if (!strcmp (arg, "-max"))
	       {
		  if (1 != sscanf (argv[i], "%lf", &user_max))
		    usage (pgm);
		  has_max = 1;
		  i++;
	       }
	     else if (!strcmp (arg, "-dx"))
	       {
		  if ((1 != sscanf (argv[i], "%lf", &dx))
		      || (dx <= 0.0))
		    usage (pgm);

		  i++;
	       }
	     else usage (pgm);
	  }
	else usage (pgm);
     }

   if (isatty (fileno(stdin)))
     usage (pgm);
   
   cindex [column - 1] = 1;
   if (-1 == (nread = JDMread_column_fdata (NULL, data, cindex, column)))
     {
	JDMmsg_error (pgm);
	return 1;
     }
   
   array = data[column - 1];

   a1 = a0 = array;
   amax = a0 + nread;
   if (has_min) 
     {
	fmin = user_min;
	while (a1 < amax)
	  {
	     float aval;
	     
	     if ((aval = *a1) >= fmin)
	       *a0++ = aval;
	     a1++;
	  }
	nread = a0 - array;
     }
   else if (nread > 0)
     {
	fmin = *a0++;
	while (a0 < amax)
	  {
	     if (*a0 < fmin) fmin = *a0;
	     a0++;
	  }
     }
   
   a1 = a0 = array;
   amax = a0 + nread;

   if (has_max) 
     {
	fmax = user_max;
	while (a1 < amax)
	  {
	     float aval;
	     
	     if ((aval = *a1) < fmax)
	       *a0++ = aval;
	     a1++;
	  }
	nread = a0 - array;
     }
   else if (nread > 0)
     {
	fmax = *a0++;
	while (a0 < amax)
	  {
	     if (*a0 > fmax) fmax = *a0;
	     a0++;
	  }
     }
   
   if (nread == 0)
     {
	fprintf (stderr, "No data in range\n");
	return 1;
     }

   num_bins = (unsigned int) ((fmax - fmin) / dx);
   if (num_bins == 0) num_bins = 1;

   if (NULL == (bin = (unsigned int *) JDMinteger_vector (num_bins)))
     {
	JDMmsg_error (pgm);
	return 1;
     }
   
   JDMhistogram_f (array, nread, bin, num_bins, fmin, fmax);
   
   i = 0;
   while (i < num_bins)
     {
	unsigned int b = bin[i];
	fprintf (stdout, "%e %u\n", fmin, bin[i]);
	i++;
	fmin += dx;
	
	if ((b == 0) && zero_flag)
	  {
	     while ((i < num_bins) && (bin[i] == 0))
	       {
		  i++;
		  fmin += dx;
	       }
	  }
     }
   
   
   return 0;
}

   
