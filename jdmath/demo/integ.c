#include <stdio.h>
#include <math.h>

int main (int argc, char **argv)
{
   double x, y;
   char line [256];
   double dx;
   double last_x, last_y, last_start_x;
   double sum, total_sum;
   unsigned int num;

   num = 0;
   last_y = last_x = last_start_x = 0.0;
   total_sum = sum = 0.0;

   dx = 0.024;

   if ((argc == 2) && (1 != sscanf (argv[1], "%lf", &dx)))
     {
	fprintf (stderr, "Usage: %s PIXELSIZE\n", argv[0]);
	return 1;
     }

   while (NULL != fgets (line, sizeof (line), stdin))
     {
	double dsum;

	if (2 != sscanf (line, "%lf %lf", &x, &y))
	  continue;

	if (num == 0)
	  {
	     last_start_x = last_x = x;
	     last_y = y;
	  }

	num++;

	if (x >= last_start_x + dx)
	  {
	     double new_x, new_y;

	     new_x = last_start_x + dx;
	     new_y = last_y + (y - last_y) * ((new_x - last_x) / (x - last_x));

	     sum += (new_x - last_x) * 0.5 * (new_y + last_y);
	     total_sum += sum;

	     fprintf (stdout, "%e\t%e\t%e\n", last_start_x + 0.5 * dx,
		      sum / dx, total_sum);

	     sum = 0.0;
	     last_x = last_start_x = new_x;
	     last_y = new_y;
	  }

	dsum = (x - last_x) * 0.5 * (y + last_y);
	sum += dsum;

	last_x = x;
	last_y = y;
     }

   total_sum += sum;
   if (last_start_x != last_x)
     fprintf (stdout, "%e\t%e\t%e\n",
	      last_start_x + 0.5 * dx,
	      sum / (last_x - last_start_x),
	      total_sum);

   return 0;
}
