#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

static char *do_realloc (char *s, unsigned int len)
{
   if (s == NULL)
     s = malloc (len);
   else
     s = realloc (s, len);
   
   if (s == NULL)
     {
	fprintf (stderr, "Out of memory.\n");
	exit (1);
     }
   return s;
}

static char *str_malloc (char *s)
{
   unsigned int len;
   char *s1;
   
   len = strlen (s) + 1;
   s1 = do_realloc (NULL, len);
   strcpy (s1, s);
   return s1;
}



typedef struct 
{
   char *line;
   double value;
}
Sort_Type;

static int cmp_fun (const void *a, const void *b)
{
   double xa, xb;
   
   xa = ((Sort_Type *) a)->value;
   xb = ((Sort_Type *) b)->value;
   
   if (xa > xb) return 1;
   if (xa == xb) return 0;
   return -1;
}   
   
static void sort_them (Sort_Type *lines, unsigned int num_lines)
{
   qsort (lines, num_lines, sizeof (Sort_Type), cmp_fun);
}

static void usage (void)
{
   fprintf (stderr, "Usage: nsort < in > out\n");
   exit (1);
}

int main (int argc, char **argv)
{
   char line [4096];
   unsigned int max_lines;
   unsigned int num_lines, i;
   FILE *fp;
   Sort_Type *lines;
   
   if (isatty (fileno(stdin)))
     {
	usage ();
     }
   
   if (argc != 1)
     usage ();
   
   fp = stdin;
   
   max_lines = num_lines = 0;
   lines = NULL;
   
   while (NULL != fgets (line, sizeof (line), fp))
     {
	double value;
	
	if (1 != sscanf (line, "%lf", &value))
	  {
	     fputs (line, stdout);
	     continue;
	  }
	
	if (max_lines == num_lines)
	  {
	     max_lines += 1024;
	     
	     lines = (Sort_Type *) do_realloc ((char *) lines, max_lines * sizeof (Sort_Type));
	  }
	
	lines [num_lines].line = str_malloc (line);
	lines [num_lines].value = value;
	num_lines++;
     }
   
   sort_them (lines, num_lines);
   
   for (i = 0; i < num_lines; i++)
     {
	fputs (lines[i].line, stdout);
	free (lines[i].line);
     }
   free (lines);
   return 0;
}

   
