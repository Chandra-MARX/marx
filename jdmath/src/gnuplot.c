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
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <stdarg.h>
#include "jdmath.h"

char *JDMgnuplot_command = "gnuplot";

typedef struct
{
   FILE *fp;			       /* returned from popen */
   unsigned char flags;
}
GnuPlot_Type;

static GnuPlot_Type GnuPlots[JDM_MAX_GNUPLOTS];

static GnuPlot_Type *get_gnuplot (unsigned int id, int new_flag)
{
   GnuPlot_Type *g;

   if ((id >= JDM_MAX_GNUPLOTS)
       || (JDMgnuplot_command == NULL))
     return NULL;

   g = GnuPlots + id;

   if (new_flag)
     {
	if (g->fp != NULL) return NULL;
     }
   else if (g->fp == NULL) return NULL;

   return g;
}

int JDMgnuplot_open (unsigned int id, char *plotfile)
{
   char buf[256];
   GnuPlot_Type *g;

   if (NULL == (g = get_gnuplot (id, 1)))
     return -1;

   if ((plotfile == NULL) || (*plotfile == 0))
     strcpy (buf, JDMgnuplot_command);
   else
     sprintf (buf, "%s %s", JDMgnuplot_command, plotfile);

   if (NULL == (g->fp = popen (buf, "w")))
     return -1;

   return 0;
}

int JDMgnuplot_close (unsigned int id)
{
   GnuPlot_Type *g;

   if (NULL == (g = get_gnuplot (id, 0)))
     return -1;

   if (g->fp != NULL) return -1;

   /* hmmm... should I return exit status?? */
   (void) pclose (g->fp);
   g->fp = NULL;

   return 0;
}

int JDMgnuplot_cmd (unsigned int id, char *cmd, ...)
{
   GnuPlot_Type *g;
   va_list ap;

   if (NULL == (g = get_gnuplot (id, 0)))
     return -1;

   va_start(ap, cmd);
   if (EOF == vfprintf(g->fp, cmd, ap))
     {
	va_end(ap);
	return -1;
     }
   va_end(ap);

   if ((EOF == fputc ('\n', g->fp))
       || (EOF == fflush (g->fp)))
     return -1;

   return 0;
}

FILE *JDMgnuplot_get_fp (unsigned int id)
{
   GnuPlot_Type *g;

   if (NULL == (g = get_gnuplot (id, 0)))
     return NULL;

   return g->fp;
}

#if 0
#include <signal.h>
int main ()
{
   char buf[256];

#ifdef SIGPIPE
   signal (SIGPIPE, SIG_IGN);
#endif

   if (-1 == gnuplot_open (0, NULL))
     {
	fprintf (stderr, "Unable to start 0\n");
	return -1;
     }

   while (NULL != fgets (buf, sizeof(buf), stdin))
     {
	if (-1 == gnuplot_cmd (0, buf))
	  {
	     fprintf (stderr, "Unable to write cmd\n");
	     return -1;
	  }
     }
   if (-1 == gnuplot_close (0))
     {
	fprintf (stderr, "Unable to close 0\n");
	return -1;
     }

   return 0;
}

#endif

