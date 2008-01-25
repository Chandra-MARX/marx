/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2002-2004 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research
    under contract SV1-61010 from the Smithsonian Institution.

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
#include "config.h"
#include "marx-feat.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#include <string.h>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#include <limits.h>

#ifndef PATH_MAX
# define PATH_MAX 1024
#endif

#define MAX_PATH_LEN PATH_MAX

#include "marx.h"
#include "_marx.h"

double _Marx_TStart_Yrs;

static void do_malloc_error (void) /*{{{*/
{
   marx_error ("Memory allocation failure.");
}

/*}}}*/

char *marx_malloc (unsigned int nbytes) /*{{{*/
{
   char *s;
   
   s = malloc (nbytes);
   if (s == NULL)
     do_malloc_error ();
   return s;
}

/*}}}*/

char *marx_calloc (unsigned int n, unsigned int size) /*{{{*/
{
   char *s;
   
   s = calloc (n, size);
   if (s == NULL)
     do_malloc_error ();
   return s;
}

/*}}}*/


char *marx_realloc (char *ptr, unsigned int size)
{
   if (ptr == NULL)
     return marx_malloc (size);
   
   ptr = realloc (ptr, size);
   if (ptr == NULL)
     do_malloc_error ();
   return ptr;
}

int marx_free (char *s) /*{{{*/
{
   if (s == NULL) return -1;
   free (s);
   return 0;
}

/*}}}*/

/* Concatenate dir and file to produce dir/file.  This routine works with or
 * without the trailing / on dir.  If file is "", then dir is simply 
 * returned.
 */
char *marx_dircat (char *dir, char *file) /*{{{*/
{
   unsigned int dirlen;
   unsigned int filelen;
   char *filename;
   
   if ((dir == NULL) && (file == NULL))
     return NULL;
   
   dirlen = filelen = 0;
   
   if (dir != NULL) dirlen = strlen (dir);
   if (file != NULL) filelen = strlen (file);
   
   if (NULL == (filename = malloc (dirlen + filelen + 2)))
     {
	do_malloc_error ();
	return NULL;
     }
   
   if (dirlen)
     {
	strcpy (filename, dir);
	/* Add final / if it is not already there. */
	if (filename[dirlen - 1] != '/')
	  {
	     filename[dirlen] = '/';
	     dirlen++;
	     filename[dirlen] = 0;
	  }
     }
   if (filelen) strcpy (filename + dirlen, file);
   return filename;
}

/*}}}*/

int marx_file_exists (char *file)
{
   struct stat st;
   int m;
   
#ifdef _S_IFDIR
# ifndef S_IFDIR
#  define S_IFDIR _S_IFDIR
# endif
#endif
   
#ifndef S_ISDIR
# ifdef S_IFDIR
#  define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
# else
#  define S_ISDIR(m) 0
# endif
#endif
   
   if (stat(file, &st) < 0) return 0;
   m = st.st_mode;
   
   if (S_ISDIR(m)) return 2;
   return 1;
}

char *marx_find_file_in_path (char *dirs, char *file, char colon)
{
   char dir_buf [MAX_PATH_LEN + 1];
   char *dir0, *dir1;
   char *filename;

   /* If file is absolute, try to find it. */
   if (file == NULL) return NULL;
   
   if (*file == '/')
     {
	if (1 == marx_file_exists (file))
	  return marx_dircat (NULL, file);
     }
   
   if (dirs == NULL)
     return NULL;

   dir0 = dirs;
   
   while (*dir0)
     {
	unsigned int len;
	
	dir1 = strchr (dir0, colon);
	if (dir1 == NULL)
	  {
	     dir1 = dir0 + strlen (dir0);
	  }
	if (dir0 == dir1)
	  {
	     dir0++;
	     continue;
	  }
	if ((len = (dir1 - dir0)) <= MAX_PATH_LEN)
	  {
	     strncpy (dir_buf, dir0, len);
	     dir_buf[len] = 0;
	     
	     filename = marx_dircat (dir_buf, file);
	     if (filename == NULL)
	       return NULL;
	     
	     if (1 == marx_file_exists (filename))
	       return filename;
	     
	     free (filename);
	  }
	
	dir0 = dir1;
	if (*dir0) dir0++;
     }
   
   return NULL;
}

int marx_set_time_years (double tyrs)
{
   _Marx_TStart_Yrs = tyrs;
   return 0;
}

Param_File_Type *marx_pf_parse_cmd_line 
  (char *file, char *mode, int argc, char **argv)
{
   Param_File_Type *pf;
   
   pf = pf_parse_cmd_line_no_exit (file, mode, argc, argv);
   if (pf != NULL)
     return pf;
#ifdef MARX_PFILE_DIR
   if (PF_Errno == PF_FILE_NOT_FOUND)
     {
	fprintf (stderr, "*****\n\
Unable to locate a %s parameter file.  Marx parameter files may be\n\
found in:\n\
 %s/\n\
*****\n",
		 file, MARX_PFILE_DIR);
     }
#endif
   exit (1);
}
