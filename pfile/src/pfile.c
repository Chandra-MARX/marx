/*
    This file is part of the MIT PFILE Parameter Library

    Copyright (C) 2002 Massachusetts Institute of Technology

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

#include <stdio.h>
#include <string.h>


#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif


#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifdef HAVE_LIMITS_H
# include <limits.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef PATH_MAX
# define MAX_PATH_LEN PATH_MAX
#else
# define MAX_PATH_LEN 1024
#endif

#include <stdarg.h>

#ifndef S_ISREG
# define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#ifndef S_ISLNK
# define S_ISLNK(m)  (((m) & S_IFMT) == S_IFLNK)
#endif

#ifndef R_OK
# define R_OK    4/* test for read permission */
# define W_OK    2/* test for write permission */
# define X_OK    1/* test for execute (search) permission */
# define F_OK    0/* test for presence of file */
#endif

#include "pfile.h"
#include "_pfile.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

/* Concatenate dir and file to produce dir/file.  This routine works with or
 * without the trailing / on dir.  If file is "", then a copy of dir is
 * returned.  If file starts with /, a copy of file is returned.
 */
static char *dircat (char *dir, char *file)
{
   unsigned int dirlen;
   unsigned int filelen;
   char *filename;
   
   if (*file == '/') return _pf_create_string (file);
   
   dirlen = strlen (dir);
   filelen = strlen (file);
   
   if (NULL == (filename = _pf_malloc (dirlen + filelen + 2)))
     return NULL;
   
   if (dirlen)
     {
	strcpy (filename, dir);
	/* Add final / if it is not already there. */
	if (filelen && (filename[dirlen - 1] != '/'))
	  {
	     filename[dirlen] = '/';
	     dirlen++;
	  }
     }
   strcpy (filename + dirlen, file);
   return filename;
}

   
/* The first element of this list represents the current working dir.
 */
typedef struct _Directory_List_Type
{
   struct _Directory_List_Type *next;
   char *path;
   
   unsigned int flags;
#define PF_SYSTEM_PATH		0x1
#define PF_USER_PATH		0x2
#define PF_VERBOSE_PATH		0x8000
}
Directory_List_Type;

static Directory_List_Type *Dir_List;

static void append_dir_type_to_list (Directory_List_Type *d)
{
   Directory_List_Type *dlast;

   d->next = NULL;
   if (Dir_List == NULL)
     {
	Dir_List = d;
	return;
     }
   
   dlast = Dir_List;
   while (dlast->next != NULL)
     dlast = dlast->next;

   dlast->next = d;
}

static void prepend_dir_type_to_list (Directory_List_Type *d)
{
   d->next = Dir_List;
   Dir_List = d;
}

static Directory_List_Type *create_dir_list_type (char *dir, unsigned int dirlen, unsigned int flags)
{
   Directory_List_Type *d;

   if ((dir == NULL)
       || (dirlen == 0))
     return NULL;

   d = (Directory_List_Type *) _pf_malloc (sizeof (Directory_List_Type));
   if (d == NULL) return NULL;
   
   if (NULL == (d->path = (char *) _pf_malloc (dirlen + 2)))
     {
	SLFREE (d);
	return NULL;
     }

   strncpy (d->path, dir, dirlen);
   if (dir[dirlen - 1] != '/') 
     {
	d->path[dirlen] = '/';
	dirlen++;
     }
   d->path[dirlen] = 0;
   
   if ((flags & PF_SYSTEM_PATH) == 0)
     {
	/* make sure this directory is writable.  */
	if (dirlen)
	  d->path[dirlen - 1] = 0;
	
	if (-1 == access (d->path, W_OK))
	  {
	     flags &= ~PF_USER_PATH;
	     flags |= PF_SYSTEM_PATH;
	  }

	if (dirlen)
	  d->path[dirlen - 1] = '/';
     }
   
   d->flags = flags;
   return d;
}

static int add_dir_to_list (char *dir, unsigned int dirlen, 
			    unsigned int flags, int append)
{
   Directory_List_Type *d;

   if (dir == NULL)
     return 0;			       /* be nice */

   if (dirlen == 0)
     return 0;

   if (NULL == (d = create_dir_list_type (dir, dirlen, flags)))
     return -1;

   if (append) append_dir_type_to_list (d);
   else prepend_dir_type_to_list (d);

   return 0;
}


/* Ok, here env has the following format (example):
 * dir1 dir2:dir3 dir4; dir5 dir6
 * where dir1-dir4 are user directories and dir5 and dir6 are system.
 */
static int add_to_search_path (char *env, int append)
{
   unsigned int flags = PF_USER_PATH;

   if (env == NULL)
     return 0;

   while (1)
     {
	char *env1, ch;

	env = _pf_skip_whitespace (env);
	if (*env == 0) return 0;
	
	/* This routine NEVER returns NULL */
	env1 = _pf_strbrk (env, " ;:");
	
	if (env1 != env)
	  {
	     if (-1 == add_dir_to_list (env, (unsigned int) (env1 - env), flags, append))
	       return -1;
	  }
	     
	ch = *env1;
	if (ch != 0)
	  {
	     env1++;
	     if (ch == ';') flags = PF_SYSTEM_PATH;
	  }
	env = env1;
     }
}

static int initialize_dirlist (void)
{
   static int initialized = 0;

   if (-1 == add_to_search_path (getenv ("PFILES"), 1))
     return -1;
   if (-1 == add_to_search_path (getenv ("UPARM"), 1))
     return -1;
   /* And finally the current directory.   It is VERY important that this
    * come last because some programs do not deal very well with the parm
    * file in the current directory.
    */
   if (-1 == add_to_search_path ("./", 1))
     return -1;

   initialized = 1;
   return 0;
}

int pf_add_to_search_path (char *path, int append)
{
   return add_to_search_path (path, append);
}


static char *extract_dir (char *filename)
{
   char *file;
   
   file = _pf_rstrchr (filename, '/');
   if (file == NULL)
     {
        return _pf_create_string ("");
     }
   file++;
   return _pf_create_nstring (filename, (unsigned int) (file - filename));
}

static char *follow_link (char *filename)
{
   char filebuf[MAX_PATH_LEN + 1];
   char *dir, *file;
   int len;
   
   if (filename == NULL) return NULL;
   len = readlink (filename, filebuf, MAX_PATH_LEN);
   
   if ((len <= 0)
       || (NULL == (dir = extract_dir (filename))))
     {
	SLFREE (filename);
	return NULL;
     }
   
   SLFREE (filename);
   file = _pf_create_nstring (filebuf, (unsigned int) len);
   if (file == NULL)
     {
	SLFREE (dir);
	return NULL;
     }
   
   filename = dircat (dir, file);
   SLFREE (dir);
   SLFREE (file);
   
   return filename;
}

/* This function finds a regular file named filename.  If not found, the
 * argument filename is freed and NULL is returned.  This means that the
 * argument must be MALLOCED!!!  If filename corresponds to a link, the
 * file pointed to by the link is returned.
 */
static char *find_regular_file (char *filename, int access_mode)
{
   struct stat st;
   unsigned max_links = 10;
   
   if (filename == NULL) return NULL;

   if ((-1 == stat (filename, &st))
       || (0 == S_ISREG(st.st_mode))
       || (-1 == access (filename, access_mode)))
     {
	PF_Errno = PF_ACCESS_ERROR;
	SLFREE (filename);
	return NULL;
     }
   
   /* Check to see whether or not this is a link. */
   while (1)
     {
	if ((max_links == 0)
	    || (-1 == lstat (filename, &st)))
	  {
	     PF_Errno = PF_ACCESS_ERROR;
	     SLFREE (filename);
	     return NULL;
	  }
	
	if (0 == S_ISLNK(st.st_mode)) break;
	  
	filename = follow_link (filename);
	if (filename == NULL) return NULL;
	max_links--;
     }
   
   return filename;
}

/* If *file is 0, just return valid directory according to access_mode. */
static char *find_a_parameter_file (char *file, unsigned int *mask,
				    int access_mode)
{
   Directory_List_Type *d;
   char *filename;
   
   /* Do not perform a search if the path is absolute or explicitly
    * refers to the current directory.
    */
   if ((*file == '/')
       || ((file[0] == '.') && (file[1] == '/')))
     {
	if (NULL == (filename = dircat ("", file)))
	  return NULL;
	
	return find_regular_file (filename, access_mode);
     }
   
   d = Dir_List;
   while (d != NULL)
     {
	if (*mask & d->flags)
	  {
	     if (NULL == (filename = dircat (d->path, file)))
	       return NULL;
	     
	     if (*mask & PF_VERBOSE_PATH)
	       fprintf (stderr, "Trying %s\n", filename);
	     
	     if (*file == 0)
	       {
		  if (-1 == access (filename, access_mode))
		    {
		       SLFREE (filename);
		       return NULL;
		    }
		  return filename;
	       }
	     
	     filename = find_regular_file (filename, access_mode);
	     if (filename != NULL)
	       {
		  *mask = d->flags;
		  return filename;
	       }
	  }
	d = d->next;
     }
   
   PF_Errno = PF_FILE_NOT_FOUND;
   return NULL;
}

static Param_File_Type *Pfile_List;

static int add_parmfile_to_list (Param_File_Type *p)
{
   Param_File_Type *plast;
   
   if (Pfile_List == NULL)
     {
	Pfile_List = p;
	return 0;
     }
   
   plast = Pfile_List;
   while (plast->next != NULL) plast = plast->next;
   plast->next = p;
   
   return 0;
}

static int remove_paramfile_from_list (Param_File_Type *p)
{
   Param_File_Type *plast;
   
   plast = Pfile_List;
   if (plast == p)
     {
	Pfile_List = plast->next;
	return 0;
     }
   
   while (plast != NULL)
     {
	if (plast->next == p)
	  {
	     plast->next = p->next;
	     return 0;
	  }
	plast = plast->next;
     }
   PF_Errno = PF_BAD_ARGUMENT;
   return -1;
}

static Param_File_Type *have_parameter_file (char *filename)
{
   Param_File_Type *p = Pfile_List;
   
   while (p != NULL)
     {
	if (!strcmp (p->input_filename, filename))
	  break;
	
	p = p->next;
     }
   return p;
}


static char *find_parameter_file (char *file, unsigned int *flags, int access_mode)
{
   char *filename = NULL;
   char *ext_list, *file_ext, *ext, *save_ext_list;
   unsigned int nth;
   
   ext_list = ".par ";	       /* .rdb files may be problematic
				* since I have no idea what they are.
				*/
   if (NULL != (save_ext_list = getenv ("PFEXTN")))
     {
	save_ext_list = _pf_strcat (ext_list, save_ext_list);
	if (save_ext_list == NULL)
	  return NULL;
	ext_list = save_ext_list;
     }
	
   nth = 0;
   while (NULL != (ext = _pf_extract_string_element (ext_list, " ", &nth)))
     {
	if (*ext != 0)
	  {
	     file_ext = _pf_strcat (file, ext);
	     SLFREE (ext);
	     if (file_ext == NULL) break;
	     
	     filename = find_a_parameter_file (file_ext, flags, access_mode);
	     SLFREE (file_ext);
	     if (filename != NULL) break;
	  }
	else SLFREE (ext);
     }
	
   if (save_ext_list != NULL) SLFREE (save_ext_list);
   
   if (filename == NULL)
     {
	filename = find_a_parameter_file (file, flags, access_mode);
     }
   return filename;
}


Param_File_Type *pf_open_parameter_file (char *file, char *openmode)
{
   Param_File_Type *p;
   char *filename;
   FILE *fp;
   unsigned int dir_flags, default_dir_flags = 0;
   int access_mode;
   
   if (file == NULL)
     {
	pf_error ("pf_open_parameter_file: NULL file passed.");
	return NULL;
     }
   
   if (-1 == initialize_dirlist ())
     return NULL;
   
   if (openmode == NULL) openmode = "rw";
   
   access_mode = 0;
   if (NULL != _pf_strchr (openmode, 'r')) access_mode |= R_OK;
   if (NULL != _pf_strchr (openmode, 'w')) access_mode |= W_OK;
   if (NULL != _pf_strchr (openmode, 'x')) access_mode |= X_OK;
   if (NULL != _pf_strchr (openmode, 'f')) access_mode |= F_OK;

   /* Verbose reporting */
   if (NULL != _pf_strchr (openmode, 'V')) default_dir_flags |= PF_VERBOSE_PATH;
     
   /* Try to find a user version.  If that fails, find a system copy and then
    * and create filename to point to a user one.  If mode has R, then look
    * for a system one.  Never write to a system copy.
    */
   
   filename = NULL;
   if (NULL == _pf_strchr (openmode, 'R'))
     {
	dir_flags = PF_USER_PATH | default_dir_flags;
	filename = find_parameter_file (file, &dir_flags, access_mode);
     }
   
   if (filename == NULL) 
     {
	dir_flags = PF_SYSTEM_PATH | default_dir_flags;
	filename = find_parameter_file (file, &dir_flags, R_OK);
	if (filename == NULL) return NULL;
     }
   
   if (NULL == (p = have_parameter_file (filename)))
     {
	if (NULL == (fp = fopen (filename, "r")))
	  {
	     pf_error ("Unable to open parameter file %s.", filename);
	     PF_Errno = PF_FILE_OPEN_ERROR;
	     SLFREE (filename);
	     return NULL;
	  }

	p = _pf_read_parm_file (filename, fp);
	fclose (fp);
	
	if (p != NULL)
	  {
	     if (-1 == add_parmfile_to_list (p))
	       {
		  _pf_free_param_file (p);
		  p = NULL;
	       }
	     else
	       {
		  char modestr[PF_MAX_LINE_LEN];
		  /* get the default mode for the file */
		  p->mode = PF_QUERY_MODE | PF_LEARN_MODE;
		  
		  
		  if ((pf_parameter_exists (p, "mode"))
		      && (0 == pf_get_string (p, "mode", modestr, sizeof(modestr))))
		    {
		       if (-1 == _pf_parse_mode (modestr, &p->mode))
			 goto close_and_return_error;
		    }
		  
		  if (NULL != _pf_strchr (openmode, 'Q'))
		    p->mode |= PF_NEVER_QUERY_MODE;

		  if (NULL != _pf_strchr (openmode, 'L'))
		    p->mode |= PF_FORCE_LEARN_MODE;
	       }
	  }
     }
   
   if (p != NULL)
     {
	p->num_references += 1;
	/* If it is system file, open a user version and switch name to
	 * it.  Note that this involves looking for a directory with
	 * RWX permission bits set.
	 */
	
	if (NULL != _pf_strchr (openmode, 'w'))
	  {
	     p->flags |= PFILE_WRITE_MODE;
	  }

	if ((dir_flags & PF_SYSTEM_PATH)
	    && (access_mode & W_OK))
	  {
	     char *new_filename;
	     char *f;
	     
	     dir_flags = PF_USER_PATH;
	     new_filename = find_a_parameter_file ("", &dir_flags, 
						   R_OK | W_OK | X_OK);
	     if (new_filename == NULL)
	       goto close_and_return_error;
	     
	     /* The new_filename is really a directory.  Add to it the
	      * file that is specified in filename.
	      */
	     f = _pf_rstrchr (filename, '/');
	     if (f == NULL) f = filename; else f++;
	     f = dircat (new_filename, f);
	     SLFREE (new_filename);
	     if (f == NULL)
	       goto close_and_return_error;
	     
	     SLFREE (filename);
	     filename = NULL;
	     
	     if (p->output_filename != NULL) SLFREE (p->output_filename);
	     p->output_filename = f;
	     
	     p->flags |= PFILE_DIRTY;
	  }
     }
   
   if (filename != NULL) SLFREE (filename);
   return p;
   
   close_and_return_error:
   if (p != NULL)
     {
	remove_paramfile_from_list (p);
	_pf_free_param_file (p);
     }
   if (filename != NULL)
     SLFREE (filename);
   
   return NULL;
}

static int write_quoted (char *str, char q, FILE *fp)
{
   if (str == NULL) return 0;
   if (q) putc (q, fp);
   str = _pf_escape_string (str);
   if (str == NULL) return -1;
   fputs (str, fp);
   SLFREE (str);
   if (q) putc (q, fp);
   return 0;
}


static int save_pfile (Param_File_Type *p)
{
   Param_Type *pf;
   FILE *fp;
   char *filename;
   
   filename = p->output_filename;
   if (filename == NULL)
     filename = p->input_filename;
   
   if (NULL == (fp = fopen (filename, "w")))
     {
	pf_error ("Unable to open %s for writing.", filename);
	PF_Errno = PF_FILE_OPEN_ERROR;
	return -1;
     }

   pf = p->pf;
   while (pf != NULL)
     {
	char buf[20], *b;
	char quote_char;

	fputs (pf->name, fp);
	
	
	if (pf->type == PF_COMMENT_TYPE)
	  {
	     putc ('\n', fp);
	     pf = pf->next;
	     continue;
	  }
	
	b = buf;
	*b++ = ',';
	if (pf->type & PF_LIST_TYPE)
	  {
	    *b++ = '*';
	  }
	
	quote_char = 0;
	switch (pf->type & 0xFF)
	  {
	   case PF_BOOLEAN_TYPE:
	     *b++ = 'b';
	     break;
	     
	   case PF_INTEGER_TYPE:
	     *b++ = 'i';
	     break;
	     
	   case PF_REAL_TYPE:
	     *b++ = 'r';
	     break;
	     
	   case PF_DOUBLE_TYPE:
	     *b++ = 'd';
	     break;
	     
	   case PF_STRING_TYPE:
	     *b++ = 's';
	     quote_char = '"';
	     break;
	     
	   case PF_FILE_TYPE:
	     *b++ = 'f';
	     if (pf->type & PF_FILE_EXISTS) *b++ = 'e';
	     if (pf->type & PF_FILE_NEXISTS) *b++ = 'n';
	     if (pf->type & PF_FILE_READABLE) *b++ = 'r';
	     if (pf->type & PF_FILE_WRITABLE) *b++ = 'w';
	     quote_char = '"';
	     break;
	     
	   default:
	     pf_error ("Type %c not supported.", pf->type & 0xFF);
	     PF_Errno = PF_NOT_IMPLEMENTED;
	     fclose (fp);
	     return -1;
	  }
	*b++ = ',';
	if (pf->mode & PF_AUTO_MODE) *b++ = 'a';
	if (pf->mode & PF_QUERY_MODE) *b++ = 'q';
	if (pf->mode & PF_HIDDEN_MODE) *b++ = 'h';
	if (pf->mode & PF_LEARN_MODE) *b++ = 'l';
	*b++ = ',';
	
	*b = 0;
	fputs (buf, fp);
	
	if (-1 == write_quoted (pf->value, quote_char, fp))
	  goto error_return;
	putc (',', fp);
	
	if (-1 == write_quoted (pf->min, quote_char, fp))
	  goto error_return;
	putc (',', fp);

	if (-1 == write_quoted (pf->max, quote_char, fp))
	  goto error_return;
	putc (',', fp);

	if (-1 == write_quoted (pf->prompt, '"', fp))
	  goto error_return;

	putc ('\n', fp);
	
	pf = pf->next;
     }
   
   fclose (fp);
   return 0;
	
   error_return:
   fclose (fp);
   return -1;
}

int pf_close_parameter_file (Param_File_Type *p)
{
   if ((p == NULL)
       || (p->num_references <= 0))
     {
	PF_Errno = PF_BAD_ARGUMENT;
	return -1;
     }
   
   if (p->num_references > 1)
     {
	p->num_references -= 1;
	return 0;
     }
   
   if ((((p->flags & PFILE_DIRTY)
	 && (p->flags & PFILE_WRITE_MODE))
	|| (p->flags & PFILE_WRITE_ON_CLOSE)) 
       && (-1 == save_pfile (p)))
     {
	pf_error ("Error saving parameter file %s.", 
		  (p->output_filename == NULL ? p->input_filename : p->output_filename));
	return -1;
     }
   
   if (-1 != remove_paramfile_from_list (p))
     _pf_free_param_file (p);
   
   return 0;
}

int pf_set_output_filename (Param_File_Type *p, char *f)
{
   if (p == NULL)
     return -1;
   
   f = _pf_create_string (f);
   if (f == NULL)
     return -1;
   
   if (p->output_filename != NULL) SLFREE (p->output_filename);
   p->output_filename = f;

   /* Force a write to the new name */
   p->flags |= PFILE_WRITE_ON_CLOSE;

   return 0;
}

char *pf_get_input_filename (Param_File_Type *p)
{
   if (p == NULL)
     return NULL;
   
   return _pf_create_string (p->input_filename);
}

char *pf_get_output_filename (Param_File_Type *p)
{
   if (p == NULL)
     return NULL;

   if (p->output_filename == NULL)
     return _pf_create_string (p->input_filename);
   else
     return _pf_create_string (p->output_filename);
}
