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
#include <stdlib.h>
#endif

#include <stdarg.h>

#include "pfile.h"
#include "_pfile.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

int PF_Errno;

void pf_error (char *msg, ...)
{
   va_list ap;
   
   va_start(ap, msg);
   (void) vfprintf(stderr, msg, ap);
   va_end(ap);
   putc ('\n', stderr);
   if (PF_Errno == 0) PF_Errno = PF_UNKNOWN_ERROR;
}


/* This routine strips a string s of surrounding quotes (if present).  It 
 * performs no backslash processing.
 * 
 * Note: There is a weakness in the definition of fields that refer to other
 * parameters (indirect types).  The weakness is that there appears to be no
 * way to form a string that starts with the two characters ')*' and not have 
 * it interpreted as an indirect type.  The weakness is in the definition as
 * expressed in the man pages for parameter files and NOT in my implementation.
 */
static char *parse_string_field (char *s, char **val, char *field)
{
   char quote;
   char ch;
   
   quote = *s;
   
   if ((quote == '\'') || (quote == '"'))
     {
	/* Skip past quote */
	s++;
	*val = s;
     }
   else 
     {
	quote = ',';
	*val = s;
     }
   
   
   while (((ch = *s) != 0) && (ch != quote))
     {
	if (ch == '\\') 
	  {
	     s++;
	     if (*s == 0)
	       {
		  pf_error ("Extra \\ character in %s field.", field);
		  PF_Errno = PF_CORRUPT_FIELD;
		  return NULL;
	       }
	  }
	s++;
     }
   
   
   if (quote == ',')
     {
	if (s == *val) *val = NULL;
	else *s = 0;
	if (ch != 0) s++;
	return s;
     }
   
   if (ch == 0)
     {
	pf_error ("Unable to find closing quote for %s field.", field);
	PF_Errno = PF_CORRUPT_FIELD;
	return NULL;
     }
   *s = 0;
   s++;		       /* skip quote. */
   /* skip over whitespace too.  Error if we do not see end of field */
   s = _pf_skip_whitespace (s);
   if (*s != 0)
     {
	if (*s != ',')
	  {
	     pf_error ("Excess garbage at end of %s field.", field);
	     PF_Errno = PF_CORRUPT_FIELD;
	     return NULL;
	  }
	s++;		       /* skip , */
     }
   return s;
}
   
   


static int chop_line (char *line, char **name, char **type, char **mode, 
		      char **value, char **min, char **max, char **prompt)
{
   char *s, ch;
   
   *name = line;
   *type = *mode = *min = *max = *value = *prompt = NULL;
   
   if (*line == '#')
     return 0;
   
   s = line;
   
   /* Find type */
   while (((ch = *s) != 0) && (ch != ','))
     s++;

   if (ch == 0)
     {
	pf_error ("%s\nUnable to find TYPE information.", line);
	PF_Errno = PF_CORRUPT_FIELD;
	return -1;
     }
   *s++ = 0;
   *type = s;
  
   /* Find mode */
   while (((ch = *s) != 0) && (ch != ','))
     s++;
   if (ch == 0)
     {
	pf_error ("%s\nUnable to find MODE information.", line);
	PF_Errno = PF_CORRUPT_FIELD;
	return -1;
     }
   *s++ = 0;
   *mode = s;
   
   while (((ch = *s) != 0) && (ch != ','))
     s++;
   if (ch == 0) return 0;	       /* no more parameters */

   *s++ = 0;			       /* terminate mode */
   
   /* Rest of the fields are optional.  We must be careful here because of 
    * the possibilty of quote characters and embedded commas. 
    */
   
   /* Value field is next. */
   s = parse_string_field (s, value, "VALUE");
   if (s == NULL) return -1;
   s = parse_string_field (s, min, "MIN");
   if (s == NULL) return -1;
   s = parse_string_field (s, max, "MAX");
   if (s == NULL) return -1;
   s = parse_string_field (s, prompt, "PROMPT");
   if (s == NULL) return -1;
   if (*s != 0)
     {
	pf_error ("Extra junk after prompt field for parameter '%s'", *name);
	return -1;
     }
   
   return 0;
}

void _pf_free_current_value (Param_Type *pf)
{
   if ((((pf->type & 0xFF) == PF_STRING_TYPE)
	|| ((pf->type & 0xFF) == PF_FILE_TYPE))
       && (pf->flags & PF_CURRENT_VALUE)
       && (pf->current_value.sval != NULL))
     {
	free (pf->current_value.sval);
	pf->current_value.sval = NULL;
     }
   pf->flags &= ~(PF_CURRENT_VALUE);
}

/* This routine free all memory associated with pf.  */
static void free_param_type (Param_Type *pf)
{
   if (pf->name != NULL) SLFREE (pf->name);
   if (pf->prompt != NULL) SLFREE (pf->prompt);
   if (pf->value != NULL) SLFREE (pf->value);
   if (pf->min != NULL) SLFREE (pf->min);
   if (pf->max != NULL) SLFREE (pf->max);
   
   if ((((pf->type & 0xFF) == PF_STRING_TYPE)
	|| ((pf->type & 0xFF) == PF_FILE_TYPE))
       && (pf->flags & PF_CURRENT_VALUE)
       && (pf->current_value.sval != NULL))
     SLFREE (pf->current_value.sval);
     
   SLFREE (pf);
}

/* The type can be b (bool), i (int), r (real), s (string), or f (file).
 * If the first character of type is a *, then it is a list.
 * The f-type has subtypes: e exists, n not-exists, r readable, w-writable
 */
static int parse_type (char *str, unsigned int *typep)
{
   unsigned int type = 0;
   
   /* There should be no spaces but we will allow it. */
   str = _pf_skip_whitespace (str);
   if (*str == 0)
     {
	pf_error ("Type information is missing.");
	PF_Errno = PF_CORRUPT_FIELD;
	return -1;
     }
   
   if (*str == '*')
     {
	type = PF_LIST_TYPE;
	str++;
     }
   
   if (0 == strcmp (str, "pset"))
     {
	/* FIXME!!!
	 * What is this?  It seems that SAO parameter files allow the type
	 * field to be "pset".  It appears to be the name of another 
	 * parameter file but their software treats it as a string.
	 * 
	 * It appears that a parameter of type pset will make the 
	 * parameters of the referenced file available to this parameter file.
	 * That is, it acts like an 'include' mechanism:
	 * 
	 *   foo,pset,h,xxx,,,etc      ==> include xxx.par
	 *   foo,pset,h,,,etc          ==> include foo.par
	 */
	str = "s";
     }

   switch (*str)
     {
      default:
	pf_error ("Field has unrecognized type. (%c)", *str);
	PF_Errno = PF_CORRUPT_FIELD;
	return -1;
	
      case 'b':
	type |= PF_BOOLEAN_TYPE;
	str++;
	break;
	
      case 'i':
	type |= PF_INTEGER_TYPE;
	str++;
	break;
	
      case 'r':
	type |= PF_REAL_TYPE;
	str++;
	break;

      case 'd':
	type |= PF_DOUBLE_TYPE;
	str++;
	break;

      case 's':
	type |= PF_STRING_TYPE;
	str++;
	break;

      case 'f':
	type |= PF_FILE_TYPE;
	str++;
	
	while (*str)
	  {
	     switch (*str)
	       {
		default:
		  pf_error ("Unrecogized subtype '%c' for 'f' type.", *str);
		  PF_Errno = PF_CORRUPT_FIELD;
		  return -1;
		  
		case 'e':
		  if (type & PF_FILE_NEXISTS)
		    {
		       pf_error ("Error in file type: both n and e specified.");
		       PF_Errno = PF_CORRUPT_FIELD;
		       return -1;
		    }
		  type |= PF_FILE_EXISTS;
		  break;
		  
		case 'n':
		  if (type & PF_FILE_EXISTS)
		    {
		       pf_error ("Error in file type: both n and e specified.");
		       PF_Errno = PF_CORRUPT_FIELD;
		       return -1;
		    }
		  type |= PF_FILE_NEXISTS;
		  break;
		  
		case 'r':
		  type |= PF_FILE_READABLE;
		  break;
		  
		case 'w':
		  type |= PF_FILE_WRITABLE;
		  break;
	       }
	     str++;
	  }
	break;
     }
   
   str = _pf_skip_whitespace (str);
   if (*str != 0)
     {
	pf_error ("Garbage at end of parameter type: '%s'", str);
	PF_Errno = PF_CORRUPT_FIELD;
	return -1;
     }
   
   *typep = type;
   return 0;
}

/* Modes are: q-query, a-auto, h-hidden, l-learn */
int _pf_parse_mode (char *str, unsigned int *modep)
{
   unsigned int mode = 0;
   
   /* There should be no spaces but we will allow it. */
   str = _pf_skip_whitespace (str);
   if (*str == 0)
     {
	pf_error ("Mode information is missing.");
	PF_Errno = PF_CORRUPT_FIELD;
	return -1;
     }
   
   while (*str)
     {
	switch (*str)
	  {
	   default:
	     pf_error ("Unrecognized mode (%c).", *str);
	     PF_Errno = PF_CORRUPT_FIELD;
	     return -1;
	     
	   case 'q':
	     mode |= PF_QUERY_MODE;
	     break;
	     
	   case 'a':
	     mode |= PF_AUTO_MODE;
	     break;
	
	   case 'h':
	     mode |= PF_HIDDEN_MODE;
	     break;
	
	   case 'l':
	     mode |= PF_LEARN_MODE;
	     break;
	  }
	str++;
	/* There should be no spaces but we will allow it. */
	str = _pf_skip_whitespace (str);
     }
   
   *modep = mode;
   return 0;
}

static Param_Type *create_param_type (char *line)
{
   Param_Type *p;
   char *name, *type, *mode, *min, *max, *value, *prompt;
   
   p = (Param_Type *) _pf_malloc (sizeof (Param_Type));
   if (p == NULL) return NULL;

   if (NULL == (p->name = _pf_create_string (line)))
     {
	goto free_and_return_error;
     }
   
   /* p->name is special.  See comment in pf.h */
   if (-1 == chop_line (line, &name, &type, &mode, &value, &min, &max, &prompt))
     {
	goto free_and_return_error;
     }
   
   p->name[strlen (name)] = 0;
   
   if (type == NULL)
     {
	p->type = PF_COMMENT_TYPE;
	return p;
     }
   
   if (-1 == parse_type (type, &p->type))
     {
	goto free_and_return_error;
     }
   
   if (-1 == _pf_parse_mode (mode, &p->mode))
     {
	goto free_and_return_error;
     }
   
   if (min != NULL)
     {
	if (NULL == (p->min = _pf_unescape_string (min)))
	  goto free_and_return_error;
     }
   
   if (max != NULL)
     {
	if (NULL == (p->max = _pf_unescape_string (max)))
	  goto free_and_return_error;
     }
   
   if (value != NULL)
     {
	if (NULL == (p->value = _pf_unescape_string (value)))
	  goto free_and_return_error;
	
	/* This is a weakness in the definition of parameter files.  There 
	 * seems to be NO WAY to have a parameter whose value begins with
	 * a ')' character.  Actually, I think there is but the existing 
	 * practice seems to defeat my solution.
	 */
	if (*p->value == ')')
	  p->flags |= PF_INDIRECT_VALUE;
     }

   if (prompt != NULL)
     {
	if (NULL == (p->prompt = _pf_unescape_string (prompt)))
	  goto free_and_return_error;
     }
   
   return p;
   
   free_and_return_error:
   
   if ((p != NULL) && (p->name != NULL))
     {
	pf_error ("Error processing parameter/line %s", p->name);
     }
   
   free_param_type (p);
   return NULL;
}

	
static Param_File_Type *create_param_file_type (char *file)
{
   Param_File_Type *p;
   
   if (NULL == (p = (Param_File_Type *) _pf_malloc (sizeof (Param_File_Type))))
     return NULL;
   
   if (NULL == (p->input_filename = _pf_create_string (file)))
     {
	SLFREE (p);
	return NULL;
     }
   
   return p;
}

void _pf_free_param_file (Param_File_Type *p)
{
   Param_Type *pf, *pf_next;
   
   if (p == NULL) return;
   pf = p->pf;
   while (pf != NULL)
     {
	pf_next = pf->next;
	free_param_type (pf);
	pf = pf_next;
     }
   if (p->input_filename != NULL) SLFREE (p->input_filename);
   if (p->output_filename != NULL) SLFREE (p->output_filename);
   SLFREE (p);
}


Param_File_Type *_pf_read_parm_file (char *file, FILE *fp)
{
   Param_File_Type *parm_file;
   Param_Type *pf, *last_pf;
   unsigned int line_num;
     
   char buf[PF_MAX_LINE_LEN];
   
   if (NULL == (parm_file = create_param_file_type (file)))
     {
	return NULL;
     }
   
   line_num = 0;
   last_pf = NULL;
   
   while (NULL != fgets (buf, sizeof(buf), fp))
     {
	unsigned int len;
	char *line;
	
	line_num++;
	line = _pf_skip_whitespace (buf);
	if (*line == 0) continue;
	
	/* Now knock off final newline */
	len = strlen (line) - 1;	       /* strlen (line) is > 0 */
	if (line[len] == '\n') line[len] = 0;

	if (NULL == (pf = create_param_type (line)))
	  {
	     _pf_free_param_file (parm_file);
	     pf_error ("Error on line %d of %s", line_num, file);
	     return NULL;
	  }
	
	pf->pfile = parm_file;
	
	if (last_pf == NULL)
	  {
	     parm_file->pf = pf;
	  }
	else last_pf->next = pf;
	
	last_pf = pf;
     }
   
   return parm_file;
}

/* If the object in the location pointed to by sval is non-NULL, it will be
 * freed.
 */
static int parse_string_according_to_type (char *str, unsigned int type,
					   int *ival, double *dval, char **sval)
{
   switch (type)
     {
      case PF_UINT_TYPE:
      case PF_INTEGER_TYPE:
	return _pf_parse_single_number (str, ival, NULL);
	
      case PF_BOOLEAN_TYPE:
	return _pf_parse_boolean (str, ival);
	
      case PF_REAL_TYPE:
      case PF_DOUBLE_TYPE:
	return _pf_parse_single_number (str, NULL, dval);
	
      case PF_STRING_TYPE:
      case PF_FILE_TYPE:
	/* This is easy */
	if (*sval != NULL) SLFREE (*sval);
	if (NULL == (*sval = _pf_create_string (str)))
	  return -1;
	
	return 0;
	
      default:
	pf_error ("type %d not implemented.", type);
	PF_Errno = PF_NOT_IMPLEMENTED;
     }
   return -1;
}


/* If type == 0, any parm matching name is returned. */
Param_Type *_pf_locate_param_by_type (Param_File_Type *p,
				      char *name, unsigned int type)
{
   Param_Type *pf;
   
   if (p == NULL) return NULL;
   pf = p->pf;
   
   while (pf != NULL)
     {
	if ((pf->name != NULL)
	    && !_pf_strcasecmp (name, pf->name))
	  {
	     if ((type == 0) 
		 || ((pf->type & 0xFF) == type))
	       return pf;
	     
	     /* FIXME!!! Make this more general! */
	     if ((type == PF_STRING_TYPE)
		 && ((pf->type & 0xFF) == PF_FILE_TYPE))
	       return pf;
	     if ((type == PF_FILE_TYPE)
		 && ((pf->type & 0xFF) == PF_STRING_TYPE))
	       return pf;

	     PF_Errno = PF_UNKNOWN_PARAMETER;
	     return NULL;
	  }
	pf = pf->next;
     }
   return pf;
}


static int do_indirect_shell_cmd (char *cmd, char **val)
{
   FILE *fp;
   char buf[4 * PF_MAX_LINE_LEN];
   char *b;
   unsigned int len;
   
   if (NULL == (fp = popen (cmd, "r")))
     {
	PF_Errno = PF_INDIRECT_ERROR;
	pf_error ("Unable to execute shell command %s", cmd);
	return -1;
     }
   
   len = fread (buf, 1, sizeof (buf)-1, fp);
   buf [len] = 0;

   (void) pclose (fp);

   /* Remove trailing new line characters */
   while (len && (buf[len-1] == '\n'))
     {
	len--;
	buf [len] = 0;
     }

   b = buf;
   while (*b != 0)
     {
	if (*b == '\n') *b = ' ';
	b++;
     }
   
   if (NULL == (*val = _pf_create_string (buf)))
     return -1;
   
   return 0;
}

static int get_indirect_object (Param_File_Type *p, char *name, 
				char **val, int what)
{
   char *file, *parm;
   Param_File_Type *new_p;
   Param_Type *pf;
   int ret;
   char *value;
   int is_file_indirect;

   if (*name != ')')
     {
	PF_Errno = PF_INDIRECT_ERROR;
	return -1;
     }
   name++;

   /* ))Shell Command */
   if (*name == ')')
     return do_indirect_shell_cmd (name + 1, val);

   /* Note: The CXC parameter file library appears to support an indirection
    * of the form ")FILE.PARM", which means to use the value of PARM from FILE.
    * (note added Oct-30-2001)
    */
   is_file_indirect = 0;
   if (*name == '*')
     {
	is_file_indirect = 1;
	name++;
     }
   if (NULL != _pf_strchr (name, '.')) /* HACK!!! Undocumented CXC parmeter file behavior. */
     is_file_indirect = 1;

   if (is_file_indirect)
     {
	if (NULL == (parm = _pf_strchr (name, '.')))
	  {
	     pf_error ("Indirection does not point to a parameter.");
	     PF_Errno = PF_CORRUPT_FIELD;
	     return -1;
	  }
	
	file = _pf_create_nstring (name, (unsigned int) (parm - name));
	if (file == NULL) return -1;
	parm++;			       /* skip . */
	
	new_p = p = pf_open_parameter_file (file, "r");
	if (p == NULL)
	  {
	     pf_error ("Unable to open indirect file %s.", file);
	     SLFREE (file);
	     PF_Errno = PF_INDIRECT_ERROR;
	     return -1;
	  }
	SLFREE (file);
	name = parm;
     }
   else new_p = NULL;
   
   if (NULL == (pf = _pf_locate_param_by_type (p, name, 0)))
     {
	if (new_p != NULL) pf_close_parameter_file (new_p);
	pf_error ("Bad indirect: %s", name);
	PF_Errno = PF_INDIRECT_ERROR;
	return -1;
     }
	
   ret = 0;
   
   if (what == 0) value = pf->value;
   else if (what == 1) value = pf->min;
   else value = pf->max;

   if (pf->flags & PF_DOING_INDIRECT)
     {
	pf_error ("Indirection loop found.");
	PF_Errno = PF_INDIRECT_ERROR;
	ret = -1;
     }
   else if (value == NULL) *val = NULL;
   else if (*value == ')')
     {
	pf->flags |= PF_DOING_INDIRECT;
	ret = get_indirect_object (p, value, val, what);
	pf->flags &= ~PF_DOING_INDIRECT;
     }
   else 
     {
	if (NULL == (*val = _pf_create_string (value)))
	  ret = -1;
     }
   
   if (new_p != NULL) pf_close_parameter_file (new_p);
   return ret;
}

int _pf_get_indirect_value (Param_File_Type *p, char *name, char **val)
{
   return get_indirect_object (p, name, val, 0);
}

static int get_indirect_min (Param_File_Type *p, char *name, char **val)
{
   return get_indirect_object (p, name, val, 1);
}

static int get_indirect_max (Param_File_Type *p, char *name, char **val)
{
   return get_indirect_object (p, name, val, 2);
}

static int enumerated_value_check (Param_Type *pf, char *list)
{
   unsigned int type;
   int ival;
   double dval;
   char *value;
   unsigned int nth = 0;
   type = pf->type & 0xFF;
   
   while (NULL != (value = _pf_extract_string_element (list, "|", &nth)))
     {
	if (*value == 0)
	  {
	     SLFREE (value);
	     continue;
	  }
	
	switch (type)
	  {
	   case PF_UINT_TYPE:
	     (void) _pf_parse_single_number (value, &ival, NULL);
	     SLFREE (value);
	     if ((unsigned int) ival == pf->current_value.uval) return 0;
	     break;

	   case PF_INTEGER_TYPE:
	     (void) _pf_parse_single_number (value, &ival, NULL);
	     SLFREE (value);
	     if (ival == pf->current_value.ival) return 0;
	     break;
	     
	
	   case PF_REAL_TYPE:
	   case PF_DOUBLE_TYPE:
	     (void) _pf_parse_single_number (value, NULL, &dval);
	     SLFREE (value);
	     if (dval == pf->current_value.dval) return 0;
	     break;
	
	   case PF_STRING_TYPE:
	   case PF_FILE_TYPE:
	     if (!strcmp (value, pf->current_value.sval))
	       {
		  SLFREE (value);
		  return 0;
	       }
	     SLFREE (value);
	     break;
	     
	   default:
	     SLFREE (value);
	     return 0;
	  }
     }

   PF_Errno = PF_RANGE_ERROR;
   return -1;
}

static int range_check (Param_File_Type *p, Param_Type *pf)
{
   char *min, *max;
   int imin, imax;
   double dmin, dmax;
   unsigned int type;
   
   
   if ((pf->flags & PF_CURRENT_VALUE) == 0)
     {
	pf_error ("Unable to check range of missing value for %s.", pf->name);
	PF_Errno = PF_BAD_ARGUMENT;
	return -1;
     }
   
   min = pf->min;
   max = pf->max;
   type = pf->type & 0xFF;
   
   if ((min != NULL) && (*min == ')'))
     if (-1 == get_indirect_min (p, pf->min, &min))
       return -1;
   
   if ((max != NULL) && (*max == ')'))
     if (-1 == get_indirect_max (p, pf->max, &max))
       return -1;
   
   if ((min == NULL) && (max == NULL))
     return 0;
   
   if ((min != NULL)
       && (_pf_strchr (min, '|') || (type == PF_STRING_TYPE)))
     return enumerated_value_check (pf, min);
   
   if ((max != NULL)
       && (_pf_strchr (max, '|') || (type == PF_STRING_TYPE)))
     return enumerated_value_check (pf, max);
   
       
   switch (type)
     {
      case PF_INTEGER_TYPE:
	if (min != NULL)
	  {
	     if (-1 == _pf_parse_single_number (min, &imin, NULL))
	       return -1;
	     if (pf->current_value.ival < imin)
	       break;
	  }
	
	if (max != NULL)
	  {
	     if (-1 == _pf_parse_single_number (max, &imax, NULL))
	       return -1;
	     if (pf->current_value.ival > imax)
	       break;
	  }
	return 0;
	
      case PF_UINT_TYPE:
	if (min != NULL)
	  {
	     if (-1 == _pf_parse_single_number (min, &imin, NULL))
	       return -1;
	     if (pf->current_value.uval < (unsigned int) imin)
	       break;
	  }
	
	if (max != NULL)
	  {
	     if (-1 == _pf_parse_single_number (max, &imax, NULL))
	       return -1;
	     if (pf->current_value.uval > (unsigned int) imax)
	       break;
	  }
	return 0;
	
      case PF_DOUBLE_TYPE:
      case PF_REAL_TYPE:
	if (min != NULL)
	  {
	     if (-1 == _pf_parse_single_number (min, NULL, &dmin))
	       return -1;
	     if (pf->current_value.dval < dmin)
	       break;
	  }
	
	if (max != NULL)
	  {
	     if (-1 == _pf_parse_single_number (max, NULL, &dmax))
	       return -1;
	     if (pf->current_value.dval > dmax)
	       break;
	  }
	return 0;
	
      case PF_FILE_TYPE:
	/* The file type should also check whether or not the other
	 * file conditions are met instead of this silly test below.
	 */
      case PF_STRING_TYPE:
	if ((min != NULL) &&
	    (strcmp (min, pf->current_value.sval) > 0))
	  break;
	
	if ((max != NULL) && 
	    (strcmp (max, pf->current_value.sval) < 0))
	  break;
	return 0;
	
      default:
	return 0;
     }

   PF_Errno = PF_RANGE_ERROR;
   return -1;
}

unsigned int _pf_get_effective_mode (Param_File_Type *p, Param_Type *pf)
{	
   unsigned int mode;
   
   
   mode = pf->mode;
   if (mode & PF_AUTO_MODE)
     mode = p->mode;

   /* Apparantly the SAO parameter file library _ignores_ the query 
    * mode on a parameter if the mode of the parameter file is h.  I am 
    * _not_ going to implement this broken feature.  Instead, I am going
    * to create a new mode flag that the application can call to set.
    */
   if (p->mode & PF_NEVER_QUERY_MODE)
     mode &= ~PF_QUERY_MODE;

   if (p->mode & PF_FORCE_LEARN_MODE)
     mode |= PF_LEARN_MODE;

   return mode;
}


/* Set the current_value field of pf.  This does not perform any range
 * checking.  Such checking should be performed by the calling routine.
 * This way, the min,max fields will not be parsed and if any parse error
 * occurs, it occurs because of the value.
 * 
 * If the PF_CURRENT_VALUE bit is not set upon return, the current_value
 * is not set and must be queried by the calling routine.  This routine
 * does not query.
 * 
 * TODO: fix this for list types
 */
static int _pf_get_current_value (Param_File_Type *p, Param_Type *pf)
{
   char *value;
   
   if (pf->flags & PF_CURRENT_VALUE) return 0;
   
   if (NULL == pf->value)
     return 0;

   if (pf->flags & PF_INDIRECT_VALUE)
     {
	if (-1 == _pf_get_indirect_value (p, pf->value, &value))
	  return -1;
	if (value == NULL) return 0;
     }
   else value = pf->value;
   
   if (-1 ==  parse_string_according_to_type (value, 
					      pf->type & 0xFF,
					      &pf->current_value.ival,
					      &pf->current_value.dval,
					      &pf->current_value.sval))
     {
	if (pf->flags & PF_INDIRECT_VALUE) SLFREE (value);
	return -1;
     }
   
   pf->flags |= PF_CURRENT_VALUE;
   if (pf->flags & PF_INDIRECT_VALUE) SLFREE (value);
   return 0;
}

static int update_value (Param_Type *pf)
{
   char buf[PF_MAX_LINE_LEN];
   char *b;
   
   if ((pf->flags & PF_CURRENT_VALUE) == 0) return 0;
   if (pf->flags & PF_INDIRECT_VALUE) return 0;
   if (pf->type & PF_LIST_TYPE) return 0;
   
   b = buf;
   switch (pf->type & 0xFF)
     {
      case PF_INTEGER_TYPE:
	sprintf (buf, "%d", pf->current_value.ival);
	break;

      case PF_UINT_TYPE:
	sprintf (buf, "%u", pf->current_value.uval);
	break;
	
      case PF_REAL_TYPE:
      case PF_DOUBLE_TYPE:
	sprintf (buf, "%.16g", pf->current_value.dval);
	break;
	
      case PF_FILE_TYPE:
      case PF_STRING_TYPE:
	b = pf->current_value.sval;
	if (strlen (b) >= PF_MAX_LINE_LEN)
	  {
	     pf_error ("String value is too long for parameter file.");
	     PF_Errno = PF_BAD_ARGUMENT;
	     return -1;
	  }
	break;
	
      case PF_BOOLEAN_TYPE:
	strcpy (buf, (pf->current_value.ival ? "yes" : "no"));
	break;
	
      default:
	pf_error ("update_value: Type %c is not supported", pf->type & 0xFF);
	PF_Errno = PF_NOT_IMPLEMENTED;
	return -1;
     }
   if (pf->value != NULL) 
     SLFREE (pf->value);
   
   if (NULL == (pf->value = _pf_create_string (b)))
     return -1;
   
   return 0;
}

int _pf_get_value (Param_File_Type *p, Param_Type *pf, unsigned int ormode)
{
   int perform_query;
   unsigned int mode;

   if (pf->type == PF_COMMENT_TYPE)
     return 0;

   mode = _pf_get_effective_mode (p, pf);
   perform_query = ((0 == (pf->flags & PF_CURRENT_VALUE))
		    && (mode & PF_QUERY_MODE));
   mode |= ormode;
   if (perform_query == 0)
     perform_query = ormode & PF_QUERY_MODE;   /* explicit requested */

   if (-1 == _pf_get_current_value (p, pf))
     {
	if ((PF_Errno != PF_NUMBER_FORMAT_BAD)
	    && (PF_Errno != PF_INDIRECT_ERROR))
	  return -1;
	/* reset error and force a query */
	pf_error ("Error occured while examining %s from %s.",
		  pf->name, p->input_filename);
	PF_Errno = 0;
	perform_query = 1;
     }
   
   if (pf->flags & PF_CURRENT_VALUE)
     {
	if (-1 == range_check (p, pf))
	  {
	     pf_error ("Range error occured while looking at parameter %s.",
		       pf->name);
	     if (PF_Errno != PF_RANGE_ERROR)
	       return -1;

	     /* Force query. */
	     perform_query = 1;
	  }
     }
   else perform_query = 1;
   
   while (perform_query)
     {
	if (-1 == _pf_query_current_value (p, pf))
	  return -1;
	
	if (0 == range_check (p, pf)) break;
	if (PF_Errno != PF_RANGE_ERROR)
	  return -1;
	pf_error ("Value is out of range.");
     }
   
   if ((mode & PF_LEARN_MODE) && (pf->flags & PF_PARAM_DIRTY))
     {
	p->flags |= PFILE_DIRTY;
	if (-1 == update_value (pf))
	  return -1;
     }
   return 0;
}


static char map_type_to_char (unsigned int type)
{
   switch (type)
     {
      case PF_UINT_TYPE:
	/* drop */
      case PF_INTEGER_TYPE:
	return 'i';
      case PF_BOOLEAN_TYPE:
	return 'b';
      case PF_FILE_TYPE:
	return 'f';
      case PF_STRING_TYPE:
	return 's';
      case PF_REAL_TYPE:
	return 'r';
      case PF_DOUBLE_TYPE:
	return 'd';
     }
   return '?';
}

static Param_Type *get_object (Param_File_Type *p,
			       char *name,
			       unsigned int type,
			       unsigned int mode)
			       
{
   Param_Type *pf;
   
   if (p == NULL) return NULL;
   if (NULL == (pf = _pf_locate_param_by_type (p, name, type)))
     {
	pf_error ("%s:\nError locating parameter named '%s' of type '%c'",
		  p->input_filename, name, map_type_to_char (type));
	
	return NULL;
     }
   
   if (-1 == _pf_get_value (p, pf, mode))
     return NULL;

   return pf;
}

int pf_get_integer (Param_File_Type *p, char *name, int *ival)
{
   Param_Type *pf;
   
   if (NULL != (pf = get_object (p, name, PF_INTEGER_TYPE, 0)))
     {
	*ival = pf->current_value.ival;
	return 0;
     }
   return -1;
}

int pf_get_uint (Param_File_Type *p, char *name, unsigned int *ival)
{
   Param_Type *pf;
   
   if (NULL != (pf = get_object (p, name, PF_INT_TYPE, 0)))
     {
	*ival = pf->current_value.uval;
	return 0;
     }
   return -1;
}

int pf_get_boolean (Param_File_Type *p, char *name, int *ival)
{
   Param_Type *pf;
   
   if (NULL != (pf = get_object (p, name, PF_BOOLEAN_TYPE, 0)))
     {
	*ival = pf->current_value.ival;
	return 0;
     }
   return -1;
}

int pf_get_string (Param_File_Type *p, char *name, char *sval, unsigned int len)
{
   Param_Type *pf;
   
   if (NULL != (pf = get_object (p, name, PF_STRING_TYPE, 0)))
     {
	strncpy (sval, pf->current_value.sval, len);
	if (len) sval[len - 1] = 0;
	return 0;
     }
   return -1;
}

int pf_get_file (Param_File_Type *p, char *name, char *sval, unsigned int len)
{
   Param_Type *pf;
   
   if (NULL != (pf = get_object (p, name, PF_FILE_TYPE, 0)))
     {
	strncpy (sval, pf->current_value.sval, len);
	if (len) sval[len-1] = 0;
	return 0;
     }
   return -1;
}

int pf_get_double (Param_File_Type *p, char *name, double *dval)
{
   Param_Type *pf;
   
   if (NULL != (pf = get_object (p, name, PF_REAL_TYPE, 0)))
     {
	*dval = pf->current_value.dval;
	return 0;
     }
   return -1;
}

int pf_parameter_exists (Param_File_Type *p, char *name)
{
   return (NULL != _pf_locate_param_by_type (p, name, 0));
}

int pf_get_type (Param_File_Type *p, char *name)
{
   Param_Type *pf = _pf_locate_param_by_type (p, name, 0);
   if (pf == NULL)
     return -1;
   
   return (int) pf->type & 0xFF;
}

