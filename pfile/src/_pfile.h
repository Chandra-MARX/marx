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
#ifndef _MITCSR_PFILE_PRIVATE_H_
#define _MITCSR_PFILE_PRIVATE_H_
struct _Param_Type
{
   char *name;
   /* This holds the name of the parameter.  The value is formed from the
    * line in the file with a 0 at the first comma.  For comments, this is
    * the comment.
    */

   unsigned int type;
   /* type is partially bitmapped.  type & 0xFF give the basic type.  The
    * other bits are flags that tend to modify the type.  Pointer types, i.e.,
    * indirect types, are not considered as modifiers.  The type given here
    * actually refers to what is listed in the type field in the parameter
    * file and is what is used to access the current_value union below.
    */
   /* See PF_INTEGER_TYPE and company in public header file */

#define PF_FILE_EXISTS		0x100
#define PF_FILE_NEXISTS		0x200
#define PF_FILE_READABLE	0x400
#define PF_FILE_WRITABLE	0x800
#define PF_LIST_TYPE		0x1000

   unsigned int mode;
#define PF_AUTO_MODE		0x1
#define PF_QUERY_MODE		0x2
#define PF_HIDDEN_MODE		0x4
#define PF_LEARN_MODE		0x8
#define PF_NEVER_QUERY_MODE	0x10   /* only attached to file's mode */
#define PF_FORCE_LEARN_MODE	0x20   /* learn even if parameter is hidden */

   unsigned int flags;
#define PF_INDIRECT_VALUE	0x01
#define PF_PARAM_DIRTY		0x10
#define PF_CURRENT_VALUE	0x20
#define PF_DOING_INDIRECT	0x1000
   /* These three character fields actually refer to string representations
    * (but lacking quoting) of what would actually be written out to a file
    * if this parameter is learned.  The value string will only be referred
    * to when the PF_CURRENT_VALUE bit is not set.  If it is set, the
    * appropriate member of the current_value member will be used.
    *
    * The min, max fields will always be used if non-null.  If min refers to
    * an enumerated value, then max will be null.  If min/max are indirect
    * types, then whether or not the value is enumerated will not be known
    * until the reference is looked at.  As a result, the min/max fields are
    * not parsed until the value for the parameter is computed.  Thus there
    * is nothing in the flags that can be examined to derive this information.
    */
   char *value, *min, *max;

   union
     {
	unsigned int uval;
	int ival;
	char *sval;
	double dval;
     }
   current_value;

   char *prompt;
   struct _Param_Type *next;
   struct _Param_File_Type *pfile;
};

struct _Param_File_Type
{
   unsigned int mode;		       /* default mode: query | learn */
   int num_references;
   unsigned int flags;
#define PFILE_DIRTY		0x1
#define PFILE_WRITE_MODE	0x2
#define PFILE_WRITE_ON_CLOSE	0x4
   char *output_filename;		       /* file to be written out. */
   char *input_filename;	       /* file read in */
   Param_Type *pf;		       /* pointers to parameters in file */
   struct _Param_File_Type *next;
};

extern char *_pf_skip_whitespace (char *);
extern char *_pf_create_string (char *);
extern char *_pf_create_nstring (char *, unsigned int);
extern char *_pf_malloc (unsigned int);
extern Param_File_Type *_pf_read_parm_file (char *, FILE *);
extern void _pf_free_param_file (Param_File_Type *);
extern void _pf_free_current_value (Param_Type *);
extern int _pf_parse_parm (Param_Type *);
extern int _pf_parse_mode (char *, unsigned int *);
extern char *_pf_strchr (char *, int);
extern char *_pf_rstrchr (char *, int);
extern char *_pf_strcat (char *, char *);
extern int _pf_parse_single_number (char *, int *, double *);
extern int _pf_query_current_value (Param_File_Type *, Param_Type *);
extern int _pf_parse_boolean (char *, int *);
extern int _pf_value_in_range (Param_Type *);
extern int _pf_get_value (Param_File_Type *, Param_Type *, unsigned int);
extern unsigned int _pf_get_effective_mode (Param_File_Type *, Param_Type *);
extern char *_pf_unescape_string (char *);
extern char *_pf_escape_string (char *);
extern char *_pf_extract_string_element (char *, char *, unsigned int *);
extern char *_pf_strbrk (char *, char *);
extern Param_Type *_pf_locate_param_by_type (Param_File_Type *, char *, 
					     unsigned int);
extern int _pf_get_indirect_value (Param_File_Type *, char *, char **);

extern int _pf_strcasecmp (char *, char *);
#endif				       /* _MITCSR_PFILE_PRIVATE_H_ */
