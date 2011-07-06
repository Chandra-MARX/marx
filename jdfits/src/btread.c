/*
    Copyright (C) 2002 MIT Center For Space Research

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

#include <stdarg.h>
#include <string.h>
#if HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdfits.h"
#include "_jdfits.h"

static int is_binary_table (JDFits_Type *f)
{
   if ((f == NULL) || (f->header == NULL))
     return 0;

   return (f->header->type == JDFITS_BINTABLE_HEADER);
}

JDFits_Type *jdfits_find_binary_table (char *file,
				       int (*fun)(void *, JDFits_Type *),
				       void *clientdata)
{
   JDFits_Type *ft;

   if (NULL == (ft = jdfits_open_file (file, JDFITS_READ_MODE)))
     return NULL;

   /* Now find the extension */
   while (1)
     {
	JDFits_Header_Type *h = ft->header;

	if (h->type == JDFITS_BINTABLE_HEADER)
	  {
	     int status;

	     if (-1 == jdfits_bintable_parse_headers (ft))
	       {
		  jdfits_close_file (ft);
		  return NULL;
	       }

	     status = (*fun) (clientdata, ft);
	     if (status == -1)
	       {
		  jdfits_close_file (ft);
		  return NULL;
	       }

	     if (status == 1)
	       return ft;
	  }

	if ((-1 == jdfits_skip_to_next_header (ft))
	    || (-1 == jdfits_read_header (ft)))
	  {
	     jdfits_error ("Unable to find proper binary table in %s", file);
	     jdfits_close_file (ft);
	     return NULL;
	  }
     }
}

static int extname_match_fun (void *vstr, JDFits_Type *ft)
{
   char *str = (char *) vstr;
   char *val;
   int status;

   if (-1 == jdfits_read_keyword_string (ft, "EXTNAME", &val))
     return 0;

   status = (0 == jdfits_strcasecmp (val, str));
   jdfits_free (val);

   return status;
}

JDFits_Type *jdfits_open_binary_table (char *file, char *extname)
{
   JDFits_Type *ft;

   if (NULL == (ft = jdfits_find_binary_table (file, extname_match_fun, (void *) extname)))
     jdfits_error ("Unable to find binary table with EXTNAME=%s in %s", extname, file);
   return ft;
}

static JDFits_Bintable_Field_Type *find_column (JDFits_Bintable_Type *fbt, char *name, unsigned int *ofsp)
{
   JDFits_Bintable_Field_Type *finfo, *finfo_max;
   unsigned int ofs;

   finfo = fbt->finfo;
   finfo_max = finfo + fbt->tfields;

   ofs = 0;
   while (finfo < finfo_max)
     {
	char *ttype;

	ttype = finfo->ttype;
	if (0 == jdfits_strcasecmp (name, ttype))
	  {
	     *ofsp = ofs;
	     return finfo;
	  }

	ofs += finfo->repeat * (finfo->size / 8);
	finfo++;
     }
   return NULL;
}

int jdfits_bintable_column_exists (JDFits_Type *f, char *column)
{
   unsigned int unused;

   if (1 != is_binary_table (f))
     {
	jdfits_error ("JDFits_Type pointer is not a binary table");
	return -1;
     }

   if (-1 == jdfits_bintable_parse_headers (f))
     return -1;

   return (NULL != find_column (f->header->ext.bintable, column, &unused));
}

JDFits_Bintable_Field_Type *_jdfits_bintable_find_column (JDFits_Bintable_Type *fbt, char *name, unsigned int *ofsp)
{
   JDFits_Bintable_Field_Type *f;

   if (NULL == (f = find_column (fbt, name, ofsp)))
     jdfits_error ("Column name %s is not in table", name);

   return f;
}

void jdfits_simple_close_btable (JDFits_BTable_Read_Type *fbtrt)
{
   if (fbtrt == NULL) return;

   if (fbtrt->ft != NULL) jdfits_close_file (fbtrt->ft);
   jdfits_free ((char *) fbtrt->row_data);
   jdfits_free ((char *) fbtrt->data_types);
   jdfits_free ((char *) fbtrt->data_offsets);
   jdfits_free ((char *) fbtrt);
}

JDFits_BTable_Read_Type *jdfits_simple_aopen_btable (char *file,
						     char *extname,
						     unsigned int ncols,
						     char **column_names)
{
   JDFits_BTable_Read_Type *fbtrt;
   JDFits_Bintable_Type *fbt;
   JDFits_Type *ft;
   JDFits_Header_Type *h;
   unsigned int i;

   if ((extname == NULL) || (ncols == 0))
     return NULL;

   fbtrt = (JDFits_BTable_Read_Type *) jdfits_malloc (sizeof (JDFits_BTable_Read_Type));
   if (fbtrt == NULL)
     return NULL;
   memset ((char *) fbtrt, 0, sizeof (JDFits_BTable_Read_Type));

   if (NULL == (fbtrt->ft = jdfits_open_binary_table (file, extname)))
     goto return_error;

   ft = fbtrt->ft;
   h = ft->header;

   /* Extension found.  Check for consistency */

   fbt = h->ext.bintable;
   if ((unsigned int) fbt->tfields < ncols)
     {
	jdfits_error ("Expecting TFIELDS >= %u", ncols);
	goto return_error;
     }

   fbtrt->row_data_len = fbt->naxis1;
   fbtrt->num_rows = fbt->naxis2;
   fbtrt->num_rows_to_read = fbt->naxis2;

   if (NULL == (fbtrt->row_data = (unsigned char *) jdfits_malloc (fbtrt->row_data_len)))
     goto return_error;
   if (NULL == (fbtrt->data_types = (unsigned int *) jdfits_malloc (ncols * sizeof(unsigned int))))
     goto return_error;
   if (NULL == (fbtrt->data_offsets = (unsigned int *) jdfits_malloc (ncols * sizeof (unsigned int))))
     goto return_error;

   fbtrt->num_data_types = ncols;

   /* Make sure that the requested columns exist and that they are
    * simple scalars.
    */

   for (i = 0; i < ncols; i++)
     {
	char *colname;
	JDFits_Bintable_Field_Type *finfo;
	unsigned int ofs;

	colname = column_names[i];

	if (NULL == (finfo = _jdfits_bintable_find_column (fbt, colname, &ofs)))
	  goto return_error;

	if (finfo->repeat != 1)
	  {
	     jdfits_error ("Column %s does not have repeat count of 1",
			   colname);
	     goto return_error;
	  }

	fbtrt->data_offsets[i] = ofs;
	fbtrt->data_types[i] = finfo->type;
     }

   if (-1 == jdfits_read_open_data (fbtrt->ft))
     goto return_error;

   return fbtrt;

   /* Get here only if something went wrong */
   return_error:

   jdfits_simple_close_btable (fbtrt);
   return NULL;
}

JDFits_BTable_Read_Type *jdfits_simple_open_btable (char *file,
						    char *extname,
						    unsigned int ncols, ...)

{
   JDFits_BTable_Read_Type *bt;
   va_list ap;
   unsigned int i;
   char **column_names;

   column_names = (char **) jdfits_malloc (ncols * sizeof (char *));
   if (column_names == NULL)
     return NULL;

   va_start (ap, ncols);
   for (i = 0; i < ncols; i++)
     column_names [i] = va_arg (ap, char *);
   va_end (ap);

   bt = jdfits_simple_aopen_btable (file, extname, ncols, column_names);

   jdfits_free ((char *) column_names);

   return bt;
}

int jdfits_simple_d_read_btable (JDFits_BTable_Read_Type *bt, double *buf)
{
   unsigned int nbytes;
   unsigned char *data;
   unsigned int *data_offsets;
   unsigned int *data_types;
   unsigned int num_data_types;
   unsigned int i;

   if (bt == NULL) return -1;
   if (bt->num_rows_to_read == 0)
     return -1;

   nbytes = bt->row_data_len;
   data = bt->row_data;

   if (nbytes != jdfits_read_bytes (bt->ft, data, nbytes))
     return -1;

   bt->num_rows_to_read -= 1;

   num_data_types = bt->num_data_types;
   data_types = bt->data_types;
   data_offsets = bt->data_offsets;

   for (i = 0; i < num_data_types; i++)
     {
	int32 i32;
	int16 i16;
	float32 f32;
	float64 f64;
	double d;
	unsigned char *data_bytes;

	data_bytes = data + data_offsets [i];

	switch (data_types[i])
	  {
	   case JDFITS_INT32_TYPE:
	     (void) jdfits_str_read_int32 (&i32, 1, data_bytes);
	     d = (double) i32;
	     break;

	   case JDFITS_INT16_TYPE:
	     (void) jdfits_str_read_int16 (&i16, 1, data_bytes);
	     d = (double) i16;
	     break;

	   case JDFITS_FLOAT32_TYPE:
	     (void) jdfits_str_read_float32 (&f32, 1, data_bytes);
	     d = (double) f32;
	     break;

	   case JDFITS_FLOAT64_TYPE:
	     (void) jdfits_str_read_float64 (&f64, 1, data_bytes);
	     d = (double) f64;
	     break;

	   default:
	     jdfits_error ("jdfits_simple_d_read_btable: data type unsupported");
	     return -1;
	  }
	buf[i] = d;
     }

   return 0;
}

void jdfits_bintable_close_rows (JDFits_Row_Type *rt)
{
   if (rt == NULL)
     return;

   if (rt->col_data != NULL)
     {
	JDFits_Col_Data_Type *cd;
	unsigned int i, num;

	num = rt->num_columns;
	cd = rt->col_data;

	for (i = 0; i < num; i++)
	  {
	     if (cd[i].data.a != NULL)
	       jdfits_free (cd[i].data.a);
	  }

	jdfits_free ((char *)cd);
     }

   if (rt->row_bytes != NULL) jdfits_free ((char *) rt->row_bytes);
   jdfits_free ((char *) rt);
}

typedef unsigned char *(*Read_Fun_Type)(void *, unsigned int, unsigned char *);

static unsigned char *read_string_bytes (char *buf, unsigned int n, unsigned char *s)
{
   memcpy (buf, (char *)s, n);
   /* An extra byte has been allocated for this */
   buf[n] = 0;
   return s + n;
}

static Read_Fun_Type get_type_read_fun (int from_type, int to_type, unsigned int *sizep)
{
   switch (to_type | 0x20)
     {
      case 'a':
	*sizep = 1;
	switch (from_type)
	  {
	   case JDFITS_STRING_TYPE:
	     return (Read_Fun_Type)read_string_bytes;
	  }
	break;

      case 's':
	*sizep = sizeof (short);
	switch (from_type)
	  {
	   case JDFITS_INT16_TYPE: return (Read_Fun_Type)jdfits_read_int16_short;
	   case JDFITS_INT32_TYPE: return (Read_Fun_Type)jdfits_read_int32_short;
	   case JDFITS_FLOAT32_TYPE: return (Read_Fun_Type)jdfits_read_float32_short;
	   case JDFITS_FLOAT64_TYPE: return (Read_Fun_Type)jdfits_read_float64_short;
	  }
	break;

      case 'i':
	*sizep = sizeof (int);
	switch (from_type)
	  {
	   case JDFITS_INT16_TYPE: return (Read_Fun_Type)jdfits_read_int16_int;
	   case JDFITS_INT32_TYPE: return (Read_Fun_Type)jdfits_read_int32_int;
	   case JDFITS_FLOAT32_TYPE: return (Read_Fun_Type)jdfits_read_float32_int;
	   case JDFITS_FLOAT64_TYPE: return (Read_Fun_Type)jdfits_read_float64_int;
	  }
	break;

      case 'l':
	*sizep = sizeof (long);
	switch (from_type)
	  {
	   case JDFITS_INT16_TYPE: return (Read_Fun_Type)jdfits_read_int16_long;
	   case JDFITS_INT32_TYPE: return (Read_Fun_Type)jdfits_read_int32_long;
	   case JDFITS_FLOAT32_TYPE: return (Read_Fun_Type)jdfits_read_float32_long;
	   case JDFITS_FLOAT64_TYPE: return (Read_Fun_Type)jdfits_read_float64_long;
	  }
	break;

      case 'f':
	*sizep = sizeof (float);
	switch (from_type)
	  {
	   case JDFITS_INT16_TYPE: return (Read_Fun_Type)jdfits_read_int16_float;
	   case JDFITS_INT32_TYPE: return (Read_Fun_Type)jdfits_read_int32_float;
	   case JDFITS_FLOAT32_TYPE: return (Read_Fun_Type)jdfits_read_float32_float;
	   case JDFITS_FLOAT64_TYPE: return (Read_Fun_Type)jdfits_read_float64_float;
	  }
	break;

      case 'd':
	*sizep = sizeof (double);
	switch (from_type)
	  {
	   case JDFITS_INT16_TYPE: return (Read_Fun_Type)jdfits_read_int16_double;
	   case JDFITS_INT32_TYPE: return (Read_Fun_Type)jdfits_read_int32_double;
	   case JDFITS_FLOAT32_TYPE: return (Read_Fun_Type)jdfits_read_float32_double;
	   case JDFITS_FLOAT64_TYPE: return (Read_Fun_Type)jdfits_read_float64_double;
	  }
	break;
     }

   jdfits_error ("Conversion from type %d to %c unsupported");
   return NULL;
}

static int get_column_name_and_type (char *s, char **np, int *tp)
{
   *tp = *s;
   s = strchr (s, ':');
   if (s == NULL)
     {
	jdfits_error ("Expecting column to be specified as TYPE:NAME");
	return -1;
     }
   *np = _jdfits_skip_whitespace (s+1);
   return 0;
}

JDFits_Row_Type *
  jdfits_bintable_aopen_rows (JDFits_Type *ft, unsigned int ncols, char **column_names)
{
   JDFits_Row_Type *rt;
   JDFits_Header_Type *h;
   JDFits_Bintable_Type *fbt;
   unsigned int i;

   if ((ft == NULL) || (ncols == 0))
     return NULL;

   if (1 != is_binary_table (ft))
     return NULL;

   rt = (JDFits_Row_Type *) jdfits_malloc (sizeof (JDFits_Row_Type));
   if (rt == NULL)
     return NULL;

   memset ((char *) rt, 0, sizeof (JDFits_Row_Type));

   h = ft->header;
   fbt = h->ext.bintable;

   if ((unsigned int) fbt->tfields < ncols)
     {
	jdfits_error ("Expecting TFIELDS >= %u", ncols);
	goto return_error;
     }

   rt->num_bytes = fbt->naxis1;
   rt->num_rows = fbt->naxis2;
   rt->num_rows_to_read = fbt->naxis2;
   rt->num_columns = ncols;

   if (NULL == (rt->row_bytes = (unsigned char *) jdfits_malloc (rt->num_bytes)))
     goto return_error;

   if (NULL == (rt->col_data = (JDFits_Col_Data_Type *) jdfits_malloc (ncols * sizeof (JDFits_Col_Data_Type))))
     goto return_error;

   memset ((char *) rt->col_data, 0, ncols * sizeof (JDFits_Col_Data_Type));

   /* Make sure that the requested columns exist. */

   for (i = 0; i < ncols; i++)
     {
	char *colname;
	JDFits_Bintable_Field_Type *finfo;
	unsigned int ofs;
	JDFits_Col_Data_Type *cd;
	int to_type;
	unsigned int size;

	cd = rt->col_data + i;

	if (-1 == get_column_name_and_type (column_names[i], &colname, &to_type))
	  goto return_error;

	if (NULL == (finfo = _jdfits_bintable_find_column (fbt, colname, &ofs)))
	  goto return_error;

	if (NULL == (cd->read_fun = get_type_read_fun (finfo->type, to_type, &size)))
	  goto return_error;

	cd->repeat = finfo->repeat;
	cd->data_offset = ofs;
	cd->data_type = to_type;

	/* allow 1 byte to null terminate strings */
	if (NULL == (cd->data.a = jdfits_malloc (1 + finfo->repeat * size)))
	  goto return_error;
     }

   if (-1 == jdfits_read_open_data (ft))
     goto return_error;

   return rt;

   /* Get here only if something went wrong */
   return_error:
   jdfits_bintable_close_rows (rt);
   return NULL;
}

int jdfits_read_next_row (JDFits_Type *f, JDFits_Row_Type *r)
{
   unsigned int nbytes;
   unsigned char *data;
   JDFits_Col_Data_Type *cd, *cd_max;

   if ((f == NULL) || (r == NULL))
     return -1;

   if (r->num_rows_to_read == 0)
     return 0;

   nbytes = r->num_bytes;
   data = r->row_bytes;

   if (nbytes != jdfits_read_bytes (f, data, nbytes))
     return -1;

   r->num_rows_to_read -= 1;

   cd = r->col_data;
   cd_max = cd + r->num_columns;

   while (cd < cd_max)
     {
	unsigned char *data_bytes;

	data_bytes = data + cd->data_offset;
	(void) cd->read_fun ((void *)cd->data.a, cd->repeat, data + cd->data_offset);

	cd++;
     }

   return 1;
}

JDFits_Row_Type *
  jdfits_bintable_open_rows (JDFits_Type *ft, unsigned int ncols, ...)
{
   JDFits_Row_Type *r;
   va_list ap;
   unsigned int i;
   char **column_names;

   column_names = (char **) jdfits_malloc (ncols * sizeof (char *));
   if (column_names == NULL)
     return NULL;

   va_start (ap, ncols);
   for (i = 0; i < ncols; i++)
     column_names [i] = va_arg (ap, char *);
   va_end (ap);

   r = jdfits_bintable_aopen_rows (ft, ncols, column_names);

   jdfits_free ((char *) column_names);

   return r;
}
