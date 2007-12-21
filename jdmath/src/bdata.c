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

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>

#include "jdmath.h"
#include "_jdmath.h"

static unsigned char Magic_Bytes[4] = 
{
   0x83, 0x43, 0x59, 0x8E
};

#define BDATA_HEADER_SIZE 32
#define DATA_TYPE_OFFSET	4
#define NROWS_OFFSET		5
#define NCOLS_OFFSET		9
#define COMMENT_LEN_OFFSET     13

static int validate_data_type (int data_type)
{
   switch (data_type)
     {
      case 'A':
      case 'I':
      case 'J':
      case 'E':
      case 'F':
	break;
	
      default:
	JDMath_Error = JDMATH_INVALID_PARAMETER;
	JDMmsg_error ("Unsupported data type.");
	return -1;
     }
   
   return 0;
}

static JDMBData_File_Type *allocate_bdata_file_type (unsigned int len)
{
   JDMBData_File_Type *h;
   
   h = (JDMBData_File_Type *) _JDMmalloc (sizeof (JDMBData_File_Type), NULL);
   if (h == NULL) return NULL;
   
   memset ((char *) h, 0, sizeof (JDMBData_File_Type));

   if (NULL == (h->comment = _JDMmalloc (len + 1, NULL)))
     {
	_JDMfree ((char *)h);
	return NULL;
     }
   *h->comment = 0;
   return h;
}


JDMBData_File_Type *JDMbdata_open_file (char *file)
{
   FILE *fp;
   unsigned char header [BDATA_HEADER_SIZE];
   unsigned int comment_len;
   JDMBData_File_Type *bf;
   int32 rc_comment[3];

   bf = NULL;

   if (NULL == (fp = fopen (file, "rb")))
     {
	JDMath_Error = JDMATH_FILE_OPEN_ERROR;
	return NULL;
     }
   
   if (BDATA_HEADER_SIZE != fread (header, 1, BDATA_HEADER_SIZE, fp))
     goto read_error;
   
   if (memcmp ((char *) header, (char *)Magic_Bytes, 4))
     {
	JDMath_Error = JDMATH_CORRUPT_FILE_ERROR;
	JDMmsg_error ("Bad magic number");
	fclose (fp);
	return NULL;
     }
   
   JDMstr_read_int32 (rc_comment, 3, header + NROWS_OFFSET);
   comment_len = (unsigned int) rc_comment [2];
   bf = allocate_bdata_file_type (comment_len);
   if (bf == NULL)
     {
	fclose (fp);
	return NULL;
     }
   
   bf->nrows = (unsigned int) rc_comment[0];
   bf->ncols = (unsigned int) rc_comment[1];
   bf->data_type = (int) header [DATA_TYPE_OFFSET];

   if (-1 == validate_data_type (bf->data_type))
     {
	_JDMfree (bf->comment);
	_JDMfree ((char *)bf);
	fclose (fp);
	return NULL;
     }
   
   if (comment_len != fread (bf->comment, 1, comment_len, fp))
     goto read_error;
   
   bf->comment[comment_len] = 0;
   bf->fp = fp;
   bf->flags = JDMBDATA_READ_MODE;

   return bf;
   
   read_error:

   JDMath_Error = JDMATH_FILE_READ_ERROR;
   fclose (fp);
   if (bf != NULL)
     {
	if (bf->comment != NULL)
	  _JDMfree(bf->comment);
	_JDMfree ((char *)bf);
     }
   return NULL;
}


static int write_header (JDMBData_File_Type *b)
{
   unsigned char header[BDATA_HEADER_SIZE];
   int32 rc_comment[3];
   FILE *fp;

   fp = b->fp;

   memset ((char *) header, 0, BDATA_HEADER_SIZE);
   memcpy ((char *) header, (char *)Magic_Bytes, 4);

   header[DATA_TYPE_OFFSET] = (unsigned char) b->data_type;
   
   rc_comment[0] = (int32) b->nrows;
   rc_comment[1] = (int32) b->ncols;
   rc_comment[2] = (int32) b->comment_len;
   JDMstr_write_int32 (rc_comment, 3, header + NROWS_OFFSET);
   
   if (BDATA_HEADER_SIZE != fwrite ((char *) header, 1, BDATA_HEADER_SIZE, fp))
     return -1;

   return 0;
}


JDMBData_File_Type *JDMbdata_create_file (char *file, 
					  int data_type, 
					  unsigned int nrows, 
					  unsigned int ncols,
					  char *comment)
{
   FILE *fp;
   JDMBData_File_Type *bf;
   unsigned int comment_len;
   
   if (-1 == validate_data_type (data_type))
     {
	return NULL;
     }
   
   if (comment == NULL)
     comment_len = 0;
   else comment_len = strlen (comment);
   
   bf = allocate_bdata_file_type (comment_len);
   if (bf == NULL)
     return NULL;
   
   fp = fopen (file, "wb");
   if (fp == NULL)
     {
	JDMath_Error = JDMATH_FILE_OPEN_ERROR;
	return NULL;
     }

   bf->nrows = nrows;
   bf->ncols = ncols;
   bf->data_type = data_type;
   bf->fp = fp;
   bf->flags = JDMBDATA_WRITE_MODE;
   bf->comment_len = comment_len;

   if (-1 == write_header (bf))
     goto write_error;
   
   if (comment_len)
     {
	if (comment_len != fwrite (comment, 1, comment_len, fp))
	  goto write_error;
	strcpy (bf->comment, comment);
     }

   return bf;
   
   write_error:
   JDMath_Error = JDMATH_FILE_WRITE_ERROR;
   fclose (fp);
   if (bf != NULL)
     {
	if (bf->comment != NULL) _JDMfree (bf->comment);
	_JDMfree ((char *)bf);
     }
   
   return NULL;
}

int JDMbdata_close_file (JDMBData_File_Type *bf)
{
   int ret;
   FILE *fp;

   if (bf == NULL)
     return -1;
   
   if (NULL == (fp = bf->fp))
     return -1;

   ret = 0;

#ifndef SEEK_SET
# define SEEK_SET	0
#endif

   if (bf->flags & JDMBDATA_WRITE_MODE)
     {
	if (-1 == FSEEK (fp, 0, SEEK_SET))
	  {
	     JDMath_Error = JDMATH_FILE_WRITE_ERROR;
	     JDMmsg_error ("Seek error");
	     ret = -1;
	  }
	else if (-1 == write_header (bf))
	  {
	     JDMath_Error = JDMATH_FILE_WRITE_ERROR;
	     JDMmsg_error ("Write error");
	     ret = -1;
	  }
     }
   
   if (EOF == fclose (bf->fp))
     {
	JDMath_Error = JDMATH_FILE_CLOSE_ERROR;
	ret = -1;
     }

   bf->fp = NULL;
   if (bf->comment != NULL) _JDMfree (bf->comment);
   _JDMfree ((char *)bf);
   
   return ret;
}
