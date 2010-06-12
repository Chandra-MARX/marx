/* -*- mode: C; mode: fold; -*- */
/*
    This file is part of MARX

    Copyright (C) 2002-2010 Massachusetts Institute of Technology

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

/*{{{ Include Files */
#include <stdio.h> 
#include <string.h>


#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>

#include "marx.h"
#include "_marx.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLFREE free
#endif

/*}}}*/

typedef struct /*{{{*/
{
   int32 num_rows;
   int32 num_cols;
   FILE *fp;
} 
/*}}}*/
Write_File_Type;

#define NBITS_LONG (8 * sizeof (long))
static Write_File_Type File_Pointers [NBITS_LONG];

static unsigned char MagicChars [4] = /*{{{*/
{
   0x83, 0x13, 0x89, 0x8D
};

/*}}}*/
#define HEADER_SIZE 32
/* Num element offset MUST follow 4 magic bytes + 16 byte column name. */
#define NUM_ROWS_OFFSET 20
#define NUM_COLS_OFFSET 24
#define RESERVED_OFFSET (NUM_COLS_OFFSET + 4)
#define HEADER_RESERVED (HEADER_SIZE - RESERVED_OFFSET)

int marx_close_write_dump_file (FILE *fp, unsigned long n)
{
   char *errfmt = "marx_close_write_dump_file: %s (errno = %d)";
   int ret = 0;
   int32 i32;
   
   if (fp == NULL)
     return -1;
   
   i32 = (int32) n;

   errno = 0;
#ifndef SEEK_SET
# define SEEK_SET	0
#endif
   if (-1 == FSEEK (fp, NUM_ROWS_OFFSET, SEEK_SET))
     {
	marx_error (errfmt, "seek error", errno);
	ret = -1;
     }
   else if (1 != JDMwrite_int32 (&i32, 1, fp))
     {
	marx_error (errfmt, "write error", errno);
	ret = -1;
     }
   
   if (EOF == fclose (fp))
     {
	marx_error (errfmt, "write error", errno);
	ret = -1;
     }
   
   if (ret == -1)
     {
#ifdef ENOSPC
	marx_error ("No space left on device.");
#endif
     }
   
   return ret;
}

static int close_write_files (int do_num) /*{{{*/
{
   unsigned int i;
   FILE *fp;
   int ret = 0;
   
   for (i = 0; i < NBITS_LONG; i++)
     {
	if (NULL != (fp = File_Pointers[i].fp))
	  {
	     File_Pointers[i].fp = NULL;
	     
	     if (do_num)
	       {
		  unsigned long num_rows = File_Pointers[i].num_rows;
		  
		  if (-1 == marx_close_write_dump_file (fp, num_rows))
		    ret = -1;
	       }
	     else if (EOF == fclose (fp))
	       ret = -1;

	     /* File_Pointers[i].num_rows = 0; */
	  }
     }
   return ret;
}

/*}}}*/

FILE *marx_create_write_dump_file (char *filename, Marx_Dump_File_Type *dft)
{
   char header[HEADER_SIZE];
   char *h, *hmax, *p;
   FILE *fp;
   unsigned int len;
   
   dft->fp = NULL;
   
   if (NULL == (fp = fopen (filename, "w+b")))
     {
	marx_error ("Unable to open %s for writing.", filename);
#ifdef EMFILE
	if (errno == EMFILE)
	  marx_error (" (Too many open files.)");
#endif
	return NULL;
     }
   
   memcpy (header, (char *) MagicChars, 4);
   header[4] = dft->type;
   
   h = header + 5;
   hmax = header + HEADER_SIZE;

   p = dft->colname;
   while (h < hmax)
     {
	if (*p != 0) *h++ = *p++;
	else *h++ = 0;
     }
   
   if ((NUM_ROWS_OFFSET != fwrite (header, 1, NUM_ROWS_OFFSET, fp))
       || (1 != JDMwrite_int32 (&dft->num_rows, 1, fp))
       || (1 != JDMwrite_int32 (&dft->num_cols, 1, fp)))
     {
	fclose (fp);
	marx_error ("Write to %s failed.", filename);
	return NULL;
     }
   
   /* Now write remaining bytes. */
   h = header + RESERVED_OFFSET;
   len = HEADER_RESERVED;
   
   if (len != fwrite (h, 1, len, fp))
     {
	marx_error ("Write to %s failed.", filename);
	fclose (fp);
	return NULL;
     }
   
   dft->fp = fp;
   return fp;
}


static int open_files (char *dir, unsigned long write_mask, int new_file) /*{{{*/
{
   char *filename = NULL;
   FILE *fp;
   char *open_mode;
   char *col_name;
   char data_type;
   unsigned int nbits;
   unsigned long mask;
   
   open_mode = "r+b";	       
   /* Note: "ab" will not work because of the system 5 behavior of writes 
    * to a file opened this way: all writes will go to end of file 
    * regardless of seek!!!
    */
   
   nbits = NBITS_LONG;
   mask = 1;
   while (nbits)
     {
	nbits--;
	switch (mask & write_mask)
	  {
	   case MARX_PI_OK:
	     filename = "b_energy.dat";
	     data_type = 'E';
	     col_name = "B_ENERGY";
	     break;

	   case MARX_ENERGY_OK:
	     filename = "energy.dat";
	     data_type = 'E';
	     col_name = "ENERGY";
	     break;
	   case MARX_TIME_OK:
	     filename = "time.dat";
	     data_type = 'E';
	     col_name = "TIME";
	     break;
	   case MARX_X1_VECTOR_OK:
	     filename = "xpos.dat";
	     data_type = 'E';
	     col_name = "XPOS";
	     break;
	   case MARX_X2_VECTOR_OK:
	     filename = "ypos.dat";
	     data_type = 'E';
	     col_name = "YPOS";
	     break;
	   case MARX_X3_VECTOR_OK:
	     filename = "zpos.dat";
	     data_type = 'E';
	     col_name = "ZPOS";
	     break;
	   case MARX_P1_VECTOR_OK:
	     filename = "xcos.dat";
	     data_type = 'E';
	     col_name = "COSX";
	     break;
	   case MARX_P2_VECTOR_OK:
	     filename = "ycos.dat";
	     data_type = 'E';
	     col_name = "COSY";
	     break;
	   case MARX_P3_VECTOR_OK:
	     filename = "zcos.dat";
	     data_type = 'E';
	     col_name = "COSZ";
	     break;
	   case MARX_PULSEHEIGHT_OK:
	     filename = "pha.dat";
	     data_type = 'I';
	     col_name = "PHA";
	     break;
	   case MARX_CCD_NUM_OK:
	     filename = "detector.dat";
	     data_type = 'A';
	     col_name = "CCDID";
	     break;
	   case MARX_Y_PIXEL_OK:
	     filename = "xpixel.dat";
	     data_type = 'E';
	     col_name = "CHIPX";
	     break;
	   case MARX_Z_PIXEL_OK:
	     filename = "ypixel.dat";
	     data_type = 'E';
	     col_name = "CHIPY";
	     break;
	   case MARX_U_PIXEL_OK:
	     filename = "hrc_u.dat";
	     data_type = 'E';
	     col_name = "U";
	     break;
	   case MARX_V_PIXEL_OK:
	     filename = "hrc_v.dat";
	     data_type = 'E';
	     col_name = "V";
	     break;
	   case MARX_MIRROR_SHELL_OK:
	     filename = "mirror.dat";
	     data_type = 'I';
	     col_name = "MIRROR";
	     break;
	   case MARX_HRC_REGION_OK:
	     filename = "hrcregion.dat";
	     data_type = 'A';
	     col_name = "HRCREGION";
	     break;
	   case MARX_ORDER_OK:
	     filename = "order.dat";
	     data_type = 'A';
	     col_name = "ORDER";
	     break;
	   case MARX_SUPPORT_ORDER1_OK:
	     filename = "ofine.dat";
	     data_type = 'A';
	     col_name = "FINE";
	     break;
	   case MARX_SUPPORT_ORDER2_OK:
	     filename = "ocoarse1.dat";
	     data_type = 'A';
	     col_name = "COARSE1";
	     break;
	   case MARX_SUPPORT_ORDER3_OK:
	     filename = "ocoarse2.dat";
	     data_type = 'A';
	     col_name = "COARSE2";
	     break;
	   case MARX_SUPPORT_ORDER4_OK:
	     filename = "ocoarse3.dat";
	     data_type = 'A';
	     col_name = "COARSE3";
	     break;
#if MARX_HAS_DITHER
	   case MARX_SKY_RA_OK:
	     filename = "sky_ra.dat";
	     data_type = 'E';
	     col_name = "RA";
	     break;
	   case MARX_SKY_DEC_OK:
	     filename = "sky_dec.dat";
	     data_type = 'E';
	     col_name = "DEC";
	     break;
#endif
	   default:
	     filename = NULL;
	     data_type = 0;
	     col_name = NULL;
	  }
	
	mask = mask << 1;

	if (filename == NULL)
	  {
	     continue;
	  }
	
	     
	if (NULL == (filename = marx_dircat (dir, filename)))
	  goto return_error;

	if (new_file)
	  {
	     Marx_Dump_File_Type dft;
	     
	     memset ((char *) &dft, 0, sizeof (Marx_Dump_File_Type));
	     dft.type = data_type;
	     strcpy (dft.colname, col_name);
	     
	     fp = marx_create_write_dump_file (filename, &dft);
	  }
	else fp = fopen (filename, open_mode);
	
	if (fp == NULL)
	  {
	     marx_error ("Unable to open %s.", filename);
#ifdef EMFILE
	if (errno == EMFILE)
	  marx_error (" (Too many open files.)");
#endif

	     goto return_error;
	  }
	
	File_Pointers[nbits].fp = fp;
	
	if (new_file == 0)
	  {
	     /* See comment above regarding r+t mode. */
	     if (-1 == FSEEK (fp, 0, SEEK_END))
	       {
		  marx_error ("seek error.");
		  goto return_error;
	       }
	  }
	
	SLFREE (filename);
	filename = NULL;
     }
   
   return 0;
   
   return_error:
   if (filename != NULL) SLFREE (filename);
   (void) close_write_files (0);
   return -1;
}

/*}}}*/

int marx_write_photons (Marx_Photon_Type *p, unsigned long write_mask, /*{{{*/
			char *dir, int open_mode, double total_time)
{
   unsigned long mask;
   unsigned int i, imax, nbits;
   Marx_Photon_Attr_Type *attr, *at;
   
#if MARX_HAS_DITHER
   if (_Marx_Dither_Mode != _MARX_DITHER_MODE_NONE)
     write_mask &= (p->history | MARX_SKY_RA_OK | MARX_SKY_DEC_OK);
   else
#endif
     write_mask &= p->history;

   if (-1 == open_files (dir, write_mask, open_mode))
     return -1;
	
   imax = p->n_photons;
   attr = p->attributes;

   for (i = 0; i < imax; i++)
     {
#if MARX_HAS_DITHER
	double ra = 0, dec = 0;
#endif
	at = attr + i;
	if (at->flags & BAD_PHOTON_MASK)
	  {
#define OUTPUT_ABSORBED	0
#if OUTPUT_ABSORBED
	     if (0 == (at->flags & OUTPUT_ABSORBED))
#endif
	       continue;
	  }
#if OUTPUT_ABSORBED
	else continue;
#endif

#if MARX_HAS_DITHER
	if (write_mask & (MARX_SKY_DEC_OK|MARX_SKY_RA_OK))
	  _marx_ray_to_sky_ra_dec (at, &ra, &dec);
#endif

	nbits = NBITS_LONG;
	mask = 1;
	while (nbits)
	  {
	     FILE *fp;
	     unsigned int one;
	     float32 float_value;
	     int16 small_int;
	     SIGNED_CHAR tiny_int;

	     nbits--;
	     fp = File_Pointers[nbits].fp;

	     switch (mask & write_mask)
	       {
		case MARX_PI_OK:
		  float_value = at->pi;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_ENERGY_OK:
		  float_value = at->energy;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;
		  
		case MARX_TIME_OK: 
		  float_value = at->arrival_time + total_time;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_X1_VECTOR_OK: 
		  float_value = at->x.x;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_X2_VECTOR_OK:
		  float_value = at->x.y;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_X3_VECTOR_OK:
		  float_value = at->x.z;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_P1_VECTOR_OK:
		  float_value = at->p.x;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_P2_VECTOR_OK:
		  float_value = at->p.y;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;
		  
		case MARX_P3_VECTOR_OK:
		  float_value = at->p.z;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_PULSEHEIGHT_OK:
		  small_int = at->pulse_height;
		  one = JDMwrite_int16 (&small_int, 1, fp);
		  break;

		case MARX_CCD_NUM_OK:
		  tiny_int = at->ccd_num;
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;
		  
		case MARX_MIRROR_SHELL_OK:
		  tiny_int = at->mirror_shell;
		  one = JDMwrite_int16 (&small_int, 1, fp);
		  break;

		case MARX_HRC_REGION_OK:
		  tiny_int = at->detector_region;
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;

		case MARX_ORDER_OK:
		  tiny_int = at->order;
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;
		  
		case MARX_SUPPORT_ORDER4_OK:
		  tiny_int = at->support_orders[3];
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;

		case MARX_SUPPORT_ORDER3_OK:
		  tiny_int = at->support_orders[2];
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;

		case MARX_SUPPORT_ORDER2_OK:
		  tiny_int = at->support_orders[1];
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;

		case MARX_SUPPORT_ORDER1_OK:
		  tiny_int = at->support_orders[0];
		  one = fwrite (&tiny_int, 1, 1, fp);
		  break;

		case MARX_Y_PIXEL_OK:
		  float_value = at->y_pixel;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_Z_PIXEL_OK:
		  float_value = at->z_pixel;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_U_PIXEL_OK:
		  float_value = at->u_pixel;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;

		case MARX_V_PIXEL_OK:
		  float_value = at->v_pixel;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;


#if MARX_HAS_DITHER
		case MARX_SKY_RA_OK:
		  float_value = ra;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;
		  
		case MARX_SKY_DEC_OK:
		  float_value = dec;
		  one = JDMwrite_float32 (&float_value, 1, fp);
		  break;
#endif
		default:
		  one = 1;
	       }
	     
	     mask = mask << 1;

	     if (one != 1)
	       {
		  marx_error ("write error.");
		  (void) close_write_files (0);
		  return -1;
	       }
	     
	     File_Pointers[nbits].num_rows += 1;
	  }
     }
   
   if (-1 == close_write_files (1))
     {
	marx_error ("write error.");
	return -1;
     }
   return 0;
}

/*}}}*/

static void dump_usage (void) /*{{{*/
{
   fprintf (stderr, "Usage: marx --dump filename [filename...]\n");
}

/*}}}*/

static Marx_Dump_File_Type *Data_File_Type_Root;

static int chain_read_data_file_type (Marx_Dump_File_Type *d) /*{{{*/
{
   Marx_Dump_File_Type *dlast;
   
   if (Data_File_Type_Root == NULL)
     {
	Data_File_Type_Root = d;
	return 0;
     }
   
   dlast = Data_File_Type_Root;
   while (dlast->next != NULL) dlast = dlast->next;
   dlast->next = d;
   return 0;
}

/*}}}*/

static int close_read_data_files (void) /*{{{*/
{
   Marx_Dump_File_Type *d, *dnext;
   int ret = 0;
   
   d = Data_File_Type_Root;
   
   while (d != NULL)
     {
	dnext = d->next;
	marx_close_read_dump_file (d);
	d = dnext;
     }
   
   Data_File_Type_Root = NULL;
   return ret;
}

/*}}}*/

static int dump_data_files (void) /*{{{*/
{
   int not_done;
   Marx_Dump_File_Type *dft;
   int32 num_rows;
   unsigned int num_cols;
   
   dft = Data_File_Type_Root;
   
   fputc ('#', stdout);
   while (dft != NULL)
     {
	if (0 == (num_cols = dft->num_cols))
	  num_cols = 1;

	while (num_cols--)
	  fprintf (stdout, "%16s", dft->colname);

	dft = dft->next;
     }
   fputc ('\n', stdout);
   
   num_rows = 0;
   while (1)
     {
	not_done = 1;

	dft = Data_File_Type_Root;

	while (dft != NULL)
	  {
	     FILE *fp = dft->fp;
	     float64 dval;
	     float32 fval;
	     int32 ival;
	     int16 sval;
	     SIGNED_CHAR cval;
	     
	     num_cols = dft->num_cols;
	     if (num_cols == 0)
	       num_cols = 1;

	     switch (dft->type)
	       {
		default:
		  marx_error ("Unknown data type: '%c'", dft->type);
		  (void) close_read_data_files ();
		  return -1;
		  
		case 'A':
		  while (num_cols && not_done)
		    {
		       if (1 == fread (&cval, 1, 1, fp))
			 fprintf (stdout, "%16d", cval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;

		case 'I':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_int16 (&sval, 1, fp))
			 fprintf (stdout, "%16d", sval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;

		case 'J':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_int32 (&ival, 1, fp))
			 fprintf (stdout, "%16ld", (long) ival);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;
		  
		case 'E':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_float32 (&fval, 1, fp))
			 fprintf (stdout, "%16e", fval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;
		  
		case 'D':
		  while (num_cols && not_done)
		    {
		       if (1 == JDMread_float64 (&dval, 1, fp))
			 fprintf (stdout, "%16e", dval);
		       else not_done = 0;
		       num_cols--;
		    }
		  break;
	       }
	     
	     if (not_done == 0)
	       {
		  if (num_rows != dft->num_rows)
		    {
		       marx_error ("Read error. (%d/%d rows)",
				  num_rows, dft->num_rows);
		       close_read_data_files ();
		       return -1;
		    }
		  break;
	       }
	     dft = dft->next;
	  }
	num_rows++;
	putc ('\n', stdout);
	if (not_done == 0) 
	  break;
     }
   
   close_read_data_files ();
   return 0;
}

/*}}}*/

Marx_Dump_File_Type *marx_open_read_dump_file (char *file) /*{{{*/
{
   FILE *fp;   
   char header[HEADER_SIZE];
   char reserved[HEADER_RESERVED];
   Marx_Dump_File_Type *dft;

   dft = (Marx_Dump_File_Type *) SLMALLOC (sizeof (Marx_Dump_File_Type));
   if (dft == NULL)
     {
	marx_error ("Memory allocation failure.");
	return NULL;
     }
   
   memset ((char *) dft, 0, sizeof (Marx_Dump_File_Type));
   
   if (NULL == (fp = fopen (file, "rb")))
     {
	marx_error ("Unable to open %s.", file);
#ifdef EMFILE
	if (errno == EMFILE)
	  marx_error (" (Too many open files.)");
#endif
	SLFREE (dft);
	return NULL;
     }
	
   if ((NUM_ROWS_OFFSET != fread (header, 1, NUM_ROWS_OFFSET, fp))
       || (1 != JDMread_int32 (&dft->num_rows, 1, fp))
       || (1 != JDMread_int32 (&dft->num_cols, 1, fp))
       || (HEADER_RESERVED != fread(reserved, 1, HEADER_RESERVED, fp)))
     {
	marx_error ("Error reading header information for %s.", file);
	fclose (fp);
	SLFREE (dft);
	return NULL;
     }
	
   if (memcmp (header, (char *)MagicChars, 4))
     {
	marx_error ("File %s does not appear to be a marx output file.  Bad Magic Number.",
		   file);
	fclose (fp);
	SLFREE (dft);
	return NULL;
     }
	
   dft->type = *(header + 4);
   strncpy (dft->colname, header + 5, 15);
   dft->colname[15] = 0;
   dft->fp = fp;
   
   return dft;
}

/*}}}*/

/* Some programs may pass NULL to this.  No error please. */
void marx_close_read_dump_file (Marx_Dump_File_Type *dft) /*{{{*/
{
   if (dft == NULL)
     return;
   
   if (dft->fp != NULL)
     (void) fclose (dft->fp);
   
   SLFREE (dft);
}

/*}}}*/

int marx_dump (int argc, char **argv) /*{{{*/
{
   char *file;
   int i;
   Marx_Dump_File_Type *dft;
   
   if (argc < 3)
     {
	dump_usage ();
	return -1;
     }
   
   for (i = 2; i < argc; i++)
     {
	file = argv[i];
	
	dft = marx_open_read_dump_file (file);
	
	if ((dft == NULL)
	    || (-1 == chain_read_data_file_type (dft)))

	  {
	     marx_error ("Unable to open %s.", file);
#ifdef EMFILE
	     if (errno == EMFILE)
	       marx_error (" (Too many open files.)");
#endif
	     close_read_data_files ();
	     if (dft != NULL) fclose (dft->fp);
	     return -1;
	  }
     }
   
   return dump_data_files ();
}

/*}}}*/
