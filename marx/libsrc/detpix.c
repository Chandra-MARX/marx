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

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#include <jdmath.h>

#include "marx.h"
#include "_marx.h"

double Marx_Focal_Length = 10065.5;    /* mm */

/*{{{ ACIS-S */

static int acis_s_to_tiled (Marx_Detector_Type *det,
			    Marx_Detector_Geometry_Type *g,
			    int chip, unsigned int x, unsigned int y,
			    unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) det;
   (void) chip;

   xf = x + g->tdet_xoff;
   if (xf < 0.0)
     xf = 0.0;

   yf = y + g->tdet_yoff;
   if (yf < 0.0) 
     yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;
   return 0;
}


/* The following numbers come from table 4 in Rev 4.0 of the coordinate memo.
 */
#define REST_OF_STRUCT_FIELDS \
   0, 0, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 0, 0, 0, 0

static Marx_Detector_Geometry_Type ACIS_S_Geom [_MARX_NUM_ACIS_S_CHIPS] =
{
     {
	4,			       /* id */
	791, 1702,
	  {0.744, -81.051, -59.170},   /* x_ll */
	  {0.353, -56.478, -59.170},   /* x_lr */
	  {0.353, -56.478,  -34.590},   /* x_ur */
	  {0.744, -81.051, -34.590},   /* x_ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	5,			       /* id */
	1833, 1702,
	  {0.348, -56.047, -59.170},   /* x_ll */
	  {0.099, -31.473, -59.170},   /* x_lr */
	  {0.099, -31.473, -34.590},    /* x_ur */
	  {0.348, -56.047, -34.590},    /* x_ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	6,			       /* id */
	2875, 1702,
	  {0.096, -31.042, -59.170},   /* x_ll */
	  {-0.011, -6.466, -59.170},   /* x_lr */
	  {-0.011, -6.466, -34.590},    /* x_ur */
	  {0.096, -31.042, -34.590},    /* x_ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	7,			       /* id */
	3917, 1702,
	  {-0.011, -6.035, -59.170},   /* x_ll */
	  {0.024, 18.541, -59.170},    /* x_lr */
	  {0.024, 18.541, -34.590},     /* x_ur */
	  {-0.011, -6.035, -34.590},    /* x_ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	8,			       /* id */
	4959, 1702,
	  {0.026, 18.972, -59.170},    /* x_ll */
	  {0.208, 43.547, -59.170},    /* x_lr */
	  {0.208, 43.547, -34.590},     /* x_ur */
	  {0.026, 18.972, -34.590},     /* x_ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	9,			       /* id */
	6001, 1702,
	  {0.208, 43.978, -59.170},    /* x_ll */
	  {0.528, 68.552, -59.170},    /* x_lr */
	  {0.528, 68.552, -34.590},     /* x_ur */
	  {0.208, 43.978, -34.590},     /* x_ul */
	REST_OF_STRUCT_FIELDS
     }
};

static Marx_Detector_Type ACIS_S_Detector = 
{
   MARX_DETECTOR_ACIS_S,   	       /* detector_type */
   1024,			       /* num x pixels */
   1024,			       /* num y pixels */
   0.024,			       /* x_pixel_size */
   0.024,			       /* y_pixel_size */
   acis_s_to_tiled,
   4, 9,			       /* first/last_chip_id */
   _MARX_NUM_ACIS_S_CHIPS,	       /* num_chips */
   ACIS_S_Geom,
   {0.0, 0.0, 237.4},		       /* stt_lsi_offset (table 15 of rev 4.2) */
   {-0.684, 0.0,-190.133},		       /* stf_stt_offset  */
   /* Note: This is the ACIS-S AS1 aimpoint */
   
   /* Focal plane system info.  See section 7.2 of Rev 4.2 coord memo */
   "AXAF-FP-1.1",		       /* fp_system_name */
   NULL,			       /* pointer to fp system info */

   0				       /* is_initialized */
};


/*}}}*/

/*{{{ ACIS-I */

/* 01
 * 23
 * Note that the (0,0) pixel is not at the lower left of all these chips.
 * For chips 0 and 2, it is at the upper left with the xpixel running down
 * and the y pixel running to the right.
 */

static int
acis_i_to_tiled (Marx_Detector_Type *det,
		 Marx_Detector_Geometry_Type *g,
		 int chip, unsigned int x, unsigned int y,
		 unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) det;
   switch (chip)
     {
      case 0:
      case 2:
	xf = g->tdet_xoff + y;
	yf = g->tdet_yoff - x;
	break;

      case 1:
      case 3:
      default:
	xf = g->tdet_xoff - y;
	yf = g->tdet_yoff + x;
	break;
     }

   if (xf < 0.0) xf = 0.0;
   if (yf < 0.0) yf = 0.0;
   
   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;

   return 0;
}


/* The following numbers come from Table 17 of the Sept 9, 1996 version
 * of Jonathon McDowell's ASC Coordinates: Rev 3.2 Draft. The
 * tdet_xoff and tdet_yoff come from the Rev 4.0 version (table 5).
 * It appears that the geometry has not changed between versions 3.2
 * and 4.0.
 */
static Marx_Detector_Geometry_Type ACIS_I_Geom [_MARX_NUM_ACIS_I_CHIPS] =
{
     {
	0,
	3061, 5131,
	  { 2.361, -26.484,  23.088},   /* ll */
	  { 1.130, -26.546,  -1.458},   /* lr */
	  {-0.100,  -2.001,  -1.458},   /* ur */
	  { 1.130,  -1.939,  23.088},   /* ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	1,
	5131, 4107,
	  { 1.130,  23.086,  -1.458},   /* ll */
	  { 2.360,  23.024,  23.088},   /* lr */
	  { 1.130,  -1.521,  23.088},   /* ur */
	  {-0.100,  -1.459,  -1.458},   /* ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	2,
	3061, 4085,
	  { 1.130, -26.546,  -1.997},   /* ll */
	  { 2.361, -26.484, -26.543},   /* lr */
	  { 1.130,  -1.939, -26.543},   /* ur */
	  {-0.100,  -2.001,  -1.997},   /* ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	3,
	5131, 3061,
	  { 2.361,  23.024, -26.543},   /* ll */
	  { 1.131,  23.086,  -1.997},   /* lr */
	  {-0.100,  -1.459,  -1.997},   /* ur */
	  { 1.130,  -1.521, -26.543},   /* ul */
	REST_OF_STRUCT_FIELDS
     },
};


static Marx_Detector_Type ACIS_I_Detector = 
{
   MARX_DETECTOR_ACIS_I,   	       /* detector_type */
   1024,			       /* num x pixels */
   1024,			       /* num y pixels */
   0.024,			       /* x_pixel_size */
   0.024,			       /* y_pixel_size */
   acis_i_to_tiled,
   0, 3,			       /* first/last_chip_id */
   _MARX_NUM_ACIS_I_CHIPS,	       /* num_chips */
   ACIS_I_Geom,
   {-0.684, 0.750, 236.552},		       /* stt_lsi_offset (table 15 of rev 4.2) */
   {-0.782, 0.0, -233.592},
   /* Note: This is the AI2 ACIS-I aimpoint */
   
   /* Focal plane system info.  See section 7.2 of Rev 4.2 coord memo */
   "AXAF-FP-1.1",		       /* fp_system_name */
   NULL,			       /* pointer to fp system info */

   0				       /* is_initialized */
};


/*}}}*/

/*{{{ HRC-S */


static int
hrc_s_to_tiled (Marx_Detector_Type *det, 
		Marx_Detector_Geometry_Type *g,
		int chip, unsigned int x, unsigned int y,
		unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) chip;
   (void) det;

#if 1
   /* For AXAF-HRC-2.7S */
   /* Note: The coordinate document that I was given appears to be messed up
    * for this coord system.  At least, I think it is, so until I am proven
    * wrong, use this:
    */
# if 1
   yf = y;
   yf = yf + g->tdet_yoff;
   xf = x + g->tdet_xoff;
# else
   /* instead of this: */
   xf = x + g->tdet_xoff;
   yf = y + g->tdet_yoff;
# endif
#else
   /* For AXAF-HRC-2.6S */
   yf = y;
   xf = -yf + g->tdet_xoff;
   yf = x + g->tdet_yoff;
#endif

   if (xf < 0.0) xf = 0.0;
   if (yf < 0.0) yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;
   return 0;
}

#if 1
/* This table is patched up during runtime by the code in hrc_s_geom.c */
static Marx_Detector_Geometry_Type HRC_S_Geom [_MARX_NUM_HRC_S_CHIPS] =
{
     {
	0, 0, 0,
	  {2.234,-144.179, -9.050},    /* ll */
	  {0.100, -49.506, -9.050},    /* lr */
	  {0.100, -49.506,  9.100},    /* ur */
	  {2.234,-144.179,  9.100},     /* ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	1, 0, 0,
	  {0.100, -45.237, -9.050},    /* ll */
	  {0.100,  52.946, -9.050},    /* lr */
	  {0.100,  52.946,  9.000},    /* ur */
	  {0.100, -45.237,  9.000},     /* ul */
	REST_OF_STRUCT_FIELDS
     },
     {
	2, 0, 0,
	  {0.100,  57.704, -9.050},    /* ll */
	  {2.589, 151.810, -9.050},    /* lr */
	  {2.589, 151.810,  9.050},    /* ur */
	  {0.100,  57.704,  9.050},     /* ul */
	REST_OF_STRUCT_FIELDS
     }
};
#else
/* These data came from Mike Juda: */
/*{{{ Mike Juda's email */

/* From juda@head-cfa.harvard.edu Mon Oct 20 14:24:15 1997
 * To: davis
 * Cc: juda@head-cfa.harvard.edu, wise
 * Subject: HRC-S active area
 * Date: Mon, 20 Oct 1997 14:24:13 -0400
 * From: "Michael Juda" <juda@head-cfa.harvard.edu>
 * 
 * File: /data/juda1/juda/asc/hrc/data/hrc-s_geometry
 * 
 * Based on the post-XRCF C-K flat field data files: p197052710.rd,
 * p197053019.rd, and p197060110.rd I have found the following locations
 * for the boundaries of the HRC-S MCP "active" areas in pixel space (as
 * well as some UV/ion shield features).
 * 
 * Segment 1:
 * ----------
 * U lower limit ~600    (edge not illuminated - determined from background)
 * U upper limit 3496
 * 
 * V lower limit  1190
 * V CsI limit    1613
 * V upper limit 16250
 * 
 * U "strip" boundary 2660
 * 
 * Segment 2:
 * ----------
 * U lower limit ~600    (edge not illuminated - determined from background)
 * U upper limit 3488
 * 
 * V lower limit 16990
 * V upper limit 32261
 * 
 * U "T" boundary  2670
 * V "T" lower    22780
 * V "T" upper    27670
 * 
 * Segment 3:
 * ----------
 * U lower limit ~600    (edge not illuminated - determined from background)
 * U upper limit 3504
 * 
 * V lower limit 32925
 * V CsI limit   47650
 * V upper limit 48110
 * 
 * U "strip" boundary 2670
 * 
 * 
 * U-axis pixel size 0.00625 mm
 * V-axis pixel size 0.006429375 mm
 * 
 *  ^ +U and +Z
 *  |
 *  |
 *  |
 *  |--------> -V and +Y
 * 
 * 
 * Corner     (U,V)         (X,Y,Z)
 * ------     -----         -------
 * HRC-S1  ( 600, 1613)  (2.589, 151.810, -9.050)  Use CsI photocathode as end
 * HRC-S1  (3496, 1613)  (2.589, 151.810,  9.050)  Use CsI photocathode as end
 * HRC-S1  (3496,16250)  (0.100,  57.704,  9.050)
 * HRC-S1  ( 600,16250)  (0.100,  57.704, -9.050)
 * 
 * HRC-S2  ( 600,16990)  (0.100,  52.946, -9.050)
 * HRC-S2  (3488,16990)  (0.100,  52.946,  9.000)
 * HRC-S2  (3488,32261)  (0.100, -45.237,  9.000)
 * HRC-S2  ( 600,32261)  (0.100, -45.237, -9.050)
 * 
 * HRC-S3  ( 600,32925)  (0.100, -49.506, -9.050)
 * HRC-S3  (3504,32925)  (0.100, -49.506,  9.100)
 * HRC-S3  (3504,47650)  (2.234,-144.179,  9.100)  Use CsI photocathode as end
 * HRC-S3  ( 600,47650)  (2.234,-144.179, -9.050)  Use CsI photocathode as end
 * 
 * Optic axis (2048,25225)  (     ,   0.000,  0.000)
 * Det origin (   0,    0)  (     , 162.181,-12.800)
 */

/*}}}*/

static Marx_Detector_Geometry_Type HRC_S_Geom [_MARX_NUM_HRC_S_CHIPS] =
{
     {
	0, 0, 0,
	  {2.234,-144.179, -9.050},    /* ll */
	  {0.100, -49.506, -9.050},    /* lr */
	  {0.100, -49.506,  9.100},    /* ur */
	  {2.234,-144.179,  9.100}     /* ul */
     },
     {
	1, 0, 0,
	  {0.100, -45.237, -9.050},    /* ll */
	  {0.100,  52.946, -9.050},    /* lr */
	  {0.100,  52.946,  9.000},    /* ur */
	  {0.100, -45.237,  9.000}     /* ul */
     },
     {
	2, 0, 0,
	  {0.100,  57.704, -9.050},    /* ll */
	  {2.589, 151.810, -9.050},    /* lr */
	  {2.589, 151.810,  9.050},    /* ur */
	  {0.100,  57.704,  9.050}     /* ul */
     }
};
#endif

static Marx_Detector_Type HRC_S_Detector = 
{
   MARX_DETECTOR_HRC_S,   	       /* detector_type */
   4096,			       /* num x pixels */
   16384,			       /* num y pixels */
   6.429e-3,			       /* x_pixel_size */
   6.429e-3,			       /* y_pixel_size */
   hrc_s_to_tiled,
   1, 3,			       /* first/last_chip_id */
   _MARX_NUM_HRC_S_CHIPS,	       /* num_chips */
   HRC_S_Geom,
   {-1.533, 1.530, -251.437},
   {-0.1, 0.0, 250.1},		       /* stf_stt_offset (table 16 of rev 4.2) */
   /* This is the HS1 aimpoint */
   
   /* Focal plane system info.  See section 7.2 of Rev 4.2 coord memo */
   "AXAF-FP-2.3",		       /* fp_system_name */
   NULL,			       /* pointer to fp system info */
   0				       /* is_initialized */
};


/*}}}*/

/*{{{ HRC-I */

/* The following numbers come from Table 20 of the Sept 9, 1996 version
 * of Jonathon McDowell's ASC Coordinates: Rev 3.2 Draft.
 * 
 * Apparantly, the HRC-I geometry is not very well known and the bore-sight
 * will not be known until it it in flight. This comes from SM.
 * The active area is roughly 92.6 mm on a side instead of the numbers that
 * JMcD has in his memo.  So, until we get a more definitive answer, I am just
 * going to fudge this.  Sigh.
 */
#define YA_HRC_HACK (92.6/105.343)
static Marx_Detector_Geometry_Type HRC_I_Geom [_MARX_NUM_HRC_I_CHIPS] =
{
     {
	0,
	0, 0,
	  {0.000,   0.000,  74.489*YA_HRC_HACK},   /* ll */
	  {0.000,  74.489*YA_HRC_HACK,   0.000},   /* lr */
	  {0.000,   0.000, -74.489*YA_HRC_HACK},   /* ur */
	  {0.000, -74.489*YA_HRC_HACK,   0.000},   /* ul */
	REST_OF_STRUCT_FIELDS
     }
};


static int
hrc_i_to_tiled (Marx_Detector_Type *det, 
		Marx_Detector_Geometry_Type *g,
		int chip, unsigned int x, unsigned int y,
		unsigned int *xp, unsigned int *yp)
{
   float xf, yf;

   (void) chip;
   (void) det;
   (void) g;

   xf = (x + g->tdet_xoff);
   yf = (y + g->tdet_yoff);

   if (xf < 0.0) xf = 0.0;
   if (yf < 0.0) yf = 0.0;

   *xp = (unsigned int) xf;
   *yp = (unsigned int) yf;
   return 0;
}


static Marx_Detector_Type HRC_I_Detector =
{
   MARX_DETECTOR_HRC_I,   	       /* detector_type */
   16384,			       /* num x pixels */
   16384,			       /* num y pixels */
   6.429e-3,			       /* x_pixel_size */
   6.429e-3,			       /* y_pixel_size */
   hrc_i_to_tiled,
   0, 0,			       /* first/last_chip_id */
   _MARX_NUM_HRC_I_CHIPS,	       /* num_chips */
   HRC_I_Geom,
   {0.15, 0.0, -126.6},		       /* stt_lsi_offset (table 15 of Rev 4.2) */
   {-0.15, 0.0, 126.6},		       /* stf_stt_offset (table 16 of Rev 4.2) */
   /* Note: This is for the HI1 aimpoint */
   
   /* Focal plane system info.  See section 7.2 of Rev 4.2 coord memo */
   "AXAF-FP-2.1",		       /* fp_system_name */
   NULL,			       /* pointer to fp system info */
   0				       /* is_initialized */
};
/*}}}*/


static Marx_FP_Coord_Type Focal_Plane_Coord_Info [] = 
{
     {
	"AXAF-FP-1.1",		       /* fp_system_name */
#if 0
	0.49190 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
#else
	0.492 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
#endif
	4096.5,			       /* fp x pixel of on axis ray */
	4096.5			       /* fp y pixel of on axis ray */
     },
     {
	"AXAF-FP-2.1",		       /* fp_system_name */
#if 0
	0.132 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
#else
	0.13175 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
#endif
	16384.5,		       /* fp x pixel of on axis ray */
	16384.5			       /* fp y pixel of on axis ray */
     },
     {
	"AXAF-FP-2.3",		       /* fp_system_name */
#if 0
	0.132 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
#else
	0.1318 * (PI/(3600.0*180.0)),   /* fp_delta_s0 (radians) */
#endif
	32768.5,		       /* fp x pixel of on axis ray */
	32768.5			       /* fp y pixel of on axis ray */
     },
     {
	NULL, 0, 0, 0
     }
};


Marx_FP_Coord_Type *marx_get_fp_system_info (char *name)
{
   Marx_FP_Coord_Type *fp;
   
   fp = Focal_Plane_Coord_Info;
   while (fp->fp_system_name != NULL)
     {
	if (0 == strcmp (name, fp->fp_system_name))
	  return fp;
	fp++;
     }
   
   marx_error ("Focal plane system %s not implemented", name);
   return NULL;
}



static void 
compute_detector_normals (Marx_Detector_Type *det)
{
   Marx_Detector_Geometry_Type *d, *dmax;
   JDMVector_Type ofs;

   d = det->geom;
   dmax = d + det->num_chips;
   
   ofs = JDMv_sum (det->stf_stt_offset, det->stt_lsi_offset);

   while (d < dmax)
     {
	d->x_lr = JDMv_sum (d->x_lr, ofs);
	d->x_ll = JDMv_sum (d->x_ll, ofs);
	d->x_ur = JDMv_sum (d->x_ur, ofs);
	d->x_ul = JDMv_sum (d->x_ul, ofs);
	
	d->xhat = JDMv_diff (d->x_lr, d->x_ll);
	d->yhat = JDMv_diff (d->x_ul, d->x_ll);
	d->xlen = JDMv_length (d->xhat);
	d->ylen = JDMv_length (d->yhat);
	JDMv_normalize (&(d->xhat));
	JDMv_normalize (&(d->yhat));
	d->normal = JDMv_cross_prod (d->xhat, d->yhat);

	/* This should be handled elsewhere... */
	d->x_pixel_size = det->x_pixel_size;
	d->y_pixel_size = det->y_pixel_size;

	d++;
     }
}

/*}}}*/

/* This is necessary because the corners defined in the pixlib parameter
 * file do not produce a detector whose size is 1024 * 0.024.
 */
static void tweak_acis_pixel_size (Marx_Detector_Type *det)
{
   Marx_Detector_Geometry_Type *d, *dmax;

   d = det->geom;
   dmax = d + det->num_chips;
   
   while (d < dmax)
     {
	d->x_pixel_size = d->xlen / 1024.0;
	d->y_pixel_size = d->ylen / 1024.0;
	
	d++;
     }
}

typedef struct
{
   char *name;
   int id;
}
Det_Name_Type;

static Det_Name_Type Detectors [] = 
{
   {"HRC-S", MARX_DETECTOR_HRC_S},
   {"HRC-I", MARX_DETECTOR_HRC_I},
   {"ACIS-S", MARX_DETECTOR_ACIS_S},
   {"ACIS-I", MARX_DETECTOR_ACIS_I},
   {NULL, 0}
};

static int patch_aimpoint (Marx_Detector_Type *d, char *aimpt_parm, char *origin_parm)
{
#if 1
   (void) aimpt_parm; (void) origin_parm;
   return _marx_caldb_patch_aimpoint (d);
#else
   char *file;
   Param_File_Type *pf;
   int status;

   if (NULL == (file = marx_make_data_file_name ("pixlib/pix_sim_table_flight.par")))
     return -1;

   /* marx_message ("\t%s\n", file); */
   pf = pf_open_parameter_file (file, "rQ");
   if (pf == NULL)
     {
	marx_error ("Unable to open %s", file);
	marx_free (file);
	return -1;
     }

   status = 0;
   if ((-1 == _marx_get_vector_parm (pf, origin_parm, &d->stt_lsi_offset))
       || (-1 == _marx_get_vector_parm (pf, aimpt_parm, &d->stf_stt_offset)))
     status = -1;

   (void) pf_close_parameter_file (pf);
   marx_free (file);
   return status;
#endif
}


Marx_Detector_Type *
marx_get_detector_info (char *detname) /*{{{*/
{
   Marx_Detector_Type *d;
   Det_Name_Type *dn;
   int det;

   if (detname == NULL)
     return NULL;

   dn = Detectors;
   det = -1;
   while (dn->name != NULL)
     {
	if (0 == strcmp (detname, dn->name))
	  {
	     det = dn->id;
	     break;
	  }
	dn++;
     }

   switch (det)
     {
      default:
	marx_error ("Detector %s is unknown", detname);
	return NULL;

      case MARX_DETECTOR_ACIS_S:
	d = &ACIS_S_Detector;
	if (d->is_initialized) return d;
	if (-1 == _marx_patch_acis_s_geom (d))
	  return NULL;
	if (-1 == patch_aimpoint (d, "aimpt-AS1", "origin-ACIS"))
	  return NULL;
	break;
	
      case MARX_DETECTOR_ACIS_I:
	d = &ACIS_I_Detector;
	if (d->is_initialized) return d;
	if (-1 == _marx_patch_acis_i_geom (d))
	  return NULL;
	if (-1 == patch_aimpoint (d, "aimpt-AI2", "origin-ACIS"))
	  return NULL;
	break;
	
      case MARX_DETECTOR_HRC_S:
	d = &HRC_S_Detector;
	if (d->is_initialized) return d;
	if (-1 == _marx_patch_hrc_s_geom (d))
	  return NULL;
	if (-1 == patch_aimpoint (d, "aimpt-HS1", "origin-HRCS"))
	  return NULL;
	break;

      case MARX_DETECTOR_HRC_I:
	d = &HRC_I_Detector;
	if (d->is_initialized) return d;

	if (-1 == _marx_patch_hrc_i_geom (d))
	  return NULL;
	if (-1 == patch_aimpoint (d, "aimpt-HI1", "origin-HRCI"))
	  return NULL;

	break;
     }
   compute_detector_normals (d);
   
   switch (det)
     {
      case MARX_DETECTOR_ACIS_I:
      case MARX_DETECTOR_ACIS_S:
	tweak_acis_pixel_size (d);
	break;
     }

   d->fp_coord_info = marx_get_fp_system_info (d->fp_system_name);

   d->is_initialized = 1;
   return d;
}

/*}}}*/

int 
marx_compute_tiled_pixel (Marx_Detector_Type *d, /*{{{*/
			      int chip, unsigned int x, unsigned int y,
			      unsigned int *xp, unsigned int *yp)
{
   Marx_Detector_Geometry_Type *g;
   
   if (d == NULL) 
     return -1;
   
   if ((chip < d->first_chip_id) || (chip > d->last_chip_id))
     {
	marx_error ("chip = %d is not appropriate for this detector", chip);
	return -1;
     }

   g = d->geom + (chip - d->first_chip_id);
   
   return (*d->tiled_pixel_map_fun) (d, g, chip, x, y, xp, yp);
}

/*}}}*/

int marx_mnc_to_fpc (Marx_FP_Coord_Type *fc, 
		     JDMVector_Type *mnc,
		     double *xp, double *yp)
{
   double mnc_x;
   double factor;

   if (fc == NULL)
     return -1;
   
   mnc_x = mnc->x;
   if (mnc_x == 0.0)
     {
	marx_error ("marx_mnc_to_fpc: mnc.x is 0");
	return -1;
     }
   
   factor = 1.0 / (fc->fp_delta_s0 * mnc_x);
   
   /* Consider a ray with a negative MNC->x coordinate, but positive y and
    * z coordinates.  We want this to map to a positive FP_X and a negative
    * FP_Y value.
    */
   *xp = fc->fp_x0 - factor * mnc->y;
   *yp = fc->fp_y0 + factor * mnc->z;
   
   return 0;
}

int _marx_get_vector_parm (Param_File_Type *pf, char *parm, JDMVector_Type *v)
{
   char buf[256];
   JDMVector_Type vv;

   if (-1 == pf_get_string (pf, parm, buf, sizeof (buf)))
     {
	marx_error ("Error getting parameter value for %s", parm);
	return -1;
     }
   
   if (3 != sscanf (buf, "(%lf%*[ ,]%lf%*[ ,]%lf)", &vv.x, &vv.y, &vv.z))
     {
	marx_error ("Parameter %s does not have the proper format", parm);
	return -1;
     }
   
   *v = vv;
   return 0;
}

