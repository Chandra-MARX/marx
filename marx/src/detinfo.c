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

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <marx.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

static char *Pgm_Name = "detinfo";

static void usage (void)
{
   fprintf (stderr, "Usage: %s [--sky] ACIS-S|ACIS-I|HRC-S|HRC-I\n", Pgm_Name);
   exit (1);
}

static void dump_physical_geometry (Marx_Detector_Type *);
static void dump_sky_geometry (Marx_Detector_Type *, char *);

int main (int argc, char **argv)
{
   char *name = NULL;
   Marx_Detector_Type *det;
   int sky = 0;

   switch (argc)
     {
      case 2:
	name = argv[1];
	break;

      case 3:
	if (0 == strcmp ("--sky", argv[1]))
	  {
	     name = argv[2];
	     sky = 1;
	     break;
	  }
	/* drop */
      default:
	usage ();
     }

   if (NULL == (det = marx_get_detector_info (name)))
     usage ();
   
   fprintf (stdout, " Detector Name: %s\n", name);
   fprintf (stdout, "  Num X Pixels: %u (per chip)\n", det->num_x_pixels);
   fprintf (stdout, "  Num Y Pixels: %u (per chip)\n", det->num_y_pixels);
   fprintf (stdout, "  X Pixel Size: % 10.4e (mm)\n", det->x_pixel_size);
   fprintf (stdout, "  Y Pixel Size: % 10.4e (mm)\n", det->y_pixel_size);
   fprintf (stdout, " First Chip id: %d\n", det->first_chip_id);
   fprintf (stdout, "  Last Chip id: %d\n", det->last_chip_id);
   fprintf (stdout, "     Num Chips: %u\n", det->num_chips);
   fprintf (stdout, "STT-LSI offset: (% 10.4e, % 10.4e, % 10.4e)\n", 
	    det->stt_lsi_offset.x,
	    det->stt_lsi_offset.y,
	    det->stt_lsi_offset.z);
   fprintf (stdout, "STF-STT offset: (% 10.4e, % 10.4e, % 10.4e)\n", 
	    det->stf_stt_offset.x,
	    det->stf_stt_offset.y,
	    det->stf_stt_offset.z);

   if (sky == 0)
     dump_physical_geometry (det);
   else
     dump_sky_geometry (det, name);

   return 0;
}

static void compute_ra_dec (Marx_Chip_To_MNC_Type *chip_mnc,
			    unsigned int x, unsigned int y, 
			    double *ra, double *dec,
			    double *tx, double *ty)
{
   JDMVector_Type p;
   JDMVector_Type p0;

   (void) marx_chip_to_mnc (chip_mnc, x, y, &p);
   marx_mnc_to_ra_dec (&p, ra, dec);
   
   /* Convert to arc-sec from radians */
   *ra = *ra * (60.0 * 180.0/PI);
   *dec = *dec * (60.0 * 180.0/PI);

   p0.x = 1.0;
   p0.y = p0.z = 0;
   
   /* Flip sign of MNC */
   p.x = -p.x;
   p.y = -p.y;
   p.z = -p.z;
   
   (void) marx_vector_to_tan_plane (&p, &p0, tx, ty);
   *tx = *tx * (60.0 * 180.0/PI);
   *ty = *ty * (60.0 * 180.0/PI);
}


static double Focal_Length = 10065.5;
   
static void dump_sky_geometry (Marx_Detector_Type *det, char *name)
{
   int i;
   unsigned int min_x_pixel, max_x_pixel;
   unsigned int min_y_pixel, max_y_pixel;

   fprintf (stdout, "\nSky RA/Dec Geometry Follows (with respect to nominal)");

   min_x_pixel = 0;
   min_y_pixel = 0;
   max_x_pixel = det->num_x_pixels;
   max_y_pixel = det->num_y_pixels;

   for (i = det->first_chip_id; i <= det->last_chip_id;	i++)
     {
	Marx_Chip_To_MNC_Type *chip_mnc;
	double ra, dec;
	double tx, ty;
	
	fprintf (stdout, "\nChip %d:\n", i);
	fprintf (stdout, " Location of corners: (RA, Dec, RA--TAN, Dec--TAN (arc-min))\n");

	if (NULL == (chip_mnc = marx_allocate_chip_to_mnc (name, i)))
	  exit (1);

	if (-1 == marx_init_chip_to_mnc (chip_mnc, Focal_Length, 0, 0, 0, 0))
	  exit (1);

	compute_ra_dec (chip_mnc, min_x_pixel, min_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t% 10.4e\t(LL)\n", ra, dec, tx, ty);

	compute_ra_dec (chip_mnc, max_x_pixel, min_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t% 10.4e\t(LR)\n", ra, dec, tx, ty);

	compute_ra_dec (chip_mnc, max_x_pixel, max_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t% 10.4e\t(UR)\n", ra, dec, tx, ty);

	compute_ra_dec (chip_mnc, min_x_pixel, max_y_pixel, &ra, &dec, &tx, &ty);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t% 10.4e\t(UL)\n", ra, dec, tx, ty);
	
	marx_free_chip_to_mnc (chip_mnc);
     }
}


   
static void dump_physical_geometry (Marx_Detector_Type *det)
{
   unsigned int i;
   Marx_Detector_Geometry_Type *g;

   fprintf (stdout, "\nPhysical Geometry Follows (STF coords at nominal aimpoint, units in mm)\n");

   g = det->geom;

   for (i = 0; i < det->num_chips; i++)
     {
	fprintf (stdout, "\nChip %d:\n", g->id);
	fprintf (stdout, " X Length: % 10.4e (pixel-size: % 10.4e)\n", g->xlen, g->x_pixel_size);
	fprintf (stdout, " Y Length: % 10.4e (pixel-size: % 10.4e)\n", g->ylen, g->y_pixel_size);
	fprintf (stdout, " Location of corners:\n");
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t(LL)\n",
		 g->x_ll.x, g->x_ll.y, g->x_ll.z);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t(LR)\n",
		 g->x_lr.x, g->x_lr.y, g->x_lr.z);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t(UR)\n",
		 g->x_ur.x, g->x_ur.y, g->x_ur.z);
	fprintf (stdout, "\t% 10.4e\t% 10.4e\t% 10.4e\t(UL)\n",
		 g->x_ul.x, g->x_ul.y, g->x_ul.z);
	
	g++;
     }
}

   
     
