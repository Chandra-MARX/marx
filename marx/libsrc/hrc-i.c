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
#include <stdlib.h>
#endif
#include <string.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"


#define NUM_FILTER_REGIONS	1

static _Marx_HRC_QE_Type MCP_QEs [_MARX_NUM_HRC_I_CHIPS];
static _Marx_HRC_QE_Type Filter_QEs [NUM_FILTER_REGIONS];

static double HRC_I_Blur;
static double X_Pixel_Size;
static double Y_Pixel_Size;
static Marx_Detector_Geometry_Type *HRC_I_Geom;

static Param_Table_Type Parm_Table [] =
{
   {"HRC-I-BlurSigma",	PF_REAL_TYPE,		&HRC_I_Blur},
   {"HRC-I-QEFile",	PF_FILE_TYPE,		&MCP_QEs[0].file},
   {"HRC-I-UVISFile",	PF_FILE_TYPE,		&Filter_QEs[0].file},
   {NULL, 0, NULL}
};


void _marx_hrc_blur_position (double *dx, double *dy, double blur)
{
   double r, theta;
	     
   if (_Marx_Det_Ideal_Flag)
     return;

   theta = (2.0 * PI) * JDMrandom ();
   r = blur * JDMgaussian_random ();
	     
   *dx += r * cos(theta);
   *dy += r * sin (theta);
   if (*dx < 0.0) *dx = 0.0;
   if (*dy < 0.0) *dy = 0.0;
}



short _marx_hrc_compute_pha (double energy)
{
   double factor_1 = 141.582;
   double factor_2 = 107.299;
   double factor_3 = 115.0;
   double width_factor = 0.424661;     /* 1/sqrt(2log(2)) */
   
   if (energy <= 0.5)
     energy = factor_1 * sqrt(energy);
   else if (energy < 2.0)
     energy = factor_2 * pow (energy, 0.1);
   else
     energy = factor_3;
   
   /* We know that delta_E/E = 1 for this detector, where delta_E is FWHM */
   energy = energy * (1.0 + width_factor * JDMgaussian_random ());
   if (energy < 0.0) energy = 0.0;

   return (short) energy;
}

static int 
apply_hrc_qe (Marx_Photon_Attr_Type *at, unsigned int mcp_id)
{
   _Marx_HRC_QE_Type *mcp;
   double qe;
   double energy;

   energy = at->energy;
   mcp = MCP_QEs + mcp_id;
   if (mcp->num_energies != 0)
     {
	qe = JDMinterpolate_f ((float) energy, 
			       mcp->energies, mcp->eff, mcp->num_energies);
	if (JDMrandom () >= qe)
	  return -1;
     }

   mcp = Filter_QEs + mcp_id;
   if (mcp->num_energies != 0)
     {
	qe = JDMinterpolate_f ((float) energy, 
			       mcp->energies, mcp->eff, mcp->num_energies);
	if (JDMrandom () >= qe)
	  return -1;
     }
   
   at->pulse_height = _marx_hrc_compute_pha (energy);

   return 0;
}

int 
_marx_hrc_i_detect (Marx_Photon_Type *pt) /*{{{*/
{
   Marx_Photon_Attr_Type *at, *attrs;
   unsigned int n_photons, i;
   unsigned int *sorted_index;

   if (pt->history & MARX_CCD_NUM_OK)
     return 0;
   
   pt->history |= (MARX_HRC_REGION_OK | MARX_PULSEHEIGHT_OK 
		   | MARX_Y_PIXEL_OK | MARX_Z_PIXEL_OK | MARX_CCD_NUM_OK);

   marx_prune_photons (pt);

   attrs = pt->attributes;
   n_photons = pt->num_sorted;
   sorted_index = pt->sorted_index;
   
   for (i = 0; i < n_photons; i++)
     {
	Marx_Detector_Geometry_Type *d;
	double dx, dy;
	
	at = attrs + sorted_index[i];
	
#if MARX_HAS_DITHER
	_marx_dither_detector (&at->dither_state);
#endif
	/* Transform ray into local system */
	_marx_transform_ray (&at->x, &at->p,
			    &_Marx_Det_XForm_Matrix);

	d = _marx_intersect_with_detector (at->x, at->p,
					   HRC_I_Geom, _MARX_NUM_HRC_I_CHIPS,
					   &at->x, &dx, &dy);
	if (d == NULL)
	  {
	     at->flags |= PHOTON_MISSED_DETECTOR;
	     at->ccd_num = -1;
	  }
	else if (-1 == apply_hrc_qe (at, d->id))
	  {
	     at->flags |= PHOTON_UNDETECTED;
	     at->ccd_num = -1;
	  }
	else
	  {
	     _marx_hrc_blur_position (&dx, &dy, HRC_I_Blur);

	     at->ccd_num = d->id;

	     (void) _marx_hrc_i_compute_pixel (dx, dy, &dx, &dy);
	     at->y_pixel = dx;
	     at->z_pixel = dy;
#if 0	     
	     at->y_pixel = dx / X_Pixel_Size;
	     at->z_pixel = dy / Y_Pixel_Size;
	     
	     /* Hack to make aimpoint agree with SMs email about 
	      * latest best guess geometry of HRC-I.
	      * Thu Nov  5 16:56:49 1998
	      * 
	      * Apparantly, the "true" value will not be known until it is
	      * in flight.
	      */
	     at->y_pixel += (7669.5 - 7201.74);
	     at->z_pixel += (7723.5 - 7201.74);
#endif
	  }
	
	/* Transform detected photon back to original system */
	_marx_transform_ray_reverse (&at->x, &at->p,
				     &_Marx_Det_XForm_Matrix);
#if MARX_HAS_DITHER
	_marx_undither_detector (&at->dither_state);
#endif
     }
   
   return 0;
}

/*}}}*/


int 
_marx_hrc_i_init (Param_File_Type *p) /*{{{*/
{
   unsigned int i;
   Marx_Detector_Type *hrc;
   
   if (-1 == pf_get_parameters (p, Parm_Table))
     return -1;
   
   if (-1 == _marx_hrc_i_geom_init (p))
     return -1;
   
   if (NULL == (hrc = marx_get_detector_info ("HRC-I")))
     return -1;

   HRC_I_Geom = hrc->geom;
   X_Pixel_Size = hrc->x_pixel_size;
   Y_Pixel_Size = hrc->y_pixel_size;

   marx_message ("Reading binary HRC-I QE/UVIS data files:\n");

   for (i = 0; i < _MARX_NUM_HRC_I_CHIPS; i++)
     {
	if (-1 == _marx_hrc_read_efficiencies (MCP_QEs + i))
	  return -1;
     }
   
   for (i = 0; i < NUM_FILTER_REGIONS; i++)
     {
	if (-1 == _marx_hrc_read_efficiencies (Filter_QEs + i))
	  return -1;
     }

   return 0;
}

