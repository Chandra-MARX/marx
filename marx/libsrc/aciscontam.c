/*
    This file is part of MARX

    Copyright (C) 2002-2009 Massachusetts Institute of Technology

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

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"

/* This model is derived from AV's May 9 2004 memo "Spatial structure in the 
 * ACIS OBF contamination"
 */
static double T_Launch = 1999.56;

#define TAU_XXX(func,a,b) \
   static double func (double t) { return ((a)*(1-exp(-(t)/(b)))); }

TAU_XXX(tau_CS, 0.63997, 2.0541)
TAU_XXX(tau_BS, 0.87532, 1.4093)
TAU_XXX(tau_TS, 0.81544, 1.5362)
TAU_XXX(tau_CI, 0.67499, 1.4042)
TAU_XXX(tau_FI, 0.96650, 1.2951)

/* These will be precomputed for the observation time */
static double Tau_CS, Tau_BS, Tau_TS, Tau_CI, Tau_FI;

static float *Energies, *Mus;
static unsigned int Num_Energies;

/* ACIS-I Functions  */

static double ArcMin_Per_Pixel = 8.0/1024.0;

/* Eq 19 */
static double contam_acis_i (double en, double x, double y, double x_0, double y_0)
{
   double gamma = 2.0;
   double dx = x - x_0;
   double dy = y - y_0;
   double r = ArcMin_Per_Pixel * sqrt(dx*dx + dy*dy);
   double a = pow (8.07, gamma);
   double b = pow (3.24, gamma);
   double z = Tau_CI + (Tau_FI - Tau_CI) * ((pow (r, gamma) - b)/(a-b));
   double mu = JDMinterpolate_f (en, Energies, Mus, Num_Energies);
   return exp (-mu*z);
}

static double contam_acis_s (double en, double x, double y)
{
   /* These come from eq 8 */
   double y_0 = 512.0;
   double alpha_1 = 5.5;
   double alpha_2 = 4.5;
   double z, mu;
   
   (void) x;

   if (y <= 512)
     z = Tau_CS + (Tau_BS - Tau_CS)*pow(fabs((y-y_0)/(64-y_0)), alpha_1);
   else
     z = Tau_CS + (Tau_TS - Tau_CS)*pow(fabs((y-y_0)/(964-y_0)), alpha_2);
   
   mu = JDMinterpolate_f (en, Energies, Mus, Num_Energies);
   return exp (-mu*z);
}

static double contam_i0 (_Marx_Acis_Chip_Type *ccd, double en, double cx, double cy)
{
   (void) ccd;
   return contam_acis_i (en, cx, cy, 1024.5, 1024.5);
}
static double contam_i1 (_Marx_Acis_Chip_Type *ccd, double en, double cx, double cy)
{
   (void) ccd;
   return contam_acis_i (en, cx, cy, 0.5, 1024.5);
}
static double contam_i2 (_Marx_Acis_Chip_Type *ccd, double en, double cx, double cy)
{
   (void) ccd;
   return contam_acis_i (en, cx, cy, 0.5, 1024.5);
}
static double contam_i3 (_Marx_Acis_Chip_Type *ccd, double en, double cx, double cy)
{
   (void) ccd;
   return contam_acis_i (en, cx, cy, 1024.5, 1024.5);
}

static double contam_s012345 (_Marx_Acis_Chip_Type *ccd, double en, double cx, double cy)
{
   (void) ccd;
   return contam_acis_s (en, cx, cy);
}

static int pre_init (Param_File_Type *p)
{
   char *file;
   double years_since_launch;

   (void) p;

   if (Mus != NULL)
     return 0;			       /* already initialized */
   
   if (NULL == (file = marx_make_data_file_name ("acis/aciscontam.dat")))
     return -1;

   marx_message ("\t%s\n", file);

   if (-1 == marx_f_read_bdat (file, &Num_Energies, 2, &Energies, &Mus))
     {
	marx_free (file);
	return -1;
     }
   marx_free (file);
   
   years_since_launch = _Marx_TStart_Yrs - T_Launch;
   
   Tau_CS = tau_CS (years_since_launch);
   Tau_BS = tau_BS (years_since_launch);
   Tau_TS = tau_TS (years_since_launch);
   Tau_CI = tau_CI (years_since_launch);
   Tau_FI = tau_FI (years_since_launch);
   
   return 0;
}

int _marx_acis_contam_init (Param_File_Type *p, _Marx_Acis_Chip_Type *ccd)
{
   if (-1 == pre_init (p))
     return -1;
   
   switch (ccd->ccd_id)
     {
      case 0:
	ccd->contam_fun = contam_i0;
	break;
      case 1:
	ccd->contam_fun = contam_i1;
	break;
      case 2:
	ccd->contam_fun = contam_i2;
	break;
      case 3:
	ccd->contam_fun = contam_i3;
	break;
      case 4:
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
	ccd->contam_fun = contam_s012345;
	break;
	
      default:
	marx_error ("_marx_acis_contam_init: unsupported ccdid: %d", ccd->ccd_id);
	return -1;
     }
   return 0;
}

