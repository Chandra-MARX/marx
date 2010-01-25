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
#include <math.h>
#include <jdmath.h>

#include "marx.h"
#include "_marx.h"

#if 0
static char *Geom_File = "pixlib/pix_corner_lsi.par";

static int patch_acis_geom_internal (Marx_Detector_Geometry_Type *g,
				     char i_or_s, 
				     unsigned int first, unsigned last,
				     Param_File_Type *pf)
{
   char parm[32];
   unsigned int i;

   for (i = first; i <= last; i++)
     {
	sprintf (parm, "ACIS-%c%u-LL", i_or_s, i);
	if (-1 == _marx_get_vector_parm (pf, parm, &g->x_ll))
	  return -1;

	sprintf (parm, "ACIS-%c%u-LR", i_or_s, i);
	if (-1 == _marx_get_vector_parm (pf, parm, &g->x_lr))
	  return -1;

	sprintf (parm, "ACIS-%c%u-UL", i_or_s, i);
	if (-1 == _marx_get_vector_parm (pf, parm, &g->x_ul))
	  return -1;

	sprintf (parm, "ACIS-%c%u-UR", i_or_s, i);
	if (-1 == _marx_get_vector_parm (pf, parm, &g->x_ur))
	  return -1;

	g++;
     }
   return 0;
}


static int patch_acis_geom (Marx_Detector_Type *d, char i_or_s,
			    unsigned int first, unsigned int last)
{
   char *file;
   Param_File_Type *pf;
   int status;

   
   if (NULL == (file = marx_make_data_file_name (Geom_File)))
     return -1;

   /* marx_message ("\t%s\n", file); */
   pf = pf_open_parameter_file (file, "rQ");
   if (pf == NULL)
     {
	marx_error ("Unable to open %s", file);
	marx_free (file);
	return -1;
     }
   
   status = patch_acis_geom_internal (d->geom, i_or_s, first, last, pf);
   (void) pf_close_parameter_file (pf);
   marx_free (file);
   return status;
}

int _marx_patch_acis_i_geom (Marx_Detector_Type *d)
{
   return patch_acis_geom (d, 'I', 0, 3);
}

int _marx_patch_acis_s_geom (Marx_Detector_Type *d)
{
   return patch_acis_geom (d, 'S', 0, 5);
}

#else
int _marx_patch_acis_i_geom (Marx_Detector_Type *d)
{
   return _marx_caldb_patch_acis_geom (d);
}

int _marx_patch_acis_s_geom (Marx_Detector_Type *d)
{
   return _marx_caldb_patch_acis_geom (d);
}
#endif
