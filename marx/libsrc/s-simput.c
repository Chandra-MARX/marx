/*
    This file is part of MARX

    Copyright (C) 2002-2015 Massachusetts Institute of Technology

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
#include <stdio.h>
#include "config.h"

#include "marx-feat.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>

#include <jdmath.h>
#include <pfile.h>

#include "marx.h"
#include "_marx.h"
#include "source.def"
// Those two files could be combined
#include "s-simput.h"
#include "simputcatalog.h"

#if MARX_HAS_DYNAMIC_LINKING
#include <dlfcn.h>

// TODO: deal with mjdref in some sensible way.

static char *Simput_handle;
static char *Simput_Source;
static  (*getARF)(int);



static double running_time = 0.;

double ra_nom, dec_nom, roll_nom;

// do something about mjdref

typedef struct
{
  double time;
  double energy;
  double dec;
  double ra;
  int lightcurve_status;
}
SIMPUT_generated_Photon_Type;

static SIMPUT_generated_Photon_Type* next_photons[];
static SourceCatalog *cat;


typedef void (*Fun_Ptr)(void);

//maybe do this later. For now, keep simple and link to SIMPUT at compile time.
static Fun_Ptr simput_dlsym (char *name, int is_required)
{
   Fun_Ptr f;

   if (Simput_handle == NULL)
     return NULL;

   f = (Fun_Ptr) dlsym (Simput_handle, name);
   if ((f == NULL) && is_required)
     {
	char *err = (char *)dlerror ();

	if (err == NULL) err = "Unknown";
	marx_error ("Unable to get symbol '%s' in SIMPUT library\nReason: %s\n",
		    name, err);
	return NULL;
     }
   return f;
}


int precompute_photon (SourceCatalog *cat, long sourcenumber, 
		       double mjdref, double prevtime)
{
  int* status;
  int photon_status;

  photon_status = getSimputPhoton(cat, cat->sources[sourcenumber], 
				  prevtime, mjdref, time, energy, ra, dec, status);
  if (*status!=0){
    marx_error ("Error generating photons in SIMPUT module.");
    return -1;
  }
  *next_photons[ii]->time = time;
  *next_photons[ii]->energy = energy;
  *next_photons[ii]->ra = ra;
  *next_photons[ii]->dec = dec;
  *next_photons[ii]->lightcurve_status = lightcurve_status;
  return 0;
}


/* would be really helpful if this was part of SIMPUT. If it was,
   I would not have to even define the ARF structure here, I could
   just pass araound a void pointer.
*/
struct ARF *get_constant_ARF(float low_energy, float high_energy, float effarea)
{
  // Make ARF for SIMPUT.
  // MARX will take care of photons lost along the path itself, 
  // so all we need here is a constant effective area equal to the opening
  // of the mirrors.
  ARF *arf;
  char *telescope = "Chandra";

  //if (NULL == (getARF = (struct ARF(*)(int*)) simput_dlsym ("getARF", 1)))
  //   return -1;

  *arf = getARF(int* const status)
  if (status != 0){
     return -1;
  };
  arf->NumberEnergyBins=1;
  if (NULL == (arf->LowEnergy = (double *)marx_malloc(arf->NumberEnergyBins*sizeof(double))))
    return -1;
  arf->LowEnergy[0] = low_energy; // keV
  if (NULL == (arf->HighEnergy = (double *)marx_malloc(arf->NumberEnergyBins*sizeof(double))))
    return -1;
  arf->HighEnergy[0] = high_energy; // keV
  if (NULL == (arf->EffArea = (double *)marx_malloc(arf->NumberEnergyBins*sizeof(double))))
    return -1;
  arf->EffArea[0] = effarea;

  if (-1 == strncopy(arf->Telescope, telescope, strlen(telescope)))
      return -1;
  return arf


// argv : first arg is SIMPUT filename
int simput_open_source (Marx_Source_Type *st)
{
  SimputSrc* const src;
  double mjdref = 52000.;
  double* const time;
  double* const energy;
  double* const ra;
  double* const dec;

  

  if (-1 == marx_get_nominal_pointing (&ra_nom, &dec_nom, &roll_nom))
    return -1;

  arf = get_constant_arf(0.3, 11., Marx_Mirror_Geometric_Area);
  cat=loadSourceCatalog(*Simput_Source, arf, status);

  if (*status!=0){
    marx_error ("Error interpreting SIMPUT catalog.");
  }

  if (NULL == (nextphotons = (SIMPUT_generated_Photon_Type *)marx_malloc((cat->nextsources+1)*sizeof(SIMPUT_generated_Photon_Type)))){
    return -1;
  }

  // Make one photon for each source and store it.
  long ii;
  int lightcurve_status
  for (ii=0; ii<cat->nextsources; ii++) {
    if (-1 == precompute_photon(cat, ii, mjdref, 0.)){
      return -1;
    }
  }
}

void simput_close_source (Marx_Source_Type *st)
{
  if (Simput_handle != NULL)
     dlclose (Simput_handle);

  Simput_handle = NULL;

  (void) free_ARF(*arf);
  (void) freeSourceCatalog(*cat, int* const status);
  marx_free (*next_photons);
  return *status;
}

static int simput_generate_ray (Marx_Photon_Attribute_Type *at)
{
  sp SIMPUT_generated_PhotonType;
  unsigned int i;
  JDMVector_Type normal, *p;
  double sigma_theta;

  double min_time = 1e40;
  long next_index = 0;
  long number_valid_photons = 0;
  JDMVector_Type p;
  double ra, dec;
  double Source_Azimuth, Source_Elevation;

  long ii;

  // From the list of pre-generated photons, find the next one
  for (ii=0; ii<cat->nextsources; ii++) {
    if ((*next_photons[ii]->time < min_time) && (*next_photons[ii]->lightcurve_status == 0)){
      min_time = *next_photons[ii]->time;
      next_index = ii;
      number_valid_photons++;
    }
  }
  if (number_valid_photons == 0){
    marx_error("None of the sources in the SIMPUT catalog has a lightcurve defined for this time.");
    return -1;
  }

  sp = *next_photons[next_index];
    
  if (-1 == precompute_photon(cat, ii, mjdref, sp->time)){
    return -1;
  }
  
  JDMVector_Type src, pnt, p;
  double ra_pnt, dec_pnt, roll_pnt;
  double az, el;
  
  src = JDMv_spherical_to_vector (1.0, PI/2.0 - sp->dec, sp->ra);

  // the next 2 lines can be moved to init
  if (-1 == marx_get_pointing (&ra_pnt, &dec_pnt, &roll_pnt))
    return -1;
  pnt = JDMv_spherical_to_vector (1.0, PI/2.0 - dec_pnt, ra_pnt);
 
  src = JDMv_rotate_unit_vector (src, pnt, -roll_pnt);
  JDMv_unit_vector_to_spherical (src, &el, &az);
  el = PI/2.0 - el;

  marx_compute_ra_dec_offsets (ra_pnt, dec_pnt, az, el, &az, &el);

  // line can be moved to init?
  p = JDMv_spherical_to_vector (1.0, 0.5*PI-zoff, yoff);

  /* Now add offsets via the proper rotations */
  p = JDMv_rotate_unit_vector (p, JDMv_vector (0, -1, 0), el);
  p = JDMv_rotate_unit_vector (p, JDMv_vector (0, 0, 1), az);

  /* Finally roll it so that this point will be invariant under roll.  That is,
   * the dither transformation will (on the average) undo this rotation.
   * See the apply_dither function.
   */
  *at->p = JDMv_rotate_unit_vector (p, JDMv_vector (1, 0, 0), roll_nom);
  *at->energy = sp->energy;
  *at->arrival_time = sp->time

  return 0;
}


int marx_select_simput_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   (void) source_id;
   st->open_source = simput_open_source;
   st->create_photons = simput_create_photons;
   st->close_source = simput_close_source;

   if (-1 == pf_get_double (p, "S-SIMPUT-Source", &Simput_Source))
     return -1;
   if (-1 == pf_get_double (p, "S-SIMPUT-Library", &Simput_handle))
     return -1;

   return 0;
}

/*}}}*/

static int create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
			   unsigned int num, unsigned int *num_created)
{
   unsigned int i;
   Marx_Photon_Attr_Type *at;
   int (*efun) (Marx_Spectrum_Type *, double *);
   double t, energy, last_time;

   at = pt->attributes;

   for (i = 0; i < num; i++)
     {
	if (-1 == simput_generate_ray (at))
	  break;

	at->flags = 0;
	at->arrival_time -=pt->start_time
	at->energy = energy;
	at++;
     }

   *num_created = i;

   // need to work through the times later...
   // not quite sure what I need to set here...
   t = 1.;
   if (t >= 0.0)
     {
	pt->history = (MARX_ENERGY_OK
		       | MARX_TIME_OK
		       | MARX_X_VECTOR_OK
		       | MARX_P_VECTOR_OK);
     }

   return 0;
}
/*}}}*/

#else
int marx_select_simput_source (Marx_Source_Type *st, Param_File_Type *p,
			     char *name, unsigned int source_id)
{
   (void) st;
   (void) p;
   (void) name;
   (void) source_id;

   marx_error ("This version of MARX does not support dynmamic linking");
   return -1;
}
#endif				       /* MARX_HAS_DYNAMIC_LINKING */
