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

// todo: mjdref
// todo: source.c: 175 makes no sense for this source
// todo: dynamic linking?
// todo: make way to get zoff, yoff: Currently hardcoded to 0.
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


#include "s-simput.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include "simput.h"

#if MARX_HAS_DYNAMIC_LINKING
#include <dlfcn.h>

// TODO: deal with mjdref in some sensible way
static double mjdref=52000.;

static char Simput_handle[PF_MAX_LINE_LEN];
static char Simput_Source[PF_MAX_LINE_LEN];
static  (*getARF)(int);
static struct ARF * const arf;


static double running_time = 0.;

double ra_nom, dec_nom, roll_nom;

// do something about mjdref


static SIMPUT_generated_Photon_Type *next_photons;
static SimputCtlg * cat;

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


int precompute_photon (SimputCtlg *cat, long sourcenumber, 
		       double mjdref, double prevtime)
{
  int* status;
  int lightcurve_status;
  double* const time;
  float* const energy;
  double* const ra;
  double* const dec;
  SimputSrc* src;
  
  src = getSimputSrc(cat, sourcenumber, status);
  if (*status!=0){
    marx_error ("Could not retrieve source from SIMPUT catalog.");
    return -1;
  }

  lightcurve_status = getSimputPhoton(cat, src, 
				  prevtime, mjdref, time, energy, ra, dec, status);
  if (*status!=0){
    marx_error ("Error generating photons in SIMPUT module.");
    return -1;
  }
  next_photons[sourcenumber].time = *time;
  next_photons[sourcenumber].energy = *energy;
  next_photons[sourcenumber].ra = *ra;
  next_photons[sourcenumber].dec = *dec;
  next_photons[sourcenumber].lightcurve_status = lightcurve_status;
  return 0;
}


/* would be really helpful if this was part of SIMPUT. If it was,
   I would not have to even define the ARF structure here, I could
   just pass araound a void pointer.
*/
int get_constant_arf(float low_energy, float high_energy, float effarea, struct ARF * const arf)
{
  // Make ARF for SIMPUT.
  // MARX will take care of photons lost along the path itself, 
  // so all we need here is a constant effective area equal to the opening
  // of the mirrors.
  int* const status;
  //if (NULL == (getARF = (struct ARF(*)(int*)) simput_dlsym ("getARF", 1)))
  //   return -1;

  arf->NumberEnergyBins=1;
  if (NULL == (arf->LowEnergy = (float *)marx_malloc(arf->NumberEnergyBins*sizeof(float)))){
    marx_error ("out of memory");
    return -1;
  }
  arf->LowEnergy[0] = low_energy; // keV
  if (NULL == (arf->HighEnergy = (float *)marx_malloc(arf->NumberEnergyBins*sizeof(float)))){
    marx_error ("out of memory");
    return -1;
  }
  arf->HighEnergy[0] = high_energy; // keV
  if (NULL == (arf->EffArea = (float *)marx_malloc(arf->NumberEnergyBins*sizeof(float)))){
    marx_error ("out of memory");
    return -1;
  }
  arf->EffArea[0] = effarea;

  strcpy(arf->ARFVersion,"1.0");          /* SPECRESP extension format version */
  strcpy(arf->Telescope, "Chandra");
  strcpy(arf->Instrument, "any");
  strcpy(arf->Detector, "any");
  strcpy(arf->Filter, "none");
  strcpy(arf->ARFExtensionName, "none");
  return 0;
}


// argv : first arg is SIMPUT filename
int simput_open_source (Marx_Source_Type *st)
{
  SimputSrc* const src;
  double mjdref = 52000.;
  double* const time;
  double* const energy;
  double* const ra;
  double* const dec;
  int* const status;
  long n_sources;

  cat = openSimputCtlg(Simput_Source, READONLY, 0, 0, 0, 0, status);
  if (*status!=0){
    marx_error ("Error interpreting SIMPUT catalog.");
    return -1;
  }
  n_sources = getSimputCtlgNSources(cat);
  if (n_sources==0){
    marx_error("No Sources found in SIMPUT catalog");
    return -1;
  }

  if (-1 == marx_get_nominal_pointing (&ra_nom, &dec_nom, &roll_nom))
    return -1;

  if (-1 == get_constant_arf(0.3, 11., Marx_Mirror_Geometric_Area, arf))
    return -1;
  setSimputARF(cat, arf);


  if (NULL == (next_photons = marx_malloc(n_sources * sizeof(*next_photons)))){
    return -1;
  }

  // Make one photon for each source and store it.
  long ii;
  int lightcurve_status;
    for (ii=1; ii<(n_sources+1); ii++) {
    if (-1 == precompute_photon(cat, ii, mjdref, 0.)){
      return -1;
    }
  }
}

static int simput_close_source (Marx_Source_Type *st)
{
  int * const status;
  /* if (Simput_handle != NULL) */
  /*    dlclose (Simput_handle); */

  /* Simput_handle = NULL; */

  //(void) free_ARF(*arf);
  (void) freeSimputCtlg(&cat, status);
  marx_free (next_photons);
  (void) st;
  return *status;
}

static int simput_generate_ray (Marx_Photon_Attr_Type *at)
{
  SIMPUT_generated_Photon_Type sp;
  unsigned int i;

  double min_time = 1e40;
  long next_index = 0;
  long number_valid_photons = 0;
  double ra, dec;
  double Source_Azimuth, Source_Elevation;
  JDMVector_Type src, pnt, p;
  double ra_pnt, dec_pnt, roll_pnt;
  double az, el;

  long ii, n_sources;

  // From the list of pre-generated photons, find the next one
  n_sources = getSimputCtlgNSources(cat);

  for (ii=1; ii < (n_sources+1); ii++) {
    if ((next_photons[ii].time < min_time) && (next_photons[ii].lightcurve_status == 0)){
      min_time = next_photons[ii].time;
      next_index = ii;
      number_valid_photons++;
    }
  }
  if (number_valid_photons == 0){
    marx_error("None of the sources in the SIMPUT catalog has a lightcurve defined for this time.");
    return -1;
  }

  sp = next_photons[next_index];
    
  // Make a new photons
  if (-1 == precompute_photon(cat, ii, mjdref, sp.time)){
    return -1;
  }
    
  // Bring the selected photon from RA/DEC to MARX system
  /* This might not be the most efficient way to do it.
   * All these spherical operations are quite expensive,
   * but for now I'll jus make it work - no pre-mature optimization.
   */
  src = JDMv_spherical_to_vector (1.0, PI/2.0 - sp.dec, sp.ra);

  // the next 2 lines can be moved to init
  if (-1 == marx_get_pointing (&ra_pnt, &dec_pnt, &roll_pnt))
    return -1;
  pnt = JDMv_spherical_to_vector (1.0, PI/2.0 - dec_pnt, ra_pnt);
 
  src = JDMv_rotate_unit_vector (src, pnt, -roll_pnt);
  JDMv_unit_vector_to_spherical (src, &el, &az);
  el = PI/2.0 - el;

  marx_compute_ra_dec_offsets (ra_pnt, dec_pnt, az, el, &az, &el);

  // line can be moved to init?
  double zoff=0;
  double yoff=0;
  p = JDMv_spherical_to_vector (1.0, 0.5*PI-zoff, yoff);

  /* Now add offsets via the proper rotations */
  p = JDMv_rotate_unit_vector (p, JDMv_vector (0, -1, 0), el);
  p = JDMv_rotate_unit_vector (p, JDMv_vector (0, 0, 1), az);

  /* Finally roll it so that this point will be invariant under roll.  That is,
   * the dither transformation will (on the average) undo this rotation.
   * See the apply_dither function.
   */
  at->p = JDMv_rotate_unit_vector (p, JDMv_vector (1, 0, 0), roll_nom);
  at->energy = sp.energy;
  at->arrival_time = sp.time;

  return 0;
}


static int simput_create_photons (Marx_Source_Type *st, Marx_Photon_Type *pt, /*{{{*/
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
	at->arrival_time -=pt->start_time;
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

int marx_select_simput_source (Marx_Source_Type *st, Param_File_Type *p, /*{{{*/
			      char *name, unsigned int source_id)
{
   (void) source_id;
   st->open_source = simput_open_source;
   st->create_photons = simput_create_photons;
   st->close_source = simput_close_source;

   if (-1 == pf_get_file (p, "S-SIMPUT-Source", Simput_Source, sizeof(Simput_Source)))
     return -1;
   if (-1 == pf_get_file (p, "S-SIMPUT-Library", Simput_handle, sizeof(Simput_handle)))
     return -1;

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
