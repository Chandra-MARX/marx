#ifndef _MARX_IXOCCD_H_INCLUDED
#define _MARX_IXOCCD_H_INCLUDED

typedef struct _IXO_CCD_QE_Type IXO_CCD_QE_Type;

#define MARX_DET_FACET_PRIVATE_DATA \
   IXO_CCD_QE_Type *qeinfo; \
   double read_noise; \
   double energy_gain; \
   double fano_factor;

#endif				       /* _MARX_IXOCCD_H_INCLUDED */
