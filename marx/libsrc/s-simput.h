/* // from cfitsio/fitsio.h */
/* #define FLEN_VALUE     71  /\* max length of a keyword value string *\/ */
/* #define FLEN_KEYWORD   72  /\* max length of a keyword (HIERARCH convention) *\/ */

/* // from heasp/heasp.h */
/* typedef struct  */
/* { */

/*   long NumberEnergyBins;                         /\* Number of response energies *\/ */

/*   float* LowEnergy; /\*NumberEnergyBins*\/         /\* Start energy of bin *\/ */
/*   float* HighEnergy; /\*NumberEnergyBins*\/        /\* End energy of bin *\/ */

/*   float* EffArea;    /\*NumberEnergyBins*\/        /\* Effective areas *\/ */

/*   char ARFVersion[FLEN_KEYWORD];                 /\* SPECRESP extension format version *\/ */
/*   char Telescope[FLEN_KEYWORD]; */
/*   char Instrument[FLEN_KEYWORD]; */
/*   char Detector[FLEN_KEYWORD]; */
/*   char Filter[FLEN_KEYWORD]; */
/*   char ARFExtensionName[FLEN_VALUE];             /\* Value of EXTNAME keyword in SPECRESP extension *\/ */

/* } */
/* ARF; */

typedef struct
{
  double time;
  double energy;
  double dec;
  double ra;
  int lightcurve_status;
}
SIMPUT_generated_Photon_Type;
