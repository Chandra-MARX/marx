#ifndef _MARX_HRC_H_INCLUDED
#define _MARX_HRC_H_INCLUDED

/* stt_lsi_offset: origin of LSI in STT system. 
 *  This comes from table 18 of JMcD's coordinate memo.
 * 
 * stf_stt_offset: origin of STT in STF at nominal aimpoint for
 * the detector.  This value comes from table 19 of JMcD's coord
 * memo. 
 */
/* The following quantities have units of mm */
#define MARX_DETECTOR_TYPE_PRIVATE_DATA \
   JDMVector_Type stt_lsi_offset; \
   JDMVector_Type stf_stt_offset;

#define MARX_DET_FACET_PRIVATE_DATA \
   float tdet_xoff; \
   float tdet_yoff;

#endif				       /* _MARX_HRC_H_INCLUDED */
