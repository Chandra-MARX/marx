#ifndef _MARX_CHANDRA_H_INCLUDED
#define _MARX_CHANDRA_H_INCLUDED

/* Chandra-specific structures */

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

#endif				       /* _MARX_CHANDRA_H_INCLUDED */
