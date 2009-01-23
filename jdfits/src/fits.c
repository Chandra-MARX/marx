/*
    Copyright (C) 2002 MIT Center For Space Research

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


#include <memory.h>
#include <ctype.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "jdfits.h"
#include "_jdfits.h"

#ifndef SLMALLOC
# define SLMALLOC malloc
# define SLCALLOC calloc
# define SLREALLOC realloc
# define SLFREE free
#endif

   
static void handle_comment (JDFits_Keyword_Type *kwt, char *p, char *pmax)
{
   while ((p < pmax) && (*p == ' ')) p++;
   if (p == pmax)
     {
	kwt->comment = NULL;
     }
   else 
     {
	kwt->comment = p;
	if (*p != '/')
	  jdfits_warning ("Keyword %s has a comment that does not start with /.",
			kwt->name);
     }
   kwt->comment_len = (int) (pmax - p);
}

   
static int parse_string (JDFits_Keyword_Type *kwt)
{
   char *pend, *pmax, *pbeg, *p1, *p2;

   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   
   pmax = kwt->comment;		       /* for unparsed cards, this
					* points to the end of the 
					* card
					*/

     
   pbeg = kwt->v.sval;	       
   if (*pbeg != '\'')
     {
	/* This is an ugly hack to handle cards without values!
	 * FIXME!!!!
	 */
	kwt->v.sval = pbeg = pmax - 1;
	*pbeg = 0;
	kwt->comment = NULL;
	return 0;
     }
   
   pbeg++;			       /* skip by opening ' */
   
   p2 = p1 = pbeg;
   while (p1 < pmax)
     {
	unsigned char ch = *p1;
	
	if (ch == '\'')
	  {
	     p1++;
	     if ((p1 == pmax) || (*p1 != '\''))
	       {
		  p1--;
		  break;
	       }
	  }
	*p2++ = ch;
	p1++;
     }
   
   if (p1 == pmax)
     {
	jdfits_warning ("Card %s lacks closing ' character.", kwt->name);
	p1--;
	if (p2 == pmax) p2--;
     }
   p2--;
   
   /* Trailing blanks not significant */
   while (*p2 == ' ') p2--;
   *(p2 + 1) = 0;
   
   kwt->v.sval = pbeg;
   kwt->type = JDFITS_STRING_TYPE;
   
   pend = p1 + 1;
   handle_comment (kwt, pend, pmax);
   return 0;
}


static int parse_float32 (JDFits_Keyword_Type *kwt)
{
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   jdfits_error ("parse_float32: Not Implemented.");
   return -1;
}

static int parse_float64 (JDFits_Keyword_Type *kwt)
{
   char *p, *pbeg, *pmax;
   float64 x;
   
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   
   pbeg = kwt->v.sval;
   
   pmax = kwt->comment;		       /* for unparsed cards, this
					* points to the end of the 
					* card
					*/
   
   /* Change [+-]ddd.dddD to [+-]ddd.dddE */
   p = pbeg;
   if ((p < pmax)
       && ((*p == '+') || (*p == '-')))
     p++;

   while ((p < pmax) && isdigit(*p))
     p++;

   if ((p < pmax) && (*p == '.'))
     p++;
   
   while ((p < pmax) && isdigit(*p))
     p++;
   
   if (*p == 'D')
     *p = 'E';

   p = pbeg;   
   if (1 != sscanf (p, "%lf", &x))
     {
	jdfits_error ("Error parsing floating value for %s.", kwt->name);
	return -1;
     }
   
   while (p < pmax)
     {
	if (*p == '/') break;
	p++;
     }
   
   kwt->type = JDFITS_FLOAT64_TYPE;
   kwt->v.dval = x;
   
   handle_comment (kwt, p, pmax);
   return 0;
}

static int parse_int16 (JDFits_Keyword_Type *kwt)
{
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   jdfits_error ("parse_int16: Not Implemented.");
   return -1;
}

static int parse_int32 (JDFits_Keyword_Type *kwt)
{
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   jdfits_error ("parse_int32: Not Implemented.");
   return -1;
}

static int parse_comment (JDFits_Keyword_Type *kwt)
{
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   jdfits_error ("parse_comment: Not Implemented.");
   return -1;
}


static int parse_int_as_hex (JDFits_Keyword_Type *kwt)
{
   char *p, *pbeg, *pmax;
   int32 i;
   int sign = 1;
   
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   
   pbeg = kwt->v.sval;
   
   pmax = kwt->comment;		       /* for unparsed cards, this
					* points to the end of the 
					* card
					*/
   
   jdfits_warning ("\
Keyword %s appeard to have an integer value expressed in HEX form.\n\
  This file does NOT conform to the fits standard-- it is NOT a fits file.",
		 kwt->name);

   p = pbeg;
   if (p < pmax)
     {
	if (*p == '+') 
	  {
	     sign = 1;
	     p++;
	  }
	else if (*p == '-')
	  {
	     sign = -1;
	     p++;
	  }
     }
   else sign = 1;
	    
   i = 0;
   while (p < pmax)
     {
	long ich = *p;
	
	if ((ich <= '9') && (ich >= '0'))
	  ich = ich - '0';
	else if (((ich = ich | 0x20) >= 'a') 
		 && (ich <= 'f'))
	  ich = 10 + (ich - 'a');
	else break;
	
	i = 16 * i + ich;
	p++;
     }
   i = i * sign;
   
   while (p < pmax)
     {
	if (*p != ' ')
	  {
	     if (*p == '/') break;
	     jdfits_warning ("Excess junk at end of integer field for %s.", kwt->name);
	     break;
	  }
	p++;
     }
   
   kwt->type = JDFITS_INT_TYPE;
   kwt->v.ival = i;
   
   handle_comment (kwt, p, pmax);
   return 0;
}

static int parse_int (JDFits_Keyword_Type *kwt)
{
   char *p, *pbeg, *pmax;
   int32 i;
   long ilong;

   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   
   pbeg = kwt->v.sval;
   
   pmax = kwt->comment;		       /* for unparsed cards, this
					* points to the end of the 
					* card
					*/
   
   /* Check for a hex value--- it is illegal but it happens :( */
   p = pbeg;
   while (p < pmax)
     {
	if ((*p == '/') || (*p == ' ')) 
	  break;

	if ((unsigned int) (*p | 0x20) >= (unsigned int) 'a') 
	  return parse_int_as_hex (kwt);
	p++;
     }
   
   p = pbeg;
   
   if (*p == '+') p++;
   if (1 != sscanf ((char *) p, "%ld", &ilong))
     {
	jdfits_error ("Error parsing integer value for %s.", kwt->name);
	return -1;
     }
   i = (int32) ilong;

   while (p < pmax)
     {
	if (!isdigit (*p) && (*p != ' ') && (*p != '-'))
	  {
	     if (*p == '/') break;
	     jdfits_warning ("Excess junk at end of integer field for %s.", kwt->name);
	     break;
	  }
	p++;
     }
   
   kwt->type = JDFITS_INT_TYPE;
   kwt->v.ival = i;
   
   handle_comment (kwt, p, pmax);
   return 0;
}

static int parse_bool (JDFits_Keyword_Type *kwt)
{
   char *p, *pmax;
   if ((kwt->type & JDFITS_UNPARSED_TYPE) == 0) return 0;
   
   p = kwt->v.sval;
   pmax = kwt->comment;
   kwt->v.ival = (*p == 'T');
   kwt->type = JDFITS_BOOL_TYPE;
   handle_comment (kwt, p + 1, pmax);
   return 0;
}

static int is_commentary_keyword (char *name)
{
   return ((0 == jdfits_strcasecmp (name, "HISTORY"))
	   || (0 == jdfits_strcasecmp (name, "COMMENT")));
}


static int parse_headers (JDFits_Header_Type *ht)
{
   unsigned char *p, ch;
   unsigned int i;
   unsigned int n = ht->num_keywords;
   int naxis;
   int bitpix;
   JDFits_Keyword_Type *kwt;
   off_t naxis_product;
   
   kwt = (JDFits_Keyword_Type *) SLCALLOC (n, sizeof (JDFits_Keyword_Type));
   
   if (kwt == NULL)
     {
	jdfits_error ("parse_headers: Memory allocation error.");
	return -1;
     }
   
   ht->keys = kwt;
   
   p = ht->header_data_buf;
   
   for (i = 0; i < n; i++)
     {
	unsigned char *p1 = p + 8, *nextp = p + JDFITS_CARD_SIZE;

	kwt->name = (char *) p;
	ch = *p1;
	
	*p1-- = 0;
	while ((p1 >= p) && (*p1 == ' '))
	  *p1-- = 0;
	
	if ((ch == '=')
	    && (0 == is_commentary_keyword (kwt->name)))
	  {
	     /* value */
	     
	     if (p1 + 1 == p)
	       {
		  jdfits_warning ("Card %d appears corrupt.", i);
	       }
	     
	     p1 = p + 9;
	     while ((p1 < nextp) && (*p1 == ' ')) p1++;
	     
	     if ((p1 == nextp) || (*p1 == '/'))
	       {
		  jdfits_warning ("Card %d (%s) has no value.", i, kwt->name);
		  kwt->type = JDFITS_UNPARSED_TYPE | JDFITS_STRING_TYPE;
	       }
	     else switch (*p1)
	       {
		case '\'':
		  kwt->type = JDFITS_STRING_TYPE | JDFITS_UNPARSED_TYPE;
		  break;
		  
		case 'F':
		case 'T': 
		  kwt->type = JDFITS_BOOL_TYPE | JDFITS_UNPARSED_TYPE;
		  if (p1 == p + 29) break;
		  
		  if (*p1 != 'F')
		    {
		       jdfits_warning ("Boolean card %d (%s) has value in wrong column.", i, kwt->name);
		       break;
		    }
		  /* Drop */
		  
		case 'A': case 'a':
		case 'B': case 'b':
		case 'C': case 'c':
		case 'D': case 'd':
		case 'E': case 'e':
		case 'f':
		  /* we can have integers, floats or decimals.  
		   * Anything else? 
		   */
		case '-': case '+': case '.':
		case '0': case '1': case '2': case '3': case '4': case '5':
		case '6': case '7': case '8': case '9':
		  kwt->type = JDFITS_NUMBER_MASK | JDFITS_UNPARSED_TYPE;
		  break;
		  
		default:
		  jdfits_warning ("Card %d (%s) appears corrupt. Assuming a string.",
				i, kwt->name);
		  kwt->type = JDFITS_STRING_TYPE | JDFITS_UNPARSED_TYPE;
	       }
	     
	     kwt->v.sval = (char *)p1;
	     kwt->comment = (char *) nextp;	       /* for unparsed cards, this
						* field points to the end
						*/
	  }
	/* END of VALUE CARDS */
	else
	  {
	     /* According to the 1990 FITS User's Guide, any thing else is a
	      * comment card
	      */
	     kwt->type = JDFITS_COMMENT_TYPE;
	     p[8] = ch;
	     kwt->v.sval = kwt->comment = (char *)p + 8;
	     kwt->comment_len = (int) ((char *)nextp - kwt->comment);
	  }
	
	p = nextp + 1;
	kwt++;

     } /* End of for loop */
   
   /* Do a minimial parse so that the data size and name of the header or 
    * extension can be determined.  According to the 1993 NOST 100-1.0 
    * definition of the FITS format, the header must begin with the keywords:
    * @ SIMPLE
    * @ BITPIX
    * @ NAXIS
    * @ NAXIS[n], n = 1..NAXIS
    * @ EXTEND (presence not required)
    * The size of the primary header is:
    * @ abs(BITPIX) * NAXIS[1] * ... * NAXIS[n]
    * If the header represents an extension, it looks like:
    * @ XTENSION
    * @ BITPIX
    * @ NAXIS
    * @ NAXIS[n], n = 1..NAXIS
    * @ (zero or more other keywords including:)
    * @ PCOUNT
    * @ GCOUNT
    * The size of the extension is: 
    * @ abs(BITPIX) * GCOUNT * (PCOUNT + NAXIS[1] * ... * NAXIS[n]) 
    */
   
   kwt = ht->keys;
   ht->name = NULL;
   
   if ((kwt->type & JDFITS_STRING_TYPE) && !strcmp (kwt->name, "XTENSION"))
     {
	ht->type = JDFITS_EXTENSION_HEADER;
	if (0 == parse_string (kwt)) 
	  {
	     ht->name = kwt->v.sval;
	     if (!strcmp ((char *) ht->name, "BINTABLE")) ht->type = JDFITS_BINTABLE_HEADER;
	  }
     }
   else if ((kwt->type & JDFITS_BOOL_TYPE) && !strcmp(kwt->name, "SIMPLE"))
     {
	ht->type = JDFITS_SIMPLE_HEADER;
	(void) parse_bool (kwt);
     }
   else
     {
	jdfits_error ("Header is not a fits header. Expected to see SIMPLE or XTENSION.");
	goto error_return;
     }
   
   if ((0 == (kwt[1].type & JDFITS_NUMBER_MASK))
       || strcmp (kwt[1].name, "BITPIX")
       || (-1 == parse_int (kwt + 1)))
     {
	jdfits_error ("BITPIX keyword missing or has a bad value.");
	goto error_return;
     }
   ht->bitpix = bitpix = kwt[1].v.ival;
   
   if ((bitpix != 8) && (bitpix != 16) && (bitpix != 32)
       && (bitpix != -32) && (bitpix != -64))
     jdfits_warning ("The BITPIX keyword has a non-standard value: %d", bitpix);
   
   if ((0 == (kwt[2].type & JDFITS_NUMBER_MASK))
       || strcmp (kwt[2].name, "NAXIS")
       || (-1 == parse_int (kwt + 2))
       || ((naxis = kwt[2].v.ival) < 0))
     {
	jdfits_error ("NAXIS keyword missing or has a bad value.");
	goto error_return;
     }
   ht->naxis = naxis;
   
   if (naxis) naxis_product = 1; else naxis_product = 0;
   
   kwt += 3;
   for (i = 1; i <= (unsigned int) naxis; i++)
     {
	char naxis_str[9];
	
	sprintf (naxis_str, "NAXIS%d", i);
	
	if ((0 == (kwt->type & JDFITS_NUMBER_MASK))
	    || strcmp (kwt->name, naxis_str)
	    || (-1 == parse_int (kwt))
	    || ((naxis_product = naxis_product * kwt->v.ival) < 0))
	  {
	     jdfits_error ("NAXIS keyword missing or has a bad value.");
	     goto error_return;
	  }
	kwt++;
     }
   
   if (naxis)
     {
	ht->kw_naxis1 = kwt - naxis;
     }
   else ht->kw_naxis1 = NULL;
   
   /* Now if this is an extension header, look for the PCOUNT and GCOUNT
    * fields.  The reference makes it seem that PCOUNT comes first.  However,
    * other documentation is not as explicit.
    */
   
   if (ht->type != JDFITS_SIMPLE_HEADER)
     {
	JDFits_Keyword_Type *kwmax = ht->keys + ht->num_keywords;
	ht->pcount = ht->gcount = -1;
	
	while (kwt < kwmax)
	  {
	     if (kwt->type & JDFITS_NUMBER_MASK)
	       {
		  if ((ht->pcount == -1) && !strcmp (kwt->name, "PCOUNT"))
		    {
		       if ((-1 == parse_int (kwt))
			   || ((ht->pcount = kwt->v.ival) < 0))
			 {
			    jdfits_error ("PCOUNT card is defective.");
			    goto error_return;
			 }
		       if (ht->gcount != -1) break;
		    }
		  else if ((ht->gcount == -1) && !strcmp (kwt->name, "GCOUNT"))
		    {
		       if ((-1 == parse_int (kwt))
			   || ((ht->gcount = kwt->v.ival) <= 0))
			 {
			    jdfits_error ("GCOUNT card is defective.");
			    goto error_return;
			 }
		       if (ht->pcount != -1) break;
		    }
	       }
	     kwt++;
	  }
	if (ht->pcount == -1)
	  {
	     jdfits_warning ("PCOUNT card is missing.  Assuming a value of zero.");
	     ht->pcount = 0;
	  }
	
	if (ht->gcount == -1)
	  {
	     jdfits_warning ("GCOUNT card is missing.  Assuming a value of one.");
	     ht->gcount = 1;
	  }
     }
   else 
     {
	ht->pcount = 0;
	ht->gcount = 1;
     }
   
   ht->size = ht->gcount * (ht->pcount + naxis_product);
   ht->size *= abs (bitpix) / 8;	       /* convert to 8 bit BYTES */
   
   return 0;
   
   error_return:
   
   SLFREE (ht->keys);
   ht->keys = NULL;
   ht->name = NULL;
   ht->type = 0;
   return -1;
}

   
	
static void free_fits_header (JDFits_Type *ft)
{
   JDFits_Header_Type *h = ft->header;
   
   if (h == NULL) return;
   if (h->free_routine != NULL) (*h->free_routine) (h);
   jdfits_free ((char *) h->header_data_buf);
   jdfits_free ((char *) h->keys);
   jdfits_free ((char *) h);
   ft->header = NULL;
}


static int jdfits_read_record (JDFits_Type *ft, unsigned char *buf)
{
   int n;
   if ((ft->mode != JDFITS_READ_MODE) || (ft->fp == NULL))
     {
	jdfits_error ("jdfits_read_record: File not open for read.");
	return -1;
     }
	
   n = fread (buf, 1, JDFITS_RECORD_SIZE, ft->fp);
   if (n == 0) return -1;
   
   if (JDFITS_RECORD_SIZE != n)
     {
	jdfits_error ("jdfits_read_record: Unable to read whole record. Corrupt file?");
	return -1;
     }
   return 0;
}

int jdfits_parse_header (JDFits_Type *ft)
{
   if (-1 == parse_headers (ft->header))
     {
	jdfits_free ((char *) ft->header->header_data_buf);
	jdfits_free ((char *) ft->header);
	ft->header = NULL;
	return -1;
     }
   return 0;
}

   
/* This function simply reads the header at the current position of the file
 * pointer.  It makes no attempt to parse it.  The current header attached
 * to the @JDFits_Type@ structure is freed.  If the calling routine wants to
 * save it, it should set the @ft->header@ field to NULL.  When the routine
 * returns, it sets the fields @ft->header@, @ft->header->num_keywords@ 
 * and @ft->header->header_data@ structure fields appropriately.
 * It returns 0 upon success and -1 upon failure.
 */
static int jdfits_read_raw_header (JDFits_Type *ft)
{
   unsigned char *buf, *bufp, *bufp_max;
   unsigned int nrecords, ncards;
   unsigned char *header_data_buf;

   /* free up the current header */
   if (ft->header != NULL) free_fits_header (ft);
   
   if ((NULL == (ft->header = (JDFits_Header_Type *) jdfits_malloc (sizeof (JDFits_Header_Type))))
       || (NULL == (buf = (unsigned char *) jdfits_malloc (JDFITS_RECORD_SIZE))))
     {
	if (ft->header != NULL) SLFREE (ft->header);
	ft->header = NULL;
	jdfits_error ("jdfits_read_header: Memory allocation failure.");
	return -1;
     }
   memset ((char *) ft->header, 0, sizeof (JDFits_Header_Type));

   nrecords = 0;
   bufp = buf;
   ncards = 0;
   
   while (1)
     {
	if (-1 == jdfits_read_record (ft, bufp))
	  {
	     if (nrecords != 0) jdfits_error ("jdfits_read_header: Error reading record.");
	     SLFREE (buf);
	     return -1;
	  }
	
	nrecords++;
	
	/* check to see if this has an end card. */
	bufp_max = bufp + JDFITS_RECORD_SIZE;
	while (bufp < bufp_max)
	  {
	     if (!strncmp ((char *) bufp, "END     ", 8)) break;	     
	     ncards++;
	     bufp += JDFITS_CARD_SIZE;
	  }
	
	if (bufp < bufp_max) 
	  {
	     bufp += 8;
	     while (bufp < bufp_max)
	       {
		  if (*bufp != ' ')
		    {
		       jdfits_warning ("\
The record following the END card is not padded with blanks as required.");
		       break;
		    }
		  bufp++;
	       }
	     break;    /* END found */
	  }
	
	
	/* Not found, read another record */
	
	if (NULL == (buf = (unsigned char *) SLREALLOC (buf, (nrecords + 1) * JDFITS_RECORD_SIZE)))
	  {
	     jdfits_error ("jdfits_read_header: Memory RE-allocation failure.");
	     return -1;
	  }
	
	bufp = buf + nrecords * JDFITS_RECORD_SIZE;
     }

   header_data_buf = (unsigned char *) jdfits_malloc (ncards * (JDFITS_CARD_SIZE+1));
   if (header_data_buf == NULL)
     {
	jdfits_free ((char *) ft->header); ft->header = NULL;
	jdfits_free ((char *) buf);
	return -1;
     }
   
   ft->header->header_data_buf = header_data_buf;
   ft->header->num_keywords = ncards;

   bufp = buf;
   while (ncards)
     {
	memcpy ((char *)header_data_buf, (char *) bufp, JDFITS_CARD_SIZE);
	bufp += JDFITS_CARD_SIZE;
	header_data_buf += JDFITS_CARD_SIZE;
	*header_data_buf++ = 0;
	ncards--;
     }
   jdfits_free ((char *) buf);
   return 0;
}

int jdfits_read_header (JDFits_Type *ft)
{
   if (-1 == jdfits_read_raw_header (ft))
     return -1;
   
   return jdfits_parse_header (ft);
}



JDFits_Type *jdfits_open_file (char *file, int mode)
{
   FILE *fp;
   JDFits_Type *ft;
   
   if (NULL == (ft = (JDFits_Type *) jdfits_malloc(sizeof (JDFits_Type))))
     {
	return NULL;
     }
   memset ((char *) ft, 0, sizeof (JDFits_Type));
   
   switch (mode)
     {
      case JDFITS_READ_MODE:
	if (NULL != (fp = fopen (file, "rb")))
	  {
	     char buf[JDFITS_RECORD_SIZE];
	     if (JDFITS_RECORD_SIZE != fread (buf, 1, JDFITS_RECORD_SIZE, fp))
	       {
		  jdfits_error ("\
Unable to read %d bytes.\n\
 File %s does not appear to be a fits file.",
			      JDFITS_RECORD_SIZE, file);
		  fclose (fp);
		  fp = NULL;
	       }
	     else if (strncmp (buf, "SIMPLE  =", 9))
	       {
		  jdfits_error ("\
File %s lacks the SIMPLE keyword.  It is NOT a fits file.",
			      file);
		  fclose (fp);
		  fp = NULL;
	       }
	     else rewind (fp);
	  }
	break;
	
      case JDFITS_WRITE_MODE:
	if (NULL == (fp = fopen (file, "wb"))) break;
	if (NULL == (ft->write_buffer = (unsigned char *)jdfits_malloc(JDFITS_RECORD_SIZE)))
	  {
	     fclose (fp);
	     fp = NULL;
	  }
	ft->write_buffer_len = 0;
	break;
	
      default:
	jdfits_error ("Unable to open %s: Invalid mode: %d.", file, mode);
	fp = NULL;
     }
   
   if (fp == NULL)
     {
	SLFREE (ft);
	return NULL;
     }
   
   
   ft->fp = fp;
   ft->mode = mode;
   
   if (mode == JDFITS_READ_MODE)
     {
	if (-1 == jdfits_read_raw_header (ft))
	  {
	     jdfits_error ("Error reading header for %s.", file);
	     jdfits_close_file (ft);
	     return NULL;
	  }
	
	/* Is this a fits file?? */
	if (strncmp ("SIMPLE  =", (char *) ft->header->header_data_buf, 8)
	    || (-1 == jdfits_parse_header (ft))
	    || ((ft->header->keys[0].type & JDFITS_BOOL_TYPE) == 0))
	  {
	     jdfits_error ("%s is not a fits file.", file);
	     jdfits_close_file (ft);
	     return NULL;
	  }
	
	if (ft->header->keys[0].v.ival == 0)
	  {
	     jdfits_warning ("Fits file %s may not conform to the standard.", file);
	  }
     }

   return ft;
}




int jdfits_close_file (JDFits_Type *ft)
{
   if (ft->fp == NULL)
     {
	jdfits_error ("File is not open.");
	return -1;
     }
   
   switch (ft->mode)
     {
      case JDFITS_READ_MODE:
	fclose (ft->fp);
	ft->fp = NULL;
	break;
	
      case JDFITS_WRITE_MODE:
	if (ft->write_buffer_len)
	  {
	     jdfits_warning ("\
Either the header or the data portion of the file was not ended.  I am\n\
 assuming this is the data portion and filling it with 0.");
	     jdfits_flush_output (ft, 0);
	  }
	
	if (EOF == fclose (ft->fp))
	  {
	     jdfits_error ("Error closing file.  Check disk space.");
	     return -1;
	  }
	if (ft->write_buffer != NULL) SLFREE (ft->write_buffer);
	ft->write_buffer = NULL;
	ft->fp = NULL;
	break;
	
      default:
	jdfits_error ("Error closing file: Unknown mode.");
	return -1;   
     }
   
   free_fits_header (ft);
   SLFREE (ft);
   return 0;
}

/* The assumption made by this routine is that the header has been read but
 * the data for the header has not.
 */
int jdfits_skip_to_next_header (JDFits_Type *ft)
{
   JDFits_Header_Type *h = ft->header;
   off_t size, dsize;
   
   if (h == NULL)
     {
	jdfits_error ("jdfits_skip_to_next_header: Current header is NULL");
	return -1;
     }
   
   if ((ft->mode & JDFITS_READ_MODE) == 0)
     {
	jdfits_error ("jdfits_skip_to_next_header: file not open for read access.");
	return -1;
     }
   
   size = h->size;
   dsize = size % JDFITS_RECORD_SIZE;
   if (dsize) size += JDFITS_RECORD_SIZE - dsize;
   
#ifndef SEEK_CUR
# define SEEK_CUR 1
#endif
   if (-1 == FSEEK (ft->fp, size, SEEK_CUR))
     {
	jdfits_error ("jdfits_skip_to_next_header: seek error.");
	return -1;
     }
   
   free_fits_header (ft);
   return 0;
}


static int parse_number (JDFits_Keyword_Type *k)
{
   char *p, ch, *pmax;
   double x, y;
   int n, ix, iy;
   p = k->v.sval;
   pmax = k->comment;
   
   ch = *p;
   if ((ch == '-') || (ch == '+')) p++;
   while ((p < pmax) && isdigit (*p)) p++;
   
   if ((*p == '.') || (*p == 'E') || (*p == 'e'))
     {
	/* assume double -- worry about complex later */
	/* This is an _ugly_ hack.  The sscanf may pick up an integer on
	 * the next row if the keyword starts with an integer and this 
	 * card contains no comment.  Use the fact that k->name + 80 points
	 * to the next card.  This is _ugly_!!!
	 * FIXME!!!!!!!
	 */
	char ugly_ch = k->name[79];	       /* last character on this card */
	k->name[79] = 0;
	n = sscanf ((char *) k->v.sval, "%lf %lf", &x ,&y);
	k->name[79] = ugly_ch;

	if (n == 2)
	  {
	jdfits_warning ("\
parse_number: complex float type not supported.  For now, I will assume\n\
the other value is a comment. (Keyword = %s)",
		      k->name);
	  }
	if (n == 0)
	  {
	     jdfits_error ("Error parsing card (%s) as float.", k->name);
	     return -1;
	  }
	k->v.dval = x;
	k->type = JDFITS_FLOAT64_TYPE;
	while ((p < pmax) && (*p != '/') && (*p != ' ')) p++;
	handle_comment (k, p, pmax);
	return 0;
     }
   
   n = sscanf ((char *) k->v.sval, "%d %d", &ix ,&iy);
   if (n == 2)
     {
	jdfits_warning ("\
parse_number: complex integer type not supported.  For now, I will assume\n\
the other value is a comment. (Keyword = %s)",
		      k->name);
     }
   return parse_int (k);
}



int jdfits_parse_key (JDFits_Keyword_Type *k, unsigned int type)
{
   int ret;
   
   if ((k->type & JDFITS_UNPARSED_TYPE) == 0) 
     {
	if (type & k->type) return 0;
	return -1;
     }
#if 1
   if (type == JDFITS_ALL_TYPES)
     type = k->type;
#else
   type = k->type;
#endif
   
   switch (type & JDFITS_ALL_TYPES)
     {
      case 0:
      case JDFITS_NUMBER_MASK:
	ret = parse_number (k);
	break;
	
      case JDFITS_INT_TYPE:
	ret = parse_int (k);
	break;
	
      case JDFITS_FLOAT32_TYPE:
	ret = parse_float32 (k);
	break;
	
      case JDFITS_INT32_TYPE:
	ret = parse_int32 (k);
	break;
	
      case JDFITS_BOOL_TYPE:
	ret = parse_bool (k);
	break;
	
      case JDFITS_FLOAT64_TYPE:
	ret = parse_float64 (k);
	break;
	
      case JDFITS_INT16_TYPE:
	ret = parse_int16 (k);
	break;
	
      case JDFITS_COMMENT_TYPE:
	ret = parse_comment (k);
	break;
	
      case JDFITS_STRING_TYPE:
	ret = parse_string (k);
	break;
	
      default:
	jdfits_error ("jdfits_parse_key: unknown type: 0x%X.", type);
	ret = -1;
     }
   return ret;
}

/* This function does not parse the keyword. */
JDFits_Keyword_Type *_jdfits_find_keyword (JDFits_Header_Type *h, char *name)
{
   JDFits_Keyword_Type *k, *kmax;
   char ch;

      
   if (h == NULL)
     {
	jdfits_error ("_jdfits_find_keyword: header is NULL");
	return NULL;
     }

   k = h->keys;
   kmax = k + h->num_keywords;
   ch = *name++;

   ch |= 0x20;
   while (k < kmax)
     {
	if ((ch == (*k->name | 0x20))
	    && (0 == jdfits_strcasecmp (name, k->name + 1)))
	  return k;
	
	k++;
     }
   return NULL;
}

/* This returns NULL if the keyword does not exist. */
JDFits_Keyword_Type *jdfits_parse_keyword (JDFits_Header_Type *h, char *name, unsigned int type)
{
   JDFits_Keyword_Type *k;
   
   if (NULL == (k = _jdfits_find_keyword (h, name)))
     return NULL;

   if (-1 == jdfits_parse_key (k, type))
     {
	jdfits_error ("jdfits_parse_keyword: Error parsing %s.", k->name);
	return NULL;
     }
   return k;
}

