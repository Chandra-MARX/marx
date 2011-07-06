#include "config.h"
#include "marx-feat.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>

#include <jdmath.h>

#include "marx.h"
#include "_marx.h"

struct _Marx_QE_Type
{
   unsigned int num_refs;
   unsigned int num_energies;
   float *energies;
   float *qes;
};

void _marx_qe_free (Marx_QE_Type *qeinfo)
{
   if (qeinfo == NULL)
     return;

   if (qeinfo->num_refs > 1)
     {
	qeinfo->num_refs--;
	return;
     }
   /* NULLs ok */
   marx_free ((char *)qeinfo->energies);
   marx_free ((char *)qeinfo->qes);
   marx_free ((char *)qeinfo);
}

static Marx_QE_Type *alloc_ccd_qe_type (unsigned int num_energies)
{
   Marx_QE_Type *qeinfo;

   if (NULL == (qeinfo = (Marx_QE_Type *)marx_malloc (sizeof(Marx_QE_Type))))
     return NULL;
   memset ((char *)qeinfo, 0, sizeof(Marx_QE_Type));

   qeinfo->num_refs = 1;
   qeinfo->num_energies = num_energies;
   if ((NULL == (qeinfo->energies = (float *)marx_malloc(num_energies*sizeof(float))))
       || (NULL == (qeinfo->qes = (float *)marx_malloc(num_energies*sizeof(float)))))
     {
	_marx_qe_free (qeinfo);
	return NULL;
     }
   return qeinfo;
}

Marx_QE_Type *_marx_qe_read_file (char *file, char *ext, char *encol, char *qecol, char *filtercol)
{
   Marx_QE_Type *qeinfo = NULL;
   JDFits_Type *ft;
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   unsigned int i, num_rows, num_cols;
   char *enbuf = NULL, *qebuf = NULL, *filter_buf=NULL;
   char *cols[3];

   ft = jdfits_open_binary_table (file, ext);
   if (ft == NULL)
     {
	marx_error ("Unable to open QE file %s[%s]", file, ext);
	return NULL;
     }

   if ((NULL == (enbuf = (char *)marx_malloc (3 + strlen(encol))))
       || (NULL == (qebuf = (char *)marx_malloc (3 + strlen(qecol))))
       || ((filtercol != NULL)
	   && (NULL == (filter_buf = (char *)marx_malloc (3 + strlen(filtercol))))))
     goto return_error;

   num_cols = 0;
   (void) sprintf (enbuf, "f:%s", encol);
   (void) sprintf (qebuf, "f:%s", qecol);
   cols[num_cols++] = enbuf;
   cols[num_cols++] = qebuf;
   if (filtercol != NULL)
     {
	(void) sprintf (filter_buf, "f:%s", filtercol);
	cols[num_cols++] = filter_buf;
     }

   r = jdfits_bintable_aopen_rows (ft, num_cols, cols);
   if (r == NULL)
     goto return_error;

   num_rows = r->num_rows;
   if (num_rows < 2)
     {
	marx_error ("Expecting more than 2 rows in %s", file);
	goto return_error;
     }

   if (NULL == (qeinfo = alloc_ccd_qe_type (num_rows)))
     goto return_error;

   c = r->col_data;
   for (i = 0; i < num_rows; i++)
     {
	float qe;

	if (1 != jdfits_read_next_row (ft, r))
	  {
	     marx_error ("File %s appears to be corrupt", file);
	     goto return_error;
	  }
	qeinfo->energies[i] = c[0].data.f[0];
	qe = c[1].data.f[0];
	if (filtercol != NULL)
	  qe *= c[2].data.f[0];

	qeinfo->qes[i] = qe;
     }

   goto free_and_return;

return_error:

   marx_error ("Error processing %s", file);
   /* NULLs ok */
   _marx_qe_free (qeinfo);
   qeinfo = NULL;

   /* drop */

free_and_return:

   jdfits_bintable_close_rows (r);
   (void) jdfits_close_file (ft);
   marx_free (qebuf);
   marx_free (enbuf);
   return qeinfo;
}

double _marx_qe_compute (Marx_QE_Type *qeinfo, double energy)
{
   if (qeinfo == NULL)
     return 0.0;

   return JDMinterpolate_f ((float) energy, qeinfo->energies, qeinfo->qes, qeinfo->num_energies);
}

void _marx_qe_inc_ref (Marx_QE_Type *qeinfo)
{
   if (qeinfo != NULL)
     qeinfo->num_refs++;
}
