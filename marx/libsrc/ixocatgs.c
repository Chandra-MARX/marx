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

#include "ixoccd.h"
#include "marx.h"
#include "_marx.h"

typedef struct
{
   unsigned int num_refs;
   int min_order, max_order;
   unsigned int num_orders;	       /* max_order-min_order+1 */
   float *energies;
   unsigned int num_energies;
   float **cum_efficiencies;	       /* [num_energies][num_orders] */
}
Grating_Eff_Type;

typedef struct _Grating_Type
{
   struct _Grating_Type *next;
   double period;
   double theta_blur;
   double dp_over_p;
   Grating_Eff_Type *geff;

   /* n, l, and d are the grating vectors.  n is normal to the
    * facet, l is in the direction of the lines, and d is in the is normal
    * to the lines.
    */
   JDMVector_Type n, l, d;

   JDMVector_Type origin;	       /* at the center of the facet */
   JDMVector_Type xhat, yhat, zhat;    /* facet basis vectors */
   double xlen, ylen;
   struct _Grating_Type *support_structure;
}
Grating_Type;

typedef struct
{
   Grating_Type *gratings;
   double xmin, xmax, ymin, ymax, zmin, zmax;   /* bounding box */
   double dispersion_angle;
}
Grating_Module_Type;

static int Use_Finite_Facets = 1;
static int Use_Support_Structure = 0;

#define SECTOR_SIZE (30.0/2)	       /* degrees */
#define MIN_SECTOR_RADIUS 200.0	       /* mm */
#define MAX_SECTOR_RADIUS 800.0	       /* mm */

static void free_grating_eff_type (Grating_Eff_Type *geff)
{
   if (geff == NULL)
     return;

   if (geff->num_refs > 1)
     {
	geff->num_refs--;
	return;
     }

   if (geff->cum_efficiencies != NULL)
     JDMfree_float_matrix (geff->cum_efficiencies, geff->num_energies);

   if (geff->energies != NULL)
     JDMfree_float_vector (geff->energies);

   marx_free ((char *) geff);
}

static Grating_Eff_Type *alloc_grating_eff_type (unsigned num_energies, int min_order, int max_order)
{
   Grating_Eff_Type *geff;

   if (NULL == (geff = (Grating_Eff_Type *)marx_malloc (sizeof (Grating_Eff_Type))))
     return NULL;

   memset ((char *)geff, 0, sizeof(Grating_Eff_Type));

   geff->num_refs = 1;
   geff->min_order = min_order;
   geff->max_order = max_order;
   geff->num_energies = num_energies;
   geff->num_orders = max_order - min_order + 1;

   geff->energies = JDMfloat_vector (num_energies);
   if (geff->energies == NULL)
     {
	free_grating_eff_type (geff);
	return NULL;
     }

   if (NULL == (geff->cum_efficiencies = JDMfloat_matrix (num_energies, geff->num_orders)))
     {
	free_grating_eff_type (geff);
	return NULL;
     }

   return geff;
}

static Grating_Module_Type The_Left_Gratings;
static Grating_Module_Type The_Right_Gratings;
static double Left_Dispersion_Angle = 2.0;
static double Right_Dispersion_Angle = -2.0;

static double Rowland_Theta = 1.5, Rowland_Phi = 4.5;
static double Rowland_Sin_Phi, Rowland_Cos_Phi, Rowland_Tan_Phi;
static double Rowland_R, Rowland_A, Rowland_C;
static double Rowland_A2, Rowland_C2, Rowland_R2;

static int diffract_photon (Grating_Type *g,
			    Marx_Photon_Attr_Type *at,
			    int order)
{
   JDMVector_Type p, x;
   JDMVector_Type n, l, d, dp;
   double factor;
   double n_lambda_over_d, dp_over_p;
   double p_n, p_d, p_l;

   n_lambda_over_d = order * (2.0 * PI * HBAR_C) / g->period / at->energy;

   x = at->x;
   p = at->p;

   n = g->n;
   l = g->l;
   d = g->d;

   /* Apply a statistical blur to represent small misalignments */
   if (g->theta_blur != 0.0)
     {
	JDMVector_Type l_tmp, d_tmp;
	double c, s, theta;

	theta = g->theta_blur * JDMgaussian_random ();
	/* Now rotate about n axis by theta. */
	c = cos(theta);
	s = sin(theta);
	l_tmp = l;
	d_tmp = d;
	l = JDMv_ax1_bx2 (c, l_tmp, s, d_tmp);
	d = JDMv_ax1_bx2 (-s, l_tmp, c, d_tmp);
     }

   p_d = n_lambda_over_d + JDMv_pdot_prod (&p, &d);
   p_l = JDMv_pdot_prod (&p, &l);
   p_n = 1.0 - p_l * p_l - p_d * p_d;

   if (p_n < 0.0)
     return -1;

   p_n = sqrt (p_n);

   p = JDMv_ax1_bx2 (p_d, d, p_n, n);
   at->p = JDMv_ax1_bx2 (1.0, p, p_l, l);

   /* Now apply dp/p blur */
   dp_over_p = g->dp_over_p * JDMgaussian_random ();

   if (dp_over_p == 0)
     return 0;

   factor = n_lambda_over_d * dp_over_p;

   dp = JDMv_ax1_bx2 (-factor, d, factor * (p_d / p_n), n);

   at->p.x += dp.x;
   at->p.y += dp.y;
   at->p.z += dp.z;

   at->p = JDMv_unit_vector (at->p);
   return 0;
}

static int newtons_quartic (double a, double b, double c, double d,
			    double t0, double *tp)
{
   unsigned int max_it = 10;
   double a2, a3, b2;
   double t;
   double eps = 1.0e-6;
   double num, den;

   a2 = 2.0 * a;
   a3 = 3.0 * a;
   b2 = 2.0 * b;

   while (1)
     {
	num = -d + t0*(t0*(b + t0*(a2 + t0*3.0)));
	den = c + t0*(b2 + t0*(a3 + t0*4.0));
	t = num / den;

	if (fabs (t - t0) < eps) break;

	max_it--;
	if (max_it == 0) return -1;

	t0 = t;
     }

   *tp = t;
   return 0;
}

static void compute_grating_vectors (Grating_Type *g, JDMVector_Type *xp)
{
   JDMVector_Type n, l, d;

   /* Here we assume that the facets are
    * perfect in the sense that the facet normal lies in the direction of
    * the origin.
    */
   /* n = -x/|x| */
   n = JDMv_unit_vector (*xp);
   n.x = -n.x; n.y = -n.y; n.z = -n.z;

   /* Choose l to be perp to both the axis of the torus and n.  That is
    *   l.n = 0 and l.a = 0.
    * ==>
    *   lx*nx + ly*ny + lz*nz = 0
    *   lx*ax + ly*ay + lz*az = 0
    *
    *   When ay = 0:
    *
    *    lz = -lx*ax/az
    * ==>
    *    lx*nx + ly*ny - lx*nz*ax/az = 0
    * ==>
    *    lx = (-ly*ny)/(nx - nz*ax/az)
    */
   l.y = -1.0;
   l.x = (-l.y*n.y)/(n.x - n.z*Rowland_Tan_Phi);
   l.z = -l.x*Rowland_Tan_Phi;

   l = JDMv_unit_vector (l);

   /* Finally rotate about l such that the angle between the normal and the
    * direction of the ray is equal to the blaze angle, which is given by
    * Rowland_Theta(=t):
    *
    *    p.n' = cos(t)
    *
    * Currently, p.n=1.
    */
   n = JDMv_rotate_unit_vector (n, l, -Rowland_Theta);
   d = JDMv_cross_prod (n, l);

   g->n = n;
   g->l = l;
   g->d = d;
}

static int intersect_torus (JDMVector_Type *xp, JDMVector_Type *pp)
{
   JDMVector_Type p, x;
   double t, a, b, c, d;
   double pA, x0A, x2, Ax0A, r2x0A;

   p = *pp;
   x = *xp;

   /* project the ray from x0 to x such that p.x = 0.
    * We have x=x0+pt ==> p.x = p.x0 + t ==> t = -p.x0
    * ==> x = x0 - (p.x0)p
    */
   t = -(p.x*x.x + p.y*x.y + p.z*x.z);
   x.x += p.x*t;
   x.y += p.y*t;
   x.z += p.z*t;

   pA = p.x*Rowland_Sin_Phi + p.z*Rowland_Cos_Phi;
   x0A = x.x*Rowland_Sin_Phi + x.z*Rowland_Cos_Phi;
   x2 = x.x*x.x + x.y*x.y + x.z*x.z;
   Ax0A = Rowland_A * x0A;
   r2x0A = Rowland_R2*x0A;

   a = -4.0*Rowland_A*pA;
   b = 2.0*x2 + 4.0*(Rowland_R2*pA*pA - Ax0A - Rowland_C2);
   c = 4.0*pA*(2*r2x0A - Rowland_A*x2);
   d = x2*(x2 - 4.0*(Ax0A + Rowland_C2)) + 4.0*r2x0A*x0A;

   t = 2.0*Rowland_A*pA - 2*Rowland_C*sqrt((1.0-pA)*(1.0+pA));

   if (-1 == newtons_quartic (a, b, c, d, t, &t))
     return -1;

   x.x += p.x*t;
   x.y += p.y*t;
   x.z += p.z*t;

   t = d + t*(c + t*(b + t*(a + t)));

   *xp = x;
   return 1;
}

static Grating_Type *intersect_facets (JDMVector_Type *xp, JDMVector_Type *pp)
{
   Grating_Type *g;
   Grating_Module_Type *module;
   double y, z, t;
   if (Use_Finite_Facets == 0)
     {
	if (1 != intersect_torus (xp, pp))
	  return NULL;

	g = (xp->y < 0) ? The_Left_Gratings.gratings : The_Right_Gratings.gratings;
	compute_grating_vectors (g, xp);

	return g;
     }

   module = (xp->y < 0) ? &The_Left_Gratings : &The_Right_Gratings;
   /* Project the ray to the plane of the bounding box closest to the mirror.
    * If the position does not fall within the bbbox there, reject it.
    * It is possible wild ray to still intersect the grating module from
    * entering the side of the bbox, but that seems extremely unlikely.
    *
    * X = X0 + pt
    * x = x0 + px*t ==> t = (x-x0)/px
    * y = y0 + py/px*(x-x0)
    * z = z0 + pz/px*(x-x0)
    */
   t = (module->xmax - xp->x)/pp->x;
   y = xp->y + pp->y*t;
   if ((y < module->ymin) || (y > module->ymax))
     return NULL;
   z = xp->z + pp->z*t;
   if ((z < module->zmin) || (z > module->zmax))
     return NULL;

   g = module->gratings;

   while (g != NULL)
     {
	JDMVector_Type normal, r;
	double pdotn;
	double rx, ry;

	/* Intersect with the plane.  Let X' be the intersection point:
	 *   X' = X0 + pt
	 * where (X'-O).n = 0 and O is the origin.
	 * ==> 0 = (X0 + pt - O).n
	 * ==> t = -(X0-O).n / p.n
	 * X'-O = X0-O + pt
	 *      = r + pt
	 *      = r - p*r.n/p.n
	 */
	normal = g->zhat;
	pdotn = JDMv_pdot_prod (pp, &normal);

	if (pdotn == 0)
	  {
	     g = g->next;
	     continue;
	  }

	r = JDMv_diff (*xp, g->origin);

	r = JDMv_pax1_bx2 (1.0, &r,
			   -1.0 * JDMv_pdot_prod (&r, &normal)/pdotn, pp);
	/* This r is defined in the plane of the facet with origin
	 * at g->origin.
	 */

	/* Now check to see if the coordinates in the plane
	 * is within the boundaries of the chip.  Here we assume the
	 * facet is rectangular.
	 */

	rx = JDMv_dot_prod (r, g->xhat);
	if (fabs(rx) > 0.5*g->xlen)
	  {
	     g = g->next;
	     continue;
	  }

	ry = JDMv_dot_prod (r, g->yhat);
	if (fabs(ry) > 0.5*g->ylen)
	  {
	     g = g->next;
	     continue;
	  }

	*xp = JDMv_sum (r, g->origin);
	return g;
     }

   /* Missed */
   return NULL;
}

static int compute_diffraction_order (double energy, Grating_Eff_Type *geff, SIGNED_CHAR *orderp)
{
   unsigned int n0, n1;
   double w0, w1, r;
   unsigned int i, num_orders;
   float *eff0, *eff1;
   double cumsum;

   n0 = JDMbinary_search_f ((float) energy, geff->energies, geff->num_energies);
   n1 = n0 + 1;
   if (n1 == geff->num_energies)
     {
	n0--; n1--;
     }
   w1 = (energy - geff->energies[n0])/(geff->energies[n1] - geff->energies[n0]);
   w0 = (1-w1);

   r = JDMrandom ();

   eff0 = geff->cum_efficiencies[n0];
   eff1 = geff->cum_efficiencies[n1];
   num_orders = geff->num_orders;

   cumsum = 0.0;

   for (i = 0; i < num_orders; i++)
     {
	double c;

	c = w0*eff0[i] + w1*eff1[i];
	if (r < c)
	  {
	     *orderp = geff->min_order + i;
	     return 0;
	  }
     }
   return -1;
}

int _marx_catgs_diffract (Marx_Photon_Type *pt)
{
   Marx_Photon_Attr_Type *photon_attributes, *at;
   unsigned int num_sorted, i, *sorted_index;

   if (pt->history & MARX_ORDER_OK)
     return 0;
   pt->history |= MARX_ORDER_OK;

   marx_prune_photons (pt);
   num_sorted = pt->num_sorted;
   photon_attributes = pt->attributes;
   sorted_index = pt->sorted_index;

   for (i = 0; i < num_sorted; i++)
     {
	double theta;
	Grating_Type *g;

	at = photon_attributes + sorted_index[i];

	/* at->flags |= PHOTON_UNDIFFRACTED; */
	at->order = 0;

	if (Use_Finite_Facets == 0)
	  {
	     if ((at->mirror_shell < 73) || (at->mirror_shell > 175))
	       continue;

	     theta = atan2 (at->x.z, at->x.y);   /* -PI <= theta <= PI */
	     theta *= 180.0/PI;

	     theta = fabs(theta);
	     if ((theta > SECTOR_SIZE) && (theta < 180.0-SECTOR_SIZE))
	       continue;
	  }

	g = intersect_facets (&at->x, &at->p);
	if (g == NULL)
	  continue;

	if (0)
	  {
	     double len;
	     JDMVector_Type dx, x1 = at->x;
	     (void) intersect_torus (&at->x, &at->p);
	     dx = JDMv_diff (x1, at->x);
	     len = JDMv_length (dx);
	     if (1)
	       {
	       }
	  }

	if (-1 == compute_diffraction_order (at->energy, g->geff, &at->order))
	  {
	     at->flags |= PHOTON_UNDIFFRACTED;
	     continue;
	  }

	if (-1 == diffract_photon (g, at, at->order))
	  {
	     at->flags |= PHOTON_UNDIFFRACTED;
	     continue;
	  }

	if (NULL == (g = g->support_structure))
	  continue;

	pt->history |= MARX_ORDER1_OK;

	if (-1 == compute_diffraction_order (at->energy, g->geff, &at->support_orders[0]))
	  {
	     at->flags = PHOTON_UNDIFFRACTED;
	     continue;
	  }

	if (-1 == diffract_photon (g, at, at->support_orders[0]))
	  {
	     at->flags = PHOTON_UNDIFFRACTED;
	     continue;
	  }
     }

   return 0;
}

static double Rowland_Distance = 19700.0;
static double Left_Grating_Period, Right_Grating_Period;   /* microns */
static char *Left_Grating_Eff_File;
static char *Right_Grating_Eff_File;
static double Facet_Size = 60.0;
static double Finite_Facet_Center_Y = 1045.0;	       /* from raytrace */
static double dP_Over_P_Sigma = 0.0;

static double Support_Period = 5.0;    /* microns */

#define MIN_GEFF_ORDER -2
#define MAX_GEFF_ORDER 16
#define MAX_GEFF_COLUMNS (MAX_GEFF_ORDER-MIN_GEFF_ORDER+2)
static char *Geff_Columns [MAX_GEFF_COLUMNS] =
{
   "f:ENERGY",
   "f:EFFm2", "f:EFFm1", "f:EFF0",
   "f:EFF1", "f:EFF2", "f:EFF3", "f:EFF4",
   "f:EFF5", "f:EFF6", "f:EFF7",
   "f:EFF8", "f:EFF9", "f:EFF10",
   "f:EFF11", "f:EFF12", "f:EFF13"
   "f:EFF14", "f:EFF15", "f:EFF16"
};

static Grating_Eff_Type *read_geff_caldb_file (char *file)
{
   int min_order, max_order;
   unsigned int num_energies, num_orders = 0;
   JDFits_Type *f;
   char hduname[16];
   JDFits_Row_Type *r;
   JDFits_Col_Data_Type *c;
   unsigned int i;
   Grating_Eff_Type *geff = NULL;

   sprintf (hduname, "IXO_GREFF");
   marx_message ("\t%s[%s]\n", file, hduname);

   if (NULL == (f = _marx_open_binary_hdu (file, hduname)))
     return NULL;

   min_order = -2; max_order = 1;
   while (max_order < MAX_GEFF_ORDER)
     {
	char colnam [16];
	int status;

	sprintf (colnam, "EFF%d", max_order+1);
	status = jdfits_bintable_column_exists (f, colnam);
	if (status == -1)
	  goto return_error;
	if (status == 0)
	  break;
	max_order++;
     }
   num_orders = max_order - min_order + 1;

   r = jdfits_bintable_aopen_rows (f, 1+num_orders, Geff_Columns);

   if (r == NULL)
     goto return_error;

   num_energies = r->num_rows;
   c = r->col_data;
   if (c[0].repeat != 1)
     {
	marx_error ("The ENERGY column must be a scalar column");
	goto return_error;
     }

   if (NULL == (geff = alloc_grating_eff_type (num_energies, min_order, max_order)))
     goto return_error;

   for (i = 0; i < num_energies; i++)
     {
	unsigned int j;
	double cumsum;
	float *cumeffs;

	if (1 != jdfits_read_next_row (f, r))
	  {
	     marx_error ("Error reading row %u", i+1);
	     goto return_error;
	  }

	geff->energies[i] = c[0].data.f[0];

	cumsum = 0.0;
	cumeffs = geff->cum_efficiencies[i];
	for (j = 0; j < num_orders; j++)
	  {
	     cumsum += c[j+1].data.f[0];
	     cumeffs[j] = cumsum;
	  }
     }
   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);

   return geff;

return_error:

   jdfits_bintable_close_rows (r);
   jdfits_close_file (f);
   free_grating_eff_type (geff);

   return NULL;
}

static Grating_Eff_Type *read_geff_file (char *file)
{
   Grating_Eff_Type *geff;

   if (NULL == (file = marx_make_data_file_name (file)))
     return NULL;

   geff = read_geff_caldb_file (file);
   marx_free (file);
   return geff;
}

static void free_gratings (Grating_Type *g)
{
   while (g != NULL)
     {
	Grating_Type *next = g->next;
	free_grating_eff_type (g->geff);
	if (g->support_structure != NULL)
	  free_gratings (g->support_structure);
	marx_free ((char *)g);
	g = next;
     }
}

static Grating_Type *alloc_grating (Grating_Eff_Type *geff, double period)
{
   Grating_Type *g;

   if (NULL == (g = (Grating_Type *) marx_malloc (sizeof (Grating_Type))))
     return NULL;
   memset ((char *) g, 0, sizeof (Grating_Type));

   g->period = period;
   g->theta_blur = 0.0;
   g->dp_over_p = dP_Over_P_Sigma;

   geff->num_refs++;
   g->geff = geff;
   return g;
}

static Param_Table_Type IXO_CATGS_Parm_Table [] =
{
   {"IXO_Rowland_Theta",	PF_REAL_TYPE,	&Rowland_Theta},
   {"IXO_Rowland_Phi",		PF_REAL_TYPE,	&Rowland_Phi},
   {"IXO_Rowland_Distance",	PF_REAL_TYPE,	&Rowland_Distance},
   {"IXO_Left_Grating_Eff_File",PF_FILE_TYPE,	&Left_Grating_Eff_File},
   {"IXO_Right_Grating_Eff_File",PF_FILE_TYPE,	&Right_Grating_Eff_File},
   {"IXO_Left_Grating_Period",	PF_REAL_TYPE,	&Left_Grating_Period},
   {"IXO_Right_Grating_Period",	PF_REAL_TYPE,	&Right_Grating_Period},
   {"IXO_Left_Dispersion_Angle",PF_REAL_TYPE,	&Left_Dispersion_Angle},
   {"IXO_Right_Dispersion_Angle",PF_REAL_TYPE,	&Right_Dispersion_Angle},
   {"IXO_Grating_Facet_Size", PF_REAL_TYPE,	&Facet_Size},
   {"IXO_Grating_Facet_Y", 	PF_REAL_TYPE,	&Finite_Facet_Center_Y},
   {"IXO_CATGS_dPoverP",	PF_REAL_TYPE,	&dP_Over_P_Sigma},
   {NULL, 0, NULL}
};

static int Variables_Inited = 0;
static void _marx_catgs_init_variables (void)
{
   if (Variables_Inited)
     return;

   Rowland_Theta *= PI/180.0;
   Rowland_Phi *= PI/180.0;

   Rowland_R = 0.5*Rowland_Distance/cos(Rowland_Theta);
   Rowland_A = Rowland_R * sin(Rowland_Theta + Rowland_Phi);
   Rowland_C = Rowland_R * cos(Rowland_Theta + Rowland_Phi);
   Rowland_R2 = Rowland_R*Rowland_R;
   Rowland_A2 = Rowland_A*Rowland_A;
   Rowland_C2 = Rowland_C*Rowland_C;

   Rowland_Sin_Phi = sin (Rowland_Phi);
   Rowland_Cos_Phi = cos (Rowland_Phi);
   Rowland_Tan_Phi = tan (Rowland_Phi);

   Use_Finite_Facets = (Facet_Size > 0);
   Left_Dispersion_Angle *= PI/180.0;
   Right_Dispersion_Angle *= PI/180.0;

   Variables_Inited = 1;
}

static int read_ixo_catgs_parms (Param_File_Type *p)
{
   if (-1 == pf_get_parameters (p, IXO_CATGS_Parm_Table))
     return -1;

   Variables_Inited = 0;
   _marx_catgs_init_variables ();

   return 0;
}

static void compute_grating_module_bbox (Grating_Module_Type *module)
{
   double
     xmin = 1e30, xmax = -1e30,
     ymin = 1e30, ymax = -1e30,
     zmin = 1e30, zmax = -1e30;
   Grating_Type *g;
   int count = 0;

   g = module->gratings;

   while (g != NULL)
     {
	double x, y, z, xc, yc, zc, dx, dy,
	  xhat_x, xhat_y, xhat_z, yhat_x, yhat_y, yhat_z;
	int sx, sy;

	count++;
	xc = g->origin.x; yc = g->origin.y; zc = g->origin.z;
	dx = 0.5*g->xlen; dy = 0.5*g->ylen;

	xhat_x = g->xhat.x; xhat_y = g->xhat.y; xhat_z = g->xhat.z;
	yhat_x = g->yhat.x; yhat_y = g->yhat.y; yhat_z = g->yhat.z;

	for (sx = -1; sx < 2; sx += 2)
	  {
	     for (sy = -1; sy < 2; sy += 2)
	       {
		  x = xc + sx*dx*xhat_x + sy*dy*yhat_x;
		  y = yc + sx*dx*xhat_y + sy*dy*yhat_y;
		  z = zc + sx*dx*xhat_z + sy*dy*yhat_z;
		  if (x < xmin) xmin = x;
		  if (y < ymin) ymin = y;
		  if (z < zmin) zmin = z;
		  if (x > xmax) xmax = x;
		  if (y > ymax) ymax = y;
		  if (z > zmax) zmax = z;
	       }
	  }
	g = g->next;
     }
   module->xmin = xmin;
   module->xmax = xmax;
   module->ymin = ymin;
   module->ymax = ymax;
   module->zmin = zmin;
   module->zmax = zmax;
   fprintf (stderr, "%d grating\n", count);
}

static void rotate_vector (JDMVector_Type *v, double c, double s)
{
   double vy, vz;
   vy = c*v->y - s*v->z;
   vz = s*v->y + c*v->z;
   v->y = vy;
   v->z = vz;
}

static void rotate_grating (Grating_Type *g, double c, double s)
{
   /* The rotation is about the x axis using the right hand rule. */
   rotate_vector (&g->n, c, s);
   rotate_vector (&g->l, c, s);
   rotate_vector (&g->d, c, s);
   rotate_vector (&g->origin, c, s);
   rotate_vector (&g->xhat, c, s);
   rotate_vector (&g->yhat, c, s);
   rotate_vector (&g->zhat, c, s);
   if (g->support_structure != NULL)
     rotate_grating (g->support_structure, c, s);
}

static void rotate_grating_module (Grating_Module_Type *module, double dispersion_angle)
{
   double s, c;
   Grating_Type *g;

   s = sin(dispersion_angle);
   c = cos(dispersion_angle);
   g = module->gratings;
   while (g != NULL)
     {
	rotate_grating (g, c, s);
	g = g->next;
     }
}

static Grating_Type *
make_finite_facet_module (double period, Grating_Eff_Type *geff,
			  double center_y, double center_z)
{
   Grating_Type *gratings, *tail;
   double min_y, min_z;
   JDMVector_Type p;
   double dy = Facet_Size, dz = Facet_Size;
   double gap = 1.0;
   double rmin = MIN_SECTOR_RADIUS, rmax = MAX_SECTOR_RADIUS;
   int num_z = 12, num_y = 8;	       /* assumed even!!! */
   int i, j;
   double s;

   /* dy /= 12; dz /= 12; num_z *= 12; num_y *= 12; */
   (void) center_z;

   s = (center_y < 0) ? -1 : 1;

   min_z = -rmax*sin(PI/180.0*SECTOR_SIZE);
   num_z = 1 + (int) ((-2*min_z)/(dz+gap));
   min_y = (center_y < 0) ? -rmax : rmin;
   num_y = 1 + (int)((rmax-rmin)/(dy+gap));

   min_z += 0.5*dz;
   min_y += 0.5*dz;

   /*   +---+---+---+---+---+---+
    *   |   |   |   |   |   |   |
    *   +---+---+---+---+---+---+
    */
#if 0
   min_y = center_y - (num_y/2)*(dy + gap) + 0.5*dy;
   min_z = center_z - (num_z/2)*(dz + gap) + 0.5*dz;
#endif

   gratings = NULL;
   tail = NULL;

   for (i = 0; i < num_y; i++)
     {
	double y = min_y + i*(dy + gap);
	for (j = 0; j < num_z; j++)
	  {
	     Grating_Type *g;
	     JDMVector_Type x, x0, x1, x2, x3, x4;
	     double r, theta;
	     double z = min_z + j*(dz + gap);

	     /* If the center falls in the sector, accept the grating */
	     theta = fabs (atan2 (z, y) * (180.0/PI));
	     if ((theta > SECTOR_SIZE) && (theta < 180.0-SECTOR_SIZE))
	       continue;
	     r = hypot (y, z);
	     if ((r < rmin) || (r > rmax))
	       continue;

	     if (NULL == (g = alloc_grating (geff, period)))
	       {
		  free_gratings (gratings);
		  return NULL;
	       }
	     if (gratings == NULL)
	       gratings = g;
	     else
	       tail->next = g;
	     tail = g;

	     /* The facet is tangent */

	     /* Intersect the center and the 4 corners with the torus.  Take
	      * the weighted mean:
	      *    X = 1/2 * (Xcntr + 1/4 * sum_corners)
	      */
	     p.x = -1; p.y = 0; p.z = 0;
	     x0.x = Rowland_Distance; x0.y = y; x0.z = z;
	     x1 = JDMv_sum (x0, JDMv_ax1_bx2 (-0.5*dy, g->xhat, -0.5*dz, g->yhat));
	     x2 = JDMv_sum (x0, JDMv_ax1_bx2 (-0.5*dy, g->xhat, 0.5*dz, g->yhat));
	     x3 = JDMv_sum (x0, JDMv_ax1_bx2 (0.5*dy, g->xhat, -0.5*dz, g->yhat));
	     x4 = JDMv_sum (x0, JDMv_ax1_bx2 (0.5*dy, g->xhat,  0.5*dz, g->yhat));

	     if ((1 != intersect_torus (&x0, &p))
		 || (1 != intersect_torus (&x1, &p))
		 || (1 != intersect_torus (&x2, &p))
		 || (1 != intersect_torus (&x3, &p))
		 || (1 != intersect_torus (&x4, &p)))
	       {
		  marx_error ("Failed to intersect the grating facet with the torus!");
		  free_gratings (gratings);
		  return NULL;
	       }

	     x.x = 0.5*x0.x + 0.125*(x1.x + x2.x + x3.x + x4.x);
	     x.y = 0.5*x0.y + 0.125*(x1.y + x2.y + x3.y + x4.y);
	     x.z = 0.5*x0.z + 0.125*(x1.z + x2.z + x3.z + x4.z);
	     compute_grating_vectors (g, &x);
	     /* The facet normal is aligned with the grating normal.  For
	      * modules on the left/right, yhat is in the dispersion direction,
	      * which correponds to the marx z direction
	      */
	     g->zhat = g->n;
	     g->yhat = g->d;
	     g->xhat = g->l;
	     g->origin = x;
	     g->xlen = dy;
	     g->ylen = dz;
	  }
     }
   return gratings;
}

static Grating_Eff_Type *get_support_efficiencies (void)
{
   Grating_Eff_Type *geff;
   unsigned int num_energies = 2;
   int order, min_order = -1, max_order = 1;
   double dc, c;

   if (NULL == (geff = alloc_grating_eff_type (num_energies, min_order, max_order)))
     return NULL;
   geff->energies[0] = 0.001;
   geff->energies[1] = 20.0;

   dc = 1.0 / (max_order - min_order + 1);
   c = 0.0;

   for (order = min_order; order <= max_order; order++)
     {
	c += dc;
	geff->cum_efficiencies[0][order-min_order] = c;
	geff->cum_efficiencies[1][order-min_order] = c;
     }

   return geff;
}

static int add_support_structure (Grating_Type *gratings)
{
   Grating_Eff_Type *geff;
   double theta = PI/2.0;

   if (Use_Support_Structure == 0)
     return 0;

   if (NULL == (geff = get_support_efficiencies ()))
     return -1;

   while (gratings != NULL)
     {
	Grating_Type *g;

	if (NULL == (g = alloc_grating (geff, Support_Period)))
	  {
	     free_grating_eff_type (geff);
	     return -1;
	  }

	g->n = gratings->n;
	g->l = JDMv_rotate_unit_vector (g->n, gratings->l, theta);
	g->d = JDMv_rotate_unit_vector (g->n, gratings->d, theta);

	gratings->support_structure = g;
	gratings = gratings->next;
     }

   free_grating_eff_type (geff);
   return 0;
}

static int make_finite_facet_gratings (Grating_Eff_Type *left_geff, Grating_Eff_Type *right_geff)
{
   double center_y = Finite_Facet_Center_Y;
   double center_z = 0.0;
   Grating_Module_Type *module;

   module = &The_Left_Gratings;
   if (NULL == (module->gratings
		= make_finite_facet_module (Left_Grating_Period, left_geff,
					    -center_y, center_z)))
     return -1;
   module->dispersion_angle = Left_Dispersion_Angle;
   rotate_grating_module (module, module->dispersion_angle);
   compute_grating_module_bbox (module);

   module = &The_Right_Gratings;
   if (NULL == (module->gratings
		= make_finite_facet_module (Right_Grating_Period, right_geff,
					    center_y, center_z)))
     return -1;
   module->dispersion_angle = Right_Dispersion_Angle;
   rotate_grating_module (module, module->dispersion_angle);
   compute_grating_module_bbox (module);

   if (-1 == add_support_structure (The_Right_Gratings.gratings))
     return -1;

   if (-1 == add_support_structure (The_Left_Gratings.gratings))
     return -1;

   return 0;
}

int _marx_catgs_init (Param_File_Type *pf)
{
   Grating_Eff_Type *left_geff, *right_geff;
   int ret;

   marx_message ("Initializing CATGS...\n");

   if (-1 == read_ixo_catgs_parms (pf))
     return -1;

   if (NULL == (left_geff = read_geff_file (Left_Grating_Eff_File)))
     return -1;
   if (NULL == (right_geff = read_geff_file (Right_Grating_Eff_File)))
     {
	free_grating_eff_type (left_geff);
	return -1;
     }

   ret = -1;

   if (Use_Finite_Facets)
     {
	if (-1 == make_finite_facet_gratings (left_geff, right_geff))
	  goto free_and_return;
     }
   else
     {
	if (NULL == (The_Left_Gratings.gratings = alloc_grating (left_geff, Left_Grating_Period)))
	  goto free_and_return;
	The_Left_Gratings.dispersion_angle = Left_Dispersion_Angle;

	if (NULL == (The_Right_Gratings.gratings = alloc_grating (right_geff, Right_Grating_Period)))
	  goto free_and_return;
	The_Right_Gratings.dispersion_angle = Right_Dispersion_Angle;
     }

   ret = 0;
   /* drop */

free_and_return:
   free_grating_eff_type (left_geff);
   free_grating_eff_type (right_geff);

   return ret;
}

#include "ixoccd.c"

