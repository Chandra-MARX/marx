#!/usr/bin/env jdl-script
require ("fits");

if (__argc != 5)
{
   () = fprintf (stderr, "Usage: %s FEF CCDID CHIPX CHIPY\n", __argv[0]);
   exit (1);
}

_debug_info = 1; _traceback=1;
static variable File = __argv[1];
static variable Chip_Id = integer (__argv[2]);
static variable Chip_X = atof (__argv[3]);
static variable Chip_Y = atof (__argv[4]);

static define compute_fef (t, p, kth)
{
   variable i, imax = 10;
   
   variable y = Double_Type[length (kth),length(p)];
   
   for (i = 1; i <= imax; i++)
     {
	variable sigma = get_struct_field (t, sprintf ("g%d_fwhm", i))[kth]
	  * 0.4246609001440095;
	variable pos = get_struct_field (t, sprintf ("g%d_pos", i))[kth];
	variable amp = get_struct_field (t, sprintf ("g%d_ampl", i))[kth];
	variable k0=0;
	vmessage ("[%d %d]: fwhm=%g, pos=%g, ampl=%g (regnum=%d)", i, k0, sigma[k0], pos[k0], amp[k0],
		  t.regnum[kth[0]]);
	
	variable j;
	
	_for (0, length (kth)-1,1)
	  {
	     j = ();
	     if (sigma[j] != 0)
	       {
		  variable x = (p - pos[j])/sigma[j];
		  y[j,*] += amp[j] * exp (-0.5*x*x)/(sqrt(2*PI) * sigma[j]);
	       }
	  }
     }

   variable e = t.energy[kth];
   variable c = t.channel[kth];
   e = interpol (p, c, e);
   return (e,y);
   %return (p,y);
}

static define read_fef (file, ccdid, x, y)
{
   variable t = fits_read_table (file);
   variable i = where ((t.chipx_lo <= x)
		       and (t.chipx_hi >= x)
		       and (t.chipy_lo <= y)
		       and (t.chipy_hi >= y)
		       and (t.ccd_id == ccdid));
   
   variable p = [-4096:4096];
   variable e;
   (e,y) = compute_fef (t, p, i);
   %plot_open ("out.ps/CPS");

   pset_ylabel ("PHA");
   pset_xlabel ("ENERGY");
   pset_x (0.0, 10.0);
   pset_points;
   plot2d (t.energy[i], t.channel[i]);
   plot_pause ();

   pset_lines;
   variable j;
   for (j = 0; j < length(i); j++)
     {
	variable en = t.energy[i[j]];
	variable yy = y[j,*];
	variable de = shift (e,1)-e;
	variable mean_en = sum(yy * e *de)/sum(yy*de);
	variable mean_chan = sum (p*yy)/sum(yy);
	variable k = where (yy == max(yy))[0];
	pset_title (sprintf ("En=%g, Mean En=%g, Peak En=%g, Mean Chan=%g=%g keV", 
			     en, mean_en, e[k],
			     mean_chan, 
			     interpol (mean_chan, t.channel[i], t.energy[i])[0]));

	pset_x (e[k]-0.4, e[k]+0.4);
	%pset_logy; pset_y (0.001,);pset_x (0,3);
	plot2d (e,yy);
	plot_pause;
     }
   plot_close;
}

read_fef (File, Chip_Id, Chip_X, Chip_Y);
