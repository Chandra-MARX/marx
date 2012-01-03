#!/usr/bin/env slsh
require ("fits");

define slsh_main ()
{
   if (__argc != 3)
     {
	() = fprintf (stderr, "Usage: %s <in-caldb-file> <out-marx-file>\n", __argv[0]);
	exit (1);
     }
   variable infile = __argv[1];
   variable outfile = __argv[2];

   variable fp = fits_open_file (outfile, "c");
   fits_create_image_hdu (fp, NULL, Char_Type, Int_Type[0]);
   fits_write_history (fp, sprintf ("This file was created using the marx script `%s %s %s`",
				    __argv[0], __argv[1], __argv[2]));

   variable keywords = struct
     {
	mission = "AXAF",
	instrum = "ACIS",
	detnam,
	ccd_id,
	marxvers = 5.0,
     };

   _for (0, 9, 1)
     {
	variable ccd = ();
	variable vfile = sprintf ("%s[AXAF_CONTAM%d]", infile, ccd+1);
	variable t = fits_read_table (vfile,
				      ["component", "n_energy", "energy", "mu",
				       "n_time", "time", "tau0", "tau1"]);
	keywords.detnam = "ACIS-$ccd"$;
	keywords.ccd_id = ccd;
	fits_write_binary_table (fp, "ACIS${ccd}_CONTAM"$, t, keywords);
     }

   fits_close_file (fp);
}
