require ("fits");

%variable File = "acisD2000-01-29fef_phaN0002.fits";
variable File = "acisD2000-01-29fef_phaN0005.fits";
variable Outfile = "test_newfef.fits";

variable Table = fits_read_table (File);
variable RegNum = Table.regnum;
variable Energy = Table.energy;

static define sort_fun (i, j)
{
   variable d = RegNum[i] - RegNum[j];
   if (d)
     return d;
   return sign (Energy[i] - Energy[j]);
}

static variable i = array_sort ([0:length(Energy)-1], &sort_fun);

foreach (get_struct_field_names (Table))
{
   static variable name = ();
   set_struct_field (Table, name, get_struct_field (Table, name)[i]);
}


fits_write_binary_table (Outfile, "FUNCTION", Table, NULL, 
			 sprintf ("This file differs from %s in that it is sorted differently",
				  File));

   
