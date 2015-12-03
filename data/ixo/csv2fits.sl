require ("readascii");
require ("setfuns");
require ("fits");

define slsh_main ()
{
   variable file = "mirror.csv";
   file = "/aluche/d1/xrtdesigns/axsio/axsio_prescription_10mFL_pbr_062211design.csv";

   variable header_start = "\"shell\",";
   variable header_start_len = strlen (header_start);
   variable comments = "";

   variable fp = fopen (file, "r");
   if (fp == NULL)
     {
	() = fprintf (stderr, "Unable to open %s\n", file);
	exit (1);
     }
   variable line;
   variable header_line = NULL;
   while (-1 != fgets (&line, fp))
     {
	line = strtrim_end (line, "\n");
	if (strncmp (line, header_start, header_start_len))
	  {
	     comments += line;
	     continue;
	  }
	header_line = line;
	break;
     }
   then
     {
	() = fprintf (stderr, "Unable to locate a line beginning with %s\n", header_start);
	exit (1);
     }

   % convert the comments to ascii
   comments = strtrans (comments, "^\\7", "?");

   variable col_names, col_data;

   header_line = strtrans (header_line, " ", "_");
   header_line = strtrans (header_line, "\"", "");
   col_names = strchop (header_line, ',', 0);
   % The column names are of the forms:
   %   name
   %   name_(unit)
   %   name1_name2_(unit)
   variable num_cols = length (col_names);
   variable units = String_Type[num_cols];
   variable i;
   _for i (0, num_cols-1, 1)
     {
	variable colname_unit = strchop (col_names[i], '(', 0);
	col_names[i] = strtrim (colname_unit[0], "_");
	if (length (colname_unit) == 1)
	  units[i] = "";
	else
	  units[i] = strtrim (colname_unit[1], "_)");
     }
   () = readascii (fp, &col_data; ncols=num_cols, delim=",");

   i = unique (col_names);
   if (length (i) != num_cols)
     {
	() = fprintf (stderr, "***Warning: Not all column names are unique.\n");
	col_data = col_data[i];
	col_names = col_names[i];
	units = units[i];
	num_cols = length (col_names);
     }

   variable tunit_names = array_map(String_Type, &sprintf, "tunit%d", [1:num_cols]);

   col_names = [col_names, "zmin", "zmax", "h_zmin", "h_zmax"];
   variable colstruct = @Struct_Type(col_names);
   variable tunit_struct = @Struct_Type(tunit_names);

   _for i (0, num_cols-1, 1)
     {
	set_struct_field (colstruct, col_names[i], col_data[i]);
	set_struct_field (tunit_struct, tunit_names[i], units[i]);
     }
   % Convert the shell column to an integer
   colstruct.shell = nint(colstruct.shell);
   variable nrows = length (colstruct.shell);

   variable flen = 10000, plen = 225, hlen = 225, gap=50;
   colstruct.zmin = Double_Type[nrows] + flen + 0.5*gap;
   colstruct.zmax = Double_Type[nrows] + flen + plen;
   colstruct.h_zmin = Double_Type[nrows] + flen - hlen;
   colstruct.h_zmax = Double_Type[nrows] + flen - 0.5*gap;

   variable outfile = "axsiomirror.fits";
   fits_write_binary_table (outfile, "IXO_MIRROR_GEOM",
			    colstruct, tunit_struct,
			    struct {comment = comments});
}
