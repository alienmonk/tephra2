# This script is designed to take the output file created by tephra2 and create a contour plot of the data.
# In order to run, the user must have the following programs installed
	# perl
	# generic mapping tools (GMT)
	# proj

#Additional files utilized by this program (in addition to the Input Files listed below):
	# perl script parse_tephra.out.pl available for download from the tephra2 tool
	# optional: replace 'gray' with a customized color pallet

# Example of input text file: Format: "latitude longitude txt_size text_angle font_number justification text"
#	-103.620 19.514 8 0 22 RB V Colima
#	-103.6091 19.5622 8 0 22 LB Nevado


# Input Files 
$T2out="Tephra2out.dat";	# Tephra2 output script
$vents = "vents.ll";		# File with vent locations, in lat long
$cities = "cities.ll";		# File containing city locations, in lat long
$contours="countours";
$topography1="16_08/ll.grd";	# .grd file, in this case derived from SRTM data

# Output Files
$tephra = "tephra2.out.utm.ll";		# Data in format for mapping with GMT
$out = "location_map.eps";		# Location Map with topography
$out2="location_map_no_topo.eps";	# Location map without topography

# Map Boundaries
$west = -104;			
$east = -102.5;
$south = 19;
$north = 20.5;
$cellsize = 0.000833333333*2;

# Map title
$Map_Label="Colima Plinian Eruption";

# Convert locations to lat long
system "perl parse_tephra.out.pl $T2out > tephra2.out.utm";
system "invproj +proj=utm +datum=WGS84 +ellps=WGS84 +zone=13 -f %.6f tephra2.out.utm > $tephra";

# Prepare extraneous files for GMT
system "makecpt -Cgray -V > topog.cpt";
system "grdgradient $topography1 -Gintensity1.grd -E45/80/.5/.2/.2/10 -Ne.1";

# Use GMT to create basemap
system "psbasemap --ANNOT_FONT_PRIMARY=2 --LABEL_FONT_SIZE=10 --ANNOT_FONT_SIZE=8 --D_FORMAT=%f -Y2i -X1i -Jm1:1500000 -R$west/$east/$south/$north -Ba0.5:'':/a0.5:'':/:.'$Map_Label':WSen -P -V -K > $out";
system "grdimage 16_08/ll.grd -Jm -R -Ctopog.cpt -Iintensity1.grd -V -O -K >> $out"; 
system "grdimage 16_09/ll.grd -Jm -R -Ctopog.cpt -Iintensity2.grd -V -O -K >> $out"; 

#Contour Map with topography
system "surface --D_FORMAT=%g $tephra -Gtephra.grd -I$cellsize -R -V";
system "pscoast --ANNOT_FONT_PRIMARY=1 --HEADER_FONT_SIZE=12 --LABEL_FONT_SIZE=10 --ANNOT_FONT_SIZE=8 -R -L-102.77/19.1/19.1/50 -T-104.6/17.5/.5i -Jm -Sblue -Di -Ia/0.25p/230 -Na/0.5p/50ta  -V -O -K >> $out";
system "grdcontour tephra.grd -C$contours -G0.5i+r0.1i -A100+f7+s4+kblack -Wa0.5p,black -Wc.18p,black -L0.1/10e6 -Jm -R$west/$east/$south/$north -K -O >> $out";
system "psxy $vents -Jm -R -Sc2p -G255/0/0 -V -O -K >> $out";
system "pstext $vents -Jm -R -S.5p,230 -V -O -K >> $out";
system "psxy $cities -Jm -R -Sc4p -G0/0/255 -W.1p,230/230/230 -V -O -K >> $out";
system "pstext $cities -Jm -S.4p,230 -R -V -O -K >> $out";

# Text to ensure whitespace boundary around map
system "pstext -R -J -G0 -O -P -N <<EOF>> $out
-103.350 18.85 10 0 14 BL kg/m^2
-104.17 18.82 1 0 14 BL .  
-102.35 18.82 1 0 14 BL . 
EOF";

# Convert eps file
system "ps2raster $out -A -Tg";

# Remove temporary files
system "rm *.eps";
system "rm tephra.grd";
system "rm topog.cpt";
system "rm intensity*";
system "rm tephra2.out*";

# Display Map
system "gm display location_map.png";
