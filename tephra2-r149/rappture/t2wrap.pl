#!/usr/bin/env perl

use Rappture;

# open the XML file containing the run parameters
$driver = Rappture::RpLibrary->new($ARGV[0]);

# Set variable paths
$exampleDir = $ENV{'TOOLDIR'} . "/examples";
$rapDir = $ENV{'TOOLDIR'} . "/rappture";
$saveDir = $ENV{'HOME'} . "/tephra2output";

#set filenames
$run_number = $driver->get("input.choice(Run_Number).current");
$save_workspace = $driver->get("input.boolean(Save_Workspace).current");
$configname="T2_Input.txt";
$outputname="T2_Output.txt";
$windname="T2_Wind.txt";
$gridname="grid.txt";
$locmapname="location_map.";
if ($save_workspace eq "yes"){
$configname="T2_Input".$run_number.".txt";
$outputname="T2_Output".$run_number.".txt";
$windname="T2_Wind".$run_number.".txt"; 
$gridname="grid".$run_number.".txt";
$locmapname = "location_map".$run_number.".";
}

$ReWind = $driver->get("input.boolean(ReWind).current");


#Find out what language to print in
$run_number = $driver->get("input.string(Run_Number).current");
$language = $driver->get("input.string(language).current");

# Obtain run parameters
$Plume_Height = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Plume_Height).current");
$Eruption_Mass = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Eruption_Mass).current");
$Max_Grainsize = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Max_Grainsize).current");
$Min_Grainsize = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Min_Grainsize).current");
$Median_Grainsize = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Median_Grainsize).current");
$STD_Grainsize = $driver->get("input.group(tabs).group(Eruption_Parameters).number(STD_Grainsize).current");
$Vent_Easting = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Vent_Easting).current");
$Vent_Northing = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Vent_Northing).current");
$Vent_Elevation = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Vent_Elevation).current");
$Eddy_Const = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Eddy_Constant).current");
$Diffusion_Coefficient = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Diffusion_Coefficient).current");
$Fall_Time_Threshold = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Fall_Time_Threshold).current");
$Lithic_Density = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Lithic_Density).current");
$Pumice_Density = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Pumice_Density).current");
$Col_Steps = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Column_Steps).current");
$Plume_Model = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Plume_Model).current");
$Plume_Ratio = $driver->get("input.group(tabs).group(Eruption_Parameters).number(Plume_Ratio).current");
#$Save_File = $driver->get("input.(Save_File).current");
#$Filename = $driver->get("input.(Filename).current");

$wind = $driver->get("input.(wind).current");
$grid = $driver->get("input.(grid).current");
$config = $driver->get("input.(config).current");
$configexist = length ($config);

# Making the grid:
$Make_grid = $driver->get("input.group(tabs).group(Grid).boolean(Grid_loaded).current");


if ($Make_grid eq 'no') {
$Elev_Landing = $driver->get("input.group(tabs).group(Grid).number(Elev_Landing).current");
$Min_East = $driver->get("input.group(tabs).group(Grid).number(Min_East).current")*1000;
$Max_East = $driver->get("input.group(tabs).group(Grid).number(Max_East).current")*1000;
$Min_North = $driver->get("input.group(tabs).group(Grid).number(Min_North).current")*1000;
$Max_North = $driver->get("input.group(tabs).group(Grid).number(Max_North).current")*1000;
$Spacing = $driver->get("input.group(tabs).group(Grid).number(Spacing).current")*1000;

$MinE = $Vent_Easting - $Min_East;
$MaxE = $Vent_Easting + $Max_East;
$MinN = $Vent_Northing - $Min_North;
$MaxN = $Vent_Northing + $Max_North;

open FILEGRID, ">$gridname";
for ($x = $MinE; $x <= $MaxE; $x += $Spacing) {
  for ($y = $MinN; $y <= $MaxN; $y += $Spacing) {
	print FILEGRID "$x $y $Elev_Landing\n";
  }
}  
close(FILEGRID);}

if ($Make_grid eq 'yes') {
open(GRID, ">$gridname");
print GRID "$grid\n";
close(GRID);}

# Obtain Plot Location Parameters
$Volcano_Name = $driver->get("input.group(tabs).group(Locations).string(Volcano_Name).current");
$Loc_Name1 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).string(Location1).current");
$Easting1 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Easting1).current");
$Northing1 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Northing1).current");
$Loc_Name2 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).string(Location2).current");
$Easting2 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Easting2).current");
$Northing2 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Northing2).current");
$Loc_Name3 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).string(Location3).current");
$Easting3 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Easting3).current");
$Northing3 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Northing3).current");
$Loc_Name4 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).string(Location4).current");
$Easting4 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Easting4).current");
$Northing4 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Northing4).current");
$Loc_Name5 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).string(Location5).current");
$Easting5 = $driver->get("input..group(tabs).group(Locations).group(Locations2Plot).number(Loc_Easting5).current");
$Northing5 = $driver->get("input.group(tabs).group(Locations).group(Locations2Plot).number(Loc_Northing5).current");

# Only plotting?
$RePlot = $driver->get("input.group(tabs).group(Plot_Refinement).boolean(RePlot).current");

#Obtain map boundary changes
$Extra_Left = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Map_Bound).number(Extra_Left).current")*1000;
$Extra_Right = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Map_Bound).number(Extra_Right).current")*1000;
$Extra_Top = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Map_Bound).number(Extra_Top).current")*1000;
$Extra_Bottom = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Map_Bound).number(Extra_Bottom).current")*1000;

#Obtain label locations
$EshiftV = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Volc_Name).current");
$EshiftL1 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc1).current");
$EshiftL2 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc2).current");
$EshiftL3 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc3).current");
$EshiftL4 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc4).current");
$EshiftL5 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc5).current");

$NshiftV = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Volc_NameN).current");
$NshiftL1 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc1N).current");
$NshiftL2 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc2N).current");
$NshiftL3 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc3N).current");
$NshiftL4 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc4N).current");
$NshiftL5 = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Plot_Labels).number(Shift_Loc5N).current");

#Reset Plot Locations that equal zero
if ($Easting1 == 0){$Easting1=$Vent_Easting};
if ($Northing1 == 0){$Northing1=$Vent_Northing}; 
if ($Easting2 == 0){$Easting2=$Vent_Easting}; 
if ($Northing2 == 0){$Northing2=$Vent_Northing};
if ($Easting3 == 0){$Easting3=$Vent_Easting}; 
if ($Northing3 == 0){$Northing3=$Vent_Northing};
if ($Easting4 == 0){$Easting4=$Vent_Easting}; 
if ($Northing4 == 0){$Northing4=$Vent_Northing} ;
if ($Easting5 == 0){$Easting5=$Vent_Easting;} 
if ($Northing5 == 0){$Northing5=$Vent_Northing};


# Obtain Contours
$Contours = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Contours).string(Contours2Plot).current");
$Cmin = $driver->get("input.group(tabs).group(Plot_Refinement).group(tabs).group(Contours).number(Cmin).current");


# Set the wind field
if ($wind eq "constant") {
    $wind = "";
    open(IN, "$exampleDir/w_constant.wind");
    while (<IN>) { $wind .= $_; }
    close(IN);}

if ($wind eq "stratified") {
    $wind = "";
    open(IN, "$exampleDir/w_stratified.wind");
    while (<IN>) { $wind .= $_; }
    close(IN);}

if (($wind eq "Random_CN") && ($ReWind eq "no")) {
$range = 1826;
$wind_directory = "$exampleDir/wind_CN1990_1995";
$random_number = int(rand($range));
$windfilename="$wind_directory/wind$random_number";
local( $/, *EXWIND ) ;
open EXWIND, "$windfilename";
$wind=<EXWIND>;
close(EXWIND);}

if (($wind eq "Random_I") && ($ReWind eq "no")) {
$range = 1826;
$wind_directory = "$exampleDir/wind_I1965_1970";
$random_number = int(rand($range));
$windfilename="$wind_directory/wind$random_number";
local( $/, *EXWIND ) ;
open EXWIND, "$windfilename";
$wind=<EXWIND>;
close(EXWIND);}

if (($wind eq "Random_P") && ($ReWind eq "no")) {
$range = 1826;
$wind_directory = "$exampleDir/wind_P2005_2010";
$random_number = int(rand($range));
$windfilename="$wind_directory/wind$random_number";
local( $/, *EXWIND ) ;
open EXWIND, "$windfilename";
$wind=<EXWIND>;
close(EXWIND);}



# Save the data to the respective files so tephra can access them
open(WIND, ">$windname");
print WIND "$wind\n";
close(WIND);

if ($ReWind eq "no"){system "cp $windname wind.tmp";}
if ($ReWind eq "yes"){system "cp wind.tmp $windname";}




# Create configuration file
$configtoprint="PLUME_HEIGHT $Plume_Height\n
ERUPTION_MASS $Eruption_Mass\n
MAX_GRAINSIZE $Max_Grainsize\n
MIN_GRAINSIZE $Min_Grainsize\n
MEDIAN_GRAINSIZE $Median_Grainsize\n
STD_GRAINSIZE $STD_Grainsize\n
VENT_EASTING $Vent_Easting\n
VENT_NORTHING $Vent_Northing\n
VENT_ELEVATION $Vent_Elevation\n
EDDY_CONST $Eddy_Const\n
DIFFUSION_COEFFICIENT $Diffusion_Coefficient\n
FALL_TIME_THRESHOLD $Fall_Time_Threshold\n
LITHIC_DENSITY $Lithic_Density\n
PUMICE_DENSITY $Pumice_Density\n
COL_STEPS $Col_Steps\n
PLUME_MODEL $Plume_Model\n
PLUME_RATIO $Plume_Ratio";


if ($configexist > 100){
open CONFIG, ">$configname";
print CONFIG "$config\n";
close(CONFIG);}
else {
open FILE1, ">$configname";
print FILE1 "$configtoprint";
close(FILE1);}


# Create Contour File
open CONT, ">contours";
print CONT "$Contours\n";
close(CONT);



#Run tephra2: output script: Tephra2out#.dat
if (($language eq "English") && ($RePlot eq "no")) {print"\nTEPHRA2 is calculating the mass accumulation of tephra at the points specified...\n\n";}
if (($language eq "Spanish") && ($RePlot eq "no")) {print"\nTEPHRA2 esta calculando la masa de tefra acumulada en los puntos especificados...\n\n";}

if ($RePlot eq "no"){
	system "tephra2 $configname $gridname $windname > Tephra2out.dat";}


###################################  Creating isopach map #############################
$tephra = "tephra2.out.temp";
if ($language eq "English"){print"\nFormatting Tephra2 output for plotting\n";}
if ($language eq "Spanish"){print"\nFormateando los resultados de Tephra2 para graficarlos\n";}

if ($RePlot eq "no"){system "perl $rapDir/parse_tephra.out.pl Tephra2out.dat > $tephra";}
system "mv Tephra2out.dat $outputname";
system "cp $rapDir/.gmt* .";



# Determine the max and min dimensions of the input grid and
# The max and min dimensions of the map where >0.1kg/m^2 of tephra accumulated
open(INP,"$tephra");
$ii=-1;
$i=0;


while (<INP>) {
 ($east, $north, $mass) = split " ", $_;
if($ii>(-1)){
	if ($mass >= $Cmin)
		{ $xx[$ii]=$east;
		  $yy[$ii]=$north;
		  $masa[$ii]=$mass;

  		  $xplot[$i]=$east;
		  $yplot[$i]=$north;
		  $i++;}
	else
		{ $xx[$ii]=$east;
		  $yy[$ii]=$north;
		  $masa[$ii]=$mass;}}
	  $ii++;
}
close(INP);


#Consider location points in determining map size
$xshift1=8000;	#Make room to write volcano name
$xshift2=5000;	#Make room for symbol to appear on map
$yshift=4000;	#Make room for symbol to appear on map
push(@xplot,$Vent_Easting,$Vent_Easting+$xshift1,$Vent_Easting-$xshift2, $Easting1, $Easting1+$xshift1,$Easting1-$xshift2,$Easting2, $Easting2+$xshift1,$Easting2-$xshift2, $Easting3,$Easting3+$xshift1,$Easting3-$xshift2, $Easting4, $Easting4+$xshift1,$Easting4-$xshift2,$Easting5, $Easting5+$xshift1,$Easting5-$xshift2);
push(@yplot,$Vent_Northing,$Vent_Northing+$yshift,$Vent_Northing-$yshift, $Northing1, $Northing1+$yshift,$Northing1-$yshift, $Northing2, $Northing2+$yshift,$Northing2-$yshift, $Northing3, $Northing3+$yshift, $Northing3-$yshift, $Northing4, $Northing4+$yshift, $Northing4-$yshift, $Northing5, $Northing5+$yshift, $Northing5-$yshift);


#Sort points in order to set map boundaries
@xsortm=sort {$a <=> $b}  (@xplot);
@ysortm=sort {$a <=> $b} (@yplot);
@xsort=sort {$a <=> $b} (@xx);
@ysort=sort {$a <=> $b} (@yy);

$xminm = $xsortm[0]-$Extra_Left;
$xmaxm = $xsortm[-1]+$Extra_Right;
$yminm = $ysortm[0]-$Extra_Bottom;
$ymaxm = $ysortm[-1]+$Extra_Top;
$xmin = $xsort[0];
$xmax = $xsort[-1];
$ymin = $ysort[0];
$ymax = $ysort[-1];
$delta_x=($xmax-$xmin)/500;


#Shift Labels for better viewing
$Vent_North_Name=$Vent_Northing+$NshiftV;
$North1=$Northing1+$NshiftL1;
$North2=$Northing2+$NshiftL2;
$North3=$Northing3+$NshiftL3;
$North4=$Northing4+$NshiftL4;
$North5=$Northing5+$NshiftL5;

$Vent_East_Name=$Vent_Easting+$EshiftV;
$East1=$Easting1+$EshiftL1;
$East2=$Easting2+$EshiftL2;
$East3=$Easting3+$EshiftL3;
$East4=$Easting4+$EshiftL4;
$East5=$Easting5+$EshiftL5;
$Volcano="Volcano";

# Use GMT to make the plot
$out="$locmapname"."eps";
$locmapname2="$locmapname"."png";
if ($language eq "English"){print"\nCreating the location map and contouring the data. . . \n";}
system "psbasemap --FONT_ANNOT_PRIMARY=12 --FONT_LABEL=10 --FONT_ANNOT_SECONDARY=8 --FORMAT_FLOAT_OUT=%.0f -JX4i -R$xminm/$xmaxm/$yminm/$ymaxm -Ba10000:'':/a10000:'':/:.'Isomass Map': -P -V -K > $out";
system "surface --FORMAT_FLOAT_OUT=%g $tephra -Gtephra.grd -I$delta_x -R$xmin/$xmax/$ymin/$ymax -V";
system "psxy -R$xminm/$xmaxm/$yminm/$ymaxm -JX4i -G0/0/0 -O -P -K -St0.25i -W0.5 -V << EOF>> $out
$Vent_Easting $Vent_Northing
EOF";
system "psxy -R$xminm/$xmaxm/$yminm/$ymaxm -JX4i -G0/0/0 -O -P -K -Sc0.1i -W0.5 -V << EOF>> $out
$Easting1 $Northing1 
$Easting2 $Northing2 
$Easting3 $Northing3 
$Easting4 $Northing4 
$Easting5 $Northing5 
EOF";
system "pstext -R$xminm/$xmaxm/$yminm/$ymaxm -JX4i -F+f12p+jBL -O -P -K << EOF>> $out
$Vent_East_Name $Vent_North_Name $Volcano_Name 
EOF";
system "pstext -R$xminm/$xmaxm/$yminm/$ymaxm -JX4i -F+f11p+jBL -O -P -K << EOF>> $out
$East1 $North1 $Loc_Name1
$East2 $North2 $Loc_Name2
$East3 $North3 $Loc_Name3
$East4 $North4 $Loc_Name4
$East5 $North5 $Loc_Name5
EOF";
system "grdcontour tephra.grd -Ccontours -Gd0.5i+r0.1i -A+f7+s8+kblack -Wa0.5p,black -Wc.18p,black -L0.1/10e6 -JX4i -O -R$xminm/$xmaxm/$yminm/$ymaxm >> $out";


###################################### Done Creating Isopach Map  ####################

#Remove unnecesary files
system "ps2raster $out -A -Tg";
system "rm *.eps";
system "rm grid.txt";
system "rm node_";
system "rm tephra.grd";
system "rm .gmt*";
system "rm contours";
#system "rm tephra2.out.utm";


#Load Output into graphical user interface (rappture)
local( $/, *OUTPUT ) ;
open OUTPUT, "$outputname";
$test=<OUTPUT>;

$driver->put("output.image(outa).about.label","Isopach Map (kg/m^2):",0);
$driver->putFile("output.image(outa).current", "$locmapname2", 1, 0);

$driver->put("output.string(outb).about.label","TEPHRA2 Output",0);
$driver->put("output.string(outb).current",$test,0);

$driver->put("output.string(outc).about.label","TEPHRA2 Model Inputs",0);
$driver->put("output.string(outc).current",$configtoprint,0);

$driver->put("output.string(outd).about.label","TEPHRA2 Wind Input",0);
$driver->put("output.string(outd).current",$wind,0);

$driver->result();


if ($save_workspace eq "yes"){

#$savedir=

if (-d $saveDir){}
else {system "mkdir $saveDir";}

system "mv $configname $saveDir/$configname"; 
system "mv $outputname $saveDir/$outputname";
system "mv $windname $saveDir/$windname";  
system "mv $gridname $saveDir/$gridname"; 
system "mv $locmapname2 $saveDir/$locmapname2"; 
}

system "rm T2_Output.txt"; 
system "rm T2_Input.txt"; 
system "rm $configname"; 
system "rm $locmapname2";
system "rm tephra2.out.utm";
system "rm $windname";

exit;
