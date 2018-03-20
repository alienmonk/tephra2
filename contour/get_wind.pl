#!/usr/bin/env perl

use URI;
use URI::Escape ('uri_escape');
use LWP::UserAgent;
use Net::FTP;
#use Rappture;


# #####################  *** Enter Parameters *** #####################

$Easting = 110.444817;		# Longitude, in degrees (NOT UTM Easting!)
$Longitude = "E"; 	# E or W
$Northing = 7.540831;
$Latitude = "S";	# N or S

$Year_Begin = 2006;	# 1948 to present	
$Month_Begin = "April";	# Jan, Feb, March, April, May, June, July, Aug, Sept, etc
$Day_Begin = 1;
$Year_End = 2006;
$Month_End = "Aug";
$Day_End = 6;

# Reformat Locations
$Long="${Easting}${Longitude}";
$Lat="${Northing}${Latitude}";

#########################################################################

print "Longitude: $Long\n";
print "Latitude: $Lat\n";
print "Start Date: Year:$Year_Begin Month:$Month_Begin Day:$Day_Begin\n";
print "End Date:   Year:$Year_End Month:$Month_End Day:$Day_End\n";



#--------------Searching for the files ------------------
# -------------------------------------------------------

# Obtaining daily Numbers for filename
#---------------
$urla = URI ->new('http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&Variable=U-wind&group=0&submit=Search');

my $responsea = LWP::UserAgent->new->get( $urla );
die "Error: ", $responsea->status_line unless $responsea->is_success;
#     print $response->content;

# Extracting the numbers
extract3($responsea->content);
	$N3=$N3temp;
    print "Number 3: $N3\n";

#------------
$urlb = URI ->new('http://www.esrl.noaa.gov/psd/cgi-bin/DataAccess.pl?DB_dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&DB_variable=u-wind&DB_statistic=Mean&DB_tid=' .
uri_escape($N3) .
'&DB_did=32&DB_vid=666');

print $urlb, "\n"; # so we can see it

my $responseb = LWP::UserAgent->new->get( $urlb );
die "Error: ", $responseb->status_line unless $responseb->is_success;
#     print $response->content;

# Extracting the numbers
extract2($responseb->content);
	$N1=$N1temp;
	$N2=$N2temp;
    print "Number 1: $N1\nNumber2: $N2\n";
#--------------

# ---------------- Geopotential Height, Daily-------------
# The URL:
$url_P = URI ->new(
'http://www.esrl.noaa.gov/psd/cgi-bin/GrADS.pl?dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&DB_did=32&file=%2FDatasets%2Fncep.reanalysis.dailyavgs%2Fpressure%2Fhgt.1948.nc+hgt.%25y4.nc+' .
uri_escape($N1) .
'&variable=hgt&DB_vid=663&DB_tid=' .
uri_escape($N2) .
'&units=m&longstat=Mean&DB_statistic=Mean&stat=&lat-begin=' .
uri_escape($Lat) .
'&lat-end=' .
uri_escape($Lat) . 
'&lon-begin='.
uri_escape($Long) . 
'&lon-end=' .
uri_escape($Long) .
'&dim0=level&level+units=millibar&level=1000.00&level=925.00&level=850.00&level=700.00&level=600.00&level=500.00&level=400.00&level=300.00&level=250.00&level=200.00&level=150.00&level=100.00&level=70.00&level=50.00&level=30.00&level=20.00&level=10.00&dim1=time&year_begin=' .
uri_escape($Year_Begin) .
'&mon_begin=' . 
uri_escape($Month_Begin) .
'&day_begin=' .
uri_escape($Day_Begin) .
'&year_end=' .
uri_escape($Year_End) .
'&mon_end=' .
uri_escape($Month_End) .
'&day_end=' .
uri_escape($Day_Begin) .
'&X=lon&Y=lat&output=file&bckgrnd=black&cint=&range1=&range2=&scale=100&submit=Create+Plot+or+Subset+of+Data');

#print $url_P, "\n"; # so we can see it



# Acquiring the web page
my $responseP = LWP::UserAgent->new->get( $url_P );
die "Error: ", $responseP->status_line unless $responseP->is_success;
#     print $response->content;

# Extracting the file name
extract($responseP->content);
	$filenameP=$tempfilename;
        print "Pressure Filename: $filenameP\n";





# ----------------- U-wind, Daily ------------------------

print "Looking up U-Wind data . . .\n";

# The URL:
$url_U = URI ->new(
'http://www.esrl.noaa.gov/psd/cgi-bin/GrADS.pl?dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&DB_did=32&file=%2FDatasets%2Fncep.reanalysis.dailyavgs%2Fpressure%2Fuwnd.1948.nc+uwnd.%25y4.nc+' .
uri_escape($N1) . 
'&variable=uwnd&DB_vid=666&DB_tid=' .
uri_escape($N2) . 
'&units=m%2Fs&longstat=Mean&DB_statistic=Mean&stat=&lat-begin=' .
uri_escape($Lat) .
'&lat-end=' .
uri_escape($Lat) . 
'&lon-begin='.
uri_escape($Long) . 
'&lon-end=' .
uri_escape($Long) .
'&dim0=level&level+units=millibar&level=1000.00&level=925.00&level=850.00&level=700.00&level=600.00&level=500.00&level=400.00&level=300.00&level=250.00&level=200.00&level=150.00&level=100.00&level=70.00&level=50.00&level=30.00&level=20.00&level=10.00&dim1=time&year_begin=' .
uri_escape($Year_Begin) .
'&mon_begin=' . 
uri_escape($Month_Begin) .
'&day_begin=' .
uri_escape($Day_Begin) .
'&year_end=' .
uri_escape($Year_End) .
'&mon_end=' .
uri_escape($Month_End) .
'&day_end=' .
uri_escape($Day_End) .
'&X=lon&Y=lat&output=file&bckgrnd=black&use_color=on&cint=&range1=&range2=&scale=100&submit=Create+Plot+or+Subset+of+Data');
 
#print $url_U, "\n"; # so we can see it

# Acquiring the web page
my $responseU = LWP::UserAgent->new->get( $url_U );
die "Error: ", $responseU->status_line unless $responseU->is_success;
#     print $response->content;

# Extracting the file name
extract($responseU->content);
	$filenameU=$tempfilename;
    print "U-Wind Filename: $filenameU\n";

# ---------------- V-Wind, Daily -------------

print"Looking up V-Wind data . . .\n";

# The URL:
$url_V = URI ->new('http://www.esrl.noaa.gov/psd/cgi-bin/GrADS.pl?dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&DB_did=32&file=%2FDatasets%2Fncep.reanalysis.dailyavgs%2Fpressure%2Fvwnd.1948.nc+vwnd.%25y4.nc+' .
uri_escape($N1) .
'&variable=vwnd&DB_vid=667&DB_tid=' .
uri_escape($N2) .
'&units=m%2Fs&longstat=Mean&DB_statistic=Mean&stat=&lat-begin=' .
uri_escape($Lat) . 
'&lat-end=' .
uri_escape($Lat) . 
'&lon-begin=' .
uri_escape($Long) . 
'&lon-end=' .
uri_escape($Long) . 
'&dim0=level&level+units=millibar&level=1000.00&level=925.00&level=850.00&level=700.00&level=600.00&level=500.00&level=400.00&level=300.00&level=250.00&level=200.00&level=150.00&level=100.00&level=70.00&level=50.00&level=30.00&level=20.00&level=10.00&dim1=time&year_begin=' .
uri_escape($Year_Begin) .
'&mon_begin=' . 
uri_escape($Month_Begin) .
'&day_begin=' .
uri_escape($Day_Begin) .
'&year_end=' .
uri_escape($Year_End) .
'&mon_end=' .
uri_escape($Month_End) .
'&day_end=' .
uri_escape($Day_End) .
'&X=lon&Y=lat&output=file&bckgrnd=black&use_color=on&cint=&range1=&range2=&scale=100&submit=Create+Plot+or+Subset+of+Data');
  
# Acquiring the web page
my $responseV = LWP::UserAgent->new->get( $url_V );
die "Error: ", $responseV->status_line unless $responseV->is_success;
#     print $response->content;

# Extracting the file name
extract($responseV->content);
	$filenameV=$tempfilename;
    print "V-Wind Filename: $filenameV\n";
# ----------------------------------------------------------------


# ----------------- FTP the files ----------------------------

print "\nEstablishing connection with ftp.cdc.noaa.gov . . .\n";

# Specify ftp site, directory, and filename
my $host="ftp.cdc.noaa.gov";
my $directory="Public/www";

# Connect to CDC site
$ftp=Net::FTP->new($host,Timeout=>240) or $newerr=1;
  push @ERRORS, "Can't ftp to $host: $!\n" if $newerr;
  myerr() if $newerr;
print "Connected\n";

# Login to the ftp server as anonymous; 
$ftp->login('anonymous') 
   or die "Can't login ($CPANhost):" . $ftp->message;  

#Switch to the appropriate directory
$ftp->cwd($directory) or $newerr=1; 
  push @ERRORS, "Can't cd  $!\n" if $newerr;
  myerr() if $newerr;
  $ftp->quit if $newerr;

# Enable binary mode
$ftp->binary();

# Get the files
print "Downloading Pressure data . . .\n";
$ftp->get($filenameP)	 or 
            warn "Couldn't get '$filenameP', skipped: $!";  

print "Downloading U-Wind data . . .\n";
$ftp->get($filenameU)	 or 
            warn "Couldn't get '$filenameU', skipped: $!"; 

print "Downloading V-Wind data . . .\n";
$ftp->get($filenameV)	 or 
            warn "Couldn't get '$filenameV', skipped: $!";
 

# Close the connection to the FTP server. 
$ftp->quit or die "Couldn't close the connection 
                      cleanly: $!";  


# ----------------------------------------------------------

# ---------------------Format data for use with tephra2 ------------

print "\nFormatting the data for use with Tephra2. . .\n";

# Use netcdf to recover the data
system "ncdump -v hgt,time -f C $filenameP > geo_height.nc.out";
system "ncdump -v uwnd,time -f C $filenameU > u-wind.nc.out";
system "ncdump -v vwnd,time -f C $filenameV > v-wind.nc.out";

system "perl combine_reanalysis_data.pl";

system "rm X*";



print "\nProgram Complete\n\n";


# ------------------------------------------------------------
# -------------------------- Subroutines ---------------------




sub extract { 
    $_[0] =~ m{<ul>(.*)</ul>}s;
    foreach my $entry (split /<li>/, $1) {
	next unless $entry =~ m{a href=ftp(.*?)}s;
	$length_back=index($entry,">FTP");
	$length_front=42;
	$length=$length_back-$length_front;
#	print "Length: $length\n";
	$tempfilename = substr($entry,42,$length);

  }
}

sub extract2 { 
    $_[0] =~ m{<em>(.*)</em>}s;
#	print "$_[0]\n";
   foreach my $entry (split /\//, $1) {
	next unless $entry =~ m{uwnd.1948.nc uwnd.%y4.nc (.*?)}s;
	$N1temp = substr($entry,26,4);
	$N2temp = substr($entry,-25,5);
  }
}

sub extract3 { 
    $_[0] =~ m{<dl>(.*)</dl>}s;
#	print "$_[0]\n";
	my $i = 1;
   foreach my $entry (split /"/, $1) {
	$i++;
	$entry2[$i]=$entry;
  }		
#print "Entry: $entry2[-2]\n";
$N3temp = substr($entry2[-2],73,5);
}




exit;
