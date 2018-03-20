# This script creates a directory of wind files for use by tephra2
# USAGE: combile_reanalysis_data.pl <configuration file>

use Math::Trig;

my $conf = "combine_reanalysis_data.conf";
print STDERR "Opening configuration file: $conf\n";
open CONF, "< $conf" or die "Can't open $conf : $!";

my %Param;
while(<CONF>) {
	if (/^$/ or /^#/) { next; }
	($key, $value) = split "=",$_;
        chomp($value);
		$Param{$key} = $value;
        print STDERR "$key=$Param{$key}\n";
} 
$level_file = $Param{LEVEL_FILE};
$uwind_file = $Param{UWIND_FILE};
$vwind_file = $Param{VWIND_FILE};

open (LEVEL, "<$level_file") || die ("$!");
open (UWIND, "<$uwind_file") || die ("$!");
open (VWIND, "<$vwind_file") || die ("$!");

@level = <LEVEL>;
@uwind = <UWIND>;
@vwind = <VWIND>;
@time_no;

my $data_ct = -1;
my $scale = 0;
my $offset = 0;
my $num_levels = 0;
foreach $lev (@level) {
	($name, $eq, $value) = split " ", $lev;
	if ($data_ct >=0) {
		chop($name);
		$level[$data_ct] = $name * $scale + $offset;
		$data_ct++;
     }
	if ($name eq "hgt") {
	  print STDERR "$name $eq\n";
	  $data_ct++;
	}
		if ($name eq "hgt:scale_factor") {
		#print STDERR "$name $eq $value\n";
		chomp($value);
		chop($value);
		print STDERR "$name $eq $value\n";
		$scale = $value;
    }
    if ($name eq "hgt:add_offset") {
      #print STDERR "$name $eq $value\n";
      chomp($value);
      chop($value);
      print STDERR "$name $eq $value\n";
      $offset = $value;
    }
    if ($name eq "level") {
      chomp($value);
      print STDERR "$name $eq $value\n";
      $num_levels = $value;
    }
}

$data_ct = -1;
$scale = 0;
$offset = 0;
foreach $uw (@uwind) {
	($name, $eq, $value) = split " ", $uw;
	if ($data_ct >=0) {
		chop($name);
		$uwind[$data_ct] = $name * $scale + $offset;
		$data_ct++;
     }
	if ($name eq "uwnd") {
	  print STDERR "$name $eq\n";
	  $data_ct++;
	}
	if ($name eq "uwnd:scale_factor") {
	    #print STDERR "$name $eq $value\n";
	    chomp($value);
		chop($value);
		print STDERR "$name $eq $value\n";
		$scale = $value;
    }
    if ($name eq "uwnd:add_offset") {
      #print STDERR "$name $eq $value\n";
      chomp($value);
      chop($value);
      print STDERR "$name $eq $value\n";
      $offset = $value;
    } 
}

$data_ct = -1;
$scale = 0;
$offset = 0;
foreach $vw (@vwind) {
	($name, $eq, $value) = split " ", $vw;
	if ($data_ct >=0) {
		chop($name);
		$vwind[$data_ct] = $name * $scale + $offset;
		$data_ct++;
     }
	if ($name eq "vwnd") {
	  print STDERR "$name $eq\n";
	  $data_ct++;
	}
		if ($name eq "vwnd:scale_factor") {
		 #print STDERR "$name $eq $value\n";
		  chomp($value);
		  chop($value);
		  print STDERR "$name $eq $value\n";
		  $scale = $value;
    }
    if ($name eq "vwnd:add_offset") {
     # print STDERR "$name $eq $value\n";
        chomp($value);
        chop($value);
        print STDERR "$name $eq $value\n";
        $offset = $value;
    } 
}
$data_ct--;
mkdir "wind_db", 0777 unless -d "wind_db";
my $file_no = 0;
my $j = $num_levels;
for ($i = 0; $i < $data_ct; $i++) {
	if ($j == $num_levels) {
		close OUT;
		$j = 0;
		my $file = "wind_db/wind$file_no";
		$file_no++;
		open OUT, "> $file" or die "Can't open $file: $!";
	}
	if ($uwind[$i] == 0) {
		if ($vwind[$i] < 0) {
			$wind_direction = 180;
		}
		elsif ($vwind > 0) {
			$wind_direction =0;
		}
	}
	else {
		$wind_direction = -180/pi *atan($vwind[$i]/$uwind[$i]);
			if ($uwind[$i] > 0) {
				$wind_direction += 90;
			}
			else {
				$wind_direction += 270;
			}
	}
	$wind_speed = sqrt($uwind[$i]*$uwind[$i]+$vwind[$i]*$vwind[$i]);
    printf OUT "%.0f %.1f %.1f\n", $level[$i], $wind_speed, $wind_direction;
    $j++;
}
$file_no--;
print STDERR "Wrote out $file_no wind files.\n";
print STDERR "Done\n";  
