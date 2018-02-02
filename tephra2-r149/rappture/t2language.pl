#!/usr/bin/env perl

# ----------------------------------------------------------------------
#  Wrapper Script for Tephra2: Language Selection
# ======================================================================
#  AUTHOR:  Leah M. Courtland
#  2010
# ======================================================================

use Rappture;

# open the XML file containing the run parameters
$driver = Rappture::RpLibrary->new($ARGV[0]);
$rapDir = $ENV{'TOOLDIR'} . "/rappture";

$language = $driver->get("input.choice(language).current");

if ($language eq "English"){print "Loading Tephra2: English User Interface. . . ";}
if ($language eq "Spanish"){print "Cargando Tephra2: Interfaz de Usuario en Español. . .";}

if ($language eq "English"){system "rappture -tool $rapDir/toolenglish.xml";}
if ($language eq "Spanish"){system "rappture -tool $rapDir/toolspanish.xml";}

$driver->result();
system "wmctrl -c Language_Select";


