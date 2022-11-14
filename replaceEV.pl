#!/usr/bin/perl 
use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Data::Dumper;


my $help = "usage : perl replace.pl -t template -d iniDat -u tau -m rhomax
    \t This script replaces the key words TAU and RHO_MAX
    \t iniside a Therminator2 Blast Wave model .ini file
    \t parameters provided as command line
    \t template  == template iniDat file for HGEANT containng keywords,
    \t iniDat    == iniDat file will be created from template iniDat
    \t tau       == freeze-out time [fm]
    \t rhoMax   == maximum transverse radius [fm]
    \n";

my $template        = "";
my $iniDat          = "";
my $seed            = 10;

my $NARGS = $#ARGV ; #  $#ARGV will be rest by getopt()!!!

my %args;
getopt('tdum', \%args);

$template   = $args{"t"} ? $args{"t"} : "" ;
$iniDat     = $args{"d"} ? $args{"d"} : "" ;
$seed       = $args{"s"} ? $args{"s"} : 10 ;

my %parameters = ("template"    => $template,
                  "iniDat"      => $iniDat,
                  "seed"     => $seed,
                 );

#print Dumper \%parameters;

if($NARGS < 3 || 
   $template   eq "" ||
   $iniDat     eq ""
  ){
    print $help;
    exit(1);
}


 
open (rFile,"$template")
    or die "Cannot open initial ini file: $template!";

open (wFile,">","$iniDat")
    or die "Cannot open modified ini file: $iniDat!";

while(<rFile>)
{
    s/SEED/$seed/;
    print wFile "$_";
}

close rFile;
close wFile;
