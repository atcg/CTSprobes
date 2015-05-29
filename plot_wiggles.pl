#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $help = 0;
my $inFile;
my $outFile;
my $blastFile;

GetOptions  ("in=s"      => \$inFile,
             "blast=s"   => \$blastFile,
             "out=s"      => \$outFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outFile or $help) {
    die "Must supply --in and --out.\n";
}

open(my $inFH, "<", $inFile) or die "Couldn't open $inFile for reading: $!\n";


while (my $line = <$inFH>) {
    my $locus;
    if ($line =~ /(.*)\.depth/) {
        $locus = $1;
    }
    
    
    
    
    
    
    
    
}



