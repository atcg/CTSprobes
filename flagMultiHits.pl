#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Bio::SearchIO;

my $inFile;
my $outFile;

GetOptions  ("in=s"           => \$inFile,
             "out=s"          => \$outFile) || die "Couldn't get options!";


my $searchIn = Bio::SearchIO->new(-file => $inFile,
                                  -format => "blasttable");

while (my $result = $searchIn->next_result()) {
    if ( scalar($result->hits) > 1) {
        print "Multiple hits found for " . $result->query_name() . "\n";
    }
    while (my $hit = $result->next_hit() ) {
        if (scalar($hit->hsps) > 1) {
            print "Multiple hsps found for " . $hit->name . "\n";
        }
    }
}
