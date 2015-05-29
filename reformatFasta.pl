#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $inFile;
my $outFile;

GetOptions  ("in=s"           => \$inFile,
             "out=s"          => \$outFile) || die "Couldn't get options!";


my $seqIn = Bio::SeqIO->new(-file => $inFile,
                            -format => 'fasta');

my $seqOut = Bio::SeqIO->new(-file => ">$outFile",
                             -format => 'fasta');


while (my $seq = $seqIn->next_seq()) {
    $seqOut->write_seq($seq);
}
