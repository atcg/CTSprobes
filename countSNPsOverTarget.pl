#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;
use Data::Dumper;
use List::Util qw(sum);

my $help = 0;
my $blastFile;
my $vcf;

GetOptions  ("blast=s"        => \$blastFile,
             "vcf=s"          => \$vcf,
             "help|man"       => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$blastFile or !$vcf or $help) {
    die "Must supply --blast <blastoutput> --vcf <vcfFile>.\n The blast file should be a blast output file\
    in tabular format (-outfmt 6) where the target sequences have been blasted against\
    the chimera-masked RBBH sequences\n";
}

# Here's what the blast file should look like
#contig203242|E15C10|OPA allCTS_:_contig203242|E15C10|OPA_:_Contig001_contig203242|E15C10|OPA    98.20   389     7       0       1       389     499     887     0.0       680
#contig189275|E15E12|OPA allCTS_:_contig189275|E15E12|OPA_:_Contig001_contig189275|E15E12|OPA    96.72   366     9       1       1       366     516     878     9e-174    606
#contig42365|E15E3|OPA   allCTS_:_contig42365|E15E3|OPA_:_Contig001_contig42365|E15E3|OPA        95.48   398     9       2       4       393     555     951     7e-180    627


# First, we want to create a list of the regions where the targets overlap with the
# assembled contig that represents it as the reciprocal best blast hit. We'll put
# these into a hash:
my %targetMatchingRegions;
my %snpsWithinTargets;

my $searchIO = Bio::SearchIO->new(-file => $blastFile,
                              -format => 'blasttable');
while (my $result = $searchIO->next_result()) {
    my $hitCounter = 0;
    my $hit = $result->next_hit(); # Just get the best scoring hit (ignore the other ones)
    my $hsp = $hit->next_hsp(); # Just get the best scoring HSP (ignore the other ones)
    my $startingBase = $hsp->start('subject');
    my $endingBase = $hsp->end('subject');
    my $smallerBase;
    my $largerBase;
    if ($startingBase > $endingBase) {
        $largerBase = $startingBase;
        $smallerBase = $endingBase;
    } else {
        $smallerBase = $startingBase;
        $largerBase = $endingBase;
    }
    my @bases = $smallerBase .. $largerBase; # This generates a sequence of integers between these two endpoints (including the endpoints)
    foreach my $base (@bases) {
        $targetMatchingRegions{$hit->name()}{$base}++;
    }
}

# print Dumper \%targetMatchingRegions;

open(my $vcfFH, "<", $vcf) or die "Couldn't open $vcf for reading: $!\n";

while (my $line = <$vcfFH>) {
    my @fields = split(/\t/, $line);
    if (exists $targetMatchingRegions{$fields[0]}{$fields[1]} ) {
        #print "Match for $fields[0]\n";
        $snpsWithinTargets{$fields[0]}++;
    } else {
        #print "No match\n";
    }
}

foreach my $target (keys %snpsWithinTargets) {
    print $target . "\t" . $snpsWithinTargets{$target} / scalar(keys %{$targetMatchingRegions{$target}}) . "\n";
}

















