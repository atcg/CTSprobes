#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $help = 0;
my $inFile;
my $outDir;
my $adapterType = "truseq";

GetOptions  ("in=s"       => \$inFile,
             "out=s"      => \$outDir,
             "adapters=s" => \$adapterType,
             "help|man"   => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outDir or $help) {
    die "Must supply --in and --out.\n";
}

unless (-d $outDir) {
    mkdir($outDir) or die "Couldn't make output directory $outDir: $!\n";
}

open(my $inFH, "<", $inFile) or die "Couldn't open $inFile for reading: $!\n";

while (my $line = <$inFH>) {
    if ($adapterType eq "truseq") { # This is the default
        if ($line =~ /(.*)\t(.*)/) {
            # $1 is the sample name
            # $2 is the index sequence
            my $indexedAdapterSeq = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" . $2 . "ATCTCGTATGCCGTCTTCTGCTTG";
            my $indexedAdapterRevComp = revcomp($indexedAdapterSeq);
            my $universalAdapterSeq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
            my $universalAdapterSeqRevComp = revcomp($universalAdapterSeq);
            open(my $outFH, ">", "$outDir/$1\.adapters") or die "Couldn't open $outDir/$1\.adapters for writing: $!\n";
            print $outFH ">$1\_indexed_adapter\n$indexedAdapterSeq\n";
            print $outFH ">$1\_indexed_adapter_revcomp\n$indexedAdapterRevComp\n";
            print $outFH ">$1\_universal_adapter\n$universalAdapterSeq\n";
            print $outFH ">$1\_universal_adapter_revcomp\n$universalAdapterSeqRevComp\n";
        }
    } elsif ($adapterType eq "itru") { # For our dual-indexed iTru5 and iTru7 smaples
        if ($line =~ /(.*)\t(.*)\t(.*)/) {
            open(my $outFH, ">", "$outDir/$1\.adapters") or die "Couldn't open $outDir/$1\.adapters for writing: $!\n";

            # $1 is the sample name
            # $2 is the iTru7 index
            # $3 is the iTru5 index
            ### TODO: write the iTru section
            
            my $i7index = revcomp($2); # These are actually revcomps of what I have in the current file (found by grepping)
            my $i5index = revcomp($3); # These are actually revcomps of what I have in the current file (found by grepping)

            my $i7seq = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" . $i7index . "ATCTCGTATGCCGTCTTCTGCTTG";
            my $i5seq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" . $i5index . "GTGTAGATCTCGGTGGTCGCCGTATCATT";
            my $i7revcomp = revcomp($i7seq);
            my $i5revcomp = revcomp($i5seq);
            
            print $outFH ">$1\_i7\n$i7seq\n";
            print $outFH ">$1\_i5\n$i5seq\n";
            print $outFH ">1_i7revcomp\n$i7revcomp\n";
            print $outFH ">1_i5revcomp\n$i5revcomp\n";
        }
    } else {
        die "I don't recognize the adapter choice you're trying to use. Please
        select --adapter truseq or --adapter itru\n";
    }
    
}


sub revcomp {
    my $dna = shift;

	# reverse the DNA sequence
    my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}




