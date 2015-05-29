#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $help = 0;
my $inFile;
my $outFile;

GetOptions  ("in=s"      => \$inFile,
             "out=s"      => \$outFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outFile or $help) {
    die "Must supply --in and --out.\n";
}

my %removeHash;
$removeHash{"allCTS_:_contig315752|FAM198A|11_:_Contig001_contig315752|FAM198A|11"}++;
$removeHash{"allCTS_:_contig315214|SNX33|6_:_Contig002_contig315214|SNX33|6"}++;
$removeHash{"allCTS_:_contig359783|DERL2|6_:_Contig001_contig359783|DERL2|6"}++;
$removeHash{"allCTS_:_contig211570|TLCD1|11_:_Contig001_contig211570|TLCD1|11"}++;
$removeHash{"allCTS_:_contig124104|BET1|8_:_Contig001_contig124104|BET1|8"}++;
$removeHash{"allCTS_:_contig217149|RYBP|12_:_Contig001_contig217149|RYBP|12"}++;
$removeHash{"allCTS_:_contig334836|MOCS2|6_:_Contig003_contig334836|MOCS2|6"}++;
$removeHash{"allCTS_:_contig319555|NARS2|6_:_Contig001_contig319555|NARS2|6"}++;
            
my $seqIn = Bio::SeqIO->new(-file => $inFile,
                            -format => 'fasta');

my $seqOut = Bio::SeqIO->new(-file => ">$outFile",
                             -format => 'fasta');


while (my $seq = $seqIn->next_seq()) {
    next if exists($removeHash{$seq->display_id});
    $seqOut->write_seq($seq);
}
