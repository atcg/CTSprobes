#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;
use Data::Dumper;

my $help = 0;
my $inFile;
my $sequences;
my $outFile;

GetOptions  ("in=s"      => \$inFile,
             "sequences=s" => \$sequences,
             "out=s"      => \$outFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$sequences or !$outFile or $help) {
    die "Must supply --in and --sequences and --out.\n";
}


# All_:_contig140241|CRMP1|6_:_Unfinished001_contig140241|CRMP1|6 All_:_contig140241|CRMP1|6_:_Unfinished001_contig140241|CRMP1|6 100.00  1233    0       0       1       1233    1       1233    0.0     2278
# All_:_contig140265|AQP1|11_:_Unfinished001_contig140265|AQP1|11 All_:_contig140265|AQP1|11_:_Unfinished001_contig140265|AQP1|11 100.00  1370    0       0       1       1370    1       1370    0.0     2531
# All_:_contig140271|TCF4|8_:_Unfinished001_contig140271|TCF4|8   All_:_contig140271|TCF4|8_:_Unfinished001_contig140271|TCF4|8   100.00  1304    0       0       1       1304    1       1304    0.0     2409
# All_:_contig140316|NEMF|11_:_Unfinished001_contig140316|NEMF|11 All_:_contig140316|NEMF|11_:_Unfinished001_contig140316|NEMF|11 100.00  1451    0       0       1       1451    1       1451    0.0     2680
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   100.00  2045    0       0       1       2045    1       2045    0.0     3777
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig162323|MLLT4|8_:_Unfinished001_contig162323|MLLT4|8 84.21   190     7       9       1       167     111     300     1e-39    163
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig71398|6-Mar|1_:_Unfinished005_contig71398|6-Mar|1   83.78   185     7       5       1       162     33      217     7e-37    154
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig214506|CAMSAP2|6_:_Unfinished001_contig214506|CAMSAP2|6     83.78   185     7       5       1       162     1332    1148    7e-37   
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig155894|RANBP3L|6_:_Unfinished001_contig155894|RANBP3L|6     83.06   183     8       5       1       160     244     426     4e-34   

my $searchIO = Bio::SearchIO->new(-file=>$inFile,
                                  -format=>'blasttable');
my $seqIn = Bio::SeqIO->new(-file=>$sequences,
                            -format=>'fasta');

my %seqHash;
while (my $seq = $seqIn->next_seq()) {
    
    $seqHash{$seq->display_id()} = $seq->seq();
}
my %maskingHash;
while (my $result = $searchIO->next_result()) {
    my $hitCounter = 0;
    while (my $hit = $result->next_hit()) {
        # Look for multiple hits. First one should be 100% and be to itself
        # For all other hits, we want to mask the regions that are blasting to those regions
        if ($hitCounter == 0) {
            $hitCounter++;
            next;
        }
        while (my $hsp = $hit->next_hsp()) {
            # Translate the query and subject names into the actual target names
            my $subjectName = $hit->name();
            my $queryName = $result->query_name();

            # Now find the offending sections:
            my $queryStart = $hsp->start('query');
            my $queryEnd = $hsp->end('query');
            my $queryString = $queryStart . ":" . $queryEnd;
            my $subjectStart = $hsp->start('subject');
            my $subjectEnd = $hsp->end('subject');
            my $subjectString = $subjectStart . ":" . $subjectEnd;
            # And store them:
            push(@{$maskingHash{$queryName}}, $queryString);
            push(@{$maskingHash{$subjectName}}, $subjectString);
        }
        # For all of these remaining hits, we want to eventually mask the overlapping sequence in both the query and the subject contigs
        $hitCounter++;   
    }
}


foreach my $target (sort keys %maskingHash) {
    # Let's individually mask every one of these regions...
    
    foreach my $badRegion (@{$maskingHash{$target}}) {
        my @fields = split(/:/, $badRegion);
        unless ($fields[0] < $fields[1]) {
            die '$fields[0] is less than $fields[1] for target ' . $target . ". We are assuming they're sorted...\n";
        }
        #print "Sequence: " . $seqHash{$target};
        #print "Badbeginning: $fields[0]\n";
        #print "Badend: fields[1]\n";
        my $badLength = $fields[1] - $fields[0] + 1; # The +1 is to include both ends
        #print "Badlength = $badLength\n";
        
        my $replacement = "N" x $badLength; # This will give, for instance, NNNN if $badLength = 4
        substr($seqHash{$target}, $fields[0], $badLength, $replacement);
        #print $seqHash{$target} . "\n";
    }
}


# print Dumper(\%maskingHash);
# print Dumper(\%seqHash);


open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";

foreach my $target (sort keys %seqHash) {
    print $outFH ">$target\n$seqHash{$target}\n";
}

