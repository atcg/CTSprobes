#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

my $help = 0;
my $assembly;
my $targets;
my $outFile;

GetOptions ("assembly=s"    => \$assembly,
            "targets=s"     => \$targets,
            "out=s"         => \$outFile,
            "help|man"      => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

#GetOptions  ("assemblybl2targ=s"      => \$assemblybl2targ,
#             "targbl2assembly=s"      => \$targbl2assembly,
#             "assemblyseqs=s"         => \$assemblySeqs,
#             "out=s"                  => \$outFile,
#             "help|man"               => \$help)           || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$targets or !$assembly or !$outFile or $help) {
    die "Must supply --assembly and --targets and --out.\nBlast commands should look like\
    blastn -db assembly -query targets.fasta -outfmt 6 -max_target_seqs 1 -out targets_bl2_assembly.blast and \
    blastn -db targets -query assembly.fasta -outfmt 6 -max_target_seqs 1 -out assembly_bl2_targets.blast";
    
}

system("makeblastdb -in $assembly -dbtype nucl -out findRBBHassembly");
system("makeblastdb -in $targets -dbtype nucl -out findRBBHtargets");
system("blastn -db findRBBHassembly -query $targets -outfmt 6 -max_target_seqs 1 -out findRBBH_targets_bl2_assembly.blast");
system("blastn -db findRBBHtargets -query $assembly -outfmt 6 -max_target_seqs 1 -out findRBBH_assembly_bl2_targets.blast");


# Make an index of all the sequences in the assembly
my %assemblyHash;
my $assemblyIn = Bio::SeqIO->new(-file => $assembly,
                                -format => 'fasta');
while (my $seq = $assemblyIn->next_seq()) {
    $assemblyHash{$seq->display_id()} = $seq;
}


# Find the reciprocal best blast hits
my $assembly2targetsBlast = Bio::SearchIO->new(-file => "findRBBH_assembly_bl2_targets.blast",
                                            -format => 'blasttable');

my %assem2targBlastResults;
while (my $result = $assembly2targetsBlast->next_result()) {
    $assem2targBlastResults{$result->query_name()} = ($result->next_hit())->name(); # This finds the name of the first hit and makes that the value for the key that is the name of the query
}



my $targ2assembBlast = Bio::SearchIO->new(-file => "findRBBH_targets_bl2_assembly.blast",
                                        -format => 'blasttable');
my %targ2assemBlastResults;
while (my $result = $targ2assembBlast->next_result()) {
    # print $result->query_name() . "\n";
    $targ2assemBlastResults{$result->query_name()} = ($result->next_hit())->name();  # This finds the name of the first hit and makes that the value for the key that is the name of the query
}

my $seqOut = Bio::SeqIO->new(-file=>">$outFile",
                            -format => 'fasta');

my $counter = 0;
foreach my $seqName (sort keys %targ2assemBlastResults) {
    # $seqName should be something like "Contig63"
    # We want to see if $assem2targBlastResults{[the name of the contig in the assembly that was best match for target]} is the same as the name of the contig in the targets that matched the assembly
    if ($assem2targBlastResults{$targ2assemBlastResults{$seqName}} eq $seqName) {
        my $assembledContigName = $1 if ($assemblyHash{$targ2assemBlastResults{$seqName}}->display_id() =~ /\_\:\_(.*)\_\:\_/);
        next unless ($assembledContigName eq $seqName);
        $counter++;

        print "Name 1: " . $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id() . ". Name 2: " . $seqName . "\n";
	# print "test: " . $targ2assemBlastResults{$seqName} . "\n";        
        my $seqDisplayName = $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id();
        my $newSeqName = $seqDisplayName . "_$seqName";
        $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id($newSeqName);
        $seqOut->write_seq($assemblyHash{$targ2assemBlastResults{$seqName}});
    }
}
print "$counter total RBB hits found for assembly\n";

unlink("findRBBH_targets_bl2_assembly.blast", "findRBBH_assembly_bl2_targets.blast", "findRBBHassembly.nsq", "findRBBHassembly.nin", "findRBBHassembly.nhr", "findRBBHtargets.nsq", "findRBBHtargets.nin", "findRBBHtargets.nhr");

