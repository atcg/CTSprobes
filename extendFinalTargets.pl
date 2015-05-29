#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::SearchIO;

my $help = 0;
my $targetsFile;
my $blastFile;
my $assemblyFile;
my $minLength = 250;
my $outFile;

GetOptions  ("targets=s"      => \$targetsFile,
             "blast=s"        => \$blastFile,
             "assembly=s"     => \$assemblyFile,
             "minlength=i"    => \$minLength,
             "out=s"          => \$outFile,
             "help|man"       => \$help) || pod2usage(2);

if (!$targetsFile or !$blastFile or !$assemblyFile or !$outFile or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}


# First let's put the assembly contigs into a hash:
my %assemblyHash;
my $assemblyIn = Bio::SeqIO->new(-file => $assemblyFile,
                                 -format => 'fasta');
while (my $seq = $assemblyIn->next_seq()) {
    $assemblyHash{$seq->display_id()} = $seq;
}

# Now make a hash of our final list of acceptable targets (depth over 5 and at least 1 SNP)
my %acceptableTargetsHash;
open(my $acceptListFH, "<", $targetsFile) or die "Couldn't open $targetsFile for reading: $!\n";
# The contig name should be in the second column. Like so:
#    "Target"	"AverageTargetDepthAcrossHSP"	"max100AverageTargetDepth"	"Length"	        "informativeSNPs"
#    "1"	    "contig00003|E19A4|OPA"	        12.613540164167	            15.3680791774806	388	2
#    "2"	    "contig01346|E20C7|OPA"	        12.7969556306204        	15.3990884098662	376	3
#    "3"	    "contig03967|E7E8|OPA"	        14.1285259092589        	17.5862664070905	410	8

#"Target"	"AverageTargetDepthAcrossHSP"	"max100AverageTargetDepth"	"SNPpercentage"
#contig00003|E19A4|OPA	12.613540164167	15.3680791774806	0.00515463917525773
#contig01346|E20C7|OPA	12.7969556306204	15.3990884098662	0.00806451612903226
#contig03967|E7E8|OPA	14.1285259092589	17.5862664070905	0.0196078431372549


while (my $line = <$acceptListFH>) {
    chomp($line);
    my @fields = split(/\t/, $line);
    $acceptableTargetsHash{$fields[0]}++;
    #print $fields[1] . "\n";
}

#foreach my $key (sort keys %acceptableTargetsHash) {
#    print "Hash key: " . $key . "\n";
#}

print "List contained a total of " . scalar(keys %acceptableTargetsHash) . " acceptable targets.\n";


# Make file for the output sequences:
open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";


# Now let's iterate over the blast file and pull out the matching sequences
my $searchIO = Bio::SearchIO->new(-file => $blastFile,
                                  -format => 'blasttable');


my $skipCounter = 0;
while (my $result = $searchIO->next_result()) {
    # We're only going to process the best HSP for the best hit (so no looping here)
    my $hit = $result->next_hit;
    my $hsp = $hit->next_hsp;
    
    # Skip the target if it isn't in our acceptable targets list (determined by depth
    # and presence of at least one SNP)
    print "Blast result query name: " . $result->query_name();
    print "Value: " . $acceptableTargetsHash{$result->query_name()} . "\n";
    unless (exists $acceptableTargetsHash{$result->query_name()}) {
        print "Key doesn't exist in line 73!\n";
        next;
    }
    
    # Skip the target if the total assembled contig length is less than $minLength
    if ($assemblyHash{$hit->name()}->length < $minLength) {
        $skipCounter++;
        next;
    }
    
    my $contigMatchStart = $hsp->start("subject");
    my $contigMatchEnd   = $hsp->end("subject");
    if ($contigMatchStart > $contigMatchEnd) {
        die "The starting base is higher than the ending base for contig " . $hit->name() . "\n"
    }
    
    # Now let's check how long the HSP is. If it's over our minimum length, we'll just
    # print the matching CTS sequence from the HSP to the file.
    if ($hsp->length("hit") >= $minLength) {
        print $outFH ">" . $hit->name . "\n";
        print $outFH $assemblyHash{$hit->name()}->subseq($contigMatchStart, $contigMatchEnd);
        print $outFH "\n";
    } else { # If the HSP sequence is NOT long enough, then we'll try to extend it until it is
        my $extraBPneeded = $minLength - $hsp->length("hit");
        
        # First check to make sure that there is at least $extraBPneeded/2 bp available on each side:
        if (  ($contigMatchStart >= ($extraBPneeded/2))  &&  (($assemblyHash{$hit->name()}->length() - $contigMatchEnd) >= ($extraBPneeded/2))  ) {
            # If we reach this block then we have enough sequence on each end
            print $outFH ">" . $hit->name . "\n";
            my $newStartingBase = $contigMatchStart - ($extraBPneeded/2);
            my $newEndingBase = $contigMatchEnd + ($extraBPneeded/2);
            print $outFH $assemblyHash{$hit->name()}->subseq($newStartingBase, $newEndingBase);
            print $outFH "\n";
        } else {
            # If there's not enough sequence on each end ($extraBPneeded/2) then we'll extend
            # as far as possible on one end and make up the difference on the other end.
            # We know we should be able to do this, because we skipped assembled contigs that
            # were less than the $minLength cutoff at the top of the while loop above
            
            # First find which end is the shorter one:
            if ($contigMatchStart > ($assemblyHash{$hit->name()}->length() - $contigMatchEnd)) {
                # If we get here, then we know the end side is the short end
                print $outFH ">" . $hit->name . "\n";
                my $newStartingBase = $assemblyHash{$hit->name()}->length() - $minLength;
                my $newEndingBase = $assemblyHash{$hit->name()}->length();
                print $outFH $assemblyHash{$hit->name()}->subseq($newStartingBase, $newEndingBase);
                print $outFH "\n";
            } else {
                # If we get here, the start side is the short end
                print $outFH ">" . $hit->name . "\n";
                my $newStartingBase = 1;
                my $newEndingBase = $minLength;
                print $outFH $assemblyHash{$hit->name()}->subseq($newStartingBase, $newEndingBase);
                print $outFH "\n";
            }
        }
    }
}

print "Skipped a total of $skipCounter records because the assembled contigs were less than $minLength bp long.\n";



#Documentation
__END__

=head1 NAME

extendFinalTargets.pl

=head1 SYNOPSIS 

perl extendFinalTargets.pl --targets <listOfTargetsFile> --assembly <assemblyFile.fasta> --blast <origTargetsBL2finalContigs.blast> --out <outputFile>

 Options:
   -targets=s       File with a list of targets in final list
   -blast=s         File that contains blast results of original targets to assembly
   -assembly=s      File with the assembly contigs that will be trimmed to final target regions
   -out=s           Output fasta file
   -minlength=i     This is the minimum desired length of target sequences
   -help|man        Prints out documentation


=head1 DESCRIPTION

This program processes the final list of targets (in our case, targets that have at
least one potentially ancestry-informative SNP within the HSP region and have
a depth of at least 5 in the maximally-scoring 100bp window within the HSP region).

It takes the assembled contigs (made from CTS sequence in our case) and pulls out the
matching CTS sequence that corresponds to the original designed target sequence. If
the assembled contig is less than $minLength, then the target sequence is extended
equally in both directions until it is long enough. That final sequence is then
dumped into the output file.


=cut