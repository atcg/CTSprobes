#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;

my $help = 0;
my $configFile;
my $reference;

GetOptions  ("config=s"      => \$configFile,
             "reference=s"   => \$reference,
             "help|man"      => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$configFile or !$reference or $help) {
    die "Must supply --config <file> --reference </path/to/bwaReference.fasta>.\n";
}
my $startingDir = getcwd();
open(my $configFH, "<", $configFile) or die "Couldn't open $configFile for reading: $!\n";

my %params;
while (my $line = <$configFH>) {
    if ($line =~ /^AdapterFolder\=(.*)/) {
        $params{"adapterFolder"} = $1;
    } elsif ($line =~ /^ReadsFolder\=(.*)/) {
        $params{"readsFolder"} = $1;
    } elsif ($line =~ /^QCoutputFolder\=(.*)/) {
        $params{"qcFolder"} = $1;
    } elsif ($line =~ /^ARCworkingFolder\=(.*)/) {
        $params{"arcWorkingFolder"} = $1;
    } elsif ($line =~ /^ARCoutputFolder\=(.*)/) {
        $params{"arcOutputFolder"} = $1;
    } elsif ($line =~ /^threads\=(.*)/) {
        $params{"threads"} = $1;
    } elsif ($line =~ /^Sample\=(.*)/) {
        push(@{$params{"samples"}}, $1);
    } elsif ($line =~ /^\s*$/) { # This is just a line that contains only whitespace (blank line)
        next;
    } elsif ($line =~ /^targets\=(.*)/) {
        $params{"targets"} = $1;
    } elsif ($line =~ /^cycles\=(.*)/) {
        $params{"cycles"} = $1;
    } elsif ($line =~ /^timeout\=(.*)/) {
        $params{"timeout"} = $1;  
    } else {
        die "Invalid config file line:\n$line\n";
    }
}
close($configFH);


foreach my $sample (@{$params{"samples"}}) {
    my $singlesBamFile = "mapping/" . $sample . ".singlesAndJoined.bam";
    my $pairedBamFile = "mapping/" . $sample . ".paired.bam";
    my $mergedBamFile = "mapping/" . $sample . ".merged.bam";
    my $reads1 = $params{"qcFolder"} . "/fastq-join/" . $sample . ".Ns.un1.fastq";
    my $reads2 = $params{"qcFolder"} . "/fastq-join/" . $sample . ".Ns.un2.fastq";
    my $readsSingles = $params{"qcFolder"} . "/fastq-join/" . $sample . ".Ns.combinedJoinedAndSingles_trimmed.fastq";
    system("bwa mem -t $params{'threads'} -M $reference $readsSingles | samtools view -@ $params{'threads'} -bS - > $singlesBamFile");
    system("bwa mem -t $params{'threads'} -M $reference $reads1 $reads2 | samtools view -@ $params{'threads'} -bS - > $pairedBamFile");
    system("samtools merge -@ $params{'threads'} $mergedBamFile $singlesBamFile $pairedBamFile");
    
    # Mark duplicates and use mpileup
    my $cleanedBam = "mapping/" . $sample . ".merged.cleaned.bam";
    my $sortedBam = "mapping/" . $sample . ".merged.cleaned.sorted.bam";
    my $markDupsBam = "mapping/" . $sample . ".merged.cleaned.sorted.markDups.bam";
    my $markDupsMetrics = "mapping/" . $sample . ".merged.cleaned.sorted.markDups.metrics";
    my $depthFile = "mapping/depths/" . $sample . ".merged.cleaned.sorted.markDups.depth";
    system("java -Xmx16g -jar /home/evan/bin/picard-tools-1.125/picard.jar CleanSam I=$mergedBamFile O=$cleanedBam");    
    system("java -Xmx16g -jar /home/evan/bin/picard-tools-1.125/picard.jar AddOrReplaceReadGroups I=$cleanedBam O=$sortedBam SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=Lib1 RGID=$sample RGSM=$sample VALIDATION_STRINGENCY=LENIENT");
    system("java -Xmx16g -jar /home/evan/bin/picard-tools-1.125/picard.jar MarkDuplicates I=$sortedBam O=$markDupsBam METRICS_FILE=$markDupsMetrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 ASSUME_SORTED=true REMOVE_DUPLICATES=false");
    print "Mapping stats for $sample\n";
    system("samtools flagstat mapping/$sample.merged.cleaned.sorted.markDups.bam");
    system("samtools depth -q 20 -Q 20 $markDupsBam > $depthFile");
}

