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
my $assembly;
my $outFilePrefix;
my $depthDir;

GetOptions  ("blast=s"        => \$blastFile,
             "assembly=s"     => \$assembly,
             "depthfiledir=s" => \$depthDir,
             "out=s"          => \$outFilePrefix,
             "help|man"       => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$blastFile or !$assembly or !$depthDir or !$outFilePrefix  or $help) {
    die "Must supply --blast <blastoutput> --depthfiledir </path/to/dir/> --assembly <fasta> --out <filePrefix>.\n The blast file should be a blast output file\
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




my @samples = (
               "01A_0801_CTS",
               "01A_0809_F1",
               "01A_0905_BTS",
               "01B_0802_CTS",
               "01B_0810_F1",
               "01B_0906_BTS",
               "01C_0803_CTS",
               "01C_0811_F1",
               "01C_0907_BTS",
               "01D_0804_CTS",
               "01D_0812_F1",
               "01D_0908_BTS",
               "01E_0805_CTS",
               "01E_0901_BTS",
               "01E_0909_BTS",
               "01F_0806_CTS",
               "01F_0902_BTS",
               "01F_0910_DSN_BTS",
               "01G_0807_F1",
               "01G_0903_BTS",
               "01G_0911_DSN_CTS",
               "01H_0808_F1",
               "01H_0904_BTS",
               "02A_0801_F1",
               "02A_0905_BTS",
               "02B_0802_F1",
               "02B_0906_CTS",
               "02C_0803_F1",
               "02C_0907_CTS",
               "02D_0804_F1",
               "02D_0908_CTS",
               "02E_0805_F1",
               "02E_0909_CTS",
               "02F_0806_BTS",
               "02F_0910_CTS",
               "02G_0807_BTS",
               "02G_0911_CTS",
               "02H_0808_BTS",
               "02H_0912_F1",
               "03A_0809_BTS",
               "03A_0905_CTS",
               "03B_0810_BTS",
               "03B_0906_F1",
               "03C_0811_BTS",
               "03C_0907_F1",
               "03D_0812_CTS",
               "03D_0908_F1",
               "03E_0901_CTS",
               "03E_0909_F1",
               "03F_0902_CTS",
               "03F_0910_F1",
               "03G_0903_CTS",
               "03G_0911_F1",
               "03H_0904_CTS",
               "03H_0912_BTS");

my %weightingsHash = (
                      "01A_0801_CTS" =>  0.672849571,
                      "01A_0809_F1" =>  0.639243628,
                      "01A_0905_BTS" =>  0.365032135,
                      "01B_0802_CTS" =>  0.893572736,
                      "01B_0810_F1" =>  0.822493028,
                      "01B_0906_BTS" =>  0.178918215,
                      "01C_0803_CTS" =>  0.71826707,
                      "01C_0811_F1" =>  0.826762193,
                      "01C_0907_BTS" =>  0.17182868,
                      "01D_0804_CTS" =>  0.956316948,
                      "01D_0812_F1" =>  0.687751098,
                      "01D_0908_BTS" =>  0.143539041,
                      "01E_0805_CTS" =>  0.791871027,
                      "01E_0901_BTS" =>  0.902430467,
                      "01E_0909_BTS" =>  0.376379171,
                      "01F_0806_CTS" =>  0.873892117,
                      "01F_0902_BTS" =>  0.999901051,
                      "01F_0910_DSN_BTS" =>  0.537961023,
                      "01G_0807_F1" =>  0.911113481,
                      "01G_0903_BTS" =>  1,
                      "01G_0911_DSN_CTS" =>  0.34160701,
                      "01H_0808_F1" =>  0.751802161,
                      "01H_0904_BTS" =>  0.733224629,
                      "02A_0801_F1" =>  0.535103229,
                      "02A_0905_BTS" =>  0.770791541,
                      "02B_0802_F1" =>  0.385594786,
                      "02B_0906_CTS" =>  0.28030874,
                      "02C_0803_F1" =>  0.616267724,
                      "02C_0907_CTS" =>  0.669865039,
                      "02D_0804_F1" =>  0.323149059,
                      "02D_0908_CTS" =>  0.286276709,
                      "02E_0805_F1" =>  0.738164993,
                      "02E_0909_CTS" =>  0.603938742,
                      "02F_0806_BTS" =>  0.410905105,
                      "02F_0910_CTS" =>  0.331357946,
                      "02G_0807_BTS" =>  0.330243715,
                      "02G_0911_CTS" =>  0.811626187,
                      "02H_0808_BTS" =>  0.290986203,
                      "02H_0912_F1" =>  0.454679297,
                      "03A_0809_BTS" =>  0.337022306,
                      "03A_0905_CTS" =>  0.294045201,
                      "03B_0810_BTS" =>  0.292356386,
                      "03B_0906_F1" =>  0.273154532,
                      "03C_0811_BTS" =>  0.281317849,
                      "03C_0907_F1" =>  0.776689419,
                      "03D_0812_CTS" =>  0.274925807,
                      "03D_0908_F1" =>  0.339460045,
                      "03E_0901_CTS" =>  0.700577468,
                      "03E_0909_F1" =>  0.728416849,
                      "03F_0902_CTS" =>  0.334234769,
                      "03F_0910_F1" =>  0.286825336,
                      "03G_0903_CTS" =>  0.555839238,
                      "03G_0911_F1" =>  0.724592837,
                      "03H_0904_CTS" =>  0.625264738,
                      "03H_0912_BTS" => 0.828501189
                );

my $outFile1 = $outFilePrefix . ".wholeHSPcounts.txt";
my $outFile2 = $outFilePrefix . ".aveDepths.txt";
open(my $outFH1, ">", $outFile1) or die "Couldn't open $outFile1 for writing: $!\n";
open(my $outFH2, ">", $outFile2) or die "Couldn't open $outFile2 for writing: $!\n";


# Print header rows
print $outFH1 "Target\taveUnder5Counts\taveOver5Counts\taveOver10Counts\taveOver15Counts\taveOver20Counts\tmax100AveUnder5Counts\tmax100AveOver5Counts\tmax100AveOver10Counts\tmax100AveOver15Counts\tmax100AveOver20Counts\n";
print $outFH2 "Target\tAverageTargetDepthAcrossHSP\tmax100AverageTargetDepth\n";


my %targetSummaryHash;
foreach my $sample (@samples) {
    my %sampleResultsHash;
    my $depthFile = $depthDir . $sample . ".merged.cleaned.sorted.markDups.depth";
    open(my $depthFH, "<", $depthFile) or die "Couldn't open $depthFile for reading: $!\n";
    
    # First we'll hash-ify all of the data for a sample.
    while (my $line = <$depthFH>) {
        my @fields = split(/\t/, $line);
        my $target = $fields[0];
        my $base = $fields[1];
        my $depth = $fields[2];
        $sampleResultsHash{$target}{$base} = $depth; # Hash of hashes--with the first key being the name of the target and the second being the position in the target
    }
    
    
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
        if (scalar(@bases) >= 100) {
            # Calculate average depth over HSP, and check if there's a stretch at least 100bp with average coverage greater than 5, 10, 20
            my $sumDepthOverHSP = 0;
            foreach my $base (@bases) {                        
                if (exists($sampleResultsHash{$hit->name()}{$base})) {
                    $sumDepthOverHSP += $sampleResultsHash{$hit->name()}{$base};
                } else {
                    # Do nothing, i.e. add 0
                }                        
            }
            my $average = $sumDepthOverHSP / scalar(@bases);

            
            my $weightedAverage = $average * $weightingsHash{$sample};
            if ($weightedAverage > 20) {
                $targetSummaryHash{$hit->name()}{'aveOver20'}++;
                $targetSummaryHash{$hit->name()}{'aveOver15'}++;
                $targetSummaryHash{$hit->name()}{'aveOver10'}++;
                $targetSummaryHash{$hit->name()}{'aveOver5'}++;
            } elsif ($weightedAverage > 15) {
                $targetSummaryHash{$hit->name()}{'aveOver15'}++;
                $targetSummaryHash{$hit->name()}{'aveOver10'}++;
                $targetSummaryHash{$hit->name()}{'aveOver5'}++;
            } elsif ($weightedAverage > 10) {
                $targetSummaryHash{$hit->name()}{'aveOver10'}++;
                $targetSummaryHash{$hit->name()}{'aveOver5'}++;
            } elsif ($weightedAverage > 5) {                        
                $targetSummaryHash{$hit->name()}{'aveOver5'}++;
            } else {
                $targetSummaryHash{$hit->name()}{'aveUnder5'}++;
            }

            push(@{$targetSummaryHash{$hit->name()}{'aveWindowDepth'}}, $weightedAverage);        

            # Need to implement finding the highest coverage for a 100bp window within the best HSP...
            # Go through each 100bp window of the HSP, and calculate the average depth. If it's higher
            # than the current depth, then update the max
            my $max100aveDepth = 0;
            foreach my $baseNum (0 .. (scalar(@bases)-100)) {
                my $maxBase = $bases[$baseNum]+99;
                my $maxStartingBase = $bases[scalar(@bases)-100];
                #print "Current target calc range: $bases[$baseNum] to $maxBase. Maximum starting base = $maxStartingBase\n";
                my $windowDepth = 0;
                my $baseCounter = 0;
                foreach my $innerBase ($bases[$baseNum] .. ($bases[$baseNum+99])) {
                    if (exists $sampleResultsHash{$hit->name()}{$innerBase}) {
                        $windowDepth += $sampleResultsHash{$hit->name()}{$innerBase};
                        $baseCounter++;
                    } else {
                        $windowDepth += 0;
                    }
                }
                my $aveWindowDepth;
                if ($baseCounter > 0) {
                    $aveWindowDepth = $windowDepth / $baseCounter;
                } else {
                    $aveWindowDepth = 0;
                }
                if ($aveWindowDepth > $max100aveDepth) {
                    $max100aveDepth = $aveWindowDepth;
                }
            }
            my $weightedMax100AveDepth = $max100aveDepth * $weightingsHash{$sample};
            print "weightedMax100AveDepth: $weightedMax100AveDepth\n";
            push(@{$targetSummaryHash{$hit->name()}{'max100AveWindowDepth'}}, $weightedMax100AveDepth);   
            
            if ($weightedMax100AveDepth > 20) {
                $targetSummaryHash{$hit->name()}{'max100AveOver20'}++;
                $targetSummaryHash{$hit->name()}{'max100AveOver15'}++;
                $targetSummaryHash{$hit->name()}{'max100AveOver10'}++;
                $targetSummaryHash{$hit->name()}{'max100AveOver5'}++;
            } elsif ($weightedMax100AveDepth > 15) {
                $targetSummaryHash{$hit->name()}{'max100AveOver15'}++;
                $targetSummaryHash{$hit->name()}{'max100AveOver10'}++;
                $targetSummaryHash{$hit->name()}{'max100AveOver5'}++;
            } elsif ($weightedMax100AveDepth > 10) {
                $targetSummaryHash{$hit->name()}{'max100AveOver10'}++;
                $targetSummaryHash{$hit->name()}{'max100AveOver5'}++;
            } elsif ($weightedMax100AveDepth > 5) {                        
                $targetSummaryHash{$hit->name()}{'max100AveOver5'}++;
            } else {
                $targetSummaryHash{$hit->name()}{'max100AveUnder5'}++;
            }
        } else {
            #Do nothing--don't process HSPs under 100bp in length
        }
    }
}

foreach my $target (sort keys %targetSummaryHash) {
    print $outFH1 "$target\t$targetSummaryHash{$target}{'aveUnder5'}\t$targetSummaryHash{$target}{'aveOver5'}\t$targetSummaryHash{$target}{'aveOver10'}\t$targetSummaryHash{$target}{'aveOver15'}\t$targetSummaryHash{$target}{'aveOver20'}\t$targetSummaryHash{$target}{'max100AveUnder5'}\t$targetSummaryHash{$target}{'max100AveOver5'}\t$targetSummaryHash{$target}{'max100AveOver10'}\t$targetSummaryHash{$target}{'max100AveOver15'}\t$targetSummaryHash{$target}{'max100AveOver20'}\n";
    print $outFH2 "$target\t" . sum(@{$targetSummaryHash{$target}{'aveWindowDepth'}})/@{$targetSummaryHash{$target}{'aveWindowDepth'}} . "\t" . sum(@{$targetSummaryHash{$target}{'max100AveWindowDepth'}})/@{$targetSummaryHash{$target}{'max100AveWindowDepth'}} . "\n";
}


