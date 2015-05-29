#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $help = 0;
my $outDir;

GetOptions ("out=s" => \$outDir,
            "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$outDir or $help) {
    die "Must supply --out <outputDirectory>.\n";
}

unless (-d $outDir) {
    mkdir $outDir;
}

my @samples = ( "01A_0801_CTS","01A_0809_F1","01A_0905_BTS","01B_0802_CTS","01B_0810_F1","01B_0906_BTS","01C_0803_CTS","01C_0811_F1","01C_0907_BTS","01D_0804_CTS","01D_0812_F1","01D_0908_BTS","01E_0805_CTS","01E_0901_BTS","01E_0909_BTS","01F_0806_CTS","01F_0902_BTS","01F_0910_DSN_BTS","01G_0807_F1","01G_0903_BTS","01G_0911_DSN_CTS","01H_0808_F1","01H_0904_BTS","02A_0801_F1","02A_0905_BTS","02B_0802_F1","02B_0906_CTS","02C_0803_F1","02C_0907_CTS","02D_0804_F1","02D_0908_CTS","02E_0805_F1","02E_0909_CTS","02F_0806_BTS","02F_0910_CTS","02G_0807_BTS","02G_0911_CTS","02H_0808_BTS","02H_0912_F1","03A_0809_BTS","03A_0905_CTS","03B_0810_BTS","03B_0906_F1","03C_0811_BTS","03C_0907_F1","03D_0812_CTS","03D_0908_F1","03E_0901_CTS","03E_0909_F1","03F_0902_CTS","03F_0910_F1","03G_0903_CTS","03G_0911_F1","03H_0904_CTS","03H_0912_BTS");

foreach my $sample (@samples) {
    my $depthFile = $sample . ".merged.cleaned.sorted.markDups.depth";
    open(my $depthFH, "<", $depthFile) or die "Couldn't open $depthFile for reading: $!\n";
    while (my $line = <$depthFH>) {
        my @fields = split(/\t/, $line);
        if ($fields[0] =~ /.*\_\:\_(.*)\_\:\_/) {
            my $locus = $1;
            $locus =~ s/\|/\-\-/g; # Replace all the |'s with --'s
            my $fileName = $outDir . "/" . $locus . ".depth";
            open(my $outFH, ">>", $fileName) or die "Couldn't open $fileName for appending: $!\n";
            print $outFH "$sample\t$line";
        }
    }
}

