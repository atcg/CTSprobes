#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;

my $help = 0;
my $configFile;

GetOptions  ("config=s"      => \$configFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$configFile or $help) {
    die "Must supply --config <file>.\n";
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

unless (-d $params{"qcFolder"}) {
    mkdir $params{"qcFolder"} or die qq{Couldn't make output directory $params{"qcFolder"}: $!\n};
}

############################################################################
# Remove all CASAVA-filtered reads:
############################################################################
my $casavaFilteredDir = $params{"qcFolder"} . "/casavaFilt";
unless (-d $casavaFilteredDir) {
    mkdir $casavaFilteredDir;
}
my @fastqIlluminaFilterCommands;
foreach my $sample (@{$params{"samples"}}) {
    my $inFileNameR1 = $params{"readsFolder"} . "/" . $sample . ".R1.fastq";
    my $inFileNameR2 = $params{"readsFolder"} . "/" . $sample . ".R2.fastq";
    my $outFileNameR1 = $casavaFilteredDir . "/" . $sample . ".Ns.R1.fastq";
    my $outFileNameR2 = $casavaFilteredDir . "/" . $sample . ".Ns.R2.fastq";
    push(@fastqIlluminaFilterCommands, "fastq_illumina_filter -vN $inFileNameR1 > $outFileNameR1");
    push(@fastqIlluminaFilterCommands, "fastq_illumina_filter -vN $inFileNameR2 > $outFileNameR2");
}
my $forkManagerCasava = new Parallel::ForkManager($params{"threads"});
foreach my $casavaCommand (@fastqIlluminaFilterCommands) {
    print "Running the following command: \n$casavaCommand\n";
    sleep 5;
    $forkManagerCasava->start and next;
    system("$casavaCommand");
    print "Finished running the following:\n$casavaCommand\n\n";
    $forkManagerCasava->finish;
}
$forkManagerCasava->wait_all_children;
print "--------------------------------------------------\n";
print "Finished running all fastq_illumina_filter (CASAVA filter) commands\n";
print "--------------------------------------------------\n\n";


############################################################################
# Run Trimmomatic on all of the samples.
############################################################################
my $trimmomaticDir = $params{"qcFolder"} . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}

my @trimmomaticCommands;
foreach my $sample (@{$params{"samples"}}) {
    my $R1 = $casavaFilteredDir . "/$sample" . ".Ns.R1.fastq";
    my $R2 = $casavaFilteredDir . "/$sample" . ".Ns.R2.fastq";
    my $adapterFile = $params{"adapterFolder"} . "/" . $sample . ".adapters";
    
    my $R1OutFilePaired = "$trimmomaticDir/$sample" . ".Ns.R1_paired_trimmed.fastq";
    my $R2OutFilePaired = "$trimmomaticDir/$sample" . ".Ns.R2_paired_trimmed.fastq";
    my $R1OutFileSingles = "$trimmomaticDir/$sample" . ".Ns.R1_singles_trimmed.fastq";
    my $R2OutFileSingles = "$trimmomaticDir/$sample" . ".Ns.R2_singles_trimmed.fastq";
    
    push(@trimmomaticCommands, "java -Xmx8G -jar ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads $params{'threads'} -phred33 $R1 $R2 $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adapterFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");
}
print "Running all Trimmomatic commands\n";
foreach my $trimmomaticCommand(@trimmomaticCommands) {
    print "Running the following command: \n$trimmomaticCommand\n";
    system($trimmomaticCommand);
}
print "Finished running all Trimmomatic commands\n";


############################################################################
# Run fastq-join on all the samples.
############################################################################
my $fastqJoinDir = $params{"qcFolder"} . "/fastq-join";
unless (-d $fastqJoinDir) {
    mkdir $fastqJoinDir or die "Couldn't make directory $fastqJoinDir: $!\n";
}

my @fastqJoinCommands;
chdir($fastqJoinDir);
foreach my $sample (@{$params{"samples"}}) {
    my $fastqJoinOutputName = $sample . ".Ns.%.fastq";
    push(@fastqJoinCommands,"fastq-join -v ' ' ../trimmomatic/$sample.Ns.R1_paired_trimmed.fastq ../trimmomatic/$sample.Ns.R2_paired_trimmed.fastq -o $fastqJoinOutputName");
}
my $forkManagerFastqJoin = new Parallel::ForkManager($params{"threads"});
foreach my $fastqJoinCommand (@fastqJoinCommands) {
    print "Running the following command: \n$fastqJoinCommand\n";
    sleep 5;    
    $forkManagerFastqJoin->start and next;
    system($fastqJoinCommand);
    print "Finished running the following command: \n$fastqJoinCommand\n";
    $forkManagerFastqJoin->finish;
}
$forkManagerFastqJoin->wait_all_children;

chdir($startingDir); # Go back to the starting directory


############################################################################
# Combine the singleton and joined read files
############################################################################
foreach my $sample (@{$params{"samples"}}) { # Don't need threading here because this is fast
    my $combinedFile = $fastqJoinDir . "/$sample.Ns.combinedJoinedAndSingles_trimmed.fastq";
    my $singlesFile1 = "$trimmomaticDir/$sample.Ns.R1_singles_trimmed.fastq"; 
    my $singlesFile2 = "$trimmomaticDir/$sample.Ns.R2_singles_trimmed.fastq";
    my $joinedFile   = "$fastqJoinDir/$sample.Ns.join.fastq";
    print "Running the following command: \ncat $joinedFile $singlesFile1 $singlesFile2 > $combinedFile\n\n";
    system("cat $joinedFile $singlesFile1 $singlesFile2 > $combinedFile");
    unlink($singlesFile1, $singlesFile2, $joinedFile);
}


############################################################################
# Generate the CTS reads files and the ARC config file
############################################################################
unless(-d "ARC") {
    mkdir "ARC";
}
my $CTSreadsCommandJoined = "cat ";
my $CTSreadsCommandUn1= "cat ";
my $CTSreadsCommandUn2 = "cat ";
foreach my $sample (@{$params{"samples"}}) {
    if ($sample =~ /\_CTS$/) {
        $CTSreadsCommandJoined = $CTSreadsCommandJoined . "$fastqJoinDir/$sample.Ns.combinedJoinedAndSingles_trimmed.fastq ";
        $CTSreadsCommandUn1 = $CTSreadsCommandUn1 . "$fastqJoinDir/$sample.Ns.un1.fastq ";
        $CTSreadsCommandUn2 = $CTSreadsCommandUn2 . "$fastqJoinDir/$sample.Ns.un2.fastq ";
    }
}
$CTSreadsCommandJoined = $CTSreadsCommandJoined . "> ARC/AllCTS.Ns.combinedJoinedAndSingles_trimmed.fastq";
$CTSreadsCommandUn1 = $CTSreadsCommandUn1 . "> ARC/AllCTS.Ns.un1.fastq";
$CTSreadsCommandUn2 = $CTSreadsCommandUn2 . "> ARC/AllCTS.Ns.un2.fastq";
print "Making the AllCTS read files\n";
print "Running the following command:\n$CTSreadsCommandJoined";
system($CTSreadsCommandJoined);
print "Running the following command:\n$CTSreadsCommandUn1";
system($CTSreadsCommandUn1);
print "Running the following command:\n$CTSreadsCommandUn2";
system($CTSreadsCommandUn2);

# Write the config file
open(my $ARCconfigFH, ">", "ARC/ARC_config.txt") or die "Couldn't open ARC/ARC_config.txt\n";
print $ARCconfigFH "## Name=value pairs:\n";
print $ARCconfigFH "## reference: contains reference sequences in fasta format\n";
print $ARCconfigFH "## numcycles: maximum number of times to try remapping\n";
print $ARCconfigFH "## mapper: the mapper to use (blat/bowtie2)\n";
print $ARCconfigFH "## assembler: the assembler to use (newbler/spades)\n";
print $ARCconfigFH "## nprocs: number of cores to use\n";
print $ARCconfigFH "## format: fasta or fasta, all must be the same\n";
print $ARCconfigFH "## verbose: control mapping/assembly log generation (True/False)\n";
print $ARCconfigFH "## urt: For Newbler, enable use read tips mode (True/False)\n";
print $ARCconfigFH "## map_against_reads: On iteration 1, skip assembly, map against mapped reads (True/False)\n";
print $ARCconfigFH "## assemblytimeout: kill assemblies and discard targets if they take longer than N minutes\n";
print $ARCconfigFH "##\n";
print $ARCconfigFH "## Columns:\n";
print $ARCconfigFH "## Sample_ID:Sample_ID\n";
print $ARCconfigFH "## FileName: path for fasta/fasta file\n";
print $ARCconfigFH "## FileType: PE1, PE2, or SE\n";
print $ARCconfigFH "## FileFormat: fasta or fasta\n";
print $ARCconfigFH "# reference=$startingDir/$params{'targets'}\n";
print $ARCconfigFH "# numcycles=$params{'cycles'}\n";
print $ARCconfigFH "# mapper=bowtie2\n";
print $ARCconfigFH "# assembler=spades\n";
print $ARCconfigFH "# nprocs=$params{'threads'}\n";
print $ARCconfigFH "# format=fastq\n";
print $ARCconfigFH "# verbose=True\n";
print $ARCconfigFH "# urt=True\n";
print $ARCconfigFH "# map_against_reads=False\n";
print $ARCconfigFH "# assemblytimeout=$params{'timeout'}\n";
print $ARCconfigFH "# bowtie2_k=5\n";
print $ARCconfigFH "# rip=True\n";
print $ARCconfigFH "# cdna=False\n";
print $ARCconfigFH "# subsample=1\n";
print $ARCconfigFH "# maskrepeats=True\n";
print $ARCconfigFH "# workingdirectory=/home/evan/ramdisk/\n"; # Make this configurable(?)
print $ARCconfigFH "Sample_ID	FileName	FileType\n";
print $ARCconfigFH "allCTS	AllCTS.Ns.un1.fastq	PE1\n";
print $ARCconfigFH "allCTS	AllCTS.Ns.un2.fastq	PE2\n";
print $ARCconfigFH "allCTS	AllCTS.Ns.combinedJoinedAndSingles_trimmed.fastq	SE\n";
close($ARCconfigFH);

# Run ARC
chdir("ARC");
print "\n\nRunning the ARC assembly pipeline.\n";
system("ARC -c ARC_config.txt > ARC_allCTS.log 2>&1");
print "Finished running ARC. All done\n";
chdir($startingDir);





