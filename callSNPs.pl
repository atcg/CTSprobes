#!/usr/bin/perl


system("bwa mem -t 30 ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta AllCTS.Ns.un1.fastq AllCTS.Ns.un2.fastq | samtools view -@ 30 -bS - > AllCTS.pe.bam");
system("bwa mem -t 30 ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta AllCTS.Ns.combinedJoinedAndSingles_trimmed.fastq | samtools view -@ 30 -bS - > AllCTS.se.bam");

system("bwa mem -t 30 ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta AllF1.Ns.un1.fastq AllF1.Ns.un2.fastq | samtools view -@ 30 -bS - > AllF1.pe.bam");
system("bwa mem -t 30 ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta AllF1.Ns.combinedJoinedAndSingles_trimmed.fastq | samtools view -@ 30 -bS - > AllF1.se.bam");



system("samtools sort -@ 30 -o AllCTS.se.sorted.bam -T blah AllCTS.se.bam");
system("samtools sort -@ 30 -o AllCTS.pe.sorted.bam -T blah2 AllCTS.pe.bam");

system("samtools sort -@ 30 -o AllF1.se.sorted.bam -T blah3 AllF1.se.bam");
system("samtools sort -@ 30 -o AllF1.pe.sorted.bam -T blah4 AllF1.pe.bam");

system("samtools merge -@ 30 AllCTS.both.bam AllCTS.pe.sorted.bam AllCTS.se.sorted.bam");
system("samtools merge -@ 30 AllF1.both.bam AllF1.pe.sorted.bam AllF1.se.sorted.bam");

system("samtools sort -@ 30 -o AllCTS.both.sorted.bam -T blah5 AllCTS.both.bam ");
system("samtools sort -@ 30 -o AllF1.both.sorted.bam -T blah6 AllF1.both.bam");

system("java -Xmx80g -jar ~/bin/picard-tools-1.125/picard.jar CleanSam I=AllCTS.both.sorted.bam O=AllCTS.cleaned.bam");
system("java -Xmx80g -jar ~/bin/picard-tools-1.125/picard.jar AddOrReplaceReadGroups I=AllCTS.cleaned.bam O=AllCTS.cleaned.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=Lib1 RGID=CTS RGSM=CTS VALIDATION_STRINGENCY=LENIENT");
system("java -Xmx40g -jar ~/bin/picard-tools-1.125/picard.jar MarkDuplicates I=AllCTS.cleaned.RG.bam O=AllCTS.cleaned.RG.MD.bam METRICS_FILE=AllCTS.cleaned.RG.MD.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false");

system("java -Xmx80g -jar ~/bin/picard-tools-1.125/picard.jar CleanSam I=AllF1.both.sorted.bam O=AllF1.cleaned.bam");
system("java -Xmx80g -jar ~/bin/picard-tools-1.125/picard.jar AddOrReplaceReadGroups I=AllF1.cleaned.bam O=AllF1.cleaned.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=Lib1 RGID=F1 RGSM=F1 VALIDATION_STRINGENCY=LENIENT");
system("java -Xmx40g -jar ~/bin/picard-tools-1.125/picard.jar MarkDuplicates I=AllF1.cleaned.RG.bam O=AllF1.cleaned.RG.MD.bam METRICS_FILE=AllF1.cleaned.RG.MD.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false");


system("java -Xmx80g -jar ~/bin/picard-tools-1.125/picard.jar MergeSamFiles SO=coordinate AS=true I=AllCTS.cleaned.RG.MD.bam I=AllF1.cleaned.RG.MD.bam O=CTSandF1.bam");
system("samtools index CTSandF1.bam");

## Realign around indels
system("java -Xmx80g -jar ~/bin/GATK.jar -T RealignerTargetCreator -R ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta -I CTSandF1.bam --minReadsAtLocus 4 -o CTSandF1.intervals");
system("java -Xmx80g -jar ~/bin/GATK.jar -T IndelRealigner -R ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta -I CTSandF1.bam -targetIntervals CTSandF1.intervals -LOD 3.0 -o CTSandF1-realigned.bam");

## Call SNPs
system("java -Xmx80g -jar ~/bin/GATK.jar -T UnifiedGenotyper -R ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta -I CTSandF1-realigned.bam -gt_mode DISCOVERY -stand_call_conf 30 -stand_emit_conf 10 -o CTSandF1-rawSNPS-Q30.vcf");

## annotate SNPs
system("java -Xmx80g -jar ~/bin/GATK.jar -T VariantAnnotator -R ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta -I CTSandF1-realigned.bam -G StandardAnnotation -V:variant,VCF CTSandF1-rawSNPS-Q30.vcf -XA SnpEff -o CTSandF1-rawSNPS-Q30-annotated.vcf");

## call indels
system("java -Xmx80g -jar ~/bin/GATK.jar -T UnifiedGenotyper -R ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta -I CTSandF1-realigned.bam -gt_mode DISCOVERY -glm INDEL -stand_call_conf 30 -stand_emit_conf 10 -o CTSandF1-inDels-Q30.vcf");

## filter SNP calls around indels
system('java -Xmx80g -jar ~/bin/GATK.jar -T VariantFiltration -R ../ARC/finished_allCTS/RBBHs.CTSonly6iter.chimeraMasked.fasta -V CTSandF1-rawSNPS-Q30-annotated.vcf --mask CTSandF1-inDels-Q30.vcf --maskExtension 5 --maskName InDel --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "BadValidation" --filterExpression "QUAL < 30.0" --filterName "LowQual" --filterExpression "QD < 5.0" --filterName "LowVQCBD" -o CTSandF1-Q30SNPs.vcf');
