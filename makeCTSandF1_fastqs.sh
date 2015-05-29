#!/bin/bash

cat ../qc/fastq-join/*F1.Ns.un1.fastq > AllF1.Ns.un1.fastq
cat ../qc/fastq-join/*F1.Ns.un2.fastq > AllF1.Ns.un2.fastq
cat ../qc/fastq-join/*F1.Ns.combinedJoinedAndSingles_trimmed.fastq > AllF1.Ns.combinedJoinedAndSingles_trimmed.fastq
cat ../qc/fastq-join/*CTS.Ns.un1.fastq > AllCTS.Ns.un1.fastq
cat ../qc/fastq-join/*CTS.Ns.un2.fastq > AllCTS.Ns.un2.fastq
cat ../qc/fastq-join/*CTS.Ns.combinedJoinedAndSingles_trimmed.fastq > AllCTS.Ns.combinedJoinedAndSingles_trimmed.fastq