#!/bin/bash

# The header of bam files from public source were different with our own bam files, which resulted in troubles in following analysis
steps. This shell script transformed the bam files to PE fastq files and move corresponding fastq files derived from a same sample to
a same directory.

for i in $(ls *.bam)
do
    samtools sort -n $i ${i:0:5}_sorted #sort queries according to query name.
    bedtools bamtofastq -i ${i:0:5}_sorted.bam -fq ${i:0:5}_R1.fq -fq2 ${i:0:5}_R2.fq
    mkdir ${i:0:5}_fqs
    mv ${i:0:5}_R1.fq ${i:0:5}_R2.fq ${i:0:5}_fqs
    rm $i
    rm ${i:0:5}_sorted.bam
done
