#!/bin/bash

# Obtain the sequence data
fastq-dump --split-files ERR5743893

# Create a directory to save FASTQC outputs
mkdir -p QC_Reports

# Run fastqc on the two fastq files
fastqc ERR5743893_1.fastq ERR5743893_2.fastq --outdir QC_Reports

# Generate MultiQC Report for all by doing:
"""cd QC_Reports"""
multiqc .

# Directory to store the results from BWA-MEM
mkdir Mapping

# Index the reference genome
bwa index MN908947.fasta

# Mapping sequence from target sample to the reference genome
bwa mem MN908947.fasta ERR5743893_1.fastq ERR5743893_2.fastq > Mapping/ERR5743893.sam

# Convert SAM file to a BAM file
samtools view -@ 20 -S -b Mapping/ERR5743893.sam > Mapping/ERR5743893.bam

# Sort bam file
samtools sort -@ 20 -o Mapping/ERR5743893.sorted.bam Mapping/ERR5743893.bam

# Index the sorted bam file
samtools index Mapping/ERR5743893.sorted.bam

# Index the reference genome
samtools faidx MN908947.fasta