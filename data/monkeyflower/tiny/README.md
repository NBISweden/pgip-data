# Tiny monkeyflower dataset

## About

The tiny data set has been used to test workflows and to generate
course web pages.

## Genotype calls

Genotype calls are for region LG4:12000000-12100000 in
M_aurantiacus_v1_splitline_ordered.fasta.

Base quality score recalibration was applied to input bam files, using
a first round of raw filtered HaplotypeCaller calls as known sites.
GATK HaplotypeCaller was then run on bqsr-files in GVCF mode. The
results were combined with GATK CombineGVCFs, and joint genotyping was
run with GATK GenotypeGVCFs.

## Files and directories

### ubam

Unmapped BAM files for region used as input to generate FASTQ files.

### fastq

FASTQ files for 37 samples. Reads have been pre-mapped to genome and
extracted for the test region, such that remapping will not generate
any unmapped reads.

### gatk-hc-bqsr

HaplotypeCaller GVCF files called on BQSR BAM files.

### gatk-combine-gvcf-bqsr

Combined GVCF files called on BQSR BAM files. Contains results for all
(n=37) and redyellow (n=10) samples.

### gatk-genotype-gvcf-bqsr

Genotype GVCF files from above. Contains results for all (n=37) and
redyellow (n=10) samples.

### ref

Reference genome files trimmed to region.

### rm

Repeat library files for region.
