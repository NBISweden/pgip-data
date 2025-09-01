# Small monkeflower dataset

## About

VCF files large enough to be used for variant filtering in exercises.
The data also contains a know locus under selection.

## Genotype calls

Genotype calls for region LG4:11000000-14000000 in
M_aurantiacus_v1_splitline_ordered.fasta. The region contains the
MaMyb2 gene described by Stakowski et al (2019).

Base quality score recalibration was applied to input bam files, using
a first round of raw filtered HaplotypeCaller calls as known sites.
GATK HaplotypeCaller was then run on bqsr-files in GVCF mode. The
results were combined with GATK CombineGVCFs, and joint genotyping was
run with GATK GenotypeGVCFs.

all.variantsites.vcf.gz - all 37 samples

redyellow.variantsites.vcf.gz - 10 samples from red and yellow
ecotypes
