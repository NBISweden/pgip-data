import pandas as pd

# Read sampleinfo and subset to
sampleinfo = pd.read_csv("sampleinfo.csv")
sampleinfo = sampleinfo.loc[sampleinfo.SampleAlias.str.startswith("PUN")]\
                       .set_index("SampleAlias")
samples = sampleinfo.index.values

REFERENCE = "ref/M_aurantiacus_v1.fasta"

rule all:
    input: "multiqc_report.html",

rule samtools_faidx:
    """Run samtools faidx on reference"""
    output: REFERENCE + ".fai",
    input: REFERENCE,
    shell: """samtools faidx {input}"""

BWA_INDEX_SUFFIX = ["amb", "ann", "bwt", "pac", "sa"]

rule bwa_index:
    """Bwa index the reference"""
    output: expand(REFERENCE + ".{suffix}", suffix=BWA_INDEX_SUFFIX),
    input: REFERENCE
    shell: """bwa index {input}"""

rule fastqc:
    """Run fastqc"""
    output: directory("fastqc/{sample}_{read}_fastqc"),
    input: "fastq/{sample}_{read}.fastq.gz",
    shell: """fastqc --extract -o fastqc {input}"""

rule bwa_mem:
    """Map reads to reference"""
    output: "bam/{sample}.bam",
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
        index=expand(REFERENCE + ".{suffix}", suffix=BWA_INDEX_SUFFIX),
    params:
        rgid=lambda wildcards: sampleinfo.loc[wildcards.sample].Run,
        ref=REFERENCE,
    threads: 4
    shell:
        """
        bwa mem -R "@RG\\tID:{params.rgid}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
        -t {threads} -M {params.ref} {input.r1} {input.r2} | \
        samtools sort - | \
        samtools view --with-header --output {output}
        """

rule qualimap_bamqc:
    """Run qualimap bamqc on a bam file"""
    output: "qualimap/{sample}_stats/genome_results.txt",
    input: "bam/{sample}.bam",
    shell:
        """
        unset DISPLAY; qualimap bamqc -bam {input} -outdir qualimap/{wildcards.sample}_stats
        """
        
rule picard_create_sequence_dictionary:
    output: REFERENCE.replace(".fasta", ".dict"),
    input: REFERENCE,
    shell: """picard CreateSequenceDictionary R={input} O={output}"""

rule picard_mark_duplicates:
    output:
        bam="md/{sample}.bam",
        bai="md/{sample}.bai",
        metrics="md/{sample}.dup_metrics.txt",
    input: "bam/{sample}.bam",
    shell:
        """
        picard MarkDuplicates --INPUT {input} --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics} --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
        """

rule gatk_haplotypecaller_raw:
    """Run GATK HaplotypeCaller in GVCF mode to generate raw variants"""
    output:
        vcf="gatk-hc-raw/{sample}.g.vcf.gz",
        tbi="gatk-hc-raw/{sample}.g.vcf.gz.tbi",
    input:
        bam="md/{sample}.bam",
        bai="md/{sample}.bai",
        ref=REFERENCE,
        dict=REFERENCE.replace(".fasta", ".dict"),
        fai=REFERENCE + ".fai",
    threads: 4
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -OVI true --emit-ref-confidence GVCF \
        --annotation FisherStrand -A QualByDepth -A MappingQuality -G StandardAnnotation \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        --native-pair-hmm-threads {threads}
        """

rule gatk_base_recalibrator:
    """Run GATK BaseRecalibrator to generate recalibration table"""
    output: "gatk-bqsr/{sample}.table",
    input:
        bam="md/{sample}.bam",
        bai="md/{sample}.bai",
        ref=REFERENCE,
        known_sites=expand("gatk-hc-raw/{sample}.g.vcf.gz", sample=samples),
    params:
        known_sites=lambda wildcards, input: ' --known-sites '.join(input.known_sites),
    shell:
        """
        gatk --java-options "-Xmx4g" BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --known-sites \
        {params.known_sites} \
        -O {output}
        """

rule gatk_apply_bqsr:
    """Run GATK ApplyBQSR to apply base quality score recalibration"""
    output:
        bam="gatk-bqsr/{sample}.bam",
        bai="gatk-bqsr/{sample}.bai",
    input:
        bam="md/{sample}.bam",
        bai="md/{sample}.bai",
        ref=REFERENCE,
        table="gatk-bqsr/{sample}.table",
    shell:
        """
        gatk --java-options "-Xmx4g" ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        --bqsr-recal-file {input.table} \
        -O {output.bam}
        """

rule gatk_haplotypecaller_bqsr:
    """Run GATK HaplotypeCaller in GVCF mode on BQSR'd BAM to generate raw variants"""
    output:
        vcf="gatk-hc-bqsr/{sample}.g.vcf.gz",
        tbi="gatk-hc-bqsr/{sample}.g.vcf.gz.tbi",
    input:
        bam="gatk-bqsr/{sample}.bam",
        bai="gatk-bqsr/{sample}.bai",
        ref=REFERENCE,
        dict=REFERENCE.replace(".fasta", ".dict"),
        fai=REFERENCE + ".fai",
    threads: 4
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -OVI true --emit-ref-confidence GVCF \
        --annotation FisherStrand -A QualByDepth -A MappingQuality -G StandardAnnotation \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        --native-pair-hmm-threads {threads}
        """

rule gatk_combine_gvcfs:
    """Combine GVCFs from all samples into a single GVCF"""
    output:
        vcf="gatk-combine-gvcfs/combined.g.vcf.gz",
        tbi="gatk-combine-gvcfs/combined.g.vcf.gz.tbi",
    input:
        samples=expand("gatk-hc-bqsr/{sample}.g.vcf.gz", sample=samples),
        ref=REFERENCE,
    params:
        variants=lambda wildcards, input: expand("-V gatk-hc-bqsr/{sample}.g.vcf.gz", sample=samples),
    shell:
        """
        gatk --java-options "-Xmx4g" CombineGVCFs \
        -R {REFERENCE} {params.variants} \
        -O {output.vcf}
        """

rule gatk_genotype_gvcfs:
    """Genotype combined GVCF to produce raw VCF"""
    output:
        vcf="gatk-genotype-gvcfs/allsites.vcf.gz",
        tbi="gatk-genotype-gvcfs/allsites.vcf.gz.tbi",
    input:
        gvcf="gatk-combine-gvcfs/combined.g.vcf.gz",
        ref=REFERENCE,
        dict=REFERENCE.replace(".fasta", ".dict"),
        fai=REFERENCE + ".fai",
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R {input.ref} \
        -V {input.gvcf} \
        -O {output.vcf} \
        --all-sites
        """

rule bcftools_stats:
    """Calculate basic stats from vcf file."""
    output: "stats/{prefix}.stats.txt",
    input:
        vcf="gatk-genotype-gvcfs/{prefix}.vcf.gz",
        tbi="gatk-genotype-gvcfs/{prefix}.vcf.gz.tbi",
    shell:
        """
        bcftools stats {input.vcf} > {output}
        """

rule multiqc:
    """Run MultiQC on all data"""
    output: "multiqc_report.html",
    input: 
        fastqc_r1=expand("fastqc/{sample}_R1_fastqc", sample=samples),
        fastqc_r2=expand("fastqc/{sample}_R2_fastqc", sample=samples),
        stats=expand("stats/{prefix}.stats.txt", prefix=["allsites"]),
        md=expand("md/{sample}.dup_metrics.txt", sample=samples),
        qualimap=expand("qualimap/{sample}_stats/genome_results.txt", sample=samples)
    shell:
        """
        multiqc -f fastqc stats md qualimap
        """
