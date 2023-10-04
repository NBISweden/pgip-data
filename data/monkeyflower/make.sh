#!/bin/bash
#
# ABOUT: download bioproject PRJNA549183 sequence data and related
# reference data, map, do variant calling and generate data sets
#
set -eo pipefail

function safe_exit {
    message=$1
    echo $message "$?"
    if [ "$?" != 0 ]; then
	echo $message
	exit
    fi
}

##############################
# bwa index
##############################
function bwa_index {
    ref=$1
    if [ ! -e ${ref}.amb ]; then
	echo "Indexing ${ref}"
	bwa index ${ref}
    else
	echo "Bwa index ${ref} exists; skipping"
    fi
    if [ "$?" != 0 ]; then
	echo "bwa index $ref failed"
	exit
    fi
}

##############################
# Make new ROI reference sequence
##############################
function make_roi_ref {
    chrom=$1
    start=$2
    end=$3
    ROI=$chrom:${start}-${end}
    bedfile=${ROI}.bed
    echo -e "${chrom}\t${start}\t${end}" > $bedfile
    ROIREF=${REF%.fasta}.${ROI}.fasta
    if [ ! -e "${ROIREF}" ]; then
	seqtk subseq ${REF} ${bedfile} > ${ROIREF}
    fi
    bwa_index ${ROIREF}
    if [ "$?" != 0 ]; then
	echo "make_roi_ref $ROI failed"
	exit
    fi
}

function run_samtools_faidx {
    ref=$1
    fai=${ref}.fai
    if [ ! -e $fai ]; then
	samtools faidx $ref
    fi
    if [ "$?" != 0 ]; then
	echo "samtools faidx $ref failed"
	exit
    fi
}

function run_picard_dict {
    ref=$1
    dict=${ref%.fasta}.dict
    if [ ! -e $dict ]; then
	picard CreateSequenceDictionary --REFERENCE $ref
    fi
    if [ "$?" != 0 ]; then
	echo "picard CreateSequenceDictionary $ref failed"
	exit
    fi
}

##############################
# qualimap
##############################
function run_qualimap {
    suffix=$1
    bedfile=$2
    if [ "${bedfile}" == "" ]; then
	gff=""
    else
	gff="-gff $bedfile"
    fi
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	bam=${id}/${id}.${suffix}
	genome_results=${bam%.bam}_stats/genome_results.txt
	if [ ! -e ${genome_results} ]; then
	    echo "Running qualimap on ${bam}"
	    qualimap bamqc -nt 12 $gff -bam $bam
	fi
    done
}

##############################
# mosdepth
##############################
function run_mosdepth {
    suffix=$1
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	bam=${id}/${id}.${suffix}
	mosdepth=${bam%.bam}.mosdepth.summary.txt
	if [ ! -e ${mosdepth} ]; then
	    echo "Running mosdepth on ${bam}"
	    mosdepth -t 14 ${bam%.bam} ${bam}
	fi
    done
}

##############################
# subset_bam_for_roi
#
# Given a region of interest (roi), subset sample bam file to roi
#
##############################
function subset_bam_for_roi {
    ROI=$1
    echo "Generate BAM files for ${ROI}"
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	inbam=${id}/${id}.sort.md.bam
	outbam=${id}/${id}.sort.md.${ROI}.bam
	echo "viewing region ${ROI}, file ${inbam}"
	if [ ! -e "$outbam" ]; then
	    samtools view -P $inbam "$ROI" -h -b -o $outbam
	    samtools index $outbam
	else
	    echo "Bam file $outbam exists; skipping!"
	fi
    done
    run_qualimap "sort.md.${ROI}.bam"
    run_mosdepth "sort.md.${ROI}.bam"
}

##############################
# make_ubam_for_roi
#
# Given a region of interest (roi), subset reads to a unmapped bam
# file
#
# NB: using uBAM format for easier storage:
# https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
#
##############################
function make_ubam_for_roi {
    ROI=$1
    echo "Generate uBAM files for ${ROI}"
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	echo "Processing id ${id}: extract reads from ROI and convert to ubam"
	pfx=${id}/${id}.sort.md.${ROI}
	bamfile=${pfx}.bam
	ubam=${pfx}.unmapped.bam
	if [ ! -e $ubam ]; then
	    picard RevertSam --INPUT ${bamfile} \
		   --OUTPUT ${ubam} \
		   --SANITIZE true \
		   --MAX_DISCARD_FRACTION 0.5 \
		   --ATTRIBUTE_TO_CLEAR XT \
		   --ATTRIBUTE_TO_CLEAR XN \
		   --ATTRIBUTE_TO_CLEAR AS \
		   --ATTRIBUTE_TO_CLEAR OC \
		   --ATTRIBUTE_TO_CLEAR OP \
		   --SORT_ORDER queryname \
		   --RESTORE_ORIGINAL_QUALITIES true \
		   --REMOVE_DUPLICATE_INFORMATION true \
		   --REMOVE_ALIGNMENT_INFORMATION true
	fi
    done
}

##############################
# remap_ubam_to_roi
#
# remap ubam reads to roi reference
##############################
function remap_ubam_to_roi {
    ROI=$1
    ref=$2
    for ((i=0; i<${#srr[@]}; i++)); do
	id=${srr[i]}.${sample[i]}
	pfx=${id}/${id}.sort.md.${ROI}
	ubam=${pfx}.unmapped.bam
	outbam=${pfx}.remap.bam
	rg="-R @RG\\tID:${srr[i]}\\tSM:${sample[i]}\\tPL:ILLUMINA"
	options="-M"
	if [ ! -e $outbam ]; then
	    picard SamToFastq --INPUT ${ubam} --INTERLEAVE true --FASTQ /dev/stdout |
		bwa mem ${rg} -p -t 14 ${options} ${ref} - |
		samtools fixmate -m - /dev/stdout |
		samtools sort - | samtools markdup - /dev/stdout |
		samtools view -h -b -o $outbam
	    samtools index $outbam
	    mosdepth ${pfx}.remap $outbam
	fi
    done
    run_qualimap "sort.md.${ROI}.remap.bam"
    run_mosdepth "sort.md.${ROI}.remap.bam"
}

##################################################
#
# M A I N
#
##################################################
# Download reference genome; bush monkeyflower
# See https://github.com/madeline-chase/Mimulus_genomic_landscape/blob/master/annotation_pipeline_december2016.txt
if [ ! -e M_aurantiacus_v1_splitline_ordered.fasta ]; then
    wget http://mimubase.org/FTP/Genomes/Maurantiacus_v1.0/M_aurantiacus_v1_splitline_ordered.fasta
fi
if [ ! -e MAUR_annotation_functional_submission.gff ]; then
    wget http://mimubase.org/FTP/Genomes/Maurantiacus_v1.0/MAUR_annotation_functional_submission.gff
fi

# Download SraRunInfo.csv for bioproject PRJNA549183; see
if [ ! -e SraRunInfo.csv ]; then
    esearch -db sra -query 'PRJNA549183' |
	efetch -format runinfo |
	awk 'BEGIN {FS=","} {if($13=="WGS" || $13 == "LibraryStrategy") print }' |
	awk 'BEGIN {FS=","} {if($30 != "Assembly Plant") print}' > SraRunInfo.csv
else
    echo "SraRunInfo.csv exists; skipping download"
fi

SRRLIST=$( cat SraRunInfo.csv | grep -v "Run" | awk 'BEGIN {FS=","} {print $1}'  | tr "\n" " " | sed -e "s/ $//")
SAMPLELIST=$(cat SraRunInfo.csv | grep -v "Run" | awk 'BEGIN {FS=","}  {print $30}'  | tr "\n" " " | sed -e "s/ $//")

read -a srr <<< $SRRLIST
read -a sample <<< $SAMPLELIST

if [ ! -e "sampleinfo.csv" ]; then
    cat SraRunInfo.csv | cut -d "," -f 1,29,30 |
	csvtk mutate - -n Sample -f SampleName -p "^[A-Z]+-(.+)" |
	csvtk mutate - -R -n Sample -f Sample -p "^([A-Z].+)_[0-9]+" |
	csvtk mutate - -R -n Sample -f Sample -p "^([^TK1][A-Z]+)[0-9]*-[0-9]+" |
	csvtk replace -f Sample -p "CLV" -r "CLV_" |
	csvtk replace -f Sample -p "GH1" -r "CLV_GH" > sraruninfo.csv
    csvtk join samples.csv sraruninfo.csv -f "Sample" > sampleinfo.csv
fi

##############################
# make bwa index
##############################
REF=M_aurantiacus_v1_splitline_ordered.fasta
bwa_index ${REF}

##############################
# Loop bioproject SRR files for whole genome resequencing and map with
# bwa
#
# See
# https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
# and https://edwards.flinders.edu.au/fastq-dump/
##############################
for ((i=0;i<${#srr[@]};i++)); do
    echo Fetching data for SRR id "${srr[i]}", sample "${sample[i]}"
    id="${srr[i]}.${sample[i]}"
    outbam="${id}/${id}.sort.md.bam"
    if [ ! -e ${outbam} ]; then
	rg="-R @RG\\tID:${srr[i]}\\tSM:${sample[i]}"
	options=""
	prefetch -p ${srr[i]}
	fasterq-dump --skip-technical -Z --split-spot ${srr[i]} |
	    bwa mem ${rg} -p -t 14 ${options} ${REF} - |
	    samtools fixmate -m - /dev/stdout |
	    samtools sort - | samtools markdup - /dev/stdout |
	    samtools view -h -b -o ${outbam}
	samtools index ${outbam}
	rm -rf ${srr[i]}
    else
	echo "Bam file ${outbam} exists; skipping"
    fi
done
run_qualimap "sort.md.bam"
run_mosdepth "sort.md.bam"

##############################
# Select ROI as region surrounding MaMyb2
##############################
start=11000000
end=14000000
chrom=LG4
ROI=${chrom}:${start}-${end}
ROIREF=${REF%.fasta}.${ROI}.fasta
make_roi_ref ${chrom} $start $end
subset_bam_for_roi $ROI
make_ubam_for_roi $ROI
remap_ubam_to_roi $ROI $ROIREF


##############################
# git repo output: subset ubam files to a 100k region. Similar to above.
##############################
start=12000000
end=12100000
chrom=LG4
ROI=${chrom}:${start}-${end}
ROIREF=${REF%.fasta}.${ROI}.fasta
make_roi_ref ${chrom} $start $end
subset_bam_for_roi $ROI
make_ubam_for_roi $ROI
remap_ubam_to_roi $ROI $ROIREF

##############################
# Variant calling workflow
##############################

# 1. filter_input_data - use process_shortreads from stacks
# NB: even when processing the full data set, no reads are removed, so
# presumably filtered reads only have been submitted to ENA

# 2. bwa-mem

# remap.bam are the files to use


# 3. picard mark duplicates
#
# NB: samtools markdup is run as piped command in combination with bwa
# mem so this step is unnecessary

# 4. initial vc with UnifiedGenotyper
#
# NB: here we use HaplotypeCaller.
function run_haplotypecaller {
    inbam=$1
    ref=$2
    options=$3
    mode=$4
    if [ "$mode" == "g.vcf" ]; then
	label="g.vcf"
	options+=("-ERC" "GVCF")
    else
	label="vcf"
    fi
    echo "Running gatk HaplotypeCaller on $inbam"
    outvcf=${inbam%.bam}.hc.${label}.gz
    if [ ! -e ${outvcf} ]; then
	gatk HaplotypeCaller -OVI true "${options[@]}" --input $inbam --output $outvcf --reference $ref
    fi
    safe_exit "gatk HaplotypeCaller $inbam failed"
}
function run_variantfiltration {
    invcf=$1
    options="$2"
    outvcf=${invcf%.vcf.gz}.filtered.vcf.gz
    echo "running variant filtering on $invcf"
    if [ ! -e ${outvcf} ]; then
	gatk VariantFiltration -OVI true --variant $invcf --output $outvcf   "${options[@]}"
    fi
    safe_exit "gatk VariantFiltration $invcf failed"
}
# Loop samples and call steps sequentially
function run_raw_variantcalling {
    ref=$1
    roi=$2
    suffix=$3
    if [ "$suffix" == "" ]; then
	suffix=sort.md.${roi}.remap.bam
    fi
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	inbam=${id}/${id}.${suffix}
	options="-A FisherStrand -A QualByDepth -A MappingQuality -G StandardAnnotation --minimum-mapping-quality 50"
	options=("-A" "FisherStrand" "-A" "QualByDepth" "-A" "MappingQuality" "-G" "StandardAnnotation" "--minimum-mapping-quality" "50")
	run_haplotypecaller $inbam $ref $options
	options=("--filter-name" "FisherStrand")
	options+=("--filter" "FS > 50.0")
	options+=("--filter-name" "QualByDepth")
	options+=("--filter" "QD < 4.0")
	options+=("--filter-name" "MappingQuality")
	options+=("--filter" "MQ < 50.0")
	run_variantfiltration ${inbam%.bam}.hc.vcf.gz "$options"
    done
}

# 5. bqsr
function run_bqsr {
    ref=$1
    roi=$2
    suffix=$3
    if [ "$suffix" == "" ]; then
	suffix=sort.md.${roi}.remap.hc.filtered.vcf.gz
    fi
    vcflist=()
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	vcf=${id}/${id}.${suffix}
	vcflist+=($vcf)
    done
    knownsites=$(echo knownsites.${roi}.vcf.gz | sed -e "s/[:-]/_/g")
    if [ ! -e $knownsites ]; then
	bcftools merge "${vcflist[@]}" -O z -o $knownsites
	tabix $knownsites
    fi
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	inbam=${id}/${id}.${suffix%.hc.filtered.vcf.gz}.bam
	recal=${inbam%.bam}.recal.table
	if [ ! -e $recal ]; then
	    gatk BaseRecalibrator -I $inbam -R $ref --known-sites $knownsites -O $recal
	fi
	safe_exit "gatk BaseRecalibrator $inbam failed"
	recalbam=${inbam%.bam}.recal.bam
	if [ ! -e $recalbam ]; then
	    gatk ApplyBQSR -bqsr $recal -I $inbam -O $recalbam
	fi
	safe_exit "gatk ApplyBQSR $inbam failed"
    done
}

# 6. unifiedgenotyper - after bqsr
function run_bqsr_variantcalling {
    ref=$1
    roi=$2
    suffix=$3
    if [ "$suffix" == "" ]; then
	suffix=sort.md.${roi}.remap.recal.bam
    fi
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	inbam=${id}/${id}.${suffix}
	options=("-A" "FisherStrand" "-A" "QualByDepth" "-A" "MappingQuality" "-G" "StandardAnnotation" "--minimum-mapping-quality" "50")
	run_haplotypecaller $inbam $ref $options "g.vcf"
	options=("--filter-name" "FisherStrand" "--filter" "FS > 60.0" "--filter-name" "QualByDepth" "--filter" "QD < 2.0" "--filter-name" "MappingQuality" "--filter" "MQ < 40.0")
	# Gvcf mode
	# run_variantfiltration ${inbam%.bam}.hc.vcf.gz $options
    done
}

# 7. unifiedgenotyper with EMIT_ALL_SITES - could compare results from
# manual coverage filtering and all sites vcf
function run_genotype_gvcfs {
    ref=$1
    roi=$2
    options=$3
    suffix=$4
    if [ "$suffix" == "" ]; then
	suffix=sort.md.${roi}.remap.recal.hc.g.vcf.gz
    fi
    vcflist=()
    for ((i=0;i<${#srr[@]};i++)); do
	id=${srr[i]}.${sample[i]}
	vcf=${id}/${id}.${suffix}
	vcflist+=("-V" $vcf)
    done
    combinevcf=$(echo "${roi}.combine.g.vcf.gz"| sed -e "s/[:\-]/_/g")
    if [ ! -e ${combinevcf} ]; then
	gatk CombineGVCFs -OVI true --output $combinevcf --reference $ref ${vcflist[@]}
    fi
    outvcf=$(echo "${roi}.allsites.vcf.gz"| sed -e "s/[:\-]/_/g")
    if [ ! -e ${outvcf} ]; then
	gatk GenotypeGVCFs -OVI true -R $ref -V $combinevcf -O $outvcf ${options[@]}
    fi
}

# Small ROI
start=12000000
end=12100000
chrom=LG4
ROI=${chrom}:${start}-${end}
ROIREF=${REF%.fasta}.${ROI}.fasta
# run_samtools_faidx $ROIREF
# run_picard_dict $ROIREF
# run_raw_variantcalling $ROIREF $ROI
# run_bqsr $ROIREF $ROI
# run_bqsr_variantcalling $ROIREF $ROI
# options=("--all-sites")
# run_genotype_gvcfs $ROIREF $ROI $options

# Large ROI
start=11000000
end=14000000
chrom=LG4
ROI=${chrom}:${start}-${end}
ROIREF=${REF%.fasta}.${ROI}.fasta
run_samtools_faidx $ROIREF
run_picard_dict $ROIREF
run_raw_variantcalling $ROIREF $ROI
run_bqsr $ROIREF $ROI
run_bqsr_variantcalling $ROIREF $ROI
options=("--all-sites")
run_genotype_gvcfs $ROIREF $ROI $options


##############################
# Phylogenetic analyses
##############################

# 1. phasing with BEAGLE

# 2. MVFtools to run raxml

# 3. tree visualisation in DensiTree

# 4. ape to assess concordance among trees

# 5. pca in plink

##############################
# Population genomics workflow
##############################

# 1. pi, dxy, fst using simon martin's tools

