#!/bin/bash

INPUT_DIRECTORY="unset"

usage() {
    echo "Usage: $0 [-d input_directory] rsync_options"
    exit 2
}

while getopts "d:h" opt; do
    case $opt in
    h)
        usage
        ;;
    d)
        INPUT_DIRECTORY=$OPTARG
        shift $((OPTIND - 1))
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        usage
        ;;
    esac
done

if [ "$INPUT_DIRECTORY" = "unset" ]; then
    echo "Error: input directory not set"
    usage
fi

OPTIONS=("$@")

echo "${OPTIONS[@]}"
echo running rsync -amv "${OPTIONS[@]}" "${INPUT_DIRECTORY}" ...

rsync -amv "${OPTIONS[@]}" \
    --include="*/" \
    --include="*.variantsites*.vcf.gz" \
    --include="*.variantsites*.vcf.gz.tbi" \
    --exclude="*.subset.vcf.gz*" \
    --exclude="*" \
    "${INPUT_DIRECTORY}/gatk-genotype-gvcf-bqsr" \
    "${INPUT_DIRECTORY}/vcftools-filter-bqsr" \
    .
