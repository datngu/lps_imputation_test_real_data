#!/bin/bash

set -e  # Exit immediately if any command fails


BAM="$1"
OUT_PREFIX="$2"
CORES="${3:-1}"  # Default to 1 if CORES is not provided

# Help message function
function display_help {
    echo "Usage: $0 BAM_FILE OUT_PREFIX [CORES]"
    echo "  BAM_FILE: Input BAM file"
    echo "  OUT_PREFIX: Output imputed vcf file prefix, chr# will be automaticly tagged"
    echo "  CORES (optional): Number of CPU cores to use"
    echo "Example: $0 input.bam input_test 4"
}

# Check for the correct number of arguments
if [ "$#" -lt 2 ]; then
    echo "Error: Insufficient arguments!"
    display_help
    exit 1
fi


samtools index $BAM

for i in {1..22}
do
    run_imputation.sh \
        ${BAM} \
        chr${i} \
        ${OUT_PREFIX}.chr${i}_imputed.vcf.gz \
        ${CORES}
done
