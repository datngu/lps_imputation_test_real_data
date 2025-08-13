#!/bin/bash

set -e  # Exit immediately if any command fails

# Help message function
function display_help {
    echo "Usage: $0 BAM_LIST CHR OUT_DIR [CORES] [REF_PREFIX]"
    echo "  BAM_LIST: Input list (.txt) of BAM files"
    echo "  CHR: Chromosome or region (e.g., chr1, chr2)"
    echo "  OUT: Output imputed vcf file"
    echo "  CORES (optional): Number of CPU cores to use"
    echo "  REF_PREFIX (optional): Prefix for reference files"
    echo "Example: $0 bam_list.txt chr1 chr1_imputed.vcf.gz 4 binary_ref"
}

# Check for the correct number of arguments
if [ "$#" -lt 3 ]; then
    echo "Error: Insufficient arguments!"
    display_help
    exit 1
fi

# Assign input arguments to variables
BAM_LIST="$1"
CHR="$2"
OUT="$3"
CORES="${4:-1}"  # Default to 1 if CORES is not provided
REF_PREFIX="${5:-binary_ref}"  # Default to 'binary_ref' if REF_PREFIX is not provided
TEM_PREFIX="${OUT%%.*}"

# ## test dir /mnt/ScratchProjects/Aqua-Faang/dat_projects/igsr_data/raw_data/ref
# ## export PATH="$PWD/bin:$PATH"
# ## ln -s /mnt/ScratchProjects/Aqua-Faang/dat_projects/lps_imputation/results/dedup_sorted_bam/HG01761_sorted_dedup* .


# BAM=HG01761_sorted_dedup.bam
# CHR=chr22
# OUT=HG01761_chr22.vcf.gz
# CORES=4
# TEM_PREFIX="${OUT%%.*}"
# REF_PREFIX='binary_ref'


## do imputation
while IFS="" read -r LINE || [ -n "$LINE" ]; 
do   
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	CHR=$(echo ${LINE} | cut -d" " -f2)
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	GLIMPSE2_phase_static --bam-list ${BAM_LIST} \
        --reference ${REF_PREFIX}_${CHR}_${REGS}_${REGE}.bin \
        --output ${TEM_PREFIX}_${CHR}_${REGS}_${REGE}_tempo.bcf \
        --keep-monomorphic-ref-sites \
        --threads ${CORES}
done < chunks.${CHR}.txt


## ligate

LST=list.${CHR}.txt
ls -1v ${TEM_PREFIX}*_tempo.bcf > ${LST}

GLIMPSE2_ligate_static --input ${LST} --output ${OUT}

rm ${TEM_PREFIX}*_tempo*
rm $LST