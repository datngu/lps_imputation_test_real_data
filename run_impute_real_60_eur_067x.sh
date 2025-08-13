#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=imp_real067   
#SBATCH --mem=4G                
#SBATCH --mail-user=n.dat@outlook.com
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/24.04.2
module load singularity/rpm


export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity


nextflow run bam2vcf.nf \
    -w 'work_bam2vcf' \
    --bam '/mnt/ScratchProjects/Aqua-Faang/dat_projects/test_real_data/results/dedup_sorted_bam/*.bam{,.bai}' \
    --ref_vcf '/mnt/ScratchProjects/Aqua-Faang/dat_projects/igsr_data/ref_for_120/*.vcf.gz{,.csi}' \
    -resume \
    -profile cluster \