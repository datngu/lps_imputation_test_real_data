#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=bwa_real_60_eur_067x 
#SBATCH --mem=4G                
#SBATCH --mail-user=n.dat@outlook.com
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm


export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity


nextflow run real_data_fastq2bam.nf \
    -w 'work_fastq2bam' \
    --bwa_idx '/mnt/users/ngda/genomes/human_igsr/*' \
    -resume \
    -profile cluster \

