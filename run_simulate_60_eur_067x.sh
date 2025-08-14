#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=simulae_60_eur_067x 
#SBATCH --mem=4G                
#SBATCH --mail-user=n.dat@outlook.com
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
#module load Nextflow/21.03
module load nextflow/20.10.0
module load singularity/rpm


export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity


nextflow run simulate_data_bam.nf \
    -w 'bam_60_eur_067x' \
    --metadata 'meta_data/eur_rep1_high_coverage_info.csv' \
    --genome '/mnt/users/ngda/genomes/human_igsr/hg38.fa' \
    -resume \
    -profile cluster \

