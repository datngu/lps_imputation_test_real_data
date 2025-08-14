#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=bwa_index   
#SBATCH --mem=16G                
#SBATCH --partition=orion
#SBATCH --mail-user=n.dat@outlook.com
#SBATCH --mail-type=ALL


img=/mnt/users/ngda/sofware/singularity/ndatth-ubuntu-22.04.img

## You can download the genome file from the 1000 Genomes Project
# wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O /mnt/users/ngda/genomes/human_igsr/hg38.fa

data_dir=/mnt/users/ngda/genomes/human_igsr
genome_file=hg38.fa

singularity exec --bind ${PWD}:/work_dir --bind ${data_dir}:/data \
    $img bwa index /data/${genome_file}


## singularity exec $img java -version