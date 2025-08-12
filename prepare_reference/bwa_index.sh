#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=bwa_index   
#SBATCH --mem=16G                
#SBATCH --partition=orion
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


img=/mnt/users/ngda/sofware/singularity/ndatth-ubuntu-22.04.img
data_dir=/mnt/users/ngda/genomes/human_igsr
genome_file=hg38.fa

singularity exec --bind ${PWD}:/work_dir --bind ${data_dir}:/data \
    $img bwa index /data/${genome_file}


singularity exec $img java -version