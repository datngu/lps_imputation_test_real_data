#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=imp_sim067   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm


export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity



nextflow run bam2vcf.nf \
    -w 'work_bam2vcf_sim_0_67x' \
    --outdir sim_results \
    --trace_dir sim_trace_dir \
    --bam '/mnt/ScratchProjects/Aqua-Faang/dat_projects/lps_simulator/results/simulated_bam/*0.67_lps.bam{,.bai}' \
    --ref_vcf '/mnt/ScratchProjects/Aqua-Faang/dat_projects/igsr_data/ref_for_120x/*.vcf.gz{,.csi}' \
    -resume \
    -profile cluster \