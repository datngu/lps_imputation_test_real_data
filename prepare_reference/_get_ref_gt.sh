#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --nodes=1                
#SBATCH --job-name=slipt_ref   
#SBATCH --mem=16G                
#SBATCH --partition=orion
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


filter_and_slipt_vcf(){
    i=$1
    dir=$2
    ref_dir=$3
    gt_dir=$4
    # wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz -O ${dir}/chr${i}_raw.vcf.gz

    bcftools view \
        --threads 4 \
        -m2 -M2 \
        -v snps \
        --exclude 'AC==1' \
        -S 2504_samples.txt ${dir}/chr${i}_raw.vcf.gz  \
        -Oz -o ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz

    ## view reference
    bcftools view \
        --threads 4 \
        -S ^unique_samples.txt ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz \
        -Oz -o ${ref_dir}/chr${i}.vcf.gz

    bcftools index -f ${ref_dir}/chr${i}.vcf.gz

    ## view eur gt
    bcftools view \
        --threads 4 \
        -S eur_samples.txt ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz \
        -Oz -o ${gt_dir}/eur_chr${i}.vcf.gz

    bcftools index -f ${gt_dir}/eur_chr${i}.vcf.gz

    ## view afr gt
    bcftools view \
        --threads 4 \
        -S afr_samples.txt ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz \
        -Oz -o ${gt_dir}/afr_chr${i}.vcf.gz

    bcftools index -f ${gt_dir}/afr_chr${i}.vcf.gz

}

mkdir -p ref_for_120
mkdir -p ground_true_120
## test chr22
for i in {22..1}
do
    filter_and_slipt_vcf $i raw_data ref_for_120 ground_true_120
done


