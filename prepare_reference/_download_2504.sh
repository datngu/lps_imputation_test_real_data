#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=download   
#SBATCH --mem=16G                
#SBATCH --partition=orion
#SBATCH --mail-type=ALL




download_vcf(){
    i=$1
    dir=$2
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz -O ${dir}/chr${i}_raw.vcf.gz

    bcftools view \
        --threads 4 \
        -m2 -M2 \
        -v snps \
        --exclude 'AC==1' ${dir}/chr${i}_raw.vcf.gz  \
        -Oz -o ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz

    ## view reference
    bcftools view \
        --threads 4 \
        -S ^unique_samples.txt ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz \
        -Oz -o ${dir}/ref_chr${i}.vcf.gz

    bcftools index -f ${dir}/ref_chr${i}.vcf.gz

    ## view eur gt
    bcftools view \
        --threads 4 \
        -S eur_samples.txt ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz \
        -Oz -o ${dir}/eur_chr${i}.vcf.gz

    bcftools index -f ${dir}/eur_chr${i}.vcf.gz

    ## view afr gt
    bcftools view \
        --threads 4 \
        -S afr_samples.txt ${dir}/chr${i}_2504_biallelic_snps_no_singleton.vcf.gz \
        -Oz -o ${dir}/afr_chr${i}.vcf.gz

    bcftools index -f ${dir}/afr_chr${i}.vcf.gz

}


## test chr22
filter_and_slipt_vcf 22 raw_data


## fix error of chr9
download_vcf 9 raw_data




# download_vcf(){
#     i=$1
#     dir=$2
#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz -O ${dir}/chr${i}_raw.vcf.gz

#     wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi -O ${dir}/chr${i}_raw.vcf.gz.tbi
# }

# mkdir raw_data
# for i in {1..11}
# do
#     download_vcf $i raw_data &
# done 

# wait

# for i in {12..22}
# do
#     download_vcf $i raw_data &
# done 

# wait
