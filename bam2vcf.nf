#!/usr/bin/env nextflow
/*
========================================================================================
                            BAM to VCF via GLIMPSE2
========================================================================================
 Author: Dat T Nguyen
 Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/



/*
 Define the default parameters
*/ 
params.bam             = "$baseDir/data/*.bam{,.bai}"
params.ref_vcf         = "$baseDir/data/*.vcf{,.csi}"
params.maps            = "$baseDir/maps/*.gmap.gz"
// vcf in should be named by the convention: chr22.vcf.gz  chr22.vcf.gz.csi
// chr1, .., chr22 will be use directly as contig name, be aware of this

params.trace_dir       = "trace_dir"
params.outdir          = "results"

nextflow.enable.dsl=2

workflow {

    BAM_ch = Channel.fromFilePairs(params.bam)
    VCF_ref_ch = Channel.fromFilePairs(params.ref_vcf)
    MAP_ch = Channel.fromPath(params.maps).collect() 
    BAM_ch.view()

    Build_index(VCF_ref_ch, MAP_ch)

    Build_index.out.view()
    
    //Impute_input_ch = BAM_ch.combine(Build_index.out)

    //IMPUTE_glimpse2(Impute_input_ch)

    //IMPUTE_result_ch = IMPUTE_glimpse2.out.groupTuple().map{it -> [it[0], it[1].flatten()]}
    //IMPUTE_result_ch.view()
    //JOIN_vcfs(IMPUTE_result_ch)
    //Build_index.out.view()
}


process Bam2Vcf {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path all_bam
    path all_bam_idx
    path all_index
    tuple val(i), val(lps_cov)

    cpus 6
    memory '32GB'
    maxForks 30

    output:
    path "*imputed.vcf.gz"

    script:
    """
    ls *${lps_cov}.bam > run_bam_list.txt


    run_imputation_bam_list.sh \
        run_bam_list.txt \
        chr${i} \
        chr${i}_${lps_cov}_imputed.vcf.gz \
        ${task.cpus}
        
    """
}


process JOIN_vcfs {

    publishDir "${params.outdir}/imputed_joined_glimpse2", mode: 'copy', overwrite: true

    input:
    tuple val(chr), path(vcfs)

    cpus 1
    memory '32GB'

    output:
    tuple path("${chr}_imputed.vcf.gz"), path("${chr}_imputed.vcf.gz.csi")


    script:
    """
    bcftools merge *${chr}.vcf.gz -Oz -o "${chr}_imputed.vcf.gz"
    bcftools index "${chr}_imputed.vcf.gz"

    """
}


process IMPUTE_glimpse2 {

    publishDir "${params.trace_dir}/imputed_glimpse2", mode: 'symlink', overwrite: true

    input:
    tuple val(sample_id), path(bam), val(chr), path(idx_files)

    cpus 4
    memory '32GB'

    output:
    tuple val(chr), path("${sample_id}.${chr}.vcf.*")


    script:
    """
    run_imputation.sh \
        ${bam[0]} \
        ${chr} \
        ${sample_id}.${chr}.vcf.gz \
        ${task.cpus}

    """
}



process Build_index {

    publishDir "${params.trace_dir}/glimpse2_idxs", mode: 'symlink', overwrite: true

    input:
    tuple val(chr), path(vcf)
    path maps

    cpus 4
    memory '32GB'

    output:
    tuple val(chr), path("${chr}_idx/*")


    script:
    """
    buid_ref.sh ${vcf[0]} ${chr} ${chr}.b38.gmap.gz ${chr}_idx ${task.cpus}
    
    """
}

