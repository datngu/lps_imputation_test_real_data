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

    Build_index(VCF_ref_ch, MAP_ch)

    BAM_all_ch = BAM_ch
        .map { sample_id, files -> files }  // keep only [bam, bai] list
        .flatten()                          // now each emitted item is a Path
        .collect()                          // collect into a single list
    
    BAM_ch.view()
    BAM_all_ch.view()

    Bam2Vcf(BAM_all_ch, Build_index.out)

}


process Bam2Vcf {

    publishDir "${params.outdir}/imputed_vcf", mode: 'copy', overwrite: true

    input:
    path all_bam
    tuple val(chr), path(idx_files)

    cpus 6
    memory '32GB'
    maxForks 30

    output:
    path "*imputed.vcf.gz"

    script:
    """
    ls *.bam > run_bam_list.txt


    run_imputation_bam_list.sh \
        run_bam_list.txt \
        ${chr} \
        ${chr}_imputed.vcf.gz \
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

