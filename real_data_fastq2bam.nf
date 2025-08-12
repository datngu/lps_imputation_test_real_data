#!/usr/bin/env nextflow
/*
========================================================================================
                                FastQ to BAM
========================================================================================
 Author: Dat T Nguyen
 Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/



/*
 Define the default parameters
*/ 
params.metadata        = "$baseDir/meta_data/eur_rep1.csv"
params.bwa_idx         = "$baseDir/bwa_index/*"
params.genome_fn       = 'hg38.fa' // name of fasta file in bwa_idx

params.trace_dir       = "trace_dir"
params.outdir          = "results"

nextflow.enable.dsl=2

workflow {

    // Create a channel of samples from the metadata
    samplesChannel = Channel.fromPath(params.metadata) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sample_id, row.read1, row.read2)}

    bwa_idx_Ch = Channel.fromPath( params.bwa_idx, checkIfExists: true ).collect()    

    Download_fastq(samplesChannel)
    BWA_mapping(bwa_idx_Ch, params.genome_fn, Download_fastq.out)
    DedupBAM(BWA_mapping.out)

}

process Download_fastq {

    publishDir "${params.trace_dir}/fastq", mode: 'symlink', overwrite: true

    input:
    tuple val(sample_id), val(read1), val(read2)

    cpus 1
    memory '8GB'

    output:

    tuple val(sample_id), path ("${sample_id}_R1.fastq.gz"), path ("${sample_id}_R2.fastq.gz")

    script:
    """
    wget ${read1} -O ${sample_id}_R1.fastq.gz
    wget ${read2} -O ${sample_id}_R2.fastq.gz

    """
}



process BWA_mapping {

    publishDir "${params.trace_dir}/bwa_map", mode: 'symlink', overwrite: true

    input:
    path genome_idxs
    val genome_fn
    tuple val(sample_id), path(read1), path(read2)

    cpus 8
    memory '32GB'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem -M -t ${task.cpus} \
        ${genome_fn} \
        ${read1} \
        ${read2} | \
        samtools view -b -h -o ${sample_id}.bam

    """
}


process DedupBAM {

    publishDir "${params.outdir}/dedup_sorted_bam", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam_file)

    cpus 4
    memory '32GB'

    output:
    tuple path ("${sample_id}_sorted_dedup.bam"), path ("${sample_id}_sorted_dedup.bam.bai")

    script:
    """
    samtools sort -O BAM --threads ${task.cpus} -o sorted.bam $bam_file 
    samtools index sorted.bam

    java -Xmx32g -jar /biotools/picard.jar MarkDuplicates \\
        I=sorted.bam \\
        O=${sample_id}_sorted_dedup.bam \\
        METRICS_FILE=marked_dup_metrics.txt \\

    samtools index ${sample_id}_sorted_dedup.bam
    
    rm sorted.bam sorted.bam.bai

    """
}
