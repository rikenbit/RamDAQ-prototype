#!/usr/bin/env nextflow

/* 
 * nextflow run download-fastq.nf -c download-fastq.config
 */

Channel
    .fromPath(params.list_run_ids)
    .splitCsv(header: true, sep: "\t")
    .map{ 
        it.Run_ID
    }
    .into{run_ids; run_ids_}

nthreads = params.nthreads
repo = params.repo

n_fastq = run_ids_.count().getVal()
println "/*\nWe are now downloading FASTQ files for ${n_fastq} samples...\n*/"

process download_fastq {
    publishDir "output_download_fastq", mode: 'move'
    errorStrategy 'ignore'

    container "yuifu/cwltool-with-bash:1.0.20180809224403"

    input:
    val run_id from run_ids

    output:
    file "*.fastq.gz"

    script:
    """
    cwltool --debug "https://raw.githubusercontent.com/pitagora-network/pitagora-cwl/master/workflows/download-fastq/download-fastq.cwl" --run_ids ${run_id} --nthreads ${nthreads} --repo ${repo}
    """
}
