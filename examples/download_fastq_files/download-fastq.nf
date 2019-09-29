#!/usr/bin/env nextflow

/* 
 * nextflow run download-fastq.nf -c download-fastq.config
 */

ymls = Channel.fromPath(params.path_yml)


process download_fastq {
    publishDir "output_download_fastq", mode: 'move'

    container "yuifu/cwltool-with-bash:1.0.20180809224403"

    input:
    file yml from ymls

    output:
    file "*.fastq.gz"

    script:
    """
    cwltool --debug "https://raw.githubusercontent.com/pitagora-network/pitagora-cwl/master/workflows/download-fastq/download-fastq.cwl" $yml
    """
}

