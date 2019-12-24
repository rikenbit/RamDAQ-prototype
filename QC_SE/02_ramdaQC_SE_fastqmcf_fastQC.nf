#!/usr/bin/env nextflow

proj_id = params.project_id

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

fastq_files = Channel
        .fromPath("output_" + params.project_id + "/**/01_fastq_files/*.fastq.gz")
        .map { [file(file(it).parent.toString().replaceAll('/01_fastq_files','')).name, it.baseName.replaceAll('.fastq', ''), it]}

fastqmcf_options = Channel
        .from(params.fastqmcf_condition)
        .map{ [it[0], it[1]] }

adapter = Channel
        .from(params.adapter)
        .map{ [it[0], file(it[1])] }

fastq_files
    .combine(fastqmcf_options)
    .combine(adapter)
    .into{fastqc_bef_conditions; fastqmcf_conditions;}

//fastqmcf_conditions.println()

process run_fastQC_bef {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/03_fastQC_beforeTrimming", mode: 'copy', overwrite: true

    container "genomicpariscentre/fastqc:0.11.5"

    input:
    val proj_id
    set run_id, fastq_name, file(fastq), option_name, option, adapter_name, file(adapterfile) from fastqc_bef_conditions

    output:
    set run_id, file("${fastq_name}_fastqc") into fastqc_bef_output

    script:
    """
    fastqc -o . --nogroup $fastq && unzip ${fastq_name}_fastqc.zip
    """
}

process collect_fastQC_bef_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"
    input:
    val proj_id
    set run_id, fastq_dir from fastqc_bef_output.groupTuple()

    output:
    file "*.txt"

    script:
    def script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_fastqc_summary.py"

    """
    python $script_path $PWD/output_$proj_id/$run_id/03_fastQC_beforeTrimming summary_beforetrim_fastQC_result
    """

}


process run_fastqmcf  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/02_fastqmcf", mode: 'copy', overwrite: true
   
    container "docker.io/myoshimura080822/fastqmcf:1.0"
  
    input:
    val proj_id
    set run_id, fastq_name, file(fastq), option_name, option, adapter_name, file(adapterfile) from fastqmcf_conditions

    output:
    set run_id, fastq_name, file("${fastq_name}_trim.fastq.gz") into fastqmcf_output

    script:
    def output_dir = "output_${proj_id}/${run_id}/02_fastqmcf"

    """
    fastq-mcf $adapterfile $fastq -o ${fastq_name}_trim.fastq $option ; gzip ${fastq_name}_trim.fastq
    """
}

process run_fastQC {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/03_fastQC", mode: 'copy', overwrite: true
    
    container "genomicpariscentre/fastqc:0.11.5"

    input:
    val proj_id
    set run_id, fastq_name, file(fastq_file) from fastqmcf_output
 
    output:
    set run_id, file("${fastq_name}_trim_fastqc") into fastqc_output
 
    script:
    """
    fastqc -o . --nogroup $fastq_file && unzip ${fastq_name}_trim_fastqc.zip
    """
}

process collect_fastQC_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"
    input:
    val proj_id
    set run_id, fastq_dir from fastqc_output.groupTuple()
    path script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_fastqc_summary.py"    

    output:
    file "*.txt"

    script:
    """
    python $script_path $PWD/output_$proj_id/$run_id/03_fastQC summary_fastQC_result
    """

}
