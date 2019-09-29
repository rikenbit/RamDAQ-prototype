#!/usr/bin/env nextflow

proj_id = params.project_id

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

fastq_files = Channel
    .fromFilePairs("output_" + params.project_id + "/**/*_{R1,R2}_trim.fastq.gz")
    .map{[file(it[1][0]).parent.toString().replaceAll('/02_fastqmcf','').split('/')[file(it[1][0]).parent.toString().replaceAll('/02_fastqmcf','').split('/').length - 1], file(it[1][0]).baseName.replaceAll('.fastq',''), file(it[1][1]).baseName.replaceAll('.fastq',''), file(it[1][0]).baseName.replaceAll('_R1_trim.fastq',''), it[1][0], it[1][1]]}

fastq_files
    .into{
        fastq_files_input
        fastq_files_to_count
    }

sailfish_options = Channel
        .from(params.sailfish_condition)
        .map{ [it[0], it[1]] }

sailfish_index = Channel
        .from(params.sailfish_index)
        .map{ [it[0], file(it[1])] }

sailfish_conditions = fastq_files_input
    .combine(sailfish_options)
    .combine(sailfish_index)

//sailfish_conditions.println()

process run_sailfish  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/06_sailfish", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/sailfish:1.0"

    input:
    val proj_id
    set run_id, fastq_L_name, fastq_R_name, fastq_name, fastq_L, fastq_R, option_name, option, index_name, index from sailfish_conditions

    output:
    set run_id, fastq_name, file("sailfish_${fastq_name}_trim") into sailfish_output
    file "**/quant.sf" into sailfish_output_to_count

    script:

    """
    sailfish quant -i ${index} ${option} -1 <(gzip -dc ${fastq_L}) -2 <(gzip -dc ${fastq_R}) -o ./sailfish_${fastq_name}_trim
    """
}

process collect_sailfish_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"

    input:
    val proj_id
    set run_id, fastq_name, sailfish_outdir from sailfish_output.groupTuple()

    output:
    file "*.txt"

    script:
    def summary_script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_sailfish_summary.py"
    def collectcounts_script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_sailfish_counts.py"

    
    """
    python $summary_script_path $PWD/output_$proj_id/$run_id/06_sailfish summary_sailfish_results 
    python $collectcounts_script_path $PWD/output_$proj_id/$run_id/06_sailfish countdata_sailfish
    """
}

// Check number of files with >0 byte file size
n_fastq = fastq_files_to_count.count {it[2].size() > 0}.getVal()
n_sailfish_output = sailfish_output_to_count.count {it.size() > 0}.getVal()

println "======== Checking the number of files ========"
if(n_fastq == n_sailfish_output) {
    println "Number of files are same:-)"
    println "    fastq: $n_fastq"
    println "    sailfish: $n_sailfish_output"
} else{
    println "!!Caution!! Number of files are different:"
    println "    fastq: $n_fastq"
    println "    sailfish: $n_sailfish_output"
}


