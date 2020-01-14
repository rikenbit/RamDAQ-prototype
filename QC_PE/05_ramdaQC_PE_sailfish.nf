#!/usr/bin/env nextflow

proj_id = params.project_id

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

fastq_files = Channel
    .fromFilePairs("output_" + params.project_id + "/**/*_{R1,R2}_trim.fastq.gz", flat: true) { file -> file.name.replaceAll(/_R1|_R2/,''    ).replaceAll('_trim', '').replaceAll('.fastq.gz', '') }
    .map{[file(file(it[1]).parent.toString().replaceAll('/02_fastqmcf','')).name, it[0], it[1], it[2]]}

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
        .map{ [it[0], it[1], file(it[1]+"*")] }

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
    set run_id, fastq_name, file(fastq_L), file(fastq_R), option_name, option, index_name, index, file(index_files)  from sailfish_conditions

    output:
    set run_id, fastq_name, file("sailfish_${fastq_name}_trim") into sailfish_output

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
    path summary_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_sailfish_summary.py"
    path collectcounts_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_sailfish_counts.py"

    output:
    file "*.txt"

    script:
    
    """
    python $summary_script_path $PWD/output_$proj_id/$run_id/06_sailfish summary_sailfish_results 
    python $collectcounts_script_path $PWD/output_$proj_id/$run_id/06_sailfish countdata_sailfish
    """
}


