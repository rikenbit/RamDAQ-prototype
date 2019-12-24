#!/usr/bin/env nextflow

proj_id = params.project_id

fastq_filelist = params.fastq_filelist

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

Channel
    .fromPath(fastq_filelist)
    .splitCsv(header: true, sep: "\t")
    .map{ 
        [file(it.Fastq1).parent.toString().split('/')[-2], it.Sample_ID + "_R1", it.Sample_ID + "_R2", file(it.Fastq1), file(it.Fastq2), file(it.Fastq1).baseName.replaceAll('.fastq',''), file(it.Fastq2).baseName.replaceAll('.fastq','')]
        //run_id, fastq_L_name, fastq_R_name, fastq_L, fastq_R, fastq_L_basename, fastq_R_basename
    }
    .into{fastq_files; fastq_files_tmp; }
        
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

//fastq_files_tmp.println()

process run_fastQC_bef {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/03_fastQC_beforeTrimming", mode: 'copy', overwrite: true

    container "genomicpariscentre/fastqc:0.11.5"

    input:
    val proj_id
    set run_id, fastq_L_name, fastq_R_name, file(fastq_L), file(fastq_R), fastq_L_basename, fastq_R_basename, option_name, option, adapter_name, file(adapterfile) from fastqc_bef_conditions

    output:
    set run_id, file("${fastq_L_basename}_fastqc"), file("${fastq_R_basename}_fastqc") into fastqc_bef_output

    script:
    """
    fastqc -o . --nogroup $fastq_L && unzip ${fastq_L_basename}_fastqc.zip
    fastqc -o . --nogroup $fastq_R && unzip ${fastq_R_basename}_fastqc.zip
    """
}

process collect_fastQC_bef_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"
    input:
    val proj_id
    set run_id, fastq_L_dir, fastq_R_dir from fastqc_bef_output.groupTuple()
    path script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_fastqc_summary.py"

    output:
    file "*.txt"

    script:
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
    set run_id, fastq_L_name, fastq_R_name, file(fastq_L), file(fastq_R), fastq_L_basename, fastq_R_basename, option_name, option, adapter_name, file(adapterfile) from fastqmcf_conditions

    output:
    set run_id, fastq_L_name, fastq_R_name, file("${fastq_L_name}_trim.fastq.gz"), file("${fastq_R_name}_trim.fastq.gz") into fastqmcf_output

    script:
    """
    fastq-mcf $adapterfile $fastq_L $fastq_R -o ${fastq_L_name}_trim.fastq -o ${fastq_R_name}_trim.fastq $option && gzip ${fastq_L_name}_trim.fastq && gzip ${fastq_R_name}_trim.fastq
    """

}


process run_fastQC {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/03_fastQC", mode: 'copy', overwrite: true
 
    container "genomicpariscentre/fastqc:0.11.5"

    input:
    val proj_id
    set run_id, fastq_L_name, fastq_R_name, file(fastq_L_trim), file(fastq_R_trim) from fastqmcf_output
 
    output:
    set run_id, file("${fastq_L_name}_trim_fastqc"), file("${fastq_R_name}_trim_fastqc") into fastqc_output
 
    script:
    """
    fastqc -o . --nogroup $fastq_L_trim && unzip ${fastq_L_name}_trim_fastqc.zip
    fastqc -o . --nogroup $fastq_R_trim && unzip ${fastq_R_name}_trim_fastqc.zip
    """
}

process collect_fastQC_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"

    input:
    val proj_id
    set run_id, fastq_L_dir, fastq_R_dir from fastqc_output.groupTuple()
    path script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_fastqc_summary.py"    

    output:
    file "*.txt"

    script:
    """
    python $script_path $PWD/output_$proj_id/$run_id/03_fastQC summary_fastQC_result
    """

}

