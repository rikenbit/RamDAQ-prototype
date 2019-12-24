#!/usr/bin/env nextflow

proj_id = params.project_id
pipeline_class = params.pipeline_class
pipeline_species = params.pipeline_species

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

hisat2_options = Channel
    .from(params.hisat2_rrna_condition)
    .map{ [it[0], it[1]] }

hisat2_strandedness = Channel
    .from(params.hisat2_strandedness)
    .map{ [it[0], it[1]] }

hisat2_index = Channel
    .from(params.hisat2_rrna_index)
    .map{ [it[0], it[1], file(it[1]+"*")] }

hisat2_conditions = fastq_files_input
    .combine(hisat2_options)
    .combine(hisat2_strandedness)
    .combine(hisat2_index)
    .combine(pipeline_class)

//hisat2_conditions.println()

collect_rrnaQC_annot_ts = params.collect_rrnaQC_annot_ts
collect_rrnaQC_annot_refseq = params.collect_rrnaQC_annot_refseq

process run_hisat2  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/04_hisat2_rrna", mode: 'copy', overwrite: true

    clusterOptions = '-S /bin/bash -l nc=8'
    container "docker.io/myoshimura080822/hisat2_set:190611"

    input:
    val proj_id
    set run_id, fastq_name, file(fastq_L), file(fastq_R), option_name, option, strandedness_name, strandedness, index_name, index, file(index_files), pipeline_class from hisat2_conditions

    output:
    set run_id, pipeline_class, fastq_name, file("*.bam") into hisat2_output, hisat2_output_forsummary
    file "*.bai"
    file "*.command.err"
    file "*.bam" into hisat2_output_to_count

    script:
    def scripts_dir = workflow.scriptFile.parent.parent + "/bamtools_scripts"

    """
    hisat2 $option -x $index -1 $fastq_L -2 $fastq_R $strandedness 2> ${fastq_name}.hisat2.command.err | samtools view -bS - | samtools sort - -o ${fastq_name}_rrna_trim.sort.bam
    samtools index ${fastq_name}_rrna_trim.sort.bam 
    """
}

featurecounts_options = Channel
        .from(params.featurecounts_rrna_options)
        .map{ [it[0], it[1]] }

featurecounts_gtfs = Channel
        .from(params.featurecounts_rrna_gtfs)
        .map{ [it[0], file(it[1])] }

featurecounts_gtfs
    .into{
        featurecounts_gtfs_input
        featurecounts_gtfs_count
    }

featurecounts_conditions = hisat2_output
    .combine(featurecounts_options)
    .combine(featurecounts_gtfs_input)
    //.println()

process run_featurecounts  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/07_featurecounts_rrna/${gtf_name}/", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/hisat2_set:190611"

    input:
    val proj_id
    set run_id, pipeline_class, bam_name, file(bam_file), option_name, option, gtf_name, file(gtf) from featurecounts_conditions

    output:
    set run_id, gtf_name, pipeline_class, bam_name, file("fcounts_${bam_name}_rrna_trim.txt"), file("fcounts_${bam_name}_rrna_trim.txt.summary") into featurecounts_output
    file "fcounts_${bam_name}_rrna.log"
    file "*_rrna_trim.txt" into fcounts_output_to_count

    script:

    """
    featureCounts $option -g gene_id -a $gtf -o ./fcounts_${bam_name}_rrna_trim.txt $bam_file >& ./fcounts_${bam_name}_rrna.log
    """
}

process collect_hisat2_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"

    input:
    val proj_id
    set run_id, pipeline_class, bam_name, file(bam_file) from hisat2_output_forsummary.groupTuple()
    path summary_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_hisat2_summary.py"   

    output:
    file "*.txt"

    script:
    def summary_script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_hisat2_summary.py"

    """
    python $summary_script_path $PWD/output_${proj_id}/${run_id}/04_hisat2_rrna summary_hisat2_rrna_results PE
    """
}

process collect_featurecounts_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"

    input:
    val proj_id
    set run_id, gtf_name, pipeline_class, bam_name, file(fcounts_file), file(fcounts_summary_file) from featurecounts_output.groupTuple(by: [0,1])
    path collectcounts_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_counts.py"

    output:
    set run_id, file("*.txt") into collect_featurecounts_results_output

    script:

    """
    python $collectcounts_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts_rrna/${gtf_name}/ mergefcounts_${gtf_name}
    """
}

process collect_highsensitivity_rrnaQC_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/highsensitivity_rrnaQC_summary/", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/r-devel:2.4"

    input:
    val proj_id
    val collect_rrnaQC_annot_ts
    val collect_rrnaQC_annot_refseq
    set run_id, file(merge_fcounts_file) from collect_featurecounts_results_output.groupTuple()
    path collectcounts_script_mouse from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_highsensitivity_rrnaQC_summary_mouse.R"
    path collectcounts_script_human from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_highsensitivity_rrnaQC_summary_human.R"
    path collect_rrnaQC_annot_ts_file from collect_rrnaQC_annot_ts
    path collect_rrnaQC_annot_refseq_file from collect_rrnaQC_annot_refseq

    output:
    file "*.txt"

    script:

    if( pipeline_species == 'mouse' )
        """
        Rscript $collectcounts_script_mouse $PWD/output_${proj_id}/${run_id} ${collect_rrnaQC_annot_ts_file} ${collect_rrnaQC_annot_refseq_file}
        """
    else if( pipeline_species == 'human' )
        """
        Rscript $collectcounts_script_human $PWD/output_${proj_id}/${run_id} ${collect_rrnaQC_annot_refseq_file}
        """
}

// Check number of files with >0 byte file size
n_fastq = fastq_files_to_count.count {it[2].size() > 0}.getVal()
n_hisat2_output = hisat2_output_to_count.count {it.size() > 0}.getVal()
n_fcounts_output = fcounts_output_to_count.count {it.size() > 0}.getVal()
n_fcounts_gtf_length = featurecounts_gtfs_count.count {it[1]}.getVal()

println "======== Checking the number of files ========"

println "Number of gtf are $n_fcounts_gtf_length..."
n_fcounts_output = n_fcounts_output/n_fcounts_gtf_length

if(n_fastq == n_hisat2_output) {
    println "Number of files are same:-)"
    println "    fastq: $n_fastq"
    println "    hisat2: $n_hisat2_output"
} else{
    println "!!Caution!! Number of files are different:"

    println "    fastq: $n_fastq"
    println "    hisat2: $n_hisat2_output"
}

if(n_fastq == n_fcounts_output) {
    println "Number of files are same:-)"
    println "    fastq: $n_fastq"
    println "    featurecounts: $n_fcounts_output"
} else{
    println "!!Caution!! Number of files are different:"
    println "    fastq: $n_fastq"
    println "    featurecounts: $n_fcounts_output"
}

