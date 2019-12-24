#!/usr/bin/env nextflow

proj_id = params.project_id

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

bam_files = Channel
        .fromPath("output_" + params.project_id + "/**/04_hisat2/*_trim.sort.bam")
        .map { [file(file(it).parent.toString().replaceAll('/04_hisat2','')).name, it.baseName.replaceAll('_trim', '').replaceAll('.sort', ''), it]}

bam_files
    .into{
        bam_files_input
        bam_files_to_count
        bam_files_test
    }

//bam_files_test.println()
        
featurecounts_options = Channel
        .from(params.featurecounts_options)
        .map{ [it[0], it[1]] }

featurecounts_metafeature = Channel
        .from(params.featurecounts_metafeatures)
        .map{ [it[0], it[1]] }

featurecounts_metafeature
    .into{
        featurecounts_metafeature_input
        featurecounts_metafeature_count
    }

featurecounts_gtfs = Channel
        .from(params.featurecounts_gtfs)
        .map{ [it[0], file(it[1])] }

featurecounts_gtfs
    .into{
        featurecounts_gtfs_input
        featurecounts_gtfs_count
    }

featurecounts_conditions = bam_files_input
    .combine(featurecounts_options)
    .combine(featurecounts_metafeature_input)
    .combine(featurecounts_gtfs_input)
    //.println()

process run_featurecounts  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/hisat2_set:190611"

    input:
    val proj_id
    set run_id, bam_name, file(bam_file), option_name, option, metafeature_name, metafeature, gtf_name, file(gtf) from featurecounts_conditions

    output:
    set run_id, gtf_name, metafeature_name, bam_name, file("fcounts_${bam_name}_trim.txt"), file("fcounts_${bam_name}_trim.txt.summary") into featurecounts_output
    file "fcounts_${bam_name}_trim.log"
    file "*_trim.txt" into fcounts_output_to_count

    script:

    """
    featureCounts $option $metafeature -a $gtf -o ./fcounts_${bam_name}_trim.txt $bam_file >& ./fcounts_${bam_name}_trim.log
    """
}

process collect_featurecounts_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"

    input:
    val proj_id
    set run_id, gtf_name, metafeature_name, bam_name, file(fcounts_file), file(fcounts_summary_file) from featurecounts_output.groupTuple(by: [0,1,2])
    path summary_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_summary.py"
    path collectcounts_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_counts.py"

    output:
    file "*.txt"

    script:

    """
    python $summary_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name} summary_featurecounts_results_${gtf_name}_${metafeature_name} 
    python $collectcounts_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name} mergefcounts_${gtf_name}_${metafeature_name} 
    """
}

// Check number of files with >0 byte file size
n_bam = bam_files_to_count.count {it[2].size() > 0}.getVal()
n_fcounts_output = fcounts_output_to_count.count {it.size() > 0}.getVal()
n_fcounts_gtf_length = featurecounts_gtfs_count.count {it[1]}.getVal()
n_fcounts_metafeature_length = featurecounts_metafeature_count.count {it[1]}.getVal()

println "======== Checking the number of files ========"

println "Number of metafeature are $n_fcounts_metafeature_length..."
n_fcounts_output = n_fcounts_output/n_fcounts_metafeature_length

println "Number of gtf are $n_fcounts_gtf_length..."
n_fcounts_output = n_fcounts_output/n_fcounts_gtf_length

if(n_bam == n_fcounts_output) {
    println "Number of files are same:-)"
    println "    bam: $n_bam"
    println "    featurecounts: $n_fcounts_output"
} else{
    println "!!Caution!! Number of files are different:"
    println "    bam: $n_bam"
    println "    featurecounts: $n_fcounts_output"
}




