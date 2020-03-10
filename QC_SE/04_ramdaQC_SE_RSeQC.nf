#!/usr/bin/env nextflow

proj_id = params.project_id

pipeline_class = params.pipeline_class

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

bam_files = Channel
    .fromFilePairs("output_" + params.project_id + "/**/04_hisat2/*.{bam,bai}", flat: true) { file -> file.name.replaceAll(/.bam|.bai$/,'').replaceAll('_trim', '') }
    .map {[file(file(it[1]).parent.toString().replaceAll('/04_hisat2','')).name, it[0], file(it[1]), file(it[2])]}

bam_files
    .into{
        bam_files_input
        bam_files_to_count
        bam_tmp
    }

//bam_tmp.println()

bam_files_sortcount = Channel
        .fromPath("output_" + params.project_id + "/**/04_hisat2/*.sort.bam")

ref_bed = Channel
        .from(params.ref_beds)
        .map{ [it[0], file(it[1])] }

RSeQC_conditions = bam_files_input
    .combine(ref_bed)
    .combine(pipeline_class)

RSeQC_conditions.into{RSeQC_conditions1; RSeQC_conditions2; RSeQC_conditions3}

process run_RSeQC_readDist  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/read_distribution", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"
    
    input:
    val proj_id
    set run_id, bam_name, file(bam_file), file(bai_file), bed_name, file(bed_file), pipeline_class from RSeQC_conditions1

    output:
    set run_id, pipeline_class, file("*_readdist.txt") into readdist_output
    file "*_readdist.txt" into readDist_output_to_count

    script:
    """
    read_distribution.py -i $bam_file -r $bed_file > ${bam_name}_readdist.txt
    """
}

process run_RSeQC_geneBC  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/gene_bodycoverage", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/julia_genebodycoverage:1.2"
    
    input:
    val proj_id
    set run_id, bam_name, file(bam_file), file(bai_file), bed_name, file(bed_file), pipeline_class from RSeQC_conditions2

    output:
    set run_id, pipeline_class, file("${bam_name}.geneBodyCoverage.txt") into genebc_output
    file "*.geneBodyCoverage.txt" into geneBC_output_to_count

    script:

    """
    julia /opt/run.jl $bam_file $bed_file $bam_name
    """
}

process run_RSeQC_inferexp  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/infer_experiment", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set run_id, bam_name, file(bam_file), file(bai_file), bed_name, file(bed_file), pipeline_class from RSeQC_conditions3

    output:
    set run_id, pipeline_class, file("*.inferexp.txt") into inferexp_output
    file "*.inferexp.txt" into inferexp_output_to_count

    when:
    bam_file =~ /.sort./

    script:

    """
    infer_experiment.py -i $bam_file -r $bed_file > ./${bam_name}.inferexp.txt
    """
}



process collect_RSeQC_summary_readDist {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"
    
    input:
    val proj_id
    set run_id, pipeline_class, file(readdist_file) from readdist_output.groupTuple()
    path readdist_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_RSeQC_ReadDist_summary.py"

    output:
    file "*.txt"

    script:

    if( pipeline_class[0] == 'stranded' )

        """
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_SE_sort
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution forward summary_RSeQC_ReadDist_results_SE_forward
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution reverse summary_RSeQC_ReadDist_results_SE_reverse
        """

    else if( pipeline_class[0] == 'unstranded' )

        """
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution unstranded summary_RSeQC_ReadDist_results_SE_sort
        """

}


process collect_RSeQC_summary_geneBC {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"
    
    input:
    val proj_id
    set run_id, pipeline_class, file(genebc_file) from genebc_output.groupTuple()
    path genebc_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_RSeQC_geneBC_summary.py"

    output:
    file "*.txt"

    script:

    if( pipeline_class[0] == 'stranded' )

        """
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_SE_sort 
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage forward summary_RSeQC_geneBC_results_SE_forward
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage reverse summary_RSeQC_geneBC_results_SE_reverse
        """

    else if( pipeline_class[0] == 'unstranded' )

        """
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage unstranded summary_RSeQC_geneBC_results_SE_sort
        """

}



process collect_RSeQC_summary_inferexp {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"
    
    input:
    val proj_id
    set run_id, pipeline_class, file(infer_file) from inferexp_output.groupTuple()
    path infer_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_inferexperiment_summary.py"

    output:
    file "*.txt"

    script:

    if( pipeline_class[0] == 'stranded' )

        """
        python $infer_script_path $PWD/output_$proj_id/$run_id/05_rseqc/infer_experiment summary_RSeQC_inferexperiment_results_SE_sort SE
        """

    else if( pipeline_class[0] == 'unstranded' )

        """
        python $infer_script_path $PWD/output_$proj_id/$run_id/05_rseqc/infer_experiment summary_RSeQC_inferexperiment_results_SE_sort SE
        """

}




// Check number of files with >0 byte file size
n_bam = bam_files_to_count.count {it[3].size() > 0}.getVal()
n_bam_sort = bam_files_sortcount.count {it.size() > 0}.getVal()
n_readDist_output = readDist_output_to_count.count {it.size() > 0}.getVal()
n_geneBC_output = geneBC_output_to_count.count {it.size() > 0}.getVal()
n_inferexp_output = inferexp_output_to_count.count {it.size() > 0}.getVal()

println "======== Checking the number of files ========"

if(n_bam == n_readDist_output) {
    println "Number of files are same:-)"
    println "    bam: $n_bam"
    println "    readDist: $n_readDist_output"
} else {
    println "!!Caution!! Number of files are different:"
    println "    bam: $n_bam"
    println "    readDist: $n_readDist_output"
}

if(n_bam == n_geneBC_output) {
    println "Number of files are same:-)"
    println "    bam: $n_bam"
    println "    geneBC: $n_geneBC_output"
} else {
    println "!!Caution!! Number of files are different:"
    println "    bam: $n_bam"
    println "    geneBC: $n_geneBC_output"
}

if(n_bam_sort == n_inferexp_output) {
    println "Number of files are same:-)"
    println "    bam: $n_bam_sort"
    println "    inferexp: $n_inferexp_output"
} else {
    println "!!Caution!! Number of files are different:"
    println "    bam: $n_bam_sort"
    println "    inferexp: $n_inferexp_output"
}
