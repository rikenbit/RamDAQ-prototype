#!/usr/bin/env nextflow

proj_id = params.project_id

pipeline_class = params.pipeline_class

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

bam_files = Channel
        .fromPath("output_" + params.project_id + "/**/04_hisat2/*.bam")
        .map { file -> tuple(file.parent.toString().replaceAll('/04_hisat2','').split('/')[file.parent.toString().replaceAll('/04_hisat2','').split('/').length - 1], file.baseName.replaceAll('_trim', ''), file.baseName.toString().replaceAll('_trim', '').split('\\.')[0], file.baseName.toString().replaceAll('_trim', '').split('\\.')[1], file) }

bam_files
    .into{
        bam_files_input
        bam_files_to_count
    }

bam_files_sortcount = Channel
        .fromPath("output_" + params.project_id + "/**/04_hisat2/*.sort.bam")
        .map { file -> tuple(file.parent.toString().replaceAll('/04_hisat2','').split('/')[file.parent.toString().replaceAll('/04_hisat2','').split('/').length - 1], file.baseName.replaceAll('_trim', ''), file.baseName.toString().replaceAll('_trim', '').split('\\.')[0], file) }

ref_bed = Channel
        .from(params.ref_beds)
        .map{ [it[0], file(it[1])] }

RSeQC_conditions = bam_files_input
    .combine(ref_bed)
    .combine(pipeline_class)

//RSeQC_conditions.println()

process run_RSeQC_readDist  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/read_distribution", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"
    
    input:
    val proj_id
    set run_id, bam_name, sample_name, strand_option, bam_file, bed_name, bed_file, pipeline_class from RSeQC_conditions

    output:
    set run_id, pipeline_class, bam_name, sample_name, strand_option, bam_file, bed_name, bed_file, file("*_readdist.txt") into readdist_output
    file "*_readdist.txt" into readDist_output_to_count

    script:

    """
    read_distribution.py -i $bam_file -r $bed_file > ${bam_name}_readdist.txt
    """
}

process run_RSeQC_geneBC  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/gene_bodycoverage", mode: 'copy', overwrite: true

    clusterOptions = '-S /bin/bash -l nc=1'
    
    container "docker.io/myoshimura080822/rseqc:1.2"   
 
    input:
    val proj_id
    set run_id, pipeline_class, bam_name, sample_name, strand_option, bam_file, bed_name, bed_file, readdist_file from readdist_output

    output:
    set run_id, pipeline_class, strand_option, bam_name, bam_file, bed_file, readdist_file, file("${bam_name}.geneBodyCoverage.txt") into genebc_output
    file "*.geneBodyCoverage.r"
    file "*.geneBodyCoverage.txt" into geneBC_output_to_count

    script:

    """
    geneBody_coverage.py -i $bam_file -r $bed_file -o ./${bam_name}
    """
}

process run_RSeQC_inferexp  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/infer_experiment", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set run_id, pipeline_class, strand_option, bam_name, bam_file, bed_file, readdist_file, geneBC_file from genebc_output

    output:
    set run_id, pipeline_class, strand_option, file("*.sort.inferexp.txt") into inferexp_output
    file "*.inferexp.txt" into inferexp_output_to_count

    when:
    bam_file =~ /.sort./

    script:

    """
    infer_experiment.py -i $bam_file -r $bed_file > ./${bam_name}.inferexp.txt
    """
}

process collect_RSeQC_summary {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/scientific_python2.7:1.0"

    input:
    val proj_id
    set run_id, pipeline_class, strand_option, infer_file from inferexp_output.groupTuple()
    
    output:
    file "*.txt"

    script:
    def readdist_script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_RSeQC_ReadDist_summary.py"
    def genebc_script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_RSeQC_geneBC_summary.py"
    def infer_script_path = workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_inferexperiment_summary.py"

    if( pipeline_class[0] == 'stranded' )
        """
         python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_PE_sort
         python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_PE_sort

         python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution forward summary_RSeQC_ReadDist_results_PE_forward
         python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage forward summary_RSeQC_geneBC_results_PE_forward

         python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution reverse summary_RSeQC_ReadDist_results_PE_reverse
         python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage reverse summary_RSeQC_geneBC_results_PE_reverse

         python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution R1 summary_RSeQC_ReadDist_results_PE_R1
         python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage R1 summary_RSeQC_geneBC_results_PE_R1

         python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution R2 summary_RSeQC_ReadDist_results_PE_R2
         python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage R2 summary_RSeQC_geneBC_results_PE_R2

         python $infer_script_path $PWD/output_$proj_id/$run_id/05_rseqc/infer_experiment summary_RSeQC_inferexperiment_results_PE_sort PE
        """
    else if( pipeline_class[0] == 'unstranded' )
        """
         python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_PE_sort
         python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_PE_sort
         python $infer_script_path $PWD/output_$proj_id/$run_id/05_rseqc/infer_experiment summary_RSeQC_inferexperiment_results_PE_sort PE
        """
}

// Check number of files with >0 byte file size
n_bam = bam_files_to_count.count {it[4].size() > 0}.getVal()
n_bam_sort = bam_files_sortcount.count {it[3].size() > 0}.getVal()
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


