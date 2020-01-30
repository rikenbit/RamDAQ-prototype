#!/usr/bin/env nextflow

proj_id = params.project_id

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
    set run_id, fastq_name, file(fastq), option_name, option, adapter_name, file(adapterfile) from fastqmcf_conditions

    output:
    set run_id, fastq_name, file("${fastq_name}_trim.fastq.gz") into fastqmcf_output

    script:
    def output_dir = "output_${proj_id}/${run_id}/02_fastqmcf"

    """
    fastq-mcf $adapterfile $fastq -o ${fastq_name}_trim.fastq $option ; gzip ${fastq_name}_trim.fastq
    """
}

fastqmcf_output
    .into{
        fastqmcf_output_tofastqc
        fastqmcf_output_tohisat2
    }

process run_fastQC {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/03_fastQC", mode: 'copy', overwrite: true
    
    container "genomicpariscentre/fastqc:0.11.5"

    input:
    val proj_id
    set run_id, fastq_name, file(fastq_file) from fastqmcf_output_tofastqc
 
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
    set run_id, file("*.txt") into fastqc_results

    script:
    """
    python $script_path $PWD/output_$proj_id/$run_id/03_fastQC summary_fastQC_result
    """
}

pipeline_class = params.pipeline_class

hisat2_options = Channel
        .from(params.hisat2_condition)
        .map{ [it[0], it[1]] }

hisat2_strandedness = Channel
        .from(params.hisat2_strandedness)
        .map{ [it[0], it[1]] }

hisat2_index = Channel
        .from(params.hisat2_index)
        .map{ [it[0], it[1], file(it[1]+"*")] }

ref_chrsize = Channel
        .from(params.ref_chrsize)
        .map{ [it[0], file(it[1])] }

hisat2_conditions = fastqmcf_output_tohisat2
    .combine(hisat2_options)
    .combine(hisat2_strandedness)
    .combine(hisat2_index)
    .combine(ref_chrsize)
    .combine(pipeline_class)

process run_hisat2 {
    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/04_hisat2", mode: 'copy', overwrite: true

    clusterOptions = '-S /bin/bash -l nc=8'
    container "docker.io/myoshimura080822/hisat2_set:190611"

    input:
    val proj_id
    set run_id, fastq_name, file(fastq), option_name, option, strandedness_name, strandedness, index_name, index, file(index_files), chrom_size, file(chrom_size_file), pipeline_class from hisat2_conditions
    path scripts_dir from workflow.scriptFile.parent.parent + "/bamtools_scripts"

    output:
    set chrom_size, pipeline_class, run_id, fastq_name, file("*.bam"), file("*.bai") into hisat2_output
    file "*.bam" into hisat2_output_to_count

    script:

    if( pipeline_class == 'stranded' )
        """
        hisat2 $option -x $index -U $fastq $strandedness | samtools view -bS - | samtools sort - -o ${fastq_name}_trim.sort.bam
        samtools index ${fastq_name}_trim.sort.bam

        bamtools filter -in ${fastq_name}_trim.sort.bam -out ${fastq_name}_trim.forward.bam -script ${scripts_dir}/bamtools_f_SE.json
        samtools index ${fastq_name}_trim.forward.bam

        bamtools filter -in ${fastq_name}_trim.sort.bam -out ${fastq_name}_trim.reverse.bam -script ${scripts_dir}/bamtools_r_SE.json
        samtools index ${fastq_name}_trim.reverse.bam
        """
    else if( pipeline_class == 'unstranded' )
        """
        hisat2 $option -x $index -U $fastq $strandedness | samtools view -bS - | samtools sort - -o ${fastq_name}_trim.sort.bam
        samtools index ${fastq_name}_trim.sort.bam
        """
}

hisat2_output
    .transpose()
    .into{
        hisat2_output_tobam2wig
        hisat2_output_tofcount
        hisat2_output_torseqc
        hisat2_output_print
    }

//hisat2_output_print.println()

ref_chrsize2 = Channel
    .from(params.ref_chrsize)
    .map{ [it[0], file(it[1])] }

bam2wig_input = hisat2_output_tobam2wig
    .combine(ref_chrsize2, by: 0)

process run_bam2wig {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/04_hisat2_bigwig", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set chrom_size, pipeline_class, run_id, fastq_name, file(bam_files), file(bai_files), file(chrom_size_file) from bam2wig_input

    output:
    file "*.bw" into bam2wig_output_to_count
    file "*.wig"

    script:

    """
    bam2wig.py -i $bam_files -s $chrom_size_file -u -o ${bam_files.baseName}
    """
}

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

featurecounts_conditions = hisat2_output_tofcount
    .combine(featurecounts_options)
    .combine(featurecounts_metafeature_input)
    .combine(featurecounts_gtfs_input)

featurecounts_conditions
    .into {
        featurecounts_conditions_fcount
        featurecounts_conditions_print
    }

//featurecounts_conditions_print.println()

process run_featurecounts  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/hisat2_set:190611"

    input:
    val proj_id
    set chrom_size, pipeline_class, run_id, bam_name, file(bam_file), file(bai_files), option_name, option, metafeature_name, metafeature, gtf_name, file(gtf) from featurecounts_conditions_fcount

    output:
    set run_id, gtf_name, metafeature_name, pipeline_class, bam_name, file("fcounts_${bam_name}_trim.txt"), file("fcounts_${bam_name}_trim.txt.summary") into featurecounts_output
    file "fcounts_${bam_name}_trim.log"
    file "*_trim.txt" into fcounts_output_to_count
    
    when:
    bam_file =~ /.sort./

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
    set run_id, gtf_name, metafeature_name, pipeline_class, bam_name, file(fcounts_file), file(fcounts_summary_file) from featurecounts_output.groupTuple(by: [0,1,2])
    path summary_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_summary.py"
    path collectcounts_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_counts.py"

    output:
    set run_id, file("summary_featurecounts_results_*.txt") into featurecounts_out01
    set run_id, file("mergefcounts_*.txt") into featurecounts_out02

    script:

    if( pipeline_class[0] == 'stranded' )
        """
        python $summary_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name} summary_featurecounts_results_${gtf_name}_${metafeature_name} stranded
        python $collectcounts_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name} mergefcounts_${gtf_name}_${metafeature_name}
        """
    else if( pipeline_class[0] == 'unstranded' )
        """
        python $summary_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name} summary_featurecounts_results_${gtf_name}_${metafeature_name} unstranded
        python $collectcounts_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${metafeature_name} mergefcounts_${gtf_name}_${metafeature_name}
        """

}

ref_bed = Channel
    .from(params.ref_beds)
    .map{ [it[0], file(it[1])] }

RSeQC_conditions = hisat2_output_torseqc
    .combine(ref_bed)

RSeQC_conditions
    .into{
        RSeQC_conditions1; 
        RSeQC_conditions2; 
        RSeQC_conditions3;
        RSeQC_conditions_print
    }

//RSeQC_conditions_print.println()

process run_RSeQC_readDist  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/read_distribution", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set chrom_size, pipeline_class, run_id, bam_name, file(bam_file), file(bai_files), bed_name, file(bed_file) from RSeQC_conditions1

    output:
    set run_id, pipeline_class, file("*_readdist.txt") into readdist_output
    file "*_readdist.txt" into readDist_output_to_count
    
    script:
    fileName = bam_file.baseName
    """
    read_distribution.py -i $bam_file -r $bed_file > ${fileName}_readdist.txt
    """
}

process run_RSeQC_geneBC  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/gene_bodycoverage", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/julia_genebodycoverage:1.1"

    input:
    val proj_id
    set chrom_size, pipeline_class, run_id, bam_name, file(bam_file), file(bai_files), bed_name, file(bed_file) from RSeQC_conditions2

    output:
    set run_id, pipeline_class, file("*.geneBodyCoverage.txt") into genebc_output
    file "*.geneBodyCoverage.txt" into geneBC_output_to_count
    
    script:
    fileName = bam_file.baseName
    """
    julia /opt/run.jl $bam_file $bed_file ${fileName}
    """
 
}

process run_RSeQC_inferexp  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/infer_experiment", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set chrom_size, pipeline_class, run_id, bam_name, file(bam_file), file(bai_files), bed_name, file(bed_file) from RSeQC_conditions3

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
    set run_id, file("*.txt") into rseqc_readdist_results

    script:

    if( pipeline_class[0] == 'stranded' )

        """
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_SE_sort
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution forward summary_RSeQC_ReadDist_results_SE_forward
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution reverse summary_RSeQC_ReadDist_results_SE_reverse
        """

    else if( pipeline_class[0] == 'unstranded' )

        """
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_SE_sort
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
    set run_id, file("*.txt") into rseqc_genebc_results

    script:

    if( pipeline_class[0] == 'stranded' )

        """
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_SE_sort
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage forward summary_RSeQC_geneBC_results_SE_forward
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage reverse summary_RSeQC_geneBC_results_SE_reverse
        """

    else if( pipeline_class[0] == 'unstranded' )

        """
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_SE_sort
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
    set run_id, file("*.txt") into rseqc_inferexp_results

    script:
        """
        python $infer_script_path $PWD/output_$proj_id/$run_id/05_rseqc/infer_experiment summary_RSeQC_inferexperiment_results_SE_sort SE
        """
}

fastqc_results
    .join(rseqc_readdist_results)
    .join(rseqc_genebc_results)
    .join(rseqc_inferexp_results)
    .into{summary_fastqcrseqc; summary_fastqcrseqc_join}

featurecounts_out01.mix(featurecounts_out02)
    .groupTuple(by: 0)
    .into{summary_fcounts; summary_fcounts_print}

summary_fastqcrseqc_join
    .join(summary_fcounts)
    .into{nbconvert_input; nbconvert_input_print}

//nbconvert_input_print.println()

process execute_nbconvert {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/jupyternotebook:2.4"

    input:
    val proj_id
    set run_id, file('*'), file('*'), file('*'), file('*'), file('*') from nbconvert_input
    path proj_dir from workflow.workDir.parent + "/output_${proj_id}"
    path notebook_path_unstranded from workflow.scriptFile.parent.parent + "/R_QCplot/RamDA-SeqQC_template_SE_unstranded_nbconvert.ipynb"
    path function_file from workflow.scriptFile.parent.parent + "/R_QCplot/00_sampleQC_function_nbconvert.R"

    output:
    file "*.html"
    file "*.ipynb"

    script:
    """
    jupyter nbconvert --to html --execute $notebook_path_unstranded --output ${run_id}_notebook_SE_unstranded.html --ExecutePreprocessor.timeout=2678400 --allow-errors --debug

    jupyter nbconvert --to notebook --execute $notebook_path_unstranded --output ${run_id}_notebook_SE_unstranded.ipynb --ExecutePreprocessor.timeout=2678400 --allow-errors --debug
    """
}









