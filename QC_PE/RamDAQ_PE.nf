#!/usr/bin/env nextflow

proj_id = params.project_id

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

// checking files exist from config setting

chk_adapter = Channel
    .from(params.adapter)
    .map{file(it[1], checkIfExists: true)}
    .ifEmpty { exit 1, "adapter file not found" }

chk_index = Channel
    .from(params.hisat2_index)
    .flatMap{file(it[1]+"*", checkIfExists: true)}
    .ifEmpty { exit 1, "hisat2 index not found" }

chk_chrsize = Channel
    .from(params.ref_chrsize)
    .map{file(it[1], checkIfExists: true)}
    .ifEmpty { exit 1, "chromsize file not found" }

chk_bedfile = Channel
    .from(params.ref_beds)
    .map{file(it[1], checkIfExists: true)}
    .ifEmpty { exit 1, "bed file not found" }

chk_gtffile = Channel
    .from(params.featurecounts_gtfs)
    .map{file(it[1], checkIfExists: true)}
    .ifEmpty { exit 1, "gtf file not found" }

Channel
    .fromPath(params.fastq_filelist)
    .splitCsv(header: true, sep: "\t")
    .map{
        [file(it.Fastq1).parent.toString().split('/')[-2], it.Sample_ID, it.Sample_ID + "_R1", it.Sample_ID + "_R2", file(it.Fastq1), file(it.Fastq2), file(it.Fastq1).baseName.replaceAll('.fastq',''), file(it.Fastq2).baseName.replaceAll('.fastq','')]
        //run_id, fastq_name, fastq_L_name, fastq_R_name, fastq_L, fastq_R, fastq_L_basename, fastq_R_basename
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

//fastqmcf_conditions.println()

process run_fastQC_bef {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/03_fastQC_beforeTrimming", mode: 'copy', overwrite: true

    container "genomicpariscentre/fastqc:0.11.5"

    input:
    val proj_id
    set run_id, fastq_name, fastq_L_name, fastq_R_name, file(fastq_L), file(fastq_R), fastq_L_basename, fastq_R_basename, option_name, option, adapter_name, file(adapterfile) from fastqc_bef_conditions

    output:
    set run_id, file("${fastq_L_basename}_fastqc"), file("${fastq_R_basename}_fastqc") into fastqc_bef_output

    script:
    """
    fastqc -o . --nogroup $fastq_L
    unzip ${fastq_L_basename}_fastqc.zip
    fastqc -o . --nogroup $fastq_R
    unzip ${fastq_R_basename}_fastqc.zip
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
    set run_id, fastq_name, fastq_L_name, fastq_R_name, file(fastq_L), file(fastq_R), fastq_L_basename, fastq_R_basename, option_name, option, adapter_name, file(adapterfile) from fastqmcf_conditions

    output:
    set run_id, fastq_L_name, fastq_R_name, file("${fastq_L_name}_trim.fastq.gz"), file("${fastq_R_name}_trim.fastq.gz") into fastqmcf_output_tofastqc
    set run_id, fastq_name, file("${fastq_L_name}_trim.fastq.gz"), file("${fastq_R_name}_trim.fastq.gz") into fastqmcf_output_tohisat2

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
    set run_id, fastq_L_name, fastq_R_name, file(fastq_L_trim), file(fastq_R_trim) from fastqmcf_output_tofastqc

    output:
    set run_id, file("${fastq_L_name}_trim_fastqc"), file("${fastq_R_name}_trim_fastqc") into fastqc_output

    script:
    """
    fastqc -o . --nogroup $fastq_L_trim
    unzip ${fastq_L_name}_trim_fastqc.zip
    fastqc -o . --nogroup $fastq_R_trim
    unzip ${fastq_R_name}_trim_fastqc.zip
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
    set run_id, fastq_name, file(fastq_L), file(fastq_R), option_name, option, strandedness_name, strandedness, index_name, index, file(index_files), chrom_size, file(chrom_size_file), pipeline_class from hisat2_conditions
    path scripts_dir from workflow.scriptFile.parent.parent + "/bamtools_scripts"

    output:
    set chrom_size, pipeline_class, run_id, fastq_name, file("*.bam"), file("*.bai") into hisat2_output
    file "*.bam" into hisat2_output_to_count

    script:
    if( pipeline_class == 'stranded' )
        """
        hisat2 $option -x $index -1 $fastq_L -2 $fastq_R $strandedness -S ${fastq_name}_trim.sam
        samtools view -bS ${fastq_name}_trim.sam | samtools sort - -o ${fastq_name}_trim.sort.bam
        samtools index ${fastq_name}_trim.sort.bam
        rm ${fastq_name}_trim.sam

        bamtools filter -in ${fastq_name}_trim.sort.bam -out ${fastq_name}_trim.forward.bam -script ${scripts_dir}/bamtools_f_PE.json
        samtools index ${fastq_name}_trim.forward.bam

        bamtools filter -in ${fastq_name}_trim.sort.bam -out ${fastq_name}_trim.reverse.bam -script ${scripts_dir}/bamtools_r_PE.json
        samtools index ${fastq_name}_trim.reverse.bam

        samtools view -bS -f 0x40 ${fastq_name}_trim.sort.bam -o ${fastq_name}_trim.R1.bam
        samtools index ${fastq_name}_trim.R1.bam

        samtools view -bS -f 0x80 ${fastq_name}_trim.sort.bam -o ${fastq_name}_trim.R2.bam
        samtools index ${fastq_name}_trim.R2.bam
        """
    else if( pipeline_class == 'unstranded' )
        """
        hisat2 $option -x $index -1 $fastq_L -2 $fastq_R $strandedness -S ${fastq_name}_trim.sam
        samtools view -bS ${fastq_name}_trim.sam | samtools sort - -o ${fastq_name}_trim.sort.bam
        samtools index ${fastq_name}_trim.sort.bam
        rm ${fastq_name}_trim.sam
        """
}

hisat2_output
    .transpose()
    .into{
        hisat2_output_tobam2wig
        hisat2_output_tofcount
        hisat2_output_torseqc
    }

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

process run_featurecounts  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${option_name}/${metafeature_name}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/hisat2_set:190611"

    input:
    val proj_id
    set chrom_size, pipeline_class, run_id, bam_name, file(bam_file), file(bai_files), option_name, option, metafeature_name, metafeature, gtf_name, file(gtf) from featurecounts_conditions_fcount

    output:
    set run_id, gtf_name, option_name, metafeature_name, pipeline_class, bam_name, file("fcounts_${bam_name}_trim.txt"), file("fcounts_${bam_name}_trim.txt.summary") into featurecounts_output
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
    set run_id, gtf_name, option_name, metafeature_name, pipeline_class, bam_name, file(fcounts_file), file(fcounts_summary_file) from featurecounts_output.groupTuple(by: [0,1,2,3])
    path summary_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_summary.py"
    path collectcounts_script_path from workflow.scriptFile.parent.parent + "/collect_output_scripts/collect_featurecounts_counts.py"

    output:
    set run_id, file("summary_featurecounts_results_*.txt") into featurecounts_out01
    set run_id, file("mergefcounts_*.txt") into featurecounts_out02 

    script:

    """
    python $summary_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${option_name}/${metafeature_name} summary_featurecounts_results_${gtf_name}_${option_name}_${metafeature_name}
    python $collectcounts_script_path $PWD/output_${proj_id}/${run_id}/07_featurecounts/${gtf_name}/${option_name}/${metafeature_name} mergefcounts_${gtf_name}_${option_name}_${metafeature_name}
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
        RSeQC_conditions3
    }

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

process run_genebodycoverage  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/05_rseqc/gene_bodycoverage", mode: 'copy', overwrite: true

    container "docker.io/yuifu/julia_genebodycoverage:1.3.1-4"

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
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_PE_sort
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution forward summary_RSeQC_ReadDist_results_PE_forward
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution reverse summary_RSeQC_ReadDist_results_PE_reverse
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution R1 summary_RSeQC_ReadDist_results_PE_R1
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution R2 summary_RSeQC_ReadDist_results_PE_R2
        """
    else if( pipeline_class[0] == 'unstranded' )

        """
        python $readdist_script_path $PWD/output_$proj_id/$run_id/05_rseqc/read_distribution sort summary_RSeQC_ReadDist_results_PE_sort
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
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_PE_sort
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage forward summary_RSeQC_geneBC_results_PE_forward
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage reverse summary_RSeQC_geneBC_results_PE_reverse
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage R1 summary_RSeQC_geneBC_results_PE_R1
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage R2 summary_RSeQC_geneBC_results_PE_R2
        """
    else if( pipeline_class[0] == 'unstranded' )

        """
        python $genebc_script_path $PWD/output_$proj_id/$run_id/05_rseqc/gene_bodycoverage sort summary_RSeQC_geneBC_results_PE_sort
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
        python $infer_script_path $PWD/output_$proj_id/$run_id/05_rseqc/infer_experiment summary_RSeQC_inferexperiment_results_PE_sort PE
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

pipeline_class_forreport = params.pipeline_class

process execute_nbconvert {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/jupyternotebook:2.4"

    input:
    val proj_id
    val pipeline_class_forreport
    set run_id, file('*'), file('*'), file('*'), file('*'), file('*') from nbconvert_input
    path proj_dir from workflow.workDir.parent + "/output_${proj_id}"
    path notebook_path_unstranded from workflow.scriptFile.parent.parent + "/R_QCplot/RamDA-SeqQC_template_PE_unstranded_nbconvert.ipynb"
    path notebook_path_stranded from workflow.scriptFile.parent.parent + "/R_QCplot/RamDA-SeqQC_template_PE_stranded_nbconvert.ipynb"
    path function_file from workflow.scriptFile.parent.parent + "/R_QCplot/00_sampleQC_function_nbconvert.R"

    output:
    file "*.html"
    file "*.ipynb"
    file "*TPM.txt"

    script:
    if( pipeline_class_forreport[0] == 'stranded' )
        """
        jupyter nbconvert --to html --execute $notebook_path_stranded --output ${run_id}_notebook_PE_stranded.html --ExecutePreprocessor.timeout=2678400 --allow-errors --debug

        jupyter nbconvert --to notebook --execute $notebook_path_stranded --output ${run_id}_notebook_PE_stranded.ipynb --ExecutePreprocessor.timeout=2678400 --allow-errors --debug
        """
    else if( pipeline_class_forreport[0] == 'unstranded' )
        """
        jupyter nbconvert --to html --execute $notebook_path_unstranded --output ${run_id}_notebook_PE_unstranded.html --ExecutePreprocessor.timeout=2678400 --allow-errors --debug

        jupyter nbconvert --to notebook --execute $notebook_path_unstranded --output ${run_id}_notebook_PE_unstranded.ipynb --ExecutePreprocessor.timeout=2678400 --allow-errors --debug
        """
}









