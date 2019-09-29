#!/usr/bin/env nextflow

proj_id = params.project_id

pipeline_class = params.pipeline_class

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

fastq_files = Channel
        .fromPath("output_" + params.project_id + "/**/*_trim.fastq.gz")
        .map { file -> tuple(file.parent.toString().replaceAll('/02_fastqmcf','').split('/')[file.parent.toString().replaceAll('/02_fastqmcf','').split('/').length - 1], file.baseName.replaceAll('_trim.fastq',''), file) }
        //.println()

fastq_files
    .into{
        fastq_files_input
        fastq_files_to_count
    }

hisat2_options = Channel
        .from(params.hisat2_condition)
        .map{ [it[0], it[1]] }

hisat2_strandedness = Channel
        .from(params.hisat2_strandedness)
        .map{ [it[0], it[1]] }

hisat2_index = Channel
        .from(params.hisat2_index)
        .map{ [it[0], file(it[1])] }

ref_chrsize = params.ref_chrsize

hisat2_conditions = fastq_files_input
    .combine(hisat2_options)
    .combine(hisat2_strandedness)
    .combine(hisat2_index)
    .combine(ref_chrsize)
    .combine(pipeline_class)

//hisat2_conditions.println()

process run_hisat2  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/04_hisat2", mode: 'copy', overwrite: true

    clusterOptions = '-S /bin/bash -l nc=8'
    container "docker.io/myoshimura080822/hisat2_set:190611"
    
    input:
    val proj_id
    set run_id, fastq_name, fastq, option_name, option, strandedness_name, strandedness, index_name, index, chrom_size, pipeline_class from hisat2_conditions

    output:
    set pipeline_class, run_id, fastq_name, file("*.bam"), chrom_size into hisat2_output
    file "*.bai"
    file "*.bam" into hisat2_output_to_count

    script:
    def scripts_dir = workflow.scriptFile.parent.parent + "/bamtools_scripts"

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

process run_bam2wig  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/04_hisat2_bigwig", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set pipeline_class, run_id, fastq_name, bam_files, chrom_size from hisat2_output

    output:
    file "*.bw" into bam2wig_output_to_count
    file "*.wig"

    script:
    
    if( pipeline_class == 'stranded' )
        """
        bam2wig.py -i ${bam_files[0]} -s $chrom_size -u -o ${bam_files[0].baseName}
        bam2wig.py -i ${bam_files[1]} -s $chrom_size -u -o ${bam_files[1].baseName}
        bam2wig.py -i ${bam_files[2]} -s $chrom_size -u -o ${bam_files[2].baseName}
        """
    else if( pipeline_class == 'unstranded' )
        """
        bam2wig.py -i $bam_files -s $chrom_size -u -o ${bam_files.baseName}
        """
}

// Check number of files with >0 byte file size
n_fastq = fastq_files_to_count.count {it[2].size() > 0}.getVal()
n_hisat2_output = hisat2_output_to_count.count {it.size() > 0}.getVal()
n_bam2wig_output = bam2wig_output_to_count.count {it.size() > 0}.getVal()

println "======== Checking the number of files ========"

if(n_fastq == n_hisat2_output) {
    println "Number of files are same:-)"
    println "    fastq: $n_fastq"
    println "    hisat2: $n_hisat2_output"
} else{
    println "!!Caution!! Number of files are different:"
    println "    fastq: $n_fastq"
    println "    hisat2: $n_hisat2_output"
}

if(n_fastq == n_bam2wig_output) {
    println "Number of files are same:-)"
    println "    fastq: $n_fastq"
    println "    bam2wig: $n_bam2wig_output"
} else{
    println "!!Caution!! Number of files are different:"
    println "    fastq: $n_fastq"
    println "    bam2wig: $n_bam2wig_output"
}











