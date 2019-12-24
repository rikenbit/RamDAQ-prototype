#!/usr/bin/env nextflow

proj_id = params.project_id

pipeline_class = params.pipeline_class

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

fastq_files = Channel
    .fromFilePairs("output_" + params.project_id + "/**/*_{R1,R2}_trim.fastq.gz", flat: true) { file -> file.name.replaceAll(/_R1|_R2/,'').replaceAll('_trim', '').replaceAll('.fastq.gz', '') }
    .map{[file(file(it[1]).parent.toString().replaceAll('/02_fastqmcf','')).name, it[0], it[1], it[2]]}

fastq_files
    .into{
        fastq_files_test
        fastq_files_input
        fastq_files_to_count
    }

//fastq_files_test.println()

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
    set run_id, fastq_name, file(fastq_L), file(fastq_R), option_name, option, strandedness_name, strandedness, index_name, index, file(index_files), chrom_size, file(chrom_size_file), pipeline_class from hisat2_conditions
    path scripts_dir from workflow.scriptFile.parent.parent + "/bamtools_scripts"

    output:
    set pipeline_class, run_id, fastq_name, file("*.bam"), file("*.bai"), file(chrom_size_file) into hisat2_output
    file "*.bai"
    file "*.bam" into hisat2_output_to_count

    script:

    if( pipeline_class == 'stranded' )
        """
        hisat2 $option -x $index -1 $fastq_L -2 $fastq_R $strandedness | samtools view -bS - | samtools sort - -o ${fastq_name}_trim.sort.bam
        samtools index ${fastq_name}_trim.sort.bam

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
        hisat2 $option -x $index -1 $fastq_L -2 $fastq_R $strandedness | samtools view -bS - | samtools sort - -o ${fastq_name}_trim.sort.bam
        samtools index ${fastq_name}_trim.sort.bam 
        """
}


process run_bam2wig  {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}/04_hisat2_bigwig", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/rseqc:1.2"

    input:
    val proj_id
    set pipeline_class, run_id, fastq_name, file(bam_files), file(bai_files), file(chrom_size) from hisat2_output

    output:
    file "*.bw" into bam2wig_output_to_count
    file "*.wig"

    script:

    if( pipeline_class == 'stranded' )
        """
        bam2wig.py -i ${bam_files[2]} -s $chrom_size -u -o ${bam_files[2].baseName}
        bam2wig.py -i ${bam_files[3]} -s $chrom_size -u -o ${bam_files[3].baseName}
        bam2wig.py -i ${bam_files[4]} -s $chrom_size -u -o ${bam_files[4].baseName}
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



