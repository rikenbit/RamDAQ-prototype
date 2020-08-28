#!/usr/bin/env nextflow

proj_id = params.project_id
bso_dir = params.bso_dir
run_ids = Channel.from(params.run_ids)
          .map{[it[0]]}

process run_bcl2fastq  {

    clusterOptions = '-S /bin/bash -l nc=8'
    
    tag { "${run_id}"}
    publishDir "output_${proj_id}/${run_id}/01_fastq_files"

    container "docker.io/myoshimura080822/bcl2fastq2:1.0"

    input:
    val proj_id
    val bso_dir
    val run_id from run_ids

    script:
    def runfolder_dir = "${bso_dir}/${run_id}"
    def output_dir = "output_${proj_id}/${run_id}/01_fastq_files" 

    """
    bcl2fastq --no-lane-splitting --runfolder-dir $runfolder_dir --interop-dir $runfolder_dir/InterOp --input-dir $runfolder_dir/Data/Intensities/BaseCalls --sample-sheet $runfolder_dir/SampleSheet.csv --output-dir $PWD/$output_dir --stats-dir $PWD/$output_dir/Stats --reports-dir $PWD/$output_dir/Reports && rm -rf $PWD/$output_dir/Undetermined*
    """
}

