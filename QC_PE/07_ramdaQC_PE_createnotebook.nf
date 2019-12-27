#!/usr/bin/env nextflow

proj_id = params.project_id
run_ids = Channel.from(params.run_ids)
          .map{[it[0]]}

process execute_nbconvert {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    container "docker.io/myoshimura080822/jupyternotebook:2.4"

    input:
    val proj_id
    val run_id from run_ids
    path proj_dir from workflow.workDir.parent + "/output_${proj_id}"
    path notebook_path_unstranded from workflow.scriptFile.parent.parent + "/R_QCplot/RamDA-SeqQC_template_PE_unstranded_nbconvert.ipynb"
    path function_file from workflow.scriptFile.parent.parent + "/R_QCplot/00_sampleQC_function_nbconvert.R"

    output:
    file "*.html"
    file "*.ipynb"

    script:
    """
    ln -s $proj_dir/${run_id}/summary_* .
    ln -s $proj_dir/${run_id}/mergefcounts_gencode_mrna_gene.txt .
    
    jupyter nbconvert --to notebook --execute $notebook_path_unstranded --output ${run_id}_notebook_PE_unstranded.ipynb --ExecutePreprocessor.timeout=2678400 --allow-errors --debug

    jupyter nbconvert --to html --execute $notebook_path_unstranded --output ${run_id}_notebook_PE_unstranded.html --ExecutePreprocessor.timeout=2678400 --allow-errors --debug
    """
}


