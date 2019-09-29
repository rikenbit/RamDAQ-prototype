#!/usr/bin/env nextflow

proj_id = params.project_id
run_ids = Channel.from(params.run_ids)
          .map{[it[0]]}

process execute_nbconvert {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_id}", mode: 'copy', overwrite: true

    input:
    val proj_id
    set run_id from run_ids

    script:
    def runfolder_dir = workflow.workDir.parent + "/output_${proj_id}/${run_id}"
    def accessory_path = workflow.scriptFile.parent.parent + "/R_QCplot"
    def notebook_path_unstranded = "/accessory/RamDA-SeqQC_template_SE_unstranded_nbconvert.ipynb"

    """
    /usr/bin/docker run --rm --user=0 -v $accessory_path:/accessory -v $runfolder_dir:/data docker.io/myoshimura080822/jupyternotebook:2.4 jupyter nbconvert --to notebook --execute $notebook_path_unstranded --output /data/${run_id}_notebook_SE_unstranded.ipynb --ExecutePreprocessor.timeout=2678400 --allow-errors --debug

    /usr/bin/docker run --rm --user=0 -v $accessory_path:/accessory -v $runfolder_dir:/data docker.io/myoshimura080822/jupyternotebook:2.4 jupyter nbconvert --to html /data/${run_id}_notebook_SE_unstranded.ipynb --output /data/${run_id}_notebook_SE_unstranded.html --ExecutePreprocessor.timeout=2678400 --allow-errors --debug
    """
}


