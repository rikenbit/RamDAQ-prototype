# How to download example FASTQ files

## Quick start

Here, we download RamDA-seq data on human neural stem cells (NSC) from [GSE125288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125288) by the following command.

```
cd $HOME/RamDAQ/examples/download_fastq_files
nextflow run download-fastq.nf -c download-fastq.config
mkdir -p $HOME/RamDAQ_example/output_RamDA_NSC/human_NSC_001
mv output_download_fastq $HOME/RamDAQ_example/output_RamDA_NSC/human_NSC_001/01_fastq_files
```

## What happens?


1. SRA files (`.sra`) listed in `job.yml` are downloaded.
2. The SRA files are coverted into (gzipped) FASTQ format by [`pfastq_dump`](https://github.com/inutano/pfastq-dump).

**Note:** The above pipeline utilize [`download-fastq.cwl`](https://github.com/pitagora-network/pitagora-cwl/tree/master/workflows/download-fastq) based on [Common Workflow Language](https://www.commonwl.org/) (CWL). We appreciate [Tazro Inutano Ohta](https://github.com/inutano), who is a great developer and has developed `download-fastq.cwl` and `pfastq_dump`.
