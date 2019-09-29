# RamDAQ
RamDAQ is a computational pipeline for quality control (QC) of RamDA-seq experiments.

RamDA-seq is a single-cell total RNA sequencing method. Publication is [here](https://doi.org/10.1038/s41467-018-02866-0). The detailed protocol is [here](https://bit.riken.jp/protocols/ramda-seq/) and [kit is available from Toyobo](https://www.toyobo-global.com/news/2019/release_105.html).

<img src="img/top-image.png" width="800" />

## Requirements
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

## Getting started
### 1. Installing Nextflow
1. Make sure 8 or later is installed on your computer by using the command: `java -version`
2. Enter the below commands in your terminal (The command creates a file nextflow in `~/bin`)
```
mkdir -p ~/bin
cd ~/bin
wget -qO- https://get.nextflow.io | bash
```

3. Run the classic Hello world by entering the following command: `~/bin/nextflow run hello`

### 2. Installing Docker
- For Mac and Windows users: Download installer from [here](https://docs.docker.com/install/#supported-platforms).
- For Linux users: Follow instructions for your platforms [here](https://docs.docker.com/install/#supported-platforms).

### 3. Clone this repository

```
cd
git clone rikenbit/RamDAQ
```

### 4. Downloading example FASTQ files
```
cd $HOME/RamDAQ/examples/download_fastq_files
~/bin/nextflow run download-fastq.nf -c download-fastq.config
mkdir -p $HOME/RamDAQ_example/output_RamDA_NSC/human_NSC_001
mv output_download_fastq $HOME/RamDAQ_example/output_RamDA_human_NSC/human_NSC_001/01_fastq_files
```

For more information, see [here](examples/download_fastq_files).

### 5. Preparing annotation files
See [Preparing human annotation files](tutorials/Preparing_human_annotation_files.md)

### 6. Modifying config file

First, copy config file to the directory.
```
cd $HOME/RamDAQ_example
cp ~/RamDAQ/RamDAQ_SE_unstranded_human.config .
```

Then, modify `RamDAQ_SE_unstranded_human.config` using your favorite as follows:

```
project_id = "RamDA_human_NSC"`
```
```
run_ids = [
            ["human_NSC_001"]
	    ]
```
```
maxreadLength = 50
minReadLength = 36
```

### 7. Running RamDAQ basic pipeline


```
cd $HOME/RamDAQ_example

~/bin/nextflow run ~/bin/RamDAQ/QC_SE/02_ramdaQC_SE_fastqmcf_fastQC.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.02_ramdaQC_SE_fastqmcf_fastQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/03_ramdaQC_SE_hisat2.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.03_ramdaQC_SE_hisat2.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/04_ramdaQC_SE_RSeQC.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.04_ramdaQC_SE_RSeQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/06_ramdaQC_SE_featurecounts.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.04_ramdaQC_SE_RSeQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/07_ramdaQC_SE_createnotebook.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.07_ramdaQC_SE_createnotebook.html
```

Finally, you can get a QC report in html format under `$HOME/RamDAQ_example/output_RamDA_human_NSC/human_NSC_001/${run_id}_notebook_SE_unstranded.html`.


## Usage
### 1. Preparing annotation files
- For human data, see [Preparing human annotation files](tutorials/Preparing_human_annotation_files.md)
- For mouse data, see [Preparing mouse annotation files](tutorials/Preparing_mouse_annotation_files.md)

### 2. Preparing input FASTQ files

- Naming convention: The extensions must be `.fastq.gz`.
- FASTQ files must be located in a single directory: `my_favorite_path/output_${project_id}/${run_id}/01_fastq_files` (*`my_favorite_path` can be any path!*)

#### Converting FASTQ files from a BCL file
See [tutorial on bcl2fastq](tutorials/bcl2fastq.md).

### 3. Modifying config file
The following section in `*.config` file should be changed.

- `project_id`: (string)
- `run_ids`: (string)
- `maxReadLength`: (int) Maximum read length to be retaiend after read trimming.
    - Note: For Illumina sequencer data, the last nucleotide of sequenced read  (In practice, if read length of your FASTQ file is 51, 76, or 101, you should set `maxReadLength` to 50, 75, or 101, respectively)
- `minReadLength`: (int) Minimum read length to be retaiend after read trimming. We recommend to use the half of `readLength`.


### 4-1. Running RamDAQ basic pipeline on single-end (SE) data

```
# Move to the directory where your FASTQ files are contained in `my_favorite_path/output_${project_id}/${run_id}/01_fastq_files`
cd my_favorite_path

# Copy a config file (The example here is in case of SE)
cp ~/RamDAQ/RamDAQ_SE_unstranded_human.config .

# Rewrite RamDAQ_SE_unstranded_human.config (See instruction above)

# Run RamDAQ pipeline
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/02_ramdaQC_SE_fastqmcf_fastQC.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.02_ramdaQC_SE_fastqmcf_fastQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/03_ramdaQC_SE_hisat2.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.03_ramdaQC_SE_hisat2.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/04_ramdaQC_SE_RSeQC.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.04_ramdaQC_SE_RSeQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/06_ramdaQC_SE_featurecounts.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.06_ramdaQC_SE_featurecounts.html
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/07_ramdaQC_SE_createnotebook.nf -c RamDAQ_SE_unstranded_human.config -resume -with-report log.07_ramdaQC_SE_createnotebook.html
```

The output files are save in `my_favorite_path/output_${project_id}`.

### 4-2. Running RamDAQ basic pipeline on paired-end (PE) data

```
# Make a directory with your favorite name
mkdir my_favorite_path

# Move into the directory
cd my_favorite_path

# Copy a config file (The example here is in case of PE)
cp ~/RamDAQ/ramdaQC_PE_unstranded.config .

# Modify ramdaQC_PE_unstranded.config (See instruction above)

# Run RamDAQ pipeline
~/bin/nextflow run ~/bin/RamDAQ/QC_PE/02_ramdaQC_PE_fastqmcf_fastQC.nf -c ramdaQC_PE_unstranded.config -resume -with-report log.2_ramdaQC_PE_fastqmcf_fastQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_PE/03_ramdaQC_PE_hisat2.nf -c ramdaQC_PE_unstranded.config -resume -with-report log.3_ramdaQC_PE_hisat2.html
~/bin/nextflow run ~/bin/RamDAQ/QC_PE/04_ramdaQC_PE_RSeQC.nf -c ramdaQC_PE_unstranded.config -resume -with-report log.4_ramdaQC_PE_RSeQC.html
~/bin/nextflow run ~/bin/RamDAQ/QC_PE/06_ramdaQC_PE_featurecounts.nf -c ramdaQC_PE_unstranded.config -resume -with-report log.6_ramdaQC_PE_featurecounts.html
~/bin/nextflow run ~/bin/RamDAQ/QC_PE/07_ramdaQC_PE_createnotebook.nf -c ramdaQC_PE_unstranded.config -resume -with-report log.7_ramdaQC_PE_createnotebook.html
```
The output files are save in `my_favorite_path/output_${project_id}`.


## What occurs in RamDAQ pipeline
### RamDAQ basic pipeline

1. Read trimming (fastq-mcf)
    - <img src="img/02_ramdaQC_PE_fastqmcf_fastQC.png" width="400" />
2. Read mapping (HISAT2)
    - <img src="img/03_ramdaQC_SE_hisat2.png" width="400" />
3. Expression level quantification (featureCounts)
    - <img src="img/06_ramdaQC_SE_featurecounts.png" width="400" />
4. Automatic reporting (Jupyter notebook)
    - <img src="img/07_ramdaQC_SE_createnotebook.png" width="400" />

### RamDAQ optional pipelines

- FASTQ file generation (bcl2fastq)
    - <img src="img/01_ramdaQC_SE_bcl2fastq.png" width="400" />
- Expression level quantification (Sailfish)
    - <img src="img/05_ramdaQC_SE_sailfish.png" width="400" />
- High sensitivity rRNA quantification (HISAT2, featureCounts)
    - <img src="img/03_ramdaQC_SE_hisat2_rrna.png" width="400" />
    - What is 'High sensitivity rRNA mapping' ?: See [here](/tutorials/High_sensitivity_rRNA_quantification.md)

## Computational time and resource usage
- Machine: Linux-x64 / 24 CPU / Memory 660 GB
- Data: RamDA-seq data (n=96) on human neural stem cells (NSC)
    - 1,999,153 reads/cell, 51 nt SE, GC% 44.56
- Computational time
    - 02_ramdaQC_SE_fastqmcf_fastQC: 5m 39s
    - 03_ramdaQC_SE_hisat2: 25m 34s
    - 04_ramdaQC_SE_RSeQC: 10h 22m 7s
    - 06_ramdaQC_SE_featurecounts: 7m 49s
 - Memory usage
    - up to ~6 GB

## Contact
- Issues
- Email: support-bit (at) riken (dot) jp

## Maintainers
- Mika Yoshimura
- Haruka Ozaki
