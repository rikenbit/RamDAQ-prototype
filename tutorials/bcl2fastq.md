# Converting FASTQ files from a BCL file
First, build a Docker container for bcl2fastq with the following command.
```
cd ~/RamDAQ/Dockerfiles/bcl2fastq
docker build -t bcl2fastq2:1.0 .
```

Then, run test.
```
$ docker run --rm bcl2fastq2:1.0 bcl2fastq --help
BCL to FASTQ file converter
bcl2fastq v2.17.1.14
Copyright (c) 2007-2015 Illumina, Inc.
...
```
Then, follow the below instructions.

```
# Make a directory with your favorite name
mkdir XXXXXX

# Move into the directory
cd XXXXXX

# Copy a config file (The example here is in case of SE)
cp ~/RamDAQ/ramdaQC_SE_unstranded.config .

# Rewrite config (with vi, vim, emacs, or your favorite editors)
vi ramdaQC_SE_unstranded.config

# Run bcl2fastq via nextflow
~/bin/nextflow run ~/bin/RamDAQ/QC_SE/01_ramdaQC_SE_bcl2fastq.nf -c ramdaQC_SE_unstranded.config -resume -with-report
```

