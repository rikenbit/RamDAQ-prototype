# Preparing human annotation files

In the following, we use [GENCODE release 32](https://www.gencodegenes.org/human/release_32.html) (based on GRCh38).

## Requirements
- Docker

## Making a directory to gather annotation files in

```
mkdir -p $HOME/annotations/human/
```

## Copying annotation files bundled with RamDAC

```
cd $HOME/annotations/human
cp $HOME/RamDAQ/annotations/all_sequencing_WTA_adopters.fa .
cp $HOME/RamDAQ/annotations/Human_GRCh38_rmk_rRNA_RefSeq_RNA45S_all.gtf .
```

## Getting Reference genome FASTA file

Here, we download reference genome sequence from GENCODE.

```
cd $HOME/annotations/human
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```
## Making 'ref_chrsize' file
This file contains chromosome size. The following command will `GRCh38.primary_assembly.genome.fa.fai`.

```
cd $HOME/annotations/human
docker run --rm -v $PWD:$PWD -w $PWD docker.io/myoshimura080822/hisat2_set:190611 \
    samtools faidx GRCh38.primary_assembly.genome.fa
```

## Building HISAT2 index

```
cd $HOME/annotations/human
mkdir hisat2_index
docker run --rm -v $PWD:$PWD -w $PWD docker.io/myoshimura080822/hisat2_set:190611 \
    hisat2-build GRCh38.primary_assembly.genome.fa hisat2_index/GRCh38.primary_assembly.genome
```

## Getting reference gene annotations (GTF)
Here, we download reference gene annotations (GTF) from GENCODE.

```
cd $HOME/annotations/human
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
gunzip gencode.v32.primary_assembly.annotation.gtf.gz
```


## Extraciting specific gene annotations

```
cd $HOME/annotations/human
less gencode.v32.primary_assembly.annotation.gtf | grep chrM > gencode.v32.primary_assembly.annotation.mt.gtf
less gencode.v32.primary_assembly.annotation.gtf | grep HIST > gencode.v32.primary_assembly.annotation.histone.gtf
```


## Converting Reference gene annotation (BED12 format)
The following commands convert the GTF file into a BED12-formatted file.

```
cd $HOME/annotations/human
docker run --rm -v $PWD:$PWD -w $PWD julia:1.1.1 \
    julia $HOME/RamDAQ/annotations/scripts/gtfToBed12.jl \
    gencode.v32.primary_assembly.annotation.gtf \
    gencode.v32.primary_assembly.annotation.gtf.bed
```

You can also convert a GFF3 file into BED12 format by the following command:
```
docker run --rm -v $PWD:$PWD -w $PWD julia:1.1.1 \
    julia $HOME/RamDAQ/annotations/scripts/gff3ToBed12.jl $gtf $bed12
```


