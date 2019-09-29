# Preparing mouse annotation files

In the following, we use [GENCODE release M23](https://www.gencodegenes.org/mouse/release_M23.html) (based on GRCm38).

## Requirements
- Docker

## Making a directory to gather annotation files in

```
mkdir -p $HOME/annotations/mouse
```

## Copying annotation files bundled with RamDAC

```
cd $HOME/annotations/human
cp $HOME/RamDAQ/annotations/all_sequencing_WTA_adopters.fa .
cp $HOME/RamDAQ/annotations/mm10_rmsk_rRNA4.5S_Refseq_45S4.5Smerge_all.gtf .
```

## Getting Reference genome FASTA file

Here, we download reference genome sequence from GENCODE.

```
cd $HOME/annotations/mouse
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
```
## Making 'ref_chrsize' file
This file contains chromosome size. The following command will `GRCm38.primary_assembly.genome.fa.fai`.

```
cd $HOME/annotations/mouse
docker run --rm -v $PWD:$PWD -w $PWD docker.io/myoshimura080822/hisat2_set:190611 \
    samtools faidx GRCm38.primary_assembly.genome.fa
```

## Building HISAT2 index

```
cd $HOME/annotations/mouse
mkdir hisat2_index
docker run --rm -v $PWD:$PWD -w $PWD docker.io/myoshimura080822/hisat2_set:190611 \
    hisat2-build GRCm38.primary_assembly.genome.fa hisat2_index/GRCm38.primary_assembly.genome
```

## Getting reference gene annotations (GTF)
Here, we download reference gene annotations (GTF) from GENCODE.

```
cd $HOME/annotations/mouse
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
gunzip gencode.vM23.primary_assembly.annotation.gtf.gz
```


## Extraciting specific gene annotations

```
cd $HOME/annotations/mouse
less gencode.vM23.primary_assembly.annotation.gtf | grep chrM > gencode.vM23.primary_assembly.annotation.mt.gtf
less gencode.vM23.primary_assembly.annotation.gtf | grep Hist > gencode.vM23.primary_assembly.annotation.histone.gtf
```


## Converting Reference gene annotation (BED12 format)
The following commands convert the GTF file into a BED12-formatted file.

```
cd $HOME/annotations/mouse
docker run --rm -v $PWD:$PWD -w $PWD julia:1.1.1 \
    julia $HOME/RamDAQ/annotations/scripts/gtfToBed12.jl \
    gencode.vM23.primary_assembly.annotation.gtf \
    gencode.vM23.primary_assembly.annotation.gtf.bed
```

You can also convert a GFF3 file into BED12 format by the following command:
```
docker run --rm -v $PWD:$PWD -w $PWD julia:1.1.1 \
    julia $HOME/RamDAQ/annotations/scripts/gff3ToBed12.jl $gtf $bed12
```


