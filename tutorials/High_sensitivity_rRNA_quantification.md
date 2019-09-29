# High sensitivity rRNA quantification

In QC of scRNA-seq data, it is important to quantify the proportions of the sources of sequenced reads: For example, reads are derived from biological sequences, artificial sequences (spike-ins), and other unknown issues (unmapped). Usualliy, procedure to quantify the proportion of read origins starts with read mapping on a reference genome. However, such procedure could underestimate rRNA-derived reads (or overestimate 'informative' biological reads) due to multi-mapping and underestimate biological reads due to difficulity in mapping spliced reads.

Thus, we devised a pipeline to classify sources of sequenced reads with high sensitivity. We classified reads into 'Informative', 'rRNA', 'tRNA', 'Spike-in', and 'Unmapped'. Note that Informative reads are derived from biological sequences other than rRNA and tRNA sequences.

## Used software in the pipeline
- HISTA2
- featureCounts

## Required files
- Reference sequences
    - Genome sequences
    - rRNA sequence in RefSeq
    - Transcriptome sequences (GENOCODE) (only mouse)
- Gene Annotations
    - GENCODE
    - RepeatMasker
    - RefSeq
