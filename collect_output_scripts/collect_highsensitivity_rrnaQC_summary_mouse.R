library(tidyverse)
library(magrittr)

args=(commandArgs(TRUE))
inputdir = args[1]
path_GENCODE_annot = args[2]
pathRefseqAnnot = args[3]
outdir = "."

result_hisat2 = sprintf("%s/summary_hisat2_rrna_results.txt", inputdir)
result_Genomic_rRNA = sprintf("%s/mergefcounts_Genomic.rRNA.txt", inputdir)
result_Genomic_tRNA = sprintf("%s/mergefcounts_Genomic.tRNA.txt", inputdir)
result_GTR = sprintf("%s/mergefcounts_G.T.R.txt", inputdir)


########################################
## mapping rate summary

dfunmap = read_tsv(result_hisat2, col_names=TRUE)
dfunmap %<>% select(-1)

outfile = sprintf("%s/hisat2_G.T.R_loose_mapped.mapping_rate.summary.txt", outdir)
write_tsv(dfunmap, outfile)

########################################
# Genomic rRNA

df_genomic_rrna = read_tsv(result_Genomic_rRNA, col_names=TRUE)
df_genomic_rrna %<>% select(-1) %>% gather(Sample_ID, count, -Geneid, -Chr, -Start, -End, -Strand, -Length)
df_genomic_rrna %<>% group_by(Sample_ID) %>% summarise(count_genomic_rRNA = sum(count)) %>% mutate(Sample_ID = str_replace_all(Sample_ID, "_rrna", ""))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_Genomic.rRNA_multiFrac.Genomic_rRNA.summary.txt", outdir)
write_tsv(df_genomic_rrna, outfile)

########################################
# Genomic tRNA

df_genomic_trna = read_tsv(result_Genomic_tRNA, col_names=TRUE)
df_genomic_trna %<>% select(-1) %>% gather(Sample_ID, count, -Geneid, -Chr, -Start, -End, -Strand, -Length)
df_genomic_trna %<>% group_by(Sample_ID) %>% summarise(count_genomic_tRNA = sum(count)) %>% mutate(Sample_ID = str_replace_all(Sample_ID, "_rrna", ""))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_Genomic.tRNA_multiFrac.Genomic_tRNA.summary.txt", outdir)
write_tsv(df_genomic_trna, outfile)

########################################
# GTR

dfGTR = read_tsv(result_GTR, col_names=TRUE)
dfGTR %<>% select(-1) %>% gather(Sample_ID, count, -Geneid, -Chr, -Start, -End, -Strand, -Length) %>%
    select(Geneid, count, Sample_ID) %>% rename(target_id = Geneid) %>% mutate(Sample_ID = str_replace_all(Sample_ID, "_rrna", ""))

### load annotation
dfGencodeAnnot = read_tsv(path_GENCODE_annot, col_names=TRUE)

#dfGencodeAnnotRrna = dfGencodeAnnot %>% filter(gene_type == "rRNA" | gene_type == "Mt_rRNA")
dfGencodeAnnotRrna = dfGencodeAnnot %>% filter(gene_type == "rRNA" | gene_type == "rRNA_pseudogene" | gene_type == "Mt_rRNA")
dfGencodeAnnotTrna = dfGencodeAnnot %>% filter(gene_type == "tRNA" | gene_type == "Mt_tRNA")

dfRefSeqAnnot = read_tsv(pathRefseqAnnot, col_names=FALSE) %>% mutate(target_id = gsub(">", "", gsub(" .+", "", X1)))

### GENCODE rRNA
dfGencodeRrna = inner_join(dfGTR, dfGencodeAnnotRrna %>% select(transcript_id) %>%
    rename(target_id=transcript_id), by="target_id") %>% group_by(Sample_ID) %>% summarise(count_GENCODE_rRNA = sum(count))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_G.T.R_multiFrac.Gencode_rRNA.summary.txt", outdir)
write_tsv(dfGencodeRrna, outfile)

### GENCODE tRNA
dfGencodeTrna = inner_join(dfGTR, dfGencodeAnnotTrna %>% select(transcript_id) %>%
    rename(target_id=transcript_id)) %>% group_by(Sample_ID) %>%
    summarise(count_GENCODE_tRNA = sum(count)) %>% mutate(Sample_ID = str_replace_all(Sample_ID, "_rrna", ""))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_G.T.R_multiFrac.Gencode_tRNA.summary.txt", outdir)
write_tsv(dfGencodeTrna, outfile)

### RefSeq rRNA
dfRefSeqRrna = inner_join(dfGTR, dfRefSeqAnnot %>% select(target_id)) %>%
    group_by(Sample_ID) %>% summarise(count_RefSeq_rRNA = sum(count)) %>% mutate(Sample_ID = str_replace_all(Sample_ID, "_rrna", ""))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_G.T.R_multiFrac.RefSeq_rRNA.summary.txt", outdir)
write_tsv(dfRefSeqRrna, outfile)

### ERCC
dfErcc = dfGTR %>% filter(grepl("ERCC-", target_id)) %>% group_by(Sample_ID) %>% summarise(count_ERCC = sum(count))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_G.T.R_multiFrac.ERCC.summary.txt", outdir)
write_tsv(dfErcc, outfile)

### SIRV
dfSirv = dfGTR %>% filter(grepl("SIRVome_isoforms", target_id)) %>% group_by(Sample_ID) %>% summarise(count_SIRV = sum(count))

outfile = sprintf("%s/featureCounts_hisat2_G.T.R_loose_mapped_G.T.R_multiFrac.SIRV.summary.txt", outdir)
write_tsv(dfSirv, outfile)


##################################################
# create summary

dfsummay =
    dfunmap %>%
    inner_join(df_genomic_rrna, by="Sample_ID") %>%
    inner_join(dfGencodeRrna, by="Sample_ID") %>%
    inner_join(dfRefSeqRrna, by="Sample_ID") %>%
    inner_join(df_genomic_trna, by="Sample_ID") %>%
    inner_join(dfGencodeTrna, by="Sample_ID") %>%
    inner_join(dfErcc, by="Sample_ID") %>%
    inner_join(dfSirv, by="Sample_ID")

dfsummay %<>%
    mutate(rRNA_read = count_genomic_rRNA + count_GENCODE_rRNA + count_RefSeq_rRNA) %>%
    mutate(tRNA_read = count_genomic_tRNA + count_GENCODE_tRNA) %>%
    mutate(mapped_read = total_reads - unmapped_read) %>%
    mutate(ERCC_read = count_ERCC) %>%
    mutate(SIRV_read = count_SIRV) %>%
    mutate(non_rRNA_tRNA_read = mapped_read - rRNA_read - tRNA_read - count_ERCC - count_SIRV)

dfsummay %<>% mutate(
    p_unmapped = unmapped_read/total_reads*100,
    p_rRNA = rRNA_read/total_reads*100,
    p_tRNA = tRNA_read/total_reads*100,
    p_non_rRNA_tRNA = non_rRNA_tRNA_read/total_reads*100,
    p_ERCC = ERCC_read/total_reads*100,
    p_SIRV = SIRV_read/total_reads*100)

outfile = sprintf("%s/featureCounts.hisat2_G.T.R_loose_mapped.summary_full.txt", outdir)
write_tsv(dfsummay, outfile)

dfsummaySimple = dfsummay %>% select(Sample_ID, total_reads, unmapped_read, mapped_read, non_rRNA_tRNA_read, rRNA_read, tRNA_read, ERCC_read, SIRV_read, matches("p_"))

outfile = sprintf("%s/featureCounts.hisat2_G.T.R_loose_mapped.summary_simple.txt", outdir)
write_tsv(dfsummaySimple, outfile)


sessionInfo()










