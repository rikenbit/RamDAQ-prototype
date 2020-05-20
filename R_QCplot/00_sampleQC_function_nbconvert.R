library("plyr")
library("dplyr")
library("tidyr")
library("stringr")
library("ggplot2")
library("reshape2")
library("matrixStats")
library("Rtsne")
library("edgeR")
library("gridExtra")
library("repr")
library("cowplot")
library("umap")

#library("factoextra")
#library("rgl")
#library("VennDiagram")
#library("pvclust")
#library("plotly")
#library("destiny")
#library("gplots")
#library("clusterSim")
#library('spatgraphs')
#library("igraph")
#library("RColorBrewer")
#library("plotrix")
#library("grid")
#library("gridExtra")
#library("lattice")
#library("caret")
#library("pheatmap")
#library("ggbiplot")
#library("entropy")
#library("tagcloud")
#library("GGally")
#library("vioplot")
#library("wesanderson")
#library("ShortRead")


trim_dataset <- function(data, pair=F, trimobj_1="_R1", trimobj_2="_R2"){
    
    if ("Sample_ID" %in% colnames(data)){
        data$name = data$Sample_ID
    }

    data$name = str_replace_all(data$name, "_trim", "")

    if (pair){
        data$readname = str_replace_all(str_replace_all(data$name, as.character(trimobj_1), ""),as.character(trimobj_2), "")
    }
    
    if ("X" %in% colnames(data)){
        data = data[,-1]    
    }
    return(data)
}


SDset <- function(x, y){
  data.frame(x=x, y=y, upper = mean(y)+3*sd(y),lower = mean(y)-3*sd(y))
}

plot_totalseq <- function(plotdata, title, outdir, nonplot=NULL){
  .e = environment()
  
  if (!is.null(nonplot)){
    plotdata = subset(plotdata, !name %in% nonplot)
  }

  plotdata = SDset(plotdata$name, plotdata$totalseq)

  g = ggplot(plotdata, aes(x=x,y=y), environment=.e) +
    geom_bar(alpha=0.7, stat="identity") +
    xlab("Sample") + ylab("TotalSeq") + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(paste(title, "", sep=" "))
  return(g)
  #ggsave(file = paste0(outdir,"/fastqc_totalseq_",title,".png"), plot=g, dpi=100, width=12, height=5)

}

plot_perGC <- function(plotdata, title, outdir, nonplot=NULL){
  .e = environment()

  if (!is.null(nonplot)){
    plotdata = subset(plotdata, !name %in% nonplot)
  }

  plotdata = SDset(plotdata$name, plotdata$perGC)

  g = ggplot(plotdata, aes(x=x,y=y), environment=.e) +
    geom_bar(alpha=0.7, stat="identity") +
    ylim(0, 100) +
    xlab("Sample")+ylab("%GC") +
    #scale_x_discrete(labels=plotdata$samplename) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(paste(title, "", sep=" "))
  return(g)
  #ggsave(file = paste0(outdir,"/fastqc_perGC_",title,".png"), plot=g, dpi=100, width=12, height=5)
}

plot_assignedgene_rate <- function(fastqcdata, readDistdata, title, outdir, ylim, pair=F, nonplot=NULL){
  
  readDistdata$name = str_replace_all(readDistdata$name, ".sort", "")

  if (!is.null(nonplot)){
    fastqcdata = subset(fastqcdata, !name %in% nonplot)
    readDistdata = subset(readDistdata, !name %in% nonplot)
  }

  if (pair){
    fastqcdata = fastqcdata %>% group_by(readname) %>% summarise(totalseq_sum = sum(totalseq))
    plotdata = dplyr::left_join(readDistdata, fastqcdata, by=c("name"="readname"))
    plotdata$assignedgene_rate = (plotdata$totalread / plotdata$totalseq_sum) *100
  } else {
    plotdata = dplyr::left_join(fastqcdata[,c("name", "totalseq")], readDistdata, by=c("name"))
    plotdata$assignedgene_rate = (plotdata$totalread / plotdata$totalseq) *100
  }

  plotdata = SDset(plotdata$name, plotdata$assignedgene_rate)
  
  g = ggplot(plotdata, aes(x=x,y=y)) +
    geom_bar(alpha=0.7, stat="identity") +
    ylim(0, as.integer(ylim)) +
    xlab("Sample")+ylab("Assigned gene rate") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(paste(title, "", sep=" "))
  #ggsave(file = paste0(outdir,"/hisat2_alinedgene_rate_",title,".png"), plot=g, dpi=100, width=12, height=5)
  return(g)
}

plot_readDist_summary <- function(plot_data, title, outdir, nonplot=NULL){

  if (!is.null(nonplot)){
    plot_data = subset(plot_data, !name %in% as.character(nonplot))
  }

  plot_data$per_intergenic = (plot_data$totaltag-plot_data$assignedtag)/plot_data$totaltag *100
  plot_data$per_cds = plot_data$cds/plot_data$totaltag *100
  plot_data$per_utr5 = plot_data$utr5/plot_data$totaltag *100
  plot_data$per_utr3 = plot_data$utr3/plot_data$totaltag *100
  plot_data$per_intron = plot_data$intron/plot_data$totaltag *100
  plot_data$per_tssup10 = plot_data$tssup10/plot_data$totaltag *100
  plot_data$per_tesdown10 = plot_data$tesdown10/plot_data$totaltag *100

  rownames(plot_data) = plot_data$name

  plot_data = plot_data[, grepl("per_",colnames(plot_data))]
  plot_data$id = rownames(plot_data)

  g = ggplot(melt(plot_data),aes(x=id,y=value,fill=variable)) +
    geom_bar(stat="identity",position="stack",colour="gray32",show.legend=T) +
    xlab("dataset") + ylab("percent") +
    ggtitle(title) +
    theme(axis.text.y=element_text(size=7), legend.text=element_text(size=8)) +
    coord_flip()

  return(g)
  #ggsave(file = paste0(outdir,"/RSeQC_readDist_",title,".png"), plot=g, dpi=100, width=10, height=30)
}

plot_geneBodyCov_heatmap <- function(rseqc_data,title,xaxis,outdir,nonplot=NULL,toMatch=NULL){

  if (!is.null(nonplot)){
    rseqc_data = subset(rseqc_data, !name %in% as.character(nonplot))
  }
  plot_data = rseqc_data

  rownames(plot_data) = rseqc_data$name
  plot_data = plot_data[order(rownames(plot_data)),]
  
  colnames(plot_data) = c("id", paste("val",formatC(1:100,width=3,flag="0"),sep=""))
  plot_data_t = (melt(plot_data,id.var = c("id")))

  if (is.null(toMatch)){
    g = ggplot(plot_data_t,aes(x=variable,y=value,group=id)) +
      geom_line(aes(color=id)) +
      xlab("Gene body percentile (5'->3')") + ylab("Coverage") +
      geom_text(data=plot_data_t[plot_data_t$variable==xaxis,],
                aes(label=id,colour=id),angle=90,vjust=0.5,hjust=1,show.legend=F) +
      ggtitle(title) +
      theme(axis.text.x=element_blank(),legend.position = "none")
    return(g)
    #ggsave(file = paste0(outdir,"/geneBC_line_",title,".png"), plot=g, dpi=100, width=9, height=6)
  }else{
    g = ggplot(plot_data_t,aes(x=variable,y=value,group=id)) +
      geom_line(aes(color=id)) +
      xlab("Gene body percentile (5'->3')") + ylab("Coverage") +
      geom_text(data=plot_data_t[plot_data_t$variable==xaxis & grepl(paste(toMatch,collapse="|"),plot_data_t$id),],
                aes(label=id,colour=id),angle=90,vjust=0.5,hjust=1,show.legend=F) +
      ggtitle(paste(title, "summary of gene body coverage", sep=" ")) +
      theme(axis.text.x=element_blank(),legend.position = "none", plot.title=element_text(size=10,face="bold"))
    return(g)
    #ggsave(file = paste0(outdir,"/geneBC_line_",title,".png"), plot=g, dpi=100, width=9, height=6)
  }
}

plot_bar_inferexp <- function(plot_data, title, outdir, nonplot=NULL){

  if (!is.null(nonplot)){
    plot_data = subset(plot_data, !name %in% as.character(nonplot))
  }
  plot_data = plot_data[,c("name", "undetermined_fraction", "sense_fraction", "antisense_fraction")]

  plot_data = melt(plot_data, id="name")
  g = ggplot(plot_data, aes(name, value, fill = variable)) +
    geom_bar(stat = "identity") +
    ggtitle(title) +
    theme(axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y=element_text(size=7), legend.text=element_text(size=8))

  return(g)
  #ggsave(file = paste0(outdir,"/RSeQC_readDist_",title,".png"), plot=g, dpi=100, width=10, height=30)
}

plot_bar_highsenst_rrna <- function(plot_data, title, outdir, nonplot=NULL){

  if (!is.null(nonplot)){
    plot_data = subset(plot_data, !name %in% as.character(nonplot))
  }

  plot_data_tmp = plot_data
  plot_data_tmp = select(plot_data, starts_with('p_'))
  plot_data = cbind(plot_data_tmp, plot_data$name)
  colnames(plot_data) = c(colnames(plot_data_tmp), "name")

  plot_data = melt(plot_data, id="name")
  g = ggplot(plot_data, aes(name, value, fill = variable)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(size=7, angle=90, hjust=1), axis.text.y=element_text(size=7), legend.text=element_text(size=8))

  return(g)
  #ggsave(file = paste0(outdir,"/RSeQC_readDist_",title,".png"), plot=g, dpi=100, width=10, height=30)
}

plot_assigned_rate <- function(plotdata, title, outdir, ylim, nonplot=NULL){
  
  if (!is.null(nonplot)){
    plotdata = subset(plotdata, !name %in% as.character(nonplot))
  } 
  plotdata = SDset(plotdata$name, plotdata$assigned_rate)
  
  g = ggplot(plotdata, aes(x=x,y=y)) +
    geom_bar(alpha=0.7, stat="identity") +
    #geom_ribbon(aes(ymin=lower, ymax=upper, group=type, fill=type), alpha=0.25) +
    ylim(0, as.integer(ylim)) +
    xlab("Sample")+ylab("Hista2 alined-rate") +
    #scale_x_discrete(labels=plotdata$x) + 
    theme(axis.text.x=element_text(size=6, angle=90, hjust=1), legend.text=element_text(size=8)) +
    ggtitle(paste(title, "alined rate of hisat2", sep=" "))
 return(g) 
 #ggsave(file = paste0(outdir,"/hisat2_alined_rate_",title,".png"), plot=g, dpi=100, width=12, height=5)
}

trim_countfc_sample <- function(data, exclude_id){

  #data = subset(data, !grepl("ERCC", Geneid))
  rownames(data) = data$Geneid
  data = data[,-c(1,2,3,4,5,6,7), drop=FALSE]
  samplename = str_replace(colnames(data), "_trim", "")
  colnames(data) = c(samplename)

  data = data[,!colnames(data) %in% exclude_id, drop=FALSE]
  return(data)
}

calc_tpm <- function(counts,len) {
  x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
}

calc_tpm_counts <- function(counts, gene_length, exclude_name){

    counts_trim = trim_countfc_sample(counts, exclude_name)
    counts_trim$Geneid = rownames(counts_trim)

    counts_tpm = dplyr::left_join(counts_trim, gene_length, by=c("Geneid"))
    rownames(counts_tpm) = counts_tpm$Geneid

    cat(paste0("before na_omit genenum: ", length(counts_tpm$length), "\n"))
    cat(paste0("after na_omit genenum: ", length(na.omit(counts_tpm)$length), "\n"))

    counts_tpm = calc_tpm(counts_tpm[,!colnames(counts_tpm) %in% c("Geneid","length"), drop=FALSE], as.numeric(counts_tpm$length))
    counts_tpm = as.data.frame(counts_tpm)
    counts_tpm_log = log10(counts_tpm+1)

    return(list(tpm = counts_tpm, tpm_log = counts_tpm_log))

}

calc_rpkm_counts <- function(counts, gene_length, exclude_name){

    counts_trim = trim_countfc_sample(counts, exclude_name)
    counts_trim$Geneid = rownames(counts_trim)

    counts_rpkm = dplyr::left_join(counts_trim, gene_length, by=c("Geneid"))
    rownames(counts_rpkm) = counts_rpkm$Geneid

    cat(paste0("before na_omit genenum: ", length(counts_rpkm$length), "\n"))
    cat(paste0("after na_omit genenum: ", length(na.omit(counts_rpkm)$length), "\n"))

    counts_rpkm = rpkm(counts_rpkm[,!colnames(counts_rpkm) %in% c("Geneid","length")], gene.length = as.numeric(counts_rpkm$length))
    counts_rpkm = as.data.frame(counts_rpkm)
    counts_rpkm_log = log10(counts_rpkm+1)

    return(list(rpkm = counts_rpkm, rpkm_log = counts_rpkm_log))

}

panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.05,"p<0.05",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}

panel.smooth<-function (x, y, col = "blue", bg = NA, pch = 18, 
  cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  # if (any(ok)) 
  #   lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
  #     col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

plot_detgene_bar <- function(counts, counts_log, title="", nonplot=NULL){

    if (!is.null(nonplot)){
        counts = counts[,!colnames(counts) %in% nonplot, drop=FALSE]
        counts_log = counts_log[,!colnames(counts_log) %in% nonplot, drop=FALSE]
    }

    counts_gene = counts[!grepl("ERCC", rownames(counts)),,drop=FALSE]
    counts_detgenenum = data.frame(colSums(counts_gene >1))
    colnames(counts_detgenenum) = c("NumOfGenes")
    counts_detgenenum$samplename = row.names(counts_detgenenum)

    if (ncol(counts_log) > 1 & ncol(counts_log) < 8){
      pairs(counts_log, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist, main=title)
    } else if(ncol(counts_log) == 1) {
      cat(paste0("No plot is shown when the number of samples is 1.", "\n"))
    } else {
      pairs(counts_log[,sample(ncol(counts_log), 8)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist, main=title)
    }

    g = ggplot(counts_detgenenum, aes(x=samplename,y=NumOfGenes)) +
        geom_bar(alpha=0.7, stat="identity") +
        xlab("Sample") + ylab("Number of detected genes (TPM>1)") +
        theme(axis.text.x=element_text(size=6, angle=90, hjust=1.0)) +
        ggtitle(title)
    return(g)

}

create_pca_tsne_umap <- function(countdata, perplexity, local_connectivity=1.0, n_neighbors=15){

  ### pca
  data_pca = prcomp(t(countdata), scale = FALSE, center = TRUE)
  data_pca_df = as.data.frame(data_pca$x)

  ### tSNE
  data_tsne = Rtsne(t(countdata), dims=2, perplexity=perplexity)
  data_tsne_df = as.data.frame(data_tsne$Y)

  ### UMAP
  data_umap = umap(t(countdata), local_connectivity=local_connectivity, n_neighbors=n_neighbors)
  data_umap_df = as.data.frame(data_umap$layout)

  return(list(pca=data_pca, pca_df=data_pca_df, tsne=data_tsne, tsne_df=data_tsne_df, umap=data_umap, umap_df=data_umap_df))

}

ggplot_2D <- function(dset, x, y, title, celltype, outname, label=F) {
  .e = environment()
  p=ggplot(dset, aes_string(x=x, y=y), environment=.e) +
    geom_point(alpha=.5, size = 3, aes(colour=as.factor(celltype))) +
    theme(legend.position="none") +
    ggtitle(title)
  return(p) 
  
}

create_pca_tsne_umap_mode <- function(countdata, perplexity, mode, local_connectivity=1.0, n_neighbors=15){

  if(ncol(countdata) == 1){
    cat("No plot is shown when the number of samples is 1.\n")
  } else { 

    if (mode=="pca"){
      ### pca
      data_pca = prcomp(t(countdata), scale = FALSE, center = TRUE)
      data_pca_df = as.data.frame(data_pca$x)
      return(list(data=data_pca, data_df=data_pca_df))
    } else if (mode=="tsne") {
      ### tSNE
      data_tsne = Rtsne(t(countdata), dims=2, perplexity=perplexity)
      data_tsne_df = as.data.frame(data_tsne$Y)
      return(list(data=data_tsne, data_df=data_tsne_df))
    } else {
      ### UMAP
      data_umap = umap(t(countdata), local_connectivity=local_connectivity, n_neighbors=n_neighbors)
      data_umap_df = as.data.frame(data_umap$layout)
      return(list(data=data_umap, data_df=data_umap_df))
    }
  }
}

