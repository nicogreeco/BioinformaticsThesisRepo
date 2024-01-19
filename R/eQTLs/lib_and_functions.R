## libraries and function for other scripts ##
#### library ####

library(stringr)
library(tidyr)
library(rtracklayer)
library(Gviz)
library(GenomicRanges)
library(ggplot2)
library(biomaRt)
library(GenomicInteractions)
library(plyranges)
#### Function ####

averagePerBin <- function(x, binsize, mcolnames=NULL){
  ### 'x': a GenomicRanges objects with non-NA seqlengths.
  ### 'binsize': a single positive integer.
  ### 'mcolnames': names of numeric metadata columns in 'x' to "average"
  ###              i.e. to propagate to the result after averaging them
  ###              on each bin.
  ### Returns a GRanges object with: (a) the same seqinfo as 'x',
  ### (b) ranges of width 'binsize' covering all the sequences in
  ### 'seqinfo(x)', and (c) the "averaged" metadata columns specified
  ### in 'mcolnames'.
  if (!is(x, "GenomicRanges"))
    stop("'x' must be a GenomicRanges object")
  if (any(is.na(seqlengths(x))))
    stop("'seqlengths(x)' contains NAs")
  bins <- IRangesList(lapply(seqlengths(x),
                             function(seqlen)
                               IRanges(breakInChunks(seqlen, binsize))))
  ans <- as(bins, "GRanges")
  seqinfo(ans) <- seqinfo(x)
  if (is.null(mcolnames))
    return(ans)
  averageMCol <- function(colname)
  {
    cvg <- coverage(x, weight=colname)
    views_list <- RleViewsList(
      lapply(names(cvg),
             function(seqname)
               Views(cvg[[seqname]], bins[[seqname]])))
    unlist(viewMeans(views_list), use.names=FALSE)
  }
  mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], averageMCol))
  ans
}
get_gene_data <- function(gene_list) {
  # a function that take as input a list of genes and return a df with gene_symbol and size
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position", "strand"), filters="hgnc_symbol", values=gene_list, mart=human)
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  return(gene_coords[,c("hgnc_symbol","ensembl_gene_id","start_position","end_position", "size", "strand")])
}
trcol <- function(rgb, alpha) {
  col <- rgb(red = rgb[1],
             green = rgb[2],
             blue = rgb[3],
             maxColorValue = 255, alpha = alpha)
  return(col)
}
load_eQTLs <- function(lnc) {
  path_result_eQTLs <- paste0("./data/GTEx/results/", lnc, "_eQTLs_results.txt")
  
  lnc_eQTLs <- read.table(file = path_result_eQTLs, sep = "\t", header = T)
  return(lnc_eQTLs)
}
create_palette <- function(lnc_eQTLs) {
  colfunc <- colorRampPalette(c("#9f9f9f", "#8b0000"))
  as.data.frame(table(lnc_eQTLs$phenotype_id))
  df_genes <- as.data.frame(table(lnc_eQTLs$phenotype_id))
  colnames(df_genes) <- c("gene", "freq")
  df_genes$gene <- sapply(df_genes$gene, function(x) sub("\\..*", "", x))
  colors <- colfunc(max(df_genes$freq))
  df_genes$color <- NA 
  for (i in 1:nrow(df_genes)) {
    n <- df_genes$freq[i]
    df_genes$color[i] <- colors[n]
  }
  return(df_genes)
}
find_gene_symbol <- function(lnc_eQTLs) {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))                                                                                                        
  ensembl_id = as.character(unique(lnc_eQTLs$phenotype_id))
  ensembl_id <- sub("\\..*", "", ensembl_id)
  genes_converted <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembl_id, mart=mart)                                                                                                                 
  return(genes_converted)  
}
plot_graph <- function(lncc, gtf) {
  lnc <- lncc
  if (lnc == "lncMB1") {
    lnc_gene <- "LOC105371730"
    lnc_ENSG <- "ENSG00000214708"
  } else if (lnc == "lncMB2") {
    lnc_gene <- "LOC100507403"
    lnc_ENSG <- "ENSG00000253123"
  } else if (lnc == "lncMB3") {
    lnc_gene <- "LOC105378520"
    lnc_ENSG <- "ENSG00000278484"
  }
  
  
  path_result_eQTLs <- paste0("./data/GTEx/results/", lnc, "_eQTLs_results.txt")
  
  lnc_eQTLs <- read.table(file = path_result_eQTLs, sep = "\t", header = T)
  
  # gtf=as.data.frame(gtf)
  
  # create a Grange obj from snps of eQTLS
  
  lnc_eQTLs_for_Grange <- lnc_eQTLs[,-c(1,5, 7:15, 18:27)]
  lnc_eQTLs_for_Grange <- lnc_eQTLs_for_Grange[!duplicated(lnc_eQTLs_for_Grange[,c(1:2, 6)]),]
  
  
  snps_Grange <- makeGRangesFromDataFrame(lnc_eQTLs_for_Grange,
                                          keep.extra.columns=T,
                                          seqinfo=NULL,
                                          seqnames.field=c("chr"),
                                          start.field="pos",
                                          end.field=c("pos"),
                                          strand.field="strand",
                                          starts.in.df.are.0based=FALSE)
  
  # using the gtf file to obtain data on the genes in the eqtls results in order to plot them with the snps
  # upload the gtf
  
  # subsetting for the genes results in eQTLs
  gene_ENSG <- as.character(unique(lnc_eQTLs$phenotype_id))
  gene_ENSG <- sub("\\..*", "", gene_ENSG)
  lnc_subset_gtf <- subset(gtf, gene_id %in% gene_ENSG | gene_id == lnc_ENSG)
  lnc_subset_gtf <- subset(lnc_subset_gtf, type == "gene")
  mcols(lnc_subset_gtf)[lnc_subset_gtf$gene_id == lnc_ENSG,"gene_name"] <- lnc
  
  palette_df <- create_palette(lnc_eQTLs)
  mcols(lnc_subset_gtf)$color <- NA
  for (i in 1:nrow(mcols(lnc_subset_gtf))) {
    id <- mcols(lnc_subset_gtf)$gene_id[i]
    col <- palette_df[palette_df$gene == id,"color"]
    if (!rlang::is_empty(col)) {
      mcols(lnc_subset_gtf)$color[i] <- col
    } else {mcols(lnc_subset_gtf)$color[i] <- "gray"}
    
  }
  
  # atrack <- AnnotationTrack(lnc_subset_gtf, name = "eQTLs")
  
  gen <- "hg38"
  chr <- paste0("chr", as.character(unique(seqnames(lnc_subset_gtf))))
  atrack <- AnnotationTrack(snps_Grange, name = "SNPs")
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  grtrack <- GeneRegionTrack(lnc_subset_gtf, genome = gen,
                             transcriptAnnotation = "symbol",
                             background.panel = "#FFFEDB",
                             background.title = "darkblue",
                             chromosome = chr, name = "eQTLs",
                             symbol = ifelse(!is.na(lnc_subset_gtf$gene_name), lnc_subset_gtf$gene_name, lnc_subset_gtf$gene_id),
                             stacking = "pack", shape = "arrow", fill=lnc_subset_gtf$color)
  
  gtrack <- GenomeAxisTrack(add53 = TRUE,
                            add35 = TRUE)
  
  plotTracks(list(itrack, gtrack, grtrack, atrack))
  
  atrack <- AnnotationTrack(snps_Grange, name = "SNPs", id = paste0("rs", snps_Grange$snp_id),   featureAnnotation = "id", fontcolor.feature = "grey", cex = 0.5, transcriptAnnotation = "id")
  distance <- (max(lnc_eQTLs$pos) - min(lnc_eQTLs$pos))/2
  
  plotTracks(list(itrack, gtrack, grtrack, atrack), fill = "#a50000", 
             from = min(lnc_eQTLs$pos) - distance, to = max(lnc_eQTLs$pos) + distance)
  
}
plot_bar_chart <- function(lnc_eQTLs, x, fill, title, fill.legend, xlab) {
  
  
  ggplot(lnc_eQTLs, aes(x = as.factor(x), fill = fill)) + geom_bar() +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 50, hjust = 1, 
                                     color = "darkblue")) +
    ggtitle(title) + labs(fill = fill.legend) +
    ylab("Count") + xlab(xlab)  
  
}
plot_scatter_plot <- function(lnc_eQTLs, x, y, color, color.name, size, size.name, title) {
  
  ggplot(lnc_eQTLs, aes(x = x, y = y, color = as.factor(color), size = size)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2)) + 
    scale_size_continuous(range = c(1, 6), limits = range(size)) +
    scale_color_discrete(name = color.name) + labs(size = size.name) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 50, hjust = 1, 
                                     color = "darkblue"),
          legend.key.size = unit(0.4, "cm")) +  
    xlab("Genes") + ylab("Effect on Expression") + ggtitle(title)
  
  
  
  
  
}
create_gene_symbol_col <- function(lnc_eQTLs) {
  genes_converted <- find_gene_symbol(lnc_eQTLs)
  lnc_eQTLs$gene_symbol <- NA
  for (i in 1:nrow(lnc_eQTLs)) {
    id <- lnc_eQTLs$phenotype_id[i]
    id <- sub("\\..*", "", id)
    if (id == lnc_ENSG) {
      lnc_eQTLs$gene_symbol[i] <- lnc
      next
    }
    symbol <- genes_converted$hgnc_symbol[genes_converted$ensembl_gene_id == id]
    if (!is.na(symbol)) {
      
      if (symbol != "") {
        lnc_eQTLs$gene_symbol[i]  <- symbol
      } else {
        lnc_eQTLs$gene_symbol[i] <- id
      }
    } else {lnc_eQTLs$gene_symbol[i] <- id}
  }
  return(lnc_eQTLs)
}

