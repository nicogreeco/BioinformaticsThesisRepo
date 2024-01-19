
library(Gviz)
library(GenomicRanges)
library(ggplot2)
library(biomaRt) 


setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs")
setwd("~/Univeristy/Tesi/R/lnc_eQTLs")

#### General Stuff  #####
setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs")
setwd("~/Univeristy/Tesi/R/lnc_eQTLs")

path_result_eQTLs <- paste0("./data/GTEx/results/", lnc, "_eQTLs_results.txt")



lnc_eQTLs <- read.table(file = path_result_eQTLs, sep = "\t", header = T)



lnc_eQTLs_smll <- lnc_eQTLs
lnc_eQTLs_smll[,c("chr", "validation_status", "function_class", "frequency", "simple_func_class")] <- NULL
lnc_eQTLs_smll$tissue <- as.factor(lnc_eQTLs_smll$tissue)
lnc_eQTLs_smll$phenotype_id <- as.factor(lnc_eQTLs_smll$phenotype_id)

table(lnc_eQTLs_smll$variant_id_prefix)
# chr17_32141005 chr17_32141547 chr17_32142076 chr17_32142404 chr17_32142595 chr17_32143009 chr17_32144276 
# 30             30              6             30             10             30             28 

summary(lnc_eQTLs_smll$pos)

table(lnc_eQTLs_smll$tissue)

table(lnc_eQTLs_smll$phenotype_id, lnc_eQTLs_smll$variant_id_prefix)



#### Plot the Chromosome ####
library(Gviz)
library(GenomicRanges)
gtf <- rtracklayer::import('C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp/data/Homo_sapiens.GRCh38.109.chr.gtf')
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
gene_ENSG <- as.character(unique(lnc_eQTLs_smll$phenotype_id))
gene_ENSG <- sub("\\..*", "", gene_ENSG)
lnc_subset_gtf <- subset(gtf, gene_id %in% gene_ENSG | gene_id == lnc_ENSG)
lnc_subset_gtf <- subset(lnc_subset_gtf, type == "gene")
mcols(lnc_subset_gtf)[lnc_subset_gtf$gene_id == lnc_ENSG,"gene_name"] <- lnc

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
                           stacking = "pack")

gtrack <- GenomeAxisTrack(add53 = TRUE,
                          add35 = TRUE)

plotTracks(list(itrack, gtrack, grtrack, atrack))

atrack <- AnnotationTrack(snps_Grange, name = "SNPs", id = paste0("rs", snps_Grange$snp_id),   featureAnnotation = "id", fontcolor.feature = "darkblue")

plotTracks(list(itrack, gtrack, grtrack, atrack),  
           from = 32100000, to = 32200000)

summary(lnc_eQTLs_smll$pos)

#### Change ENSG to Symbol ####


library("biomaRt")                                                                                                                   
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))                                                                                                        
ensembl_id = as.character(unique(lnc_eQTLs$phenotype_id))
ensembl_id <- sub("\\..*", "", ensembl_id)
genes_converted <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembl_id, mart=mart)                                                                                                                 
return(genes_converted)  



#### Function for total plot ####


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
plot_graph <- function(lnc_eQTLs, lnc) {
  
  if (lnc == "lncMB1") {
    lnc_gene <- "LOC105371730"
    lnc_ENSG <- "ENSG00000214708"
    chr <- 17
  } else if (lnc == "lncMB2") {
    lnc_gene <- "LOC100507403"
    lnc_ENSG <- "ENSG00000253123"
    chr <- 8
  } else if (lnc == "lncMB3") {
    lnc_gene <- "LOC105378520"
    lnc_ENSG <- "ENSG00000278484"
    chr <- 10
  }
  
  
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
  #lnc_subset_gtf <- subset(gtf, gene_id %in% gene_ENSG | gene_id == lnc_ENSG)
  #lnc_subset_gtf <- subset(lnc_subset_gtf, type == "gene")
  #mcols(lnc_subset_gtf)[lnc_subset_gtf$gene_id == lnc_ENSG,"gene_name"] <- lnc
  
  #palette_df <- create_palette(lnc_eQTLs)
  #mcols(lnc_subset_gtf)$color <- NA
  #for (i in 1:nrow(mcols(lnc_subset_gtf))) {
  #  id <- mcols(lnc_subset_gtf)$gene_id[i]
  #  col <- palette_df[palette_df$gene == id,"color"]
  #  if (!rlang::is_empty(col)) {
  #    mcols(lnc_subset_gtf)$color[i] <- col
  #  } else {mcols(lnc_subset_gtf)$color[i] <- "gray"}
    
  #}
  
  # atrack <- AnnotationTrack(lnc_subset_gtf, name = "eQTLs")
  
  gen <- "hg38"
  chr <- paste0("chr", chr)
  atrack <- AnnotationTrack(snps_Grange, name = "SNPs")
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  biomTrack <- BiomartGeneRegionTrack(genome = "hg38", name = "Genes", transcriptAnnotation = "symbol", chromosome = chr,
                                      filters = list(ensembl_gene_id = c(gene_ENSG, lnc_ENSG)), collapseTranscripts = "meta",
                                      shape = "arrow", background.title = "#974810", fill = "#cd8f72", collapseTranscripts = "meta",cex.title = 0.6)
  
  gtrack <- GenomeAxisTrack(add53 = TRUE,
                            add35 = TRUE)
  
  plotTracks(list(itrack, gtrack, biomTrack, atrack))
  
  atrack <- AnnotationTrack(snps_Grange, name = "SNPs", id = paste0("rs", snps_Grange$snp_id),   featureAnnotation = "id", fontcolor.feature = "grey", cex = 0.5, transcriptAnnotation = "id")
  distance <- (max(lnc_eQTLs$pos) - min(lnc_eQTLs$pos))/2
  
  plotTracks(list(itrack, gtrack, biomTrack, atrack), fill = "#a50000", 
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
    xlab("Genes") + ylab("Effect on Expression") + ggtitle(title) +  ylim( -0.8,0.4) 
  
  
  
  
  
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


#### other graphs ####
library(ggplot2)
find_gene_symbol <- function(lnc_eQTLs) {
  library("biomaRt")                                                                                                                   
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))                                                                                                        
  ensembl_id = as.character(unique(lnc_eQTLs$phenotype_id))
  ensembl_id <- sub("\\..*", "", ensembl_id)
  genes_converted <- getBM(filters = "ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembl_id, mart=mart)                                                                                                                 
  return(genes_converted)  
}

lnc_eQTLs <- read.table(file = path_result_eQTLs, sep = "\t", header = T)
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

print("first Graph")
ggplot(lnc_eQTLs, aes(x = as.factor(gene_symbol), fill = tissue)) + geom_bar() +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1, 
                                   color = "darkblue")) +
  ggtitle("eQTLs per Gene") + labs(fill = "Tissue of eQTL association") +
  ylab("Count") + xlab("Gene") + ylim(0, 100)

print("second Graph")
ggplot(lnc_eQTLs, aes(x = as.factor(snp_id), fill = gene_symbol)) + geom_bar() +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1, 
                                   color = "darkblue")) +
  ggtitle("Gene association per SNP") + labs(fill = "Genes") +
  ylab("eQTLs number") + xlab("SNP id")

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

#### Complete Analysis ####

lnc <- "lncMB1"
sub_MB1 <- T # option to subset the analisys for the important 5 snps
source(file = './source/lib_and_functions.R')

pdf(file = paste0("./data/Graphs/",lnc,"_eQTLS_Graphs_Chromosome.pdf"),
    width = 10, height = 7,
    paper = "USr")


# setting the lnc_gene obj and lnc_ENSG with gene symbol and ensg id
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
# upload the eqtls result
lnc_eQTLs <- load_eQTLs(lnc)

# add a gene symbol col from biomart
lnc_eQTLs <- create_gene_symbol_col(lnc_eQTLs)

### subsetting everithing only on the 5 similiar snps
if (sub_MB1 == T) {
  summary(lnc_eQTLs$pos)
  snpss <- unique(lnc_eQTLs$snp_id)
  lowest_row <- lnc_eQTLs[which.min(lnc_eQTLs$pos), "snp_id"]
  highest_row <- lnc_eQTLs[which.max(lnc_eQTLs$pos), "snp_id"]
  elements_to_remove <- c(highest_row, lowest_row, 376459993)
  snpss <- snpss[-which(snpss %in% elements_to_remove)]
  lnc_eQTLs <- subset(lnc_eQTLs, lnc_eQTLs$snp_id %in% snpss)
  
}

### now i save some txt files, needed for the GTex onlyne database wt graph###
unique <- unique(lnc_eQTLs[,c(1:17,29)])[,c(2:6,16:18)]
unique[order(unique$gene_symbol, unique$pos),]
print <- c()
print <- apply( unique[ , c(6,7) ] , 1 , paste , collapse = "," )


write.table(unique(print),
            file = paste0("./data/GTEx/results/", lnc, "_eQTL_ENSG,Variant_id.txt"),
            sep="\t", row.names=FALSE, col.names = F, quote = F)

write.table(unique(lnc_eQTLs$variant_id),
            file = paste0("./data/GTEx/results/", lnc, "_eQTL_Variant_id.txt"),
            sep="\t", row.names=FALSE, col.names = F, quote = F)

### --- PLOTS -- ###

plot_bar_chart(lnc_eQTLs, x = lnc_eQTLs$gene_symbol, fill = lnc_eQTLs$tissue, 
               title = "eQTLs per Gene", fill.legend = "Tissue of eQTL association", xlab = "Gene" )

plot_bar_chart(lnc_eQTLs, x = lnc_eQTLs$snp_id, fill = lnc_eQTLs$gene_symbol, 
               title = "Gene association per SNP", fill.legend = "Genes", xlab = "SNP Id" )

plot_bar_chart(lnc_eQTLs, x = lnc_eQTLs$tissue, fill = as.factor(lnc_eQTLs$snp_id), 
               title = "SNPs association per Tissue", fill.legend = "SNPs id", xlab = "Tissue" )

# gtf <- rtracklayer::import('C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp/data/Homo_sapiens.GRCh38.109.chr.gtf')

## chr plot
plot_graph(lnc_eQTLs, lnc)

### study the snp - creating a df with only data of the snps
snp_subset <- lnc_eQTLs[!duplicated(lnc_eQTLs[,1]),c(2:4,6:7,9,12,14)]
snp_subset <- snp_subset[order(snp_subset$pos),]

pdf(file = paste0("./data/Graphs/",lnc,"_eQTLS_Graphs_MB1_3SNPs.pdf"),
    width = 7, height = 4.3,
    paper = "USr")

plot_scatter_plot(lnc_eQTLs, x = lnc_eQTLs$gene_symbol, y = lnc_eQTLs$slope,
                  color = as.factor(lnc_eQTLs$snp_id), color.name = "Snp_id", 
                  size = -log10(lnc_eQTLs$pval_beta), size.name = "P Value",
                  title = "Effect of eQTLs")


for (gene in unique(lnc_eQTLs$gene_symbol)) {
  title <- paste0("Effect of eQTLs for ", gene)
  lnc_eQTLs_COPRS <- lnc_eQTLs[lnc_eQTLs$gene_symbol == gene,]
  print(plot_scatter_plot(lnc_eQTLs_COPRS, x = lnc_eQTLs_COPRS$gene_symbol, y = lnc_eQTLs_COPRS$slope,
                  color = as.factor(lnc_eQTLs_COPRS$tissue), color.name = "Tissue of Association", 
                  size = -log10(lnc_eQTLs_COPRS$pval_beta), size.name = "-log10(adj p Value)",
                  title = title))
}

### This is same as the previous for loop but 
## to plot the plots in a row in the same page

scatter_plots <- list()

for (gene in c("UTP6", "COPRS")) {
  title <- paste0("Effect of eQTLs for ", gene)
  lnc_eQTLs_COPRS <- lnc_eQTLs[lnc_eQTLs$gene_symbol == gene,]
  
  # Generate the scatter plot using ggplot
  scatter_plot <- ggplot(lnc_eQTLs_COPRS, aes(x = tissue, y = slope, color = as.factor(snp_id), size = 2)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
    geom_errorbar(aes(ymin = slope - slope_se, ymax = slope + slope_se), width = 0.2, size = 0.5, position = position_dodge(width = 0.7)) +
    # scale_size_continuous(range = c(1, 6), limits = range(-log10(lnc_eQTLs_COPRS$pval_beta)), guide = FALSE) +
    scale_color_discrete(name = "SNP Category") +
    #labs(size = "-log10(adj p Value)") +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 50, hjust = 1, color = "darkblue"),
      legend.key.size = unit(0.4, "cm")
    ) +
    xlab("Tissue of Association") + ylab("Effect on Expression") +
    ggtitle(title) + ylim(-0.6, 0.2)
  
  scatter_plots[[gene]] <- scatter_plot
}

# Arrange the scatter plots in a row matrix using grid.arrange
grid.arrange(grobs = scatter_plots, nrow = 1)

### --- PLOT SOME TABLES --- ###

library(gt)
snp_subset[,c("clinical_significance","snp_id")]
table <- gt(snp_subset[,c("clinical_significance","snp_id")])
print(table)
snp_subset[,c("extra_simplified_function_class","snp_id")]
table <- gt(snp_subset[,c("extra_simplified_function_class","snp_id")])
print(table)

dev.off()

