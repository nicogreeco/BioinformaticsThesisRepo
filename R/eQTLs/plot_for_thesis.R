

lnc <- "lncMB1"
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
source(file = './source/lib_and_functions.R')

get_gene_data <- function(gene_list) {
  # a function that take as input a list of genes and return a df with gene_symbol and size
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position", "strand"), filters="ensembl_gene_id", values=gene_list, mart=human)
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  return(gene_coords[,c("hgnc_symbol","ensembl_gene_id","start_position","end_position", "size", "strand")])
}


# upload the eqtls result
lnc_eQTLs <- load_eQTLs(lnc)

summary(lnc_eQTLs$pos)
snpss <- unique(lnc_eQTLs$snp_id)
lowest_row <- lnc_eQTLs[which.min(lnc_eQTLs$pos), "snp_id"]
highest_row <- lnc_eQTLs[which.max(lnc_eQTLs$pos), "snp_id"]
elements_to_remove <- c(highest_row, lowest_row, 376459993, 79777546)
snpss <- snpss[-which(snpss %in% elements_to_remove)]
lnc_eQTLs <- subset(lnc_eQTLs, lnc_eQTLs$snp_id %in% snpss)

# add a gene symbol col from biomart
lnc_eQTLs <- create_gene_symbol_col(lnc_eQTLs)

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

# ENSG of genes that will be in the plot
gene_ENSG <- as.character(unique(lnc_eQTLs$phenotype_id))
gene_ENSG <- sub("\\..*", "", gene_ENSG)

### Create the GenomicInteractions obj ###
## Create GRange of target genes
TargetGene <- lnc_eQTLs[,c(16:18,29)]
# here i try to recover from ensemble the star and end position of every gene
TargetGene$ENSG <- sapply(strsplit(TargetGene[["phenotype_id"]], "\\."), function(x) x[1])
table <- get_gene_data(unique(TargetGene$ENSG)) # a function that take as input a list of genes and return a df with gene info

# add to each row the start and end of the gene
TargetGene$start <- NA
TargetGene$end <- NA
for (i in 1:nrow(TargetGene)) {
  id <- TargetGene$ENSG[i]
  start <- table$start_position[table$ensembl_gene_id == id]
  end <- table$end_position[table$ensembl_gene_id == id]
  TargetGene$start[i] <- start
  TargetGene$end[i] <- end
}

TargetGene$chr <- paste0("chr", chr)

# now i create the GRange obj
TargetGene <- makeGRangesFromDataFrame(TargetGene,
                                       keep.extra.columns=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chr",
                                       start.field="start",
                                       end.field="end",
                                       strand.field="*",
                                       starts.in.df.are.0based=FALSE)
## Create GRange for SNPs
SNPs <- lnc_eQTLs[,c(2:4,6,23,27)]
SNPs$chr <- paste0("chr", chr)
# Enancher <- Enancher[!duplicated(Enancher),]
SNPs <- makeGRangesFromDataFrame(SNPs,
                                     keep.extra.columns=TRUE,
                                     seqinfo=NULL,
                                     seqnames.field="chr",
                                     start.field="pos",
                                     end.field="pos",
                                     strand.field="*",
                                     starts.in.df.are.0based=FALSE)

## function GInteraction() to form a GInteraction obj is from library "InteractionSet"
## function GenomicInteractions() to form GenomicInteractions obj is from library "GenomicInteractions"
gi <-  GenomicInteractions(TargetGene, SNPs, 
                           #counts = as.integer(-log10(SNPs$pval_beta))
                           counts = abs(SNPs$slope)
                           )
mcols(gi) <- mcols(gi)[,c("counts",  "anchor1.variant_id" , "anchor1.gene_symbol", "anchor1.ENSG", "anchor2.snp_id", "anchor2.slope")]
colnames(mcols(gi)) <- c("counts", "variant_id" , "gene_symbol", "ENSG", "snp_id", "slope")



interactions.track <- InteractionTrack(name='Test', gi, chromosome=paste0("chr", chr))
displayPars(interactions.track) <- list(plot.anchors= F, plot.outside = TRUE,
                                        min.height = 1.5,
                                        min.distance = 1,
                                        background.title = "lightblue3")
gen <- "hg38"
# chr <- paste0("chr", chr)
atrack <- AnnotationTrack(snps_Grange, name = "SNPs")
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

biomTrack <- BiomartGeneRegionTrack(genome = "hg38", name = "Genes", transcriptAnnotation = "symbol", chromosome = chr,
                                    filters = list(ensembl_gene_id = c(gene_ENSG, lnc_ENSG)), collapseTranscripts = "meta",
                                    shape = "arrow", background.title = "darkblue", fill = "darkgray", collapseTranscripts = "meta",cex.title = 0.6)

gtrack <- GenomeAxisTrack(add53 = TRUE,
                          add35 = TRUE)

plotTracks(list(itrack, gtrack, interactions.track, biomTrack, atrack), sizes=c(1,1,3,1,0.5))

atrack <- AnnotationTrack(snps_Grange, name = "SNPs", id = paste0("rs", snps_Grange$snp_id),   featureAnnotation = "id", fontcolor.feature = "grey", cex = 0.5, transcriptAnnotation = "id")
distance <- (max(lnc_eQTLs$pos) - min(lnc_eQTLs$pos))/2

plotTracks(list(itrack, gtrack, biomTrack, atrack), 
           from = min(lnc_eQTLs$pos) - distance, to = max(lnc_eQTLs$pos) + distance)

