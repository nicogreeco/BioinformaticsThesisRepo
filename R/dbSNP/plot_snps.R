

## ## ## ## 


#### inport of file ####
lnc <- "lncMB3"

setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp")
source(file = "C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs/source/lib_and_functions.R")

c <- c("lncMB1", "lncMB2", "lncMB3")
for (lnc in c){
  if (lnc == "lncMB1") {
    lnc_gene <- "LOC105371730"
    lnc_ENSG <- "ENSG00000214708"
    ch_len <- 	83257441
    chr <- "chr17"
  } else if (lnc == "lncMB2") {
    lnc_gene <- "LOC100507403"
    lnc_ENSG <- "ENSG00000253123"
    ch_len <- 	145138636
    chr <- 8
  } else if (lnc == "lncMB3") {
    lnc_gene <- "LOC105378520"
    lnc_ENSG <- "ENSG00000278484"
    ch_len <- 	133797422
    chr <- 10
  }
  
  lnc_snps <- as.data.frame(read.table(paste0("./data/snp_processed/simplified_function_class/",lnc, "_1Mb_Gn_fc_snv.txt"), sep = "\t", header = T))
  dim(lnc_snps)
  lnc_snps <- lnc_snps[!is.na(lnc_snps$GnomAD),]
  dim(lnc_snps)
  
  #### from LOC ID to ENSG ####
  
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  genes_in_snps <- unique(unlist(lapply(unique(lnc_snps$gene), function(x) unlist(strsplit(x, ";")))))
  
  # in this way i converd the genes in symbol
  gene_list <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "entrezgene_id"), 
                       filters="hgnc_symbol",
                       values=genes_in_snps, 
                       mart=human)
  
  # here the loc to ensembl. loc is just the entrez id with a loc in front
  gene_list <- rbind(gene_list, getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "entrezgene_id"), 
                                          filters="entrezgene_id",
                                          values= gsub("^LOC", "", genes_in_snps), 
                                          mart=human))
  
  #### Track construction ####
  
  fontsize = 22
  
  biomTrack_all <- BiomartGeneRegionTrack(genome = "hg38", chromosome = chr, 
                                          collapseTranscripts = "meta", 
                                          background.title = rgb(209,150,226, max = 255), fill = rgb(202,211,250, max = 255),
                                          shape = "arrow", filters = list(ensembl_gene_id = gene_list$ensembl_gene_id), 
                                          transcriptAnnotation = "symbol", min.height = 2.5, collapse = T, fontsize=fontsize, fontsize.group=fontsize, stackHeight=0.75)
  
  # Grange_SNPs <-  makeGRangesFromDataFrame(dplyr::sample_n(lnc_snps[lnc_snps$variant_type == "snv",c(1:6,9,11,13)], 30000),
  
  Grange_SNPs <-  makeGRangesFromDataFrame(lnc_snps[lnc_snps$variant_type == "snv", c(1:6,9,11,13)],
                                           keep.extra.columns=T,
                                           seqinfo= Seqinfo(seqnames= c(as.character(unique(lnc_snps$chr))),
                                                            seqlengths= c(ch_len),
                                                            genome= "hg38"),
                                           seqnames.field="chr",
                                           start.field="pos",
                                           end.field="pos",
                                           strand.field="*",
                                           starts.in.df.are.0based=FALSE)
  
  # for a graph of snp density per bin (in terms of snps on 100 nucleotides)
  mcols(Grange_SNPs)$density <- 1
  snp_per_bin <- averagePerBin(Grange_SNPs, 15000, mcolnames = "density")
  # snp_per_bin <- subset(snp_per_bin, start(ranges(snp_per_bin) > min(lnc_snps$pos))
  snp_per_bin <- snp_per_bin[snp_per_bin$density != 0]
  mcols(snp_per_bin)$density <- (mcols(snp_per_bin)$density)*100
  mean_snps <- nrow(lnc_snps[lnc_snps$variant_type == "snv",])/(max(lnc_snps$pos) - min(lnc_snps$pos))
  
  # for a graph of snp density per bin, by avareging the frequency (in terms of snps on 100 nucleotides)
  mcols(Grange_SNPs)$freq <- unlist(lapply(mcols(Grange_SNPs)$GnomAD, function(x) as.numeric(str_extract(x, "(?<=:)[0-9.]+"))))
  snp_per_bin_freq <- averagePerBin(Grange_SNPs, 15000, mcolnames = "freq")
  # snp_per_bin <- subset(snp_per_bin, start(ranges(snp_per_bin) > min(lnc_snps$pos))
  snp_per_bin_freq <- snp_per_bin_freq[snp_per_bin_freq$freq > 0.00001]
  
  #SNP_track <- AnnotationTrack(Grange_SNPs, name = "SNPs", background.title = "coral3", fill = "coral2")
  
  SNP_track <- DataTrack(snp_per_bin, name = "SNPs Density", background.title = "orange3", fill = "#f8dec6", fontsize=fontsize, fontsize.group=fontsize)
  SNP_track_freq <- DataTrack(snp_per_bin_freq, name = "SNPs Density", background.title = "orange3", fill = "#f8dec6", fontsize=fontsize, fontsize.group=fontsize)
  
    
  itrack <- IdeogramTrack(genome = "hg38",
                          chromosome = chr)
  gtrack <- GenomeAxisTrack(add53 = TRUE,
                            add35 = TRUE)
  #### Plot ####
  
  # plot 1
  #save_path <- paste0(lnc, ".png")
  
  #png(filename = save_path,
  #    width = 1700, height = 619, units = "px")
  
  plotTracks(list(itrack, gtrack, SNP_track, biomTrack_all), type="histogram",
             from = min(lnc_snps$pos), to = max(lnc_snps$pos))
  
  #dev.off()
  
  # plot 2
  #save_path <- paste0(lnc, "_freq.png")
  
  #png(filename = save_path,
  #    width = 1700, height = 619, units = "px")
  
  plotTracks(list(itrack, gtrack, SNP_track_freq, biomTrack_all), type="histogram",
             from = min(lnc_snps$pos), to = max(lnc_snps$pos))
  
  #dev.off()
}
  
  
  
  