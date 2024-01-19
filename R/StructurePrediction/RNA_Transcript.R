

setwd("~/Univeristy/Tesi")
source(file = '~/Univeristy/Tesi/R/lnc_eQTLs/source/lib_and_functions.R')



#### for spliced transcript ####
# Read the contents of the file
# Extract the first row as a character string
file_path <- "./lncMB1_structure_prediction.vienna.txt"
file_content <- readLines(file_path)
first_row <- file_content[1]

lnc_gene <- "LOC105371730"
lnc_ENSG <- "ENSG00000214708"
chr <- 17

gtf <- rtracklayer::import('C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp/data/Homo_sapiens.GRCh38.109.chr.gtf')
gtf_lnc <- subset(gtf, gene_id == lnc_ENSG & type == "exon")
length(gtf_lnc) 

exons <- as.data.frame(ranges(gtf_lnc))

setwd("~/Univeristy/Tesi/R/lnc_eQTLs")
lnc_eQTLs <- load_eQTLs("lncMB1")
snps <- unique(lnc_eQTLs[,c(3,6,17)])
five_id <- as.list(table(lnc_eQTLs$snp_id))
five_id <- as.numeric(names(five_id[five_id>25]))

snps <- subset(snps, snp_id %in% five_id)

snps$pos_in_transc <- NA
for (snp in snps$snp_id) {
  row <- snps[snps$snp_id == snp,]
  if (row$pos > exons$start[1] & row$pos < exons$end[1]) {
    snps[snps$snp_id == snp, "pos_in_transc"] <- exons$end[1] - row$pos
  } else if (row$pos > exons$start[2] & row$pos < exons$end[2]) {
    snps[snps$snp_id == snp, "pos_in_transc"] <- exons$end[2] - row$pos + 111
  }
}


change_character <- function(sequence, index, new_character) {
  # Convert the sequence into a character vector
  sequence_vector <- strsplit(sequence, "")[[1]]
  
  # Change the character at the specified index
  sequence_vector[index] <- new_character
  
  # Join the character vector back into a string
  new_sequence <- paste(sequence_vector, collapse = "")
  
  return(new_sequence)
}
get_character <- function(sequence, index) {
  
  sequence_vector <- strsplit(sequence, "")[[1]]
  
  return(sequence_vector[index])
}

chr17_32141547_T_C_b38 <- change_character(first_row, 1277, "C")
chr17_32142404_C_A_b38 <- change_character(first_row, 420, "A")


first_row

# get exon seq

sequence_vector <- strsplit(first_row, "")[[1]]
exon1 <- sequence_vector[1:exons$width[1]]
exon2 <- sequence_vector[exons$width[1]:length(sequence_vector)]

exon1 <- paste(exon1, collapse = "")
exon2 <- paste(exon2, collapse = "")

#### for unspliced ####
# Read the contents of the file
# Extract the first row as a character string
file_path <- "./lncMB1_total_transcript.txt"
file_content <- readLines(file_path)
total_trans <- file_content[1]

### intron
sequence_vector <- strsplit(total_trans, "")[[1]]
intron <- sequence_vector[112:421] # --------- 112-421 are the coordinates of intron in the unspliced transcript
intron <- paste(intron, collapse = "")

lnc_gene <- "LOC105371730"
lnc_ENSG <- "ENSG00000214708"
chr <- 17

gtf <- rtracklayer::import('C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp/data/Homo_sapiens.GRCh38.109.chr.gtf')
gtf_lnc <- subset(gtf, gene_id == lnc_ENSG & type == "exon")
length(gtf_lnc) 

exons <- as.data.frame(ranges(gtf_lnc))

setwd("~/Univeristy/Tesi/R/lnc_eQTLs")
lnc_eQTLs <- load_eQTLs("lncMB1")
snps <- unique(lnc_eQTLs[,c(3,6,17)])
five_id <- as.list(table(lnc_eQTLs$snp_id))
five_id <- as.numeric(names(five_id[five_id>25]))

snps <- subset(snps, snp_id %in% five_id)

snps$pos_in_total_transc <- NA
for (snp in snps$snp_id) {
  row <- snps[snps$snp_id == snp,]
  if (row$pos > exons$start[2] & row$pos < exons$end[1]) {
    snps[snps$snp_id == snp, "pos_in_total_transc"] <- exons$end[1] - row$pos
  }
}


change_character <- function(sequence, index, new_character) {
  # Convert the sequence into a character vector
  sequence_vector <- strsplit(sequence, "")[[1]]
  
  # Change the character at the specified index
  sequence_vector[index] <- new_character
  
  # Join the character vector back into a string
  new_sequence <- paste(sequence_vector, collapse = "")
  
  return(new_sequence)
}
get_character <- function(sequence, index) {
  
  sequence_vector <- strsplit(sequence, "")[[1]]
  
  return(sequence_vector[index])
}

chr17_32141547_T_C_b38 <- change_character(first_row, 1277, "C")
chr17_32142404_C_A_b38 <- change_character(first_row, 420, "A")



#### PLOT #####
gen <- "hg38"
chr <- paste0("chr", chr)
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
grtrack <- GeneRegionTrack(gtf_lnc, genome = gen,
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue",
                           chromosome = chr, name = "eQTLs",
                           symbol = gtf_lnc$gene_id,
                           stacking = "pack")

gtrack <- GenomeAxisTrack(add53 = TRUE,
                          add35 = TRUE)

plotTracks(list(gtrack, grtrack))

