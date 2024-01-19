setwd("~/Univeristy/Tesi/R")
setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R")

lnc_snps <- read.table("./data/dbSNPs/snp_lncMB3_1Mb.txt", sep = "\t", header = T)
lnc_snps <- read.table("./data/dbSNPs/snp_lncMB1.txt", sep = "\t", header = T)

head(lnc_snps[,])
str(lnc_snps)
colnames(lnc_snps)
dim(lnc_snps)
#### extract Gnomad ####
unique(lnc_snps$validation_status)  
# "by-frequency" , "by-cluster" , "by-alfa"
# and their combination

unique(lnc_snps$clinical_significance)


unique(lnc_snps$gene)
# all the snps are genes = "RHOT1;LOC105371730"



head(lnc_snps$frequency)
# frequencies are in the format:
# C:0.270144:1353:1000Genomes
# C is the minor allele, 0.27 is MAF, 1353 is the sample size
# and 1000Genomes is the source 




# this is a function that ill use in the gnomad_col_freq
# given the content of the frequency column for a row (SNP) in the dbSNPS
# dataset, it returns a matrix containing for each row one study reporting
# frequency of the snp and columns with header "nucleotide","freq","sample","study"
#     nuc freq      sample  study
# [1] G   0.111649  559     1000Genomes    
# [2] G   0.075765  292     ALSPAC         
from_freq_to_matrix <- function(col_freq) {
  vec_freq = unlist(strsplit(col_freq, '\\|'))
  my_mat <- matrix(nrow = length(vec_freq), ncol = 4)
  for (i in 1:length(vec_freq)) {
    my_mat[i,] <- unlist(strsplit(vec_freq[i], ":"))
  }
  colnames(my_mat) <- c("nucleotide","freq","sample","study")
  return(as.data.frame(my_mat))
}

# this function instead extract the freq of the Gnomad study
# for each row and add it to a new column wit it.
# for each row it create the matrix freq_mat with the function above 
# than it check the index of the row in freq_mat where
# is present the study GnomAD and uses it to extrat the freq and 
# save it in the new col of the dbdataset
# 

Gnomad_col_freq <- function(snps_df) {
  snps_df$GnomAD <- NA
  for (i in 1:nrow(snps_df)){
    if (snps_df[i,"frequency"] != "") {
      freq_mat <- from_freq_to_matrix(snps_df[i,"frequency"])
      Gnomad_idx <- match("GnomAD", freq_mat$study)
      if (!is.na(freq_mat[Gnomad_idx, "nucleotide"])) {
        snps_df[i,"GnomAD"] <- paste(freq_mat[Gnomad_idx, "nucleotide"], freq_mat[Gnomad_idx, "freq"], sep = ":")
      }
    }
  }
  return(snps_df)
}

# to remove 0s values and
# replace it with NA
for (i in seq_len(nrow(lnc_snps))) {
  if (is.na(lnc_snps[i, "GnomAD"]) || lnc_snps[i, "GnomAD"] == 0) {
    lnc_snps[i, "GnomAD"] <- NA
  }
}


lnc_snps <- Gnomad_col_freq(lnc_snps)


write.table(lnc_snps, file = "./data/snp_processed/lncMB3_1Mb_GnomAD.txt", sep = "\t", row.names = FALSE)






##### match qith GTEx eqtls df #####
##### This was done on 27-03 before 

lnc_snps <- lnc_snps[,-c(6,7,8,9,10)]
# removing "clinical_significance", "validation_status" 
# and "function_class" 

add_coordinates <- function(snps_df) {
  
  # add cordinates in the format chr8:37328275
  snps_df$coord <- paste0("chr", snps_df$chr, ":",  snps_df$pos)
  
  # add variant_id in the format chr8_37331645_C_G_b38 as in GTEx eQTLs database
  snps_df$newvariation <- gsub(",.*", "", snps_df$variation) # here i mantain only firt variation ie G>C,G,T --> G>C
  snps_df$variant_id <- paste0("chr", snps_df$chr, "_", snps_df$pos, "_", 
                               gsub(pattern = "[>]", replacement = "_", snps_df$newvariation), "_b38")
  # In the above code, we used the gsub() function to replace 
  # the > character in the variation column with an underscore (_). 
  # We then used the paste0() function to concatenate values from chr, pos, variation,
  # and the suffix string "_b38" into a new column called new_col.
  snps_df$newvariation <- NULL
  return(snps_df)
}

lnc_snps <- add_coordinates(lnc_snps)


