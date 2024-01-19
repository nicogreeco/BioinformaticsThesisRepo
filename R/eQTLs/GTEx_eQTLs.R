
#### - data inport and processing - ####
setwd("~/Univeristy/Tesi/R/lnc_eQTLs")


path <- "./data/GTEx/GTEx_Analysis_v8_eQTL_EUR/Brain_Amygdala.v8.EUR.signif_pairs.txt.gz"
eQTLs <- read.table(gzfile(path), sep = "\t", header = T)
lnc_snp <- read.table(file = "~/Univeristy/Tesi/R/db_snp/data/dbSNPs/snp_lncMB3.txt", sep = "\t", header = T) 
lnc_snp_mb <- read.table(file = "~/Univeristy/Tesi/R/db_snp/data/dbSNPs/snp_lncMB3_1Mb.txt", sep = "\t", header = T) 


lnc_snp_mb <- subset(lnc_snp_mb, grepl("LOC105378520", gene) )
dim(lnc_snp_mb)
dim(lnc_snp)
## [1] 1956   10 --> obtaining snps from the 1mb ds i obtain ore tham the 
# ds from exact sequence of lncMB1

#eQTLs[3, "variant_id"] <- "chr17_32141547_C_G_b38"
#eQTLs[157, "variant_id"] <- "chr17_32143810_C_G_b38"

add_coordinates <- function(snps_df) {
  
  # add cordinates in the format chr8:37328275
  snps_df$coord <- paste0("chr", snps_df$chr, ":",  snps_df$pos)
  
  # add cordinates in the format chr8_37328275
  snps_df$variant_id_prefix <- paste0("chr", snps_df$chr, "_", snps_df$pos)
  
  return(snps_df)
}

lnc_snp_mb <- add_coordinates(lnc_snp_mb)
lnc_snp_mb <- lnc_snp_mb[,-c(6,7,8,10)]

lnc_snp_mb <- lnc_snp_mb[!duplicated(snp__),]

eQTLs$variant_id_prefix <- sapply(eQTLs$variant_id,
                      function(x)  sub("^(chr[^_]+_[^_]+).*", "\\1", x))




#### 



#### - find eQTLS from GTEx - ####

snp__ <- lnc_snp_mb
eqtls <- eQTLs


extract_rows <- function(snp_db, eQTL) {
  
  # Extract prefixes of eQTL ids that match snp_db ids
  match_ids <- unique(sapply(snp_db$variant_id, function(x) grep(x, eQTL$variant_id, value = TRUE)))
  print(length(match_ids))
  df3 <- subset(eQTL, variant_id %in% match_ids)
  
  # Join the two data frames on the id columns
  df4 <- merge(snp_db, df3, by = "variant_id")
  
  return(df4)
}


lnc_eqtls <- merge(snp__, eqtls, by = "variant_id_prefix")

#### - final function ####

add_coordinates <- function(snps_df) {
  
  # add cordinates in the format chr8:37328275
  snps_df$coord <- paste0("chr", snps_df$chr, ":",  snps_df$pos)
  
  # add cordinates in the format chr8_37328275
  snps_df$variant_id_prefix <- paste0("chr", snps_df$chr, "_", snps_df$pos)
  
  return(snps_df)
}

## --------------------------------------------------------------------## 
###                         compete function                           ##
## - it take as input the lncRNA name in lnc and a T/F in Only_Brain - ##
## - Only_Brain == T will find eQTLs only for Brain tissues          - ##
## - default is Only_Brain = F. The function returns a df of all     - ##
## - rows of GTEx eQTLs whose id was in the lnc snp dataset, merged  - ##
## - with it. It will also save the file                             - ##
find_eQTLs <- function(lnc, Only_Brain = F, QTLs = "e") {
  library(progress)
  if (lnc == "lncMB1") {
    lnc_gene <- "LOC105371730"
  } else if (lnc == "lncMB2") {
    lnc_gene <- "LOC100507403"
  } else if (lnc == "lncMB3") {
    lnc_gene <- "LOC105378520"
  }
  
  
  ## path with all the eQTLs file
  if (QTLs == "e") {
    path_to_folder <- "C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs/data/GTEx/GTEx_Analysis_v8_eQTL_EUR"
  } else if (QTLs == "s") {
    path_to_folder <- "C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs/data/GTEx/GTEx_Analysis_v8_sQTL_EUR"
  }
  
  # obtain the file list i will iterate on
  # if Only_Brain == T it will only save name of eQTLs on Brain tissue thus iterate only on those
  # else it will iterate on all eQTLs
  if (Only_Brain == T) {
    file_list <- list.files(path = path_to_folder, pattern =  "^Brain.*\\.signif_pairs\\.txt\\.gz$")
  } else {
    file_list <- list.files(path = path_to_folder, pattern =  "\\.signif_pairs\\.txt\\.gz$")
  }
  
  
  
  ## importing the snps df      # CHANGED THE DIRECTORY FOR THE DESKTOP
  path_snp <- paste0("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp/data/snp_processed/simplified_function_class/",lnc, "_1Mb_GnomAD_function_class.txt")
  lnc_snp_mb <- read.table(file = path_snp , sep = "\t", header = T)
  ## subset the lnc df with only lnMB
  lnc_snp_mb <- subset(lnc_snp_mb, grepl(lnc_gene, gene))
  
  
  ## processing the data
  # add coordinates to snps df in the format of the eQTLs file "chr17_32143810_C_G_b38"
  # but just until the snps id

  lnc_snp_mb <- add_coordinates(lnc_snp_mb)
  
  # removing duplicates
  lnc_snp_mb <- lnc_snp_mb[!duplicated(lnc_snp_mb),]
  
  
  print(paste0("The snps annotated on the ", lnc , " are ", nrow(lnc_snp_mb)))
  
  all_matches <- data.frame()
  pb <- progress_bar$new(total = length(file_list))
  for (file in file_list) {
    cat("\n")
    cat("\n")  
    cat("\n")
    cat("\n")
    ## saving the tissue name
    tissue <- sub("\\..*", "", file)
    print(paste0("GTEx dataset containing eQTLs for ", tissue ))
    
    ## upload the data
    
    path_eQ <- paste0(path_to_folder, "/", file)
    eQTLs <- read.table(gzfile(path_eQ), sep = "\t", header = T)
    print(paste0("GTEx dataset contains ", nrow(eQTLs), " rows"))
    
    eQTLs <- eQTLs[!duplicated(eQTLs),]
    
      # forat of the snps id from "chr17_32143810_C_G_b38" to "chr17_32143810"
      # as the one we created 
    eQTLs$variant_id_prefix <- sapply(eQTLs$variant_id,
                                      function(x)  sub("^(chr[^_]+_[^_]+).*", "\\1", x))
    
      # "^" means to start at the beginning of the string.
      # "chr" is matching the literal string "chr".
      # "[^_]+" means to match one or more of any character except underscore.
      # "_" is matching the literal underscore character.
      # This pattern is repeated with "[^_]+" to match the next substring after the underscore.
      # Finally, ".*" matches any characters that follow this point.
      
      # The replacement string is "\\1". This means that the pattern that was matched should 
      # be replaced with the first subpattern, which was captured by enclosing it in parentheses. 
      # In the case of the pattern used here, this will be the substring from the beginning of the 
      # variant_id up to (but not including) the second underscore.
    
    eQTLs$tissue <- tissue
    
    ## find matching eQTLs
   
    temp <- merge(lnc_snp_mb, eQTLs, by = "variant_id_prefix")
    
    
    # add the 'tissue' column to the temporary data frame
    # temp$tissue <- tissue
    # append the new matching rows to the all_matches data frame using rbind()
    all_matches <- rbind(all_matches, temp)
    
    print(paste0("The eqtls for ", lnc , " are ", nrow(temp)))
  
    pb$tick()
  }
  
  
  if (Only_Brain == T) {
    if (QTLs == "e") {
      write.table(all_matches, file = paste0("./data/GTEx/results/",lnc, "_brain_eQTLs_results.txt"), sep = "\t", row.names = FALSE)
    } else if (QTLs == "s") {
      write.table(all_matches, file = paste0("./data/GTEx/results/",lnc, "_brain_sQTLs_results.txt"), sep = "\t", row.names = FALSE)
    }
  } else {
    if (QTLs == "e") {
      write.table(all_matches, file = paste0("./data/GTEx/results/",lnc, "_eQTLs_results.txt"), sep = "\t", row.names = FALSE)
    } else if (QTLs == "s") {
      write.table(all_matches, file = paste0("./data/GTEx/results/",lnc, "_sQTLs_results.txt"), sep = "\t", row.names = FALSE)
    }
  }
  
  
  return(all_matches)
}


setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs")
setwd("~/Univeristy/Tesi/R/lnc_eQTLs")

c <- c("lncMB1", "lncMB2", "lncMB3")
for (lnc in c) {
  find_eQTLs(lnc, QTLs = "s")
}



#### - trans ####
setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs")
lnc <- "lncMB1"
if (lnc == "lncMB1") {
  lnc_gene <- "LOC105371730"
} else if (lnc == "lncMB2") {
  lnc_gene <- "LOC100507403"
} else if (lnc == "lncMB3") {
  lnc_gene <- "LOC105378520"
}
path_snp <- paste0("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp/data/snp_processed/simplified_function_class/",lnc, "_1Mb_GnomAD_function_class.txt")
lnc_snp_mb <- read.table(file = path_snp , sep = "\t", header = T)
lnc_snp_mb <- add_coordinates(lnc_snp_mb)

path_eQ <- "C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/lnc_eQTLs/data/GTEx/GTEx_Analysis_v8_trans_eGenes_fdr05.txt"
eQTLs <- read.table(gzfile(path_eQ), sep = "\t", header = T)
eQTLs <- eQTLs[!duplicated(eQTLs),]

eQTLs$variant_id_prefix <- sapply(eQTLs$variant_id,
                                  function(x)  sub("^(chr[^_]+_[^_]+).*", "\\1", x))

temp <- merge(lnc_snp_mb, eQTLs, by = "variant_id_prefix")


# NO TRANS eQTLS ON MBs
 ## here just wandering wheater all snps within lnc locus are annotated for it as gene. It is
lnc_snp_mb_lnc <- subset(lnc_snp_mb, grepl(lnc_gene, gene))
lnc_snp_mb_lnc_2 <- subset(lnc_snp_mb, summ$Min. < pos & summ$Max. > pos)
dimsumm <- as.list(summary(lnc_snp_mb_lnc$pos))
