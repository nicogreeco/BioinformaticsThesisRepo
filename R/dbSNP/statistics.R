
### lncMB1 --> LOC105371730
### lncMB2 --> LOC100507403
### lncMB3 --> LOC105378520
# lnc_snps_GnomAD
lnc <- "lncMB1"
lnc_gene <- "LOC105371730"
stat <- TRUE
## --------------------------------------------------------------------------- ##
#### importing the data and processing ##### 
## --------------------------------------------------------------------------- ##

### inport
#setwd("~/Univeristy/Tesi/R")
setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R/db_snp")

# lnc_snps_GnomAD <- read.table(paste0("./data/snp_processed/",lnc, "_1Mb_GnomAD.txt"), sep = "\t", header = T)
# lnc_snps_GnomAD <- read.table("./data/dbSNPs/snp_lncMB1.txt", sep = "\t", header = T)
lnc_snps_GnomAD <- read.table(paste0("./data/snp_processed/simplified_function_class/",lnc, "_1Mb_Gn_fc_snv.txt"), sep = "\t", header = T)


library(ggplot2)
library(progress)

### processing

# funtion to split GnomAD into two and assign the freq to Gnomad and the nuc to GnomAD_Nucleotide
add_GnomAD_values <- function(x) {
  # Split the string by the colon separator
  split_string <- strsplit(x, split=":")[[1]]
  # Extract the nucleotide and numeric parts
  nucleotide <- split_string[1]
  numeric_value <- as.numeric(split_string[2])
  # Assign the nucleotide and numeric values to 'GnomAD_Nucleotide' and 'GnomAD' columns respectively
  return(c(GnomAD_Nucleotide=nucleotide, GnomAD=numeric_value))
}

# Use sapply() to apply the function to 'GnomAD' column of the dataframe
lnc_snps_GnomAD[, c('GnomAD_Nucleotide', 'GnomAD')] <- t(sapply(lnc_snps_GnomAD$GnomAD, add_GnomAD_values))
lnc_snps_GnomAD$GnomAD <-as.numeric(lnc_snps_GnomAD$GnomAD)

# substitution of lnc
# Split the contents of the 'gene' column in 'lnc_snps_GnomAD' into separate rows using ; as delimiter
lnc_snps_GnomAD <- tidyr::separate_rows(lnc_snps_GnomAD, gene, sep = ";")
lnc_snps_GnomAD$gene[lnc_snps_GnomAD$gene == lnc_gene] <- lnc

# function to get the order of the genes based on the mean position of snps

get_order_genes <- function(db_snp) {
  # initialize a df that will contain the avarage snps position per gene
  a <- unique(db_snp$gene)
  b <- integer(length(a))
  pos <- data.frame(gene = a, position = b)
  
  for (gen in a) {
    subset <- subset(db_snp, gene == gen)
    avarage_pos <- sum(subset$pos, na.rm = TRUE)/nrow(subset)
    pos[pos$gene == gen, "position"] <- avarage_pos
  }
  
  pos <- pos[pos$gene != "",]
  pos <- pos[order(pos$position),]
  return(pos$gene)
}
pos <- get_order_genes(lnc_snps_GnomAD)
# substitute empty sttrrings with na  in Genes
lnc_snps_GnomAD$gene[lnc_snps_GnomAD$gene == ""] <- NA


lnc_snps_GnomAD$gene <- factor(x = lnc_snps_GnomAD$gene, levels = pos)



pdf(file = paste0("./data/Graphs/",lnc,".pdf"),
    width = 10, height = 7,
    paper = "USr")

## --------------------------------------------------------------------------- ##
#### -FUNCTIONS for exon and gene length #### NOT FINISHED 
###
### COMMENTA CON GPT
# genes exon length for function normalization
find_total_exon_length <- function(gene_df) {
  gene_df <- gene_df[order(gene_df$start, gene_df$end), ]
  last_processed_nc <- 0
  total_width <- 0
  for (i in 1:nrow(gene_df)) {
    row <- gene_df[i,]
    if (row$end > last_processed_nc) {
      if (row$start > last_processed_nc) {
        total_width <- total_width + row$width
        last_processed_nc <- row$end
      } else { 
        total_width <- total_width + (row$end - last_processed_nc)
        last_processed_nc <- row$end
      }
      
    }
    
  }
  return(total_width)
}

find_exon_length_genes <- function(genes,gtf) {
  lnc_subset_gtf <- subset(gtf, gene_name %in% genes)
  genes <- unique(lnc_subset_gtf$gene_name)
  a <- genes
  b <- integer(length(a))
  result <- data.frame(gene = a, width = b)
  for (gen in genes) {
    gene_df <- subset(lnc_subset_gtf, gene_name == gen)
    gene_df <- gene_df[gene_df$type == "exon",]
    gene_df <- subset(gene_df, transcript_biotype != "nonsense_mediated_decay")
    width <- find_total_exon_length(gene_df)
    result[result$gene == gen, "width"] <- width
  }
  return(result)
}

find_gene_length <- function(genes,gtf) {
  lnc_subset_gtf <- subset(gtf, gene_name %in% genes)
  genes <- unique(lnc_subset_gtf$gene_name)
  a <- genes
  b <- integer(length(a))
  result <- data.frame(gene = a, width = b)
  for (gen in genes) {
    gene_df <- subset(lnc_subset_gtf, gene_name == gen)
    width <- gene_df[gene_df$type == "gene",]$width
    if (length(width) > 1) {
      width <- max(width)
    }
    result[result$gene == gen, "width"] <- width
  }
  return(result)
}

#### exon and gene length #### NOT FINISHED
library(GenomicRanges)

gtf <- rtracklayer::import('./data/Homo_sapiens.GRCh38.109.chr.gtf')
gtf <- as.data.frame(gtf)
genes <- unique(lnc_snps_GnomAD$gene)
genes <- genes[!is.na(genes)]

# functions to find from gtf files the cumulative length of exons per each gene
# and the gene length
genes_length <- find_gene_length(genes,gtf)
exons_length <- find_exon_length_genes(genes,gtf)


#### -FUNCTIONS for genes frequencies ##### 
## --------------------------------------------------------------------------- ##

# sto a fa una funzione per contare i geni 
# in lnc$gene anche con quelli insieme
# "ZNF703;LOC124901932;LOC124901933"
# unique_gen <- function(vect_genes) {
#  genes <- c()
#  for (i in vect_genes) {
#    if (grepl(";", i)) {
#      genes <- append(genes, unlist(strsplit(i, ";")))
#    } else {
#      genes <- append(genes, i)
#    }
#    
#  }
#  return(table(genes))
#}

# it seems to do not work but i dont want to debug
# same function but using lapply
# to make it much shorter
unique_gen <- function(vect_genes) {
  genes <- unlist(lapply(vect_genes, function(x) unlist(strsplit(x, ";"))))
  return(as.data.frame(table(genes)))
}




## --------------------------------------------------------------------------- ##
#### table of the GENE FREQUENCIES ##### 
## --------------------------------------------------------------------------- ##

#genes_frequ <- unique_gen(lnc_snps_GnomAD$gene)
#(genes_frequ)


genes_frequ <- table(lnc_snps_GnomAD$gene)
genes_frequ <- as.data.frame(genes_frequ)
colnames(genes_frequ) <- c("gene", "Freq")
# saving the path in a variable and using it for the write table
path <- paste0("./data/Statistic_Results/", lnc, "/GenesFrequencies_", lnc, ".txt")

# create the file if it dose not already exist
if (!file.exists(path)) {
  write.table(genes_frequ, file = path, sep="\t", row.names=FALSE, col.names = TRUE)
}


## plotting the snps per gene
#levels(genes_frequ$genes) <- c(levels(genes_frequ$genes), lnc)
#genes_frequ$genes[genes_frequ$genes == lnc_gene] <- lnc
#genes_frequ$genes <- droplevels(genes_frequ$genes)

### NOT FINISHED 
normaized_gene_snps <- data.frame(gene = genes_frequ$gene, mut_rate = integer(length(genes_frequ$gene)))
for (gen in genes_frequ$gene) {
  if (gen %in% genes_length$gene) {
    snps <- as.integer(genes_frequ[genes_frequ$gene == gen, "Freq"])
    length <- as.integer(genes_length[genes_length$gene == gen, "width"])
    normaized_gene_snps[normaized_gene_snps$gene == gen, "mut_rate"] <- snps/length
   
  }
  
}

### Some graphs
breaks <- pos
colours <- ifelse(breaks == lnc, "red", "black")
faces <- ifelse(breaks == lnc, "bold", "plain")


ggplot(genes_frequ, aes(x = factor(gene, level = pos), y = Freq, fill = "red")) +
  geom_bar(stat = "identity", color = "black", fill = "pink") +
  labs(x = "Genes", y = "Number of SNPs", title = "SNPs per Genes") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1, 
                                   color = colours,
                                   face = faces))



ggplot(normaized_gene_snps, aes(x = factor(gene, level = pos), y = mut_rate, fill = "red")) +
  geom_bar(stat = "identity", color = "black", fill = "pink") +
  labs(x = "Genes", y = "Number of SNPs", title = "SNPs per Genes") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1, 
                                   color = colours,
                                   face = faces)) + ylim(c(0,1))

# creating a new df sub_lnc_snps only containing some col
# use it mainly for graphs. 
sub_lnc_snps <- lnc_snps_GnomAD[,c("snp_id","clinical_significance","gene","GnomAD")]


# tidyr::separate_rows separate values in the specified column according to the sep
# spified and create onerow for each values of it
# f.e. "benign;likely-benign" is splitted into "benign" and "likely-benign" and two rows 
# for the same snp are created each with one of the two
sub_lnc_snps <- tidyr::separate_rows(sub_lnc_snps, clinical_significance, sep = ";")


# empty as NA and chategorical variable as factors
sub_lnc_snps[sub_lnc_snps == ""] <- NA
sub_lnc_snps$clinical_significance <- as.factor(sub_lnc_snps$clinical_significance)
sub_lnc_snps$gene <- as.factor(sub_lnc_snps$gene)
sub_lnc_snps <- subset(sub_lnc_snps, grepl("pathogenic", clinical_significance))
sub_lnc_snps <- subset(sub_lnc_snps, clinical_significance != "conflicting-interpretations-of-pathogenicity")



# plot the snps per genes colured by clinc sign


breaks <- levels(as.factor(sub_lnc_snps$gene))
colours <- ifelse(breaks == lnc, "red", "black")
faces <- ifelse(breaks == lnc, "bold", "plain")

ggplot(sub_lnc_snps, aes(x = as.factor(gene), fill = clinical_significance)) + geom_bar() +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1, 
                                   color = colours,
                                   face = faces)) +
  ggtitle("SNPs Pathogenicity per Gene") + labs(fill = "Clinical Significance class") +
  ylab("Count") + xlab("Gene")






 
## --------------------------------------------------------------------------- ##
#### -FUNCTIONS of the variation type #####
## --------------------------------------------------------------------------- ##


# this is a function that simply return the table() function as df 
find_freq_variant_type <- function(dbSNPs) {
  
  table_vt <- table(dbSNPs$variant_type)
  
  table_vt <- as.data.frame(table_vt)
  
  colnames(table_vt) <- c("variant_type", "freq")
  
  return(table_vt)
  
}


 
## --------------------------------------------------------------------------- ##
#### table of the VARIATION TYPE #####
## --------------------------------------------------------------------------- ##

if (stat == TRUE) {
  variant_type <- find_freq_variant_type(lnc_snps_GnomAD)
  (variant_type)
  
  # path for save the file
  path <- paste0("./data/Statistic_Results/", lnc, "/VariantType_", lnc, ".txt")
  
  # save the file if it dose not exist
  if (!file.exists(path)) {
    write.table(variant_type, path , sep="\t", row.names=FALSE, col.names = TRUE)
  }
  
  ggplot(variant_type, aes(x = variant_type, y = freq, fill = "red", )) +
    geom_bar(stat="identity", color = "black", fill="pink") +
    labs(x = "Variant Type",
         y = "Number of SNPs",
         title = "Variation by Type") +
    scale_y_continuous(trans = "log10") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  
  
  
}


### --This is a little code to plot a graph of the variation type 

lnc <- "lncMB3"
path <- paste0("./data/Statistic_Results/", lnc, "/VariantType_", lnc, ".txt")
variant_type <- read.table(path, sep = "\t", header = T)
row.names(variant_type) <- variant_type$variant_type
variant_type$variant_type <- NULL
colnames(variant_type) <- lnc 
lncMB1 <- variant_type
lncMB2 <- variant_type
lncMB3 <- variant_type
variant_type <- cbind(lncMB1, lncMB2, lncMB3)

prop2 <- prop.table(as.matrix(variant_type), margin = 2)
barplot(prop2)

library("reshape2")
data <- variant_type
data_long <- as.data.frame(data)    # Reshape data from wide to long
data_long$subgroup <- rownames(data_long)
data_long <- melt(data_long, id.vars = "subgroup")
data_long                           # Printing long data

ggplot(data_long,            # Create ggplot2 plot scaled to 1.00
       aes(x = variable,
           y = value,
           fill = subgroup)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())

## --------------------------------------------------------------------------- ##
#### -FUNCTIONS for clinical significance #####
## --------------------------------------------------------------------------- ##

## obtain the set of unique factor of clinical_significance column
# a <- unique(lnc_snps_GnomAD$clinical_significance)
# b <- unlist(lapply(a, function(x) unlist(strsplit(x, ";"))))
# unique(b) 
# here is the unique set of values:
# [1] "risk-factor"                                  "benign"                                      
# [3] "likely-benign"                                "uncertain-significance"                      
# [5] "benign-likely-benign"                         "conflicting-interpretations-of-pathogenicity"
# [7] "pathogenic"                                   "likely-pathogenic"                           
# [9] "pathogenic-likely-pathogenic"                 "not-provided"
# [11] "drug-response" 


# substitute empty strings with NA values
lnc_snps_GnomAD$clinical_significance[lnc_snps_GnomAD$clinical_significance == ""] <- NA


## function to count the clinical_significance frewuencies
unique_clin_sign <- function(clinical_significance) {
  clin_sign <- unlist(lapply(clinical_significance, function(x) unlist(strsplit(x, ";"))))
  clin_sign <- factor(x = clin_sign, levels = c("risk-factor","benign","likely-benign",
                  "uncertain-significance","benign-likely-benign","conflicting-interpretations-of-pathogenicity",
                  "pathogenic","likely-pathogenic","not-provided","pathogenic-likely-pathogenic", "drug-response"))
  # this is needed to split values of some rows. clinc sign data for some snps are recorded in a semicolumn separated string:
  # "benign-likely-benign;uncertain-significance;likely-benign"
  # I firstly separate them using lapply, than unlist to create a loneg vector containing each occurrence
  # furthermore i save clin_sign as a factor, specifing the levels (ie values it can take) in order to obtain
  # a value zero if a clinical_significance status is not present
  clin_sign_table <- as.data.frame(table(clin_sign)) # tranform it in df
  colnames(clin_sign_table) <- c("clinical_significance", "freq") # change col names
  return(clin_sign_table)
  # I finally create a table with the vectro
}




### imma create a function to calculate the clinical_significance for each genes
clinic_var_for_gene <- function(dbSNP) {
  
  
  
  
  ## find the unique list of genes in the dbsnp dtaframe:
  genes <- unlist(lapply(dbSNP$gene, function(x) unlist(strsplit(as.character(x), ";"))))
  genes <- unique(genes)  #obtain the set of unique genes that appear in the df
  
  ## find the unique list of clinical_significance status: 
  clin_sign <- unlist(lapply(dbSNP$clinical_significance, function(x) unlist(strsplit(x, ";"))))
  clin_sign <- unique(clin_sign)
  
  ## initialize the df that will be returned 
  clin_per_gene <- data.frame(matrix(ncol = length(clin_sign) + 1, nrow = 0)) 
  colnames(clin_per_gene) <- c("gene", clin_sign) # set column names
  # the dataframe clin_per_gene will have the first column for the genes
  # and the other for the severeal clinical_significance status:
  # in each of them there will be saved the freq of that status in the gene of that row
  
  # to create a processing bar for the function
  # library(progress)
  # pb <- progress_bar$new(total = length(genes), format = "(:spin) [:bar] :percent")
  
  ## now i count the freq of clinical_significance status for each gene
  for (i in genes) {
    # select for each gene i
    # only the rows in wich i appear in gene col
    # even if in semicolumn separated string
    subset_gene <- subset(dbSNP, grepl(i, gene)) 
    subset_gene_table <- as.data.frame(t(unique_clin_sign(subset_gene$clinical_significance)))
    colnames(subset_gene_table) <- subset_gene_table[1,]
    subset_gene_table <- subset_gene_table[-1,] # im creating a df one row and same structure of clin_per_gene
    subset_gene_table$gene <- i                 # with data of freq of clinical_significance status for i gene
    # than ill add this to the clin_per_gene, one gene at time
    clin_per_gene <- rbind(clin_per_gene, subset_gene_table)
    # pb$inc() # increment the progress bar
  }                                              
  rownames(clin_per_gene)<- c()
  return(clin_per_gene)
  
  
}



## --------------------------------------------------------------------------- ##
#### table of the CLINICAL SIGNIFICANCE #####
## --------------------------------------------------------------------------- ##

if (stat == TRUE) {
  
  # df with frequencies for each clinical_significance status
  clin_sign_freq <- unique_clin_sign(lnc_snps_GnomAD$clinical_significance)
  (clin_sign_freq)
  
  # path for save the file
  # path <- paste0("./data/Statistic_Results/", lnc, "/ClinicalSignificanceFreq_", lnc, ".txt")
  
  # save the file if it dose not exist
  # if (!file.exists(path)) {
  #  write.table(clin_sign_freq, path , sep="\t", row.names=FALSE, col.names = TRUE)
  #}
  
  na_clin_sign <- sum(is.na(lnc_snps_GnomAD$clinical_significance))
  (na_clin_sign)
  
  
  aaaaa <- clinic_var_for_gene(lnc_snps_GnomAD)
  
  
}


# number of na in clinical_significance 
# lncMB2 -- 770183
# lncMB1 -- 751542
# lncMB3 -- 739995


# path for save the file
path <- paste0("./data/Statistic_Results/", lnc, "/ClinicalSignificanceForGene_", lnc, ".txt")

# save the file if it dose not exist
 if (!file.exists(path)) {
  write.table(aaaaa, path , sep="\t", row.names=FALSE, col.names = TRUE)
}


## --------------------------------------------------------------------------- ##
#### -FUNCTIONS for mutation freq #####
## --------------------------------------------------------------------------- ##



## split the variation column into two coloum one for reference and one for mutated
split_variation <- function(db_snp) {
  ## split the variation column into two xoloum
  db_snp$reference <- NA
  db_snp$mutant <- NA
  pb <- progress_bar$new(total = nrow(db_snp))
  for (row in 1:nrow(db_snp)) {
    # Extract the reference nucleotide using base R - substring() function
    db_snp$reference[row] <- substring(db_snp$variation[row], 1, 1)
    # Extract the mutant nucleotide using base R - substring() function
    db_snp$mutant[row] <- substring(db_snp$variation[row], 3, 3)
    pb$tick()
  }
  return(db_snp)
}

# same function but using apply so much more efficient
split_variation <- function(db_snp) {
  db_snp$reference <- apply(db_snp, 1, function(x) substr(x["variation"], 1, 1))
  db_snp$mutant <- apply(db_snp, 1, function(x) substr(x["variation"], 3, 3))
  return(db_snp)
}




## --------------------------------------------------------------------------- ##
#### table of MUTATION FREQUENCIES #####
## --------------------------------------------------------------------------- ##



## here im workin only on rows with non NA Gnomad and variant_type = snv
lnc_snps <- subset(lnc_snps_GnomAD, variant_type == "snv")
lnc_snps <- subset(lnc_snps, !is.na(GnomAD))

## split the variation column into two coloum one for reference and one for mutated
lnc_snps <- split_variation(lnc_snps)



freq_table <- round(as.matrix(prop.table(table(lnc_snps$reference, lnc_snps$mutant))),4)
names(dimnames(freq_table)) = c("Reference", "Mutant")




## --------------------------------------------------------------------------- ##
#### processing for graph of FUNCTION CLASS #####
## --------------------------------------------------------------------------- ##

sub_lnc_snps <- lnc_snps_GnomAD


# Split the contents of the 'extra_simplified_function_class' column in 'sub_lnc_snps' into separate rows
sub_lnc_snps <- tidyr::separate_rows(sub_lnc_snps, extra_simplified_function_class, sep = ";")
sub_lnc_snps <- subset(sub_lnc_snps, extra_simplified_function_class != "" & extra_simplified_function_class != "intron_variant" & extra_simplified_function_class != "intragenic" )
sub_lnc_snps <- subset(sub_lnc_snps, !is.na(extra_simplified_function_class))


# sub_lnc_snps$gene <- factor(sub_lnc_snps$gene, levels = pos)


breaks <- pos
colours <- ifelse(breaks == lnc, "red", "black")
faces <- ifelse(breaks == lnc, "bold", "plain")


ggplot(sub_lnc_snps, aes(x = gene, fill = extra_simplified_function_class)) + geom_bar() +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1, 
          color = colours,
          face = faces)) +  labs(fill = "Variation Class from SO ") +
  ggtitle("SNPs Function CLass per Gene")



#### graphics statistics #####


##----- histigrams with diff susbet of the df and bin width------##



if (stat == TRUE) {
  
  # ggplot(lnc_snps_GnomAD, aes(x=GnomAD)) + 
  #  geom_histogram(binwidth=.005, colour="black", fill="pink")
  
  ggplot(subset(lnc_snps_GnomAD,GnomAD <= 0.01) , aes(x=GnomAD)) + 
    geom_histogram(binwidth=.0002, colour="black", fill="pink") + 
    ggtitle("Distribution of GnomAD frequencies < 1% ")
  
  # ggplot(subset(lnc_snps_GnomAD,GnomAD <= 0.005) , aes(x=GnomAD)) + 
  #  geom_histogram(binwidth=.0001, colour="black", fill="pink")
  
  # ggplot(subset(lnc_snps_GnomAD,GnomAD <= 0.0005) , aes(x=GnomAD)) + 
  #   geom_histogram(binwidth=.00001, colour="black", fill="pink")
  
  # subset of snp with frew > 0.0001 (0.01%)
  ggplot(subset(lnc_snps_GnomAD,GnomAD >= 0.0001) , aes(x=GnomAD)) + 
    geom_histogram(binwidth=.01, colour="black", fill="pink") +
    ylim(0, 1000) +
    xlim(0, 0.6) + ggtitle("Distribution of GnomAD frequencies > 0.01% ")
  
  
  na_GnomA = sum(is.na(lnc_snps_GnomAD$GnomA))
  # number of na in GnomAD
  # lncMB2 -- 412584
  # lncMB1 -- 411433
  # lncMB3 -- 393494
  
  
  ratio_small = nrow(subset(lnc_snps_GnomAD,GnomAD <= 0.0001))/(nrow(lnc_snps_GnomAD) - na_GnomA)
  # ratio of snps with freq lower than 0.0001 (0.01%) GnomAD freq over the totsl non na snps
  # lncMB2 -- 0.8785488
  # lncMB1 -- 0.8679929
  # lncMB3 -- 0.8551886
  
  
  ratio_big = nrow(subset(lnc_snps_GnomAD,GnomAD >= 0.01))/(nrow(lnc_snps_GnomAD) - na_GnomA)
  # ratio of snps with freq higher than 0.01 (1%) GnomAD freq over the totsl non na snps
  # lncMB2 -- 0.0239943
  # lncMB1 -- 0.0339884
  # lncMB3 -- 0.0469724
  
}


##---------------------------------------------------------------##


dev.off()







