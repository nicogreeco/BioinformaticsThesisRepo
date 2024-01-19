setwd("C:/Users/nicco/OneDrive/Documenti/Univeristy/Tesi/R")
library(progress)


## --------------------------------------------------------------------------- ##
#### function class substitution ##### 
## --------------------------------------------------------------------------- ##

# this two function is to first simplify the classes
# i have looked on sequence ontology website the meaning of each of 
# the annotation class in the db and subsitute them with less specific
# ontology class
change_class <- function(vect) {
  funtion_classes <- unlist(strsplit(vect, ";"))
  vec <- c()
  for (class in funtion_classes) {
    if (class == "2KB_upstream_variant") {
      vec <- append(vec, "intragenic") #1
    } else if (class == "3_prime_UTR_variant") {
      vec <- append(vec, "exon_variant") #2
    } else if (class == "5_prime_UTR_variant") {
      vec <- append(vec, "exon_variant") #3
    } else if (class == "500B_downstream_variant") {
      vec <- append(vec, "intragenic") #4
    } else if (class == "coding_sequence_variant") {
      vec <- append(vec, "coding_sequence_variant") #5
    } else if (class == "downstream_transcript_variant") {
      vec <- append(vec, "intragenic") #6
    } else if (class == "frameshift_variant") {
      vec <- append(vec, "frameshift_variant") #7
    } else if (class == "genic_downstream_transcript_variant") {
      vec <- append(vec, "gene_variant") #8
    } else if (class == "genic_upstream_transcript_variant") {
      vec <- append(vec, "gene_variant") #9
    } else if (class == "inframe_deletion") {
      vec <- append(vec, "inframe_variant") #10
    } else if (class == "inframe_indel") {
      vec <- append(vec, "inframe_variant") #11
    } else if (class == "inframe_insertion") {
      vec <- append(vec, "inframe_variant") #12
    } else if (class == "initiator_codon_variant") {
      vec <- append(vec, "coding_sequence_variant") #13
    } else if (class == "intron_variant") {
      vec <- append(vec, "intron_variant") #14
    } else if (class == "missense_variant") {
      vec <- append(vec, "nonsynonimus_variant") #15
    } else if (class == "non_coding_transcript_variant") {
      vec <- append(vec, "non_coding_transcript_variant") #16
    } else if (class == "splice_acceptor_variant") {
      vec <- append(vec, "intron_variant") #17
    } else if (class == "splice_donor_variant") {
      vec <- append(vec, "intron_variant") #18
    } else if (class == "stop_gained") {
      vec <- append(vec, "nonsynonimus_variant") #19
    } else if (class == "stop_lost") {
      vec <- append(vec, "nonsynonimus_variant") #20
    } else if (class == "synonymous_variant") {
      vec <- append(vec, "coding_sequence_variant") #21
    } else if (class == "terminator_codon_variant") {
      vec <- append(vec, "coding_sequence_variant") #22
    } else if (class == "upstream_transcript_variant") {
      vec <- append(vec, "intragenic") #23
    }
  }
  return(unique(vec))
}
# i apply it to every row in order to have a simplified value
simplify_function_class <- function(db_snp) {
  db_snp$simple_func_class <- NA
  pb <- progress_bar$new(total = nrow(db_snp))
  for (row in 1:nrow(db_snp)) {
    db_snp$simple_func_class[row] <- paste(change_class(db_snp$function_class[row]), collapse = ";")
    pb$tick()
  }
  return(db_snp)
}


# there is still a lot of redundancy:
#[1] "non_coding_transcript_variant" "coding_sequence_variant"       "exon_variant"                 
#[4] "gene_variant"                  "nonsynonimus_variant"  
# coding_seq. is a subset of gene variant end also nonsyn.
#  i would like to mantain only the more specific


# i used igraph to model the hierarchy of the ontology classes
# as a directed tree that can be seen if printed
create_graph <- function() {
  library(igraph)
  g <- make_empty_graph(directed = TRUE)
  g <- add_vertices(g, nv = 9)
  V(g)$name <- c("intragenic", "exon_variant", "coding_sequence_variant", "frameshift_variant", "gene_variant", "inframe_variant", "intron_variant", "nonsynonimus_variant", "non_coding_transcript_variant")
  g <- add_edges(g, edges = c(
    "gene_variant", "non_coding_transcript_variant",
    "gene_variant", "intron_variant",
    "gene_variant", "exon_variant",
    "exon_variant", "coding_sequence_variant",
    "coding_sequence_variant", "frameshift_variant",
    "coding_sequence_variant", "inframe_variant",
    "coding_sequence_variant", "nonsynonimus_variant"))
  return(g)
}

plot(g)


# I want to usee the tree to remove the more general class
# and mantain only the specific ones, avoiding redundancy

# Define a recursive function to traverse the graph
# given an node start (a functional class) it expand to its childs
# check wheter they are in the lst (wich is the content of funtional class variable)
# if not it expand further in tree, if true it remove the start bc 
# we have found a more specific funcitional class
# ex: ("gene_variant", "exon_variant"), 
# we want to remove gene variant and mantain only the more specific one
traverse_vertices <- function(g, start, lst) {
  neighbors <- neighbors(g, start, mode = "out")
  if (length(neighbors) == 0) {
    return(FALSE)
  } else {
    for (neighbor in neighbors) {
      if (V(g)$name[neighbor] %in% lst) {
        # print(paste("Vertex", V(g)$name[neighbor], "matched the condition"))
        return(TRUE)
      }
      # print(paste("Traversing to vertex", V(g)$name[neighbor]))
      if (traverse_vertices(g, neighbor, lst)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
}


extra_simplified_function_class <- function(db_snp) {
  db_snp$extra_simplified_function_class <- NA
  pb <- progress_bar$new(total = nrow(db_snp))
  for (row in 1:nrow(db_snp)) {
    simple_func_class <- unlist(strsplit(db_snp$simple_func_class[row], ";"))
    temporary <- simple_func_class
    for (class in simple_func_class) {
      if (traverse_vertices(g, class, simple_func_class) == TRUE) {
        temporary <- temporary[temporary != class]
      }
    }
    db_snp$extra_simplified_function_class[row] <- paste(sort(temporary), collapse = ";")
    pb$tick()
  }
  return(db_snp)
}




c <- c("lncMB1","lncMB2","lncMB3")
g <- create_graph()
for (lnc in c) {
  print(paste0("Im processing:", lnc))
  lnc_snps_GnomAD <- read.table(paste0("./data/snp_processed/",lnc, "_1Mb_GnomAD.txt"), sep = "\t", header = T)
  print(paste0(lnc, ": read - ", format(Sys.time(), "%X")))
  lnc_snps_GnomAD <- simplify_function_class(lnc_snps_GnomAD)
  print(paste0(lnc, ": first processing step - ", format(Sys.time(), "%X")))
  lnc_snps_GnomAD <- extra_simplified_function_class(lnc_snps_GnomAD)
  print(paste0(lnc, ": second processing step - ", format(Sys.time(), "%X")))
  write.table(lnc_snps_GnomAD, file = paste0("./data/snp_processed/simplified_function_class/",lnc, "_1Mb_GnomAD_function_class.txt"), sep = "\t", row.names = FALSE)
  print(paste0(lnc, " has been saved"))
  
}

# --------------------------------------------------------------------------- ##
#### snv substitution ##### 
## --------------------------------------------------------------------------- ##


## now i want to mantain only the row with GnomAD

Gnomad_no_NA <- (subset(lnc_snps_GnomAD, !is.na(GnomAD)))
# 348148 rows
dim(subset(Gnomad_no_NA, nchar(variation) > 3))
# 67341 rows with more than one nuc eg T>G,A,C



## now i check if the nucleotide whose gnomad frewuencies is refered to
## is present in the variation variable

## i created a new column containing t or f according to 
## its presence in the variation variable

Gnomad_no_NA$checkGnomADTotal <- NA
pb <- progress_bar$new(total = nrow(Gnomad_no_NA))
for (row in 1:nrow(Gnomad_no_NA)) {
  nuc <- unlist(strsplit(Gnomad_no_NA$GnomAD[row], ":"))[1]
  var <- Gnomad_no_NA$variation[row]
  grepl(nuc,var)
  Gnomad_no_NA$checkGnomADTotal[row] <- grepl(nuc,var)
  # print(row)
  pb$tick()
}

Gnomad_no_NA <- subset(Gnomad_no_NA, variant_type == "snv")
# considering only snv i get 312959 rows
# and all f them were true
# with all the variant_type instead there were 2815 FALSE and 345333 TRUE

# i do thge same as before but cheking the presence of the nucleotide
# of Gnomad in the first 3 character of variation eg T>A of T>A,G,C

Gnomad_no_NA$checkGnomAD <- NA
pb <- progress_bar$new(total = nrow(Gnomad_no_NA))
for (row in 1:nrow(Gnomad_no_NA)) {
  nuc <- unlist(strsplit(Gnomad_no_NA$GnomAD[row], ":"))[1]
  var <- Gnomad_no_NA$variation[row]
  var <- substr(var, 1, 3)
  grepl(nuc,var)
  Gnomad_no_NA$checkGnomAD[row] <- grepl(nuc,var)
  # print(row)
  pb$tick()
}

check_nuc <- function(db_snp, row) {
  nuc <- unlist(strsplit(db_snp$GnomAD[row], ":"))[1]
  var <- db_snp$variation[row]
  var <- substr(var, 1, 3)
  sol <- grepl(nuc,var)
  return(sol)
}

#  for the ones that do not have the gnomad nucleotide in the first three
# i substitute it as third character gnomad nucleotide = C : T>A,G,C --> T>C

pb <- progress_bar$new(total = nrow(Gnomad_no_NA))
for (row in 1:nrow(Gnomad_no_NA)) {
  if (nchar(Gnomad_no_NA$variation[row]) > 3) {
    # print(row)
    var <- Gnomad_no_NA$variation[row]
    if (Gnomad_no_NA$checkGnomAD[row] == FALSE) {
      nuc <- unlist(strsplit(Gnomad_no_NA$GnomAD[row], ":"))[1]
      Gnomad_no_NA$variation[row] <- paste0(substr(var, 1, 2), nuc)
    } else {
      Gnomad_no_NA$variation[row] <- substr(var, 1, 3)
    }
  } 
  pb$tick()
}


## funtion to perform the whole process of substitution of the variation with the Gnomad
## one automatically 
change_nucl <- function(db_snp) {
  pb <- progress_bar$new(total = nrow(db_snp))
  for (row in 1:nrow(db_snp)) {
    if (nchar(db_snp$variation[row]) > 3) {
      var <- db_snp$variation[row]
      if(!is.na(check_nuc(db_snp, row))) {
        if(check_nuc(db_snp, row) == FALSE) { 
          nuc <- unlist(strsplit(db_snp$GnomAD[row], ":"))[1]
          db_snp$variation[row] <- paste0(substr(var, 1, 2), nuc)
        } else {
          db_snp$variation[row] <- substr(var, 1, 3)
        }
      }
    }
    pb$tick()
  }
  return(db_snp)
}


## for loop to do it on all the three dataset
c <- c("lncMB1","lncMB2","lncMB3")
for (lnc in c) {
  print(paste0("Im processing:", lnc))
  lnc_snps_GnomAD <- read.table(paste0("./data/snp_processed/simplified_function_class/",lnc, "_1Mb_GnomAD_function_class.txt"), sep = "\t", header = T)
  print(paste0(lnc, ": read - ", format(Sys.time(), "%X")))
  lnc_snps_GnomAD <- change_nucl(lnc_snps_GnomAD)
  print(paste0(lnc, ": first processing step - ", format(Sys.time(), "%X")))
  write.table(lnc_snps_GnomAD, file = paste0("./data/snp_processed/simplified_function_class/",lnc, "_1Mb_Gn_fc_snv.txt"), sep = "\t", row.names = FALSE)
  print(paste0(lnc, " has been saved"))
  
}









