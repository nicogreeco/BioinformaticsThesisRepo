# Script Description

- **GTEx_eQTLs.R**: This R script imports and processes data from the GTEx project's eQTL dataset. 
It reads in the relevant files containing eQTL data and the snp annotations for the specific lncRNAs of interest. 
The script then matches eQTLs with the snp annotations based on their variant IDs and creates a merged dataset for further analysis.

- **study_eQTLs.R**: This script provides functions for loading eQTLs data for specific lncRNAs, finding corresponding gene symbols using BioMart, and plotting relevant genomic features. It includes functions for creating bar charts and scatter plots to visualize the associations between eQTLs and genes, as well as processing the results to provide useful summary statistics and tables.

- **plot_for_thesis.R**: This script generates a genomic visualization of eQTLs for the specific lncRNA being analyzed. It extracts genetic data related to eQTLs, SNPs and target genes within the lncRNA loci, prepares the data by creating GenomeRanges objects, and utilizes packages like GenomicInteractions and Gviz for plotting the genomic interactions between SNPs and target genes.
