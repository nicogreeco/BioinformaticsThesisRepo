# Script Description

- **plot_snps.R**: This script generates genomic visualizations by plotting the SNP density and distribution around the genomic coordinates of three specific lncRNAs. 
It begins by importing necessary files and setting parameters. The code then processes genetic data, filtering, organizing, and extracting relevant SNPs within a specified interval. 
It converts gene IDs to gene symbols using pre-existing databases. The script constructs visualization tracks, including ideograms, genome axes, and data tracks representing SNP density. 
Finally, it generates plots for each lncRNA locus, saving them as PNG files.

- **processing_data_1.R**: The script filters and organizes the data, ensuring that relevant information regarding genomic coordinates and variant IDs is added. 
Additionally, this script prepares the dataset for further analysis by matching it with GTEx eQTLs data, generating a refined dataset for downstream analyses and saving it for future use.

- **processing_data_2.R**: This script performs several steps to process and refine genetic variation data downloaded from the dbSNP database.
  It includes simplifying and further categorizing the function classes of the variants, creating a hierarchy graph to remove redundant classes, and substituting nucleotides in the variations with the corresponding nucleotides from the GnomAD database.
  The resulting processed data is then saved for further analysis.

- **statistics.R**: In this script, the SNPs data is processed and analyzed to obtain important statistical insights.
  The script imports the data, organizes it, and extracts relevant information such as gene frequencies, variation types, and clinical significance frequencies.
  It also includes functions to calculate mutation frequencies and create histograms to visualize the distribution of GnomAD frequencies.
  The results are saved in separate files and some graphs are also plotted for visualization purposes.

