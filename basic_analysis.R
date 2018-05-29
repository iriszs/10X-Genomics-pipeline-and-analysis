# Only for computer at the office
.libPaths(.libPaths()[2:4])
# Does adding not work? -> open disk first.
.libPaths("/media/imgorter/firstDisk/Iris/R")
.libPaths()

#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
#source("https://bioconductor.org/biocLite.R")

# Libraries
library(SingleCellExperiment)
library(scater)
library(cellrangerRkit)
library(Seurat)
library(dplyr)


#################################################################
#                   Argument Processing                         #
#################################################################


# Get the commandline arguments
args = commandArgs(trailingOnly=TRUE)

# First argument is the output path
output_path <- args[1]

# Second argument is the run ID
id <- args[2]

# Third argument is the species (mouse or human)
species <- args[3]

# split the outputpath on /
splitted_output_path <- unlist(strsplit(output_path, "/"))
# Remove the data map and the species map from the path (which are the last two)
new_output_path <- splitted_output_path[0:(length(splitted_output_path)-2)]
# Put the arguments back together
new_output_path <- trimws(paste(new_output_path, collapse = '/'))
# Create a new R analysis directory where the pdf's will be saved.
dir.create(paste0(new_output_path, "/R_analysis/", sep=""), showWarnings = FALSE)

# Quit when there are less than 3 arguments
if (length(args) < 3) {
  stop("Three arguments must be supplied!", call. = FALSE)
}


#################################################################
#                       Create SCE object                       #
#################################################################


# Cell barcodes output file from Cell Ranger
cellbarcodes <- read.table(paste0(output_path, "barcodes.tsv"))
# Gene file output from Cell Ranger
genenames <- read.table(paste0(output_path, "genes.tsv"))
# Matrix file output from Cell Ranger
molecules <- Matrix::readMM(paste0(output_path, "matrix.mtx"))

# set the rownames of the molecules matrix
rownames(molecules) <- genenames[,1]
# Set the colnames of the molecules matrix
colnames(molecules) <- paste(id, cellbarcodes[,1], sep="_")

#Convert molecules as a matrix 
molecules <- as.matrix(molecules)

# Create a Single Cell Experiment object using the molecules matrix
sceset <- SingleCellExperiment(assays = list(counts=as.matrix(molecules)))


#################################################################
#                 Quality Control & Filtering                   #
#################################################################


# Find genes that are not expressed in any cell
keep_feature <- rowSums(counts(sceset) > 0) > 0

# not expressed cells: all cells minus the cells that are above 0
not_expressed_cells  <- length(keep_feature) - length(which(keep_feature))

# Notify user how many genes will be removed with this filtering step
print(paste(not_expressed_cells, "genes that are not expressed in any cell have been removed"))

# Remove genes that are not expressed in any cell
sceset <- sceset[keep_feature, ]

# Open up a PDF to save all QC plots
pdf(paste0(new_output_path, "/R_analysis/QC.pdf"))

containsSpike <- spikeNames(sceset)
print(containsSpike)

if (is.null(containsSpike)) {
  print("Your dataset does not contain any spikes. Skipping spike-ins quality check........")
  
} else{
  isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
  plotPhenoData(sceset, aes_string(x = "total_features", y = "pct_counts_ERCC"))
  # Calculate the quality metrics
  sceset <- calculateQCMetrics(sceset, feature_controls = list(ERCC = isSpike(sceset, "ERCC")))
  
  # ERCC filter might not be possible, it has to be defined by eye? Our data does not have metadata.
  #filter_by_ERCC <- sceset$batch != "NA19098.r2" & umi$pct_counts_ERCC < 10
  #table(filter_by_ERCC)
}


# If the species is human, find the human mitochondrial genes
if (species == "Human") {
  isSpike(sceset, "MT") <- rownames(sceset) %in% 
    c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
      "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
      "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
      "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
      "ENSG00000198840")
  # Calculate QC for human mitochondrial genes
  sceset <- calculateQCMetrics(sceset, feature_controls = list(MT = isSpike(sceset, "MT")))
  
  
# If the species is mouse, find the mouse mitochondrial genes (TODO: implement)
} else if (species == "Mouse") {
  isSpike(sceset, "mt") <- rownames(sceset) %in%
    c("ENSMUSG")
  # Calculate QC for mouse mitochondrial genes
  sceset <- calculateQCMetrics(sceset, feature_controls = list(MT = isSpike(sceset, "mt")))
  

# If there is a species that is not mouse or human, give a notification to the user.
} else {
  print("This species is not supported (yet) to find mitochondrial genes")
}

# Calculate the QC for the mitocho

# Plot ratio between ERCC spike-ins and endogenous RNA
# Is a main title possible? would be: Percentage of counts in MT genes
plotPhenoData(sceset, aes_string(x = "total_features", y = "pct_counts_MT"))

# Plot histogram of the library size
hist(sceset$total_counts, breaks = 100, main = paste("Histogram of library sizes for all cells"))
# TODO: at 25%?, place the abline, or always 2500?
abline(v = 2500, col = "red")

# MT filter based on? 
#filter_by_MT <- sceset$pct_counts_MT < 10
#table(filter_by_MT)

# Define the filter to keep cells that have more than 2500 counts in a cell
filter_by_total_counts <- (sceset$total_counts > 2500)
table(filter_by_total_counts)

# Plot histogram of the detected genes
hist(sceset$total_features, breaks = 100, main = paste("Histogram of the number of unique detected genes in all cells"))
# TODO: at 75%? place the abline or always 1500?
abline(v = 1500, col = "red")

# Define the filter to remove cells that have over 7000 genes in a cell
filter_by_expr_features <- (sceset$total_features > 7000)
table(filter_by_expr_features)

# Cell filter based on the previous quality controls.
# TODO: implement different sceset$use based on spike ins or not
#sceset$use <- (filter_by_expr_features & filter_by_total_counts & filter_by_ERCC & filter_by_MT)
sceset$use <- (filter_by_expr_features & filter_by_total_counts)

plotQC(sceset, type = "highest-expression")

#close the PDF
dev.off()

#filter the genes that are 
filter_genes <- apply(counts(sceset[ , colData(sceset)$sceset]), 1, function(x) length(x[x > 1]) >= 2)
rowData(sceset)$use <- filter_genes

table(filter_genes)

# Decide here to filter or not to filter


#Confounding factors
#Normalization
#Biologic analysis