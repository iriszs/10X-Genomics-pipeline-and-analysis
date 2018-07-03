# Only for computer at the office
#.libPaths(.libPaths()[2:4])
# Does adding not work? -> open disk first.
#.libPaths("/media/imgorter/firstDisk/Iris/R")
#.libPaths()

#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
source("https://bioconductor.org/biocLite.R")


# Functions to test if package is already installed or not.

# Package test for the packages that need to be installed using install.packages()
pkgTest <- function(x){
  if (!require(x,character.only = TRUE)){
    install.packages(x, repos = "http://cran.us.r-project.org")
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# Package test for the packages that need to be installed using bioconductor
bioPkgTest <- function(x){
  if (!require(x,character.only = TRUE)){
    source("https://bioconductor.org/biocLite.R")
    biocLite(x)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}


# Libraries

pkgTest("devtools")

if (!require(scRNA.seq.funcs)){
  devtools::install_github("hemberg-lab/scRNA.seq.funcs")
}

bioPkgTest("SingleCellExperiment")
bioPkgTest("scater")
bioPkgTest("scran")
bioPkgTest("limma")

# Options

options(stringsAsFactors = FALSE)
set.seed(1234567)



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

anno <- args[4]
# split the outputpath on /
splitted_output_path <- unlist(strsplit(output_path, "/"))
# Remove the data map and the species map from the path (which are the last two)
new_output_path <- splitted_output_path[0:(length(splitted_output_path)-2)]
# Put the arguments back together
new_output_path <- trimws(paste(new_output_path, collapse = '/'))
# Create a new R analysis directory where the pdf's will be saved.
dir.create(paste0(new_output_path, "/R_analysis/", sep=""), showWarnings = FALSE)

# Quit when there are less than 3 arguments
if (length(args) < 4) {
  stop("Four arguments must be supplied!", call. = FALSE)
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
rownames(molecules) <- make.names(genenames[,2], unique = TRUE)
# Set the colnames of the molecules matrix
colnames(molecules) <- paste(id, cellbarcodes[,1], sep="_")

#Convert molecules as a matrix 
molecules <- as.matrix(molecules)

# If annotation file is not present.
if (anno == ""){
  print("No annotation file.")
  # Create a Single Cell Experiment object using the molecules matrix
  sceset <- SingleCellExperiment(assays = list(counts=as.matrix(molecules)))
} else{
  # If annotation file is present
  annotations <- read.table(paste(output_path, anno, sep= ""), sep = "\t", header = TRUE)
  # Create a Single Cell Experiment object using the molecules matrix
  sceset <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotations)
  
  
}

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


if (is.null(containsSpike)) {
  print("Your dataset does not contain any spikes. Skipping spike-ins quality check........")
  
} else{
  isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
  plotPhenoData(sceset, aes_string(x = "total_features", y = "pct_counts_ERCC")) + ggtitle("ERCC Spike-in counts")
  # Calculate the quality metrics
  sceset <- calculateQCMetrics(sceset, feature_controls = list(ERCC = isSpike(sceset, "ERCC")))
  
}


# If the species is human, find the human mitochondrial genes
if (species == "Human") {
  isSpike(sceset, "MT") <- grep(pattern = "^MT\\.", rownames(sceset), value = TRUE)
  # Calculate QC for human mitochondrial genes
  sceset <- calculateQCMetrics(sceset, feature_controls = list(MT = isSpike(sceset, "MT")))
  
  # If the species is mouse, find the mouse mitochondrial genes
} else if (species == "Mouse") {
  isSpike(sceset, "mt") <- grep(pattern = "^mt\\.", rownames(sceset), value = TRUE)
  # Calculate QC for mouse mitochondrial genes
  sceset <- calculateQCMetrics(sceset, feature_controls = list(MT = isSpike(sceset, "mt")))
  
  # If there is a species that is not mouse or human, give a notification to the user.
} else {
  print("This species is not supported (yet) to find mitochondrial genes")
}

# Plot ratio between ERCC spike-ins and endogenous RNA
# Is a main title possible? would be: Percentage of counts in MT genes
plotPhenoData(sceset, aes_string(x = "total_features", y = "pct_counts_MT")) + ggtitle("Percentage of counts in MT genes")

# Plot histogram of the library size
hist(sceset$total_counts, breaks = 100, main = paste("Histogram of library sizes for all cells"))

# Plot histogram of the detected genes
hist(sceset$total_features, breaks = 100, main = paste("Histogram of the number of unique detected genes in all cells"))

# Filters
filter_by_total_counts <- (sceset$total_counts > 25000)
table(filter_by_total_counts)

filter_by_expr_features <- (sceset$total_features > 7000)
table(filter_by_expr_features)

# 
sceset$use <- (filter_by_expr_features & filter_by_total_counts)

# Plot the genes that have the highest expression
plotQC(sceset, type = "highest-expression") + ggtitle("Genes with the highest expression")

# Close the PDF
dev.off()


#################################################################
#                         Normalization                         #
#################################################################

# Remove genes that are 'undetectable'
filter_genes <- apply(counts(sceset[ , colData(sceset)$use]), 1, function(x) length(x[x > 1]) >= 2)
rowData(sceset)$use <- filter_genes

table(filter_genes)
dim(sceset[rowData(sceset)$use, colData(sceset)$use])

# Lognormalize the raw counts
assay(sceset, "logcounts_raw") <- log2(counts(sceset) + 1)
reducedDim(sceset) <- NULL

# Only use the data that matches the quality control requirements
sceset.qc <- sceset[rowData(sceset)$use, colData(sceset)$use]
endog_genes <- !rowData(sceset.qc)$is_feature_control

# Open PDF for the normalization comparison plots
pdf(paste0(new_output_path, "/R_analysis/Normalization.pdf"))


# CPM normalization


assay(sceset, "logcounts_raw") <- log2(counts(sceset) + 1)

logcounts(sceset) <- log2(calculateCPM(sceset, use.size.factors = FALSE) + 1)

# Compare raw (log normalized) with CPM using a relative log expression plot

plotRLE(sceset[endog_genes, ], exprs_mats = list(Raw = sceset@assays$data$logcounts_raw, CPM = sceset@assays$data$logcounts), exprs_logged = c(TRUE, TRUE))


# SCRAN normalization


qclust <- quickCluster(sceset.qc)
sceset.qc <- computeSumFactors(sceset.qc, clusters = qclust)
#Can return warning, might want to suppress if possible
sceset.qc <- normalize(sceset.qc)

# Compare raw (log normalized) with SCRAN using a relative log expression plot

plotRLE(sceset[endog_genes, ], exprs_mats = list(Raw = sceset@assays$data$logcounts_raw, scran = sceset@assays$data$logcounts), exprs_logged = c(TRUE, TRUE))


# Close PDF
dev.off()
