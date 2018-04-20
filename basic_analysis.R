### combination script with Seurat, scater and CellrangerRkit###
#Scater for determination of clusters? see heimberg website

.libPaths(.libPaths()[2:4])
.libPaths("/media/imgorter/firstDisk/Iris/R")
.libPaths()

#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
#source("https://bioconductor.org/biocLite.R")
library(cellrangerRkit)
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

output_path <- args[1]

id <- args[2]

species <- args[3]


if (length(args)==0) {
  stop("At least one argument must be supplied!", call.=FALSE)
} else if (length(args)==1) {
  # default id
  args[2] = "Study ID unknown"
}

data.data <- Read10X(args[1]) 

data <- CreateSeuratObject(raw.data = data.data, min.cells = 3, min.genes = 100, project = args[2])

#if data is human, ^MT-, if data is mouse ^mt-
#maybe figure that out with HG path and MM path
#what to do when cellranger is run with both

if (species == "Human") {
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = data@data), value = TRUE)
} else {
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = data@data), value = TRUE)
}  

percent.mito <- Matrix::colSums(data@raw.data[mito.genes, ])/Matrix::colSums(data@raw.data)

data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = data, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = data, gene1 = "nUMI", gene2 = "nGene")

data <- FilterCells(object = data, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

percentage_removed_cells <- (length(data@data) / length(data@raw.data)) * 100

if (percentage_removed_cells >= 70) {
  print("Removed over 70% of your cells!")
} else {
  print(paste0("Removed ", percentage_removed_cells, "% of your cells"))
}

