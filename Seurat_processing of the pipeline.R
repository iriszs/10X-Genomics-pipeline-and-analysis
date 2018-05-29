#################################################################
#                      Seurat processing                        #
#################################################################


data.data <- Read10X(args[1]) 

data <- CreateSeuratObject(raw.data = data.data, min.cells = 3, min.genes = 100, project = args[2])

#if data is human, ^MT-, if data is mouse ^mt-
#maybe figure that out with HG path and MM path
#what to do when cellranger is run with both

if (species == "Human") {
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = data@data), value = TRUE)
} else if (species == "Mouse") {
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = data@data), value = TRUE)
}  else {
  stop("This pipeline does not support " + species, call. = FALSE)
}

# Calculate the mitochondrial genes percentage
percent.mito <- Matrix::colSums(data@raw.data[mito.genes, ])/Matrix::colSums(data@raw.data)
# Create PDF of the violin plot containing the mitochondrial information 
pdf(paste0(new_output_path, "/R_analysis/meta_Data_violin.pdf"))
# Add meta data to the to the Seurat object 
data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# Create PDF containing GenePlot with the mitochondrial genes and number of genes in comparison with the amount of UMI's
pdf(paste0(new_output_path, "/R_analysis/mito_gene_umi.pdf"))
par(mfrow = c(1, 2))
GenePlot(object = data, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = data, gene1 = "nUMI", gene2 = "nGene")
dev.off()

# Filter cells based on the number of genes and the mitochondrial genes percentage.
data <- FilterCells(object = data, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# Calculate the amount of cells that have been removed 
percentage_removed_cells <- (length(data@data) / length(data@raw.data)) * 100

# Give warning when over 70 % of the cells are removed with the filter.
if (percentage_removed_cells >= 70) {
  print("Removed over 70% of your cells!")
} else {
  print(paste0("Removed ", percentage_removed_cells, "% of your cells. This might influence the results."))
}
