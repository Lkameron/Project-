# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(readxl)
library(dplyr)

# Set a seed for reproducibility
set.seed(123)

# Define the data directory and check if it exists
data_dir <- "D:/Users/project data/CELLRANGER_ANALYSIS/raw_gene_bc_matrices/mm10"
if (!dir.exists(data_dir)) {
  stop("Data directory not found!")
}

# Load data
seurat_data <- Seurat::Read10X(data.dir = data_dir)

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = seurat_data, project = "raw_Matrix_RTNdata")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Data quality control by mitochondrial content
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
VlnPlot(seurat_object, features = "percent.mt", pt.size = 0) +
  geom_vline(xintercept = 5, color = "red", linetype = "dashed") +
  ggtitle("Mitochondrial Content")
# Find cells with higher mitochondrial content, using a lower threshold
higher_mito_cells <- subset(seurat_object, subset = percent.mt > 3)  # Adjust as needed
# Check the number of cells found
num_higher_mito_cells <- ncol(higher_mito_cells)
print(paste("Number of cells with mitochondrial content above 3%:", num_higher_mito_cells))
# Example of subsetting based on a specific threshold
filtered_seurat_object <- subset(seurat_object, subset = percent.mt < 5)  # This filters out cells with mitochondrial content over 5%
# Visualize the distribution of mitochondrial content
VlnPlot(seurat_object, features = "percent.mt", pt.size = 0.1) +
  ggtitle("Distribution of Mitochondrial Content")

# Visualize the distribution of nFeature_RNA and nCount_RNA
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0.1) +
  ggtitle("Quality Control Metrics")

# Check the number of cells in the new subsetted object
num_filtered_cells <- ncol(filtered_seurat_object)
print(paste("Number of cells after filtering:", num_filtered_cells))
# Save the filtered Seurat object
saveRDS(filtered_seurat_object, file = "filtered_seurat_object.rds")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Subsetting based on quality metrics
filtered_seurat_object <- subset(filtered_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Data normalization and identification of variable features
filtered_seurat_object <- NormalizeData(filtered_seurat_object) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(vars.to.regress = "percent.mt")

# Perform PCA and visualize
filtered_seurat_object <- RunPCA(filtered_seurat_object, features = VariableFeatures(object = filtered_seurat_object))
DimPlot(filtered_seurat_object, reduction = "pca", label = TRUE)

# Clustering
filtered_seurat_object <- FindNeighbors(filtered_seurat_object, dims = 1:10) %>%
  FindClusters(resolution = 0.5)
# Run UMAP for visualization
filtered_seurat_object <- RunUMAP(filtered_seurat_object, dims = 1:10)
DimPlot(filtered_seurat_object, reduction = "umap", label = TRUE)

# Find all markers for each cluster
cluster_markers <- FindAllMarkers(filtered_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Plotting function for gene markers
plot_gene_markers <- function(genes) {
  FeaturePlot(filtered_seurat_object, features = genes, min.cutoff = "q10", max.cutoff = "q90") +
    ggtitle(paste("Expression of", paste(genes, collapse=", "), "across Clusters"))
}

# Perform t-SNE
filtered_seurat_object <- RunTSNE(filtered_seurat_object, dims = 1:10)
# Visualize t-SNE
DimPlot(filtered_seurat_object, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("t-SNE Plot")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Identify variable features using the "vst" method
filtered_seurat_object <- FindVariableFeatures(filtered_seurat_object, selection.method = "vst", nfeatures = 2000)

# Identify variable features using the "vst" method
filtered_seurat_object <- FindVariableFeatures(filtered_seurat_object, selection.method = "vst", nfeatures = 2000)
# Retrieve variable features and their dispersion
variable_features <- VariableFeatures(filtered_seurat_object)
# Get the dispersion metrics
dispersion_metrics <- HVFInfo(filtered_seurat_object)  # HVF stands for Highly Variable Features
# Variable feature plot showing mean expression and dispersion
VariableFeaturePlot(filtered_seurat_object) +
  ggtitle("Variable Features")
# Plot variable features
VariableFeaturePlot(filtered_seurat_object) +
  ggtitle("Variable Features by Average Expression and Dispersion")
# Get the top 2,000 variable features
variable_features <- VariableFeatures(filtered_seurat_object)
# Use variable features for clustering
filtered_seurat_object <- FindNeighbors(filtered_seurat_object, dims = 1:10)  # Dimensions based on PCA
filtered_seurat_object <- FindClusters(filtered_seurat_object, resolution = 0.5)
# Run UMAP on the Seurat object
filtered_seurat_object <- RunUMAP(filtered_seurat_object, dims = 1:10)  # Adjust dimensions as needed
# Plot UMAP with cluster labels
DimPlot(filtered_seurat_object, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of Clusters")
# Visualize UMAP with clustering information
DimPlot(filtered_seurat_object, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Cluster Relationships")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

FeaturePlot(filtered_seurat_object, features = c("S100b", "Gfap", "Aldh1l1", "Slc1a3", "Slc1a2", "Aqp4", "Fabp7")) # vasculature associated astrocytes 
FeaturePlot(filtered_seurat_object, features = c("Mbp", "Olig2", "Olig1", "Cntn1", "Mog")) # oligodendrocytes 
FeaturePlot(filtered_seurat_object, features = c("Csf1r", "P2ry12", "Tmem119", "Trem2", "Siglech")) #Microglia 
FeaturePlot(filtered_seurat_object, features = c("Neurod1", "Gap43", "Tubb3", "Sox11")) # neuron
FeaturePlot(filtered_seurat_object, features = c("Cldn5", "Klf2", "Vwa1", "Egfl7", "Esam", "Cd34", "Kdr", "Eng", "Slc2a1", "Nampt", "Cd93")) # Endothelial 
FeaturePlot(filtered_seurat_object, features = c("Rgs5", "Myl9", "Pdgfrb", "Abcc9", "Mcam", "Cspg4", "Pcdh18")) # Pericyte cell 
FeaturePlot(filtered_seurat_object, features = c("Kcnj10", "Kcna1", "Kcnc1", "Kcnq2", "Kcnq3")) #potassium channels
FeaturePlot(filtered_seurat_object, features = c("Scn1a", "Scn3a")) #sodium channels voltage dependent
FeaturePlot(filtered_seurat_object, features = c("Clcn2", "Clcn3", "Clcn4")) #chloride channel 
FeaturePlot(filtered_seurat_object, features = c("Gria1", "Gria2", "Gria3", "Gria4")) #AMPA receptors 
FeaturePlot(filtered_seurat_object, features = c("Grin1", "Grin2a", "Grin2b")) #NMDA receptors 
FeaturePlot(filtered_seurat_object, features = c("Htr2c", "Htr2a", "Htr1a", "Tph2", "Slc6a4", "Maoa")) # serotonin receptors genes  
FeaturePlot(filtered_seurat_object, features = c("Gabbr2", "Gad1", "Gad2", "Slc32a1", "Gabra1", "Gabra2", "Gabra3")) # GABA receptor genes 

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

