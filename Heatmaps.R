# Heat maps - all cell marker genes 
# Loading necessary library
library(Seurat)
library(dplyr)# Ensure dplyr is loaded for any data manipulation

# Combine all gene markers into a single vector
gene_markers <- c("S100b", "Gfap", "Aldh1l1", "Slc1a3", "Slc1a2", "Aqp4", "Fabp7", 
                  "Mbp", "Olig2", "Olig1", "Cntn1", "Mog",
                  "Csf1r", "P2ry12", "Tmem119", "Trem2", "Siglech",
                  "Neurod1", "Gap43", "Tubb3", "Sox11",
                  "Cldn5", "Klf2", "Vwa1", "Egfl7", "Esam", "Cd34", "Kdr", "Eng", "Slc2a1", "Nampt", "Cd93",
                  "Rgs5", "Myl9", "Pdgfrb", "Abcc9", "Mcam", "Cspg4", "Pcdh18")

# Check if all genes are present in the dataset
valid_genes <- gene_markers[gene_markers %in% rownames(filtered_seurat_object)]

# Creating the heatmap for valid markers
DoHeatmap(filtered_seurat_object, features = valid_genes, disp.min = -2, disp.max = 2) + 
  ggtitle("Heatmap of Marker Gene Expression Across Clusters") 

#--------------------------------------------------------------------------------------------------------------------------------------


