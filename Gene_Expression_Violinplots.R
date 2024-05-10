# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("S100b", "Gfap", "Aldh1l1", "Slc1a3", "Slc1a2", "Aqp4", "Fabp7") # vasculature associated astrocytes 
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")


# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Mbp", "Olig2", "Olig1", "Cntn1", "Mog") # oligodendrocytes 
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")



# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Csf1r", "P2ry12", "Tmem119", "Trem2", "Siglech") #Microglia 
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")



# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Neurod1", "Gap43", "Tubb3", "Sox11") # neuron 
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")



# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Cldn5", "Klf2", "Vwa1", "Egfl7", "Esam", "Cd34", "Kdr", "Eng", "Slc2a1", "Nampt", "Cd93") # Endothelial
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")




# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Rgs5", "Myl9", "Pdgfrb", "Abcc9", "Mcam", "Cspg4", "Pcdh18") # Pericyte cell
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")




# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Kcnj10", "Kcna1", "Kcnc1", "Kcnq2", "Kcnq3") #potassium channels
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")






# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Clcn2", "Clcn3", "Clcn4") #chloride channel
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")




# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Gria1", "Gria2", "Gria3", "Gria4") #AMPA receptors
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")




# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Grin1", "Grin2a", "Grin2b") #NMDA receptors 
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")





# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Htr2c", "Htr2a", "Htr1a", "Tph2", "Slc6a4", "Maoa") # serotonin receptors genes  
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")





# FeaturePlot to visualize spatial gene expression----- protoplasmic astrocytes 
feature_genes <- c("Gabbr2", "Gad1", "Gad2", "Slc32a1", "Gabra1", "Gabra2", "Gabra3") # GABA receptor genes 
FeaturePlot(filtered_seurat_object, features = feature_genes, min.cutoff = "q10", max.cutoff = "q90") + 
  ggtitle("Expression of Marker Genes across Clusters")
# VlnPlot to visualize distribution of gene expression across clusters
VlnPlot(filtered_seurat_object, features = feature_genes, pt.size = 0) + 
  ggtitle("Distribution of Marker Gene Expression")

