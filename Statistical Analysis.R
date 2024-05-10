
# Define the gene list
gene_list <- c("S100b", "Gfap", "Aldh1l1", "Slc1a3", "Slc1a2", "Aqp4", "Fabp7",
               "Mbp", "Olig2", "Olig1", "Cntn1", "Mog",
               "Csf1r", "P2ry12", "Tmem119", "Trem2", "Siglech",
               "Neurod1", "Gap43", "Tubb3", "Sox11",
               "Cldn5", "Klf2", "Vwa1", "Egfl7", "Esam", "Cd34", "Kdr", "Eng", "Slc2a1", "Nampt", "Cd93",
               "Rgs5", "Myl9", "Pdgfrb", "Abcc9", "Mcam", "Cspg4", "Pcdh18",
               "Kcnj10", "Kcna1", "Kcnc1", "Kcnq2", "Kcnq3",
               "Scn1a", "Scn3a",
               "Clcn2", "Clcn3", "Clcn4",
               "Gria1", "Gria2", "Gria3", "Gria4",
               "Grin1", "Grin2a", "Grin2b",
               "Htr2c", "Htr2a", "Htr1a", "Tph2", "Slc6a4", "Maoa",
               "Gabbr2", "Gad1", "Gad2", "Slc32a1", "Gabra1", "Gabra2", "Gabra3")

# Perform differential expression for each gene, assuming comparison between two hypothetical clusters, e.g., 0 and 1
results <- lapply(gene_list, function(gene) {
  FindMarkers(filtered_seurat_object, ident.1 = 0, ident.2 = 1, features = gene)
})

# Collect results into a data frame
results_df <- do.call(rbind, results)

# View results
print(results_df)