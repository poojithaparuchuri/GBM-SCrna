# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(devtools)
library(presto)

########################### previous code
saveRDS(filtered_gbm, file = "rds/filt_gbm_plots_singlecell.rds") 


##############################################################################
##############################################################################
######################Endothelial cell analysis###############################
#############################################################################

#extract only endothelial cell data 
#endothelial_cluster <- subset(filtered_gbm, idents = "Endothelial")


#saveRDS(endothelial_cluster, file = "rds/endothelial_cluster.rds")
#endothelial_cluster <- readRDS("endothelial_cluster.rds")

# Check column names in the metadata
colnames(endothelial_cluster@meta.data)

# View the structure of the metadata
head(endothelial_cluster@meta.data)


### analysis 

# Identification of highly variable features
endothelial_cluster <- FindVariableFeatures(endothelial_cluster, selection.method = "vst", nfeatures = 2000)

# Run PCA using only variable features
endothelial_cluster <- RunPCA(endothelial_cluster, features = VariableFeatures(object = endothelial_cluster))

# Cluster cells
endothelial_cluster <- FindNeighbors(endothelial_cluster, dims = 1:10)
endothelial_cluster <- FindClusters(endothelial_cluster, resolution = 0.5)

# UMAP
endothelial_cluster <- RunUMAP(endothelial_cluster, dims = 1:10)

#########################
#Plotting UMAPS
#########################
###############################################################################

#####################Plot UMAP by seurat clusters#####################
endothelial_cluster_umap <- DimPlot(endothelial_cluster, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 2)
endothelial_cluster_umap
ggsave("plots/endothelial_clusters_umap.png", plot = endothelial_cluster_umap, width = 15, height = 8) 

# Plot UMAP by tumor location
EC_location_cluster <- DimPlot(endothelial_cluster, reduction = "umap", group.by = "sample_type", label = TRUE, pt.size = 2, cols = c("blue", "red"))
EC_location_cluster
ggsave("plots/EC_location_cluster.png", plot = EC_location_cluster, width = 15, height = 8) 

#####################UMAP by patient cluster#####################
#copy cluster in case it fails
EC_patient <- endothelial_cluster
# Add a column to distinguish between patient samples
EC_patient$patient <- ifelse(grepl("T", EC_patient$orig.ident), "patient", "peripheral")

# Add another column to identify the specific patient
EC_patient$patient_id <- ifelse(grepl("R1", EC_patient$orig.ident), "1",
                                     ifelse(grepl("R2", EC_patient$orig.ident), "2",
                                            ifelse(grepl("R3", EC_patient$orig.ident), "3",
                                                   ifelse(grepl("R4", EC_patient$orig.ident), "4", NA))))

# Check column names in the metadata
colnames(EC_patient@meta.data)
# View the structure of the metadata
head(EC_patient@meta.data)

# Plot UMAP by sample origin 
EC_patient_cluster <- DimPlot(EC_patient, reduction = "umap", group.by = "patient_id", label = TRUE, pt.size = 2)
EC_patient_cluster
ggsave("plots/EC_patient_cluster.png", plot = EC_patient_cluster, width = 15, height = 8)
#########################################################################

########################## some figure 2 analysis #######################


####Relative contribution of endothelial cells from sample origin type 
# Calculate the percentage of endothelial cells in each sample origin type
# 'sample_type' is the column in the metadata, and 'seurat_clusters' is the cluster assignment column
endothelial_cluster$seurat_clusters <- as.factor(endothelial_cluster$seurat_clusters)

# Create a dataset
plot_data <- data.frame(sample_type = endothelial_cluster$sample_type, cluster = endothelial_cluster$seurat_clusters)

# Stacked + percent

png("plots/location_cluster_stacked_bar_graph.png", width = 400, height = 800)
ggplot(plot_data, aes(fill = sample_type, y = ..count../sum(..count..), x = cluster)) +
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 1)) +  # Adjust limits to 0-1
  labs(title = "Stacked Bar Graph of Sample Type within Seurat Clusters",
       x = "Seurat Clusters", y = "Fraction of Cells") +
  theme_minimal()
dev.off()

# Extract the data as a CSV file
write.csv(plot_data, "data_output/location_cluster_stacked_bar_graph_data.csv", row.names = FALSE)


# Relative contribution of endothelial cells from individual patient.
# 'patient_id' is the column in the metadata, and 'seurat_clusters' is the cluster assignment column
EC_patient$seurat_clusters <- as.factor(EC_patient$seurat_clusters)

# Create a dataset
plot_data <- data.frame(patient_id = EC_patient$patient_id, cluster = EC_patient$seurat_clusters)

# Stacked + percent

png("plots/patient_clusters_stacked_bar_graph.png", width = 400, height = 600)
ggplot(plot_data, aes(fill = patient_id, y = ..count../sum(..count..), x = cluster)) +
  geom_bar(position = "fill", stat = "count") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 1)) +  # Adjust limits to 0-1
  labs(title = "Stacked Bar Graph of Sample Type within Seurat Clusters",
       x = "Seurat Clusters", y = "Fraction of Cells") +
  theme_minimal()
dev.off()

####### Heatmap -- Gene expression levels of top 50 marker genes in different endothelial subclusters.

# Find markers
markers <- FindAllMarkers(endothelial_cluster, only.pos = TRUE, min.pct = 0.25)

# Filter top 50 marker genes
top50_genes <- head(markers$gene, 50)

# DoHeatmap for top 50 genes
DoHeatmap(endothelial_cluster, features = top50_genes, cells = 1:500, size = 4, angle = 90) + NoLegend()
ggsave("heatmap.png", plot = heatmap, , width = 15, height = 8)



#Heatmap showing top 10 enriched GO terms in different endothelial subclusters based on top 50 marker genes.





 