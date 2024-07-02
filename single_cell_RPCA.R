#RPCA 
# Single-cell RNA-seq analysis - QC
#https://hbctraining.github.io/scRNA-seq_online/lessons/03_SC_quality_control-setup.html

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(devtools)
library(presto)

# Load the single-cell seq data from all patients
# T = tumor core
# N = peripheral
R1_T.data <- Read10X(data.dir = "raw_counts_matrix/R1_T/")
R1_N.data <- Read10X(data.dir = "raw_counts_matrix/R1_N/")
R2_T.data <- Read10X(data.dir = "raw_counts_matrix/R2_T/")
R2_N.data <- Read10X(data.dir = "raw_counts_matrix/R2_N/")
R3_T.data <- Read10X(data.dir = "raw_counts_matrix/R3_T/")
R3_N.data <- Read10X(data.dir = "raw_counts_matrix/R3_N/")
R4_T.data <- Read10X(data.dir = "raw_counts_matrix/R4_T/")
R4_N.data <- Read10X(data.dir = "raw_counts_matrix/R4_N/")

# Create Seurat objects for each dataset
T1 <- CreateSeuratObject(counts = R1_T.data, project = "R1_T", min.features = 200)
N1 <- CreateSeuratObject(counts = R1_N.data, project = "R1_N", min.features = 200)
T2 <- CreateSeuratObject(counts = R2_T.data, project = "R2_T", min.features = 200)
N2 <- CreateSeuratObject(counts = R2_N.data, project = "R2_N", min.features = 200)
T3 <- CreateSeuratObject(counts = R3_T.data, project = "R3_T", min.features = 200)
N3 <- CreateSeuratObject(counts = R3_N.data, project = "R3_N", min.features = 200)
T4 <- CreateSeuratObject(counts = R4_T.data, project = "R4_T", min.features = 200)
N4 <- CreateSeuratObject(counts = R4_N.data, project = "R4_N", min.features = 200)

# Normalize data 
T1 <- NormalizeData(T1, normalization.method = "LogNormalize", scale.factor = 10000)
N1 <- NormalizeData(N1, normalization.method = "LogNormalize", scale.factor = 10000)
T2 <- NormalizeData(T2, normalization.method = "LogNormalize", scale.factor = 10000)
N2 <- NormalizeData(N2, normalization.method = "LogNormalize", scale.factor = 10000)
T3 <- NormalizeData(T3, normalization.method = "LogNormalize", scale.factor = 10000)
N3 <- NormalizeData(N3, normalization.method = "LogNormalize", scale.factor = 10000)
T4 <- NormalizeData(T4, normalization.method = "LogNormalize", scale.factor = 10000)
N4 <- NormalizeData(N4, normalization.method = "LogNormalize", scale.factor = 10000)

# Create a list of Seurat objects for each dataset
tumor_list <- list(T1, T2, T3, T4)
normal_list <- list(N1, N2, N3, N4)

# Run PCA on each dataset
tumor_list <- lapply(tumor_list, function(seurat_obj) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
  return(seurat_obj)
})


normal_list <- lapply(normal_list, function(seurat_obj) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
  return(seurat_obj)
})


# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = c(tumor_list, normal_list), dims = 1:20)

# Integrate datasets
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20)

#save
saveRDS(integrated_data, file = "rds/integrated_data.rds")

####################Run the standard workflow for visualization and clustering 
#scale data
integrated_data <- ScaleData(integrated_data, verbose = FALSE)
# Run PCA on integrated data
integrated_data <- RunPCA(integrated_data, features = VariableFeatures(integrated_data), verbose = FALSE)
# Run UMAP 
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)
#save plots for now
plot_umap_rpc <-DimPlot(integrated_data, reduction = "umap", pt.size = 0.1, raster = FALSE)
ggsave("plot_umap_rpc.png", plot = plot_umap_rpc, width = 15, height = 15)
#find neighbors 
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
#find clusters
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Check column names in the metadata
colnames(integrated_data@meta.data)

# View the structure of the metadata
head(integrated_data@meta.data)

##############now add a column to distinguis between tumor and peripheral

# Add a column to distinguish between tumor and peripheral samples
integrated_data$sample_type <- ifelse(integrated_data$orig.ident %in% c("R1_T", "R2_T", "R3_T", "R4_T"), "tumor", "peripheral")
# Check column names in the metadata
colnames(integrated_data@meta.data)
# View the structure of the metadata
head(integrated_data@meta.data)


# Run UMAP on integrated data
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)

# Plot UMAP with colors based on sample type
DimPlot(integrated_data, reduction = "umap", group.by = "sample_type", label = TRUE)


#visualization 
p1 <- DimPlot(integrated_data, reduction = "umap", group.by = "sample_type")
p2 <- DimPlot(integrated_data, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 

view(p1)


# Find markers enriched in clusters, show only positive markers
markers_df <- FindAllMarkers(integrated_data, only.pos = TRUE)
markers_df %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Assuming you've already obtained 'markers_df'
# Filter markers based on avg_log2FC and select specific clusters
selected_clusters <- c(0, 1, 2, 8, 10, 12, 14, 15, 18)
filtered_markers_df <- markers_df %>%
  filter(cluster %in% selected_clusters, avg_log2FC > 1)

# Use the top marker genes for the selected clusters
marker_genes <- filtered_markers_df$gene

#filter out the clusters that are not needed into a new seurat object
filtered_gbm <- subset(integrated_data, idents = selected_clusters)

# Check the order of cluster names
cluster_order <- levels(Idents(filtered_gbm))
print(cluster_order)

##################### Rename cluster identities in Seurat object
# Assuming you have a vector of new cluster names called new_cluster_names
new_cluster_names <- c("Macrophage", "Microglia", "Dendritic Cell", "Neutrophil", "B Cell", "Glia/Neuronal cell", "Mural Cell", "Endothelial", "T cell")

# Assign new cluster names and associate them with correct levels
names(new_cluster_names) <- levels(filtered_gbm)
filtered_gbm <- RenameIdents(filtered_gbm, new_cluster_names)

cluster_order <- levels(Idents(filtered_gbm))
print(cluster_order)


# Plot UMAP with new cluster names
gbm_umap <- DimPlot(filtered_gbm, reduction = "umap", pt.size = 0.05) + NoLegend()
ggsave("umap_plot_all_cells.png", plot = gbm_umap, width = 15, height = 12)

####################
#find all marker genes again with filtered clusters seurat to create a dot plot
####################
# Find markers enriched in clusters, show only positive markers
markers_df_filt <- FindAllMarkers(filtered_gbm, only.pos = TRUE)
markers_df_filt %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#filt top 20 markers
markers_df_filt %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

#set 20 markers to new markers variable 
markers <- top20$gene
markers

######the seurat object already has the new cluster names, so we do not need to assign the names, just plot 


# Assuming 'merged_gbm' has been updated with new cluster names
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
# Assuming 'merged_gbm' has been updated with new cluster names
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
#singlecell_dotplot <- DotPlot(filtered_gbm, features = marker_genes, group.by = 'cluster') + RotatedAxis()

# Define all genes in one list
all_genes <- c("CLDN5", "VWF", "CD34", "APOC1", "CD163", "F13A1",
               "IL1R2", "CD3D", "CD3E", "GZMK", "IGHG1", "IGHG3", "CD79A",
               "HLA-DQA1", "FABP7", "PTPRZ1", "RGS5", "PDGFRB", "NOTCH3")

# Find the top 5 markers in each cluster
markers_df_filt %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5_markers

# Extract the gene names
top5_genes <- top5_markers$gene

# DotPlot with filtered data and top 5 genes per cluster

singlecell_dotplot <- DotPlot(filtered_gbm, features = top5_genes, cluster.idents = TRUE, cols = "RdYlBu", dot.scale = 18) +
  RotatedAxis()

singlecell_dotplot
ggsave("singlecell_dotplot_top5.png", plot = singlecell_dotplot, width = 25, height = 14)



#save to csv
write.csv(top5_markers, "data_output/top5genes_percluster.csv", row.names = FALSE)
write.csv(top20, "data_output/top20genes_percluster.csv", row.names = FALSE)
write.csv(top50, "data_output/top50genes_percluster.csv", row.names = FALSE)


########################### save project before exiting 

saveRDS(filtered_gbm, file = "rds/filt_gbm_plots_singlecell.rds") 


########################extract only endothelial cell data 
endothelial_cluster <- subset(filtered_gbm, idents = "Endothelial")
