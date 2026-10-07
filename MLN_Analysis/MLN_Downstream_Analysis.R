library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# Ignore local RAM throttle 
mem.maxVSize(vsize = Inf)

# Merge samples
#sample_601 <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Pooled/seurat_pool601_merged.rds")
#sample_602 <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Pooled/seurat_pool602_merged.rds")
#sample_603 <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Pooled/seurat_pool603_merged.rds")

# Merge them into one Seurat object
#merged <- merge(
 # sample_601,
 # y = list(sample_602, sample_603),
 # add.cell.ids = c("601", "602", "603"), 
 # project = "MLN"
#)

#saveRDS(merged, "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Pooled/seurat_merged_demux.rds")

# Load merged data
dat <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Pooled/seurat_merged_demux.rds")

head(dat)
#write.csv(table(dat$orig.ident), "cells_per_biological_sample.csv")
#write.csv(rownames(dat), "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed//Pre_QC_genelist.csv")

# Batch corrected object
#dat <- readRDS("BatchCorrected_PostQC_dat")
dim(dat)
# Cell counts before QC
cell_counts <- as.data.frame(table(dat$orig.ident))
#write.csv(cell_counts, "pre_QC_Sample_Cell_Counts.csv")

# Check layers and assays
head(dat)
Layers(dat)
dat@assays
unique(dat@meta.data$orig.ident)


# Look at ADT count matrix
adt_counts <- GetAssayData(dat, assay = "ADT", layer = "counts")
head(adt_counts)
#write.csv(rownames(adt_counts), "ADT_markers_pre_qc.csv")

# RNA Filtering QC column
DefaultAssay(dat) <- "RNA"
dat[["percent.mt"]] <- PercentageFeatureSet(
  dat,
  #pattern = "(?i)^mt-",
  pattern = "mt-",
  assay="RNA"
)

# Ribosmal filtering QC column
dat[["percent.ribo"]] <- PercentageFeatureSet(
  dat,
  pattern = "^Rpl|^Rps",
  assay = "RNA"
)

head(dat)
# Visualize distribution before QC
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
             ncol = 1, pt.size = 0)&
  theme(
    axis.text.x = element_text(size = 35, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 40),
    axis.title = element_text(size = 40),
    plot.title = element_text(size = 44),
    strip.text = element_text(size = 42)
  )
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/MLN_NoQC.png", plot = p,
      # width = 30,
      # height = 50,
      # dpi = 300, limitsize = FALSE)

# Apply QC filter: keep cells with 200–2500 genes and <5% mito
dat <- subset(dat, subset = nFeature_RNA > 200 & 
                #nFeature_RNA < 2500 & 
                percent.mt < 5)

# join layers
dat <- JoinLayers(dat, assay = "RNA")

# Remove lowly expressed genes
DefaultAssay(dat) <- "RNA"
min_pct <- 0.005
min_cells <- ceiling(ncol(dat) * min_pct)
rna_counts <- GetAssayData(dat, assay = "RNA", layer = "counts")
n_cells_expressed <- Matrix::rowSums(rna_counts > 0)
keep_genes <- names(n_cells_expressed)[n_cells_expressed >= min_cells]

cat("Cells:", ncol(dat), "\n")
cat("Min cells required (0.5%):", min_cells, "\n")
cat("RNA genes before:", nrow(rna_counts), "\n")
cat("RNA genes kept:", length(keep_genes), "\n")

dat[["RNA"]] <- subset(dat[["RNA"]], features = keep_genes)
rm(rna_counts, n_cells_expressed, keep_genes)


# Visualize distribution After QC
p_filtered <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
             ncol = 1, pt.size = 0)  & theme(
               axis.text.x = element_text(size = 35, angle = 35, hjust = 1),
               axis.text.y = element_text(size = 40),
               axis.title = element_text(size = 40),
               plot.title = element_text(size = 44),
               strip.text = element_text(size = 42)
             )
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed//Jing_Figures/MLN_MT_AfterQC.png", plot = p_filtered,
      # width = 30,
      # height = 50,
     #  dpi = 300, limitsize = FALSE)


# Just percent.mt
p <- VlnPlot(dat, features = c("percent.mt"), 
             ncol = 1, pt.size = 0) +theme_minimal()+
  theme(axis.text.x = element_text(size = 15, angle = 35, hjust = 1))
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed//Jing_Figures/MLN_MTOnly_AfterQC.png", plot = p,
      # width = 25,
      # height = 10,
       #dpi = 300, limitsize = FALSE)

# Number of Genes before Mitochondrial filter
dat_with_MT <- dim(dat)

# Remove MT genes from RNA assay
rna_feats <- rownames(dat[["RNA"]])
mt_genes <- grep("(?i)^mt-", rna_feats, value = TRUE) 
length(mt_genes)
if (length(mt_genes) > 0) {
  keep_feats <- setdiff(rna_feats, mt_genes)
  # Subset just the RNA assay's features
  dat[["RNA"]] <- subset(dat[["RNA"]], features = keep_feats)
}
dim(dat)

if ("RNA" %in% names(Assays(dat))) {
  cat("RNA features:", nrow(dat[["RNA"]]), 
      "| Cells:", ncol(dat[["RNA"]]), "\n")
}
if ("ADT" %in% names(Assays(dat))) {
  cat("ADT features:", nrow(dat[["ADT"]]), 
      "| Cells:", ncol(dat[["ADT"]]), "\n")
}
dat_with_MT
dim(dat)



# Normalize data and run PCA
DefaultAssay(dat) <- 'RNA'
dat<- NormalizeData(dat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca')

# Visualize ribosomal percent on PCA
FeaturePlot(dat, features = "percent.ribo", reduction = "pca")

# Check correlation with PCs
pc_embeddings <- Embeddings(dat, "pca")
# Correlation with percent ribo
cor(pc_embeddings[,1], dat$percent.ribo)
cor(pc_embeddings[,2], dat$percent.ribo)
cor(pc_embeddings[,3], dat$percent.ribo)
VizDimLoadings(dat, dims = 1:2, reduction = "pca")
# Check correlation with percent mt
cor(pc_embeddings[,1], dat$percent.mt)
cor(pc_embeddings[,2], dat$percent.mt)
cor(pc_embeddings[,3], dat$percent.mt)

# Correlation with library size
cor(pc_embeddings[,1], dat$nCount_RNA)
cor(pc_embeddings[,2], dat$nCount_RNA)
cor(pc_embeddings[,3], dat$nCount_RNA)
cor(pc_embeddings[,4], dat$nCount_RNA)


# Re-normalize after adjusting for ribosomal correlation
dat <- NormalizeData(dat) %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c("percent.ribo", "nCount_RNA")) %>%
  RunPCA()

# Visualize ribosomal percent on PCA
FeaturePlot(dat, features = "percent.ribo", reduction = "pca")

# Check correlation with PCs
pc_embeddings <- Embeddings(dat, "pca")
# Correlation with percent ribo
cor(pc_embeddings[,1], dat$percent.ribo)
cor(pc_embeddings[,2], dat$percent.ribo)
cor(pc_embeddings[,3], dat$percent.ribo)
VizDimLoadings(dat, dims = 1:2, reduction = "pca")
# Correlation with library size
cor(pc_embeddings[,1], dat$nCount_RNA)
cor(pc_embeddings[,2], dat$nCount_RNA)
cor(pc_embeddings[,3], dat$nCount_RNA)
cor(pc_embeddings[,4], dat$nCount_RNA)

# ##########Doublet removal ###########

# re-normalize/scale/calculate PCs


# Process ADT Assay
DefaultAssay(dat) <- 'ADT'
VariableFeatures(dat) <- rownames(dat[["ADT"]])
dat <- NormalizeData(dat, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

# Save Processed/Cleaned Object
#saveRDS(dat, "post_qc_no_cluster.rds)

# View Reductions
Reductions(dat)
names(dat@reductions)
dat[["pca"]]
Embeddings(dat, "pca")[1:5, 1:5]

# Generate initial UMAP
dat <- FindMultiModalNeighbors(
  dat,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:15, 1:18),
  modality.weight.name = "RNA.weight",
  k.nn = 20
)

# UMAP + clustering
dat <- RunUMAP(dat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

# Keep clusters with >= 200 cells
clust_counts <- table(Idents(dat))
keep_clusters <- names(clust_counts)[clust_counts >= 200]
dat <- subset(dat, idents = keep_clusters)

# Plots
head(dat@meta.data)
p_by_sample <- DimPlot(dat, reduction = "wnn.umap", group.by = "pool_id", label = FALSE) +
  ggtitle("WNN UMAP, colored by sample_id") +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/MLN_NoBatch_bySample.png", plot = p_by_sample)

p_by_transcript <- DimPlot(dat, reduction = "wnn.umap", label = TRUE) +
  ggtitle("WNN UMAP, colored by Transcript") +
  theme(plot.title = element_text(hjust = 0.5))  
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed//Jing_Figures/MLN_NoBatch_byTranscript.png", plot = p_by_transcript)

# drop unused factors
dat$seurat_clusters <- droplevels(dat$seurat_clusters)

# Table of post QC cluster and sample cell counts
#write.csv(table(dat$seurat_clusters), "post_QC_Cluster_CellCounts.csv")
#write.csv(table(dat$orig.ident), "post_QC_Cluster_SampleCounts.csv")
#head(dat)
#DimPlot(dat)
#saveRDS(dat, "Original_dat_clustered.rds")







### BATCH CORRECTION WITH HARMONY ###
mem.maxVSize(vsize = Inf)

dat <- readRDS("Original_dat_clustered.rds")

set.seed(4538)

# Harmony on RNA PCs
dat <- RunHarmony(
  object = dat,
  group.by.vars = "orig.ident",
  reduction = "pca",
  assay.use = "RNA",
  dims.use = 1:15,
  reduction.save = "harmony",
  verbose = FALSE
)

# Harmony on ADT PCs
dat <- RunHarmony(
  object = dat,
  group.by.vars = "orig.ident",
  reduction = "apca",
  assay.use = "ADT",
  dims.use = 1:18,
  reduction.save = "aharmony",
  verbose = FALSE
)

# Use harmonized reductions for WNN
dat <- FindMultiModalNeighbors(
  dat,
  reduction.list = list("harmony", "aharmony"),
  dims.list = list(1:15, 1:18),
  modality.weight.name = "RNA.weight",
  k.nn = 20
)

# UMAP + clustering
dat <- RunUMAP(dat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


# Keep clusters with >= 200 cells
clust_counts <- table(Idents(dat))
keep_clusters <- names(clust_counts)[clust_counts >= 200]
dat <- subset(dat, idents = keep_clusters)

# Plots
head(dat@meta.data)
p_by_sample <- DimPlot(dat, reduction = "wnn.umap", group.by = "pool_id", label = FALSE) +
  ggtitle("WNN UMAP After Batch, colored by sample_id") +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/MLN_AfterBatch_bySample.png", plot = p_by_sample)

p_by_transcript <- DimPlot(dat, reduction = "wnn.umap", label = FALSE) +
  ggtitle("WNN UMAP After Batch, colored by Transcript") +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/MLN_AfterBatch_byTranscript.png", plot = p_by_transcript)

#saveRDS(dat, "BatchCorrected_PostQC_dat.rds")

print("Batch Correction Complete")




# Feature Plots
#dat <- readRDS("BatchCorrected_PostQC_dat.rds")
dat <- readRDS("Original_dat_clustered.rds")
p <- DimPlot(dat, label=TRUE)
ggsave(plot=p, filename="annotated_UMAP.png")
# ADT Markers
DefaultAssay(dat) <- "ADT"
#write.csv(rownames(dat), "ADT_Genes.csv")
#ADT_marker_list <- list("MHCI", "CD45","CD3","CD4","CD8a","CD8b",
                        #"TCRb","TCRgd","CD19","CD11c","CD16/32","NK1.1","CD14")
# Not present:"MHCI", "TCRb","TCRgd", "CD16/32", "CD14"

# MHCI 
#TCRb
#TCRgd
#CD16/32
#CD14
ADT_plot <- FeaturePlot(dat,
  features = c("Ms.CD45","Ms.CD3","Ms.CD4","Ms.CD8a","Ms.CD8b",
               "Ms.CD19","Ms.CD11c","Ms.NK.1.1"),
  reduction = "wnn.umap",
  min.cutoff = "q10",
  max.cutoff = "q99",
  cols = c("lightgrey","darkgreen"),
  order=TRUE,
  combine = TRUE
)
ADT_plot
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/ADT_Markers_NoBatch.png", plot = ADT_plot)







# ADT Markers
DefaultAssay(dat) <- "RNA"
#write.csv(rownames(dat), "RNA_Genes.csv")

# MHCI: H2-K1, H2-D1
#TCRb: Trbc1, Trbc2
#TCRgd" Trdc
#CD16/32: Fcgr3, Fcgr2b
#CD14: NOT FOUND

RNA_plot <- FeaturePlot(dat,
                        features = c("Fcgr3", "Fcgr2b", "Trdc", "Trbc1", "Trbc2", "H2-K1", "H2-D1"),
                        reduction = "wnn.umap",
                        min.cutoff = "q10",
                        max.cutoff = "q99",
                        cols = c("lightgrey","blue"),
                        order=TRUE,
                        combine = TRUE
)
RNA_plot
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/RNA_Markers_NoBatch.png", plot = RNA_plot)



cell_counts <- as.data.frame(table(dat$seurat_clusters))
colnames(cell_counts)[colnames(cell_counts) == "Var1"] <- "Cluster"
cell_counts <- cell_counts[cell_counts$Freq > 0,]
write.csv(cell_counts, "post_QC_Cluster_Cell_Counts.csv")

cell_counts <- as.data.frame(table(dat$orig.ident))
head(dat)
dim(dat)
#write.csv(cell_counts, "post_QC_Sample_Cell_Counts.csv")



### DC Subset Feature Plots

DefaultAssay(dat) <- "ADT"
# Not Present: MHCII,CD64,CX3CR1, PDCA-1
ADT_plot <- FeaturePlot(dat,
                        features = c( "Ms.CD11c", "Ms.CD103", "Ms.CD8a",
                                      "Ms.CD8a", "HuMs.CD11b", "Ms.CD103",
                                      "Ms.CD103","Ms.Siglec-H"),
                        reduction = "wnn.umap",
                        min.cutoff = "q10",
                        max.cutoff = "q99",
                        cols = c("lightgrey","darkgreen"),
                        order=TRUE,
                        combine = TRUE
)
ADT_plot
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/DC_subset_ADT_Markers.png", plot = ADT_plot)



DefaultAssay(dat) <- "RNA"
# Not Present: Xcr1, Clec9a, Il23a
RNA_plot <- FeaturePlot(dat,
                        features = c(features <- c(
                          "Itgax", "H2-Ab1","Itgam","Cd8a", "Clec9a", "Itgae", "Irf4", "Klf4",
                          "Itgam", "Cx3cr1", "Tbx21", "Itgae", "Batf3", "Irf8",
                          "Siglech", "Bst2", "Tcf4"
                        )
                        ),
                        reduction = "wnn.umap",
                        min.cutoff = "q10",
                        max.cutoff = "q99",
                        cols = c("lightgrey","blue"),
                        order=TRUE,
                        combine = TRUE
)
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/Jing_Figures/DC_subset_RNA_Markers.png", plot = RNA_plot)


print("Script Complete")
