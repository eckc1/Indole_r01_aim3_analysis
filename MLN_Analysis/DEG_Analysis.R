library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(Matrix)

dat <- readRDS("../Original_dat_clustered.rds")
DefaultAssay(dat) <- "RNA"

cluster_use <- "36"

out_dir <- "/Users/eckco/Desktop/Kuhn_Lab/Jing_Data/MLN/Demultiplexed/CL36_Analysis/cluster36_pseudobulk_volcano"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Subset to cluster 36
Idents(dat) <- "seurat_clusters"
cl36 <- subset(dat, idents = cluster_use)
dim(cl36)
table(cl36$orig.ident)
table(cl36$pool_id)
sort(table(cl36$orig.ident))

# Raw counts
counts <- GetAssayData(cl36, assay = "RNA", layer = "counts")


# Build pseudobulk counts
# one pseudobulk column per orig.ident
pb_ids <- cl36$orig.ident
pb_levels <- unique(pb_ids)

pb_counts <- sapply(pb_levels, function(x) {
  Matrix::rowSums(counts[, pb_ids == x, drop = FALSE])
})

pb_counts <- as.matrix(pb_counts)
mode(pb_counts) <- "integer"

# Metadata for DESeq2
meta_pb <- data.frame(
  pb_id = colnames(pb_counts),
  stringsAsFactors = FALSE
)

meta_pb$treatment <- sub("^pool([0-9]+)_.*$", "\\1", meta_pb$pb_id)
meta_pb$replicate <- sub("^pool[0-9]+_", "", meta_pb$pb_id)

rownames(meta_pb) <- meta_pb$pb_id

# Volcano plot function
make_volcano <- function(res_df, title, out_file,
                         lfc_thr = 0.5,
                         padj_thr = 0.05,
                         top_n = 10) {
  
  res_df$gene <- rownames(res_df)
  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
  res_df$neglog10_padj <- -log10(res_df$padj + 1e-300)
  
  res_df <- res_df %>%
    mutate(
      category = case_when(
        padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Down",
        TRUE ~ "NotSig"
      )
    )
  
  top_up <- res_df %>%
    filter(category == "Up") %>%
    arrange(padj, desc(log2FoldChange)) %>%
    slice_head(n = top_n)
  
  top_down <- res_df %>%
    filter(category == "Down") %>%
    arrange(padj, log2FoldChange) %>%
    slice_head(n = top_n)
  
  top_genes <- bind_rows(top_up, top_down)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = category), alpha = 0.75, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey70")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2
    ) +
    theme_classic() +
    labs(
      title = title,
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    )
  
  ggsave(out_file, p, width = 8, height = 6, dpi = 300)
  p
}


# Pairwise DESeq2
run_pb_deseq <- function(count_mat, meta_df, group1, group2, out_prefix, out_dir) {
  
  keep_samples <- meta_df$treatment %in% c(group1, group2)
  
  sub_meta <- meta_df[keep_samples, , drop = FALSE]
  sub_counts <- count_mat[, rownames(sub_meta), drop = FALSE]
  
  sub_meta$treatment <- factor(sub_meta$treatment, levels = c(group1, group2))
  
  # filter low-count genes
  keep_genes <- rowSums(sub_counts) >= 10
  sub_counts <- sub_counts[keep_genes, , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = sub_counts,
    colData = sub_meta,
    design = ~ treatment
  )
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("treatment", group2, group1))
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj), ]
  
  write.csv(
    res_df,
    file.path(out_dir, paste0(out_prefix, "_DESeq2_results.csv")),
    row.names = TRUE
  )
  
  p <- make_volcano(
    res_df,
    title = paste("Cluster 36 pseudobulk:", group2, "vs", group1),
    out_file = file.path(out_dir, paste0(out_prefix, "_volcano.png")),
    lfc_thr = 0.5,
    padj_thr = 0.05,
    top_n = 10
  )
  
  list(res = res_df, plot = p)
}


# Run comparisons
res_601_602 <- run_pb_deseq(
  pb_counts, meta_pb,
  group1 = "601", group2 = "602",
  out_prefix = "cluster36_602_vs_601",
  out_dir = out_dir
)

res_601_603 <- run_pb_deseq(
  pb_counts, meta_pb,
  group1 = "601", group2 = "603",
  out_prefix = "cluster36_603_vs_601",
  out_dir = out_dir
)

res_602_603 <- run_pb_deseq(
  pb_counts, meta_pb,
  group1 = "602", group2 = "603",
  out_prefix = "cluster36_603_vs_602",
  out_dir = out_dir
)

# Show plots in R
res_601_602$plot
res_601_603$plot
res_602_603$plot

head(dat)
dim(dat)




