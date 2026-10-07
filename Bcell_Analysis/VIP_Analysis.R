library(harmony)
library(patchwork)
library(dplyr)
library(Seurat)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(purrr)
library(tibble)
library(Matrix)
library(mixOmics)

# Ignore local RAM throttle
mem.maxVSize(vsize = Inf)


# Load merged data

dat <- readRDS("dat_clean_strict.rds")

DefaultAssay(dat) <- "RNA"

# Basic checks
head(dat)
Reductions(dat)
length(unique(dat$seurat_clusters))
table(dat@meta.data$marker_pass_strict)
table(dat$celltype_broad, dat$marker_pass_strict)
table(dat$pool_id)
table(dat$pool_id, dat$celltype_broad)


# output directory
vip_outdir <- "GEX_Results/pathway_Share/VIP_pool_id_results"

if (!dir.exists(vip_outdir)) {
  dir.create(vip_outdir, recursive = TRUE)
}

# Save broad cell-type proportions by pool
celltype_prop <- as.data.frame(
  prop.table(table(dat$pool_id, dat$celltype_broad), margin = 1)
)

colnames(celltype_prop) <- c("pool_id", "celltype_broad", "proportion")

write.csv(
  celltype_prop,
  file = file.path(vip_outdir, "celltype_broad_proportions_by_pool.csv"),
  row.names = FALSE
)

p_celltype_prop <- ggplot(
  celltype_prop,
  aes(x = pool_id, y = proportion, fill = celltype_broad)
) +
  geom_col(position = "fill") +
  theme_bw() +
  labs(
    title = "celltype_broad proportions by pool_id",
    x = "pool_id",
    y = "Proportion",
    fill = "celltype_broad"
  )

ggsave(
  filename = file.path(vip_outdir, "celltype_broad_proportions_by_pool.png"),
  plot = p_celltype_prop,
  width = 8,
  height = 6,
  dpi = 300
)


# Function to run VIP analysis
run_pool_vip <- function(
    seurat_obj,
    analysis_name,
    outdir,
    top_n_genes = 2000,
    top_n_vip = 30,
    min_cells_total = 20,
    min_pb_samples_total = 4,
    min_pb_samples_per_pool = 2
) {
  
  message("============================================================")
  message("Running VIP analysis: ", analysis_name)
  message("============================================================")
  
  analysis_name_safe <- gsub("[^A-Za-z0-9_\\-]", "_", analysis_name)
  analysis_outdir <- file.path(outdir, analysis_name_safe)
  
  if (!dir.exists(analysis_outdir)) {
    dir.create(analysis_outdir, recursive = TRUE)
  }
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  # marker checks
  if (!"pool_id" %in% colnames(seurat_obj@meta.data)) {
    stop("pool_id is not present in metadata.")
  }
  if (!"marker_pass_strict" %in% colnames(seurat_obj@meta.data)) {
    stop("marker_pass_strict is not present in metadata.")
  }
  
  n_cells <- ncol(seurat_obj)
  
  message("Cells in object: ", n_cells)
  message("pool_id counts:")
  print(table(seurat_obj$pool_id))
  
  if ("celltype_broad" %in% colnames(seurat_obj@meta.data)) {
    message("celltype_broad counts:")
    print(table(seurat_obj$celltype_broad))
  }
  
  if (n_cells < min_cells_total) {
    warning("Skipping ", analysis_name, ": fewer than ", min_cells_total, " cells.")
    return(NULL)
  }
  
  # Save cell counts
  cell_counts <- seurat_obj@meta.data %>%
    as.data.frame() %>%
    dplyr::count(pool_id, celltype_broad, name = "n_cells")
  
  write.csv(
    cell_counts,
    file = file.path(analysis_outdir, paste0("cell_counts_", analysis_name_safe, ".csv")),
    row.names = FALSE
  )
  
  # Define pseudobulk replicate unit
  if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$vip_sample_unit <- seurat_obj$orig.ident
  } else if ("demux_id" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$vip_sample_unit <- seurat_obj$demux_id
  } else {
    stop("Need either orig.ident or demux_id to define pseudobulk replicate units.")
  }
  
  seurat_obj$vip_pb_group <- paste(
    seurat_obj$pool_id,
    seurat_obj$vip_sample_unit,
    sep = "__"
  )
  
  pb_design <- seurat_obj@meta.data %>%
    as.data.frame() %>%
    dplyr::distinct(vip_pb_group, pool_id, vip_sample_unit) %>%
    dplyr::arrange(pool_id, vip_sample_unit)
  
  message("Pseudobulk design:")
  print(pb_design)
  print(table(pb_design$pool_id))
  
  write.csv(
    pb_design,
    file = file.path(analysis_outdir, paste0("pseudobulk_design_", analysis_name_safe, ".csv")),
    row.names = FALSE
  )
  
  if (nrow(pb_design) < min_pb_samples_total) {
    warning(
      "Skipping ", analysis_name,
      ": fewer than ", min_pb_samples_total,
      " total pseudobulk samples."
    )
    return(NULL)
  }
  
  if (length(unique(pb_design$pool_id)) < 2) {
    warning("Skipping ", analysis_name, ": fewer than 2 pool_id groups.")
    return(NULL)
  }
  
  pb_per_pool <- table(pb_design$pool_id)
  
  if (any(pb_per_pool < min_pb_samples_per_pool)) {
    warning(
      "Skipping ", analysis_name,
      ": at least one pool_id has fewer than ",
      min_pb_samples_per_pool,
      " pseudobulk replicates."
    )
    return(NULL)
  }
  

  # Create pseudobulk count matrix 
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  groups <- factor(seurat_obj$vip_pb_group)
  names(groups) <- colnames(seurat_obj)
  
  # Sparse design matrix: cells x pseudobulk groups
  group_mat <- Matrix::sparse.model.matrix(~ 0 + groups)
  colnames(group_mat) <- levels(groups)
  
  # Pseudobulk counts: genes x pseudobulk groups
  # avoid sparse -> dense conversion of the full single-cell matrix
  pb_counts <- counts %*% group_mat
  
  # Convert after aggregation only
  pb_counts <- as.matrix(pb_counts)
  
  # Match metadata to pseudobulk columns
  pb_meta <- seurat_obj@meta.data %>%
    as.data.frame() %>%
    dplyr::distinct(vip_pb_group, pool_id, vip_sample_unit) %>%
    dplyr::filter(vip_pb_group %in% colnames(pb_counts)) %>%
    dplyr::arrange(match(vip_pb_group, colnames(pb_counts)))
  
  stopifnot(all(pb_meta$vip_pb_group == colnames(pb_counts)))
  
  write.csv(
    pb_meta,
    file = file.path(analysis_outdir, paste0("pseudobulk_metadata_", analysis_name_safe, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    pb_counts,
    file = file.path(analysis_outdir, paste0("pseudobulk_counts_", analysis_name_safe, ".csv"))
  )
  
  # Filter and normalize genes
  dge <- DGEList(counts = pb_counts)
  keep <- filterByExpr(dge, group = pb_meta$pool_id)
  
  if (sum(keep) < 10) {
    warning("Skipping ", analysis_name, ": fewer than 10 genes passed filterByExpr.")
    return(NULL)
  }
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  logCPM <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Remove zero-variance genes
  gene_var <- apply(logCPM, 1, var)
  logCPM <- logCPM[gene_var > 0, , drop = FALSE]
  
  if (nrow(logCPM) < 10) {
    warning("Skipping ", analysis_name, ": fewer than 10 non-zero-variance genes.")
    return(NULL)
  }
  
  # Keep top variable genes to stabilize PLS-DA
  top_n <- min(top_n_genes, nrow(logCPM))
  top_genes <- names(sort(apply(logCPM, 1, var), decreasing = TRUE))[1:top_n]
  
  X <- t(logCPM[top_genes, , drop = FALSE])
  Y <- factor(pb_meta$pool_id)
  
  message("Pseudobulk samples: ", nrow(X))
  message("Genes used: ", ncol(X))
  message("Pool groups:")
  print(table(Y))
  
  write.csv(
    logCPM,
    file = file.path(analysis_outdir, paste0("pseudobulk_logCPM_", analysis_name_safe, ".csv"))
  )
  
  # Run PLS-DA
  ncomp_use <- min(2, length(levels(Y)) - 1, nrow(X) - 1)
  
  if (ncomp_use < 1) {
    warning("Skipping ", analysis_name, ": not enough components available.")
    return(NULL)
  }
  
  plsda_fit <- plsda(
    X = X,
    Y = Y,
    ncomp = ncomp_use
  )
  
  saveRDS(
    plsda_fit,
    file = file.path(analysis_outdir, paste0("plsda_fit_", analysis_name_safe, ".rds"))
  )
  
  # Extract VIP scores
  vip_mat <- mixOmics::vip(plsda_fit)
  
  vip_df <- as.data.frame(vip_mat) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      VIP_max = apply(dplyr::across(where(is.numeric)), 1, max)
    ) %>%
    dplyr::arrange(desc(VIP_max))
  
  write.csv(
    vip_df,
    file = file.path(analysis_outdir, paste0("VIP_pool_id_PLSDA_genes_", analysis_name_safe, ".csv")),
    row.names = FALSE
  )
  
  # Save loadings
  loadings_df <- as.data.frame(plsda_fit$loadings$X) %>%
    tibble::rownames_to_column("gene")
  
  write.csv(
    loadings_df,
    file = file.path(analysis_outdir, paste0("VIP_pool_id_PLSDA_loadings_", analysis_name_safe, ".csv")),
    row.names = FALSE
  )

  # Plot top VIP genes
  top_vip <- vip_df %>%
    dplyr::slice_head(n = top_n_vip) %>%
    dplyr::mutate(gene = factor(gene, levels = rev(gene)))
  
  p_vip <- ggplot(top_vip, aes(x = gene, y = VIP_max)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste0("Top VIP genes distinguishing pool_id groups: ", analysis_name),
      x = "Gene",
      y = "VIP score"
    )
  
  print(p_vip)
  
  ggsave(
    filename = file.path(
      analysis_outdir,
      paste0("VIP_pool_id_top", top_n_vip, "_genes_", analysis_name_safe, ".png")
    ),
    plot = p_vip,
    width = 8,
    height = 10,
    dpi = 300
  )

  # PLS-DA sample score plot
  scores <- as.data.frame(plsda_fit$variates$X) %>%
    tibble::rownames_to_column("vip_pb_group") %>%
    dplyr::left_join(pb_meta, by = "vip_pb_group")
  
  write.csv(
    scores,
    file = file.path(analysis_outdir, paste0("PLSDA_scores_", analysis_name_safe, ".csv")),
    row.names = FALSE
  )
  
  if (ncomp_use >= 2) {
    
    p_scores <- ggplot(
      scores,
      aes(x = comp1, y = comp2, color = pool_id, label = vip_sample_unit)
    ) +
      geom_point(size = 4) +
      geom_text_repel(size = 3, max.overlaps = 20) +
      theme_bw() +
      labs(
        title = paste0("PLS-DA of pseudobulk samples by pool_id: ", analysis_name),
        x = "Component 1",
        y = "Component 2",
        color = "pool_id"
      )
    
    print(p_scores)
    
    ggsave(
      filename = file.path(
        analysis_outdir,
        paste0("VIP_pool_id_PLSDA_scores_", analysis_name_safe, ".png")
      ),
      plot = p_scores,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  # Top VIP heatmap
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    
    top_vip_genes <- vip_df %>%
      dplyr::slice_head(n = top_n_vip) %>%
      dplyr::pull(gene)
    
    top_expr <- logCPM[top_vip_genes, , drop = FALSE]
    top_expr_scaled <- t(scale(t(top_expr)))
    
    # Annotation must be a regular data.frame with rownames matching columns of heatmap matrix
    annot_col <- pb_meta %>%
      as.data.frame() %>%
      tibble::remove_rownames() %>%
      dplyr::select(vip_pb_group, pool_id, vip_sample_unit)
    
    annot_col <- as.data.frame(annot_col)
    rownames(annot_col) <- annot_col$vip_pb_group
    annot_col$vip_pb_group <- NULL
    
    # Make sure annotation rows are in the same order as heatmap columns
    annot_col <- annot_col[colnames(top_expr_scaled), , drop = FALSE]
    
    pheatmap::pheatmap(
      top_expr_scaled,
      annotation_col = annot_col,
      main = paste0("Top VIP genes: ", analysis_name),
      filename = file.path(
        analysis_outdir,
        paste0("Top_VIP_genes_heatmap_", analysis_name_safe, ".png")
      ),
      width = 8,
      height = 10
    )
    
  } else {
    message("Package pheatmap not installed; skipping top VIP heatmap.")
  }
  message("Finished VIP analysis: ", analysis_name)
  message("Saved to: ", analysis_outdir)
  
  return(
    list(
      analysis_name = analysis_name,
      analysis_name_safe = analysis_name_safe,
      outdir = analysis_outdir,
      pb_meta = pb_meta,
      logCPM = logCPM,
      vip_df = vip_df,
      loadings_df = loadings_df,
      scores = scores,
      plsda_fit = plsda_fit
    )
  )
}


# Run 1: all strict marker-pass cells
dat_vip_all <- subset(
  dat,
  subset = marker_pass_strict == TRUE
)

vip_results <- list()

vip_results[["all_strict"]] <- run_pool_vip(
  seurat_obj = dat_vip_all,
  analysis_name = "all_strict",
  outdir = vip_outdir
)


# Run 2: every cell type in celltype_broad
celltypes_to_run <- sort(unique(as.character(dat$celltype_broad)))
celltypes_to_run <- celltypes_to_run[!is.na(celltypes_to_run)]

message("Cell types to run:")
print(celltypes_to_run)

for (ct in celltypes_to_run) {
  
  message("------------------------------------------------------------")
  message("Preparing cell type: ", ct)
  message("------------------------------------------------------------")
  
  dat_vip_ct <- subset(
    dat,
    subset = marker_pass_strict == TRUE & celltype_broad == ct
  )
  
  message("Cell type table for ", ct, ":")
  print(table(dat_vip_ct$celltype_broad))
  
  message("Pool table for ", ct, ":")
  print(table(dat_vip_ct$pool_id))
  
  ct_safe <- gsub("[^A-Za-z0-9_\\-]", "_", ct)
  
  vip_results[[ct_safe]] <- tryCatch(
    {
      run_pool_vip(
        seurat_obj = dat_vip_ct,
        analysis_name = ct_safe,
        outdir = vip_outdir
      )
    },
    error = function(e) {
      warning("VIP analysis failed for cell type ", ct, ": ", e$message)
      return(NULL)
    }
  )
}

# Combine VIP scores across analyses
valid_results <- vip_results[!sapply(vip_results, is.null)]

if (length(valid_results) > 0) {
  
  vip_combined <- purrr::imap_dfr(
    valid_results,
    function(res, nm) {
      res$vip_df %>%
        dplyr::mutate(analysis = nm) %>%
        dplyr::select(analysis, gene, everything())
    }
  )
  
  write.csv(
    vip_combined,
    file = file.path(vip_outdir, "VIP_scores_all_analyses_long.csv"),
    row.names = FALSE
  )
  
  # Wide table of VIP_max values
  vip_wide <- vip_combined %>%
    dplyr::select(analysis, gene, VIP_max) %>%
    tidyr::pivot_wider(
      names_from = analysis,
      values_from = VIP_max
    )
  
  write.csv(
    vip_wide,
    file = file.path(vip_outdir, "VIP_scores_all_analyses_wide.csv"),
    row.names = FALSE
  )
}


# Shared top VIP genes compared to all_strict
if (!is.null(vip_results[["all_strict"]])) {
  
  all_top30 <- vip_results[["all_strict"]]$vip_df %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::pull(gene)
  
  shared_summary <- purrr::imap_dfr(
    valid_results,
    function(res, nm) {
      
      this_top30 <- res$vip_df %>%
        dplyr::slice_head(n = 30) %>%
        dplyr::pull(gene)
      
      shared <- intersect(all_top30, this_top30)
      
      tibble(
        analysis = nm,
        n_shared_with_all_strict_top30 = length(shared),
        shared_genes = paste(shared, collapse = ";")
      )
    }
  )
  
  write.csv(
    shared_summary,
    file = file.path(vip_outdir, "Shared_top30_VIP_genes_vs_all_strict.csv"),
    row.names = FALSE
  )
}


# Save run summary
run_summary <- purrr::imap_dfr(
  vip_results,
  function(res, nm) {
    
    if (is.null(res)) {
      return(
        tibble(
          analysis = nm,
          status = "failed_or_skipped",
          n_pseudobulk_samples = NA_integer_,
          n_genes = NA_integer_,
          top_gene = NA_character_,
          top_VIP = NA_real_
        )
      )
    }
    
    tibble(
      analysis = nm,
      status = "completed",
      n_pseudobulk_samples = nrow(res$pb_meta),
      n_genes = nrow(res$vip_df),
      top_gene = res$vip_df$gene[1],
      top_VIP = res$vip_df$VIP_max[1]
    )
  }
)

write.csv(
  run_summary,
  file = file.path(vip_outdir, "VIP_run_summary.csv"),
  row.names = FALSE
)

message("============================================================")
message("All VIP analyses complete.")
message("Results saved to: ", vip_outdir)
message("============================================================")







