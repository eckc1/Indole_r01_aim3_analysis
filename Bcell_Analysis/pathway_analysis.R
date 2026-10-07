# PSEUDOBULK DESEQ2 DEG + PATHWAY ANALYSIS
# RUNS INDEPENDENTLY FOR EVERY celltype_broad

library(Seurat)
library(Matrix)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(stringr)
library(patchwork)


# Load Data

dat <- readRDS("dat_clean_strict.rds")

# Basic checks
head(dat)
Reductions(dat)
length(unique(dat$seurat_clusters))
table(dat@meta.data$marker_pass_strict)
table(dat$celltype_broad, dat$marker_pass_strict)
table(dat$pool_id)
table(dat$demux_id)
table(dat$orig.ident)
table(dat$pool_id, dat$celltype_broad)



DefaultAssay(dat) <- "RNA"


celltype_col <- "celltype_broad"

# Treatment / condition
group_col <- "pool_id"


# Remove  601 (TS + Veh) group
dat <- subset(
  dat,
  subset = pool_id != "601"
)

# Verify
table(dat$pool_id)

# Adjust treatment group labels

group_label_map <- c(
  #"601" = "TS + Veh",
  "602" = "TL + Veh",
  "603" = "TL + Ind"
)

get_group_label <- function(x) {
  
  x <- as.character(x)
  
  labelled <- group_label_map[x]
  
  labelled <- ifelse(
    is.na(labelled),
    x,
    labelled
  )
  
  return(labelled)
}

# Create new metadata column with readable treatment labels
dat$group_label <- unname(
  get_group_label(dat$pool_id)
)

# UMAPs: USE wnn reduction


# Check labels
table(dat$pool_id, dat$group_label)

# UMAPs: USE wnn reduction
DimPlot(
  dat,
  reduction = "wnn.umap",
  group.by = "group_label"
)

DimPlot(
  dat,
  reduction = "wnn.umap",
  label = TRUE
)

DimPlot(
  dat,
  reduction = "wnn.umap",
  group.by = "celltype_broad"
)

#########################################
# Biological replicate / mouse
replicate_col <- "demux_id"


# DEG thresholds
raw_p_cutoff <- 0.05
fdr_cutoff   <- 0.05
logfc_cutoff <- 0.5


# Root output directory
root_output_dir <- "GEX_Results/pathway_Share"

dir.create(
  root_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Check metadata
required_metadata <- c(
  celltype_col,
  group_col,
  replicate_col
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(dat@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}


# Display Cell Types
cell_types <- sort(
  unique(
    as.character(
      dat@meta.data[[celltype_col]]
    )
  )
)

cell_types <- cell_types[
  !is.na(cell_types) &
    cell_types != ""
]


cat("\n")
cat("CELL TYPES FOUND\n")
print(cell_types)
cat("\nNumber of cell types:", length(cell_types), "\n")


# Double check filenames

safe_name <- function(x) {
  
  x <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  x <- gsub(
    "^_|_$",
    "",
    x
  )
  
  return(x)
}



# Pathway Databases

species <- "Mus musculus"
db_species <- "HS"


cat("\n")
cat("========================================\n")
cat("Loading MSIGDB Databases\n")
cat("========================================\n")

# Hallmark

hallmark_sets <- msigdbr(
  db_species = db_species,
  species = species,
  collection = "H"
) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()



# KEGG

kegg_sets <- msigdbr(
  db_species = db_species,
  species = species,
  collection = "C2",
  subcollection = "CP:KEGG_MEDICUS"
) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()


# GO Biological Process

go_bp_sets <- msigdbr(
  db_species = db_species,
  species = species,
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()



# GO Cellular Component

go_cc_sets <- msigdbr(
  db_species = db_species,
  species = species,
  collection = "C5",
  subcollection = "GO:CC"
) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()



# GO Molecular Function

go_mf_sets <- msigdbr(
  db_species = db_species,
  species = species,
  collection = "C5",
  subcollection = "GO:MF"
) %>%
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()


go_sets <- bind_rows(
  go_bp_sets,
  go_cc_sets,
  go_mf_sets
) %>%
  distinct()


cat("\nHallmark genes:",
    length(unique(hallmark_sets$gene_symbol)),
    "\n")

cat("KEGG genes:",
    length(unique(kegg_sets$gene_symbol)),
    "\n")

cat("GO genes:",
    length(unique(go_sets$gene_symbol)),
    "\n")


# Clean Pathway Names

clean_pathway_names <- function(x) {
  
  x <- gsub(
    "^HALLMARK_",
    "",
    x
  )
  
  x <- gsub(
    "^KEGG_MEDICUS_REFERENCE_",
    "",
    x
  )
  
  x <- gsub(
    "^KEGG_MEDICUS_PATHOGEN_",
    "",
    x
  )
  
  x <- gsub(
    "^KEGG_MEDICUS_ENV_FACTOR_",
    "",
    x
  )
  
  x <- gsub(
    "^KEGG_",
    "",
    x
  )
  
  x <- gsub(
    "^GOBP_",
    "",
    x
  )
  
  x <- gsub(
    "^GOCC_",
    "",
    x
  )
  
  x <- gsub(
    "^GOMF_",
    "",
    x
  )
  
  x <- gsub(
    "_",
    " ",
    x
  )
  
  x <- stringr::str_wrap(
    x,
    width = 42
  )
  
  return(x)
}



# Raw P value volcano

make_raw_p_volcano <- function(
    results_df,
    comparison_name,
    output_dir
) {
  
  df <- results_df %>%
    mutate(
      
      plot_pvalue = ifelse(
        is.na(pvalue),
        1,
        pvalue
      ),
      
      plot_pvalue = pmax(
        plot_pvalue,
        .Machine$double.xmin
      ),
      
      neg_log10_pvalue =
        -log10(plot_pvalue),
      
      significance = case_when(
        
        !is.na(pvalue) &
          pvalue < raw_p_cutoff &
          log2FoldChange >= logfc_cutoff ~
          "Upregulated",
        
        !is.na(pvalue) &
          pvalue < raw_p_cutoff &
          log2FoldChange <= -logfc_cutoff ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      )
    )
  
  
  # Top 10 upregulated by raw P-value
  top_up <- df %>%
    filter(
      significance == "Upregulated"
    ) %>%
    arrange(
      pvalue,
      desc(log2FoldChange)
    ) %>%
    slice_head(
      n = 10
    )
  
  
  # Top 10 downregulated by raw P-value
  top_down <- df %>%
    filter(
      significance == "Downregulated"
    ) %>%
    arrange(
      pvalue,
      log2FoldChange
    ) %>%
    slice_head(
      n = 10
    )
  
  
  top_labels <- bind_rows(
    top_up,
    top_down
  )
  
  
  p <- ggplot(
    df,
    aes(
      x = log2FoldChange,
      y = neg_log10_pvalue,
      color = significance
    )
  ) +
    
    geom_point(
      alpha = 0.70,
      size = 1.5
    ) +
    
    geom_vline(
      xintercept = c(
        -logfc_cutoff,
        logfc_cutoff
      ),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    
    geom_hline(
      yintercept =
        -log10(raw_p_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(
        label = gene
      ),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "gray70"
      )
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(
        "DESeq2 Volcano: ",
        comparison_name
      ),
      
      subtitle = paste0(
        "Raw P < ",
        raw_p_cutoff,
        " and |log2FC| >= ",
        logfc_cutoff,
        "; labels = top 10 up/down by raw P-value"
      ),
      
      x = "log2 fold-change",
      
      y = "-log10(raw P-value)",
      
      color = "Result"
    )
  
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Volcano_",
        safe_name(comparison_name),
        "_RawP.png"
      )
    ),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  
  return(p)
}



# FDR Volcano

make_fdr_volcano <- function(
    results_df,
    comparison_name,
    output_dir
) {
  
  df <- results_df %>%
    mutate(
      
      plot_padj = ifelse(
        is.na(padj),
        1,
        padj
      ),
      
      plot_padj = pmax(
        plot_padj,
        .Machine$double.xmin
      ),
      
      neg_log10_FDR =
        -log10(plot_padj),
      
      significance = case_when(
        
        !is.na(padj) &
          padj < fdr_cutoff &
          log2FoldChange >= logfc_cutoff ~
          "Upregulated",
        
        !is.na(padj) &
          padj < fdr_cutoff &
          log2FoldChange <= -logfc_cutoff ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      )
    )
  
  
  # Top 10 upregulated by FDR
  top_up <- df %>%
    filter(
      significance == "Upregulated"
    ) %>%
    arrange(
      padj,
      desc(log2FoldChange)
    ) %>%
    slice_head(
      n = 10
    )
  
  
  # Top 10 downregulated by FDR
  top_down <- df %>%
    filter(
      significance == "Downregulated"
    ) %>%
    arrange(
      padj,
      log2FoldChange
    ) %>%
    slice_head(
      n = 10
    )
  
  
  top_labels <- bind_rows(
    top_up,
    top_down
  )
  
  
  p <- ggplot(
    df,
    aes(
      x = log2FoldChange,
      y = neg_log10_FDR,
      color = significance
    )
  ) +
    
    geom_point(
      alpha = 0.70,
      size = 1.5
    ) +
    
    geom_vline(
      xintercept = c(
        -logfc_cutoff,
        logfc_cutoff
      ),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    
    geom_hline(
      yintercept =
        -log10(fdr_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(
        label = gene
      ),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "gray70"
      )
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(
        "DESeq2 Volcano: ",
        comparison_name
      ),
      
      subtitle = paste0(
        "FDR < ",
        fdr_cutoff,
        " and |log2FC| >= ",
        logfc_cutoff,
        "; labels = top 10 up/down by FDR"
      ),
      
      x = "log2 fold-change",
      
      y = "-log10(FDR)",
      
      color = "Result"
    )
  
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Volcano_",
        safe_name(comparison_name),
        "_FDR0.05.png"
      )
    ),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  
  return(p)
}



# GSEA Barplot

make_gsea_barplot <- function(
    results_df,
    database_name,
    comparison_name,
    result_dir,
    fdr_cutoff = 0.25,
    top_n = 10
) {
  
  plot_df <- results_df %>%
    filter(
      !is.na(NES),
      !is.na(p.adjust),
      p.adjust < fdr_cutoff
    )
  
  
  if (nrow(plot_df) == 0) {
    
    cat(
      "No significant pathways available for GSEA bar plot.\n"
    )
    
    return(NULL)
  }
  
  
  positive_df <- plot_df %>%
    filter(
      NES > 0
    ) %>%
    arrange(
      desc(NES)
    ) %>%
    slice_head(
      n = top_n
    )
  
  
  negative_df <- plot_df %>%
    filter(
      NES < 0
    ) %>%
    arrange(
      NES
    ) %>%
    slice_head(
      n = top_n
    )
  
  
  plot_df <- bind_rows(
    negative_df,
    positive_df
  )
  
  
  if (nrow(plot_df) == 0) {
    
    return(NULL)
  }
  
  
  plot_df <- plot_df %>%
    mutate(
      
      Direction = ifelse(
        NES > 0,
        paste0("Enriched in ", group2_label),
        paste0("Enriched in ", group1_label)
      ),
      
      Pathway = clean_pathway_names(
        Description
      )
    ) %>%
    arrange(
      NES
    ) %>%
    mutate(
      Pathway = factor(
        Pathway,
        levels = unique(Pathway)
      )
    )
  
  fill_colors <- setNames(
    c("#E41A1C", "#377EB8"),
    c(
      paste0("Enriched in ", group2_label),
      paste0("Enriched in ", group1_label)
    )
  )
  p <- ggplot(
    plot_df,
    aes(
      x = NES,
      y = Pathway,
      fill = Direction
    )
  ) +
    
    geom_col(
      width = 0.7
    ) +
    scale_fill_manual(
      values = fill_colors,
      drop = FALSE
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    
    labs(
      title = paste0(
        database_name,
        " GSEA"
      ),
      subtitle = paste0(
        comparison_name,
        "\nFDR < ",
        fdr_cutoff
      ),
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = NULL
    ) +
    
    theme_classic(
      base_size = 9
    ) +
    
    theme(
      legend.position = "top",
      
      legend.text = element_text(
        size = 7
      ),
      
      plot.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        size = 8
      ),
      
      axis.title.x = element_text(
        size = 9
      ),
      
      axis.text.x = element_text(
        size = 7
      ),
      
      axis.text.y = element_text(
        size = 6.5,
        lineheight = 0.9
      ),
      
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
    )
  
  
  plot_height <- max(
    7,
    0.55 * nrow(plot_df) + 2
  )
  
  
  cutoff_label <- gsub(
    "\\.",
    "_",
    sprintf("%.2f", fdr_cutoff)
  )
  
  ggsave(
    filename = file.path(
      result_dir,
      paste0(
        database_name,
        "_GSEA_FDR_",
        cutoff_label,
        "_positive_negative_barplot.png"
      )
    ),
    plot = p,
    width = 12,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
  
  
  write.csv(
    plot_df,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_GSEA_FDR_",
        cutoff_label,
        "_positive_negative_barplot_data.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  return(p)
}



# ORA Barplot

make_ora_barplot <- function(
    results_df,
    database_name,
    comparison_name,
    direction,
    result_dir,
    fdr_cutoff = 0.05,
    top_n = 20
) {
  
  plot_df <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < fdr_cutoff
    )
  
  
  if (nrow(plot_df) == 0) {
    
    cat(
      "No significant pathways available for ORA bar plot.\n"
    )
    
    return(NULL)
  }
  
  
  plot_df <- plot_df %>%
    separate(
      GeneRatio,
      into = c(
        "GeneRatio_num",
        "GeneRatio_den"
      ),
      sep = "/",
      remove = FALSE,
      convert = TRUE
    ) %>%
    
    mutate(
      
      GeneRatio_numeric =
        GeneRatio_num / GeneRatio_den,
      
      Pathway = clean_pathway_names(
        Description
      ),
      
      minus_log10_FDR =
        -log10(p.adjust)
    ) %>%
    
    arrange(
      desc(GeneRatio_numeric)
    ) %>%
    
    slice_head(
      n = top_n
    ) %>%
    
    arrange(
      GeneRatio_numeric
    ) %>%
    
    mutate(
      Pathway = factor(
        Pathway,
        levels = unique(Pathway)
      )
    )
  
  
  p <- ggplot(
    plot_df,
    aes(
      x = GeneRatio_numeric,
      y = Pathway
    )
  ) +
    
    geom_col(
      width = 0.65
    ) +
    
    labs(
      title = paste0(
        database_name,
        " ORA"
      ),
      
      subtitle = paste0(
        comparison_name,
        " - ",
        direction,
        "\nFDR < ",
        fdr_cutoff
      ),
      
      x = "Gene Ratio",
      y = NULL
    ) +
    
    theme_classic(
      base_size = 9
    ) +
    
    theme(
      
      plot.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        size = 8
      ),
      
      axis.title.x = element_text(
        size = 9
      ),
      
      axis.text.x = element_text(
        size = 7
      ),
      
      axis.text.y = element_text(
        size = 6.5,
        lineheight = 0.9
      ),
      
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
    )
  
  
  plot_height <- max(
    8,
    0.65 * nrow(plot_df) + 2
  )
  
  
  cutoff_label <- gsub(
    "\\.",
    "_",
    sprintf("%.2f", fdr_cutoff)
  )
  
  ggsave(
    filename = file.path(
      result_dir,
      paste0(
        database_name,
        "_",
        direction,
        "_ORA_FDR_",
        cutoff_label,
        "_barplot.png"
      )
    ),
    plot = p,
    width = 13,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
  
  
  write.csv(
    plot_df,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_",
        direction,
        "_ORA_FDR_",
        cutoff_label,
        "_barplot_data.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  return(p)
}



# GSEA Function

run_gsea <- function(
    deg_df,
    gene_sets,
    database_name,
    comparison_name,
    output_dir
) {
  
  cat("\n")
  cat("========================================\n")
  cat("GSEA:", database_name, "\n")
  cat("Comparison:", comparison_name, "\n")
  cat("========================================\n")
  
  
  # DESeq2 Wald statistic is used for ranking.
  # This incorporates both effect size and uncertainty and is
  # generally preferable to ranking only by log2FoldChange.
  
  ranked_df <- deg_df %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(stat),
      is.finite(stat)
    ) %>%
    
    group_by(
      gene
    ) %>%
    
    slice_max(
      order_by = abs(stat),
      n = 1,
      with_ties = FALSE
    ) %>%
    
    ungroup() %>%
    
    arrange(
      desc(stat)
    )
  
  
  gene_list <- ranked_df$stat
  
  names(gene_list) <- ranked_df$gene
  
  
  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )
  
  
  pathway_genes <- unique(
    gene_sets$gene_symbol
  )
  
  
  overlapping_genes <- intersect(
    names(gene_list),
    pathway_genes
  )
  
  
  cat(
    "Genes in ranked DESeq2 table:",
    length(gene_list),
    "\n"
  )
  
  
  cat(
    "Genes matching",
    database_name,
    ":",
    length(overlapping_genes),
    "\n"
  )
  
  
  if (length(gene_list) > 0) {
    
    cat(
      "Percent mapped:",
      round(
        100 *
          length(overlapping_genes) /
          length(gene_list),
        2
      ),
      "%\n"
    )
  }
  
  
  if (length(overlapping_genes) < 10) {
    
    cat(
      "Too few genes overlap with pathway database. Skipping.\n"
    )
    
    return(NULL)
  }
  
  
  gene_list <- gene_list[
    names(gene_list) %in% overlapping_genes
  ]
  
  
  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )
  
  
  gsea_result <- tryCatch(
    {
      
      GSEA(
        geneList = gene_list,
        TERM2GENE = gene_sets,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
    },
    
    error = function(e) {
      
      message(
        "GSEA failed: ",
        e$message
      )
      
      return(NULL)
    }
  )
  
  
  if (is.null(gsea_result)) {
    
    return(NULL)
  }
  
  
  results_df <- as.data.frame(
    gsea_result
  )
  
  
  safe_comparison <- safe_name(
    comparison_name
  )
  
  
  result_dir <- file.path(
    output_dir,
    database_name,
    safe_comparison
  )
  
  
  dir.create(
    result_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # All pathways
  write.csv(
    results_df,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_GSEA_all_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # FDR < 0.25
  sig_025 <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.25
    )
  
  
  write.csv(
    sig_025,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_GSEA_FDR_0.25.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # FDR < 0.10
  sig_010 <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.10
    )
  
  
  write.csv(
    sig_010,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_GSEA_FDR_0.10.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # FDR < 0.05
  sig_005 <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    )
  
  
  write.csv(
    sig_005,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_GSEA_FDR_0.05.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  cat(
    "GSEA pathways returned:",
    nrow(results_df),
    "\n"
  )
  
  
  cat(
    "FDR < 0.25:",
    nrow(sig_025),
    "\n"
  )
  
  
  cat(
    "FDR < 0.10:",
    nrow(sig_010),
    "\n"
  )
  
  
  cat(
    "FDR < 0.05:",
    nrow(sig_005),
    "\n"
  )
  
  
  make_gsea_barplot(
    results_df = results_df,
    database_name = database_name,
    comparison_name = comparison_name,
    result_dir = result_dir,
    fdr_cutoff = 0.25,
    top_n = 10
  )
  
  make_gsea_barplot(
    results_df = results_df,
    database_name = database_name,
    comparison_name = comparison_name,
    result_dir = result_dir,
    fdr_cutoff = 0.10,
    top_n = 10
  )
  
  make_gsea_barplot(
    results_df = results_df,
    database_name = database_name,
    comparison_name = comparison_name,
    result_dir = result_dir,
    fdr_cutoff = 0.05,
    top_n = 10
  )
  
  

  # GSEA dotplot
  
  if (nrow(results_df) > 0) {
    
    dot_plot <- dotplot(
      gsea_result,
      showCategory = 20,
      split = ".sign"
    ) +
      
      facet_grid(
        . ~ .sign
      ) +
      
      scale_y_discrete(
        labels = function(x) {
          
          stringr::str_wrap(
            clean_pathway_names(x),
            width = 42
          )
        }
      ) +
      
      labs(
        title = paste0(
          database_name,
          " GSEA"
        ),
        subtitle = comparison_name
      ) +
      
      theme_bw(
        base_size = 8
      )
    
    
    dot_height <- max(
      9,
      0.5 *
        min(
          20,
          nrow(results_df)
        ) +
        3
    )
    
    
    ggsave(
      filename = file.path(
        result_dir,
        paste0(
          database_name,
          "_GSEA_dotplot.png"
        )
      ),
      plot = dot_plot,
      width = 14,
      height = dot_height,
      dpi = 300,
      limitsize = FALSE
    )
  }
  
  
  # Individual enrichment plots for top significant pathways
  
  if (nrow(sig_025) > 0) {
    
    top_pathways <- sig_025 %>%
      arrange(
        p.adjust
      ) %>%
      slice_head(
        n = 10
      ) %>%
      pull(
        ID
      )
    
    
    for (pathway in top_pathways) {
      
      safe_pathway <- safe_name(
        pathway
      )
      
      
      pathway_plot <- gseaplot2(
        gsea_result,
        geneSetID = pathway,
        title = clean_pathway_names(pathway),
        pvalue_table = TRUE
      )
      
      
      ggsave(
        filename = file.path(
          result_dir,
          paste0(
            safe_pathway,
            "_enrichment_plot.png"
          )
        ),
        plot = pathway_plot,
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }
  
  
  return(gsea_result)
}



# ORA Function

run_ora <- function(
    genes,
    background_genes,
    gene_sets,
    database_name,
    direction,
    comparison_name,
    output_dir
) {
  
  genes <- unique(
    genes
  )
  
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  
  background_genes <- unique(
    background_genes
  )
  
  
  background_genes <- background_genes[
    !is.na(background_genes) &
      background_genes != ""
  ]
  
  
  pathway_genes <- unique(
    gene_sets$gene_symbol
  )
  
  
  mapped_genes <- intersect(
    genes,
    pathway_genes
  )
  
  
  mapped_background <- intersect(
    background_genes,
    pathway_genes
  )
  
  
  cat("\n")
  cat("========================================\n")
  cat("ORA:", database_name, "\n")
  cat("Comparison:", comparison_name, "\n")
  cat("Direction:", direction, "\n")
  cat("========================================\n")
  
  
  cat(
    "Input genes:",
    length(genes),
    "\n"
  )
  
  
  cat(
    "Genes mapped:",
    length(mapped_genes),
    "\n"
  )
  
  
  cat(
    "Background genes mapped:",
    length(mapped_background),
    "\n"
  )
  
  
  if (length(mapped_genes) < 3) {
    
    cat(
      "Fewer than 3 genes map to pathway database. Skipping.\n"
    )
    
    return(NULL)
  }
  
  
  result_dir <- file.path(
    output_dir,
    database_name,
    safe_name(comparison_name),
    direction
  )
  
  
  dir.create(
    result_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ora_result <- tryCatch(
    {
      
      enricher(
        gene = mapped_genes,
        universe = mapped_background,
        TERM2GENE = gene_sets,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        minGSSize = 5,
        maxGSSize = 500
      )
    },
    
    error = function(e) {
      
      message(
        "ORA failed: ",
        e$message
      )
      
      return(NULL)
    }
  )
  
  
  if (is.null(ora_result)) {
    
    return(NULL)
  }
  
  
  results_df <- as.data.frame(
    ora_result
  )
  
  
  write.csv(
    results_df,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_",
        direction,
        "_ORA_all_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  
  # Save ORA results at FDR < 0.25, 0.10, and 0.05
  
  sig_025 <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.25
    )
  
  write.csv(
    sig_025,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_",
        direction,
        "_ORA_FDR_0.25.csv"
      )
    ),
    row.names = FALSE
  )
  
  sig_010 <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.10
    )
  
  write.csv(
    sig_010,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_",
        direction,
        "_ORA_FDR_0.10.csv"
      )
    ),
    row.names = FALSE
  )
  
  sig_005 <- results_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    )
  
  write.csv(
    sig_005,
    file.path(
      result_dir,
      paste0(
        database_name,
        "_",
        direction,
        "_ORA_FDR_0.05.csv"
      )
    ),
    row.names = FALSE
  )
  
  cat(
    "ORA pathways returned:",
    nrow(results_df),
    "\n"
  )
  
  cat(
    "FDR < 0.25:",
    nrow(sig_025),
    "\n"
  )
  
  cat(
    "FDR < 0.10:",
    nrow(sig_010),
    "\n"
  )
  
  cat(
    "FDR < 0.05:",
    nrow(sig_005),
    "\n"
  )
  

  # Save ORA barplots at all three FDR thresholds
  make_ora_barplot(
    results_df = results_df,
    database_name = database_name,
    comparison_name = comparison_name,
    direction = direction,
    result_dir = result_dir,
    fdr_cutoff = 0.25,
    top_n = 20
  )
  
  make_ora_barplot(
    results_df = results_df,
    database_name = database_name,
    comparison_name = comparison_name,
    direction = direction,
    result_dir = result_dir,
    fdr_cutoff = 0.10,
    top_n = 20
  )
  
  make_ora_barplot(
    results_df = results_df,
    database_name = database_name,
    comparison_name = comparison_name,
    direction = direction,
    result_dir = result_dir,
    fdr_cutoff = 0.05,
    top_n = 20
  )
  
  

  # ORA dotplot
  
  if (nrow(results_df) > 0) {
    
    categories_to_show <- min(
      20,
      nrow(results_df)
    )
    
    
    ora_plot <- dotplot(
      ora_result,
      showCategory = categories_to_show
    ) +
      
      scale_y_discrete(
        labels = function(x) {
          
          stringr::str_wrap(
            clean_pathway_names(x),
            width = 42
          )
        }
      ) +
      
      labs(
        title = paste0(
          database_name,
          " ORA"
        ),
        
        subtitle = paste0(
          comparison_name,
          " - ",
          direction
        )
      ) +
      
      theme_bw(
        base_size = 8
      )
    
    
    dot_height <- max(
      9,
      0.7 * categories_to_show + 2
    )
    
    
    ggsave(
      filename = file.path(
        result_dir,
        paste0(
          database_name,
          "_",
          direction,
          "_ORA_dotplot.png"
        )
      ),
      plot = ora_plot,
      width = 14,
      height = dot_height,
      dpi = 300,
      limitsize = FALSE
    )
  }
  
  
  return(ora_result)
}



# MAIN LOOP THROUGH CELL TYPES

for (cell_type in cell_types) {
  
  cat("\n\n")
  cat("############################################################\n")
  cat("CELL TYPE:", cell_type, "\n")
  cat("############################################################\n")
  
  
  safe_cell_type <- safe_name(
    cell_type
  )
  
  
  # Cell-type output directory
  
  celltype_output_dir <- file.path(
    root_output_dir,
    safe_cell_type
  )
  
  
  deg_output_dir <- file.path(
    celltype_output_dir,
    "DEG"
  )
  
  
  raw_volcano_dir <- file.path(
    celltype_output_dir,
    "Volcano_Plots",
    "Raw_Pvalue"
  )
  
  
  fdr_volcano_dir <- file.path(
    celltype_output_dir,
    "Volcano_Plots",
    "FDR_0.05"
  )
  
  
  pathway_output_dir <- file.path(
    celltype_output_dir,
    "Pathway_Analysis"
  )
  
  
  gsea_output_dir <- file.path(
    pathway_output_dir,
    "GSEA"
  )
  
  
  ora_output_dir <- file.path(
    pathway_output_dir,
    "ORA"
  )
  
  
  dir.create(
    deg_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    raw_volcano_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    fdr_volcano_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    gsea_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    ora_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # Subset to current cell type
  
  cells_keep <- rownames(
    dat@meta.data
  )[
    dat@meta.data[[celltype_col]] == cell_type &
      !is.na(dat@meta.data[[celltype_col]])
  ]
  
  
  if (length(cells_keep) == 0) {
    
    message(
      "No cells found for ",
      cell_type,
      ". Skipping."
    )
    
    next
  }
  
  
  cell_dat <- subset(
    dat,
    cells = cells_keep
  )
  
  
  cat(
    "\nNumber of cells:",
    ncol(cell_dat),
    "\n"
  )
  
  
  cat(
    "\nCells per treatment:\n"
  )
  
  print(
    table(
      cell_dat@meta.data[[group_col]]
    )
  )
  
  
  cat(
    "\nCells per treatment and biological replicate:\n"
  )
  
  print(
    table(
      cell_dat@meta.data[[group_col]],
      cell_dat@meta.data[[replicate_col]]
    )
  )
  
  
 
  # Create Pseudobulk IDs
  
  cell_dat$treatment <- as.character(
    cell_dat@meta.data[[group_col]]
  )
  
  
  cell_dat$biological_replicate <- as.character(
    cell_dat@meta.data[[replicate_col]]
  )
  
  
  # True pseudobulk = treatment + mouse replicate
  cell_dat$pseudobulk_id <- paste0(
    "treat.",
    cell_dat$treatment,
    ".mouse.",
    cell_dat$biological_replicate
  )
  
  

  # Pseudobulk metadata
  
  pb_meta <- cell_dat@meta.data %>%
    dplyr::select(
      pseudobulk_id,
      treatment,
      biological_replicate
    ) %>%
    distinct() %>%
    arrange(
      treatment,
      biological_replicate
    )
  
  
  cat(
    "\nPseudobulk samples:\n"
  )
  
  print(pb_meta)
  
  
  cat(
    "\nPseudobulk samples per treatment:\n"
  )
  
  print(
    table(
      pb_meta$treatment
    )
  )
  
  

  # Get Raw Counts

  counts <- LayerData(
    object = cell_dat,
    assay = "RNA",
    layer = "counts"
  )
  
  
  # Make sure metadata order matches count matrix
  cell_meta <- cell_dat@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  
  pseudobulk_factor <- factor(
    cell_meta$pseudobulk_id,
    levels = pb_meta$pseudobulk_id
  )
  
  
  # Sparse cell x pseudobulk aggregation matrix
  pb_design <- sparse.model.matrix(
    ~ 0 + pseudobulk_factor
  )
  
  
  colnames(pb_design) <- levels(
    pseudobulk_factor
  )
  
  
  # genes x pseudobulk samples
  pb_counts <- counts %*% pb_design
  
  
  pb_counts <- as(
    pb_counts,
    "dgCMatrix"
  )
  
  
  # Match metadata to count columns
  
  pb_meta <- as.data.frame(
    pb_meta
  )
  
  
  rownames(pb_meta) <- pb_meta$pseudobulk_id
  
  
  pb_meta <- pb_meta[
    colnames(pb_counts),
    ,
    drop = FALSE
  ]
  
  

  # Save pseudobulk data
  write.csv(
    pb_meta,
    file = file.path(
      celltype_output_dir,
      paste0(
        safe_cell_type,
        "_pseudobulk_metadata.csv"
      )
    ),
    row.names = TRUE
  )
  
  
  write.csv(
    as.data.frame(
      as.matrix(pb_counts)
    ) %>%
      rownames_to_column(
        "gene"
      ),
    
    file = file.path(
      celltype_output_dir,
      paste0(
        safe_cell_type,
        "_pseudobulk_counts.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # Treatment Comparisons
  
  treatment_ids <- sort(
    unique(
      pb_meta$treatment
    )
  )
  
  
  if (length(treatment_ids) < 2) {
    
    message(
      "Fewer than two treatments present for ",
      cell_type,
      ". Skipping DEG analysis."
    )
    
    next
  }
  
  
  comparisons <- combn(
    treatment_ids,
    2,
    simplify = FALSE
  )
  
  
  celltype_all_results <- list()
  
  

  # Run Each Pairwise comparison

  
  for (comp in comparisons) {
    
    group1 <- as.character(comp[1])
    group2 <- as.character(comp[2])
    
    group1_label <- get_group_label(group1)
    group2_label <- get_group_label(group2)
    
    comparison_name <- paste0(
      group2_label,
      " vs ",
      group1_label
    )
    
    safe_comparison <- safe_name(
      comparison_name
    )
    
    
    cat("\n")
    cat("------------------------------------------------------------\n")
    cat("DESeq2 comparison:", comparison_name, "\n")
    cat("------------------------------------------------------------\n")
    
    
    # ---------------------------------------------------------
    # Select pseudobulk samples for comparison
    # ---------------------------------------------------------
    
    samples_keep <- pb_meta$treatment %in%
      c(
        group1,
        group2
      )
    
    
    meta_sub <- pb_meta[
      samples_keep,
      ,
      drop = FALSE
    ]
    
    
    counts_sub <- pb_counts[
      ,
      rownames(meta_sub),
      drop = FALSE
    ]
    
    
    meta_sub$treatment <- factor(
      meta_sub$treatment,
      levels = c(
        group1,
        group2
      )
    )
    
    
    cat(
      "\nPseudobulk samples used:\n"
    )
    
    print(meta_sub)
    
    
    cat(
      "\nReplicates per treatment:\n"
    )
    
    print(
      table(
        meta_sub$treatment
      )
    )
    
    
    
    # Require at least 2 biological replicates per treatment
  
    replicate_table <- table(
      meta_sub$treatment
    )
    
    
    if (
      length(replicate_table) < 2 ||
      any(replicate_table < 2)
    ) {
      
      message(
        "Skipping ",
        comparison_name,
        " for ",
        cell_type,
        ": fewer than 2 pseudobulk replicates in one or more groups."
      )
      
      next
    }
    
    
   # Gene filter
    
    counts_sub_matrix <- as.matrix(
      counts_sub
    )
    
    
    storage.mode(
      counts_sub_matrix
    ) <- "integer"
    
    
    # Keep genes with at least 10 total counts across samples
    keep_genes <- rowSums(
      counts_sub_matrix
    ) >= 10
    
    
    counts_sub_matrix <- counts_sub_matrix[
      keep_genes,
      ,
      drop = FALSE
    ]
    
    
    cat(
      "\nGenes remaining after count filtering:",
      nrow(counts_sub_matrix),
      "\n"
    )
    
    
    if (nrow(counts_sub_matrix) == 0) {
      
      message(
        "No genes remain after filtering. Skipping comparison."
      )
      
      next
    }
    
    

    # DESEQ2
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts_sub_matrix,
      colData = meta_sub,
      design = ~ treatment
    )
    
    
    dds <- DESeq(
      dds
    )
    
    
    deseq_res <- results(
      dds,
      contrast = c(
        "treatment",
        group2,
        group1
      ),
      alpha = fdr_cutoff
    )
    
    
    res_df <- as.data.frame(
      deseq_res
    ) %>%
      
      rownames_to_column(
        "gene"
      ) %>%
      
      mutate(
        cell_type = cell_type,
        comparison = comparison_name,
        group1 = group1,
        group2 = group2
      ) %>%
      
      arrange(
        pvalue
      )
    
    
    celltype_all_results[[
      comparison_name
    ]] <- res_df
    
    

    # Save DESEQ2 results
    
    write.csv(
      res_df,
      file.path(
        deg_output_dir,
        paste0(
          safe_comparison,
          "_DESeq2_ALL_results.csv"
        )
      ),
      row.names = FALSE
    )
    
    
    
    # Raw P < 0.05
    
    
    raw_sig <- res_df %>%
      filter(
        !is.na(pvalue),
        pvalue < raw_p_cutoff
      ) %>%
      arrange(
        pvalue
      )
    
    
    write.csv(
      raw_sig,
      file.path(
        deg_output_dir,
        paste0(
          safe_comparison,
          "_DESeq2_RawP_less_than_0.05.csv"
        )
      ),
      row.names = FALSE
    )
    
    
    
    # FDR < 0.05
    
    
    fdr_sig <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff
      ) %>%
      arrange(
        padj
      )
    
    
    write.csv(
      fdr_sig,
      file.path(
        deg_output_dir,
        paste0(
          safe_comparison,
          "_DESeq2_FDR_less_than_0.05.csv"
        )
      ),
      row.names = FALSE
    )
    
    
    
    # FDR < 0.05 + logFC threshold
    
    
    fdr_logfc_sig <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff,
        abs(log2FoldChange) >= logfc_cutoff
      ) %>%
      arrange(
        padj
      )
    
    
    write.csv(
      fdr_logfc_sig,
      file.path(
        deg_output_dir,
        paste0(
          safe_comparison,
          "_DESeq2_FDR0.05_log2FC",
          logfc_cutoff,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    
    
    # Save top raw genes
    
    
    top_raw_up <- res_df %>%
      filter(
        !is.na(pvalue),
        pvalue < raw_p_cutoff,
        log2FoldChange >= logfc_cutoff
      ) %>%
      arrange(
        pvalue,
        desc(log2FoldChange)
      ) %>%
      slice_head(
        n = 10
      ) %>%
      mutate(
        direction = "Upregulated"
      )
    
    
    top_raw_down <- res_df %>%
      filter(
        !is.na(pvalue),
        pvalue < raw_p_cutoff,
        log2FoldChange <= -logfc_cutoff
      ) %>%
      arrange(
        pvalue,
        log2FoldChange
      ) %>%
      slice_head(
        n = 10
      ) %>%
      mutate(
        direction = "Downregulated"
      )
    
    
    write.csv(
      bind_rows(
        top_raw_up,
        top_raw_down
      ),
      
      file.path(
        deg_output_dir,
        paste0(
          safe_comparison,
          "_Top10_Up_Down_RawP.csv"
        )
      ),
      
      row.names = FALSE
    )
    
    
    
    # Save top FDR genes
    
    
    top_fdr_up <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff,
        log2FoldChange >= logfc_cutoff
      ) %>%
      arrange(
        padj,
        desc(log2FoldChange)
      ) %>%
      slice_head(
        n = 10
      ) %>%
      mutate(
        direction = "Upregulated"
      )
    
    
    top_fdr_down <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff,
        log2FoldChange <= -logfc_cutoff
      ) %>%
      arrange(
        padj,
        log2FoldChange
      ) %>%
      slice_head(
        n = 10
      ) %>%
      mutate(
        direction = "Downregulated"
      )
    
    
    write.csv(
      bind_rows(
        top_fdr_up,
        top_fdr_down
      ),
      
      file.path(
        deg_output_dir,
        paste0(
          safe_comparison,
          "_Top10_Up_Down_FDR0.05.csv"
        )
      ),
      
      row.names = FALSE
    )
    
    
    
    # Volcano Plots
    
    
    make_raw_p_volcano(
      results_df = res_df,
      comparison_name = comparison_name,
      output_dir = raw_volcano_dir
    )
    
    
    make_fdr_volcano(
      results_df = res_df,
      comparison_name = comparison_name,
      output_dir = fdr_volcano_dir
    )
    
    

    # Pathway analysis from DESEQ2 table
    
    cat("\n")
    cat("Starting pathway analysis for:", comparison_name, "\n")
    
    
    # Background = genes successfully tested by DESeq2
    background_genes <- res_df %>%
      filter(
        !is.na(pvalue)
      ) %>%
      pull(
        gene
      ) %>%
      unique()
    
    
 
    # GSEA
    # Uses ALL DESeq2 genes ranked by Wald statistic
    
    
    run_gsea(
      deg_df = res_df,
      gene_sets = hallmark_sets,
      database_name = "Hallmark",
      comparison_name = comparison_name,
      output_dir = gsea_output_dir
    )
    
    
    run_gsea(
      deg_df = res_df,
      gene_sets = kegg_sets,
      database_name = "KEGG",
      comparison_name = comparison_name,
      output_dir = gsea_output_dir
    )
    
    
    run_gsea(
      deg_df = res_df,
      gene_sets = go_sets,
      database_name = "GO",
      comparison_name = comparison_name,
      output_dir = gsea_output_dir
    )
    
    
    
    # ORA
    # Uses DESeq2 FDR < 0.05 genes
    
    
    significant_genes <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff
      ) %>%
      pull(
        gene
      ) %>%
      unique()
    
    
    upregulated_genes <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff,
        log2FoldChange > 0
      ) %>%
      pull(
        gene
      ) %>%
      unique()
    
    
    downregulated_genes <- res_df %>%
      filter(
        !is.na(padj),
        padj < fdr_cutoff,
        log2FoldChange < 0
      ) %>%
      pull(
        gene
      ) %>%
      unique()
    
    
    for (db in c(
      "Hallmark",
      "KEGG",
      "GO"
    )) {
      
      gene_sets <- switch(
        db,
        Hallmark = hallmark_sets,
        KEGG = kegg_sets,
        GO = go_sets
      )
      
      
      
      # All significant genes
      
      run_ora(
        genes = significant_genes,
        background_genes = background_genes,
        gene_sets = gene_sets,
        database_name = db,
        direction = "All_Significant",
        comparison_name = comparison_name,
        output_dir = ora_output_dir
      )
      
      
      
      # Upregulated genes
      
      
      run_ora(
        genes = upregulated_genes,
        background_genes = background_genes,
        gene_sets = gene_sets,
        database_name = db,
        direction = "Upregulated",
        comparison_name = comparison_name,
        output_dir = ora_output_dir
      )
      
      
      
      # Downregulated genes
     
      
      run_ora(
        genes = downregulated_genes,
        background_genes = background_genes,
        gene_sets = gene_sets,
        database_name = db,
        direction = "Downregulated",
        comparison_name = comparison_name,
        output_dir = ora_output_dir
      )
    }
    
    
    cat(
      "\nCompleted:",
      cell_type,
      "-",
      comparison_name,
      "\n"
    )
  }
  
  

  # Save cobmined DEG table for cell type

  
  if (length(celltype_all_results) > 0) {
    
    combined_celltype_results <- bind_rows(
      celltype_all_results
    )
    
    
    write.csv(
      combined_celltype_results,
      file.path(
        celltype_output_dir,
        paste0(
          safe_cell_type,
          "_ALL_DESeq2_comparisons.csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  
  cat("\n")
  cat("============================================================\n")
  cat("COMPLETED CELL TYPE:", cell_type, "\n")
  cat("============================================================\n")
}



cat("\n\n")
cat("############################################################\n")
cat("ALL GEX ANALYSES COMPLETE\n")
cat("############################################################\n")

cat(
  "\nResults saved to:\n",
  root_output_dir,
  "\n"
)