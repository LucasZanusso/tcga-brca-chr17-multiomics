# =============================================================================
# 05_rnaseq.R — RNA-seq Differential Expression Analysis
# Project : TCGA-BRCA Multi-Omics (chr17 focus)
# Input   : data/rna/tumor/*.tsv  — STAR augmented gene counts (GDC open access)
#           data/rna/normal/*.tsv — same format, solid tissue normal (4 cases)
# Output  : figures/rna/
#
# Design  : ~condition (Tumor vs Normal), unpaired
#           9 tumors vs 4 normals — unbalanced but valid for DESeq2
#
# STAR augmented counts format (9 columns):
#   gene_id | gene_name | gene_type | unstranded | stranded_first |
#   stranded_second | tpm_unstranded | fpkm_unstranded | fpkm_uq_unstranded
#
#   → Use 'unstranded' as count input for DESeq2
#   → Use 'tpm_unstranded' for expression visualization
#   → Skip first 4 rows (summary stats: N_unmapped, N_multimapping, etc.)
#   → Strip Ensembl version from gene_id (ENSG00000000003.15 → ENSG00000000003)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggrepel)       # non-overlapping labels on volcano plot
  library(pheatmap)      # heatmap
})

# Bioconductor — install once with:
# BiocManager::install(c("DESeq2", "apeglm", "EnhancedVolcano",
#                        "org.Hs.eg.db", "AnnotationDbi"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
RNA_TUMOR_DIR  <- "data/rna/tumor"
RNA_NORMAL_DIR <- "data/rna/normal"
FIGURES_DIR    <- "figures/rna"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# Significance thresholds
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1.0     # |log2FoldChange| > 1 → ~2x expression change

# Genes of interest — chr17 focus + key cancer genes
GENES_OF_INTEREST <- c("ERBB2", "TP53", "BRCA1",
                        "EGFR", "MKI67", "ESR1", "PGR")

# Metadata — same as 04_cnv_analysis.R for consistent coloring
METADATA <- tribble(
  ~case_id,        ~pam50,           ~her2,
  "TCGA-B6-A0I9",  "HER2-enriched",  "Positive",
  "TCGA-BH-A18R",  "HER2-enriched",  "Positive",
  "TCGA-BH-A0DZ",  "HER2-enriched",  "Positive",
  "TCGA-A1-A0SK",  "Basal-like",     "Negative",
  "TCGA-A2-A0CM",  "Basal-like",     "Negative",
  "TCGA-BH-A18V",  "Basal-like",     "Negative",
  "TCGA-BH-A18Q",  "Basal-like",     "Negative",
  "TCGA-A2-A0CU",  "Luminal A",      "Negative",
  "TCGA-AR-A0TR",  "Luminal A",      "Negative"
)

SUBTYPE_COLORS <- c(
  "HER2-enriched" = "#E84855",
  "Basal-like"    = "#3B4CC0",
  "Luminal A"     = "#59A14F",
  "Normal"        = "#B8B8B8"
)

SAMPLE_ORDER <- METADATA |>
  arrange(factor(pam50, levels = c("HER2-enriched", "Basal-like", "Luminal A"))) |>
  pull(case_id)

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------
message("=== Loading STAR gene count files ===")

extract_case_id <- function(filepath) {
  basename(filepath) |> str_extract("^TCGA-[A-Z0-9]+-[A-Z0-9]+")
}

load_star_file <- function(filepath, tissue_type) {
  case_id <- extract_case_id(filepath)

  # Skip first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
  raw <- read_tsv(filepath,
                  comment    = "#",           # ignora linha do GENCODE model
                  col_types  = cols(.default = col_character()),
                  show_col_types = FALSE) |>
    filter(!str_starts(gene_id, "N_")) |>   # remove N_unmapped etc.
    mutate(
      gene_id_clean = str_remove(gene_id, "\\.\\d+$"),
      case_id       = case_id,
      tissue_type   = tissue_type,
      unstranded    = as.integer(unstranded),
      tpm           = as.numeric(tpm_unstranded)
    ) |>
    dplyr::select(gene_id = gene_id_clean, gene_name, gene_type,
                  unstranded, tpm, case_id, tissue_type)
}



tumor_files  <- list.files(RNA_TUMOR_DIR,  pattern = "\\.tsv$", full.names = TRUE)
normal_files <- list.files(RNA_NORMAL_DIR, pattern = "\\.tsv$", full.names = TRUE)

if (length(tumor_files) == 0)
  stop("No tumor RNA-seq files found in ", RNA_TUMOR_DIR,
       ". Run 00_download.sh first.")

counts_long <- bind_rows(
  map_dfr(tumor_files,  load_star_file, tissue_type = "Tumor"),
  map_dfr(normal_files, load_star_file, tissue_type = "Normal")
)

message("  Tumor samples  : ", n_distinct(counts_long$case_id[counts_long$tissue_type == "Tumor"]))
message("  Normal samples : ", n_distinct(counts_long$case_id[counts_long$tissue_type == "Normal"]))
message("  Genes loaded   : ", n_distinct(counts_long$gene_id))

# ---------------------------------------------------------------------------
# 2. BUILD COUNT MATRIX FOR DESEQ2
#    Rows = genes, columns = samples
#    Keep only protein-coding genes (reduces noise, speeds up analysis)
# ---------------------------------------------------------------------------
message("\n=== Building count matrix ===")

counts_coding <- counts_long |>
  filter(gene_type == "protein_coding")

count_matrix <- counts_coding |>
  select(gene_id, case_id, tissue_type, unstranded) |>
  mutate(sample_id = paste(case_id, tissue_type, sep = "_")) |>
  select(gene_id, sample_id, unstranded) |>
  pivot_wider(names_from = sample_id, values_from = unstranded) |>
  column_to_rownames("gene_id") |>
  as.matrix()

# Remove genes with zero counts across all samples
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
message("  Protein-coding genes after filtering: ", nrow(count_matrix))

# Sample metadata for DESeq2
sample_meta <- tibble(
  sample_id    = colnames(count_matrix),
  case_id      = str_extract(sample_id, "^TCGA-[A-Z0-9]+-[A-Z0-9]+"),
  tissue_type  = str_extract(sample_id, "(Tumor|Normal)$"),
  condition    = factor(tissue_type, levels = c("Normal", "Tumor"))
) |>
  left_join(METADATA, by = "case_id") |>
  mutate(pam50 = replace_na(pam50, "Normal")) |>
  column_to_rownames("sample_id")

# Ensure column order matches
count_matrix <- count_matrix[, rownames(sample_meta)]

# ---------------------------------------------------------------------------
# 3. DESEQ2 — differential expression
# ---------------------------------------------------------------------------
message("\n=== Running DESeq2 ===")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_meta,
  design    = ~ condition     # Tumor vs Normal, unpaired
)

# Pre-filter: keep genes with >= 10 counts in at least 4 samples
keep <- rowSums(counts(dds) >= 10) >= 4
dds  <- dds[keep, ]
message("  Genes after pre-filtering (>=10 counts in >=4 samples): ", sum(keep))

dds <- DESeq(dds, quiet = TRUE)

# Shrink LFC with apeglm (more accurate than standard MLE for small N)
res_shrunk <- lfcShrink(dds,
                        coef  = "condition_Tumor_vs_Normal",
                        type  = "apeglm",
                        quiet = TRUE)

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  # Add gene symbols from Ensembl IDs
  mutate(
    gene_name = mapIds(org.Hs.eg.db,
                       keys    = gene_id,
                       column  = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first"),
    sig = case_when(
      padj < PADJ_THRESHOLD & log2FoldChange >  LFC_THRESHOLD ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up", "Down", "NS"))
  ) |>
  arrange(padj)

n_up   <- sum(res_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(res_df$sig == "Down", na.rm = TRUE)
message("  DEGs (padj < ", PADJ_THRESHOLD, ", |LFC| > ", LFC_THRESHOLD, "):")
message("    Up-regulated   : ", n_up)
message("    Down-regulated : ", n_down)

write_tsv(res_df, file.path(FIGURES_DIR, "deseq2_results.tsv"))
message("  Saved: figures/rna/deseq2_results.tsv")

# ---------------------------------------------------------------------------
# 4. PCA — sample-level quality check
# ---------------------------------------------------------------------------
message("\n=== Plot 1: PCA ===")

# VST normalisation for PCA and heatmap (variance-stabilising transformation)
vst_data <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vst_data,
                    intgroup = c("condition", "pam50"),
                    returnData = TRUE)

pct_var <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data,
                aes(x = PC1, y = PC2,
                    color = pam50, shape = condition, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 2.8, max.overlaps = 20,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(values = SUBTYPE_COLORS, name = "PAM50 / type") +
  scale_shape_manual(values = c(Normal = 1, Tumor = 16), name = "Tissue") +
  labs(
    title    = "PCA — RNA-seq Expression",
    subtitle = "VST-normalised counts · TCGA-BRCA · 9 tumors + 4 normals",
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(FIGURES_DIR, "01_pca.pdf"),
       p_pca, width = 8, height = 6)
message("  Saved: figures/rna/01_pca.pdf")

# ---------------------------------------------------------------------------
# 5. VOLCANO PLOT
# ---------------------------------------------------------------------------
message("\n=== Plot 2: Volcano plot ===")

# Label top 15 DEGs + genes of interest
top_genes <- res_df |>
  filter(sig != "NS") |>
  slice_min(padj, n = 15) |>
  pull(gene_name)

label_genes <- union(top_genes, GENES_OF_INTEREST)

plot_data <- res_df |>
  filter(!is.na(padj), !is.na(log2FoldChange)) |>
  mutate(
    neg_log10_padj = -log10(padj),
    label = if_else(gene_name %in% label_genes, gene_name, NA_character_)
  )

p_volcano <- ggplot(plot_data,
                    aes(x = log2FoldChange, y = neg_log10_padj,
                        color = sig, label = label)) +
  geom_point(aes(size = sig), alpha = 0.6) +
  geom_text_repel(
    color = "black", size = 2.8, fontface = "bold",
    max.overlaps = 25, box.padding = 0.4,
    na.rm = TRUE
  ) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_color_manual(
    values = c(Up = "#B40426", Down = "#4393C3", NS = "grey70"),
    labels = c(Up   = paste0("Up (n=", n_up, ")"),
               Down = paste0("Down (n=", n_down, ")"),
               NS   = "Not significant"),
    name = NULL
  ) +
  scale_size_manual(values = c(Up = 1.8, Down = 1.8, NS = 0.8),
                    guide = "none") +
  labs(
    title    = "Differential Expression — Tumor vs Normal",
    subtitle = "DESeq2 + apeglm shrinkage · TCGA-BRCA · 9 tumors vs 4 normals",
    x        = "log2 Fold Change (Tumor / Normal)",
    y        = "-log10(adjusted p-value)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(FIGURES_DIR, "02_volcano.pdf"),
       p_volcano, width = 9, height = 7)
message("  Saved: figures/rna/02_volcano.pdf")

# ---------------------------------------------------------------------------
# 6. HEATMAP — top 50 DEGs
# ---------------------------------------------------------------------------
message("\n=== Plot 3: Heatmap top DEGs ===")

top50_genes <- res_df |>
  filter(sig != "NS") |>
  slice_min(padj, n = 50) |>
  pull(gene_id)

# VST matrix, subset to top 50 genes
vst_mat <- assay(vst_data)
vst_top <- vst_mat[rownames(vst_mat) %in% top50_genes, ]

# Replace Ensembl IDs with gene symbols for row labels
rowname_map <- res_df |>
  filter(gene_id %in% rownames(vst_top)) |>
  select(gene_id, gene_name) |>
  deframe()
rownames(vst_top) <- rowname_map[rownames(vst_top)]

# Scale by row (z-score) for readable heatmap
vst_scaled <- t(scale(t(vst_top)))

# Column annotations
col_annot <- sample_meta |>
  select(condition, pam50) |>
  rename(Tissue = condition, PAM50 = pam50)

annot_colors <- list(
  Tissue = c(Tumor = "#E84855", Normal = "#B8B8B8"),
  PAM50  = SUBTYPE_COLORS
)

pdf(file.path(FIGURES_DIR, "03_heatmap_top50.pdf"), width = 10, height = 12)
pheatmap(
  vst_scaled,
  annotation_col  = col_annot,
  annotation_colors = annot_colors,
  color           = colorRampPalette(c("#4393C3", "white", "#B40426"))(100),
  breaks          = seq(-3, 3, length.out = 101),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 7,
  fontsize_col    = 8,
  main            = "Top 50 DEGs — Tumor vs Normal (z-score VST)"
)
dev.off()
message("  Saved: figures/rna/03_heatmap_top50.pdf")

# ---------------------------------------------------------------------------
# 7. GENES OF INTEREST — ERBB2, TP53, BRCA1 expression per sample
# ---------------------------------------------------------------------------
message("\n=== Plot 4: Expression of ERBB2, TP53, BRCA1 ===")

# Use TPM for visualisation (more interpretable than raw counts)
tpm_interest <- counts_long |>
  filter(gene_name %in% GENES_OF_INTEREST[1:3],  # ERBB2, TP53, BRCA1
         gene_type == "protein_coding") |>
  left_join(METADATA, by = "case_id") |>
  mutate(
    pam50   = replace_na(pam50, "Normal"),
    case_id = if_else(tissue_type == "Normal",
                      paste0(case_id, "\n(N)"), case_id),
    gene_name = factor(gene_name, levels = c("ERBB2", "TP53", "BRCA1"))
  )

# Also pull DESeq2 results for these genes to annotate significance
sig_labels <- res_df |>
  filter(gene_name %in% c("ERBB2", "TP53", "BRCA1")) |>
  select(gene_name, log2FoldChange, padj, sig) |>
  mutate(
    label = sprintf("LFC = %.2f\npadj = %.3f", log2FoldChange, padj)
  )

p_genes <- ggplot(tpm_interest,
                  aes(x = reorder(case_id, tpm),
                      y = log1p(tpm),
                      fill = pam50)) +
  geom_col(width = 0.75) +
  geom_text(data = sig_labels,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.3, size = 2.5, color = "grey30") +
  scale_fill_manual(values = SUBTYPE_COLORS, name = "PAM50") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 1) +
  labs(
    title    = "Expression of ERBB2, TP53, BRCA1",
    subtitle = "log1p(TPM) · TCGA-BRCA · tumor + normal samples",
    x        = NULL,
    y        = "log1p(TPM)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 7),
    strip.background   = element_rect(fill = "#2C3E50"),
    strip.text         = element_text(color = "white", face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

ggsave(file.path(FIGURES_DIR, "04_genes_of_interest.pdf"),
       p_genes, width = 10, height = 9)
message("  Saved: figures/rna/04_genes_of_interest.pdf")

# ---------------------------------------------------------------------------
# 8. EXPORT — normalised expression table for 07_integration.R
# ---------------------------------------------------------------------------
message("\n=== Exporting normalised expression table ===")

# VST matrix — all samples, all filtered genes
vst_export <- as.data.frame(assay(vst_data)) |>
  rownames_to_column("gene_id") |>
  left_join(res_df |> select(gene_id, gene_name) |> distinct(),
            by = "gene_id") |>
  relocate(gene_name, .after = gene_id)

write_tsv(vst_export,
          file.path(FIGURES_DIR, "vst_expression_matrix.tsv"))

# DESeq2 results — full table already saved above
message("  Saved: figures/rna/vst_expression_matrix.tsv")

# ---------------------------------------------------------------------------
# 9. SESSION INFO
# ---------------------------------------------------------------------------
sink(file.path(FIGURES_DIR, "session_info.txt"))
sessionInfo()
sink()

message("\n=== 05_rnaseq.R complete ===")
message("Outputs in figures/rna/:")
message("  01_pca.pdf                   — sample clustering, VST counts")
message("  02_volcano.pdf               — DEG overview, Tumor vs Normal")
message("  03_heatmap_top50.pdf         — top 50 DEGs, z-score heatmap")
message("  04_genes_of_interest.pdf     — ERBB2 / TP53 / BRCA1 expression")
message("  deseq2_results.tsv           — full DESeq2 results table")
message("  vst_expression_matrix.tsv    — input for 07_integration.R")
