# =============================================================================
# 06_methylation.R — DNA Methylation Analysis
# Project : TCGA-BRCA Multi-Omics (chr17 focus)
# Input   : data/methylation/tumor/*.txt   — beta values (GDC SeSAMe pipeline)
#           data/methylation/normal/*.txt  — same format, solid tissue normal
# Output  : figures/methylation/
#
# Data format (GDC Methylation Beta Value, no header):
#   col1: probe ID (cgXXXXXXXX — Illumina EPIC array)
#   col2: beta value (0–1, proportion of methylated molecules)
#
# Analysis:
#   1. Load + build beta matrix
#   2. Filter probes (remove NAs, low-variance, non-CpG)
#   3. M-value transformation for differential methylation (limma)
#   4. Differential methylation: Tumor vs Normal (limma, 4 paired normals)
#   5. Chr17 methylation overview — ERBB2, TP53, BRCA1 loci
#   6. Heatmap — top differentially methylated probes
#   7. Export for 07_integration.R
#
# Why M-values for statistics:
#   Beta values are bounded [0,1] — violates normality assumptions.
#   M-value = log2(beta / (1 - beta)) is unbounded and approximately normal,
#   making it suitable for linear models (limma).
#   Beta values are kept for visualization (more interpretable biologically).
# =============================================================================
setwd('/home/lucas/projeto_brca')
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
})

# Bioconductor — install once with:
# BiocManager::install(c("limma", "minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg38"))
suppressPackageStartupMessages({
  library(limma)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
METH_TUMOR_DIR  <- "data/methylation/tumor"
METH_NORMAL_DIR <- "data/methylation/normal"
FIGURES_DIR     <- "figures/methylation"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# Significance thresholds
ADJP_THRESHOLD  <- 0.05
DELTA_THRESHOLD <- 0.2    # |delta beta| > 0.2 — meaningful methylation change

# Metadata — consistent with 04 and 05
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

# Chr17 gene loci for annotation
GENE_COORDS <- tribble(
  ~gene,   ~chr,  ~start,   ~end,
  "ERBB2", "chr17", 39687914, 39730426,
  "TP53",  "chr17",  7668402,  7687538,
  "BRCA1", "chr17", 43044295, 43125483
)

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------
message("=== Loading methylation beta value files ===")

extract_case_id <- function(filepath) {
  basename(filepath) |> str_extract("^TCGA-[A-Z0-9]+-[A-Z0-9]+")
}

load_beta_file <- function(filepath, tissue_type) {
  case_id <- extract_case_id(filepath)
  read_tsv(filepath,
           col_names      = c("probe_id", "beta"),
           col_types      = cols(probe_id = col_character(),
                                 beta     = col_double()),
           show_col_types = FALSE) |>
    mutate(case_id     = case_id,
           tissue_type = tissue_type)
}

tumor_files  <- list.files(METH_TUMOR_DIR,  pattern = "\\.txt$", full.names = TRUE)
normal_files <- list.files(METH_NORMAL_DIR, pattern = "\\.txt$", full.names = TRUE)

if (length(tumor_files) == 0)
  stop("No tumor methylation files found in ", METH_TUMOR_DIR,
       ". Run 00_download.sh first.")

message("  Loading ", length(tumor_files), " tumor files...")
beta_tumor  <- map_dfr(tumor_files,  load_beta_file, tissue_type = "Tumor")

message("  Loading ", length(normal_files), " normal files...")
beta_normal <- map_dfr(normal_files, load_beta_file, tissue_type = "Normal")

beta_all <- bind_rows(beta_tumor, beta_normal)

message("  Tumor samples  : ", n_distinct(beta_tumor$case_id))
message("  Normal samples : ", n_distinct(beta_normal$case_id))
message("  Total probes   : ", n_distinct(beta_all$probe_id))

# ---------------------------------------------------------------------------
# 2. BUILD BETA MATRIX
#    Rows = probes, columns = samples
# ---------------------------------------------------------------------------
message("\n=== Building beta matrix ===")

beta_matrix <- beta_all |>
  mutate(sample_id = paste(case_id, tissue_type, sep = "_")) |>
  select(probe_id, sample_id, beta) |>
  pivot_wider(names_from = sample_id, values_from = beta) |>
  column_to_rownames("probe_id") |>
  as.matrix()

message("  Matrix dimensions: ", nrow(beta_matrix), " probes x ",
        ncol(beta_matrix), " samples")

# ---------------------------------------------------------------------------
# 3. FILTER PROBES
#    Remove: probes with any NA, non-CpG probes (ch. prefix),
#            probes on sex chromosomes, low-variance probes
# ---------------------------------------------------------------------------
message("\n=== Filtering probes ===")

# Remove probes with any NA
keep_complete <- complete.cases(beta_matrix)
beta_matrix   <- beta_matrix[keep_complete, ]
message("  After removing NAs       : ", nrow(beta_matrix), " probes")

# Keep only CpG probes (cg prefix — removes ch. probes which are non-CpG)
keep_cg     <- str_starts(rownames(beta_matrix), "cg")
beta_matrix <- beta_matrix[keep_cg, ]
message("  After keeping cg probes  : ", nrow(beta_matrix), " probes")

# Remove low-variance probes (bottom 25% by IQR)
probe_iqr   <- apply(beta_matrix, 1, IQR)
keep_var    <- probe_iqr > quantile(probe_iqr, 0.25)
beta_matrix <- beta_matrix[keep_var, ]
message("  After IQR filtering      : ", nrow(beta_matrix), " probes")

# ---------------------------------------------------------------------------
# 4. ANNOTATE PROBES — add chr, position, gene from EPIC array manifest
# ---------------------------------------------------------------------------
message("\n=== Annotating probes ===")

epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) |>
  as.data.frame() |>
  rownames_to_column("probe_id") |>
  select(probe_id, chr, pos, strand, UCSC_RefGene_Name, Relation_to_Island) |>
  as_tibble()

# Probes on chr17 — for locus-level analysis
probes_chr17 <- epic_anno |>
  filter(chr == "chr17",
         probe_id %in% rownames(beta_matrix))

message("  chr17 probes after filtering: ", nrow(probes_chr17))

# ---------------------------------------------------------------------------
# 5. M-VALUE TRANSFORMATION
#    M = log2(beta / (1 - beta))
#    Used for limma differential methylation — beta used for visualization
# ---------------------------------------------------------------------------
message("\n=== Computing M-values ===")

# Clip beta values away from 0 and 1 to avoid infinite M-values
beta_clipped <- pmax(pmin(beta_matrix, 0.999), 0.001)
m_matrix     <- log2(beta_clipped / (1 - beta_clipped))

# ---------------------------------------------------------------------------
# 6. DIFFERENTIAL METHYLATION — limma (Tumor vs Normal)
#    Uses all 9 tumors and 4 normals (unpaired design, same as DESeq2)
# ---------------------------------------------------------------------------
message("\n=== Differential methylation with limma ===")

# Sample metadata
sample_meta <- tibble(
  sample_id   = colnames(m_matrix),
  case_id     = str_extract(sample_id, "^TCGA-[A-Z0-9]+-[A-Z0-9]+"),
  tissue_type = str_extract(sample_id, "(Tumor|Normal)$"),
  condition   = factor(tissue_type, levels = c("Normal", "Tumor"))
) |>
  left_join(METADATA, by = "case_id") |>
  mutate(pam50 = replace_na(pam50, "Normal"))

# Design matrix — Tumor vs Normal
design <- model.matrix(~ condition, data = sample_meta)

# Fit linear model
fit  <- lmFit(m_matrix, design)
fit  <- eBayes(fit)

# Extract results for Tumor vs Normal contrast
dmp_results <- topTable(fit,
                        coef      = "conditionTumor",
                        number    = Inf,
                        sort.by   = "p") |>
  rownames_to_column("probe_id") |>
  as_tibble() |>
  left_join(epic_anno, by = "probe_id") |>
  rename(log2FC_M = logFC, adjP = adj.P.Val)

# Calculate mean beta per group for delta beta
mean_beta_tumor  <- rowMeans(beta_matrix[, sample_meta$tissue_type == "Tumor"],
                             na.rm = TRUE)
mean_beta_normal <- rowMeans(beta_matrix[, sample_meta$tissue_type == "Normal"],
                             na.rm = TRUE)

delta_beta_df <- tibble(
  probe_id    = names(mean_beta_tumor),
  beta_tumor  = mean_beta_tumor,
  beta_normal = mean_beta_normal,
  delta_beta  = mean_beta_tumor - mean_beta_normal
)

dmp_results <- dmp_results |>
  left_join(delta_beta_df, by = "probe_id") |>
  mutate(
    sig = case_when(
      adjP < ADJP_THRESHOLD & delta_beta >  DELTA_THRESHOLD ~ "Hypermethylated",
      adjP < ADJP_THRESHOLD & delta_beta < -DELTA_THRESHOLD ~ "Hypomethylated",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Hypermethylated", "Hypomethylated", "NS"))
  )

n_hyper <- sum(dmp_results$sig == "Hypermethylated", na.rm = TRUE)
n_hypo  <- sum(dmp_results$sig == "Hypomethylated",  na.rm = TRUE)
message("  Differentially methylated probes:")
message("    Hypermethylated: ", n_hyper)
message("    Hypomethylated : ", n_hypo)

write_tsv(dmp_results, file.path(FIGURES_DIR, "dmp_results.tsv"))
message("  Saved: figures/methylation/dmp_results.tsv")

# ---------------------------------------------------------------------------
# 7. VOLCANO PLOT — delta beta vs -log10(adjP)
# ---------------------------------------------------------------------------
message("\n=== Plot 1: Volcano plot ===")

# Label top probes near genes of interest
top_probes <- dmp_results |>
  filter(sig != "NS") |>
  slice_min(adjP, n = 20) |>
  pull(probe_id)

genes_interest_probes <- dmp_results |>
  filter(str_detect(UCSC_RefGene_Name,
                    "ERBB2|TP53|BRCA1")) |>
  pull(probe_id)

label_probes <- union(top_probes, genes_interest_probes)

plot_data <- dmp_results |>
  filter(!is.na(adjP)) |>
  mutate(
    neg_log10_adjp = -log10(adjP),
    label = if_else(probe_id %in% label_probes,
                    paste0(probe_id, "\n", UCSC_RefGene_Name),
                    NA_character_)
  )

p_volcano <- ggplot(plot_data,
                    aes(x = delta_beta, y = neg_log10_adjp,
                        color = sig, label = label)) +
  geom_point(aes(size = sig), alpha = 0.5) +
  geom_text_repel(color = "black", size = 2.2, max.overlaps = 15,
                  box.padding = 0.4, na.rm = TRUE) +
  geom_vline(xintercept = c(-DELTA_THRESHOLD, DELTA_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(ADJP_THRESHOLD),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_color_manual(
    values = c(Hypermethylated = "#B40426",
               Hypomethylated  = "#4393C3",
               NS              = "grey70"),
    labels = c(Hypermethylated = paste0("Hypermethylated (n=", n_hyper, ")"),
               Hypomethylated  = paste0("Hypomethylated (n=", n_hypo, ")"),
               NS              = "Not significant"),
    name = NULL
  ) +
  scale_size_manual(
    values = c(Hypermethylated = 1.5, Hypomethylated = 1.5, NS = 0.6),
    guide  = "none"
  ) +
  labs(
    title    = "Differential Methylation - Tumor vs Normal",
    subtitle = "limma · TCGA-BRCA · 9 tumors vs 4 normals",
    x        = "Delta Beta (Tumor - Normal)",
    y        = "-log10(adjusted p-value)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(FIGURES_DIR, "01_volcano.pdf"),
       p_volcano, width = 9, height = 7)
message("  Saved: figures/methylation/01_volcano.pdf")

# ---------------------------------------------------------------------------
# 8. CHR17 METHYLATION OVERVIEW — beta values along chr17
# ---------------------------------------------------------------------------
message("\n=== Plot 2: chr17 methylation overview ===")

chr17_beta <- beta_all |>
  filter(probe_id %in% probes_chr17$probe_id) |>
  left_join(probes_chr17 |> select(probe_id, pos), by = "probe_id") |>
  left_join(METADATA, by = "case_id") |>
  mutate(
    pam50   = replace_na(pam50, "Normal"),
    group   = if_else(tissue_type == "Normal", "Normal", pam50),
    group   = factor(group, levels = c("HER2-enriched", "Basal-like",
                                       "Luminal A", "Normal"))
  )

# Mean beta per probe per group for cleaner visualization
chr17_mean <- chr17_beta |>
  group_by(probe_id, pos, group) |>
  summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop")

p_chr17 <- ggplot(chr17_mean,
                  aes(x = pos / 1e6, y = mean_beta, color = group)) +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE,
              linewidth = 0.8, na.rm = TRUE) +
  # Gene markers
  geom_vline(data = GENE_COORDS,
             aes(xintercept = (start + end) / 2 / 1e6),
             color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(data = GENE_COORDS,
            aes(x = (start + end) / 2 / 1e6, y = 1.02, label = gene),
            inherit.aes = FALSE, size = 2.8, fontface = "bold", vjust = 0) +
  scale_color_manual(values = SUBTYPE_COLORS, name = "Group") +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(
    title    = "chr17 Methylation Profile",
    subtitle = "Mean beta value per probe · LOESS smoothing · TCGA-BRCA",
    x        = "Position (chr17)",
    y        = "Beta value (methylation level)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

ggsave(file.path(FIGURES_DIR, "02_chr17_methylation.pdf"),
       p_chr17, width = 13, height = 5)
message("  Saved: figures/methylation/02_chr17_methylation.pdf")

# ---------------------------------------------------------------------------
# 9. HEATMAP — top 50 differentially methylated probes
# ---------------------------------------------------------------------------
message("\n=== Plot 3: Heatmap top DMPs ===")

top50_probes <- dmp_results |>
  filter(sig != "NS") |>
  slice_min(adjP, n = 50) |>
  pull(probe_id)

beta_top <- beta_matrix[rownames(beta_matrix) %in% top50_probes, ]

# Row labels: probe + gene name
row_labels <- dmp_results |>
  filter(probe_id %in% rownames(beta_top)) |>
  mutate(label = if_else(UCSC_RefGene_Name != "",
                         paste0(probe_id, " (", UCSC_RefGene_Name, ")"),
                         probe_id)) |>
  select(probe_id, label) |>
  deframe()
rownames(beta_top) <- row_labels[rownames(beta_top)]

# Column annotations
col_annot <- sample_meta |>
  select(sample_id, condition, pam50) |>
  mutate(pam50 = replace_na(pam50, "Normal")) |>
  column_to_rownames("sample_id") |>
  rename(Tissue = condition, PAM50 = pam50)

annot_colors <- list(
  Tissue = c(Tumor = "#E84855", Normal = "#B8B8B8"),
  PAM50  = SUBTYPE_COLORS
)

pdf(file.path(FIGURES_DIR, "03_heatmap_top50.pdf"), width = 11, height = 13)
pheatmap(
  beta_top,
  annotation_col    = col_annot,
  annotation_colors = annot_colors,
  color             = colorRampPalette(c("#4393C3", "white", "#B40426"))(100),
  breaks            = seq(0, 1, length.out = 101),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 6,
  fontsize_col      = 8,
  main              = "Top 50 DMPs - Tumor vs Normal (beta value)"
)
dev.off()
message("  Saved: figures/methylation/03_heatmap_top50.pdf")

# ---------------------------------------------------------------------------
# 10. GENES OF INTEREST — ERBB2, TP53, BRCA1 probe-level beta values
# ---------------------------------------------------------------------------
message("\n=== Plot 4: ERBB2, TP53, BRCA1 methylation ===")

genes_probes <- dmp_results |>
  filter(str_detect(UCSC_RefGene_Name, "ERBB2|TP53|BRCA1"),
         probe_id %in% rownames(beta_matrix)) |>
  mutate(target_gene = case_when(
    str_detect(UCSC_RefGene_Name, "ERBB2") ~ "ERBB2",
    str_detect(UCSC_RefGene_Name, "TP53")  ~ "TP53",
    str_detect(UCSC_RefGene_Name, "BRCA1") ~ "BRCA1"
  ))

if (nrow(genes_probes) > 0) {
  beta_genes <- beta_all |>
    filter(probe_id %in% genes_probes$probe_id) |>
    left_join(genes_probes |> select(probe_id, target_gene, pos,
                                     Relation_to_Island, sig),
              by = "probe_id") |>
    left_join(METADATA, by = "case_id") |>
    mutate(
      pam50 = replace_na(pam50, "Normal"),
      group = if_else(tissue_type == "Normal", "Normal", pam50),
      group = factor(group, levels = c("HER2-enriched", "Basal-like",
                                       "Luminal A", "Normal")),
      target_gene = factor(target_gene, levels = c("ERBB2", "TP53", "BRCA1"))
    )

  p_genes <- ggplot(beta_genes,
                    aes(x = pos / 1e6, y = beta, color = group)) +
    geom_point(aes(shape = sig), size = 1.5, alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE,
                linewidth = 0.8, na.rm = TRUE) +
    scale_color_manual(values = SUBTYPE_COLORS, name = "Group") +
    scale_shape_manual(
      values = c(Hypermethylated = 17, Hypomethylated = 25, NS = 16),
      name   = "DMP status"
    ) +
    scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
    facet_wrap(~ target_gene, scales = "free_x", ncol = 1) +
    labs(
      title    = "Methylation at ERBB2, TP53, BRCA1 loci",
      subtitle = "Beta value per probe · shape = DMP significance",
      x        = "Genomic position (Mb)",
      y        = "Beta value"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "#2C3E50"),
      strip.text       = element_text(color = "white", face = "bold"),
      legend.position  = "right"
    )

  ggsave(file.path(FIGURES_DIR, "04_genes_of_interest.pdf"),
         p_genes, width = 10, height = 10)
  message("  Saved: figures/methylation/04_genes_of_interest.pdf")
} else {
  message("  No probes found for ERBB2/TP53/BRCA1 after filtering — skipping plot 4")
}

# ---------------------------------------------------------------------------
# 11. EXPORT — gene-level mean beta for 07_integration.R
# ---------------------------------------------------------------------------
message("\n=== Exporting tables ===")

# Mean beta per gene of interest per sample
gene_meth_export <- beta_all |>
  filter(probe_id %in% genes_probes$probe_id) |>
  left_join(genes_probes |> select(probe_id, target_gene), by = "probe_id") |>
  group_by(case_id, tissue_type, target_gene) |>
  summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop") |>
  filter(tissue_type == "Tumor") |>
  select(-tissue_type) |>
  pivot_wider(names_from = target_gene,
              values_from = mean_beta,
              names_prefix = "meth_") |>
  left_join(METADATA, by = "case_id")

write_tsv(gene_meth_export,
          file.path(FIGURES_DIR, "gene_methylation_summary.tsv"))
message("  Saved: figures/methylation/gene_methylation_summary.tsv")

# ---------------------------------------------------------------------------
# 12. SESSION INFO
# ---------------------------------------------------------------------------
sink(file.path(FIGURES_DIR, "session_info.txt"))
sessionInfo()
sink()

message("\n=== 06_methylation.R complete ===")
message("Outputs in figures/methylation/:")
message("  01_volcano.pdf               - DMP overview, Tumor vs Normal")
message("  02_chr17_methylation.pdf     - chr17 beta value profile")
message("  03_heatmap_top50.pdf         - top 50 DMPs, beta value heatmap")
message("  04_genes_of_interest.pdf     - ERBB2 / TP53 / BRCA1 probe-level")
message("  dmp_results.tsv              - full limma results table")
message("  gene_methylation_summary.tsv - input for 07_integration.R")
