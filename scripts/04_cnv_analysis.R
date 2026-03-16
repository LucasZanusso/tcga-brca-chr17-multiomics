# =============================================================================
# 04_cnv_analysis.R — CNV Analysis from GDC Masked Copy Number Segments
# Project : TCGA-BRCA Multi-Omics (chr17 focus)
# Input   : data/cnv/tumor/*.txt  — Masked Copy Number Segment (GDC open access)
#           data/cnv/normal/*.txt — same format, blood-derived normal
#           data/raw/master_sample_table.tsv
# Output  : figures/cnv/
#
# Data format (GDC Masked Copy Number Segment):
#   GDC_Aliquot | Chromosome | Start | End | Num_Probes | Segment_Mean
#   Segment_Mean = log2(copy_ratio) — positive = gain, negative = loss
#
# File naming convention (set by 00_download.sh):
#   {case_id}.{original_filename}.txt
#   Example: TCGA-A1-A0SK.SHAWM_p_TCGAb72...nocnv_grch38.seg.v2.txt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

# Bioconductor — install once with:
# BiocManager::install(c("Gviz", "GenomicRanges",
#   "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
CNV_TUMOR_DIR  <- "data/cnv/tumor"
CNV_NORMAL_DIR <- "data/cnv/normal"
SAMPLE_TABLE   <- "data/raw/master_sample_table.tsv"
FIGURES_DIR    <- "figures/cnv"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# PAM50 subtype and HER2 status for each case (from Supplementary Table 1)
# Used for coloring and grouping throughout the analysis
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

# Subtype color palette — reused in 05_rnaseq.R and 07_integration.R
SUBTYPE_COLORS <- c(
  "HER2-enriched" = "#E84855",
  "Basal-like"    = "#3B4CC0",
  "Luminal A"     = "#59A14F"
)

# Genes of interest — hg38 coordinates
GENE_COORDS <- tribble(
  ~gene,   ~chr,    ~start,   ~end,     ~plot_start, ~plot_end,
  "ERBB2", "chr17", 39687914, 39730426, 39200000,    40200000,
  "TP53",  "chr17",  7668402,  7687538,  7200000,     8200000,
  "BRCA1", "chr17", 43044295, 43125483, 42500000,    43700000
)

# Segment_Mean thresholds for CN state classification
GAIN_THRESHOLD <- 0.3    # log2 ratio > 0.3  → gain
LOSS_THRESHOLD <- -0.3   # log2 ratio < -0.3 → loss
AMP_THRESHOLD  <- 1.0    # log2 ratio > 1.0  → high-level amplification

CN_COLORS <- c(
  Amplification = "#B40426",
  Gain          = "#F4A582",
  Neutral       = "#F7F7F7",
  Loss          = "#4393C3"
)

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------
message("=== Loading CNV segment files ===")

# Extract case_id from filename — 00_download.sh prefixes files as {case_id}.{name}
extract_case_id <- function(filepath) {
  basename(filepath) |> str_extract("^TCGA-[A-Z0-9]+-[A-Z0-9]+")
}

# Load one segment file and tag with case_id and tissue type
load_seg_file <- function(filepath, tissue_type) {
  case_id <- extract_case_id(filepath)
  read_tsv(filepath,
           col_types = cols(
             GDC_Aliquot  = col_character(),
             Chromosome   = col_character(),
             Start        = col_double(),
             End          = col_double(),
             Num_Probes   = col_double(),
             Segment_Mean = col_double()
           ),
           show_col_types = FALSE) |>
    mutate(
      case_id     = case_id,
      tissue_type = tissue_type,
      # Standardise to UCSC chromosome names (chr1, chr17, etc.)
      Chromosome  = if_else(str_starts(Chromosome, "chr"),
                            Chromosome, paste0("chr", Chromosome))
    )
}

tumor_files  <- list.files(CNV_TUMOR_DIR,  pattern = "\\.txt$", full.names = TRUE)
normal_files <- list.files(CNV_NORMAL_DIR, pattern = "\\.txt$", full.names = TRUE)

if (length(tumor_files) == 0)
  stop("No tumor CNV files found in ", CNV_TUMOR_DIR,
       ". Run 00_download.sh first.")

seg_tumor  <- map_dfr(tumor_files,  load_seg_file, tissue_type = "Tumor")
seg_normal <- map_dfr(normal_files, load_seg_file, tissue_type = "Normal")
seg_all    <- bind_rows(seg_tumor, seg_normal) |>
  left_join(METADATA, by = "case_id")

message("  Tumor samples  : ", n_distinct(seg_tumor$case_id))
message("  Normal samples : ", n_distinct(seg_normal$case_id))
message("  Total segments : ", nrow(seg_all))

# ---------------------------------------------------------------------------
# 2. CLASSIFY SEGMENTS
# ---------------------------------------------------------------------------
seg_all <- seg_all |>
  mutate(
    cn_state = case_when(
      Segment_Mean >= AMP_THRESHOLD  ~ "Amplification",
      Segment_Mean >= GAIN_THRESHOLD ~ "Gain",
      Segment_Mean <= LOSS_THRESHOLD ~ "Loss",
      TRUE                           ~ "Neutral"
    ),
    cn_state = factor(cn_state,
                      levels = c("Amplification", "Gain", "Neutral", "Loss"))
  )

# Sample order: HER2-enriched → Basal-like → Luminal A
SAMPLE_ORDER <- METADATA |>
  arrange(factor(pam50, levels = c("HER2-enriched", "Basal-like", "Luminal A"))) |>
  pull(case_id)

# ---------------------------------------------------------------------------
# 3. CHROMOSOME 17 OVERVIEW
# ---------------------------------------------------------------------------
message("\n=== Plot 1: chr17 overview ===")

chr17_tumor <- seg_all |>
  filter(Chromosome == "chr17", tissue_type == "Tumor") |>
  mutate(case_id = factor(case_id, levels = SAMPLE_ORDER))

p_chr17 <- ggplot(chr17_tumor,
                  aes(x     = Start / 1e6, xend = End / 1e6,
                      y     = Segment_Mean, yend = Segment_Mean,
                      color = cn_state)) +
  geom_segment(linewidth = 2.2, alpha = 0.9) +
  geom_hline(yintercept = 0,
             linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = c(GAIN_THRESHOLD, LOSS_THRESHOLD),
             linetype = "dotted", color = "grey70", linewidth = 0.3) +
  # Gene position markers
  geom_vline(data = GENE_COORDS,
             aes(xintercept = (start + end) / 2 / 1e6),
             color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(data = GENE_COORDS,
            aes(x = (start + end) / 2 / 1e6, y = 2.3, label = gene),
            inherit.aes = FALSE, size = 2.8, fontface = "bold", vjust = 0) +
  scale_color_manual(values = CN_COLORS, name = "CN state") +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  coord_cartesian(ylim = c(-2.5, 2.7)) +
  facet_wrap(~ case_id, ncol = 1, strip.position = "right") +
  labs(
    title    = "Chromosome 17 — Copy Number Segments",
    subtitle = "TCGA-BRCA · 9 Primary Tumors · Masked CNV Segment (SNP array, hg38)",
    x        = "Position (chr17)",
    y        = "Segment Mean (log2 ratio)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background   = element_rect(fill = "#2C3E50"),
    strip.text         = element_text(color = "white", size = 7),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom"
  )

ggsave(file.path(FIGURES_DIR, "01_chr17_overview.pdf"),
       p_chr17, width = 12,
       height = 1.5 * n_distinct(chr17_tumor$case_id) + 2)
message("  Saved: figures/cnv/01_chr17_overview.pdf")

# ---------------------------------------------------------------------------
# 4. GENE-LEVEL HEATMAP
#    Weighted mean Segment_Mean across each gene body
# ---------------------------------------------------------------------------
message("\n=== Plot 2: gene-level heatmap ===")

gene_cnv <- seg_tumor |>
  left_join(METADATA, by = "case_id") |>
  cross_join(GENE_COORDS |> dplyr::rename(g_start = start, g_end = end)) |>
  filter(Chromosome == chr, Start <= g_end, End >= g_start) |>
  mutate(
    overlap = pmin(End, g_end) - pmax(Start, g_start),
    weight  = overlap / (g_end - g_start)
  ) |>
  group_by(case_id, pam50, gene) |>
  summarise(mean_log2 = weighted.mean(Segment_Mean, w = weight),
            .groups = "drop") |>
  mutate(
    gene    = factor(gene, levels = c("ERBB2", "TP53", "BRCA1")),
    case_id = factor(case_id, levels = SAMPLE_ORDER)
  )

p_heatmap <- ggplot(gene_cnv,
                    aes(x = gene, y = case_id, fill = mean_log2)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = round(mean_log2, 2)), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(
    low      = "#4393C3",
    mid      = "white",
    high     = "#B40426",
    midpoint = 0,
    limits   = c(-2, 2),
    oob      = scales::squish,
    name     = "log2 ratio"
  ) +
  labs(
    title    = "Gene-level Copy Number — ERBB2 · TP53 · BRCA1",
    subtitle = "Weighted mean log2 ratio across gene body · TCGA-BRCA",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9),
    panel.grid  = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "02_gene_heatmap.pdf"),
       p_heatmap, width = 6,
       height = max(4, n_distinct(gene_cnv$case_id) * 0.55 + 1.5))
message("  Saved: figures/cnv/02_gene_heatmap.pdf")

# ---------------------------------------------------------------------------
# 5. GVIZ LOCUS ZOOM
# ---------------------------------------------------------------------------
message("\n=== Plot 3: Gviz locus zoom (ERBB2, TP53, BRCA1) ===")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

make_locus_plot <- function(gene_name, chr,
                            gene_start, gene_end,
                            plot_start, plot_end,
                            segments_df) {

  gaxis <- GenomeAxisTrack(col = "black", fontsize = 9)

  gene_track <- GeneRegionTrack(
    txdb,
    chromosome           = chr,
    start                = plot_start,
    end                  = plot_end,
    name                 = "Genes",
    showId               = TRUE,
    transcriptAnnotation = "symbol",
    collapseTranscripts  = "meta",
    fill                 = "#2C3E50",
    col                  = NA,
    fontsize.group       = 9
  )

  # One track per sample, colored by PAM50
  sample_tracks <- map(SAMPLE_ORDER, function(cid) {
    color <- SUBTYPE_COLORS[METADATA$pam50[METADATA$case_id == cid]]
    segs  <- segments_df |>
      filter(case_id == cid, Chromosome == chr,
             tissue_type == "Tumor",
             Start <= plot_end, End >= plot_start) |>
      mutate(Start = pmax(Start, plot_start),
             End   = pmin(End,   plot_end))

    if (nrow(segs) == 0) return(NULL)

    gr <- GRanges(seqnames = chr,
                  ranges   = IRanges(segs$Start, segs$End),
                  score    = segs$Segment_Mean)

    DataTrack(
      range            = gr,
      data             = "score",
      name             = str_remove(cid, "TCGA-"),
      type             = "histogram",
      fill             = color,
      col              = color,
      ylim             = c(-2.5, 2.5),
      yTicksAt         = c(-2, -1, 0, 1, 2),
      baseline         = 0,
      col.baseline     = "grey40",
      lwd.baseline     = 0.8,
      fontsize         = 8,
      background.title = "#2C3E50",
      fontcolor.title  = "white"
    )
  }) |> compact()

  highlight <- HighlightTrack(
    trackList    = list(gene_track),
    start        = gene_start,
    end          = gene_end,
    chromosome   = chr,
    col          = "firebrick",
    fill         = "#FFE4E1",
    inBackground = TRUE
  )

  pdf_path <- file.path(FIGURES_DIR,
                        sprintf("03_locus_%s.pdf", gene_name))
  pdf(pdf_path, width = 11,
      height = 2.5 + 1.2 * length(sample_tracks))
  plotTracks(
    c(list(gaxis), sample_tracks, list(highlight)),
    chromosome = chr,
    from       = plot_start,
    to         = plot_end,
    main       = sprintf("%s locus — chr17 · TCGA-BRCA (hg38)", gene_name),
    cex.main   = 1.1
  )
  dev.off()
  message("  Saved: figures/cnv/03_locus_", gene_name, ".pdf")
}

pwalk(GENE_COORDS, function(gene, chr, start, end, plot_start, plot_end) {
  tryCatch(
    make_locus_plot(gene, chr, start, end, plot_start, plot_end, seg_all),
    error = function(e)
      message("  WARNING: Gviz failed for ", gene, " — ", e$message)
  )
})

# ---------------------------------------------------------------------------
# 6. CN STATE BARPLOT — proportion of genome per state per sample
# ---------------------------------------------------------------------------
message("\n=== Plot 4: CN state distribution ===")

cn_summary <- seg_tumor |>
  left_join(METADATA, by = "case_id") |>
  mutate(
    seg_length = End - Start,
    cn_state   = case_when(
      Segment_Mean >= AMP_THRESHOLD  ~ "Amplification",
      Segment_Mean >= GAIN_THRESHOLD ~ "Gain",
      Segment_Mean <= LOSS_THRESHOLD ~ "Loss",
      TRUE                           ~ "Neutral"
    ),
    cn_state = factor(cn_state,
                      levels = c("Amplification", "Gain", "Neutral", "Loss")),
    case_id  = factor(case_id, levels = SAMPLE_ORDER)
  ) |>
  group_by(case_id, pam50, cn_state) |>
  summarise(total_bp = sum(seg_length), .groups = "drop") |>
  group_by(case_id) |>
  mutate(pct = total_bp / sum(total_bp) * 100) |>
  ungroup()

p_barplot <- ggplot(cn_summary,
                    aes(x = case_id, y = pct, fill = cn_state)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = CN_COLORS, name = "CN state") +
  scale_y_continuous(labels = scales::label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "Genome-wide CN State Distribution per Sample",
    subtitle = "TCGA-BRCA · 9 Primary Tumors",
    x        = NULL,
    y        = "% of genome"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

ggsave(file.path(FIGURES_DIR, "04_cn_state_barplot.pdf"),
       p_barplot, width = 9, height = 5)
message("  Saved: figures/cnv/04_cn_state_barplot.pdf")

# ---------------------------------------------------------------------------
# 7. EXPORT — gene-level table for 07_integration.R
# ---------------------------------------------------------------------------
message("\n=== Exporting tables ===")

gene_cnv_export <- gene_cnv |>
  dplyr::select(case_id, gene, mean_log2) |>
  pivot_wider(names_from = gene, values_from = mean_log2) |>
  left_join(METADATA, by = "case_id")

write_tsv(gene_cnv_export,
          file.path(FIGURES_DIR, "gene_cnv_summary.tsv"))
message("  Saved: figures/cnv/gene_cnv_summary.tsv")

# ---------------------------------------------------------------------------
# 8. SESSION INFO
# ---------------------------------------------------------------------------
sink(file.path(FIGURES_DIR, "session_info.txt"))
sessionInfo()
sink()

message("\n=== 04_cnv_analysis.R complete ===")
message("Outputs in figures/cnv/:")
message("  01_chr17_overview.pdf     — genome-wide chr17 view, all 9 tumors")
message("  02_gene_heatmap.pdf       — ERBB2 / TP53 / BRCA1 across samples")
message("  03_locus_ERBB2.pdf        — Gviz locus zoom with gene models")
message("  03_locus_TP53.pdf")
message("  03_locus_BRCA1.pdf")
message("  04_cn_state_barplot.pdf   — genome-wide CN state per sample")
message("  gene_cnv_summary.tsv      — input table for 07_integration.R")
