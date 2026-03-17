# Integrative Multi-Omics Analysis of Chromosome 17 in TCGA Breast Cancer (Ongoing)

Integrative analysis of **Copy Number Variation**, **RNA-seq expression**, and **DNA Methylation** across three independent but connected projects, using open-access derived data from TCGA-BRCA.

Biological focus: **ERBB2**, **TP53**, and **BRCA1** on chromosome 17.

> Wet lab → dry lab transition project. Adapted from the [SARS-CoV-2 Genomic Surveillance Pipeline](../sars-cov2-surveillance).  
> Reference: *Comprehensive molecular portraits of human breast tumours* — TCGA Network, Nature 2012.

---

## Project structure

```
tcga-brca-multiomics/
├── data/
│   ├── raw/
│   │   ├── manifest_clean.txt        # 44-file GDC manifest (filtered)
│   │   ├── master_sample_table.tsv   # maps file_id → case, data type, tissue
│   │   └── downloads/                # temporary gdc-client output
│   ├── cnv/
│   │   ├── tumor/                    # Masked Copy Number Segment — Primary Tumor
│   │   └── normal/                   # Masked Copy Number Segment — Normal
│   ├── rna/
│   │   ├── tumor/                    # STAR gene counts TSV — Primary Tumor
│   │   └── normal/                   # STAR gene counts TSV — Normal (4 cases)
│   └── methylation/
│       ├── tumor/                    # Methylation Beta Value TXT — Primary Tumor
│       └── normal/                   # Methylation Beta Value TXT — Normal (4 cases)
├── results/
│   ├── cnv/                          # CNV analysis outputs
│   ├── rna/                          # DESeq2 outputs
│   └── methylation/                  # minfi / ChAMP outputs
├── figures/                          # Final plots (PDF/PNG)
├── scripts/
│   ├── 00_download.sh                # GDC download + file organization
│   ├── 01_qc.sh                      # FastQC + MultiQC (recycled)
│   ├── 02_alignment.sh               # BWA-MEM + samtools markdup *
│   ├── 03_cnvkit.sh                  # CNVkit batch/call/export *
│   ├── 04_cnv_analysis.R             # CNV visualization + stats
│   ├── 05_rnaseq.R                   # DESeq2 normalization + DEG  [future]
│   ├── 06_methylation.R              # minfi / ChAMP + DMR analysis [future]
│   └── 07_integration.R              # Multi-omics integration      [future]
├── environment.yml                   # Conda environment (all 3 projects)
└── README.md
```

> \* Scripts `02_alignment.sh` and `03_cnvkit.sh` are included for completeness and
> future use with controlled-access FASTQs (e.g. dbGaP / PhD cluster).
> The current open-access workflow starts at `00_download.sh` and proceeds
> directly to the analysis scripts (`04`, `05`, `06`).

---

## Pipeline overview

### Current workflow — open-access derived data

```
GDC Portal (open access)
  │
  ├── Masked Copy Number Segment (.txt)     ← pre-computed by GDC/TCGA
  ├── Gene Expression Quantification (.tsv) ← STAR counts, GRCh38
  └── Methylation Beta Value (.txt)         ← Illumina 450k array
        │
        ▼ 00_download.sh
  gdc-client → organized into data/cnv | rna | methylation
        │
        ├──▶ 04_cnv_analysis.R
        │    chr17 segments · ERBB2/TP53/BRCA1 heatmap · Gviz locus plots
        │
        ├──▶ 05_rnaseq.R              [future]
        │    DESeq2 normalization · DEG tumor vs normal · volcano plots
        │
        ├──▶ 06_methylation.R         [future]
        │    minfi QC · beta/M-value · DMR detection · methylation heatmap
        │
        └──▶ 07_integration.R         [future]
             CNV × expression correlation · methylation × expression
             multi-omics heatmap · MOFA+ factor analysis
```

### Future workflow — controlled-access FASTQs (dbGaP)

```
FASTQ (SRA / GDC controlled access)
  │
  ▼ 01_qc.sh          fastp → FastQC → MultiQC
  ▼ 02_alignment.sh   BWA-MEM → samtools markdup → .dedup.bam
  ▼ 03_cnvkit.sh      CNVkit batch → .cnr / .cns / .call.cns
  ▼ 04_cnv_analysis.R (same as above)
```

---

## Sample cohort

9 cases selected from the TCGA-BRCA 2012 paper (Supplementary Table 1).
Selection criteria: all three data types available in open access + biological
diversity across expected chr17 alteration patterns.

| Case ID | PAM50 subtype | HER2 status | Expected chr17 event |
|---------|--------------|-------------|----------------------|
| TCGA-B6-A0I9 | HER2-enriched | Positive | ERBB2 amplification |
| TCGA-BH-A18R | HER2-enriched | Positive | ERBB2 amplification |
| TCGA-BH-A0DZ | HER2-enriched | Positive | ERBB2 amplification |
| TCGA-A1-A0SK | Basal-like | Negative | TP53 deletion / BRCA1 loss |
| TCGA-A2-A0CM | Basal-like | Negative | TP53 deletion / BRCA1 loss |
| TCGA-BH-A18V | Basal-like | Negative | TP53 deletion / BRCA1 loss |
| TCGA-BH-A18Q | Basal-like | Negative | TP53 deletion / BRCA1 loss |
| TCGA-A2-A0CU | Luminal A | Negative | Biological control |
| TCGA-AR-A0TR | Luminal A | Negative | Biological control |

### Data availability per case

| Case ID | CNV tumor | CNV normal | RNA tumor | RNA normal | Meth tumor | Meth normal |
|---------|-----------|------------|-----------|------------|------------|-------------|
| TCGA-A1-A0SK | ✅ | ✅ | ✅ | — | ✅ | — |
| TCGA-A2-A0CM | ✅ | ✅ | ✅ | — | ✅ | — |
| TCGA-A2-A0CU | ✅ | ✅ | ✅ | — | ✅ | — |
| TCGA-AR-A0TR | ✅ | ✅ | ✅ | — | ✅ | — |
| TCGA-B6-A0I9 | ✅ | ✅ | ✅ | — | ✅ | — |
| TCGA-BH-A0DZ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| TCGA-BH-A18Q | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| TCGA-BH-A18R | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| TCGA-BH-A18V | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

> Normal for CNV = Peripheral Blood (sample suffix `-10A`).  
> Normal for RNA-seq and methylation = Solid Tissue Normal (`-11A`),
> available only for 4 cases. For the remaining 5, tumor-only analysis applies.

### Manifest filtering

The original GDC cart contained 48 files. Four were removed before download:

- **3 files** from `TCGA-BH-A18V-06A` (metastatic sample) — kept only Primary Tumor `-01A`
- **1 file** from `TCGA-BH-A0DZ-11A` (duplicate CNV normal, solid tissue) — kept Peripheral Blood `-10A`

---

## Key genes

| Gene | Chr17 locus (hg38) | Expected alteration in BRCA |
|------|--------------------|-----------------------------|
| ERBB2 | 17q12 (39.7 Mb) | Amplification (~15–20% of tumors) |
| TP53 | 17p13.1 (7.7 Mb) | Deletion / LOH (~35%) |
| BRCA1 | 17q21.31 (43.0 Mb) | Deletion / LOH (~15%) |

---

## TCGA data access

| Data type | Format | Access | Download |
|-----------|--------|--------|----------|
| Masked Copy Number Segment | TXT | Open ✅ | `gdc-client` |
| Gene Expression Quantification | TSV | Open ✅ | `gdc-client` |
| Methylation Beta Value | TXT | Open ✅ | `gdc-client` |
| WES / RNA-seq BAM (unaligned) | BAM | Controlled (dbGaP) | `gdc-client` + token |

All files were obtained from the [GDC Data Portal](https://portal.gdc.cancer.gov/)
under project **TCGA-BRCA**, open access tier.

---

## Quickstart

```bash
# 1. Create and activate conda environment
conda env create -f environment.yml
conda activate tcga-brca

# 2. Place manifest and sample table in data/raw/
cp manifest_clean.txt master_sample_table.tsv data/raw/

# 3. Download and organize files (~300 MB total)
bash scripts/00_download.sh

# 4. CNV analysis
Rscript scripts/04_cnv_analysis.R

# 5. RNA-seq analysis  [future]
Rscript scripts/05_rnaseq.R

# 6. Methylation analysis  [future]
Rscript scripts/06_methylation.R

# 7. Multi-omics integration  [future]
Rscript scripts/07_integration.R
```

---

## Environment

All dependencies are managed via Conda. See `environment.yml` for the full
specification covering all three projects plus the future integration module.

```bash
# Create
conda env create -f environment.yml

# Export locked versions after creation (for full reproducibility)
conda env export > environment_locked.yml

# Recreate on a new machine
conda env create -f environment_locked.yml
```

---

## Notes

- **Reference genome**: all GDC open-access files are aligned to **GRCh38/hg38**. Coordinates in the TCGA 2012 paper use hg19 — do not mix liftover positions with these files.
- **CNV data type**: `Masked Copy Number Segment` is derived from SNP array (Affymetrix GenomeWideSNP_6), not WES. It provides genome-wide segmented log2 ratios suitable for gene-level CNV analysis without requiring alignment.
- **RNA-seq counts**: files use the **STAR augmented** format with four count columns (unstranded, stranded_first, stranded_second, TPM). Use `unstranded` counts as input to DESeq2.
- **Scalability**: scripts iterate over files in `data/` directories. Adding more samples requires only updating the manifest and re-running — no script changes needed.

---
## Project Status
Active development

- ✅ CNV analysis completed
- ✅ RNA-seq analysis in completed
- ✅ Methylation analysis in completed
- 🔄 Multi-omics integration planned
- 🔄 Results discussion planned

---
## Related projects

- [SARS-CoV-2 Genomic Surveillance](../sars-cov2-surveillance) — variant calling pipeline (BWA-MEM → BCFtools → SnpEff) that originated the QC and alignment scripts reused here.
- Future: RNA-seq differential expression (`05_rnaseq.R`)
- Future: DNA methylation analysis (`06_methylation.R`)
- Future: Multi-omics integration (`07_integration.R`) — connects all three projects
