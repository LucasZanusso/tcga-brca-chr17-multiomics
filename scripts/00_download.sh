#!/usr/bin/env bash
# =============================================================================
# 00_download.sh — Download TCGA-BRCA multi-omics data via GDC Data Transfer Tool
# Project : CNV + RNA-seq + Methylation — TCGA-BRCA (9 cases, chr17 focus)
#
# Input files (place in data/raw/ before running):
#   manifest_clean.txt     — 44-file manifest (filtered from GDC cart):
#                            removed 1 metastatic sample (TCGA-BH-A18V-06A)
#                            removed 1 duplicate CNV normal (TCGA-BH-A0DZ-11A,
#                            kept peripheral blood -10A over solid tissue -11A)
#   master_sample_table.tsv — maps file_id to case_id, data type, tissue type
#
# Output structure:
#   data/cnv/tumor|normal/
#   data/rna/tumor|normal/
#   data/methylation/tumor|normal/
#
# Prerequisite:
#   conda install -c bioconda gdc-client
#   or download from: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
MANIFEST="data/raw/manifest_clean.txt"
SAMPLE_SHEET="data/raw/master_sample_table.tsv"
GDC_CLIENT="gdc-client"
THREADS=4                           # parallel connections for gdc-client
DOWNLOAD_DIR="data/raw/downloads"   # temporary download destination

# Final directories per project
CNV_TUMOR_DIR="data/cnv/tumor"
CNV_NORMAL_DIR="data/cnv/normal"
RNA_TUMOR_DIR="data/rna/tumor"
RNA_NORMAL_DIR="data/rna/normal"
METH_TUMOR_DIR="data/methylation/tumor"
METH_NORMAL_DIR="data/methylation/normal"

mkdir -p "$DOWNLOAD_DIR" \
         "$CNV_TUMOR_DIR"  "$CNV_NORMAL_DIR" \
         "$RNA_TUMOR_DIR"  "$RNA_NORMAL_DIR" \
         "$METH_TUMOR_DIR" "$METH_NORMAL_DIR"

# ---------------------------------------------------------------------------
# 1. DOWNLOAD via gdc-client
#    -m : manifest file
#    -d : output directory
#    -n : number of parallel connections
#    MD5 verification is enabled by default; remove --no-verify-md5 to keep it
# ---------------------------------------------------------------------------
N_FILES=$(tail -n +2 "$MANIFEST" | wc -l)
echo "=== [$(date +%T)] Starting download of ${N_FILES} files ==="
echo "    Temporary destination: ${DOWNLOAD_DIR}"

$GDC_CLIENT download \
    -m "$MANIFEST" \
    -d "$DOWNLOAD_DIR" \
    -n "$THREADS"

echo "=== Download complete ==="

# ---------------------------------------------------------------------------
# 2. ORGANIZE files into project directories
#    gdc-client creates one subfolder per file_id — we flatten and rename
#    using master_sample_table.tsv to route each file to the right place
#    Final filename format: {case_id}.{original_filename}
#    Example: TCGA-A1-A0SK.BONZE_p_TCGAb56...nocnv_grch38.seg.v2.txt
# ---------------------------------------------------------------------------
echo ""
echo "=== Organizing files ==="

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r file_id filename data_category data_type case_id sample_id tissue_type; do

    # gdc-client places each file under downloads/{file_id}/{filename}
    src="${DOWNLOAD_DIR}/${file_id}/${filename}"
    if [[ ! -f "$src" ]]; then
        echo "  WARNING: file not found: ${src}"
        continue
    fi

    # Route to correct directory based on data category and tissue type
    case "$data_category" in
        "Copy Number Variation")
            dest_dir=$([[ "$tissue_type" == "Tumor" ]] && echo "$CNV_TUMOR_DIR" || echo "$CNV_NORMAL_DIR")
            ;;
        "Transcriptome Profiling")
            dest_dir=$([[ "$tissue_type" == "Tumor" ]] && echo "$RNA_TUMOR_DIR" || echo "$RNA_NORMAL_DIR")
            ;;
        "DNA Methylation")
            dest_dir=$([[ "$tissue_type" == "Tumor" ]] && echo "$METH_TUMOR_DIR" || echo "$METH_NORMAL_DIR")
            ;;
        *)
            echo "  WARNING: unknown data category '${data_category}', skipping."
            continue
            ;;
    esac

    dest_file="${dest_dir}/${case_id}.${filename}"
    cp "$src" "$dest_file"
    echo "  OK  ${case_id} | ${data_category} | ${tissue_type} → $(basename "$dest_file")"
done

# ---------------------------------------------------------------------------
# 3. SUMMARY — count files per directory
# ---------------------------------------------------------------------------
echo ""
echo "=== Summary ==="
for dir in "$CNV_TUMOR_DIR"  "$CNV_NORMAL_DIR" \
           "$RNA_TUMOR_DIR"  "$RNA_NORMAL_DIR" \
           "$METH_TUMOR_DIR" "$METH_NORMAL_DIR"; do
    n=$(find "$dir" -type f | wc -l)
    echo "  ${dir}: ${n} file(s)"
done

echo ""
echo "=== 00_download.sh complete ==="
echo "Next step: run analysis scripts per project (cnv / rna / methylation)"
