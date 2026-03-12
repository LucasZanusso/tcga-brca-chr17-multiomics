# Projeto-

cnv-analysis-tcga-brca/
├── config.sh                  ← única fonte de verdade
├── scripts/
│   ├── 00_download.sh         ← lê SAMPLES de config.sh
│   ├── 02_alignment.sh        ← aceita sample como $1 ou itera config.sh
│   ├── 03_cnvkit.sh
│   ├── 04_cnv_analysis.R
│   ├── 05_rnaseq.sh           ← futuro
│   ├── 06_deseq2.R            ← futuro
│   └── 07_integration.R       ← futuro
└── slurm/
    ├── 02_alignment.slurm     ← array job pronto, só descomentar
    └── 03_cnvkit.slurm
