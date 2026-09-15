# Coffea arabica × Xylella fastidiosa — Transcriptomic Analysis

Differential expression and functional enrichment analysis comparing two *Coffea arabica* cultivars (Catuai and CR95) under *Xylella fastidiosa* infection versus saline control. The pipeline covers quality control through HISAT2 alignment (external, see report), PCA, differential expression (DESeq2), functional enrichment (hypergeometric test and CAMERA), coexpression network and module detection, and promoter motif analysis (MEME/TomTom).

Full methodological detail (thresholds, parameters, experimental design) is described in the accompanying report; this repository provides the exact R code used to produce those results.

## Repository structure
coffea-xylella-transcriptomics/
├── README.md
└── R codes/
    ├── 01_PCA.R
    ├── 02_DESeq2_enrichment_network_modules.R
    ├── 03_CAMERA.R
    └── 04_MEME_promoter_motifs.R

## Requirements
•R (>=4.2).

•CRAN packages: ggplot2, pheatmap, reshape2, ggVennDiagram, dplyr.

•Bioconductor packages: DESeq2, edgeR, ashr, clusterProfiler, limma, Go.db, AnnotationDbi, igraph, rtracklayer, Biostrings, Rsamtools.

•Bioconductor package memes (requires a local MEME Suite installation; see script 04 for details).

•rentrez(CRAN) for optional FASTA retrieval from NCBI.

## Input data (not included in this repository)
The following files are requires in the working directory and are not distributed here due to size/privacy
| File | Description | Used in |
|---|---|---|
|`Matrix_hisat2.txt` | Gene count matrix (genes x samples)| 01, 02, 03 |
| `metadata.txt` | Sample metadata (cultivar,treatment)| 02, 03 |
| `fullAnnotation_clean.txt`| Functional annotation (eggNOG-mapper output; gene_id, GOs columns) | 02, 03 |
| `GCF_036785885.1_..._genomic.fna`| C. arabica reference genome (NCBI RefSeq) | 04 |
| `genomic.gff` | Reference genome annotation (GFF3) | 04 |

## Execution order
The scripts are numbered in the order they should be run. Each script assumes the working directory contains the input files listed above:
1. **01_PCA.R** - Loads the full count matrix (both cultivars), filters low-expression genes per cultivar x treatment group, applies a variance-stabilizing transformation, and performs PCA to explore sample clustering prior to splitting the analysis by cultivar.
2. **02_DESeq2_enrichment_network_modules.R** - Main pipeline: per-cultivar differential expression (DESeq2), heatmaps, coexpression networks and Louvain module detection, Venn diagram of shared/exclusive DEGs, and hypergeometric functional enrichment (clusterProfiler) with GO annotation.
3. **03_CAMERA.R** - Complementary gene-set enrichment via CAMERA (limma/voom), evaluating the full expression ranking without a fold-change cutoff.
4. **04_MEME_promoter_motifs.R** - Extracts promoter sequences (1000 bp upstream / 200 bp downstream of TSS) for a set of genes of interest, runs de novo motif discovery with MEME, and compares the resulting motifs against JASPAR_plants via TomTom

## Key parameters and thresholds
| Parameter | Value | Used in |
|---|---|---|
| Significance (padj) | < 0.05 | DESeq2, hypergeometric test |
| Fold-change threshold | \|log2FC\| > 1 | Hypergeometric test |
| Minimum gene set size | ≥ 5 genes | CAMERA |
| Correlation (network edges) | \|r\| ≥ 0.8 | Coexpression network | 
| Minimum module size | ≥ 3 genes | Module detection (Louvain) |
| CAMERA significance (FDR) | <0.05 | CAMERA |

## Outputs
Each script writes its results (tables as .txt/.csv, figures as .pdf) to the working directory. No outputs/ folder is version-controlled in this repository; regenerate them by running the scripts against the input data described above. 
