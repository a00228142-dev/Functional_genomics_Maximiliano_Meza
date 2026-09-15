# Coffea arabica × Xylella fastidiosa — Transcriptomic Analysis

Differential expression and functional enrichment analysis comparing two *Coffea arabica* cultivars (Catuai and CR95) under *Xylella fastidiosa* infection versus saline control. The pipeline covers quality control through HISAT2 alignment (external, see report), PCA, differential expression (DESeq2), functional enrichment (hypergeometric test and CAMERA), coexpression network and module detection, and promoter motif analysis (MEME/TomTom).

Full methodological detail (thresholds, parameters, experimental design) is described in the accompanying report; this repository provides the exact R code used to produce those results.

## Repository structure
