#############################################
# PCA ANALYSIS - Coffea arabica (Catuai/CR95)
# Corrected script: paths adjusted and unnecessary
# packages removed (only DESeq2 and ggplot2 are used)
#############################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(DESeq2)
library(ggplot2)

# --- ADJUST THIS PATH TO YOUR OWN WORKING FOLDER ---
setwd("~/Desktop/(7mo) SÉPTIMO SEMESTRE/FUNCTIONAL GENOMICS AND SYNTHETIC BIOLOGY")
outpathcount <- "~/Desktop/(7mo) SÉPTIMO SEMESTRE/FUNCTIONAL GENOMICS AND SYNTHETIC BIOLOGY/"

#############################################
# DATA LOADING
#############################################

counts <- read.delim("Matrix_hisat2.txt", stringsAsFactors = FALSE, row.names = 1)
counts <- counts[rowSums(counts) > 0, ]  # Remove all genes that have no expression

# Fixes the mistyped sample name in the matrix (Cax7 -> Ca7x),
# consistent with the main script. Without this, that sample
# would not find a match in the metadata and would break the
# per-group filtering.
colnames(counts)[colnames(counts) == "Cax7"] <- "Ca7x"

sample.id <- colnames(counts)

# Cultivar metadata per sample
metadata <- data.frame(
  sample    = c("Ca1x", "Ca2x", "Ca3x", "Ca4x", "Ca5x", "Ca6x", "Ca7x", "Ca10x",
                "CaC1", "CaC2", "CaC3", "CaC6", "CaC7", "CaC9", "CaC10"),
  cultivar  = c("Catuai", "Catuai", "Catuai", "Catuai", "Catuai", "CR95", "CR95", "CR95",
                "Catuai", "Catuai", "Catuai", "CR95", "CR95", "CR95", "CR95"),
  stringsAsFactors = FALSE
)

match.idx <- match(sample.id, metadata$sample)

# Early validation: if any sample doesn't match, stop with a
# clear message instead of letting the error appear several steps later.
if (anyNA(match.idx)) {
  stop("Samples with no match in metadata: ",
       paste(sample.id[is.na(match.idx)], collapse = ", "))
}

# Cultivar identifier ("A" = Catuai, "B" = CR95)
cultivar.prefix <- ifelse(metadata$cultivar[match.idx] == "Catuai", "A", "B")
colnames(counts) <- paste0(cultivar.prefix, "_", sample.id)

cultivar <- as.factor(sapply(strsplit(colnames(counts), split = "_"), `[[`, 1))
trt <- as.factor(sub(".*([xC])\\d*$", "\\1", sapply(strsplit(colnames(counts), split = "_"), `[[`, 2)))

names(cultivar) <- colnames(counts)
names(trt) <- colnames(counts)

coldata <- data.frame(cultivar = cultivar, trt = trt)
rownames(coldata) <- colnames(counts)

#############################################
# FILTERING AND NORMALIZATION
#############################################

dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts)), colData = coldata, design = ~1)

cpm <- fpm(dds, robust = FALSE)
mat.eval <- cpm > 1

aux.eval <- paste0(cultivar, "_", trt)
aux.eval <- as.character(aux.eval)

keep <- array(FALSE, nrow(mat.eval))
cutoff <- ceiling(table(aux.eval) / 2)

for (i in unique(aux.eval)) {
  # drop = FALSE prevents R from collapsing the selection to a vector when
  # a group has only one sample, which would break rowSums().
  keep <- keep | (rowSums(mat.eval[, which(aux.eval == i), drop = FALSE]) >= cutoff[i])
}

cat("Proportion of genes retained after filtering:", round(sum(keep) / length(keep), 3), "\n")

dds.1 <- dds[keep, ]
dds.1 <- estimateSizeFactors(dds.1)
vst.1 <- vst(dds.1)

# Standardized matrix (z-score per gene)
std.mat <- t(scale(t(assay(vst.1))))

#############################################
# PCA
#############################################

pca <- princomp(std.mat)
barras <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)

cat("Variance explained by the first 5 components (%):\n")
print(round(barras[1:5], 1))

barplot(barras, las = 2, main = "Variance explained by PC", ylab = "% variance")

# Save numeric PCA results
write.csv(pca$scores, file = paste0(outpathcount, "PCA_scores.csv"))
write.csv(pca$loadings, file = paste0(outpathcount, "PCA_loadings.csv"))

#############################################
# EXPLORATORY VISUALIZATION (base R) - PC1 to PC5, all combinations
#############################################

col.trt <- as.character(trt)
col.trt[col.trt == "C"] <- "blue"
col.trt[col.trt == "x"] <- "orange"

pdf(paste0(outpathcount, "PCA_plots_loadings.pdf"), width = 8, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in 1:4) {
  for (j in (i + 1):5) {
    plot(pca$loadings[, c(i, j)], main = "PCA Loadings",
         xlab = paste0("PC ", i, " (", round(barras[i], 1), "%)"),
         ylab = paste0("PC ", j, " (", round(barras[j], 1), "%)"),
         type = "n")
    text(pca$loadings[, c(i, j)], as.character(cultivar), col = col.trt)
    legend("bottom", legend = c("C", "x"), fill = c("blue", "orange"), bty = "n")
  }
}
dev.off()

#############################################
# FINAL PLOT FOR THE REPORT (PC1 vs PC2, ggplot2)
#############################################

# IMPORTANT: princomp() was run on std.mat (genes in rows,
# samples in columns), so SAMPLES are represented in
# pca$loadings (15 rows), not in pca$scores (33,583 rows = genes).
pca_df <- data.frame(
  PC1 = pca$loadings[, 1],
  PC2 = pca$loadings[, 2],
  cultivar = as.character(cultivar),
  treatment = ifelse(as.character(trt) == "C", "saline", "xylella")
)

# Rename cultivar for the report (A -> Catuai, B -> CR95)
pca_df$cultivar <- ifelse(pca_df$cultivar == "A", "Catuai", "CR95")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment, shape = cultivar)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c(saline = "steelblue", xylella = "firebrick")) +
  labs(
    title = "PCA of gene expression (VST-scaled)",
    x = paste0("PC1 (", round(barras[1], 1), "%)"),
    y = paste0("PC2 (", round(barras[2], 1), "%)"),
    color = "Treatment", shape = "Cultivar"
  ) +
  theme_bw()

print(p_pca)
ggsave(paste0(outpathcount, "PCA_PC1_PC2.pdf"), plot = p_pca, width = 7, height = 6)
