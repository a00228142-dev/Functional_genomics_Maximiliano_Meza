setwd("~/Downloads/7mo semsestre /Reto_Coffe_R")

# NOTE: the original script was missing the trailing "/" here, which made
# every output file land in Diff_class/ directly (e.g. "expresCAMERA_Catuai.txt")
# instead of inside the "expres" folder. Added below.
outpathcount <- "~/Downloads/7mo semsestre /Reto_Coffe_R"
dir.create(outpathcount, showWarnings = FALSE)

options(repos = c(CRAN = "https://cloud.r-project.org"))

# install.packages() and BiocManager::install() both accept a vector of
# package names directly -- no loop needed. If you've already installed
# these once, feel free to comment these two lines out on later runs.
install.packages(c("ggplot2", "BiocManager"))

BiocManager::install(c("edgeR", "limma", "GO.db", "AnnotationDbi"))
library(edgeR)
library(limma)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

counts <- read.table("Matrix_hisat2.txt", header = TRUE, row.names = 1, sep = "\t",
                     check.names = FALSE)
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Fix the one mismatched sample name (Cax7 in the matrix vs Ca7x in the
# metadata) rather than hardcoding every column name -- safer, since it
# doesn't assume the matrix's column order.
colnames(counts)[colnames(counts) == "Cax7"] <- "Ca7x"

counts <- counts[, match(metadata$sample, colnames(counts))]
counts <- round(counts)
annotation <- read.table("fullAnnotation_clean.txt", header = TRUE, sep = "\t",
                         quote = "\"", comment.char = "", fill = TRUE,
                         stringsAsFactors = FALSE, na.strings = c("-", "NA", ""))

has_go <- !is.na(annotation$GOs)
gene2go <- strsplit(annotation$GOs[has_go], ",")
names(gene2go) <- annotation$gene_id[has_go]
term2gene <- split(rep(names(gene2go), lengths(gene2go)), unlist(gene2go))


length(term2gene)

min_set_size <- 5

metadata_Catuai <- metadata[metadata$cultivar == "Catuai", ]
counts_Catuai   <- counts[, metadata_Catuai$sample]

treatment_Catuai <- factor(metadata_Catuai$treatment, levels = c("saline", "xylella"))
treatment_Catuai
## [1] xylella xylella xylella xylella xylella saline  saline  saline 
## Levels: saline xylella

dge_Catuai <- DGEList(counts = counts_Catuai, group = treatment_Catuai)
keep_Catuai <- filterByExpr(dge_Catuai)
dge_Catuai <- dge_Catuai[keep_Catuai, , keep.lib.sizes = FALSE]
dge_Catuai <- calcNormFactors(dge_Catuai)

design_Catuai <- model.matrix(~treatment_Catuai)
v_Catuai <- voom(dge_Catuai, design_Catuai, plot = FALSE)

term2gene_Catuai <- lapply(term2gene, function(g) intersect(g, rownames(v_Catuai$E)))
term2gene_Catuai <- term2gene_Catuai[lengths(term2gene_Catuai) >= min_set_size]
idx_Catuai <- lapply(term2gene_Catuai, function(g) match(g, rownames(v_Catuai$E)))

res_Catuai <- camera(v_Catuai, index = idx_Catuai, design = design_Catuai, contrast = 2)
res_Catuai$GOID <- rownames(res_Catuai)

term_info_Catuai <- suppressMessages(
  AnnotationDbi::select(GO.db, keys = res_Catuai$GOID, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
)
res_Catuai <- merge(res_Catuai, term_info_Catuai, by = "GOID", all.x = TRUE)
res_Catuai$TERM[is.na(res_Catuai$TERM)] <- res_Catuai$GOID[is.na(res_Catuai$TERM)]
res_Catuai$ONTOLOGY[is.na(res_Catuai$ONTOLOGY)] <- "Unknown"
res_Catuai <- res_Catuai[order(res_Catuai$PValue), ]

write.table(res_Catuai, file = paste0(outpathcount, "CAMERA_Catuai.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

nrow(res_Catuai)                      # total GO terms tested

sum(res_Catuai$FDR < 0.05)            # significant at FDR < 0.05

sig_Catuai <- res_Catuai[res_Catuai$FDR < 0.05, ]
top_Catuai <- head(sig_Catuai[order(sig_Catuai$FDR), ], 15)
top_Catuai$TERM <- factor(top_Catuai$TERM, levels = rev(top_Catuai$TERM))

p_bar_Catuai <- ggplot(top_Catuai, aes(x = TERM, y = -log10(FDR), fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  labs(title = "Catuai", x = "GO term", y = "-log10(FDR)") +
  theme_bw()

p_bar_Catuai

ggsave(paste0(outpathcount, "barplot_Catuai.pdf"), p_bar_Catuai, width = 8, height = 6)

p_dot_Catuai <- ggplot(top_Catuai, aes(x = -log10(FDR), y = TERM, size = NGenes, color = Direction)) +
  geom_point() +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  labs(title = "Catuai", x = "-log10(FDR)", y = "GO term", size = "Genes in set") +
  theme_bw()

p_dot_Catuai

###################################
#CR95
metadata_CR95 <- metadata[metadata$cultivar == "CR95", ]
counts_CR95   <- counts[, metadata_CR95$sample]

treatment_CR95 <- factor(metadata_CR95$treatment, levels = c("saline", "xylella"))
treatment_CR95

dge_CR95 <- DGEList(counts = counts_CR95, group = treatment_CR95)
keep_CR95 <- filterByExpr(dge_CR95)
dge_CR95 <- dge_CR95[keep_CR95, , keep.lib.sizes = FALSE]
dge_CR95 <- calcNormFactors(dge_CR95)

design_CR95 <- model.matrix(~treatment_CR95)
v_CR95 <- voom(dge_CR95, design_CR95, plot = FALSE)

term2gene_CR95 <- lapply(term2gene, function(g) intersect(g, rownames(v_CR95$E)))
term2gene_CR95 <- term2gene_CR95[lengths(term2gene_CR95) >= min_set_size]
idx_CR95 <- lapply(term2gene_CR95, function(g) match(g, rownames(v_CR95$E)))

res_CR95 <- camera(v_CR95, index = idx_CR95, design = design_CR95, contrast = 2)
res_CR95$GOID <- rownames(res_CR95)

term_info_CR95 <- suppressMessages(
  AnnotationDbi::select(GO.db, keys = res_CR95$GOID, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
)
res_CR95 <- merge(res_CR95, term_info_CR95, by = "GOID", all.x = TRUE)
res_CR95$TERM[is.na(res_CR95$TERM)] <- res_CR95$GOID[is.na(res_CR95$TERM)]
res_CR95$ONTOLOGY[is.na(res_CR95$ONTOLOGY)] <- "Unknown"
res_CR95 <- res_CR95[order(res_CR95$PValue), ]

write.table(res_CR95, file = paste0(outpathcount, "CAMERA_CR95.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

nrow(res_CR95)                      # total GO terms tested
## [1] 6738
sum(res_CR95$FDR < 0.05)            # significant at FDR < 0.05

sig_CR95 <- res_CR95[res_CR95$FDR < 0.05, ]
top_CR95 <- head(sig_CR95[order(sig_CR95$FDR), ], 15)
top_CR95$TERM <- factor(top_CR95$TERM, levels = rev(top_CR95$TERM))

p_bar_CR95 <- ggplot(top_CR95, aes(x = TERM, y = -log10(FDR), fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  labs(title = "CR95", x = "GO term", y = "-log10(FDR)") +
  theme_bw()

p_bar_CR95

ggsave(paste0(outpathcount, "barplot_CR95.pdf"), p_bar_CR95, width = 8, height = 6)

p_dot_CR95 <- ggplot(top_CR95, aes(x = -log10(FDR), y = TERM, size = NGenes, color = Direction)) +
  geom_point() +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  labs(title = "CR95", x = "-log10(FDR)", y = "GO term", size = "Genes in set") +
  theme_bw()

p_dot_CR95

#############################

#############################
top_Catuai$cultivar <- "Catuai"
top_CR95$cultivar   <- "CR95"

comparison_terms <- unique(c(top_Catuai$GOID, top_CR95$GOID))

comparison_Catuai <- res_Catuai[res_Catuai$GOID %in% comparison_terms, ]
comparison_Catuai$cultivar <- "Catuai"
comparison_CR95 <- res_CR95[res_CR95$GOID %in% comparison_terms, ]
comparison_CR95$cultivar <- "CR95"

comparison_df <- rbind(comparison_Catuai, comparison_CR95)

p_comparison <- ggplot(comparison_df, aes(x = cultivar, y = TERM, size = NGenes,
                                          color = -log10(FDR), shape = Direction)) +
  geom_point() +
  scale_color_viridis_c() +
  labs(title = "GO terms compared across cultivars", x = "Cultivar", y = "GO term",
       size = "Genes in set", color = "-log10(FDR)") +
  theme_bw()

p_comparison

ggsave(paste0(outpathcount, "comparison_dotplot.pdf"), p_comparison, width = 9, height = 7)

composition_Catuai <- data.frame(cultivar = "Catuai", ONTOLOGY = sig_Catuai$ONTOLOGY,
                                 Direction = sig_Catuai$Direction)
composition_CR95 <- data.frame(cultivar = "CR95", ONTOLOGY = sig_CR95$ONTOLOGY,
                               Direction = sig_CR95$Direction)
composition_df <- rbind(composition_Catuai, composition_CR95)

p_composition <- ggplot(composition_df, aes(x = cultivar, fill = ONTOLOGY)) +
  geom_bar(position = "stack") +
  facet_wrap(~Direction) +
  labs(title = "Significant GO terms by ontology and direction",
       x = "Cultivar", y = "Number of significant terms") +
  theme_bw()

p_composition

ggsave(paste0(outpathcount, "ontology_composition.pdf"), p_composition, width = 8, height = 5)

overview_Catuai <- data.frame(cultivar = "Catuai", NGenes = res_Catuai$NGenes,
                              FDR = res_Catuai$FDR, Direction = res_Catuai$Direction)
overview_CR95 <- data.frame(cultivar = "CR95", NGenes = res_CR95$NGenes,
                            FDR = res_CR95$FDR, Direction = res_CR95$Direction)
overview_df <- rbind(overview_Catuai, overview_CR95)

p_overview <- ggplot(overview_df, aes(x = NGenes, y = -log10(FDR), color = Direction)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  facet_wrap(~cultivar) +
  labs(title = "All tested GO terms: set size vs. significance",
       x = "Genes in GO term", y = "-log10(FDR)") +
  theme_bw()

p_overview

ggsave(paste0(outpathcount, "overview_plot.pdf"), p_overview, width = 9, height = 5)
