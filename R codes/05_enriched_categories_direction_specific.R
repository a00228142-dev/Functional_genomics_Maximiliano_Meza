setwd("~/Desktop/(7mo) SÉPTIMO SEMESTRE/FUNCTIONAL GENOMICS AND SYNTHETIC BIOLOGY")
getwd()

outpathcount <- "~/Desktop/(7mo) SÉPTIMO SEMESTRE/FUNCTIONAL GENOMICS AND SYNTHETIC BIOLOGY/"
dir.create(outpathcount, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager:: install("edgeR")
if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")

library (pheatmap)
library(edgeR)
library(ggplot2)
library(ashr)

metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
class (metadata)

counts_hisat <- read.table ("Matrix_hisat2.txt", header = TRUE, row.names = 1,
                            sep = "\t", check.names = FALSE)

#colnames(counts_hisat)[colnames(counts_hisat) == "Cax7"] <- "Ca7x"

stopifnot(all(metadata$sample %in% colnames (counts_hisat)))
counts_hisat <- counts_hisat [, match(metadata$sample, colnames(counts_hisat))]
counts_hisat <- round(counts_hisat)

cultivars <- c("Catuai", "CR95")
exists("cultivars")
cultivars

#############     DESEQ2 analysis   ##########################
# CR95: Control vs xylella
# Catuai: Control vs xylella

if (!requireNamespace ("DESeq2", quietly = TRUE)) BiocManager::install ("DESeq2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages ("pheatmap")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages ("igraph")
if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

library(DESeq2)
library(pheatmap)
library(igraph)
library(ashr)
library(ggplot2)
library(reshape2)
library(edgeR) #Used for filterByExpr()

#DATA
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cultivars <- c("Catuai", "CR95")
sig_cutoff <- 0.05

metadata$cultivar <- factor(metadata$cultivar, levels = c("Catuai", "CR95"))
metadata$treatment <- factor(metadata$treatment, levels = c("saline", "xylella"))

table(metadata$cultivar)
table(metadata$treatment)

counts_hisat <- read.table("Matrix_hisat2.txt", header = TRUE, row.names = 1,
                           sep = "\t", check.names = FALSE)
colnames(counts_hisat)[colnames(counts_hisat) == "Cax7"] <- "Ca7x"
stopifnot(all(metadata$sample %in% colnames(counts_hisat)))
counts_hisat <- counts_hisat[, match(metadata$sample, colnames(counts_hisat))]
counts_hisat <- round(counts_hisat)

cultivar_labels <- c(Catuai = "Catuai", CR95 = "CR95")

#DESEQ2 PER CULTIVAR

# --- define your cutoffs ONCE, at the top ---
sig_cutoff <- 0.05
lfc_cutoff <- 1   # log2FC threshold; 1 = 2-fold change. Set to 0 if you want padj-only.

de_results_deseq <- list()
dds_objects <- list()
vsd_objects <- list()

for (cv in cultivars) {
  
  cat("=== DESeq2 for", cv, "===\n")
  
  metadata_cv <- metadata[metadata$cultivar == cv, ]
  counts_cv <- counts_hisat[, match(metadata_cv$sample, colnames(counts_hisat))]
  metadata_cv$treatment <- factor(metadata_cv$treatment, levels = c("saline", "xylella"))
  
  dge_cv <- DGEList(counts = counts_cv, group = metadata_cv$treatment)
  keep_cv <- filterByExpr(dge_cv)
  counts_cv_filt <- counts_cv[keep_cv, ]
  
  dds_cv <- DESeqDataSetFromMatrix(countData = counts_cv_filt,
                                   colData = metadata_cv,
                                   design = ~treatment)
  dds_cv <- DESeq(dds_cv)
  
  res_cv <- lfcShrink(dds_cv, contrast = c("treatment", "xylella", "saline"), type = "ashr")
  res_df <- as.data.frame(res_cv)
  res_df$gene <- rownames(res_df)
  
  # --- SINGLE, CENTRALIZED DEFINITION OF "SIGNIFICANT / DE GENE" ---
  # This column is what every downstream step (heatmap, network, barplot, venn)
  # will use, so all outputs stay consistent with each other and with sig_cutoff/lfc_cutoff.
  res_df$sig <- !is.na(res_df$padj) &
    res_df$padj < sig_cutoff &
    abs(res_df$log2FoldChange) > lfc_cutoff
  
  res_df <- res_df[order(res_df$padj), ]
  
  de_results_deseq[[cv]] <- res_df
  dds_objects[[cv]] <- dds_cv
  vsd_objects[[cv]] <- vst(dds_cv, blind = FALSE)
  
  write.table(res_df, file = paste0(outpathcount, "DEgenes_DESeq2_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  titulo_cv <- paste0(cultivar_labels[[cv]], ":control vs xylella")
  
  pdf(file = paste0(outpathcount, "dispersion_", cv, ".pdf"), width = 7, height = 6)
  plotDispEsts(dds_cv, main = paste("Dispersion -", titulo_cv))
  dev.off()
  
  # NOTE: fixed — this was calling plotDispEsts() twice before (bug).
  pdf(file = paste0(outpathcount, "Maplot_", cv, ".pdf"), width = 7, height = 6)
  plotMA(res_cv, main = paste("MA plot -", titulo_cv), ylim = c(-5, 5))
  dev.off()
}

# Count of DE genes per cultivar using the centralized "sig" column
# (padj < sig_cutoff AND |log2FC| > lfc_cutoff)
sapply(de_results_deseq, function(df) sum(df$sig, na.rm = TRUE))


#HEATMAP AND COEXPRESSION NETWORK, PER CULTIVAR

# IMPORTANTE: crear la carpeta de salida ANTES del bucle
dir.create(outpathcount, showWarnings = FALSE, recursive = TRUE)

network_objects <- list()

for (cv in cultivars) {
  res_df <- de_results_deseq[[cv]]
  sig_cv <- res_df[res_df$sig, ]
  
  cat(cv, "-- genes significativos:", nrow(sig_cv), "\n")
  
  if (nrow(sig_cv) == 0) {
    cat(cv, "-- sin genes significativos, se omite\n")
    next
  }
  
  title_cv <- paste0(cultivar_labels[[cv]], ": control vs xylella")
  
  vsd_mat_cv <- assay(vsd_objects[[cv]])
  mat_cv <- vsd_mat_cv[rownames(vsd_mat_cv) %in% sig_cv$gene, , drop = FALSE]
  
  # Quitar genes con varianza cero (si no, scale() genera NaN y pheatmap falla)
  row_var <- apply(mat_cv, 1, var)
  if (any(row_var == 0)) {
    cat(cv, "-- se remueven", sum(row_var == 0), "genes con varianza cero\n")
  }
  mat_cv <- mat_cv[row_var > 0, , drop = FALSE]
  sig_cv <- sig_cv[sig_cv$gene %in% rownames(mat_cv), ]
  
  if (nrow(mat_cv) == 0) {
    cat(cv, "-- no quedan genes tras filtrar varianza cero, se omite\n")
    next
  }
  
  mat_scaled_cv <- t(scale(t(mat_cv)))
  
  metadata_cv <- metadata[metadata$cultivar == cv, ]
  ann_col_cv <- data.frame(
    treatment = metadata_cv$treatment[match(colnames(mat_scaled_cv), metadata_cv$sample)]
  )
  rownames(ann_col_cv) <- colnames(mat_scaled_cv)
  
  while (!is.null(dev.list())) dev.off()
  
  # --- HEATMAP ---
  tryCatch({
    pdf(file = file.path(outpathcount, paste0("HEEATMAP_", cv, ".pdf")), width = 8, height = 10)
    cat("Title for this interaction:", title_cv, "\n")
    pheatmap(mat_scaled_cv,
             annotation_col = ann_col_cv,
             show_rownames = TRUE,
             fontsize_row = 6,
             main = paste0("Significant genes (padj < ", sig_cutoff,
                           ", |log2FC| > ", lfc_cutoff, ") - ", title_cv))
    dev.off()
    cat(cv, "-- heatmap generado correctamente\n")
  }, error = function(e) {
    while (!is.null(dev.list())) dev.off()
    message(cv, " -- ERROR generando heatmap: ", conditionMessage(e))
  })
  
  # --- RED DE COEXPRESIÓN ---
  cor_genes_cv <- cor(t(mat_cv), method = "pearson")
  threshold <- 0.8
  adj_cv <- abs(cor_genes_cv) >= threshold
  diag(adj_cv) <- FALSE
  
  g_cv <- igraph::simplify(graph_from_adjacency_matrix(adj_cv, mode = "undirected", diag = FALSE))
  V(g_cv)$direction <- ifelse(sig_cv$log2FoldChange[match(V(g_cv)$name, sig_cv$gene)] > 0, "Up", "Down")
  V(g_cv)$color <- ifelse(V(g_cv)$direction == "Up", "deeppink", "darkseagreen")
  
  network_objects[[cv]] <- g_cv
  
  while (!is.null(dev.list())) dev.off()
  
  tryCatch({
    pdf(file = file.path(outpathcount, paste0("coexpression_network_", cv, ".pdf")), width = 9, height = 9)
    plot(g_cv, vertex.label.cex = 0.5, vertex.size = 8,
         main = paste("Coexpression network -", title_cv))
    dev.off()
    cat(cv, "-- red de coexpresion generada correctamente\n")
  }, error = function(e) {
    while (!is.null(dev.list())) dev.off()
    message(cv, " -- ERROR generando red: ", conditionMessage(e))
  })
  
  cat(cv, "-- network connections:", ecount(g_cv), "\n")
}

dev.list()
while (!is.null(dev.list())) dev.off()

#COMPARATIVE BARPLOT: UP/DOWN REGULATED PER CULTIVAR
gene_counts <- do.call(rbind, lapply(cultivars, function(cv) {
  sig_cv <- de_results_deseq[[cv]][de_results_deseq[[cv]]$sig, ]
  data.frame(
    Cultivar = cultivar_labels[[cv]],
    Upregulated = sum(sig_cv$log2FoldChange > 0),
    Downregulated = sum(sig_cv$log2FoldChange < 0)
  )
}))

gene_counts

gene_counts_long <- reshape2::melt(gene_counts, id.vars = "Cultivar",
                                   variable.name = "Regulation",
                                   value.name = "Count")

p_updown <- ggplot(gene_counts_long, aes(x = Cultivar, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3) +
  scale_fill_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue")) +
  labs(title = paste0("Differentially expressed genes (padj < ", sig_cutoff,
                      ", |log2FC| > ", lfc_cutoff, ")"),
       x = "Cultivar", y = "Number of genes") +
  theme_minimal()

print(p_updown)
ggsave(paste0(outpathcount, "genes_up_downregulated_per_cultivar.pdf"), plot = p_updown,
       width = 8, height = 6)

#VENN DIAGRAM AND SHARED GENES
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

# Significant genes per cultivar (using the same centralized "sig" column)
sig_genes_list <- lapply(cultivars, function(cv) {
  res_df <- de_results_deseq[[cv]]
  res_df$gene[res_df$sig]
})
names(sig_genes_list) <- cultivar_labels[cultivars]

#ENRICHMENT ANALYSIS 

#GO ANNOTATION
if(!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if(!requireNamespace("limma", quietly = TRUE)) BiocManager:: install("limma")
if(!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db")
if(!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
library(clusterProfiler)
library(limma)
library(GO.db)
library(AnnotationDbi)

annotation <-read.table("fullAnnotation_clean.txt", header = TRUE, sep =  "\t",
                        quote = "\"", comment.char = "", fill = TRUE,
                        stringsAsFactors = FALSE, na.strings = c("-", "NA", ""))

has_go <- !is.na(annotation$GOs) & !is.na(annotation$gene_id)
gene2go <-strsplit(annotation$GOs[has_go], ",")
names(gene2go) <- annotation$gene_id[has_go]
term2gene <- split(rep(names(gene2go), lengths(gene2go)), unlist(gene2go))
term2gene <- lapply(term2gene, unique)

term2gene_df <- data.frame(
  term = rep(names(term2gene), lengths(term2gene)),
  gene = unlist(term2gene)
)

"GOID MAP"
all_go_ids <- unique(term2gene_df$term)
go_names <- suppressMessages(
  AnnotationDbi::select(GO.db, keys = all_go_ids, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
)
term2name_df <- go_names[, c("GOID", "TERM")]
colnames(term2name_df) <- c("term", "name")

log2fc_cutoff <- 1
min_set_size <- 5

hyper_results <- list()
camera_results2 <- list()

hyper_results <- list()

for (cv in cultivars) {
  
  cat("=== Hypergeometric (UP/DOWN) withFC for", cv, "===\n")
  
  res_df <- de_results_deseq[[cv]]
  
  universe_genes <- res_df$gene[!is.na(res_df$padj)]
  
  genes_up_fc <- res_df$gene[
    !is.na(res_df$padj) &
      res_df$padj < sig_cutoff &
      res_df$log2FoldChange >= log2fc_cutoff
  ]
  
  genes_down_fc <- res_df$gene[
    !is.na(res_df$padj) &
      res_df$padj < sig_cutoff &
      res_df$log2FoldChange <= -log2fc_cutoff
  ]
  
  enrich_up_fc <- if (length(genes_up_fc) >= 5) {
    enricher(
      gene = genes_up_fc,
      universe = universe_genes,
      TERM2GENE = term2gene_df,
      TERM2NAME = term2name_df,
      pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1
    )
  } else NULL
  
  enrich_down_fc <- if (length(genes_down_fc) >= 5) {
    enricher(
      gene = genes_down_fc,
      universe = universe_genes,
      TERM2GENE = term2gene_df,
      TERM2NAME = term2name_df,
      pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1
    )
  } else NULL
  
  up_file   <- file.path(outpathcount, paste0("enrich_hyper_withFC_UP_", cv, ".txt"))
  down_file <- file.path(outpathcount, paste0("enrich_hyper_withFC_DOWN_", cv, ".txt"))
  
  write.table(as.data.frame(enrich_up_fc), up_file,
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(as.data.frame(enrich_down_fc), down_file,
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  hyper_results[[cv]] <- list(withFC_UP = enrich_up_fc, withFC_DOWN = enrich_down_fc)
  
  cat(cv, "-- UP withFC:", length(genes_up_fc), "genes,",
      nrow(as.data.frame(enrich_up_fc)), "terms\n")
  cat(cv, "-- DOWN withFC:", length(genes_down_fc), "genes,",
      nrow(as.data.frame(enrich_down_fc)), "terms\n")
}


for (cv in cultivars) {
  cat("\n=== Breakdown (withFC) ", cv, "===\n")
  
  for (dir in c("UP", "DOWN")) {
    
    enrich_obj <- hyper_results[[cv]][[paste0("withFC_", dir)]]
    
    if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat(cv, "--", dir, ": sin términos enriquecidos, se omite\n")
      next
    }
    
    enrich_df <- as.data.frame(enrich_obj)
    enrich_df <- merge(enrich_df, go_names, by.x = "ID", by.y = "GOID", all.x = TRUE)
    
    cat("-- Counting by ontology (", dir, "):\n", sep = "")
    print(table(enrich_df$ONTOLOGY))
    
    cat("-- Top 10 BP (", dir, "):\n", sep = "")
    bp_df <- enrich_df[enrich_df$ONTOLOGY == "BP", ]
    bp_df <- bp_df[order(bp_df$p.adjust), c("ID", "Description", "p.adjust", "Count")]
    print(head(bp_df, 10))
    
    out_file <- file.path(outpathcount, paste0("enrich_hyper_withFC_", dir, "_withOntology_", cv, ".txt"))
    write.table(enrich_df, out_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
}
#CAMERA ANALYSIS

for (cv in cultivars) {
  cat ("===CAMERA for", cv, "===\n")
  
  metadata_cv <- metadata[metadata$cultivar == cv, ]
  counts_cv <- counts_hisat[, match(metadata_cv$sample, colnames(counts_hisat))]
  treatment <- factor(metadata_cv$treatment, levels = c("saline", "xylella"))
  
  dge <- DGEList(counts = counts_cv, group = treatment)
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~treatment)
  v <- voom(dge, design, plot = FALSE)
  
  term2gene_cv <- lapply(term2gene, function(g) intersect(g, rownames(v$E)))
  term2gene_cv <- term2gene_cv[lengths(term2gene_cv) >= min_set_size]
  idx <- lapply(term2gene_cv, function(g) match(g, rownames(v$E)))
  
  res <- camera(v, index = idx, design = design, contrast = 2)
  res$GOID <- rownames(res)
  
  term_info <- suppressMessages(
    AnnotationDbi::select(GO.db, keys = res$GOID, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
  )
  res <- merge(res, term_info, by = "GOID", all.x = TRUE)
  res$TERM[is.na(res$TERM)] <- res$GOID[is.na(res$TERM)]
  res$ONTOLOGY[is.na(res$ONTOLOGY)] <- "Unknown"
  res <- res[order(res$PValue), ]
  
  camera_results2[[cv]] <- res
  
  write.table(res, file = paste0(outpathcount, "CAMERA_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep ="\t")
  cat(cv, "-- GO terms tested: ", nrow(res), "--FDR <0.05", sum(res$FDR < sig_cutoff), "\n")
}

#TOP CATEGORIES - CAMERA
for (cv in cultivars) {
  cat("\n=== CAMERA BREAKDOWN", cv, "===\n")
  
  camera_df <- camera_results2[[cv]]
  
  cat("-- Counting by ontology (CAMERA): \n")
  print(table(camera_df$ONTOLOGY))
  
  cat("--Top 10 BP (CAMERA): \n")
  camera_bp <-camera_df[camera_df$ONTOLOGY == "BP", ]
  print(head(camera_bp[order(camera_bp$FDR), c("GOID", "TERM", "FDR", "NGenes")], 10))
}

#SIGNIFICATIVE MODULES SEARCH

module_results <- list()
min_module_size <- 3

for(cv in cultivars) {
  
  cat("\n=== Modules for", cv, "===\n")
  
  g_cv <- network_objects[[cv]]
  
  if (is.null(g_cv) || ecount(g_cv) == 0) {
    cat(cv, "-- red sin conexiones, no se pueden detectar módulos\n")
    next
  }
  
  comm_cv <- igraph::cluster_louvain(g_cv)
  
  V(g_cv)$module <- membership(comm_cv)
  
  module_sizes <- sizes(comm_cv)
  cat("Número total de módulos detectados:", length(module_sizes), "\n")
  print(module_sizes)
  
  sig_modules <- names(module_sizes)[module_sizes >= min_module_size]
  cat("Módulos con >=", min_module_size, "genes:", length(sig_modules), "\n")
  
  cat("Modularidad de la red:", modularity(comm_cv), "\n")
  
  module_results[[cv]] <- list(graph = g_cv, communities = comm_cv, sizes = module_sizes)
  
  module_df <- data.frame(gene = V(g_cv)$name, module = V(g_cv)$module)
  module_df <- module_df[module_df$module %in% sig_modules, ]
  module_df <- module_df[order(module_df$module), ]
  
  write.table(module_df, file = paste0(outpathcount, "modules_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  pdf(file = paste0(outpathcount, "network_modules_", cv, ".pdf"), width = 9, height = 9)
  plot(comm_cv, g_cv,
       vertex.label.cex = 0.5, vertex.size = 8,
       main = paste("Coexpression modules -", cultivar_labels[[cv]]))
  dev.off()
}


#PLOTS 

#HYPERGEOMETRIC ANALYSIS PLOT (WITH FC — UP/DOWN)

fc_types <- c("withFC_UP", "withFC_DOWN")
fc_labels <- c(withFC_UP   = "upregulated",
               withFC_DOWN = "downregulated")

for (cv in cultivars) {
  for (fc_type in fc_types) {
    
    enrich_obj <- hyper_results[[cv]][[fc_type]]
    
    if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat(cv, "-", fc_type, "-- sin resultados para graficar\n")
      next
    }
    
    titulo_cv <- paste0(cultivar_labels[[cv]], ": enriched categories (", fc_labels[[fc_type]], ")")
    
    # DOTPLOT
    p_dot <- dotplot(enrich_obj, showCategory = 10) + ggtitle(titulo_cv)
    print(p_dot)
    ggsave(paste0(outpathcount, "/dotplot_enrich_", cv, "_", fc_type, ".pdf"),
           plot = p_dot, width = 8, height = 8)
    
    # BARPLOT
    p_bar <- barplot(enrich_obj, showCategory = 10) + ggtitle(titulo_cv)
    print(p_bar)
    ggsave(paste0(outpathcount, "/barplot_enrich_", cv, "_", fc_type, ".pdf"),
           plot = p_bar, width = 8, height = 8)
  }
}

#CAMERA PLOT
for (cv in cultivars) {
  
  camera_df <- camera_results2[[cv]]
  camera_bp <- camera_df[camera_df$ONTOLOGY == "BP", ]
  top_camera <- head(camera_bp[order(camera_bp$FDR), ], 15)
  top_camera$TERM <- factor(top_camera$TERM, levels = rev(top_camera$TERM))
  
  titulo_cv <- paste0(cultivar_labels[[cv]], ": CAMERA - top categorías BP")
  
  p_camera <- ggplot(top_camera, aes(x = TERM, y = -log10(FDR))) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    labs(title = titulo_cv, x = NULL, y = "-log10(FDR)") +
    theme_minimal()
  
  print(p_camera)
  ggsave(paste0(outpathcount, "camera_barplot_", cv, ".pdf"), plot = p_camera, width = 9, height = 7)
}

#MODULES GRAPHIC
for (cv in cultivars) {
  g_cv <- module_results[[cv]]$graph
  comm_cv <- module_results[[cv]]$communities
  
  if (is.null(g_cv)) next
  
  set.seed(42)
  layout_cv <- layout_with_fr(g_cv)
  
  pdf(file = paste0(outpathcount, "network_modules_clean_", cv, ".pdf"), width = 8, height = 8)
  plot(comm_cv, g_cv,
       layout = layout_cv,
       vertex.label = NA,
       vertex.size = 5,
       main = paste("Coexpression modules -", cultivar_labels[[cv]]))
  dev.off()
}

#SEARCHING SIGNIFICANT GENES
ethylene_go <- names(Term(GOTERM))[grepl("ethylene", Term(GOTERM), ignore.case = TRUE)]

# Genes de CR95, significativos, anotados a ese término GO -- ordenados por padj
genes_ethylene_CR95 <- de_results_deseq$CR95[
  de_results_deseq$CR95$gene %in% unique(unlist(term2gene[ethylene_go])) &
    !is.na(de_results_deseq$CR95$padj) & de_results_deseq$CR95$padj < sig_cutoff, ]

genes_ethylene_CR95 <- genes_ethylene_CR95[order(genes_ethylene_CR95$padj), ]
genes_ethylene_CR95

