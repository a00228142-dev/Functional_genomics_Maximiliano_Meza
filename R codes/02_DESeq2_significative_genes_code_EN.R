setwd("~/Desktop/(7mo) SÉPTIMO SEMESTRE/FUNCTIONAL GENOMICS AND SYNTHETIC BIOLOGY")
getwd()

outpathcount <- "~/Desktop/(7mo) SÉPTIMO SEMESTRE/FUNCTIONAL GENOMICS AND SYNTHETIC BIOLOGY/"
dir.create(outpathcount, showWarnings = FALSE)

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

colnames(counts_hisat)[colnames(counts_hisat) == "Cax7"] <- "Ca7x"

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

de_results_deseq <- list()
dds_objects <- list()
vsd_objects <- list()

for (cv in cultivars) {
  
  cat("=== DESeq2 for", cv, "===\n")
  
  metadata_cv <- metadata[metadata$cultivar == cv, ]
  counts_cv <- counts_hisat[, match (metadata_cv$sample, colnames(counts_hisat))]
  metadata_cv$treatment <- factor (metadata_cv$treatment, levels = c("saline", "xylella"))
  
  dge_cv <- DGEList(counts = counts_cv, group = metadata_cv$treatment)
  keep_cv <- filterByExpr(dge_cv)
  counts_cv_filt <- counts_cv[keep_cv, ]
  
  dds_cv <- DESeqDataSetFromMatrix(countData = counts_cv_filt,
                                   colData = metadata_cv,
                                   design = ~treatment)
  dds_cv <- DESeq(dds_cv)
  
  res_cv <- lfcShrink (dds_cv, contrast = c("treatment", "xylella", "saline"), type = "ashr")
  res_df <- as.data.frame(res_cv)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj), ]
  
  de_results_deseq[[cv]] <- res_df
  dds_objects[[cv]] <- dds_cv
  vsd_objects[[cv]] <- vst(dds_cv, blind = FALSE)
  
  write.table(res_df, file = paste0(outpathcount, "DEgenes_DESeq2_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  title_cv <- paste0(cultivar_labels[[cv]], ":control vs xylella")
  
  pdf(file= paste0(outpathcount,"dispersion_", cv, ".pdf"), width = 7, height =6)
  plotDispEsts(dds_cv, main = paste("Dispersion -", title_cv))
  dev.off()
  
  pdf(file = paste0(outpathcount, "Maplot_", cv, ".pdf"), width = 7, height =6)
  plotDispEsts(dds_cv, main = paste("Dispersion -", title_cv))
  dev.off()
}

sapply(de_results_deseq, function(df) sum(df$padj < sig_cutoff, na.rm = TRUE))


#HEATMAP AND COEXPRESSION NETWORK, PER CULTIVAR
network_objects <-list()

for (cv in cultivars) {
  res_df <- de_results_deseq[[cv]]
  sig_cv <- res_df[!is.na(res_df$padj) & res_df$padj < sig_cutoff, ]
  
  if(nrow(sig_cv) == 0) next
  
  title_cv <- paste0(cultivar_labels[[cv]], ": control vs xylella")
  
  vsd_mat_cv <- assay(vsd_objects[[cv]])
  mat_cv <-vsd_mat_cv[rownames(vsd_mat_cv) %in% sig_cv$gene, , drop = FALSE]
  mat_scaled_cv <- t(scale(t(mat_cv)))
  
  metadata_cv <- metadata[metadata$cultivar == cv, ]
  ann_col_cv <- data.frame(
    treatment = metadata_cv$treatment[match(colnames(mat_scaled_cv), metadata_cv$sample)]
  )
  rownames(ann_col_cv) <- colnames(mat_scaled_cv)
  
  while (!is.null(dev.list())) dev.off()
  
  pdf(file = paste0(outpathcount, "heatmap_", cv, ".pdf"), width = 8, height = 10)
  cat("Title for this interaction:", title_cv, "\n")
  pheatmap(mat_scaled_cv,
           annotation_col = ann_col_cv,
           show_rownames = TRUE,
           fontsize_row = 6,
           main = paste("Significant genes (padj < 0.05) -", title_cv))
  dev.off()
  
  cor_genes_cv <- cor(t(mat_cv), method = "pearson")
  threshold <- 0.8
  adj_cv <- abs(cor_genes_cv) >= threshold
  diag(adj_cv) <- FALSE
  
  g_cv <- igraph::simplify(graph_from_adjacency_matrix(adj_cv, mode = "undirected", diag = FALSE))
  V(g_cv)$direction <- ifelse(sig_cv$log2FoldChange[match(V(g_cv)$name, sig_cv$gene)] > 0, "Up", "Down")
  V(g_cv)$color <- ifelse(V(g_cv)$direction == "Up", "deeppink", "darkseagreen")
  
  network_objects[[cv]] <- g_cv
  
  while(!is.null(dev.list())) dev.off()
  
  pdf(file = paste0(outpathcount, "coexpression_network_", cv, ".pdf"), width = 9, height = 9)
  plot(g_cv,
       vertex.label = NA,
       vertex.size = 5,
       main = paste("Coexpression network -", title_cv))
  legend("bottomright", legend = c("Upregulated", "Downregulated"),
         pt.bg = c("deeppink", "darkseagreen"), pch = 21, pt.cex = 1.5,
         bty = "n")
  dev.off()
  
  cat(cv, "-- network connections:", ecount(g_cv), "\n")
  
  # Supplementary table: identity of each gene in the network
  network_gene_table <- data.frame(
    gene = V(g_cv)$name,
    direction = V(g_cv)$direction,
    log2FoldChange = sig_cv$log2FoldChange[match(V(g_cv)$name, sig_cv$gene)],
    degree = igraph::degree(g_cv)
  )
  network_gene_table <- network_gene_table[order(-network_gene_table$degree), ]
  
  write.table(network_gene_table,
              file = paste0(outpathcount, "coexpression_network_genes_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
dev.list()
while (!is.null(dev.list())) dev.off()


#COMPARATIVE BARPLOT: UP/DOWN REGULATED PER CULTIVAR
gene_counts <- do.call (rbind, lapply(cultivars, function (cv) {
  sig_cv <- de_results_deseq[[cv]][!is.na(de_results_deseq[[cv]]$padj) &
                                     de_results_deseq[[cv]]$padj < sig_cutoff, ]
  data.frame(
    Cultivar = cultivar_labels[[cv]],
    Upregulated = sum(sig_cv$log2FoldChange > 0),
    Downregulated = sum(sig_cv$log2FoldChange <0)
  )
}))

gene_counts

gene_counts_long <- reshape2::melt(gene_counts, id.vars = "Cultivar",
                                   variable.name = "Regulation",
                                   value.name = "Count")

p_updown <- ggplot(gene_counts_long, aes(x = Cultivar, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label=Count), position = position_dodge(width=0.9), vjust= -0.3) +
  scale_fill_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue")) +
  labs(title = "Differentially expressed genes (Upregulated vs Downregulated", 
       x = "Cultivar", y = "Number of genes") +
  theme_minimal()

print(p_updown)
ggsave(paste0(outpathcount, "genes_up_downregulated_per_cultivar.pdf"), plot = p_updown,
       width = 8, height =6)

#VENN DIAGRAM AND SHARED GENES 
if(!requireNamespace("ggVennDiagram", quietly = TRUE)) install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

#Significant genes per cultivar
sig_genes_list <- lapply(cultivars, function(cv) {
  res_df <- de_results_deseq[[cv]]
  res_df$gene[!is.na(res_df$padj) & res_df$padj < sig_cutoff]
})
names(sig_genes_list) <- cultivar_labels[cultivars]

#VENN DIAGRAM
p_venn <- ggVennDiagram(sig_genes_list, label = "count", set_size = 4) +
  scale_fill_gradient(low = "darkturquoise", high = "chocolate") +
  scale_x_continuous(expand = expansion (mult = 0.15)) +
  labs(title = "Differentially expressed genes (padj<0.05)") +
  theme(legend.position = "right")

print(p_venn)

ggsave(paste0(outpathcount, "venn_DEGs_cultivars.pdf"), plot = p_venn,
       width = 8, height = 7)

shared_genes <- Reduce(intersect, sig_genes_list)
cat("Differentially expressed genes shared between", paste(names(sig_genes_list), collapse = " and "),
":", length(shared_genes), "\n")

#Genes exclusive per cultivar
exclusive_genes <- lapply(names(sig_genes_list), function(nm) {
  setdiff(sig_genes_list[[nm]], unlist(sig_genes_list[names(sig_genes_list) != nm]))
})
names(exclusive_genes) <-names(sig_genes_list)

#Save exclusive genes per cultivar in .txt
for(nm in names(exclusive_genes)) {
  write.table(data.frame(gene = exclusive_genes[[nm]]),
              file = paste0(outpathcount, "DEgenes_exclusive_", nm, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

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

#HYPERGEOMETRIC ANALYSIS
for(cv in cultivars) {
  
  cat("=== Hypergeometric for", cv, "===\n")
  
  res_df <- de_results_deseq[[cv]]
  universe_genes <- res_df$gene[!is.na(res_df$padj)]
  
  #1. Using foldchange filter
  genes_with_fc <- res_df$gene[!is.na(res_df$padj) & res_df$padj < sig_cutoff &
                                 abs(res_df$log2FoldChange) > log2fc_cutoff]
  
  enrich_with_fc <- enricher(gene = genes_with_fc,
                             universe = universe_genes,
                             TERM2GENE = term2gene_df,
                             TERM2NAME = term2name_df,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 1, qvalueCutoff = 1)
  
  write.table(as.data.frame(enrich_with_fc),
              file = paste0(outpathcount, "enrich_hyper_withFC_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  #2. Without using filter foldchange
  genes_no_fc <- res_df$gene[!is.na(res_df$padj) & res_df$padj < sig_cutoff]
  
  enrich_no_fc <- enricher(gene = genes_no_fc,
                           universe = universe_genes,
                           TERM2GENE = term2gene_df,
                           TERM2NAME= term2name_df,
                           pAdjustMethod = "BH",
                           pvalueCutoff = 1, qvalueCutoff = 1)
  write.table(as.data.frame(enrich_no_fc),
              file = paste0(outpathcount, "enrich_hyper_noFC_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  hyper_results[[cv]] <- list(withFC = enrich_with_fc, noFC = enrich_no_fc)
  
  cat(cv, "-- With FoldChange:", length(genes_with_fc), "entry genes,",
      nrow(as.data.frame(enrich_with_fc)), "enriched terms\n")
  cat(cv, "-- Without FoldChange:", length(genes_no_fc), "entry genes,",
      nrow(as.data.frame(enrich_no_fc)), "enriched terms\n")
}

names(hyper_results)
names(hyper_results[["Catuai"]])

#CATUAI p.adjust filter
sig_with_fc_Catuai <- as.data.frame(hyper_results[["Catuai"]]$withFC)
sum(sig_with_fc_Catuai$p.adjust < 0.05, na.rm = TRUE)

sig_no_fc_Catuai <- as.data.frame(hyper_results[["Catuai"]]$noFC)
sum(sig_no_fc_Catuai$p.adjust <0.05, na.rm = TRUE)

#CR95 p.adjust filter
sig_with_fc_CR95 <- as.data.frame(hyper_results[["CR95"]]$withFC)
sum(sig_with_fc_CR95$p.adjust < 0.05, na.rm = TRUE)

sig_no_fc_CR95 <- as.data.frame(hyper_results[["CR95"]]$noFC)
sum(sig_no_fc_CR95$p.adjust <0.05, na.rm = TRUE)


#TOP CATEGORIES HYPERGEOMETRIC ANALYSIS
hyper_results <- list()
for (cv in cultivars) {
  cat("\n === Breakdown", cv, "===\n")
  
  enrich_df <- as.data.frame(hyper_results[[cv]]$noFC)
  enrich_df <- merge(enrich_df, go_names, by.x = "ID", by.y = "GOID", all.x = TRUE)
  
  cat("-- Counting by ontology (hypergeometric without FC): \n")
  print(table(enrich_df$ONTOLOGY))
  
  cat("-- Top 10 BP (hypergeometric without FC: \n")
  print(head(enrich_df[enrich_df$ONTOLOGY == "BP", ][order(enrich_df[enrich_df$ONTOLOGY == "BP", ]$p.adjust),
                                                     c("ID", "Description", "p.adjust", "Count")], 10))
  
  write.table(enrich_df, file = paste0(outpathcount, "enrich_hyper_noFC_withOntology_", cv, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

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

#SIGNIFICANT MODULES SEARCH

module_results <- list()
min_module_size <- 3

for(cv in cultivars) {
  
  cat("\n=== Modules for", cv, "===\n")
  
  g_cv <- network_objects[[cv]]
  
  if (is.null(g_cv) || ecount(g_cv) == 0) {
    cat(cv, "-- network has no connections, modules cannot be detected\n")
    next
  }
  
  comm_cv <- igraph::cluster_louvain(g_cv)
  
  V(g_cv)$module <- membership(comm_cv)
  
  module_sizes <- sizes(comm_cv)
  cat("Total number of modules detected:", length(module_sizes), "\n")
  print(module_sizes)
  
  sig_modules <- names(module_sizes)[module_sizes >= min_module_size]
  cat("Modules with >=", min_module_size, "genes:", length(sig_modules), "\n")
  
  cat("Network modularity:", modularity(comm_cv), "\n")
  
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

#HYPERGEOMETRIC ANALYSIS PLOT (WITH AND WITHOUT FC)

fc_types <- c("withFC", "noFC")
fc_labels <- c(withFC = "with FoldChange", noFC = "without FoldChange")

for (cv in cultivars) {
  for (fc_type in fc_types) {
    
    enrich_obj <- hyper_results[[cv]][[fc_type]]
    
    if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      cat(cv, "-", fc_type, "-- no results to plot\n")
      next
    }
    
    title_cv <- paste0(cultivar_labels[[cv]], ": enriched categories (", fc_labels[[fc_type]], ")")
    
    # DOTPLOT
    p_dot <- dotplot(enrich_obj, showCategory = 15) + ggtitle(title_cv)
    print(p_dot)
    ggsave(paste0(outpathcount, "dotplot_enrich_", cv, "_", fc_type, ".pdf"),
           plot = p_dot, width = 8, height = 8)
    
    # BARPLOT
    p_bar <- barplot(enrich_obj, showCategory = 15) + ggtitle(title_cv)
    print(p_bar)
    ggsave(paste0(outpathcount, "barplot_enrich_", cv, "_", fc_type, ".pdf"),
           plot = p_bar, width = 8, height = 8)
  }
}

#CAMERA PLOT
for (cv in cultivars) {
  
  camera_df <- camera_results2[[cv]]
  camera_bp <- camera_df[camera_df$ONTOLOGY == "BP", ]
  top_camera <- head(camera_bp[order(camera_bp$FDR), ], 15)
  top_camera$TERM <- factor(top_camera$TERM, levels = rev(top_camera$TERM))
  
  title_cv <- paste0(cultivar_labels[[cv]], ": CAMERA - top BP categories")
  
  p_camera <- ggplot(top_camera, aes(x = TERM, y = -log10(FDR))) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    labs(title = title_cv, x = NULL, y = "-log10(FDR)") +
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

# CR95 genes, significant, annotated to that GO term -- sorted by padj
genes_ethylene_CR95 <- de_results_deseq$CR95[
  de_results_deseq$CR95$gene %in% unique(unlist(term2gene[ethylene_go])) &
    !is.na(de_results_deseq$CR95$padj) & de_results_deseq$CR95$padj < sig_cutoff, ]

genes_ethylene_CR95 <- genes_ethylene_CR95[order(genes_ethylene_CR95$padj), ]
genes_ethylene_CR95

#MULTIFASTA CODE TO MAKE A BLAST ANALYSIS (CR95)
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
library(rentrez)

gene_ids <- genes_ethylene_CR95$gene

fasta_sequences <- character(0)

for (g in gene_ids) {
  search_res <- entrez_search(db = "gene", term = paste0(g, "[sym] AND Coffea arabica[orgn]"))
  if (length(search_res$ids) == 0) { cat("Not found:", g, "\n"); next }
  
  link_res <- entrez_link(dbfrom = "gene", id = search_res$ids[1], db = "nuccore")
  nuc_ids  <- link_res$links$gene_nuccore_refseqrna
  if (is.null(nuc_ids)) { cat("No linked RefSeq:", g, "\n"); next }
  
  seq_fasta <- entrez_fetch(db = "nuccore", id = nuc_ids[1], rettype = "fasta")
  fasta_sequences <- c(fasta_sequences, seq_fasta)
  
  Sys.sleep(0.4)
}

writeLines(fasta_sequences, paste0(outpathcount, "genes_ethylene_CR95.fasta"))
length(fasta_sequences)  # confirms how many sequences were actually retrieved (should be close to 10)

#MULTIFASTA CODE TO MAKE A BLAST ANALYSIS (UPREGULATED WITH FOLDCHANGE) (Catuai & CR95)
photo_go <- names(Term(GOTERM))[grepl("photosynthesis", Term(GOTERM), ignore.case = TRUE)]
photo_genes <- unique(unlist(term2gene[photo_go]))

for (cv in cultivars) {
  res_photo <- de_results_deseq[[cv]][de_results_deseq[[cv]]$gene %in% photo_genes, ]
  
  n_total <- nrow(res_photo)
  n_sig_up   <- sum(!is.na(res_photo$padj) & res_photo$padj < sig_cutoff & res_photo$log2FoldChange > 0)
  n_sig_down <- sum(!is.na(res_photo$padj) & res_photo$padj < sig_cutoff & res_photo$log2FoldChange < 0)
  n_no_sig   <- n_total - n_sig_up - n_sig_down
  
  cat(cv, "-- total photosynthesis genes:", n_total,
      "| significant UP:", n_sig_up,
      "| significant DOWN:", n_sig_down,
      "| no significant change:", n_no_sig, "\n")
}

thylakoid_go <- names(Term(GOTERM))[grepl("thylakoid|photosynthetic membrane|photosynthesis|chlorophyll|plastid translation",
                                          Term(GOTERM), ignore.case = TRUE)]
Term(GOTERM)[thylakoid_go]  # check which terms were captured

thylakoid_genes <- unique(unlist(term2gene[thylakoid_go]))

sig_thylakoid_Catuai <- de_results_deseq$Catuai[
  de_results_deseq$Catuai$gene %in% thylakoid_genes &
    !is.na(de_results_deseq$Catuai$padj) & de_results_deseq$Catuai$padj < sig_cutoff, ]

cat("Total significant genes in photosynthetic machinery categories:", nrow(sig_thylakoid_Catuai), "\n")
cat("Upregulated:", sum(sig_thylakoid_Catuai$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_thylakoid_Catuai$log2FoldChange < 0), "\n")

hist(sig_thylakoid_Catuai$log2FoldChange, breaks = 30,
     main = "log2FC distribution -- photosynthetic machinery genes (Catuai)")


###########
teammate_genes <- c("LOC113738600", "LOC113741493", "LOC113716386", "LOC113722862",
                    "LOC113716505", "LOC113712194", "LOC140003830", "LOC113729884",
                    "LOC113693854", "LOC113709139", "LOC113698731")

check_concordance <- function(cv) {
  res <- de_results_deseq[[cv]]
  res_sub <- res[res$gene %in% teammate_genes, c("gene", "log2FoldChange", "padj")]
  res_sub$significant_in_my_analysis <- !is.na(res_sub$padj) & res_sub$padj < sig_cutoff
  res_sub$direction <- ifelse(res_sub$log2FoldChange > 0, "Up", "Down")
  res_sub[order(res_sub$padj), ]
}

concordance_Catuai <- check_concordance("Catuai")
concordance_CR95   <- check_concordance("CR95")

cat("=== Catuai ===\n"); concordance_Catuai
cat("=== CR95 ===\n");   concordance_CR95

cat("Missing in my Catuai:", setdiff(teammate_genes, de_results_deseq$Catuai$gene), "\n")
cat("Missing in my CR95:",   setdiff(teammate_genes, de_results_deseq$CR95$gene), "\n")

cat("Match as significant in MY Catuai:", sum(concordance_Catuai$significant_in_my_analysis, na.rm = TRUE),
    "of", nrow(concordance_Catuai), "found\n")
cat("Match as significant in MY CR95:", sum(concordance_CR95$significant_in_my_analysis, na.rm = TRUE),
    "of", nrow(concordance_CR95), "found\n")

homeostasis_go <- names(Term(GOTERM))[grepl("homeostasis", Term(GOTERM), ignore.case = TRUE)]
homeostasis_genes <- unique(unlist(term2gene[homeostasis_go]))

data.frame(
  gene = teammate_genes,
  in_my_homeostasis_category = teammate_genes %in% homeostasis_genes
)

#MULTIFASTA FILE HOMEOSTASIS CATUAI
# ============================================================
# ============================================================
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
library(rentrez)

homeostasis_gene_ids <- c(
  "LOC113693854", "LOC113698731", "LOC113709139", "LOC113712194",
  "LOC113716386", "LOC113716505", "LOC113722862", "LOC113729884",
  "LOC113738600", "LOC113741493", "LOC140003830"
)

fetch_fasta <- function(gene_ids, label) {
  fasta_sequences <- character(0)
  
  for (g in gene_ids) {
    search_res <- entrez_search(db = "gene", term = paste0(g, "[sym] AND Coffea arabica[orgn]"))
    if (length(search_res$ids) == 0) { cat("Not found:", g, "\n"); next }
    
    link_res <- entrez_link(dbfrom = "gene", id = search_res$ids[1], db = "nuccore")
    nuc_ids  <- link_res$links$gene_nuccore_refseqrna
    if (is.null(nuc_ids)) { cat("No linked RefSeq:", g, "\n"); next }
    
    seq_fasta <- entrez_fetch(db = "nuccore", id = nuc_ids[1], rettype = "fasta")
    fasta_sequences <- c(fasta_sequences, seq_fasta)
    
    Sys.sleep(0.4)
  }
  
  writeLines(fasta_sequences, paste0(outpathcount, "genes_homeostasis_", label, ".fasta"))
  cat(label, "-- sequences retrieved:", length(fasta_sequences), "of", length(gene_ids), "\n")
}

fetch_fasta(homeostasis_gene_ids, "Catuai")

#MULTIFASTA FILE CR95
# ============================================================
# SIGNIFICANT HOMEOSTASIS GENES IN CR95
# (uses the same homeostasis GO annotation already loaded)
# ============================================================

# If these are no longer in your environment from the previous session:
if (!exists("homeostasis_go")) {
  homeostasis_go <- names(Term(GOTERM))[grepl("homeostasis", Term(GOTERM), ignore.case = TRUE)]
}
if (!exists("homeostasis_genes")) {
  homeostasis_genes <- unique(unlist(term2gene[homeostasis_go]))
}

genes_homeostasis_CR95 <- de_results_deseq$CR95[
  de_results_deseq$CR95$gene %in% homeostasis_genes &
    !is.na(de_results_deseq$CR95$padj) & de_results_deseq$CR95$padj < sig_cutoff, ]

genes_homeostasis_CR95 <- genes_homeostasis_CR95[order(genes_homeostasis_CR95$padj), ]

cat("CR95 -- significant homeostasis genes:", nrow(genes_homeostasis_CR95), "\n")
genes_homeostasis_CR95

#HEATMAP Analysis
for (cv in cultivars) {
  res_photo <- de_results_deseq[[cv]][de_results_deseq[[cv]]$gene %in% photo_genes, ]
  
  n_total <- nrow(res_photo)
  n_sig_up   <- sum(!is.na(res_photo$padj) & res_photo$padj < sig_cutoff & res_photo$log2FoldChange > 0)
  n_sig_down <- sum(!is.na(res_photo$padj) & res_photo$padj < sig_cutoff & res_photo$log2FoldChange < 0)
  n_no_sig   <- n_total - n_sig_up - n_sig_down
  
  cat(cv, "-- total photosynthesis genes:", n_total,
      "| significant UP:", n_sig_up,
      "| significant DOWN:", n_sig_down,
      "| no significant change:", n_no_sig, "\n")
}

cat("Upregulated:", sum(sig_thylakoid_Catuai$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_thylakoid_Catuai$log2FoldChange < 0), "\n")
