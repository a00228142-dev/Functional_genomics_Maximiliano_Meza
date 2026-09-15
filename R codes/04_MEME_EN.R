# ============================================================
# Promoter extraction + motif analysis with MEME
# ============================================================

# --- 1. Required libraries ---
# BiocManager::install(c("rtracklayer", "Biostrings", "memes", "universalmotif"))

library(rtracklayer)
library(Biostrings)
library(memes)
library(universalmotif)

# --- 2. File paths (adjust if necessary) ---
genome_fasta   <- "GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna"
annotation_gff <- "genomic.gff"

file.exists(genome_fasta)     # should be TRUE
file.exists(annotation_gff)   # should be TRUE

# --- 3. Genes of interest ---
genes_interes <- c(
  "LOC113695805", "LOC140009270", "LOC113690379", "LOC113690378",
  "LOC113699104", "LOC113691088", "LOC113691272", "LOC113688572"
)

# --- 4. Import the GFF and keep only rows of type "gene" ---
gff <- import(annotation_gff, format = "gff3")
genes_gff <- gff[gff$type == "gene"]

# Check which column carries the gene ID (run this ONCE to confirm)
head(mcols(genes_gff))

# --- 5. Filter your genes of interest ---
# Adjust this line based on what you see in the head() above:

# Option A: the "gene" column has the clean ID (e.g. "LOC113695805")
genes_gr_filtrados <- genes_gff[genes_gff$gene %in% genes_interes]

# Option B (use only if A finds nothing): ID comes as "gene-LOC113695805"
# genes_gr_filtrados <- genes_gff[sub("gene-", "", genes_gff$ID) %in% genes_interes]

cat("Genes found:", length(genes_gr_filtrados), "of", length(genes_interes), "\n")
print(genes_gr_filtrados$gene)   # change to $ID if you used Option B

# --- 6. Define the promoter region (automatically respects +/- strand) ---
upstream_bp   <- 1000
downstream_bp <- 200

promotores_gr <- promoters(genes_gr_filtrados,
                           upstream = upstream_bp,
                           downstream = downstream_bp)
library(Rsamtools)

# --- 7. Load the genome as an indexed FaFile ---
# This creates a .fai index the first time (may take a while with large genomes)
if (!file.exists(paste0(genome_fasta, ".fai"))) {
  indexFa(genome_fasta)
}

fa <- FaFile(genome_fasta)

# --- Make sure seqnames match between the GFF and the genome ---
seqlevels(promotores_gr)          # chromosome/scaffold names in your GRanges
seqnames(scanFaIndex(fa))         # chromosome/scaffold names in the indexed FASTA

# If they don't match exactly, seqlevels(promotores_gr) must be adjusted before continuing

# --- Adjust seqinfo with the actual lengths from the indexed genome ---
fa_seqinfo <- seqinfo(fa)
seqlevels(promotores_gr) <- seqlevels(promotores_gr)[seqlevels(promotores_gr) %in% seqnames(fa_seqinfo)]
seqlengths(promotores_gr) <- seqlengths(fa_seqinfo)[seqlevels(promotores_gr)]
promotores_gr <- trim(promotores_gr)

# --- Extract the promoter sequences ---
promotores_seq <- getSeq(fa, promotores_gr)
names(promotores_seq) <- genes_gr_filtrados$gene   # or $ID, depending on which one you used

promotores_seq
writeXStringSet(promotores_seq, filepath = "promotores_8genes.fasta")

options(meme_bin = "/opt/miniconda3/envs/meme_env/bin/")
check_meme_install()
# --- 8. Run MEME de novo ---
check_meme_install()

resultados_meme <- runMeme(
  input   = promotores_seq,
  db      = NULL,
  nmotifs = 3,
  minw    = 6,
  maxw    = 15,
  mod     = "zoops",
  outdir  = "meme_output_8genes",
  parse_genomic_coord = FALSE
)

print(resultados_meme)

library(universalmotif)

# Convert from data.frame to a list of universalmotif objects
motivos_um <- to_list(resultados_meme)

# Now actually visualize
view_motifs(motivos_um)

jaspar_db <- "JASPAR_plants.meme"
file.exists(jaspar_db)   # should be TRUE

resultados_tomtom <- runTomTom(
  input    = resultados_meme,
  database = jaspar_db
)

print(resultados_tomtom)

library(dplyr)
install.packages("dplyr")
library(dplyr)
# View the top matches with their similarity score
resultados_tomtom %>%
  select(name, best_match_name, best_match_altname, best_match_pval, best_match_qval)
