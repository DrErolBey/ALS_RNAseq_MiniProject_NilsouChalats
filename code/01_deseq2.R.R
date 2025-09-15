# 01_deseq2.R  DESeq2 pipeline
# Author: Nilsou Chalats
# Data: GSE234297 (3 ALS vs 3 CTRL)
# Inputs: GSE234297_gene_raw_counts.txt
# Outputs: report/*.csv, report/*.pdf
# Deps: readr, dplyr, DESeq2, apeglm, ggplot2, AnnotationDbi, org.Hs.eg.db

root <- "ALS_RNAseq_MiniProject_NilsouChalats"

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(DESeq2); library(apeglm); library(ggplot2)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

# --- I/O ---
counts_file <- "GSE234297_gene_raw_counts.txt"
outdir <- "report"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Read counts ---
cts  <- read_tsv(counts_file)
gene <- cts[[1]]
mat  <- as.matrix(cts[ , -1]); rownames(mat) <- gene; storage.mode(mat) <- "integer"

stopifnot(is.integer(mat), !anyNA(mat), all(mat >= 0), length(unique(rownames(mat))) == nrow(mat))

# --- Select 3v3 ALS vs CTRL by column names ---
all_cols  <- colnames(mat)
als_cols  <- grep("Case|ALS",              all_cols, ignore.case=TRUE, value=TRUE)
ctrl_cols <- grep("Control|CTRL|Healthy",  all_cols, ignore.case=TRUE, value=TRUE)
n <- min(3, length(als_cols), length(ctrl_cols))
if (n < 2) stop("Insufficient samples. At least 2 ALS and 2 CTRL are required.")
keep_cols <- c(als_cols[seq_len(n)], ctrl_cols[seq_len(n)])
mat <- mat[, keep_cols]

# --- Metadata ---
condition <- ifelse(grepl("Case|ALS", keep_cols, ignore.case=TRUE), "ALS", "CTRL")
coldata <- data.frame(condition=factor(condition, levels=c("CTRL","ALS")), row.names=keep_cols)

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData=mat, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, , drop=FALSE]
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

# --- Results table ---
res_tbl <- as.data.frame(res)
res_tbl$gene <- rownames(res)
res_tbl <- arrange(res_tbl, padj)
write_csv(res_tbl, file.path(outdir, sprintf("deseq2_results_%sv%s.csv", n, n)))
write_csv(arrange(res_tbl, padj), file.path(outdir, "deseq2_results_full.csv"))

# --- QC: PCA & MA ---
vsd  <- vst(dds)
pdat <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
pv   <- round(100*attr(pdat, "percentVar"))
p <- ggplot(pdat, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", pv[1], "%")) +
  ylab(paste0("PC2: ", pv[2], "%")) +
  theme_bw()
ggsave(file.path(outdir, "pca_plot.pdf"), p, width=6, height=4)
pdf(file.path(outdir, "ma_plot.pdf")); plotMA(res); dev.off()

# --- Up/Down lists (correct sign + FDR) ---
ok <- res_tbl %>%
  select(gene, log2FoldChange, padj, baseMean) %>%
  filter(!is.na(padj), !is.na(log2FoldChange))
up <- ok %>% filter(padj < 0.05, log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>% slice_head(n=10)
down <- ok %>% filter(padj < 0.05, log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>% slice_head(n=10)
write_csv(up,   file.path(outdir, "top_up_genes.csv"))
write_csv(down, file.path(outdir, "top_down_genes.csv"))
stopifnot(all(up$log2FoldChange > 0), all(down$log2FoldChange < 0))

# --- Volcano consistent with up/down ---
vtbl <- res_tbl %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(dir = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))
v <- ggplot(vtbl, aes(x=log2FoldChange, y=-log10(padj), color=dir)) +
  geom_point(size=1.6, alpha=0.7) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()
ggsave(file.path(outdir, "volcano_plot.pdf"), v, width=6, height=4)

# --- SYMBOL mapping (ENTREZID -> SYMBOL) ---
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys=as.character(res_tbl$gene),
                             keytype="ENTREZID",
                             columns="SYMBOL")
res_sym <- res_tbl %>%
  mutate(gene = as.character(gene)) %>%
  left_join(map, by=c("gene"="ENTREZID")) %>%
  relocate(SYMBOL, .before=gene)
write_csv(arrange(res_sym, padj), file.path(outdir, "deseq2_results_full_symbol.csv"))

up_sym <- read_csv(file.path(outdir, "top_up_genes.csv"),
                   col_types=cols(gene=col_character(), .default=col_double())) %>%
  left_join(map, by=c("gene"="ENTREZID")) %>% relocate(SYMBOL, .before=gene)
down_sym <- read_csv(file.path(outdir, "top_down_genes.csv"),
                     col_types=cols(gene=col_character(), .default=col_double())) %>%
  left_join(map, by=c("gene"="ENTREZID")) %>% relocate(SYMBOL, .before=gene)
write_csv(up_sym,   file.path(outdir, "top_up_genes_symbol.csv"))
write_csv(down_sym, file.path(outdir, "top_down_genes_symbol.csv"))

# --- Optional: GO enrichment (assumes packages already installed) ---
if (requireNamespace("clusterProfiler", quietly=TRUE)) {
  library(clusterProfiler)
  sig <- subset(res_tbl, !is.na(padj) & padj < 0.05)$gene
  eg  <- bitr(sig, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)
  ego <- enrichGO(eg$SYMBOL, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP")
  write_csv(as.data.frame(ego), file.path(outdir, "go_bp.csv"))
}

# --- Copy to portfolio (run after volcano/go exist) ---
root <- "ALS_RNAseq_MiniProject_NilsouChalats"
dir.create(file.path(root,"figs"),    showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(root,"results"), showWarnings=FALSE, recursive=TRUE)

file.copy(file.path(outdir, "pca_plot.pdf"),
          file.path(root, "figs", "PCA_3v3_GSE234297.pdf"), overwrite=TRUE)
file.copy(file.path(outdir, "ma_plot.pdf"),
          file.path(root, "figs", "MAplot_3v3_GSE234297.pdf"), overwrite=TRUE)
file.copy(file.path(outdir, "volcano_plot.pdf"),
          file.path(root, "figs", "Volcano_3v3_GSE234297.pdf"), overwrite=TRUE)

file.copy(file.path(outdir, "deseq2_results_full.csv"),
          file.path(root, "results", "DESeq2_results_3v3_GSE234297.csv"), overwrite=TRUE)
if (file.exists(file.path(outdir,"go_bp.csv"))) {
  file.copy(file.path(outdir, "go_bp.csv"),
            file.path(root,  "results", "go_bp.csv"), overwrite=TRUE)
}
