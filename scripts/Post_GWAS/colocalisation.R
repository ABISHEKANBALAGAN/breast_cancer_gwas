#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(coloc)
})

# ----------------------------
# Config (edit as needed)
# ----------------------------
gwas_file  <- "results/harmonized_GC_fixed2.tsv"
eqtl_file  <- "GTEX/Breast_Mammary_Tissue.v8.eqtl_signifpairs.txt"

# GWAS design (set appropriately)
GWAS_IS_CASE_CONTROL <- TRUE
N_GWAS        <- 12898 + 388549
CASE_FRAC     <- 12898 / (12898 + 388549)

# eQTL design
N_EQTL <- 369

# Coloc filters
MIN_SNPS_PER_GENE <- 5

# Output
out_dir <- "results/coloc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Load GWAS
# ----------------------------
gwas <- fread(gwas_file)
cat("GWAS columns:\n"); print(colnames(gwas))

req_gwas <- c("variant_id","BETA","SE","MAF","P_GC")
missing_gwas <- setdiff(req_gwas, colnames(gwas))
if (length(missing_gwas)) {
  stop("GWAS file is missing required columns: ", paste(missing_gwas, collapse=", "))
}

gwas_std <- gwas %>%
  transmute(
    snp   = variant_id,
    beta  = as.numeric(BETA),
    se    = as.numeric(SE),
    maf_gwas = suppressWarnings(as.numeric(MAF)),
    p_gc = as.numeric(P_GC)
  ) %>%
  filter(!is.na(beta), !is.na(se), se > 0) %>%
  filter(!is.na(maf_gwas), maf_gwas > 0 & maf_gwas < 1) 

cat("GWAS rows after QC:", nrow(gwas_std), "\n")

# ----------------------------
# 2) Load eQTL
# ----------------------------
eqtl_cols <- c("variant_id","gene_id","maf","slope","slope_se")
eqtl <- fread(eqtl_file, select = eqtl_cols)
cat("eQTL columns:\n"); print(colnames(eqtl))

eqtl_std <- eqtl %>%
  transmute(
    snp        = variant_id,
    gene_id    = gene_id,
    beta_eqtl  = as.numeric(slope),
    se_eqtl    = as.numeric(slope_se),
    maf_eqtl   = suppressWarnings(as.numeric(maf))
  ) %>%
  filter(!is.na(beta_eqtl), !is.na(se_eqtl), se_eqtl > 0) %>%
  filter(!is.na(maf_eqtl), maf_eqtl > 0 & maf_eqtl < 1)

cat("eQTL rows after QC:", nrow(eqtl_std), "\n")

# ----------------------------
# 3) Merge by SNP
# ----------------------------
merged <- inner_join(gwas_std, eqtl_std, by = "snp", relationship = "many-to-many")
cat("Overlap SNP rows (all genes):", nrow(merged), "\n")

ngenes_any <- merged %>% summarise(n = n_distinct(gene_id)) %>% pull(n)
cat("Genes with any overlap:", ngenes_any, "\n")

# ----------------------------
# 4) Run coloc per gene
# ----------------------------
run_coloc_for_gene <- function(df_gene) {
  if (nrow(df_gene) < MIN_SNPS_PER_GENE) return(NULL)

  if (GWAS_IS_CASE_CONTROL) {
    ds1 <- list(
      beta   = df_gene$beta,
      varbeta= df_gene$se^2,
      snp    = df_gene$snp,
      MAF    = df_gene$maf_gwas,
      N      = N_GWAS,
      type   = "cc",
      s      = CASE_FRAC
    )
  } else {
    ds1 <- list(
      beta   = df_gene$beta,
      varbeta= df_gene$se^2,
      snp    = df_gene$snp,
      MAF    = df_gene$maf_gwas,
      N      = N_GWAS,
      type   = "quant",
      sdY    = 1
    )
  }

  ds2 <- list(
    beta   = df_gene$beta_eqtl,
    varbeta= df_gene$se_eqtl^2,
    snp    = df_gene$snp,
    MAF    = df_gene$maf_eqtl,
    N      = N_EQTL,
    type   = "quant",
    sdY    = 1
  )

  out <- try(coloc.abf(ds1, ds2), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

genes <- unique(merged$gene_id)
cat("Attempting coloc on", length(genes), "genes with MIN_SNPS_PER_GENE =", MIN_SNPS_PER_GENE, "\n")

summ_list <- list()
i <- 0
for (g in genes) {
  df_g <- merged %>% filter(gene_id == g)
  res  <- run_coloc_for_gene(df_g)
  if (is.null(res)) next

  # Transpose summary to get a 1-row data.frame
  summ <- as.data.frame(t(res$summary))
  summ$gene_id <- g
  summ$nsnps   <- nrow(df_g)
  summ_list[[length(summ_list) + 1]] <- summ

  # SNP-level results
  snp_res <- res$results
  snp_res$gene_id <- g
  out_prefix <- file.path(out_dir, paste0("coloc_", g))
  fwrite(snp_res, paste0(out_prefix, "_snp.tsv"), sep = "\t")

  i <- i + 1
  if (i %% 20 == 0) cat("Processed", i, "genes ...\n")
}

if (length(summ_list) == 0) {
  cat("No genes passed filters / produced results. Try lowering MIN_SNPS_PER_GENE or checking MAF.\n")
  quit(status = 0)
}

summary_all <- bind_rows(summ_list) %>% relocate(gene_id, nsnps)

# Safe sort/filter
pp_h4_col <- grep("^PP\\.H4", colnames(summary_all), value = TRUE)
if (length(pp_h4_col) == 1) {
  summary_all <- summary_all %>% arrange(desc(!!sym(pp_h4_col)))
  top_genes <- summary_all %>% filter(!!sym(pp_h4_col) > 0.7)
} else {
  warning("PP.H4 column not found; skipping top-gene filter")
  top_genes <- data.frame()
}

fwrite(summary_all, file.path(out_dir, "coloc_summary.tsv"), sep = "\t")
cat("Saved summary to", file.path(out_dir, "coloc_summary.tsv"), "\n")

cat("Top genes by PP.H4:\n")
print(head(summary_all, 10))
if (nrow(top_genes) > 0) {
  cat("Genes with PP.H4 > 0.7:\n")
  print(top_genes)
}
