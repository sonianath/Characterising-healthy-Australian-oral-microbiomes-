# Genus-level differential abundance with MaAsLin 3: one model per
# exposure, adjusted for age and gender. Animal and plant protein density
# are negatively correlated, so they are not entered jointly.

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(phyloseq); library(maaslin3)
})

ps <- readRDS(file.path(PATHS$results, "phyloseq.rds"))$ps

genus_ps  <- tax_glom(ps, taxrank = "Genus")
genus_mat <- as(otu_table(genus_ps), "matrix")
if (taxa_are_rows(genus_ps)) genus_mat <- t(genus_mat)      # samples x genera
colnames(genus_mat) <- make.unique(as.character(tax_table(genus_ps)[, "Genus"]))
genus_df <- as.data.frame(genus_mat)

meta <- data.frame(sample_data(ps), check.names = FALSE)[rownames(genus_df), ]

for (expo in EXPOSURES) {
  out_dir <- file.path(PATHS$results, "maaslin3", expo)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  maaslin3(
    input_data     = genus_df,
    input_metadata = meta,
    output         = out_dir,
    formula        = as.formula(paste("~", expo, "+ age + gender")),
    normalization  = "TSS",
    transform      = "LOG",
    augment        = TRUE,
    standardize    = TRUE,
    median_comparison_abundance  = TRUE,
    median_comparison_prevalence = FALSE,
    max_significance  = 0.1,
    min_prevalence    = 0.30,
    min_abundance     = 0,
    cores             = 1,
    plot_summary_plot = FALSE,
    plot_associations = FALSE
  )
}

da <- do.call(rbind, lapply(EXPOSURES, function(expo) {
  d <- read.delim(file.path(PATHS$results, "maaslin3", expo, "all_results.tsv"),
                  stringsAsFactors = FALSE)
  d <- d[d$metadata == expo & (is.na(d$error) | d$error == ""), ]
  d$exposure <- expo
  d[, c("feature", "exposure", "model", "coef", "stderr",
        "pval_individual", "qval_individual", "pval_joint", "qval_joint")]
}))

write.csv(da, file.path(PATHS$results, "maaslin3_all_results.csv"),
          row.names = FALSE)
