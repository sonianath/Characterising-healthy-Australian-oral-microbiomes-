# Imports the QIIME 2 artefacts, attaches participant data and creates
# the three phyloseq objects used downstream:
#   ps       all ASVs, raw counts (MaAsLin 3, genus correlations, networks)
#   ps_filt  ASVs in >=10% of samples with >=100 reads (PERMANOVA, PCA)
#   ps_rar   rarefied to 8,000 reads (alpha diversity)

setwd("/Users/a1799090/healthy_australian/protein-source-oral-microbiome")
list.files()    # should show R, data, run_all.R, README.md
source("R/00_setup.R")
suppressPackageStartupMessages({
  library(qiime2R); library(phyloseq); library(microViz)
})

meta <- read_participants()

ps <- qza_to_phyloseq(features = file.path(PATHS$qiime2, "table.qza"),
                      tree     = file.path(PATHS$qiime2, "rooted-tree.qza"),
                      taxonomy = file.path(PATHS$qiime2, "taxonomy.qza"))
ps <- prune_samples(sample_names(ps) %in% rownames(meta), ps)
sample_data(ps) <- sample_data(meta[sample_names(ps), ])
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
phyloseq_validate(ps)

ps <- ps |>
  tax_fix(min_length = 4, unknowns = c(""), sep = " ",
          anon_unique = TRUE, suffix_rank = "classified") |>
  tax_fix(unknowns = c("uncultured"))

otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
keep <- rowMeans(otu > 0) >= MIN_PREVALENCE & rowSums(otu) >= MIN_TOTAL_READS
ps_filt <- prune_taxa(keep, ps)
ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)

ps_rar <- rarefy_even_depth(ps, sample.size = RAREFY_DEPTH, rngseed = 42,
                            replace = FALSE, verbose = FALSE)

cat("Samples:", nsamples(ps), " ASVs:", ntaxa(ps),
    " filtered ASVs:", ntaxa(ps_filt),
    " rarefied samples:", nsamples(ps_rar), "\n")

saveRDS(list(ps = ps, ps_filt = ps_filt, ps_rar = ps_rar),
        file.path(PATHS$results, "phyloseq.rds"))
