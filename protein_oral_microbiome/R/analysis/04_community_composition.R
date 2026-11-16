# Community composition on Aitchison distance:
#   PERMANOVA for each exposure (ASV level, filtered), sequential terms
#   PCA of CLR abundances (filtered ASVs) and PC1 vs animal protein density
#   Constrained partial ordination (all ASVs), conditioned on age and gender
#   Spearman correlations of CLR genus abundances with each exposure

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(phyloseq); library(microViz); library(vegan); library(dplyr)
})

obj     <- readRDS(file.path(PATHS$results, "phyloseq.rds"))
ps      <- obj$ps
ps_filt <- obj$ps_filt

## PERMANOVA ------------------------------------------------------------
X <- as(otu_table(ps_filt), "matrix")
if (taxa_are_rows(ps_filt)) X <- t(X)
clr  <- clr_transform(X)
Dmat <- as.matrix(dist(clr, method = "euclidean"))
meta <- data.frame(sample_data(ps_filt), check.names = FALSE)

run_permanova <- function(expo) {
  m <- meta[complete.cases(meta[, c(PERMANOVA_COVARS, expo)]), , drop = FALSE]
  D <- as.dist(Dmat[rownames(m), rownames(m)])
  f <- as.formula(paste("D ~", paste(c(PERMANOVA_COVARS, expo), collapse = " + ")))
  set.seed(1234)
  ad <- as.data.frame(adonis2(f, data = m, permutations = 999, by = "terms"))
  ad$term <- rownames(ad)
  ad <- ad[ad$term %in% c(PERMANOVA_COVARS, expo),
           c("term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
  ad$exposure <- expo
  ad$n <- nrow(m)
  ad
}
permanova_all <- do.call(rbind, lapply(EXPOSURES, run_permanova))
permanova <- permanova_all[permanova_all$term == permanova_all$exposure, ]
permanova$q <- p.adjust(permanova$`Pr(>F)`, method = "BH")
print(permanova, row.names = FALSE)

## PCA of CLR abundances -----------------------------------------------
pca <- prcomp(clr, center = TRUE, scale. = FALSE)
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_scores <- data.frame(sample_id = rownames(pca$x), pca$x[, 1:4]) |>
  left_join(meta[, c("sample_id", "protein_group", PRIMARY)], by = "sample_id")
pc1_rho <- cor.test(pca_scores$PC1, pca_scores$animal_pct_energy,
                    method = "spearman", exact = FALSE)

## Constrained partial ordination ---------------------------------------
# NB: the "Post_glucose_pH" constraint is the salivary pH FALL after the
# glucose challenge (saliva_ph_fall), as in the published figure.
CONSTRAINTS <- c("Animal_protein", "Plant_protein", "Filled_teeth",
                 "Missing_teeth", "Post_glucose_pH", "Physical_activity")
ord <- ps |>
  tax_transform("identity") |>
  dist_calc("aitchison") |>
  ps_mutate(Animal_protein    = animal_pct_energy,
            Plant_protein     = plant_pct_energy,
            Filled_teeth      = as.numeric(filled_teeth),
            Missing_teeth     = as.numeric(missing_teeth),
            Post_glucose_pH   = saliva_ph_fall,
            Physical_activity = met_min_week,
            age_s             = as.numeric(scale(age)),
            sex_s             = as.numeric(gender == "Female")) |>
  ord_calc(constraints = CONSTRAINTS, conditions = c("age_s", "sex_s"))

## Genus-level Spearman correlations -----------------------------------
cor_test_spearman <- function(formula, data, ...) {
  mf <- model.frame(formula, data)
  cor.test(mf[[1]], mf[[2]], method = "spearman", exact = FALSE)
}
top_taxa <- tax_top(ps, 20, by = max, rank = "Genus")
genus_cor <- ps |>
  tax_model(trans = "clr", rank = "Genus", taxa = top_taxa,
            variables = as.list(EXPOSURES), type = cor_test_spearman,
            return_psx = FALSE, verbose = FALSE) |>
  tax_models2stats(rank = "Genus") |>
  group_by(term) |>
  mutate(q = p.adjust(p.value, method = "BH")) |>
  ungroup()

saveRDS(list(permanova = permanova, permanova_all = permanova_all,
             pca_scores = pca_scores, pca_var = pca_var, pc1_rho = pc1_rho,
             ord = ord, genus_cor = genus_cor,
             n_asv = ntaxa(ps_filt), n_samples = nsamples(ps_filt)),
        file.path(PATHS$results, "community_composition.rds"))

