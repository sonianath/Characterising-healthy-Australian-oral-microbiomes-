# Alpha diversity on rarefied data: Spearman correlations with each
# exposure, Kruskal-Wallis and pairwise Wilcoxon tests across protein
# source groups, and standardised linear models adjusted for age and gender.
# BH correction is applied within each exposure (or index) family.

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(phyloseq); library(picante); library(dplyr); library(tidyr)
  library(rstatix)
})

ps_rar <- readRDS(file.path(PATHS$results, "phyloseq.rds"))$ps_rar

alpha <- estimate_richness(ps_rar, measures = c("Observed", "Chao1",
                                                "Shannon", "Simpson"))
alpha$sample_id <- sample_names(ps_rar)
otu <- as(otu_table(ps_rar), "matrix")
if (taxa_are_rows(ps_rar)) otu <- t(otu)
alpha$FaithPD <- pd(otu, phy_tree(ps_rar), include.root = TRUE)$PD
alpha <- merge(alpha[, c("sample_id", METRICS)],
               data.frame(sample_data(ps_rar)), by = "sample_id")

## Spearman correlations -----------------------------------------------
spearman <- expand.grid(exposure = EXPOSURES, metric = METRICS,
                        stringsAsFactors = FALSE)
for (i in seq_len(nrow(spearman))) {
  ct <- cor.test(alpha[[spearman$metric[i]]], alpha[[spearman$exposure[i]]],
                 method = "spearman", exact = FALSE)
  spearman$rho[i] <- unname(ct$estimate)
  spearman$p[i]   <- ct$p.value
}
spearman <- spearman |>
  group_by(exposure) |>
  mutate(q = p.adjust(p, method = "BH")) |>
  ungroup()

## Protein source groups -----------------------------------------------
long <- alpha |>
  select(sample_id, protein_group, all_of(METRICS)) |>
  pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value")

group_summary <- long |>
  group_by(metric) |>
  group_modify(\(d, k) {
    s <- d |> group_by(protein_group) |>
      summarise(median = median(value), iqr = IQR(value), n = n(),
                .groups = "drop")
    s$kw_p <- kruskal.test(value ~ protein_group, data = d)$p.value
    s
  }) |>
  ungroup()

pairwise <- long |>
  group_by(metric) |>
  wilcox_test(value ~ protein_group, p.adjust.method = "BH") |>
  ungroup() |>
  as.data.frame()

## Adjusted linear models (standardised) --------------------------------
regression <- do.call(rbind, lapply(EXPOSURES, function(e) {
  do.call(rbind, lapply(METRICS, function(m) {
    d <- alpha
    d$.y <- as.numeric(scale(d[[m]]))
    d$.x <- as.numeric(scale(d[[e]]))
    s <- summary(lm(.y ~ .x + age + gender, data = d))$coefficients[".x", ]
    data.frame(exposure = e, metric = m, beta = s[["Estimate"]],
               se = s[["Std. Error"]], p = s[["Pr(>|t|)"]])
  }))
})) |>
  group_by(exposure) |>
  mutate(q = p.adjust(p, method = "BH")) |>
  ungroup()

cat("Alpha diversity: n =", nrow(alpha), "\n")
print(table(alpha$protein_group))

saveRDS(list(values = alpha, spearman = spearman,
             group_summary = group_summary, pairwise = pairwise,
             regression = regression),
        file.path(PATHS$results, "alpha_diversity.rds"))
