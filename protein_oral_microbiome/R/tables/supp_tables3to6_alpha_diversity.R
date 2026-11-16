# Supplementary Tables 3-6. Alpha diversity (rarefied to 8,000 reads).
#   S3  Spearman correlations with each exposure
#   S4  Median (IQR) by protein source group, Kruskal-Wallis
#   S5  Pairwise Wilcoxon tests for indices with Kruskal-Wallis p < 0.05
#   S6  Standardised linear regression adjusted for age and gender

source("R/00_setup.R")
suppressPackageStartupMessages(library(tidyr))

a <- readRDS(file.path(PATHS$results, "alpha_diversity.rds"))

to_wide <- function(d, value_cols) {
  d$metric <- factor(d$metric, METRICS, M_LAB[METRICS])
  d$exposure <- factor(d$exposure, EXPOSURES, EXPO_LAB[EXPOSURES])
  d <- d[order(d$exposure, d$metric), c("metric", "exposure", value_cols)]
  pivot_wider(d, id_cols = metric, names_from = exposure,
              values_from = all_of(value_cols), names_vary = "slowest",
              names_glue = "{exposure} {.value}") |>
    dplyr::rename(`Diversity index` = metric)
}

## S3 ---------------------------------------------------------------------
s3 <- a$spearman
s3$rho <- fmt_num(s3$rho)
s3$p   <- fmt_p(s3$p)
s3$q   <- fmt_p(s3$q)
write_table(to_wide(s3, c("rho", "p", "q")), "supp_table3_alpha_spearman")

## S4 ---------------------------------------------------------------------
g <- a$group_summary
g$value <- sprintf("%.2f (%.2f)", g$median, g$iqr)
s4 <- pivot_wider(g, id_cols = c(metric, kw_p), names_from = protein_group,
                  values_from = value)
s4 <- s4[match(METRICS, s4$metric), c("metric", GROUP_LEVELS, "kw_p")]
s4$metric <- M_LAB[s4$metric]
s4$kw_p   <- fmt_p(s4$kw_p)
names(s4) <- c("Diversity index", GROUP_LEVELS, "Kruskal-Wallis p")
write_table(s4, "supp_table4_alpha_groups")

## S5 ---------------------------------------------------------------------
sig <- unique(g$metric[g$kw_p < 0.05])
pw  <- a$pairwise[a$pairwise$metric %in% sig, ]
pw  <- pw[order(match(pw$metric, METRICS)), ]
s5 <- data.frame(
  `Diversity index` = M_LAB[pw$metric],
  Comparison = paste(pw$group1, "vs", tolower(pw$group2)),
  n = paste(pw$n1, "/", pw$n2),
  W = pw$statistic,
  p = fmt_p(pw$p),
  q = fmt_p(pw$p.adj),
  check.names = FALSE)
write_table(s5, "supp_table5_alpha_pairwise")

## S6 ---------------------------------------------------------------------
s6 <- a$regression
s6$`beta (SE)` <- sprintf("%s (%s)", fmt_num(s6$beta), fmt_num(s6$se))
s6$p <- fmt_p(s6$p)
s6$q <- fmt_p(s6$q)
write_table(to_wide(s6, c("beta (SE)", "p", "q")), "supp_table6_alpha_regression")
