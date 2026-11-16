# Supplementary Table 1. Genus-level features associated with each protein
# exposure in MaAsLin 3 (nominal p < 0.05), adjusted for age and gender.

source("R/00_setup.R")

da <- read.csv(file.path(PATHS$results, "maaslin3_all_results.csv"))
da <- da[da$pval_individual < 0.05, ]
da <- da[order(match(da$exposure, EXPOSURES), da$pval_individual), ]

supp1 <- data.frame(
  Exposure    = EXPO_LAB[da$exposure],
  Feature     = da$feature,
  Model       = tools::toTitleCase(da$model),
  `beta (SE)` = sprintf("%s (%s)", fmt_num(da$coef), fmt_num(da$stderr)),
  p           = fmt_p(da$pval_individual),
  q           = fmt_p(da$qval_individual),
  `q (joint)` = fmt_p(da$qval_joint),
  check.names = FALSE)

write_table(supp1, "supp_table1_maaslin3")
