# Table 2. PERMANOVA of plaque community composition (Aitchison distance)
# for each protein exposure, adjusted for age, gender and toothbrushing.

source("R/00_setup.R")

res <- readRDS(file.path(PATHS$results, "community_composition.rds"))$permanova
res <- res[match(EXPOSURES, res$exposure), ]

table2 <- data.frame(
  Exposure = EXPO_LAB[res$exposure],
  n        = res$n,
  Df       = res$Df,
  `R2 (%)` = sprintf("%.2f", 100 * res$R2),
  F        = sprintf("%.2f", res$F),
  p        = fmt_p(res$`Pr(>F)`),
  q        = fmt_p(res$q),
  check.names = FALSE)

write_table(table2, "table2_permanova")
