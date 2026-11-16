# Supplementary Figure 2. Alpha diversity by protein source group.
# Brackets show pairwise Wilcoxon tests with BH-adjusted q < 0.05.

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(rstatix)
  library(ggpubr)
})

alpha <- readRDS(file.path(PATHS$results, "alpha_diversity.rds"))

long <- alpha$values |>
  select(sample_id, protein_group, all_of(METRICS)) |>
  pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, METRICS, M_LAB[METRICS]))

stat <- long |>
  group_by(metric) |>
  wilcox_test(value ~ protein_group, p.adjust.method = "BH") |>
  add_xy_position(scales = "free_y", step.increase = 0.12) |>
  filter(p.adj < 0.05)

fig_s2 <- ggplot(long, aes(protein_group, value)) +
  geom_boxplot(aes(fill = protein_group, colour = protein_group), width = 0.55,
               alpha = 0.35, outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(aes(colour = protein_group), width = 0.16, size = 1.4, alpha = 0.75) +
  stat_pvalue_manual(stat, label = "p.adj.signif", tip.length = 0.01,
                     bracket.size = 0.35, hide.ns = TRUE) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = PAL) +
  scale_colour_manual(values = PAL) +
  scale_x_discrete(labels = c("Reference", "High\nanimal", "High\nplant")) +
  labs(title = "Alpha diversity by protein source group",
       subtitle = sprintf("Rarefied to %s reads (n = %d)",
                          format(RAREFY_DEPTH, big.mark = ","), nrow(alpha$values)),
       x = NULL, y = "Alpha diversity",
       caption = "Boxes show median and IQR. Pairwise Wilcoxon, BH-adjusted (* q < 0.05, ** q < 0.01).") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey40", hjust = 0, size = 8))

ggsave(file.path(PATHS$figures, "supp_figure2_alpha_groups.png"), fig_s2,
       width = 10, height = 7, dpi = 300, bg = "white")
