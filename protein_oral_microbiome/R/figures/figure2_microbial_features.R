# Figure 2. Microbial features associated with dietary protein source.
#   A  MaAsLin 3 abundance model, animal protein density (volcano)
#   B  Eigenvector centrality in the high animal vs high plant networks
#   C  Alpha diversity against animal and plant protein density

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(cowplot)
})

## A ---------------------------------------------------------------------
volcano <- read.csv(file.path(PATHS$results, "maaslin3_all_results.csv")) |>
  filter(exposure == "animal_pct_energy", model == "abundance") |>
  mutate(status = case_when(
    qval_individual < 0.1  & coef > 0 ~ "Positive (q < 0.1)",
    qval_individual < 0.1  & coef < 0 ~ "Negative (q < 0.1)",
    pval_individual < 0.05 & coef > 0 ~ "Positive (p < 0.05)",
    pval_individual < 0.05 & coef < 0 ~ "Negative (p < 0.05)",
    TRUE                              ~ "Not significant"),
    status = factor(status, levels = c("Positive (q < 0.1)", "Negative (q < 0.1)",
                                       "Positive (p < 0.05)", "Negative (p < 0.05)",
                                       "Not significant")))

enrich <- volcano |>
  summarise(n = n(), sig = sum(pval_individual < 0.05),
            fdr = sum(qval_individual < 0.1), expected = 0.05 * n(),
            binom_p = binom.test(sig, n, 0.05, alternative = "greater")$p.value) |>
  mutate(lab = sprintf("%d/%d at p < 0.05 (%.1f expected, p = %.3f)\n%d at q < 0.1",
                       sig, n, expected, binom_p, fdr))
xmax <- max(abs(volcano$coef), na.rm = TRUE) * 1.18

pA <- ggplot(volcano, aes(coef, -log10(pval_individual))) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, linewidth = 0.45,
             colour = "grey55") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey65") +
  geom_point(aes(colour = status, size = status, alpha = status)) +
  geom_text_repel(data = \(d) filter(d, pval_individual < 0.05),
                  aes(label = feature), size = 3.1, max.overlaps = 25,
                  min.segment.length = 0.2, box.padding = 0.45,
                  segment.colour = "grey60", segment.size = 0.3) +
  annotate("text", x = Inf, y = Inf, label = enrich$lab, hjust = 1.04,
           vjust = 1.25, size = 2.9, colour = "grey35", lineheight = 1.15) +
  scale_colour_manual(values = c("Positive (q < 0.1)"  = "#A50026",
                                 "Negative (q < 0.1)"  = "#313695",
                                 "Positive (p < 0.05)" = "#F4A582",
                                 "Negative (p < 0.05)" = "#92C5DE",
                                 "Not significant"     = "grey80"),
                      name = NULL, drop = FALSE) +
  scale_size_manual(values = c(3.8, 3.8, 2.6, 2.6, 1.7), guide = "none") +
  scale_alpha_manual(values = c(1, 1, 0.9, 0.9, 0.65), guide = "none") +
  scale_x_continuous(limits = c(-xmax, xmax)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.20))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 3.2))) +
  labs(x = "Coefficient (log2 fold change per SD of animal protein density)",
       y = expression(-log[10] * "(p)"),
       title = "Genera associated with animal protein density",
       subtitle = "MaAsLin 3 abundance model; TSS + log2, adjusted for age and gender") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

## B ---------------------------------------------------------------------
centr <- readRDS(file.path(PATHS$results, "networks.rds"))$centrality
pB <- centr |>
  mutate(label = ifelse(rank(-abs(delta)) <= 12, taxon, NA)) |>
  ggplot(aes(plant, animal, colour = delta)) +
  geom_abline(slope = 1, linetype = 2, colour = "grey60") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  scale_colour_gradient2(low = "#4F8A5B", mid = "grey85", high = "#C0504D",
                         midpoint = 0, name = "Difference") +
  coord_equal() +
  labs(x = "Eigenvector centrality — high plant",
       y = "Eigenvector centrality — high animal",
       title = "Taxon centrality in each network",
       caption = "Above the diagonal: more central in the high animal network.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey40", hjust = 0, size = 8))

## C ---------------------------------------------------------------------
alpha <- readRDS(file.path(PATHS$results, "alpha_diversity.rds"))
src_lab <- c(animal_pct_energy = "Animal protein", plant_pct_energy = "Plant protein")

alpha_long <- alpha$values |>
  select(sample_id, all_of(PRIMARY), all_of(METRICS)) |>
  pivot_longer(all_of(METRICS), names_to = "metric", values_to = "value") |>
  pivot_longer(all_of(PRIMARY), names_to = "exposure", values_to = "intake") |>
  mutate(metric = factor(metric, METRICS, M_LAB[METRICS]),
         exposure = factor(exposure, PRIMARY, src_lab[PRIMARY]))

ann <- alpha$spearman |>
  filter(exposure %in% PRIMARY) |>
  mutate(metric = factor(metric, METRICS, M_LAB[METRICS]),
         exposure = factor(exposure, PRIMARY, src_lab[PRIMARY]),
         lab = sprintf("rho = %+.2f, p = %.3f", rho, p),
         face = ifelse(p < 0.05, "bold", "plain"))

pC <- ggplot(alpha_long, aes(intake, value)) +
  geom_point(aes(colour = exposure), size = 1.7, alpha = 0.65) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "grey15", linewidth = 0.8) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab, fontface = face),
            hjust = -0.05, vjust = 1.6, size = 2.9, colour = "grey25",
            inherit.aes = FALSE) +
  facet_grid(metric ~ exposure, scales = "free") +
  scale_colour_manual(values = SRC_PAL) +
  labs(title = "Alpha diversity across animal and plant protein intake",
       x = "% of total energy", y = "Alpha diversity",
       caption = "LOESS fit; Spearman rho per panel, bold if p < 0.05.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey40", hjust = 0, size = 8))

## Assemble --------------------------------------------------------------
figure2 <- plot_grid(
  plot_grid(pA, pB, ncol = 1, labels = c("A", "B"), label_size = 16),
  plot_grid(pC, ncol = 1, labels = "C", label_size = 16),
  ncol = 2, rel_widths = c(1.15, 1))

ggsave(file.path(PATHS$figures, "figure2.png"), figure2,
       width = 14, height = 13, dpi = 300, bg = "white")
ggsave(file.path(PATHS$figures, "figure2.pdf"), figure2,
       width = 14, height = 13, device = cairo_pdf, bg = "white")
