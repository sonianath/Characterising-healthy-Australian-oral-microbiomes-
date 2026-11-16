# Figure 1. Plaque community composition and dietary protein source.
#   A  Constrained partial ordination (Aitchison, ASV level)
#   B  PC1 of CLR abundances against animal protein density
#   C  PERMANOVA R2 for each exposure
#   D  Spearman correlations of the top 20 genera with each exposure

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(microViz); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(cowplot)
})

comp <- readRDS(file.path(PATHS$results, "community_composition.rds"))

## A ---------------------------------------------------------------------
pA <- comp$ord |>
  ord_plot(colour = "protein_group", size = 3.2, alpha = 0.6,
           auto_caption = 7, constraint_vec_length = 4.2,
           constraint_vec_style = vec_constraint(1.2, colour = "grey25"),
           constraint_lab_style = constraint_lab_style(
             max_angle = 80, size = 3.2, aspect_ratio = 0.9, colour = "grey10")) +
  stat_ellipse(aes(colour = protein_group), linewidth = 0.4) +
  scale_colour_manual(values = PAL, name = "Protein source group") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.2, colour = "grey50") +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.2, colour = "grey50") +
  coord_fixed(ratio = 0.9, clip = "off") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
        axis.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA))

## B ---------------------------------------------------------------------
rho <- comp$pc1_rho
pB <- ggplot(comp$pca_scores, aes(animal_pct_energy, PC1)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3, colour = "grey75") +
  geom_point(aes(colour = protein_group), size = 2.3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20",
              fill = "grey85", linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("rho = %+.2f, p = %.3f", rho$estimate, rho$p.value),
           hjust = -0.08, vjust = 1.6, size = 3.4, fontface = "bold") +
  scale_colour_manual(values = PAL, guide = "none") +
  labs(x = "Animal protein (% energy)",
       y = sprintf("PC1 (%.1f%%)", comp$pca_var[1]),
       title = "PC1 across animal protein density") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

## C ---------------------------------------------------------------------
pC <- comp$permanova |>
  mutate(source = ifelse(grepl("animal", exposure), "Animal protein",
                         "Plant protein"),
         lab = factor(EXPO_LAB[exposure], levels = rev(EXPO_LAB[EXPOSURES])),
         txt = sprintf("R² = %.2f%%   p = %.3f", 100 * R2, `Pr(>F)`),
         sig = `Pr(>F)` < 0.05) |>
  ggplot(aes(R2, lab, colour = source)) +
  geom_segment(aes(x = 0, xend = R2, yend = lab), linewidth = 1.6, alpha = 0.55) +
  geom_point(size = 5) +
  geom_text(aes(label = txt, fontface = ifelse(sig, "bold", "plain")),
            hjust = -0.12, size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = SRC_PAL, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.45))) +
  labs(x = expression("PERMANOVA " * R^2), y = NULL,
       title = "Community composition by protein exposure",
       subtitle = "Aitchison distance, ASV level; adjusted for age, gender and toothbrushing") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

## D ---------------------------------------------------------------------
cors <- comp$genus_cor
taxa_order <- cors |>
  select(term, taxon, estimate) |>
  pivot_wider(id_cols = taxon, names_from = term, values_from = estimate) |>
  column_to_rownames("taxon") |> as.matrix() |>
  dist() |> hclust(method = "ward.D2")
taxa_order <- taxa_order$labels[taxa_order$order]

pD <- cors |>
  mutate(taxon = factor(taxon, levels = taxa_order),
         term  = factor(term, levels = EXPOSURES, labels = EXPO_LAB[EXPOSURES])) |>
  ggplot(aes(term, taxon)) +
  geom_raster(aes(fill = estimate)) +
  geom_point(data = \(x) filter(x, p.value < 0.05), shape = 8, size = 3) +
  geom_point(data = \(x) filter(x, q < 0.05), shape = 16, size = 3) +
  scale_fill_distiller(palette = "BrBG", limits = c(-0.45, 0.45), oob = squish) +
  labs(x = NULL, y = NULL, fill = "Spearman\ncorrelation",
       caption = "Asterisk indicates p < 0.05 (unadjusted)\nFilled circle indicates FDR q < 0.05") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

## Assemble --------------------------------------------------------------
figure1 <- plot_grid(
  plot_grid(pA, pB, ncol = 1, labels = c("A", "B"), label_size = 16),
  plot_grid(pC, pD, ncol = 1, rel_heights = c(1, 2),
            labels = c("C", "D"), label_size = 16),
  ncol = 2, rel_widths = c(1, 1.05))

ggsave(file.path(PATHS$figures, "figure1.png"), figure1,
       width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(PATHS$figures, "figure1.pdf"), figure1,
       width = 14, height = 12, device = cairo_pdf, bg = "white")
