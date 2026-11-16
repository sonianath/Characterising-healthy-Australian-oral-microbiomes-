# Supplementary Figure 1. Macronutrient intake against the Acceptable
# Macronutrient Distribution Ranges for Australian adults (NHMRC, 2006):
# protein 15-25%, fat 20-35%, carbohydrate 45-65% of energy.
# Energy factors (FSANZ Schedule 11): protein 17, fat 37, carbohydrate 17 kJ/g.

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

dat <- read_participants()

kj_per_g <- c(Protein = 17, Fat = 37, Carbohydrate = 17)
amdr <- data.frame(macro = names(kj_per_g), lo = c(15, 20, 45), hi = c(25, 35, 65))

pct_e <- data.frame(
  participant_id = rep(dat$participant_id, 3),
  macro = rep(names(kj_per_g), each = nrow(dat)),
  pct   = 100 * c(dat$protein_g      * kj_per_g["Protein"],
                  dat$fat_g          * kj_per_g["Fat"],
                  dat$carbohydrate_g * kj_per_g["Carbohydrate"]) /
          rep(dat$energy_kj, 3))
pct_e <- merge(pct_e, amdr, by = "macro")
pct_e$status <- factor(ifelse(pct_e$pct < pct_e$lo, "Below range",
                       ifelse(pct_e$pct > pct_e$hi, "Above range", "Within range")),
                       levels = c("Below range", "Within range", "Above range"))

macro_lab <- c(Protein = "Protein\n(AMDR 15–25%)",
               Fat = "Fat\n(AMDR 20–35%)",
               Carbohydrate = "Carbohydrate\n(AMDR 45–65%)")
pct_e$macro_f <- factor(macro_lab[pct_e$macro], levels = rev(macro_lab))
amdr$macro_f  <- factor(macro_lab[amdr$macro],  levels = rev(macro_lab))

amdr_tab <- as.data.frame(table(macro = pct_e$macro, status = pct_e$status))
amdr_tab$pct <- 100 * amdr_tab$Freq / nrow(dat)
amdr_tab <- amdr_tab[order(amdr_tab$macro, amdr_tab$status), ]
write_table(transform(amdr_tab, pct = round(pct, 1)), "amdr_classification")
print(round(tapply(pct_e$pct, pct_e$macro, mean), 1))

status_col <- c("Below range" = "#3B75AF", "Within range" = "#2A9D8F",
                "Above range" = "#E07A3F")
amdr_tab$macro_f <- factor(macro_lab[as.character(amdr_tab$macro)],
                           levels = rev(macro_lab))
amdr_tab$label <- ifelse(amdr_tab$pct >= 8,
                         sprintf("%.0f%%\n(n = %d)", amdr_tab$pct, amdr_tab$Freq), "")

pA <- ggplot(amdr_tab, aes(pct, macro_f, fill = status)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.4,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            colour = "white", size = 3, lineheight = 0.9, fontface = "bold") +
  scale_fill_manual(values = status_col, drop = FALSE, name = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100.5),
                     breaks = seq(0, 100, 25), labels = \(x) paste0(x, "%")) +
  labs(x = "Participants", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black"))

set.seed(1)
pB <- ggplot(pct_e, aes(pct, macro_f)) +
  geom_tile(data = amdr, aes(x = (lo + hi) / 2, y = macro_f, width = hi - lo),
            height = 0.84, inherit.aes = FALSE, fill = "grey90") +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, colour = "grey35") +
  geom_jitter(aes(colour = status), height = 0.15, width = 0, size = 1.3,
              alpha = 0.8) +
  scale_colour_manual(values = status_col, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, 10),
                     labels = \(x) paste0(x, "%")) +
  labs(x = "Percentage of total energy intake", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = "black"))

fig_s1 <- pA / pB + plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.15), guides = "collect") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))

ggsave(file.path(PATHS$figures, "supp_figure1_amdr.png"), fig_s1,
       width = 180, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(PATHS$figures, "supp_figure1_amdr.pdf"), fig_s1,
       width = 180, height = 150, units = "mm", device = cairo_pdf, bg = "white")
