######### TABLE 1 STUDY PARTICIPANT CHARACTERISTIC ##################################
########## K10 SCORES AND DEMOGRAPHIC VARIABES #############
####################################################################
################# WRITTEN BY SONIA NATH #################
################# DATE 15-FEB-2026 #################
### STEPS FOR ANALYISIS ##
#1. CREATING A PHYLOSEQ OBJECT AND SUMMARISING PHYLOSEQ OBJECT
#3. ALPHA DIVERSITY AFTER RAREFACTION (LINEAR REGRESSION MODELS AND GROUP COMPARISON)
#4. BETA DIVERSITY (AITCHISON'S DISTANCE) AND RUNNING PERMANOVA USING ADONIS BY TERMS AGE, SEX, SMOKING, PHYSICAL EXERCISE
#5. PARTIAL CORDINATE ORDINATION PLOT
#6. DIFFERENTIAL ABUNDANCE ANALYSIS

#Clear existing data and graphics
rm(list=ls())
graphics.off()

#SETTING THE DIRECTORY
setwd("/K10_OMT/R_analysis")
####################################################################
# 1) CHECKING THE PHYLOSEQ OBJECT
####################################################################
#Library the packages
library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.
library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(plotly) # A package to create interactive web graphics of use in 3D plots
library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.
library(philr) # This package provides functions for the analysis of compositional data 
library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step
library(adespatial) # Tools for the multiscale spatial analysis of multivariate data
library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks
library(qiime2R) # A package for importing qiime artifacts into an R session
library(MicrobeR) # Data visualization
library(microbiome) # Data analysis and visualization
library(grid) # support data visualization
library(gridExtra)  # support data visualization
library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.
library(png) # Figure download
library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'
library(ggpubr) # publication quality figures, based on ggplot2
library(RColorBrewer) # nice color options
library(DT) #for interactive tables
library(reshape2)
library(scales)
library(data.table)
library(Biostrings)
library(readr)
library(dplyr)
library(tibble)


#A) CREATING A PHYLOSEQ OBJECT
physeq <- qza_to_phyloseq(
  features  = "table.qza",          # FeatureTable[Frequency]
  tree      = "rooted-tree.qza",    # Phylogeny[Rooted]
  taxonomy  = "taxonomy.qza",       # FeatureData[Taxonomy]
  metadata  = "metadata.tsv"     # QIIME 2 metadata (CHECK first col = SampleID)
)

ps <- physeq #creating a copy-> using ps for analysis

# B) SUMMARISING THE PHYLOSEQ OBJECT
summarize_phyloseq(physeq) #gives the brief summary of the samples and the data
print(physeq) #prints the same summary
summary(sample_sums(physeq)) #fives the min, max, median, mean for the reads
sample_names(physeq)[1:5] #overview of sample names 
rank_names(physeq)  #Rank names
sample_variables(physeq)  # the variables in the metadata

otu_table(physeq)[1:5, 1:5] ## the feature table 
# -----------------------------
# B) Rarefy to 3000 reads
# -----------------------------
# If any samples < 3000, dropping them first
ps_qc <- prune_samples(sample_sums(ps) >= 3000, ps)
ps_qc <- prune_taxa(taxa_sums(ps_qc) > 0, ps_qc)
summary(sample_sums(ps_qc)) 

ps_rar<- rarefy_even_depth(
  ps_qc,
  sample.size = 8676,
  rngseed = 123,
  replace = FALSE,
  verbose = FALSE
)

summarize_phyloseq(ps_rar)
print(ps_rar)
summary(sample_sums(ps_rar))
# -----------------------------
# C) Alpha diversity table + metadata
# -----------------------------
alpha_measures <- c("Observed", "Shannon", "Simpson")

richness <- estimate_richness(ps_rar, measures = alpha_measures)

meta <- data.frame(sample_data(ps_rar), check.names = FALSE)

alpha_diversity <- cbind(meta, richness)

# Save files 
write.table(richness,
            file = "richness.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(alpha_diversity,
            file = "alpha_diversity.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# -----------------------------
# D) Quick distribution checks
# -----------------------------
par(mfrow = c(1, 3))
hist(alpha_diversity$Shannon,  main = "Shannon",  xlab = "")
hist(alpha_diversity$Observed, main = "Observed", xlab = "")
hist(alpha_diversity$Simpson,  main = "Simpson",  xlab = "")
par(mfrow = c(1, 1))

# Shapiro-Wilk 
shapiro.test(alpha_diversity$Shannon)
shapiro.test(alpha_diversity$Observed)
shapiro.test(alpha_diversity$Simpson)


# ============================================================
# Robust LM (HC3) for alpha diversity ~ K10 (continuous)
# + plots for Observed, Shannon, Simpson
# ============================================================

library(dplyr)
library(broom)
library(lmtest)
library(sandwich)
library(ggplot2)

# -----------------------------
# 1) Prepare analysis dataset
# -----------------------------
dat_k10 <- alpha_diversity %>%
  mutate(
    gender.factor = factor(gender.factor),
    age_s = as.numeric(scale(age)),
    k10_s = as.numeric(scale(k10score))
  ) %>%
  # complete cases for everything used in models
  filter(complete.cases(Observed, Shannon, Simpson, k10_s, age_s, gender.factor))


dat_k10 <- alpha_diversity %>%
  mutate(
    gender.factor = factor(gender.factor),
    smoking_cat   = factor(smoking_cat),
    age_s = as.numeric(scale(as.numeric(as.character(age)))),
    k10_s = as.numeric(scale(as.numeric(as.character(k10score))))
  ) %>%
  filter(complete.cases(
    Observed, Shannon, Simpson,
    k10_s, age_s, gender.factor,
    smoking_cat
  ))%>%
  # complete cases for everything used in models
  filter(complete.cases(Observed, Shannon, Simpson, k10_s, age_s, gender.factor))


# -----------------------------
# 2) Robust regression function (HC3)
# -----------------------------
fit_one_k10 <- function(outcome, data) {
  fml <- as.formula(paste(outcome, "~ k10_s + age_s + gender.factor + smoking_cat"))
  
  m <- lm(fml, data = data)
  
  rob <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type = "HC3"))
  
  broom::tidy(rob) %>%
    filter(term == "k10_s") %>%
    mutate(
      outcome = outcome,
      exposure = "k10score (z)",
      n = nrow(data)
    ) %>%
    select(outcome, exposure, estimate, std.error, statistic, p.value, n)
}

# -----------------------------
# 3) Run models for all 3 indices
# -----------------------------
alpha_measures <- c("Observed", "Shannon", "Simpson")

results_k10 <- bind_rows(lapply(alpha_measures, fit_one_k10, data = dat_k10)) %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value)

results_k10

# Save results
write.table(results_k10,
            file = "alpha_lm_robust_results_K10_HC3.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# -----------------------------
# 4) Plots: alpha diversity vs K10 (z)
# -----------------------------
plot_one <- function(df, yvar, ylab, file_stub) {
  p <- ggplot(df, aes(x = k10_s, y = .data[[yvar]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    labs(
      title = paste0(ylab, " vs K10 score"),
      x = "K10 score (z-score)",
      y = ylab
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(paste0(file_stub, ".png"), p, width = 6.5, height = 4.5, dpi = 300, bg = "white")
  ggsave(paste0(file_stub, ".svg"), p, width = 6.5, height = 4.5, device = "svg", bg = "white")
  
  p
}

p_obs <- plot_one(dat_k10, "Observed", "Observed features", "observed_vs_k10_z")
p_sha <- plot_one(dat_k10, "Shannon",  "Shannon diversity", "shannon_vs_k10_z")
p_sim <- plot_one(dat_k10, "Simpson",  "Simpson index",     "simpson_vs_k10_z")

p_obs; p_sha; p_sim

# Combine: 3 plots in one row (change to (p_obs / p_sha / p_sim) for one column)
library(patchwork)
p_alphaline <- p_obs | p_sha | p_sim
p_alphaline
# Save as SVG
ggsave(
  filename = "alpha_diversity_k10_alphaline.svg",
  plot     = p_alphaline,
  width    = 15,   # adjust as you like
  height   = 5,
  bg       = "white"
)

################################################
############ ALPHA DIVERSITY PLOTS###############
############### by K10 categories ###############
################################################
# ============================================================
# F) Alpha diversity plots by K10 categories (No distress vs Distress)
#    (stable: use alpha_diversity long format instead of plot_richness)
# ============================================================

comps <- list(c("No distress", "Distress"))

pal <- c("No distress" = "#009E73",
         "Distress"    = "#D55E00")

symnum.args <- list(
  cutpoints = c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1),
  symbols   = c("****", "***", "**", "*", "ns")
)

dat_cat <- alpha_diversity %>%
  filter(!is.na(k10_comb)) %>%
  mutate(
    k10_comb = factor(as.character(k10_comb),
                      levels = c("No distress", "Distress"),
                      ordered = TRUE)
  ) %>%
  pivot_longer(cols = all_of(alpha_measures),
               names_to = "measure",
               values_to = "alpha") %>%
  filter(!is.na(alpha))

p_k10 <- ggplot(dat_cat, aes(x = k10_comb, y = alpha, color = k10_comb, fill = k10_comb)) +
  geom_boxplot(color = "black", alpha = 0.65, width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.7) +
  facet_wrap(~ measure, scales = "free_y") +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(
    x = "Psychological distress (K10)",
    y = "Alpha diversity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  ) +
  # With 2 groups, Wilcoxon is the natural choice; runs per facet
  stat_compare_means(
    comparisons   = comps,
    method        = "wilcox.test",
    label         = "p.signif",
    symnum.args   = symnum.args,
    step.increase = 0.08
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label  = "p.format",
    label.y.npc = 0.98
  )

print(p_k10)

ggsave("alpha_diversity_by_k10_comb.png", p_k10, width = 8.5, height = 4.8, dpi = 300, bg = "white")
ggsave("alpha_diversity_by_k10_comb.svg", p_k10, width = 8.5, height = 4.8, device = "svg", bg = "white")


# -----------------------------
#  FIGURE 1
# -----------------------------
p_alphaline
p_k10

library(patchwork)

# Stack: A on top, B on bottom
Figure_1 <- ( p_k10/p_alphaline ) +
  plot_annotation(tag_levels = "A")  
# View
Figure_1 

# Save (edit size if you want)
ggsave("Figure_1.svg", Figure_1, width = 8.5, height = 12, bg = "white")
ggsave("Figure_1.png", Figure_1, width = 8.5, height = 12, dpi = 300, bg = "white")

####################################################
############### Setting up for Microviz package###########
###################################################
library(microbiome)
library(ggplot2)
library(dplyr)
library(vegan)
library(microViz)
##Fixing the tables
ps1<- ps
tax_table(ps1) %>% head(3)
phyloseq_validate(ps1)
#tax_fix_interactive(ps1)

ps1 %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

ps3 <- tax_fix(ps1)  
ps3 %>% tax_fix(unknowns = c("uncultured"))
ps3 <- ps3 %>% tax_fix(unknowns = c("uncultured"))

#####################  #####################  ##################### 
########### Beta diversity ordination plot (Aitchison + NMDS) ####
#####################  #####################  ##################### 
aitchison_dists <- ps3 %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc("aitchison")

perm_res <- aitchison_dists %>%
  dist_permanova(
    seed        = 1234,
    n_processes = 2,
    n_perms     = 999,
    variables   = "k10_comb"  
  )
perm_res

################################################
############ PERMANOVA ###############
############### by K10 categories ###############
################################################
library(phyloseq)
library(vegan)

# ---- inputs ----
covars    <- c("age", "gender", "smoking", "total_met_minutes", "Energy..exc.fibre", "Water")
expo_vars <- c("k10score")

# ---- 0) Start from ASV-level phyloseq object (NO FILTERING) ----
ps_asv <- ps

# ---- 1) Build samples x taxa matrix (NO FILTERING) ----
X <- as(otu_table(ps_asv), "matrix")
if (taxa_are_rows(ps_asv)) X <- t(X)    # samples x taxa
if (is.null(rownames(X))) rownames(X) <- sample_names(ps_asv)

# ---- 2) Aitchison distance = Euclidean on CLR ----
X <- X + 1
X <- X / rowSums(X)

logX <- log(X)
clrX <- logX - rowMeans(logX)

dist_ait_asv <- dist(clrX, method = "euclidean")
Dmat <- as.matrix(dist_ait_asv)

# ---- 3) Metadata ----
meta <- data.frame(sample_data(ps_asv), check.names = FALSE)
rownames(meta) <- sample_names(ps_asv)

# ---- Coerce covariates safely ----
meta$gender <- factor(meta$gender)

# smoking: 1 = No smoking, 2 = Smoking
meta$smoking <- factor(
  meta$smoking,
  levels = c(1, 2),
  labels = c("No smoking", "Smoking")
)

meta$age <- suppressWarnings(as.numeric(as.character(meta$age)))
meta$total_met_minutes <- suppressWarnings(as.numeric(as.character(meta$total_met_minutes)))

# dietary covariates (continuous)
meta$`Energy..exc.fibre` <- suppressWarnings(as.numeric(as.character(meta$`Energy..exc.fibre`)))
meta$water               <- suppressWarnings(as.numeric(as.character(meta$water)))

# ---- helper: run one model and return tidy rows (covars + exposure) ----
run_perm_tidy_asv <- function(expo, meta, Dmat, covars,
                              permutations = 999, seed = 1234, kept_ASVs = NA_integer_) {
  
  m0 <- meta
  
  # Ensure k10score is numeric (even if stored as character/factor)
  if (expo %in% c("k10score")) {
    m0[[expo]] <- suppressWarnings(as.numeric(as.character(m0[[expo]])))
  } else {
    if (is.factor(m0[[expo]]) || is.character(m0[[expo]])) {
      m0[[expo]] <- factor(m0[[expo]])
    } else {
      m0[[expo]] <- suppressWarnings(as.numeric(m0[[expo]]))
    }
  }
  
  keep <- complete.cases(m0[, c(covars, expo)])
  m2   <- m0[keep, , drop = FALSE]
  
  if (nrow(m2) < 5) {
    out <- data.frame(
      exposure  = expo,
      n         = nrow(m2),
      kept_ASVs = kept_ASVs,
      term      = c(covars, expo),
      Df        = NA, SumOfSqs = NA, R2 = NA, F = NA, `Pr(>F)` = NA,
      check.names = FALSE
    )
    return(out)
  }
  
  # subset distance to those samples 
  sam <- rownames(m2)
  D2  <- as.dist(Dmat[sam, sam])
  
  # fit model: covars first, exposure last (sequential)
  fml <- as.formula(paste("D2 ~", paste(c(covars, expo), collapse = " + ")))
  set.seed(seed)
  ad <- adonis2(fml, data = m2, permutations = permutations, by = "terms")
  
  tab <- as.data.frame(ad)
  tab$term <- rownames(tab)
  
  tab <- tab[tab$term %in% c(covars, expo),
             c("term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)"),
             drop = FALSE]
  
  tab$exposure  <- expo
  tab$n         <- nrow(m2)
  tab$kept_ASVs <- kept_ASVs
  
  tab[, c("exposure", "n", "kept_ASVs", "term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
}

# ---- 4) Run exposure(s) + bind ----
kept_ASVs <- ntaxa(ps_asv)

tidy_all_asv <- do.call(
  rbind,
  lapply(seq_along(expo_vars), function(i) {
    run_perm_tidy_asv(
      expo         = expo_vars[i],
      meta         = meta,
      Dmat         = Dmat,
      covars       = covars,
      permutations = 999,
      seed         = 1234 + i - 1,
      kept_ASVs    = kept_ASVs
    )
  })
)

# Exposure-only summary
tidy_expo_only_asv <- tidy_all_asv[tidy_all_asv$term == tidy_all_asv$exposure, ]
tidy_expo_only_asv <- tidy_expo_only_asv[order(tidy_expo_only_asv$`Pr(>F)`), ]

# (Optional) BH-FDR across exposures if you have >1 exposure
# tidy_expo_only_asv$p_adj_BH <- p.adjust(tidy_expo_only_asv$`Pr(>F)`, method = "BH")

# ---- 5) View ----
tidy_all_asv
tidy_expo_only_asv

# ---- 6) Save ----
write.csv(tidy_all_asv,
          "permanova_ASV_Aitchison_tidy_all_terms_by_terms.csv",
          row.names = FALSE)

write.csv(tidy_expo_only_asv,
          "permanova_ASV_Aitchison_tidy_exposures_only_by_terms.csv",
          row.names = FALSE)

################################################
############ RDA ANALYSIS ###############
############### ###############
################################################
library(phyloseq)
library(microViz)
library(dplyr)
library(ggplot2)

# -----------------------------
# 1) Prepare metadata
# -----------------------------
ps_base <- physeq %>%
  ps_mutate(
    # ---- numeric base ----
    age_num   = suppressWarnings(as.numeric(as.character(age))),
    k10_score = suppressWarnings(as.numeric(as.character(k10score))),
    met       = suppressWarnings(as.numeric(as.character(total_met_minutes))),
    flow      = suppressWarnings(as.numeric(as.character(saliva_flow_rate))),
    glu       = suppressWarnings(as.numeric(as.character(post.glucose_2))),
    filled    = suppressWarnings(as.numeric(as.character(filled))),
    smoking      = suppressWarnings(as.numeric(as.character(smoking))),
    carb      = suppressWarnings(as.numeric(as.character(Carbohydrate))),
    water     = suppressWarnings(as.numeric(as.character(Water))),
    energy    = suppressWarnings(as.numeric(as.character(Energy..exc.fibre))),
    
    # for plotting (already exists)
    gender = factor(as.character(gender)),
    
    # ---- condition ----
    age_s = as.numeric(scale(age_num)),
    
    # ---- z-scored continuous constraints (vectors) ----
    k10_z    = as.numeric(scale(k10_score)),
    met_z    = as.numeric(scale(met)),
    flow_z   = as.numeric(scale(flow)),
    glu_z    = as.numeric(scale(glu)),
    filled_z = as.numeric(scale(filled)),
    carb_z   = as.numeric(scale(carb)),
    water_z  = as.numeric(scale(water)),
    energy_z = as.numeric(scale(energy))
  ) %>%
  ps_filter(
    !is.na(gender),
    !is.na(age_s),
    !is.na(smoking),
    !is.na(k10_z), !is.na(met_z), !is.na(flow_z), !is.na(glu_z),
    !is.na(filled_z), !is.na(carb_z), !is.na(water_z), !is.na(energy_z)
  )

# -----------------------------
# 2) Aitchison distance
# -----------------------------
ait_full <- ps_base %>%
  tax_transform("identity") %>%
  dist_calc("aitchison")

# -----------------------------
# 3) PERMANOVA (by terms) + partial dbRDA
#    Order matters when by="terms": conditions first, then constraints
# -----------------------------
perm_full <- ait_full %>%
  dist_permanova(
    variables = c(
      # conditions first (partialled out / adjusted)
      "age_s", "smoking","met_z", 
      # continuous constraints after
      "k10_z", "flow_z", "glu_z",
      "filled_z", "carb_z", "water_z", "energy_z"
    ),
    n_perms = 999, seed = 321, by = "terms"
  )

ord_full <- perm_full %>%
  ord_calc(
    # constraints = continuous vectors on plot
    constraints = c(
      "k10_z", "flow_z", "glu_z",
      "filled_z", "carb_z", "water_z", "energy_z"
    ),
    # conditions = adjusted/partialled out
    conditions  = c("age_s", "smoking", "met_z")
  )

# -----------------------------
# 4) Plot
# -----------------------------
p_full <- ord_full %>%
  ord_plot(
    colour = "gender.factor",
    size   = 3.2,
    alpha  = 0.6,
    auto_caption = 7,
    constraint_vec_length = 4.2,
    constraint_vec_style  = vec_constraint(1.2, colour = "grey25"),
    constraint_lab_style  = constraint_lab_style(
      max_angle = 80, size = 3.2, aspect_ratio = 0.9, colour = "grey10"
    )
  ) +
  stat_ellipse(aes(colour = gender.factor), linewidth = 0.4, linetype = 1) +
  scale_colour_brewer(palette = "Set2", name = "Age group") +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.2, colour = "grey50") +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.2, colour = "grey50") +
  coord_fixed(ratio = 0.9, clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL)

print(p_full)

ggsave("ord_full_partial_dbRDA_agegroup.png", p_full, width = 8.5, height = 6.5, dpi = 300, bg = "white")
ggsave("ord_full_partial_dbRDA_agegroup.svg", p_full, width = 8.5, height = 6.5, bg = "white")

################################################
############ BETA DIVERSITY PLOTS###############
############### by K10 categories ###############
################################################

library(phyloseq)
library(vegan)
library(ggplot2)

# =========================
# Inputs
# =========================
ps_obj <- ps3
# Same colors as above
pal <- c(
  "No distress" = "#009E73",
  "Distress"    = "#D55E00"
)

# =========================
# 1) Keep samples with non-missing k10_comb and enforce level order
# =========================
ps_plot <- prune_samples(!is.na(sample_data(ps_obj)$k10_comb), ps_obj)

sample_data(ps_plot)$k10_comb <- factor(
  as.character(sample_data(ps_plot)$k10_comb),
  levels = c("No distress", "Distress"),
  ordered = TRUE
)
sample_data(ps_plot)$k10_comb <- droplevels(sample_data(ps_plot)$k10_comb)

# =========================
# 2) Build samples x taxa matrix
# =========================
X <- as(otu_table(ps_plot), "matrix")
if (taxa_are_rows(ps_plot)) X <- t(X)  # samples x taxa
if (is.null(rownames(X))) rownames(X) <- sample_names(ps_plot)

# Safety: remove any samples with zero total counts (rare but can happen)
rs <- rowSums(X)
if (any(rs == 0)) {
  keep_samp <- names(rs)[rs > 0]
  ps_plot <- prune_samples(keep_samp, ps_plot)
  X <- as(otu_table(ps_plot), "matrix")
  if (taxa_are_rows(ps_plot)) X <- t(X)
}

# =========================
# 3) Aitchison distance: Euclidean on CLR
#    (pseudocount +1 to handle zeros)
# =========================
X <- X + 1
X <- X / rowSums(X)

logX <- log(X)
clrX <- logX - rowMeans(logX)  # CLR per sample

dist_ait <- dist(clrX, method = "euclidean")

# =========================
# 4) PCoA (cmdscale) + % variance
# =========================
pcoa <- cmdscale(dist_ait, k = 2, eig = TRUE)

eig <- pcoa$eig
eig_pos <- eig[eig > 0]
pct <- eig_pos / sum(eig_pos) * 100

ord_df <- data.frame(
  sample_id = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)

meta <- data.frame(sample_data(ps_plot), check.names = FALSE)
meta$sample_id <- rownames(meta)

plot_df <- merge(ord_df, meta, by = "sample_id", all.x = TRUE)

# =========================
# 5) Plot 
# =========================
p_beta_k10 <- ggplot(plot_df, aes(x = PC1, y = PC2, color = k10_comb)) +
  geom_point(size = 2.4, alpha = 0.85) +
  stat_ellipse(aes(group = k10_comb, color = k10_comb),
               type = "t", linetype = "dashed", linewidth = 0.9, level = 0.95) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(
    title = "PCoA (Aitchison distance) by K10 distress group",
    x = sprintf("PC1 (%.1f%%)", pct[1]),
    y = sprintf("PC2 (%.1f%%)", pct[2]),
    color = "K10 group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

print(p_beta_k10)

ggsave("beta_pcoa_aitchison_k10_comb.png", p_beta_k10, width = 7, height = 5, dpi = 300, bg = "white")
ggsave("beta_pcoa_aitchison_k10_comb.svg", p_beta_k10, width = 7, height = 5, device = "svg", bg = "white")

#####################  #####################  ##################### 
################# Differential Abundance Analysis ###############
#####################  #####################  ##################### 
library(maaslin3)

set.seed(123)
df_input_metadata$k10_comb <- factor(df_input_metadata$k10_comb)
df_input_metadata$k10_comb <- 
  factor(df_input_metadata$k10_comb, levels = c('No distress','Distress' ))
fit_out <- maaslin3(
  input_data     = df_input_data,
  input_metadata = df_input_metadata,
  output         = "mas_k10",
  formula        = '~ k10_comb + age + gender.factor + smoking',
  normalization  = "NONE",
  transform      = "LOG",
  warn_prevalence  = FALSE,         
  evaluate_only    = "abundance",    
  augment          = TRUE,
  standardize      = TRUE,
  min_prevalence = 0,
  median_comparison_abundance  = FALSE,
  median_comparison_prevalence = FALSE,
  max_significance = 0.1,
  max_pngs         = 50,
  cores = 1
)

# -----------------------------
# Plot
# -----------------------------
library(ggplot2)
library(dplyr)
library(patchwork)

results <- read.csv("mas3_k10_results.csv")

scale_x_continuous(limits = c(-3500, 0))

da_plot <- ggplot(results_clean, aes(x = coef, y = feature_clean, fill = fill_col)) +
  geom_col() +
  geom_errorbarh(
    aes(xmin = coef - 1.96 * stderr, xmax = coef + 1.96 * stderr),
    height = 0.3,
    colour = "grey30",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = sig,
        x = coef + sign(coef) * 1.96 * stderr + sign(coef) * 50),
    size = 3.5
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_continuous(limits = c(-3500, 0)) +
  scale_fill_manual(values = c(
    "FDR q<0.1" = "#c0392b",
    "FDR q≥0.1" = "#5b9bd5"
  )) +
  labs(
    x = "CLR coefficient",
    y = NULL,
    title = "Differential Abundance Analysis",
    subtitle = "Negative = lower in Distress group",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
da_plot
####################
#####FIGURE 2######
###################
p_full
p_beta_k10
da_plot

library(patchwork)

# Left = C (da_plot spans full height)
# Right = A (p_full) above B (p_beta_k10), with A taller than B
Figure_2 <- (p_full / p_beta_k10 + plot_layout(heights = c(1.5, 1))) | da_plot

# Tag order in this layout is: left plot, then top-right, then bottom-right
Figure_2 <- Figure_2  +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(plot.tag = element_text(face = "bold"))

Figure_2 

ggsave("Figure_2.svg", Figure_2, width = 12, height = 7, bg = "white")


