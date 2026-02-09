######### MICROBIOME ANALYSIS ##################################
########## ORAL MICROBIOME ANALYSIS: SUGAR TYPES IN FRUITS AND PROCESSED FOODS #############
####################################################################
################# WRITTEN BY SONIA NATH #################
################# DATE 28-JAN-2026 #################


### STEPS FOR ANALYISIS ##
#1. CREATING A PHYLOSEQ OBJECT AND SUMMARISING PHYLOSEQ OBJECT
#3. ALPHA DIVERSITY AFTER RAREFACTION (LINEAR REGRESSION MODELS AND GROUP COMPARISON)
#4. BETA DIVERSITY (AITCHISON'S DISTANCE) AND RUNNING PERMANOVA USING ADONIS BY TERMS AGE, SEX, TOOTHBRUSHING FREQUENCY
#5. PARTIAL CORDINATE ORDINATION PLOT
#6. DIFFERENTIAL ABUNDANCE ANALYSIS
###################################################
rm(list=ls())
graphics.off()
##################################
#SETTING THE DIRECTORY
setwd("~/Diet_sugar/data/R_analysis/")

########################################
# PARTCIPANT SUGAR CONSUMPTION FRUITS VS. PROCESSED
################################################
library(dplyr)

df_sugar <- dat %>%
  select(fruit_total_sugars_g_day, proc_total_sugars_g_day)

# Paired t-test (mean difference)
tt <- t.test(df_sugar$fruit_total_sugars_g_day,
             df_sugar$proc_total_sugars_g_day,
             paired = TRUE)

# Wilcoxon signed-rank (median difference, robust)
wt <- wilcox.test(df_sugar$fruit_total_sugars_g_day,
                  df_sugar$proc_total_sugars_g_day,
                  paired = TRUE, exact = FALSE)

# Effect size for paired t-test (Cohen's dz)
diffs <- with(df_sugar, fruit_total_sugars_g_day - proc_total_sugars_g_day)
cohen_dz <- mean(diffs) / sd(diffs)

tt
wt
cohen_dz

library(ggplot2)
tot_sugar<- ggplot(df_sugar, aes(x = proc_total_sugars_g_day, y = fruit_total_sugars_g_day)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Processed total sugars (g/day)",
       y = "Fruit total sugars (g/day)",
       title = "Within-person comparison: fruit vs processed sugars")
ggsave("tot_sugar_person.png", tot_sugar, width = 7, height = 5, dpi = 300, bg = "white")



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


#A) CREATING A PHYLOSEQ OBHECT
physeq <- qza_to_phyloseq(
  features  = "table.qza",          # FeatureTable[Frequency]
  tree      = "rooted-tree.qza",    # Phylogeny[Rooted]
  taxonomy  = "taxonomy.qza",       # FeatureData[Taxonomy]
  metadata  = "metadata_filled.tsv"     # QIIME 2 metadata (CHECK first col = SampleID)
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

####################################################################
### 2) ALPHA RAREFACTION PLOTS####
##### TRIMMING TO 3000 READS MIN.########
####################################################################
ps <- prune_samples(sample_sums(ps) >= 3000, ps) # losing only one sample. 92 in total
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# ASV matrix (samples x ASVs)
asv_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) asv_mat <- t(asv_mat)

depth_grid <- seq(100, 3000, by = 100)

rare_df <- lapply(seq_len(nrow(asv_mat)), function(i){
  counts <- asv_mat[i, ]
  tibble(
    sample   = rownames(asv_mat)[i],
    depth    = depth_grid,
    Observed = sapply(depth_grid, function(m) vegan::rarefy(counts, sample = m))
  )
}) |> bind_rows()

ggplot(rare_df, aes(depth, Observed, group = sample)) +
  geom_line(alpha = 0.5) +
  geom_vline(xintercept = 3000, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 3000)) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Expected richness (Observed ASVs)",
    title = "Alpha rarefaction curves (Observed ASVs)"
  ) +
  theme_minimal()


# 1) sequencing depth per sample
depth_df <- data.frame(
  sample = sample_names(ps),
  reads  = sample_sums(ps)
)

summary(depth_df$reads)

ggplot(depth_df, aes(x = reads)) +
  geom_histogram(bins = 30) +
  scale_x_log10() +
  labs(x = "Reads per sample (log10)", y = "Number of samples")

####################################################################
##### 3) LINEAR REGRESSION MODELS USING ALPHA DIVERSITY METRICS
##### A) OBSERVED FEATURES B) SHANNON'S IDEX C) SIMPSON'S INDEX ####
#####################################################################

library(phyloseq)
library(vegan)
library(dplyr)
library(broom)
library(lmtest)
library(sandwich)

set.seed(123)

# -----------------------------
# A) Inputs
# -----------------------------

covars <- c("age", "gender.factor", "toothbrsh_freq")

sugar_cols <- c(
  "fruit_total_sugars_g_day", "fruit_free_sugars_g_day",
  "fruit_glucose_g_day", "fruit_fructose_g_day", "fruit_sucrose_g_day",
  "proc_glucose_g_day", "proc_fructose_g_day", "proc_sucrose_g_day",
  "proc_total_sugars_g_day", "proc_added_sugars_g_day", "proc_free_sugars_g_day"
)

alpha_measures <- c("Observed", "Shannon", "Simpson")

# -----------------------------
# B) Rarefy to 3000 reads
# -----------------------------
# If any samples < 3000, dropping them first
ps_qc <- prune_samples(sample_sums(ps) >= 3000, ps)
ps_qc <- prune_taxa(taxa_sums(ps_qc) > 0, ps_qc)

ps_rar3000 <- rarefy_even_depth(
  ps_qc,
  sample.size = 3000,
  rngseed = 123,
  replace = FALSE,
  verbose = FALSE
)

# -----------------------------
# C) Alpha diversity table + metadata
# -----------------------------
richness <- estimate_richness(ps_rar3000, measures = alpha_measures)

meta <- data.frame(sample_data(ps_rar3000), check.names = FALSE)

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

# -----------------------------
# E) Prepare analysis dataset
#    - factor coding
#    - z-score sugars + age
#    -  median split groups (kept, but not used in lm below)
# -----------------------------
dat <- alpha_diversity %>%
  mutate(
    gender.factor  = factor(gender.factor),
    toothbrsh_freq = factor(toothbrsh_freq),
    age_s = as.numeric(scale(age))
  ) %>%
  mutate(
    across(all_of(sugar_cols), ~ as.numeric(.), .names = "{.col}")  # ensure numeric
  ) %>%
  mutate(
    # z-scored sugars for comparability across exposures
    across(all_of(sugar_cols), ~ as.numeric(scale(.)), .names = "{.col}_z")
  ) %>%
  mutate(
    # optional median-based groups on ORIGINAL (not z) sugars
    across(
      all_of(sugar_cols),
      ~ ifelse(is.na(.), NA_character_,
               ifelse(. > median(., na.rm = TRUE), "High", "Low")),
      .names = "{.col}_grp"
    )
  ) %>%
  mutate(
    across(ends_with("_grp"), ~ factor(., levels = c("Low", "High")))
  )

# keep complete cases for covars (exposure handled in loop)
dat_base <- dat %>%
  filter(complete.cases(across(all_of(c("age_s", "gender.factor", "toothbrsh_freq")))))

# -----------------------------
# 5) Run robust linear regression for ALL sugars x ALL alpha outcomes
#    Model: alpha ~ sugar_z + age_s + gender + brushing
#    Robust SE: HC3
# -----------------------------
fit_one <- function(outcome, sugar_z_col, data) {
  d <- data %>% filter(!is.na(.data[[outcome]]), !is.na(.data[[sugar_z_col]]))
  fml <- as.formula(paste(outcome, "~", sugar_z_col, "+ age_s + gender.factor + toothbrsh_freq"))
  m <- lm(fml, data = d)
  
  rob <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type = "HC3"))
  
  # extract the sugar term only
  term_row <- broom::tidy(rob) %>%
    filter(term == sugar_z_col) %>%
    mutate(
      outcome = outcome,
      sugar = gsub("_z$", "", sugar_z_col),
      n = nrow(d)
    )
  
  term_row
}

results <- bind_rows(lapply(alpha_measures, function(out){
  bind_rows(lapply(paste0(sugar_cols, "_z"), function(sz){
    fit_one(outcome = out, sugar_z_col = sz, data = dat_base)
  }))
})) %>%
  mutate(
    sugar_type = ifelse(grepl("^fruit_", sugar), "fruit", "processed")
  ) %>%
  group_by(outcome) %>%
  mutate(p_adj_BH = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(outcome, sugar_type, p.value)

results

# Save results
write.table(results,
            file = "alpha_lm_robust_results_HC3.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# -----------------------------
# ALPHA DIVERSITY PLOTS (using z-scored predictors - Line plots)
# Data frame: alpha_diversity
# -----------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# 1) Ensure the derived predictors exist (create if missing)
alpha_diversity <- alpha_diversity %>%
  mutate(
    # factors (optional)
    gender.factor  = factor(gender.factor),
    toothbrsh_freq = factor(toothbrsh_freq),
    
    # z-scored sugar predictors (create these if you haven't already)
    fruit_sugar_tot  = if (!"fruit_sugar_tot"  %in% names(.)) as.numeric(scale(fruit_total_sugars_g_day)) else fruit_sugar_tot,
    fruit_sugar_free = if (!"fruit_sugar_free" %in% names(.)) as.numeric(scale(fruit_free_sugars_g_day))  else fruit_sugar_free,
    fruit_sugar_glu  = if (!"fruit_sugar_glu"  %in% names(.)) as.numeric(scale(fruit_glucose_g_day))      else fruit_sugar_glu,
    fruit_sugar_fruc = if (!"fruit_sugar_fruc" %in% names(.)) as.numeric(scale(fruit_fructose_g_day))     else fruit_sugar_fruc,
    fruit_sugar_suc  = if (!"fruit_sugar_suc"  %in% names(.)) as.numeric(scale(fruit_sucrose_g_day))      else fruit_sugar_suc,
    
    proc_sugar_glu   = if (!"proc_sugar_glu"   %in% names(.)) as.numeric(scale(proc_glucose_g_day))       else proc_sugar_glu,
    proc_sugar_fruc  = if (!"proc_sugar_fruc"  %in% names(.)) as.numeric(scale(proc_fructose_g_day))      else proc_sugar_fruc,
    proc_sugar_suc   = if (!"proc_sugar_suc"   %in% names(.)) as.numeric(scale(proc_sucrose_g_day))       else proc_sugar_suc,
    proc_sugar_tot   = if (!"proc_sugar_tot"   %in% names(.)) as.numeric(scale(proc_total_sugars_g_day))  else proc_sugar_tot,
    proc_sugar_add   = if (!"proc_sugar_add"   %in% names(.)) as.numeric(scale(proc_added_sugars_g_day))  else proc_sugar_add,
    proc_sugar_free  = if (!"proc_sugar_free"  %in% names(.)) as.numeric(scale(proc_free_sugars_g_day))   else proc_sugar_free,
    
    age_s            = if (!"age_s"            %in% names(.)) as.numeric(scale(age))                       else age_s
  )

# 2) Predictors to plot (these are the ones you said you added)
predictors <- c(
  "fruit_sugar_tot","fruit_sugar_free","fruit_sugar_glu","fruit_sugar_fruc","fruit_sugar_suc",
  "proc_sugar_glu","proc_sugar_fruc","proc_sugar_suc","proc_sugar_tot","proc_sugar_add","proc_sugar_free",
  "age_s"
)


var_labs <- c(
  fruit_sugar_tot  = "Fruit sugars (total, z)",
  fruit_sugar_free = "Fruit sugars (free, z)",
  fruit_sugar_glu  = "Fruit glucose (z)",
  fruit_sugar_fruc = "Fruit fructose (z)",
  fruit_sugar_suc  = "Fruit sucrose (z)",
  proc_sugar_glu   = "Processed glucose (z)",
  proc_sugar_fruc  = "Processed fructose (z)",
  proc_sugar_suc   = "Processed sucrose (z)",
  proc_sugar_tot   = "Processed sugars (total, z)",
  proc_sugar_add   = "Processed sugars (added, z)",
  proc_sugar_free  = "Processed sugars (free, z)",
  age_s            = "Age (z)"
)

# 3) Long format for faceting
plot_data <- alpha_diversity %>%
  dplyr::select(dplyr::all_of(c("Observed", "Shannon", "Simpson", predictors))) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(predictors),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(variable = factor(variable, levels = predictors))


# 4) Plot helper
plot_facets <- function(df, yvar, ylab, title, file_stub) {
  p <- ggplot(df, aes(x = value, y = .data[[yvar]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, colour = "red") +
    facet_wrap(~ variable, scales = "free_x", labeller = as_labeller(var_labs)) +
    labs(title = title, x = "Predictor value", y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      strip.background  = element_rect(fill = "white", colour = NA)
    )
  
  ggsave(paste0(file_stub, ".png"), p, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(paste0(file_stub, ".svg"), p, width = 12, height = 8, device = "svg", bg = "white")
  p
}

# 5) Make + save plots
p_obs <- plot_facets(plot_data, "Observed", "Observed features",
                     "Observed features vs predictors (z-scored, bivariate LM)", "observed_vs_predictors_z")

p_sha <- plot_facets(plot_data, "Shannon", "Shannon diversity",
                     "Shannon diversity vs predictors (z-scored, bivariate LM)", "shannon_vs_predictors_z")

p_sim <- plot_facets(plot_data, "Simpson", "Simpson index",
                     "Simpson index vs predictors (z-scored, bivariate LM)", "simpson_vs_predictors_z")

p_obs; p_sha; p_sim # used for supplementary

# ==================================================
# FOREST PLOTS FROM THE LINEAR REGRESSION MODELS
# ==================================================
#of adjusted effects (HC3 robust SE)
# Facet: Outcome (rows) x Fruit/Processed (cols)
# Color: pastel red (fruit), pastel blue (processed)
# Star: p <= 0.05
# Data frame: alpha_diversity
# Predictors: fruit_sugar_*, proc_sugar_* (z-scored) + covars age_s, gender.factor, toothbrsh_freq
# =========================
library(dplyr)
library(ggplot2)
library(lmtest)
library(sandwich)
library(broom)
library(forcats)

# --- make sure factors are factors ---
alpha_diversity <- alpha_diversity %>%
  mutate(
    gender.factor  = factor(gender.factor),
    toothbrsh_freq = factor(toothbrsh_freq)
  )

# --- predictors and covariates (match your added columns) ---
sugars <- c(
  "fruit_sugar_tot","fruit_sugar_free","fruit_sugar_glu","fruit_sugar_fruc","fruit_sugar_suc",
  "proc_sugar_glu","proc_sugar_fruc","proc_sugar_suc","proc_sugar_tot","proc_sugar_add","proc_sugar_free"
)
covars <- c("age_s","gender.factor","toothbrsh_freq")
outcomes <- c("Observed","Shannon","Simpson")


sugar_labels <- c(
  fruit_sugar_tot  = "Fruit total",   fruit_sugar_free = "Fruit free",
  fruit_sugar_glu  = "Fruit glucose", fruit_sugar_fruc = "Fruit fructose", fruit_sugar_suc = "Fruit sucrose",
  proc_sugar_glu   = "Processed glucose", proc_sugar_fruc = "Processed fructose", proc_sugar_suc = "Processed sucrose",
  proc_sugar_tot   = "Processed total",   proc_sugar_add = "Processed added",     proc_sugar_free = "Processed free"
)

# --- fit all models and extract the sugar term with robust HC3 SE ---
fit_one <- function(y, x){
  
  dat <- alpha_diversity %>%
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(c(y, x, covars)))))
  
  fml <- stats::as.formula(paste(y, "~", paste(c(x, covars), collapse = " + ")))
  m <- stats::lm(fml, data = dat)
  
  ct <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type = "HC3"))
  td <- broom::tidy(ct) %>% dplyr::filter(term == x)
  
  td %>%
    dplyr::mutate(
      outcome = y,
      sugar = x,
      sugar_label = sugar_labels[x],
      group = ifelse(grepl("^fruit_", x), "Fruit", "Processed"),
      n = nrow(dat),
      ci_low  = estimate - 1.96 * std.error,
      ci_high = estimate + 1.96 * std.error,
      sig = dplyr::case_when(
        p.value <= 0.05 ~ "*",
        p.value <= 0.10 ~ ".",
        TRUE ~ ""
      )
    )
}  

res <- res %>%
  dplyr::group_by(outcome, group) %>%
  dplyr::mutate(
    # panel-specific x-max (match your axis limits)
    x_max = dplyr::case_when(
      outcome == "Observed" ~ 14,
      TRUE ~ 0.5
    ),
    # scale-aware offset for where to place symbols
    panel_rng = max(ci_high, na.rm = TRUE) - min(ci_low, na.rm = TRUE),
    
    # symbol x-position: put '*' further right than '.'
    sig_x = dplyr::case_when(
      sig == "*" ~ ci_high + 0.06 * panel_rng,   # further right
      sig == "." ~ ci_high + 0.06 * panel_rng,   # slightly right
      TRUE ~ NA_real_
    ),
    
    # cap so the symbol stays inside axis
    sig_x = pmin(sig_x, x_max - 0.03),
    
    # make '.' bigger than '*'
    sig_size = dplyr::case_when(
      sig == "." ~ 9,    # bigger dot
      sig == "*" ~ 7,    # big star
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::ungroup()


# =========================
# Plot 1: Observed only  (x: -14..14, ticks every 2)
# Plot 2: Shannon+Simpson (x: -5..5, ticks every 1)
# Then stack with patchwork 
# =========================

library(ggplot2)
library(dplyr)

# --- pastel colours ---
cols <- c(Fruit = "#F4A6A6", Processed = "#A6C8F4")

# --- common theme ---
base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

# --- safety: keep stars inside panel ranges ---
res_plot <- res %>%
  mutate(
    x_max = case_when(
      outcome == "Observed" ~ 10,
      TRUE ~ 0.5
    ),
    star_x = pmin(star_x, x_max - 0.3)
  )

# =========================
# PLOT 1: Observed only
# =========================
p_obs <- res_plot %>%
  filter(outcome == "Observed") %>%
  ggplot(aes(x = estimate, y = sugar_label, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.7) +
  geom_point(size = 3.6) +  # bigger dots
  geom_text(aes(x = star_x, label = sig), colour = "black", size = 7, vjust = 0.35) +
  facet_grid(. ~ group, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = cols) +
  scale_x_continuous(
    limits = c(-14, 8),
    breaks = seq(-14, 8, by = 2),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  labs(x = "Adjusted β (per 1 SD higher sugar)", y = NULL) +
  base_theme

# =========================
# PLOT 2: Shannon + Simpson
# =========================
p_ss <- res_plot %>%
  filter(outcome %in% c("Shannon", "Simpson")) %>%
  ggplot(aes(x = estimate, y = sugar_label, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.7) +
  geom_point(size = 3.6) +
  geom_text(aes(x = star_x, label = sig), colour = "black", size = 5, vjust = 0.35) +
  facet_grid(outcome ~ group, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = cols) +
  scale_x_continuous(
    limits = c(-0.3, 0.3),
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  labs(x = "Adjusted β (per 1 SD higher sugar)", y = NULL) +
  base_theme

p_obs
p_ss

# =========================
# STACK THEM
# =========================
library(patchwork)

p_stacked <- p_obs / p_ss + plot_layout(heights = c(1, 2))
p_stacked

ggsave("alpha_forest_stacked.png", p_stacked, width = 11, height = 11, dpi = 300, bg = "white")
ggsave("alpha_forest_stacked.svg", p_stacked, width = 11, height = 11, device = "svg", bg = "white")



# ==================================================
# PREPARING PHYLOSEQ OBJECT FOR MICROVIZ PACKAGE
# A) CONTRAINED PARTIAL ANALYSIS  B) CORRELATION PLOTS
# ==================================================
library(microbiome)
library(ggplot2)
library(dplyr)
library(vegan)
library(microViz)
##Fixing the tables
ps1<- physeq
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

# ========================================
# CONTRAINED PARTIAL ORDINATION
# ========================================

# Packages
library(phyloseq)
library(microViz)
library(dplyr)
library(ggplot2)
vars <- c(
  "fruit_total_sugars_g_day", "fruit_free_sugars_g_day",
  "fruit_glucose_g_day", "fruit_fructose_g_day", "fruit_sucrose_g_day",
  "proc_glucose_g_day", "proc_fructose_g_day", "proc_sucrose_g_day",
  "proc_total_sugars_g_day", "proc_added_sugars_g_day", "proc_free_sugars_g_day"
)
# --- 1) Aitchison distance at ASV level ---
ait <- ps3 %>%
  tax_transform("identity") %>%
  dist_calc("aitchison")

# --- 2) Prepare metadata: z-score sugars; coerce covariates ---

ait <- ait %>%
  ps_mutate(
    fruit_tot_sug = as.numeric((fruit_total_sugars_g_day)),
    
    
    proc_tot_sug   = as.numeric((proc_total_sugars_g_day)),
    
    post_glucose_salivary_pH   = as.numeric((ph_fall)),  
    Filled_teeth         = as.numeric((filled)),
    Missing_teeth        = as.numeric((missing)),
    physical_activity =  as.numeric((total_met_minutes)),               
    
    age_s          = as.numeric(age_s)
  )


perm_sug <- ait %>%
  dist_permanova(
    variables   = c("fruit_tot_sug",
                    "proc_tot_sug",
                    "Filled_teeth", "Missing_teeth", "post_glucose_salivary_pH",
                    "physical_activity",
                    "age_s"),
    n_perms     = 999, seed = 321, by = "margin"
  )

ord_sug <- perm_sug %>%
  ord_calc(constraints = c("fruit_tot_sug",
                           "proc_tot_sug",
                           "Filled_teeth", "Missing_teeth", "post_glucose_salivary_pH",
                           "physical_activity"),
           conditions  = c("age_s"))


p_cpo <- ord_sug %>%
  ord_plot(
    colour = "country_born.factor",
    size   = 3.2,
    alpha  = 0.6,
    auto_caption = 7,
    # shorter arrows look cleaner; adjust as you like
    constraint_vec_length = 4.2,
    constraint_vec_style  = vec_constraint(1.2, colour = "grey25"),
    constraint_lab_style  = constraint_lab_style(
      max_angle = 80, size = 3.2, aspect_ratio = 0.9, colour = "grey10"
    )
  ) +
  # thin ellipses by group
  stat_ellipse(aes(colour = country_born.factor), linewidth = 0.4, linetype = 1) +
  scale_color_brewer(palette = "Set1", name = "Country of birth") +
  # zero line for reference
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.2, colour = "grey50") +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.2, colour = "grey50") +
  # aspect ratio & padding
  coord_fixed(ratio = 0.9, clip = "off") +
  # cleaner theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.margin = margin(2, 2, 2, 2),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(
    x = NULL, y = NULL,
    title = NULL,
    subtitle = NULL
  )

print(p_cpo)


# 4) Save
ggsave("ord_sug.png", p_cpa, width = 8.5, height = 6.5, dpi = 300, bg = "white")


#==================================================
##.       Correlation plot
# Using spearman's correlation
# comparing atop 20 taxa at Genus rank
#==================================================
library(microViz)
library(dplyr)

# sugars
vars <- c(
  "fruit_total_sugars_g_day","fruit_free_sugars_g_day",
  "fruit_glucose_g_day","fruit_fructose_g_day","fruit_sucrose_g_day",
  "proc_glucose_g_day","proc_fructose_g_day","proc_sucrose_g_day",
  "proc_total_sugars_g_day","proc_added_sugars_g_day"
)

# IMPORTANT: tax_model LISTing of variable names
vars_list <- as.list(vars)

# Spearman correlation test that works with formula + data
cor_test_spearman <- function(formula, data, ...) {
  mf <- model.frame(formula, data)
  stats::cor.test(mf[[1]], mf[[2]], method = "spearman", exact = FALSE)
}

# top taxa from your raw-count object (non-negative)
top_taxa <- tax_top(ps_counts, 20, by = max, rank = "Genus")

correlations_df <- ps_counts %>%
  tax_model(
    trans = "clr",
    rank = "Genus",
    taxa = top_taxa,
    variables = vars_list,
    type = cor_test_spearman,
    return_psx = FALSE,
    verbose = FALSE
  ) %>%
  tax_models2stats(rank = "Genus") %>%
  mutate(p.FDR = p.adjust(p.value, method = "fdr"))


taxa_hclust <- correlations_df %>%
  select(term, taxon, estimate) %>%
  pivot_wider(id_cols = taxon, names_from = term, values_from = estimate) %>%
  column_to_rownames("taxon") %>%
  as.matrix() %>%
  dist(method = "euclidean") %>%
  hclust(method = "ward.D2")

taxa_order <- taxa_hclust$labels[taxa_hclust$order]

cor_plot <- correlations_df %>%
  mutate(
    taxon = factor(taxon, levels = taxa_order),
    term  = factor(term, levels = vars)   # keep your sugar order
  ) %>%
  ggplot(aes(x = term, y = taxon)) +
  geom_raster(aes(fill = estimate)) +
  geom_point(
    data = function(x) dplyr::filter(x, p.value < 0.05),
    shape = 8, size = 3, colour = "black"   # asterisk-like
  ) +
  geom_point(
    data = function(x) dplyr::filter(x, p.FDR < 0.05),
    shape = 16, size = 3, colour = "black"  # filled circle
  ) +
  scale_y_discrete(limits = taxa_order) +
  scale_fill_distiller(
    palette = "BrBG",
    limits  = c(-0.45, 0.45),
    oob     = scales::squish
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL, y = NULL,
    fill = "Spearman\ncorrelation",
    caption = paste(
      "Asterisk indicates p < 0.05 (not FDR adjusted)",
      "Filled circle indicates FDR p < 0.05",
      sep = "\n"
    )
  )

cor_plot

# ========================= ========================= =========================
# ASV-level PERMANOVA sensitivity analysis (FILTER rare ASVs) + tidy tables + CSV export
# - Filter ASVs: prevalence >= 10% samples AND total reads >= 100 (edit as needed)
# - Distance: Aitchison (Euclidean on CLR)
# - PERMANOVA: vegan::adonis2, by = "terms" (sequential; exposure tested after confounders)
# - Output:
#     1) tidy_all_asv   : covariates + exposure rows for each model
#     2) tidy_expo_only : exposure-only summary (plus BH-FDR across exposures)
# ========================= ========================= =========================

library(phyloseq)
library(vegan)

# ---- inputs ----
covars <- c("age", "gender.factor", "toothbrsh_freq")

sugar_vars <- c(
  "fruit_total_sugars_g_day", "fruit_free_sugars_g_day",
  "fruit_glucose_g_day", "fruit_fructose_g_day", "fruit_sucrose_g_day",
  "proc_glucose_g_day", "proc_fructose_g_day", "proc_sucrose_g_day",
  "proc_total_sugars_g_day", "proc_added_sugars_g_day", "proc_free_sugars_g_day"
)

# ---- rare-ASV filter thresholds ----
min_prev_prop   <- 0.10   # >=10% of samples (93 samples => >=10 samples)
min_total_reads <- 100    # total reads across all samples per ASV

# ---- 0) Start from ASV-level phyloseq object ----
ps_asv <- ps3

# ---- 1) Filter rare ASVs ----
X0 <- as(otu_table(ps_asv), "matrix")
if (!taxa_are_rows(ps_asv)) X0 <- t(X0)  # ensure taxa x samples

n_samp <- nsamples(ps_asv)
prev_counts <- rowSums(X0 > 0)
prev_prop   <- prev_counts / n_samp
tot_reads   <- rowSums(X0)

keep_taxa <- (prev_prop >= min_prev_prop) & (tot_reads >= min_total_reads)
ps_filt <- prune_taxa(keep_taxa, ps_asv)
ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)  # safety

cat("Original ASVs:", ntaxa(ps_asv), "\n")
cat("Kept ASVs:", ntaxa(ps_filt), "\n")
cat("Samples:", nsamples(ps_filt), "\n")

# ---- 2) Build samples x taxa matrix (filtered) ----
X <- as(otu_table(ps_filt), "matrix")
if (taxa_are_rows(ps_filt)) X <- t(X)    # now samples x taxa

# ---- 3) Aitchison distance = Euclidean on CLR ----
X <- X + 1
X <- X / rowSums(X)
clrX <- log(X) - rowMeans(log(X))
dist_ait_asv_filt <- dist(clrX, method = "euclidean")
Dmat <- as.matrix(dist_ait_asv_filt)

# ---- 4) Metadata as plain data.frame ----
meta <- data.frame(sample_data(ps_filt), check.names = FALSE)
meta$gender.factor  <- factor(meta$gender.factor)
meta$toothbrsh_freq <- factor(meta$toothbrsh_freq)

# ---- helper: run one model and return tidy rows (covars + exposure) ----
run_perm_tidy_asv <- function(expo, permutations = 999, seed = 1234) {
  
  # coerce exposure to numeric if needed
  if (is.factor(meta[[expo]]) || is.character(meta[[expo]])) {
    meta[[expo]] <<- as.numeric(as.character(meta[[expo]]))
  }
  
  # complete cases for covars + exposure
  keep <- complete.cases(meta[, c(covars, expo)])
  m2   <- meta[keep, , drop = FALSE]
  
  # subset distance to those samples
  D2 <- as.dist(Dmat[rownames(m2), rownames(m2)])
  
  # fit model: covars first, exposure last (sequential)
  fml <- as.formula(paste("D2 ~", paste(c(covars, expo), collapse = " + ")))
  set.seed(seed)
  ad <- adonis2(fml, data = m2, permutations = permutations, by = "terms")
  
  # tidy output: keep only covars + exposure rows
  tab <- as.data.frame(ad)
  tab$term <- rownames(tab)
  tab <- tab[tab$term %in% c(covars, expo),
             c("term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)"),
             drop = FALSE]
  
  tab$exposure <- expo
  tab$n        <- nrow(m2)
  tab$kept_ASVs <- ntaxa(ps_filt)
  
  tab[, c("exposure", "n", "kept_ASVs", "term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
}

# ---- 5) Run ALL sugar exposures + bind ----
tidy_all_asv <- do.call(rbind, lapply(sugar_vars, run_perm_tidy_asv))

# Exposure-only summary (usually what you report)
tidy_expo_only_asv <- tidy_all_asv[tidy_all_asv$term == tidy_all_asv$exposure, ]
tidy_expo_only_asv <- tidy_expo_only_asv[order(tidy_expo_only_asv$`Pr(>F)`), ]

# BH-FDR across exposures (exposure rows only)
tidy_expo_only_asv$p_adj_BH <- p.adjust(tidy_expo_only_asv$`Pr(>F)`, method = "BH")

# ---- 6) View tables ----
tidy_all_asv
tidy_expo_only_asv

# ---- 7) Save as CSV ----
write.csv(tidy_all_asv,
          "permanova_ASV_filtered_tidy_all_terms_by_terms.csv",
          row.names = FALSE)

write.csv(tidy_expo_only_asv,
          "permanova_ASV_filtered_tidy_exposures_only_by_terms.csv",
          row.names = FALSE)

# ========================= ========================= =========================
# CREATING BETA DIVERSITY PERMANOVA PLOTS WITH R2 VALUES
## ========================= ========================= =========================
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

df <- read_csv("permanova_ASV_filtered_tidy_exposures_only_by_terms.csv", show_col_types = FALSE)

# p-value column
if ("Pr(>F)" %in% names(df)) {
  df <- df %>% rename(p = `Pr(>F)`)
} else if (!"p" %in% names(df)) {
  stop("Couldn't find p-value column (expected 'Pr(>F)' or 'p').")
}

df_plot <- df %>%
  mutate(
    group = case_when(
      grepl("^fruit_", exposure) ~ "Fruit sugars",
      grepl("^proc_",  exposure) ~ "Processed sugars",
      TRUE ~ NA_character_
    ),
    label = gsub("^(fruit_|proc_)", "", exposure),
    star  = ifelse(!is.na(p) & p < 0.05, "*", "")
  ) %>%
  filter(!is.na(group)) %>%
  arrange(factor(group, levels = c("Fruit sugars","Processed sugars")), desc(R2)) %>%
  mutate(
    label2 = paste0(label, " (", ifelse(group=="Fruit sugars","fruit","proc"), ")"),
    x = factor(label2, levels = label2)
  )

# divider between groups
n_fruit <- sum(df_plot$group == "Fruit sugars")
divider_x <- n_fruit + 0.5

# ---- Create a gradient fill within each group based on R2 ----
# (light -> darker pastel)
fruit_ramp <- colour_ramp(c("#FBE3E5", "#E76F7B"))  # pastel red gradient
proc_ramp  <- colour_ramp(c("#E3F0FF", "#5B9DFF"))  # pastel blue gradient

df_plot <- df_plot %>%
  group_by(group) %>%
  mutate(r2_scaled = rescale(R2, to = c(0, 1), from = range(R2, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(
    fill_col = ifelse(group == "Fruit sugars",
                      fruit_ramp(r2_scaled),
                      proc_ramp(r2_scaled))
  )

p_beta <- ggplot(df_plot, aes(x = x, y = R2)) +
  geom_col(aes(fill = fill_col), width = 0.75, color = "grey35", linewidth = 0.2) +
  geom_text(aes(label = star), vjust = -0.35, size = 5) +
  geom_vline(xintercept = divider_x, linetype = "dashed", linewidth = 0.6, color = "grey35") +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = expression("PERMANOVA " * R^2 * " (ASV-level filtered; by = 'terms')"),
    title = expression("ASV-level PERMANOVA " * R^2 * " by sugar exposure"),
    subtitle = "* indicates p < 0.05; bar shade is a within-group R² gradient (fruit = red, processed = blue)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

p_beta

ggsave("permanova_ASV_filtered_R2_pastel_red_blue_gradient_divider.png", p,
       width = 13, height = 5, dpi = 300)

# ========================= ========================= =========================
#                   Differential Abundance Analysis
# Negative Binomial Model
# ========================= ========================= =========================
library(Maaslin2)
df_input_data = read.table(file = "table.txt", header = TRUE, sep = "\t",
                           row.names = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]
df_input_metadata = read.table(file= "metadata_filled.tsv", header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "fruit_total_sugar",
                    fixed_effects = c("fruit_total_sugars_g_day", "age","gender.factor"))


fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "fruit_free_sugar",
                    fixed_effects = c("fruit_free_sugars_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "fruit_glucose",
                    fixed_effects = c("fruit_glucose_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "fruit_fructose",
                    fixed_effects = c("fruit_fructose_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "fruit_sucrose",
                    fixed_effects = c("fruit_sucrose_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "proc_total_sugar",
                    fixed_effects = c("proc_total_sugars_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "proc_added_sugar",
                    fixed_effects = c("proc_added_sugars_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "proc_free_sugar",
                    fixed_effects = c("proc_free_sugars_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "proc_glucose",
                    fixed_effects = c("proc_glucose_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "proc_fructose",
                    fixed_effects = c("proc_fructose_g_day","age", "gender.factor"))

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "proc_sucrose",
                    fixed_effects = c("proc_sucrose_g_day","age", "gender.factor"))

# ========================= ========================= ========================= ##################### 
#=============== FOREST PLOT USING MAASLIN OUTPUT FOR FRUITS AND PROCESSED FOOODS ===================
# ========================= ========================= ========================= 

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# ---- 1) Read & detect columns ----
dat <- read_csv("fruitsugar_da.csv", show_col_types = FALSE)
nm  <- names(dat)

# heuristics to find needed columns
col_feature <- nm[str_detect(tolower(nm), "^feat|tax|genus|otu|asv")][1]
col_coef    <- nm[str_detect(tolower(nm), "^coef|beta|effect|estimate$")][1]
col_group   <- nm[str_detect(tolower(nm), "value|group|contrast|condition|topsugar|sample")][1]

stopifnot(!is.na(col_feature), !is.na(col_coef), !is.na(col_group))

dat <- dat %>%
  rename(
    feature = all_of(col_feature),
    coef    = all_of(col_coef),
    value   = all_of(col_group)
  )

# order taxa within groups (handle 'coef' name clash explicitly)
dat <- dat %>%
  mutate(feature = as.character(.data$feature)) %>%
  arrange(.data$value, .data$coef) %>%
  mutate(feature = factor(.data$feature, levels = unique(rev(.data$feature))))

# ---- 2) Palette adapting to observed groups ----
grps <- unique(dat$value) %>% as.character()
pal  <- c(
  "Total_sugar" = "#D55E00",
  "Glucose"     = "#0046FF",
  "Fructose"    = "purple",
  "Sucrose"     = "#41A67E"
)
pal_use <- setNames(ifelse(grps %in% names(pal), pal[grps], rep("#555555", length(grps))), grps)
shp_use <- setNames(seq(16, 16 + length(grps) - 1), grps)

# ---- 3) Plot (single panel, no facets, no CIs) ----
fruit_plot <- ggplot(dat, aes(x = coef, y = feature)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(aes(color = value, shape = value),
             size = 4.5,
             position = position_jitter(height = 0.10, width = 0)) +
  labs(x = "Coefficient", y = "Taxa (Genus)", color = NULL, shape = NULL) +
  scale_color_manual(values = pal_use) +
  scale_shape_manual(values = shp_use) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    legend.position = "right",
    plot.title = element_blank()
  )

print(fruit_plot)

#  save with white background
ggsave("fruitsugar_da_plot.png", da_plot, width = 9, height = 6, dpi = 300, bg = "white")
ggsave("fruitsugar_da_plot.svg", da_plot, width = 9, height = 6, bg = "white")


# ---- 1) Read & detect columns ----
dat <- read_csv("procsugar_da.csv", show_col_types = FALSE)
nm  <- names(dat)

# heuristics to find needed columns
col_feature <- nm[str_detect(tolower(nm), "^feat|tax|genus|otu|asv")][1]
col_coef    <- nm[str_detect(tolower(nm), "^coef|beta|effect|estimate$")][1]
col_group   <- nm[str_detect(tolower(nm), "value|group|contrast|condition|topsugar|sample")][1]

stopifnot(!is.na(col_feature), !is.na(col_coef), !is.na(col_group))

dat <- dat %>%
  rename(
    feature = all_of(col_feature),
    coef    = all_of(col_coef),
    value   = all_of(col_group)
  )

# order taxa within groups (handle 'coef' name clash explicitly)
dat <- dat %>%
  mutate(feature = as.character(.data$feature)) %>%
  arrange(.data$value, .data$coef) %>%
  mutate(feature = factor(.data$feature, levels = unique(rev(.data$feature))))

# ---- 2) Palette adapting to observed groups ----
grps <- unique(dat$value) %>% as.character()
pal  <- c(
  "Total_sugar" = "#D55E00",
  "Added_sugar" = "#432323",
  "Glucose"     = "#0046FF",
  "Fructose"    = "purple",
  "Sucrose"     = "#41A67E"
)
pal_use <- setNames(ifelse(grps %in% names(pal), pal[grps], rep("#555555", length(grps))), grps)
shp_use <- setNames(seq(16, 16 + length(grps) - 1), grps)

# ---- 3) Plot (single panel, no facets, no CIs) ----
proc_plot <- ggplot(dat, aes(x = coef, y = feature)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(aes(color = value, shape = value),
             size = 4.5,
             position = position_jitter(height = 0.10, width = 0)) +
  labs(x = "Coefficient", y = "Taxa (Genus)", color = NULL, shape = NULL) +
  scale_color_manual(values = pal_use) +
  scale_shape_manual(values = shp_use) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    legend.position = "right",
    plot.title = element_blank()
  )

print(proc_plot)


# save with white background
ggsave("fruitsugar_da_plot.png", da_plot, width = 9, height = 6, dpi = 300, bg = "white")
ggsave("fruitsugar_da_plot.svg", da_plot, width = 9, height = 6, bg = "white")


#=============================================
 # =============== FIGURE 1 ===============
#  =============================================
#A. p_cpo- Contrained partial analysis
  #B. p_beta- Beta diversity using eucleadean distance
  #C. cor_plot

library(cowplot)
library(ggplot2)

# ---- Right panel: two stacked plots (A, B) ----
left_panel <- plot_grid(
  p_cpo, p_beta,
  ncol = 1,
  rel_heights = c(1.6, 1),   
  labels = c("A", "B"),
  label_size = 16
)

# ---- Left panel: center cor_plot (C) vertically ----

right_panel <- plot_grid(
  NULL, cor_plot, NULL,
  ncol = 1,
  rel_heights = c(0.5, 3, 0.5),  
  labels = c("", "C", ""),
  label_size = 16
)

# ---- Final layout: left + right ----
final_panel <- plot_grid(
  left_panel, right_panel,
  ncol = 2,
  rel_widths = c(1, 1.2)  
)

print(final_panel)

# ---- Save ----
ggsave(
  "figure1_combined.svg",
  final_panel,
  width = 14, height = 10,
  device = "svg",
  bg = "white"
)

ggsave(
  "figure1_combined.png",
  final_panel,
  width = 14, height = 10,
  dpi = 300,
  bg = "white"
)

#=============================================
#  =============== FIGURE 2 ===============
#  =============================================
#A. Differential abundance analysis- fruits- fruit_plot
  #B. Differential abundance analysis- fruits - proc_plot
  #C. Linear regression alpha diversity - p_stacked

library(cowplot)
library(ggplot2)

# ---- Left panel: A bigger than B ----
left_panel2 <- plot_grid(
  fruit_plot, proc_plot,
  ncol = 1,
  rel_heights = c(1.6, 1),
  labels = c(
    "A. Differential abundance estimates for fruit-derived sugars",
    "B. Differential abundance estimates for processed food-derived sugars"
  ),
  label_size = 13,
  label_fontface = "bold",
  label_x = 0,   # left align
  label_y = 1,   # top
  hjust = 0,
  vjust = 1
)

# ---- Right panel: one plot spanning full height ----
right_panel2 <- plot_grid(
  p_stacked,
  ncol = 1,
  labels = c("C. Linear regression estimates for dietary sugars and alpha diversity"),
  label_size = 13,
  label_fontface = "bold",
  label_x = 0,
  label_y = 1,
  hjust = 0,
  vjust = 1
)

# ---- Final layout: left + right ----
figure2 <- plot_grid(
  left_panel2, right_panel2,
  ncol = 2,
  rel_widths = c(1, 1)
)

print(figure2)

# ---- Save ----
ggsave(
  "figure2_combined.svg",
  figure2,
  width = 14, height = 10,
  device = "svg",
  bg = "white"
)

# (You had final_panel2 here, which wasn't defined—using figure2 instead)
ggsave(
  "figure_two_panel.png",
  figure2,
  width = 14, height = 10,
  dpi = 300,
  bg = "white"
)
