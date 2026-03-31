######### TABLE 1 STUDY PARTICIPANT CHARACTERISTIC ##################################
########## K10 SCORES AND DEMOGRAPHIC VARIABES #############
####################################################################
################# WRITTEN BY SONIA NATH #################
################# DATE 15-FEB-2026 #################


#Clear existing data and graphics
rm(list=ls())
graphics.off()

#SETTING THE DIRECTORY
setwd("/Users/a1799090/K10_OMT/R_analysis")

# ---------------------------
# 0) Packages
# ---------------------------
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# ---------------------------
# 1) Read data
# ---------------------------
data <- read_tsv("metadata_filled.tsv", show_col_types = FALSE)

# ---------------------------
# 2) Create variables needed for the table
# ---------------------------
data <- data %>%
  mutate(
    # Keep mf (filled + missing) WITHOUT keeping filled/missing as rows in the table
    # If one of them is missing, treat it as 0; if both missing, mf = NA
    mf = dplyr::if_else(
      is.na(filled) & is.na(missing),
      NA_real_,
      dplyr::coalesce(filled, 0) + dplyr::coalesce(missing, 0)
    ),
    
    # K10 3-category
    k10_comb1 = case_when(
      is.na(k10score) ~ NA_character_,
      k10score >= 10 & k10score <= 15 ~ "Low",
      k10score >= 16 & k10score <= 21 ~ "Moderate",
      k10score >= 22 & k10score <= 50 ~ "High/Very high",
      TRUE ~ NA_character_
    ),
    k10_comb1 = factor(k10_comb1,
                      levels = c("Low", "Moderate", "High/Very high"),
                      ordered = TRUE),
    
    # Age group
    age_group = case_when(
      is.na(age) ~ NA_character_,
      age >= 18 & age <= 25 ~ "18-25",
      age >= 26 & age <= 35 ~ "26-35",
      age > 35 ~ ">35",
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("18-25", "26-35", ">35")),
    
    # Education (3-level)
    edu_lvl3 = case_when(
      is.na(edu_lvl) ~ NA_character_,
      edu_lvl %in% 1:5 ~ "Certificate / School Level",
      edu_lvl %in% 6:7 ~ "Undergraduate / Diploma Level",
      edu_lvl %in% 8:9 ~ "Postgraduate / Graduate Level",
      TRUE ~ NA_character_
    ),
    edu_lvl3 = factor(
      edu_lvl3,
      levels = c("Certificate / School Level",
                 "Undergraduate / Diploma Level",
                 "Postgraduate / Graduate Level")
    ),
    
    # Employment
    employment = recode(emp_stat.factor,
                        "Working in some capacity" = "Employed",
                        "Not currently working"    = "Not working",
                        .default = NA_character_),
    employment = factor(employment, levels = c("Employed", "Not working")),
    
    # Physical activity
    physical_activity = case_when(
      is.na(total_met_minutes) ~ NA_character_,
      total_met_minutes < 600 ~ "Minimal activity (<600 met min)",
      total_met_minutes >= 600 & total_met_minutes <= 1500 ~ "Moderate activity (600-1500 met min)",
      total_met_minutes > 1500 ~ "High activity (> 1500 met min)",
      TRUE ~ NA_character_
    ),
    physical_activity = factor(
      physical_activity,
      levels = c("Minimal activity (<600 met min)",
                 "Moderate activity (600-1500 met min)",
                 "High activity (> 1500 met min)")
    ),
    
    # Alcohol
    alcohol_cat = case_when(
      alcohol_per_week.factor == "less than 1 drink" ~ "Occasional/non-alcohol consumers",
      alcohol_per_week.factor %in% c("1-10 Alcoholic drink", "11-19 Alcoholic drinks") ~ ">= 1 standard drink/week",
      TRUE ~ NA_character_
    ),
    alcohol_cat = factor(alcohol_cat,
                         levels = c(">= 1 standard drink/week", "Occasional/non-alcohol consumers")),
    
    # Smoking
    smoking_cat = case_when(
      smoking.factor == "Yes" ~ "Smoker",
      smoking.factor == "No"  ~ "Non-smoker",
      TRUE ~ NA_character_
    ),
    smoking_cat = factor(smoking_cat, levels = c("Smoker", "Non-smoker")),
    
    # Toothbrushing
    toothbrushing = case_when(
      toothbrsh_freq.factor == "Once"  ~ "Once",
      toothbrsh_freq.factor == "Twice" ~ "Twice or more",
      TRUE ~ NA_character_
    ),
    toothbrushing = factor(toothbrushing, levels = c("Once", "Twice or more")),
    
    # Flossing
    flossing = case_when(
      times_flossing.factor == "Once daily" ~ "Once or more everyday",
      times_flossing.factor == "Occasionally" ~ "Occasionally/never",
      TRUE ~ NA_character_
    ),
    flossing = factor(flossing, levels = c("Once or more everyday", "Occasionally/never")),
    
    # Last dental visit
    last_dental_visit = case_when(
      last_visit_dent.factor == "Less than 12 months" ~ "Less than 12 months",
      !is.na(last_visit_dent.factor) ~ "More than 12 months",
      TRUE ~ NA_character_
    ),
    last_dental_visit = factor(last_dental_visit,
                               levels = c("Less than 12 months", "More than 12 months"))
  )

# ---------------------------
# 3) Helper functions
# ---------------------------
fmt_mean_sd <- function(x, digits = 2) {
  if (all(is.na(x))) return(NA_character_)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), m, s)
}

fmt_n_pct_ci <- function(x, n, digits = 1) {
  if (is.na(n) || n == 0) return(NA_character_)
  pct <- 100 * x / n
  ci  <- binom.test(x, n)$conf.int
  sprintf(paste0("%d (%.", digits, "f%%; %.1f–%.1f)"),
          x, pct, 100 * ci[1], 100 * ci[2])
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  format.pval(p, digits = 3, eps = 0.001)
}

# overall association var vs k10 (chi-square or fisher)
p_assoc_overall <- function(df, var, outcome = "k10_comb1") {
  tmp <- df %>%
    select(all_of(c(var, outcome))) %>%
    filter(!is.na(.data[[var]]), !is.na(.data[[outcome]]))
  
  tab <- table(tmp[[var]], tmp[[outcome]])
  if (min(dim(tab)) < 2) return(NA_real_)
  
  chi <- suppressWarnings(chisq.test(tab))
  if (any(chi$expected < 5)) suppressWarnings(fisher.test(tab)$p.value) else chi$p.value
}

# pairwise p-values for categorical var across k10 groups (BH adjusted)
pairwise_p_cat <- function(df, var, outcome = "k10_comb1") {
  tmp <- df %>%
    select(all_of(c(var, outcome))) %>%
    filter(!is.na(.data[[var]]), !is.na(.data[[outcome]]))
  
  lev <- levels(tmp[[outcome]])
  prs <- combn(lev, 2, simplify = FALSE)
  
  raw <- map_dbl(prs, function(pp) {
    d2 <- tmp %>% filter(.data[[outcome]] %in% pp)
    tab <- table(d2[[var]], d2[[outcome]])
    if (min(dim(tab)) < 2) return(NA_real_)
    chi <- suppressWarnings(chisq.test(tab))
    if (any(chi$expected < 5)) suppressWarnings(fisher.test(tab)$p.value) else chi$p.value
  })
  
  adj <- p.adjust(raw, method = "BH")
  names(adj) <- map_chr(prs, ~ paste(.x[1], "vs", .x[2]))
  adj
}

# pairwise p-values for continuous var across k10 groups (Wilcoxon; BH adjusted)
pairwise_p_cont <- function(df, var, outcome = "k10_comb1") {
  tmp <- df %>%
    select(all_of(c(var, outcome))) %>%
    filter(!is.na(.data[[var]]), !is.na(.data[[outcome]]))
  
  lev <- levels(tmp[[outcome]])
  prs <- combn(lev, 2, simplify = FALSE)
  
  raw <- map_dbl(prs, function(pp) {
    d2 <- tmp %>% filter(.data[[outcome]] %in% pp)
    tryCatch(
      wilcox.test(d2[[var]] ~ d2[[outcome]], exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  })
  
  adj <- p.adjust(raw, method = "BH")
  names(adj) <- map_chr(prs, ~ paste(.x[1], "vs", .x[2]))
  adj
}

pairwise_to_string <- function(pvec) {
  if (all(is.na(pvec))) return(NA_character_)
  paste(paste0(names(pvec), "=", vapply(pvec, fmt_p, character(1))), collapse = "; ")
}

# ---------------------------
# 4) Categorical blocks (Overall column + 95% CI + overall & pairwise p-rows)
# ---------------------------
make_block <- function(df, var, label, outcome = "k10_comb1",
                       k10_levels = levels(df[[outcome]]),
                       cont1 = "k10score") {
  
  tmp <- df %>%
    select(all_of(c(var, outcome, cont1))) %>%
    filter(!is.na(.data[[var]]), !is.na(.data[[outcome]]))
  
  totalN <- nrow(tmp)
  
  counts <- tmp %>%
    count(level = .data[[var]], !!sym(outcome), name = "n") %>%
    group_by(level) %>%
    mutate(N = sum(n)) %>%
    ungroup() %>%
    complete(level, !!sym(outcome) := k10_levels, fill = list(n = 0)) %>%
    group_by(level) %>%
    mutate(N = max(N)) %>%
    ungroup() %>%
    mutate(cell = map2_chr(n, N, fmt_n_pct_ci)) %>%
    select(level, N, !!sym(outcome), cell) %>%
    pivot_wider(names_from = !!sym(outcome), values_from = cell)
  
  # Overall column: subgroup size out of total non-missing sample (with 95% CI)
  overall_col <- tmp %>%
    count(level = .data[[var]], name = "N_level") %>%
    mutate(Overall = map2_chr(N_level, totalN, fmt_n_pct_ci)) %>%
    select(level, Overall)
  
  conts <- tmp %>%
    group_by(level = .data[[var]]) %>%
    summarise(`K10 score` = fmt_mean_sd(.data[[cont1]]), .groups = "drop")
  
  p_overall <- p_assoc_overall(df, var, outcome)
  p_pair <- pairwise_p_cat(df, var, outcome)
  
  out <- counts %>%
    left_join(overall_col, by = "level") %>%
    left_join(conts, by = "level") %>%
    mutate(
      Category = label,
      Subgroups = as.character(level),
      p_value = NA_character_
    ) %>%
    select(Category, Subgroups, N, `K10 score`, Overall, all_of(k10_levels), p_value)
  
  na_k10 <- as.list(setNames(rep(NA_character_, length(k10_levels)), k10_levels))
  
  p_row_overall <- as_tibble(c(
    list(Category = label, Subgroups = "p-value (overall)", N = NA_integer_,
         `K10 score` = NA_character_, Overall = NA_character_),
    na_k10,
    list(p_value = fmt_p(p_overall))
  ))
  
  p_row_pair <- as_tibble(c(
    list(Category = label, Subgroups = "p-value (pairwise, BH)", N = NA_integer_,
         `K10 score` = NA_character_, Overall = NA_character_),
    na_k10,
    list(p_value = pairwise_to_string(p_pair))
  ))
  
  bind_rows(out, p_row_overall, p_row_pair)
}

data <- data %>%
  mutate(
    # ensure these are numeric for mean/sd + mf
    filled  = suppressWarnings(as.numeric(as.character(filled))),
    missing = suppressWarnings(as.numeric(as.character(missing))),
    
    # mf = filled + missing (if both missing -> NA)
    mf = dplyr::if_else(
      is.na(filled) & is.na(missing),
      NA_real_,
      dplyr::coalesce(filled, 0) + dplyr::coalesce(missing, 0)
    ),

  )
# ---------------------------
# 5) Continuous-variable rows (Overall + group means + overall & pairwise tests)
# ---------------------------
make_cont_row <- function(df, var, label,
                          outcome = "k10_comb1",
                          k10_levels = levels(df[[outcome]]),
                          k10_var = "k10score") {
  
  df2 <- df %>% filter(!is.na(.data[[outcome]]), !is.na(.data[[var]]))
  
  # group means
  by_k10 <- df2 %>%
    group_by(.data[[outcome]]) %>%
    summarise(cell = fmt_mean_sd(.data[[var]]), .groups = "drop") %>%
    rename(k10 = .data[[outcome]])
  
  cells <- setNames(rep(NA_character_, length(k10_levels)), k10_levels)
  cells[as.character(by_k10$k10)] <- by_k10$cell
  
  # overall
  overall_var <- fmt_mean_sd(df2[[var]])
  k10_overall <- fmt_mean_sd(df2[[k10_var]])
  
  # overall + pairwise tests
  p_overall <- tryCatch(kruskal.test(df2[[var]] ~ df2[[outcome]])$p.value,
                        error = function(e) NA_real_)
  p_pair <- pairwise_p_cont(df2, var, outcome)
  
  main_row <- tibble(
    Category = label,
    Subgroups = "",
    N = nrow(df2),
    `K10 score` = k10_overall,
    Overall = overall_var,
    Low = cells["Low"],
    Moderate = cells["Moderate"],
    `High/Very high` = cells["High/Very high"],
    p_value = fmt_p(p_overall)
  )
  
  pair_row <- tibble(
    Category = label,
    Subgroups = "p-value (pairwise, BH)",
    N = NA_integer_,
    `K10 score` = NA_character_,
    Overall = NA_character_,
    Low = NA_character_,
    Moderate = NA_character_,
    `High/Very high` = NA_character_,
    p_value = pairwise_to_string(p_pair)
  )
  
  bind_rows(main_row, pair_row)
}

# ---------------------------
# 6) Build the table
# ---------------------------

# Continuous rows
cont_blocks <- list(
  make_cont_row(data, "saliva_flow_rate",  "Saliva flow rate"),
  make_cont_row(data, "baseline_pH",       "Baseline pH"),
  make_cont_row(data, "Energy..exc.fibre", "Energy (exc fibre)"),
  make_cont_row(data, "Water",             "Water"),
  make_cont_row(data, "Carbohydrate",      "Carbohydrate"),
  make_cont_row(data, "Sugars",            "Sugars"),
  make_cont_row(data, "Protein",           "Protein"),
  make_cont_row(data, "Fat",               "Fat"),
  make_cont_row(data, "Dietary.Fibre",     "Dietary fibre"),
  make_cont_row(data, "filled",            "Filled"),
  make_cont_row(data, "missing",           "Missing"),
  make_cont_row(data, "post.glucose_2",    "Post glucose (2h)")
)

# Categorical blocks (with CI + overall and pairwise)
cat_blocks <- list(
  make_block(data, "age_group",           "Age group"),
  make_block(data, "gender.factor",       "Gender"),
  make_block(data, "country_born.factor", "Country of birth"),
  make_block(data, "edu_lvl3",            "Level of education"),
  make_block(data, "employment",          "Employment"),
  make_block(data, "physical_activity",   "Physical activity"),
  make_block(data, "alcohol_cat",         "Alcohol consumption"),
  make_block(data, "smoking_cat",         "Smoking"),
  make_block(data, "toothbrushing",       "Tooth brushing frequency"),
  make_block(data, "flossing",            "Flossing frequency"),
  make_block(data, "last_dental_visit",   "Last dental visit")
)

descriptive_tbl <- bind_rows(
  bind_rows(cont_blocks),
  bind_rows(cat_blocks)
)

# ---------------------------
# 7) Save output
# ---------------------------
write.csv(descriptive_tbl, "descriptive_table_by_k10_comb1new.csv", row.names = FALSE)

descriptive_tbl


readr::write_tsv(data, "metadata_filled.tsv")
