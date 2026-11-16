# Table 1. Participant characteristics and daily energy and protein
# intake by subgroup (mean ± SD). Welch t-test for two groups, one-way
# ANOVA for three.

source("R/00_setup.R")

dat <- read_participants()

dat$age_group <- cut(dat$age, c(-Inf, 25, 35, Inf),
                     labels = c("18–25", "26–35", ">35"))
dat$ph_group  <- ifelse(dat$saliva_ph_post_glucose < 5.5, "<5.5", "≥5.5")
dat$filled_group  <- ifelse(dat$filled_teeth == 0, "None",
                            ifelse(dat$filled_teeth <= 2, "1–2", "≥3"))
dat$missing_group <- ifelse(dat$missing_teeth == 0, "None", "≥1")

blocks <- list(
  list("age_group",         "Age (years)",          c("18–25", "26–35", ">35")),
  list("gender",            "Gender",               c("Male", "Female")),
  list("country_of_birth",  "Country of birth",     c("Australia", "Other")),
  list("education",         "Highest education",    c("Tertiary", "Secondary or less")),
  list("employment",        "Employment",           c("Employed", "Not working")),
  list("physical_activity", "Physical activity (MET-min/week)",
       c(Low = "Low (<600)", Moderate = "Moderate (600–1500)",
         High = "High (>1500)")),
  list("alcohol",           "Alcohol consumption",
       c(">=1 standard drink/week" = "≥1 standard drink/week",
         "Occasional or non-drinker" = "Occasional or non-drinker")),
  list("smoking",           "Smoking",              c("Current smoker", "Non-smoker")),
  list("toothbrushing",     "Toothbrushing frequency",
       c("Once daily", "Twice daily or more")),
  list("flossing",          "Flossing frequency",
       c("Once daily or more", "Occasionally or never")),
  list("last_dental_visit", "Last dental visit",
       c("<12 months" = "<12 months", ">=12 months" = "≥12 months")),
  list("ph_group",          "Salivary pH after glucose challenge",
       c("<5.5", "≥5.5")),
  list("filled_group",      "Filled teeth",         c("None", "1–2", "≥3")),
  list("missing_group",     "Missing teeth",        c("None", "≥1"))
)

outcomes <- c("Energy (kJ/day)"        = "energy_kj",
              "Total protein (g/day)"  = "protein_g",
              "Animal protein (g/day)" = "animal_protein_g",
              "Plant protein (g/day)"  = "plant_protein_g")

summarise_block <- function(var, label, levels) {
  if (is.null(names(levels))) names(levels) <- levels
  g <- factor(dat[[var]], names(levels), labels = levels)
  stopifnot(!anyNA(g))
  n <- table(g)

  cells <- sapply(outcomes, function(y)
    tapply(dat[[y]], g, function(x) sprintf("%.1f ± %.1f", mean(x), sd(x))))
  p <- sapply(outcomes, function(y) {
    if (nlevels(g) == 2) t.test(dat[[y]] ~ g)$p.value
    else summary(aov(dat[[y]] ~ g))[[1]][["Pr(>F)"]][1]
  })

  body <- data.frame(Characteristic = c(label, rep("", nlevels(g) - 1)),
                     Category = levels(g),
                     `n (%)` = sprintf("%d (%.1f)", n, 100 * n / nrow(dat)),
                     matrix(cells, nrow = nlevels(g),
                            dimnames = list(NULL, names(outcomes))),
                     check.names = FALSE)
  prow <- body[1, ]
  prow[] <- ""
  prow$Category <- "p value"
  prow[names(outcomes)] <- ifelse(p < 0.001, "<0.001", sprintf("%.2f", p))
  rbind(body, prow)
}

table1 <- do.call(rbind, lapply(blocks, \(b) summarise_block(b[[1]], b[[2]], b[[3]])))
write_table(table1, "table1_participant_characteristics")

## Whole-sample values quoted in the Results -----------------------------
ms <- function(x) sprintf("%.1f (SD %.1f)", mean(x), sd(x))
cat("n =", nrow(dat), "\n")
for (v in outcomes) cat(sprintf("  %-17s %s\n", v, ms(dat[[v]])))
cat(sprintf("  animal vs plant density, r = %.2f\n",
            cor(dat$animal_pct_energy, dat$plant_pct_energy)))
