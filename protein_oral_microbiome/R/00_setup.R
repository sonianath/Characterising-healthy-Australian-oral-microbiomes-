# Shared paths, constants and helpers. Sourced by every other script.
# Run all scripts from the repository root.

invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))

PATHS <- list(
  participants = "data/participant_data.csv",
  qiime2       = "data/qiime2",
  results      = "results",
  tables       = "output/tables",
  figures      = "output/figures"
)
for (d in PATHS[c("results", "tables", "figures")]) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

EXPOSURES <- c("animal_pct_energy", "plant_pct_energy",
               "animal_protein_g", "plant_protein_g")
EXPO_LAB  <- c(animal_pct_energy = "Animal protein (%E)",
               plant_pct_energy  = "Plant protein (%E)",
               animal_protein_g  = "Animal protein (g/day)",
               plant_protein_g   = "Plant protein (g/day)")
PRIMARY   <- c("animal_pct_energy", "plant_pct_energy")

# PERMANOVA covariates, entered sequentially before the exposure
PERMANOVA_COVARS <- c("age", "gender", "toothbrushing")

METRICS <- c("Observed", "Chao1", "FaithPD", "Shannon", "Simpson")
M_LAB   <- c(Observed = "Observed ASVs", Chao1 = "Chao1",
             FaithPD = "Faith's PD", Shannon = "Shannon", Simpson = "Simpson")

GROUP_LEVELS <- c("Reference", "High animal", "High plant")
PAL     <- c("Reference" = "#8C8C8C", "High animal" = "#C0504D",
             "High plant" = "#4F8A5B")
SRC_PAL <- c("Animal protein" = "#C0504D", "Plant protein" = "#4F8A5B")

RAREFY_DEPTH    <- 8000
MIN_PREVALENCE  <- 0.10   # ASV filter for compositional analyses
MIN_TOTAL_READS <- 100

read_participants <- function() {
  d <- read.csv(PATHS$participants, stringsAsFactors = FALSE)
  d$protein_group <- factor(d$protein_group, GROUP_LEVELS)
  d$gender        <- factor(d$gender, c("Male", "Female"))
  d$toothbrushing <- factor(d$toothbrushing)
  rownames(d) <- d$sample_id
  d
}

clr_transform <- function(counts) {      # samples x taxa
  x <- counts + 1
  x <- x / rowSums(x)
  log(x) - rowMeans(log(x))
}

fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001",
                              formatC(p, format = "f", digits = digits)))
}
fmt_num <- function(x, digits = 2) {
  x <- round(x, digits)
  x[x == 0] <- 0
  sub("^-", "−", formatC(x, format = "f", digits = digits))
}

write_table <- function(x, name) {
  write.csv(x, file.path(PATHS$tables, paste0(name, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  invisible(x)
}
