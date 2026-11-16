# Builds data/participant_data.csv from the raw study files in data-raw/.
# The raw files are not part of the public repository; the output file is.
#
# 93 participants with plaque sequencing data
#  - 4 with no FFQ food-level record
#  - 1 implausible dietary report (participant 140: 17,395 kJ/day,
#      417 g/day animal protein)
#  = 88

source("R/00_setup.R")

meta   <- read.delim("data-raw/metadata_filled.tsv", check.names = FALSE,
                     stringsAsFactors = FALSE)
animal <- read.csv("data-raw/protein_animal_g.csv", check.names = FALSE)
plant  <- read.csv("data-raw/protein_plant_g.csv",  check.names = FALSE)

dat <- merge(meta, animal[, c("participant_id", "total_animal_protein")],
             by = "participant_id")
dat <- merge(dat, plant[, c("participant_id", "total_plant_protein")],
             by = "participant_id")
stopifnot(nrow(dat) == 89)

## Protein density (% of energy, 17 kJ/g) ------------------------------
dat$animal_pct_energy <- 100 * dat$total_animal_protein * 17 / dat$Energy..exc.fibre
dat$plant_pct_energy  <- 100 * dat$total_plant_protein  * 17 / dat$Energy..exc.fibre

## Protein source groups ------------------------------------------------
# Top tertile of animal or plant protein density; everyone else is the
# reference group. Cut points were set on the 89 participants with FFQ
# data, before participant 140 was excluded. Participants in the top
# tertile of both are assigned to the source with the higher percentile rank.
cut_animal <- quantile(dat$animal_pct_energy, 2/3)
cut_plant  <- quantile(dat$plant_pct_energy,  2/3)
prank <- function(x) sapply(x, function(v) mean(x <= v))

hi_animal <- dat$animal_pct_energy >= cut_animal
hi_plant  <- dat$plant_pct_energy  >= cut_plant
both      <- hi_animal & hi_plant

dat$protein_group <- "Reference"
dat$protein_group[hi_animal & !hi_plant] <- "High animal"
dat$protein_group[hi_plant & !hi_animal] <- "High plant"
dat$protein_group[both] <- ifelse(
  prank(dat$animal_pct_energy)[both] >= prank(dat$plant_pct_energy)[both],
  "High animal", "High plant")

cat(sprintf("Tertile cut points: animal %.2f %%E, plant %.2f %%E\n",
            cut_animal, cut_plant))

## Exclude the implausible dietary report -------------------------------
dat <- dat[dat$participant_id != 140, ]
stopifnot(nrow(dat) == 88)

## Analysis variables ---------------------------------------------------
tertiary <- c(
  "Bachelor or Honours degree level",
  "Postgraduate degree level-University Masters degree or PhD",
  "Advanced diploma/Diploma level/Associate degree-CAE; Teachers College or Nursing; University diploma",
  "Graduate diploma/Graduate certificate level- Graduate specialisation after Bachelor degree",
  "Certificate level 1, 2, 3, 4- Trade certificate/apprenticeship/vocational (eg TAFE, hairdressing)"
)

out <- data.frame(
  sample_id         = dat$sampleid,
  participant_id    = dat$participant_id,
  age               = dat$age,
  gender            = dat$gender.factor,
  country_of_birth  = ifelse(dat$country_born.factor == "Australia",
                             "Australia", "Other"),
  education         = ifelse(dat$edu_lvl.factor %in% tertiary,
                             "Tertiary", "Secondary or less"),
  employment        = ifelse(dat$emp_stat.factor == "Working in some capacity",
                             "Employed", "Not working"),
  physical_activity = c("Minimal activity" = "Low",
                        "Moderate activity" = "Moderate",
                        "High activity" = "High")[dat$physical_activity_status],
  met_min_week      = dat$total_met_minutes,
  alcohol           = ifelse(dat$alcohol_per_week.factor == "less than 1 drink",
                             "Occasional or non-drinker",
                             ">=1 standard drink/week"),
  # "Rather not say" (n = 1) is grouped with non-smokers
  smoking           = ifelse(dat$smoking.factor == "Yes",
                             "Current smoker", "Non-smoker"),
  toothbrushing     = ifelse(dat$toothbrsh_freq.factor == "Once",
                             "Once daily", "Twice daily or more"),
  flossing          = ifelse(dat$times_flossing.factor == "Once daily",
                             "Once daily or more", "Occasionally or never"),
  last_dental_visit = ifelse(dat$last_visit_dent.factor == "Less than 12 months",
                             "<12 months", ">=12 months"),
  filled_teeth      = dat$filled,
  missing_teeth     = dat$missing,
  saliva_ph_post_glucose = dat$post.glucose_2,
  saliva_ph_fall    = dat$ph_fall,
  energy_kj         = dat$Energy..exc.fibre,
  protein_g         = dat$Protein,
  fat_g             = dat$Fat,
  carbohydrate_g    = dat$Carbohydrate,
  animal_protein_g  = dat$total_animal_protein,
  plant_protein_g   = dat$total_plant_protein,
  animal_pct_energy = dat$animal_pct_energy,
  plant_pct_energy  = dat$plant_pct_energy,
  protein_group     = dat$protein_group,
  stringsAsFactors  = FALSE
)
out <- out[order(out$participant_id), ]
stopifnot(!anyNA(out))

print(table(out$protein_group))
write.csv(out, PATHS$participants, row.names = FALSE)
