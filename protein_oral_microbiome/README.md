# Dietary protein source and the supragingival plaque microbiome

Analysis code for:

> Nath S, Jamieson L, Zilm P, Weyrich LS. Dietary protein sources modulate healthy Australian adult oral microbiomes. *Frontiers in Cellular and Infection Microbiology* (submitted).

Cross-sectional analysis of 88 healthy, caries-free Australian adults linking animal and plant protein intake (food frequency questionnaire) to the supragingival plaque microbiome (16S rRNA gene, V4).

## Repository structure

```
├── data/
│   ├── participant_data.csv        analysis dataset (n = 88)
│   ├── data_dictionary.csv         variable descriptions
│   └── qiime2/                     place table.qza, rooted-tree.qza, taxonomy.qza here
├── R/
│   ├── 00_setup.R                  paths, constants, helpers
│   ├── 01_prepare_participant_data.R   builds data/participant_data.csv from raw files
│   ├── analysis/
│   │   ├── 02_build_phyloseq.R
│   │   ├── 03_alpha_diversity.R
│   │   ├── 04_community_composition.R
│   │   ├── 05_differential_abundance.R
│   │   └── 06_cooccurrence_networks.R
│   ├── tables/                     Table 1, Table 2, Supplementary Tables 1-6
│   └── figures/                    Figures 1-2, Supplementary Figures 1-2
└── run_all.R
```

Analysis scripts write intermediate results to `results/`; table and figure scripts read them and write to `output/tables/` and `output/figures/`.

| Output | Script |
|---|---|
| Table 1 | `R/tables/table1_participant_characteristics.R` |
| Table 2 | `R/tables/table2_permanova.R` |
| Supplementary Table 1 | `R/tables/supp_table1_maaslin3.R` |
| Supplementary Table 2 | `R/tables/supp_table2_network_centrality.R` |
| Supplementary Tables 3-6 | `R/tables/supp_tables3to6_alpha_diversity.R` |
| Figure 1 | `R/figures/figure1_community_composition.R` |
| Figure 2 | `R/figures/figure2_microbial_features.R` |
| Supplementary Figure 1 | `R/figures/supp_figure1_amdr.R` |
| Supplementary Figure 2 | `R/figures/supp_figure2_alpha_groups.R` |

## Data

- **Sequence data:** NCBI BioProject [PRJNA1311717](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1311717). Processing in QIIME 2 (2024.11; DADA2, SILVA 138) is described in Supplementary Text S2.
- **Participant data:** `data/participant_data.csv` holds only the variables used in the analyses, with categorical variables collapsed to the levels reported. See `data/data_dictionary.csv`.

Analytic sample: 93 participants with plaque sequencing data, minus 4 without FFQ food-level records and 1 implausible dietary report = 88. Protein source groups are the top tertile of animal or plant protein density (% energy); cut points were set before the implausible report was excluded (see `R/01_prepare_participant_data.R`).


## Software

R 4.4.2. CRAN: vegan, picante, igraph, psych, dplyr, tidyr, tibble, ggplot2, ggrepel, ggpubr, rstatix, cowplot, patchwork, scales. Bioconductor: phyloseq, maaslin3. GitHub: [qiime2R](https://github.com/jbisanz/qiime2R), [microViz](https://github.com/david-barnett/microViz).

```r
install.packages(c("vegan", "picante", "igraph", "psych", "dplyr", "tidyr",
                   "tibble", "ggplot2", "ggrepel", "ggpubr", "rstatix",
                   "cowplot", "patchwork", "scales", "remotes", "BiocManager"))
BiocManager::install(c("phyloseq", "maaslin3"))
remotes::install_github("jbisanz/qiime2R")
install.packages("microViz", repos = c("https://david-barnett.r-universe.dev",
                                       "https://cloud.r-project.org"))
```

## Contact

Sonia Nath, Adelaide Dental School, The University of Adelaide.
