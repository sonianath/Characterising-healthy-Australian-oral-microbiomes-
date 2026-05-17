######### BIOFILM- EXCELLIGENCE ANALYSIS ##################################
########## OMT SALIVA EXCELLIGENCE DATA  #############
####################################################################
################# WRITTEN BY SONIA NATH #################
################# DATE 05-MAY-2026 #################

#Clear existing data and graphics
rm(list=ls())
graphics.off()

#SETTING THE DIRECTORY
setwd("/Users/a1799090/Excelligence/")

###################################################
# A) CLEANING INSPECTING AND FORMATTING THE DATA
#1. IMPORT AND INSPECT THE DATA
###################################################
library(tidyverse)

# File 1: maximum biofilm and time to maximum
biofilm_max <- read_csv("biofilm.csv", show_col_types = FALSE)

# File 2: full time-course data
biofilm_full <- read_csv("biofilm_full.csv", show_col_types = FALSE)

glimpse(biofilm_max)
glimpse(biofilm_full)

# Check missing values
colSums(is.na(biofilm))

###################################################
#2. CONVERT THE DATASET FROM WIDE TO LONG FORMAT
###################################################
sample_cols <- setdiff(names(biofilm_full), c("Time_Hour", "Mean_growth"))

biofilm_long <- biofilm_full %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample_number",
    values_to = "cell_index"
  ) %>%
  mutate(
    sample_number = as.integer(sample_number)
  ) %>%
  arrange(sample_number, Time_Hour)

head(biofilm_long)

##################################
##3. KINETIC CURVE ANALYSIS
#PLOT ALL BIOFILM GROWTH CURVES
###################################################
growth_curves<- ggplot(biofilm_long, aes(x = Time_Hour, y = cell_index, group = sample_number)) +
  geom_line(alpha = 0.4) +
  labs(
    title = "Biofilm formation curves over 12 hours",
    x = "Time, hours",
    y = "Cell index"
  ) +
  theme_minimal()
growth_curves

##################################
#4. PLOT THE AVERAGE BIOFILM CURVE
##################################
library(tidyverse)

biofilm_summary <- biofilm_long %>%
  group_by(Time_Hour) %>%
  summarise(
    mean_cell_index = mean(cell_index, na.rm = TRUE),
    sd_cell_index = sd(cell_index, na.rm = TRUE),
    n = sum(!is.na(cell_index)),
    se_cell_index = sd_cell_index / sqrt(n),
    lower_ci = mean_cell_index - 1.96 * se_cell_index,
    upper_ci = mean_cell_index + 1.96 * se_cell_index,
    .groups = "drop"
  )

biofilm_curve <- ggplot(biofilm_summary, aes(x = Time_Hour, y = mean_cell_index)) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    fill = "#A6CEE3",   # pastel pink
    alpha = 0.45
  ) +
  geom_line(
    colour = "#1F78B4",   # pastel mauve
    linewidth = 1.4
  ) +
  geom_point(
    colour = "black",
    fill = "#1F78B4",
    size = 2.8,
    shape = 21,
    stroke = 0.6
  ) +
  labs(
    subtitle = "Mean cell index with 95% confidence interval",
    x = "Time (hours)",
    y = "Mean cell index"
  ) +
  scale_x_continuous(
    breaks = seq(0, 12, by = 2),
    limits = c(0, 12),
    expand = c(0.01, 0.01)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      colour = "#5C4B51"
    ),
    plot.subtitle = element_text(
      size = 12.5,
      hjust = 0.5,
      colour = "#7A6A70",
      margin = margin(b = 12)
    ),
    axis.title = element_text(
      face = "bold",
      size = 13,
      colour = "#4F4F4F"
    ),
    axis.text = element_text(
      size = 11,
      colour = "#5A5A5A"
    ),
    axis.line = element_line(
      colour = "#8C8C8C",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      colour = "#8C8C8C",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.4
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    )
  )
biofilm_curve #figure 1A

#######################################################################################################################
#B. Derive new biofilm features from the full curve
#5. Calculate AUC, maximum growth rate, time to peak, and decline after peak
#The area under the curve, or AUC, is very useful because it measures the total biofilm burden over time, not only the maximum value.
#######################################################################################################################
# Trapezoidal AUC function
trapz_auc <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

# AUC within a time interval
auc_between <- function(t, y, lower, upper) {
  keep <- t >= lower & t <= upper
  
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  
  trapz_auc(t[keep], y[keep])
}

# First time that a sample crosses a threshold
first_crossing_time <- function(t, y, threshold) {
  idx <- which(y >= threshold)
  
  if (length(idx) == 0) {
    return(NA_real_)
  } else {
    return(t[min(idx)])
  }
}

biofilm_features <- biofilm_long %>%
  group_by(sample_number) %>%
  group_modify(~ {
    
    d <- .x %>%
      arrange(Time_Hour)
    
    t <- d$Time_Hour
    y <- d$cell_index
    
    max_y <- max(y, na.rm = TRUE)
    max_index <- which.max(y)
    
    slopes <- diff(y) / diff(t)
    slope_time <- (head(t, -1) + tail(t, -1)) / 2
    
    tibble(
      maxcellindex_from_full = max_y,
      time_to_max = t[max_index],
      
      auc_0_12 = trapz_auc(t, y),
      auc_0_3 = auc_between(t, y, 0, 3),
      auc_3_6 = auc_between(t, y, 3, 6),
      auc_6_12 = auc_between(t, y, 6, 12),
      
      max_growth_rate = max(slopes, na.rm = TRUE),
      time_of_max_growth_rate = slope_time[which.max(slopes)],
      
      time_to_half_max = first_crossing_time(t, y, 0.5 * max_y),
      time_to_10_percent_max = first_crossing_time(t, y, 0.1 * max_y),
      time_to_cellindex_0_05 = first_crossing_time(t, y, 0.05),
      
      end_cellindex = y[length(y)],
      decline_after_peak = max_y - y[length(y)],
      proportional_decline_after_peak = (max_y - y[length(y)]) / abs(max_y)
    )
  }) %>%
  ungroup()

biofilm_features

# Merging the files
biofilm_all_features <- biofilm_features %>%
  left_join(biofilm_max, by = "sample_number")

head(biofilm_all_features)

biofilm_all_features %>%
  summarise(
    max_difference_in_max_biofilm =
      max(abs(maxcellindex_from_full - maxcellindex), na.rm = TRUE),
    
    max_difference_in_time_to_max =
      max(abs(time_to_max - Time_max_biofilm), na.rm = TRUE)
  )

###################################################
#6.Descriptive statistics
###################################################
biofilm_all_features %>%
  summarise(
    n = n(),
    
    mean_max_biofilm = mean(maxcellindex_from_full, na.rm = TRUE),
    sd_max_biofilm = sd(maxcellindex_from_full, na.rm = TRUE),
    median_max_biofilm = median(maxcellindex_from_full, na.rm = TRUE),
    
    mean_time_to_max = mean(time_to_max, na.rm = TRUE),
    sd_time_to_max = sd(time_to_max, na.rm = TRUE),
    median_time_to_max = median(time_to_max, na.rm = TRUE),
    
    mean_auc = mean(auc_0_12, na.rm = TRUE),
    sd_auc = sd(auc_0_12, na.rm = TRUE),
    median_auc = median(auc_0_12, na.rm = TRUE),
    
    mean_max_growth_rate = mean(max_growth_rate, na.rm = TRUE),
    sd_max_growth_rate = sd(max_growth_rate, na.rm = TRUE)
  )

###################################################
#7. classifying samples by biofilm forming pattern
#Low, moderate, and high biofilm capacity: 12 IN EACH GROUP
####################################################################
biofilm_all_features <- biofilm_all_features %>%
  mutate(
    biofilm_capacity_group = ntile(maxcellindex_from_full, 3),
    biofilm_capacity_group = case_when(
      biofilm_capacity_group == 1 ~ "Low biofilm capacity",
      biofilm_capacity_group == 2 ~ "Moderate biofilm capacity",
      biofilm_capacity_group == 3 ~ "High biofilm capacity"
    )
  )

biofilm_all_features %>%
  count(biofilm_capacity_group) # 12 in each group

########################################################################
#8. Fast, medium, and slow biofilm forming samples
#This uses time to maximum biofilm.
########################################################################
biofilm_all_features <- biofilm_all_features %>%
  mutate(
    speed_group = ntile(time_to_max, 3),
    speed_group = case_when(
      speed_group == 1 ~ "Fast peak",
      speed_group == 2 ~ "Medium peak",
      speed_group == 3 ~ "Slow peak"
    )
  )

biofilm_all_features %>%
  count(speed_group) #12 in each group

################################################################################################
#9. Plot curves based on high, moderate, and low biofilm forming groups
################################################################################################
biofilm_long_grouped <- biofilm_long %>%
  left_join(
    biofilm_all_features %>%
      dplyr::select(sample_number, biofilm_capacity_group, speed_group),
    by = "sample_number"
  )

ggplot(
  biofilm_long_grouped,
  aes(x = Time_Hour, y = cell_index, group = sample_number)
) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ biofilm_capacity_group) +
  labs(
    x = "Time, hours",
    y = "Cell index"
  ) +
  theme_minimal()


library(dplyr)
library(ggplot2)
library(ggpubr)

# Arrange the groups in the correct order
biofilm_all_features <- biofilm_all_features %>%
  mutate(
    biofilm_capacity_group = factor(
      biofilm_capacity_group,
      levels = c(
        "Low biofilm capacity",
        "Moderate biofilm capacity",
        "High biofilm capacity"
      )
    )
  )

# Pairwise comparisons
my_comparisons <- list(
  c("Low biofilm capacity", "Moderate biofilm capacity"),
  c("Low biofilm capacity", "High biofilm capacity"),
  c("Moderate biofilm capacity", "High biofilm capacity")
)

max_biofilm <- ggplot(
  biofilm_all_features,
  aes(x = biofilm_capacity_group, y = time_to_max, fill = biofilm_capacity_group)
) +
  geom_boxplot(
    alpha = 0.8,
    width = 0.65,
    colour = "grey25",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.75,
    size = 2,
    colour = "black"
  ) +
  scale_fill_manual(values = c(
    "Low biofilm capacity" = "#F4C2C2",
    "Moderate biofilm capacity" = "#CDEAC0",
    "High biofilm capacity" = "#BFD7EA"
  )) +
  labs(
    subtitle = "Box plots with Kruskal–Wallis and pairwise Wilcoxon tests",
    x = "Biofilm capacity group",
    y = "Time to maximum biofilm (hours)"
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(biofilm_all_features$time_to_max, na.rm = TRUE) + 1.2
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "#F0F0F0", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

max_biofilm # INCLUDE IN SUPPLEMENTARY FIGURE 1-PLOTTING TIME AND MAXIMUM BIOFILM FORMING GROUPS

## ANOTHER PLOT
biofilm_long_grouped <- biofilm_long_grouped %>%
  mutate(
    biofilm_capacity_group = factor(
      biofilm_capacity_group,
      levels = c(
        "Low biofilm capacity",
        "Moderate biofilm capacity",
        "High biofilm capacity"
      )
    ),
    time_bin = cut(
      Time_Hour,
      breaks = seq(0, 12, by = 1),
      include.lowest = TRUE,
      right = FALSE
    )
  )


# Make sure the groups are ordered
biofilm_long_grouped <- biofilm_long_grouped %>%
  mutate(
    biofilm_capacity_group = factor(
      biofilm_capacity_group,
      levels = c(
        "Low biofilm capacity",
        "Moderate biofilm capacity",
        "High biofilm capacity"
      )
    )
  )

biofilm_groups<- ggplot(
  biofilm_long_grouped,
  aes(x = time_bin, y = cell_index, fill = biofilm_capacity_group)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.8,
    width = 0.65,
    colour = "grey25",
    outlier.shape = NA
  ) +
  geom_point(
    aes(group = biofilm_capacity_group),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.30,
    size = 1,
    colour = "black"
  ) +
  stat_compare_means(
    aes(group = biofilm_capacity_group),
    method = "kruskal.test",
    label = "p.signif",
    hide.ns = FALSE,
    size = 3
  ) +
  scale_fill_manual(values = c(
    "Low biofilm capacity" = "#F4C2C2",
    "Moderate biofilm capacity" = "#CDEAC0",
    "High biofilm capacity" = "#BFD7EA"
  )) +
  labs(
    x = "Time (hours)",
    y = "Cell index",
    fill = "Biofilm capacity group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )
biofilm_groups #FIGURE 1D

######################################################
# 10. Creating box plots based on speed groups.
######################################################
library(dplyr)
library(ggplot2)
library(ggpubr)

biofilm_all_features <- biofilm_all_features %>%
  mutate(
    speed_group = ntile(time_to_max, 3),
    speed_group = case_when(
      speed_group == 1 ~ "Fast peak",
      speed_group == 2 ~ "Medium peak",
      speed_group == 3 ~ "Slow peak"
    )
  ) 


# Arrange the groups in the order you want
biofilm_all_features <- biofilm_all_features %>%
  mutate(
    speed_group = factor(
      speed_group,
      levels = c("Fast peak", "Medium peak", "Slow peak")
    )
  )

# Pairwise comparisons
my_comparisons <- list(
  c("Fast peak", "Medium peak"),
  c("Fast peak", "Slow peak"),
  c("Medium peak", "Slow peak")
)

speed_groups <- ggplot(
  biofilm_all_features,
  aes(x = speed_group, y = maxcellindex_from_full, fill = speed_group)
) +
  geom_boxplot(
    alpha = 0.85,
    width = 0.62,
    colour = "grey25",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.75,
    size = 2.2,
    colour = "black"
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(biofilm_all_features$maxcellindex_from_full, na.rm = TRUE) + 0.06
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",   # change to "p.format" if you want exact p-values
    hide.ns = FALSE
  ) +
  scale_fill_manual(values = c(
    "Fast peak" = "#F4C2C2",
    "Medium peak" = "#CDEAC0",
    "Slow peak" = "#BFD7EA"
  )) +
  labs(
    x = "Speed group",
    y = "Maximum cell index"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      colour = "grey40"
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text = element_text(
      size = 11,
      colour = "grey20"
    ),
    legend.position = "none",
    panel.grid.major = element_line(
      colour = "#E6E6E6",
      linewidth = 0.45
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 20, 15, 15)
  )

speed_groups #SUPPLEMENTARY FIGURE 2

####################################
#11. AUC plots
####################################
library(dplyr)
biofilm_all_features %>%
  arrange(desc(auc_0_12)) %>%
  dplyr::select(
    sample_number,
    maxcellindex_from_full,
    time_to_max,
    auc_0_12,
    max_growth_rate,
    biofilm_capacity_group,
    speed_group
  )
biofilm_auc_ranked <- biofilm_all_features %>%
  arrange(desc(auc_0_12)) %>%
  dplyr::select(
    sample_number,
    maxcellindex_from_full,
    time_to_max,
    auc_0_12,
    max_growth_rate,
    biofilm_capacity_group,
    speed_group
  )

biofilm_auc_ranked

##Correlation scratter plots
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rlang)

make_biofilm_scatter <- function(data, xvar, yvar, title, subtitle,
                                 xlab, ylab,
                                 point_fill = "#BFD7EA",
                                 line_colour = "#5B8DB8",
                                 ribbon_fill = "#E8F1F8") {
  
  x_sym <- sym(xvar)
  y_sym <- sym(yvar)
  
  # Spearman correlation
  cor_res <- cor.test(
    data[[xvar]],
    data[[yvar]],
    method = "spearman"
  )
  
  cor_label <- paste0(
    "Spearman's rho = ", round(cor_res$estimate, 2),
    "\n",
    "p = ", format.pval(cor_res$p.value, digits = 3, eps = .001)
  )
  
  ggplot(data, aes(x = !!x_sym, y = !!y_sym)) +
    geom_point(
      shape = 21,
      size = 3.2,
      stroke = 0.5,
      colour = "black",
      fill = point_fill,
      alpha = 0.9
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      colour = line_colour,
      fill = ribbon_fill,
      linewidth = 1.3
    ) +
    stat_regline_equation(
      label.x.npc = "left",
      label.y.npc = 0.98,
      size = 4.2
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = cor_label,
      hjust = 1.1,
      vjust = 2.8,
      size = 4.3,
      colour = "#4F4F4F",
      fontface = "bold"
    ) +
    labs(
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, colour = "grey40"),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, colour = "grey20"),
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 20, 15, 15)
    )
}

#Maximum biofilm vs time to max biofilm
plot_time <- make_biofilm_scatter(
  data = biofilm_all_features,
  xvar = "time_to_max",
  yvar = "maxcellindex_from_full",
  subtitle = "Linear regression line with Spearman correlation",
  xlab = "Time to maximum biofilm (hours)",
  ylab = "Maximum cell index",
  point_fill = "#BFD7EA",
  line_colour = "#5B8DB8",
  ribbon_fill = "#E8F1F8"
)

plot_time #Figure 1B

#Maximum biofilm vs total biofilm burden
plot_auc <- make_biofilm_scatter(
  data = biofilm_all_features,
  xvar = "maxcellindex_from_full",
  yvar = "auc_0_12",
  subtitle = "Linear regression line with Spearman correlation",
  xlab = "Maximum cell index",
  ylab = "AUC (0–12 hours)",
  point_fill = "#A6CEE3",
  line_colour = "#1F78B4",
  ribbon_fill = "#D9ECF5"
)

plot_auc #Figure 1C

########################################################################
#12. Creating table with all the kinetics results
########################################################################
final_biofilm_results <- biofilm_all_features %>%
  arrange(desc(maxcellindex_from_full)) %>%
  dplyr::select(
    sample_number,
    maxcellindex_from_full,
    time_to_max,
    auc_0_12,
    auc_0_3,
    auc_3_6,
    auc_6_12,
    max_growth_rate,
    time_of_max_growth_rate,
    time_to_half_max,
    end_cellindex,
    decline_after_peak,
    biofilm_capacity_group,
    speed_group,
    cluster
  )

write_csv(final_biofilm_results, "biofilm_kinetic_features.csv")

################################################
#13. FIGURE 1
################################################
library(ggplot2)
library(patchwork)

layout <- "
ABC
DDD
"

figure1 <- (
  biofilm_curve + plot_time + plot_auc + biofilm_groups
) +
  plot_layout(
    design = layout,
    heights = c(1, 1.2)
  ) +
  plot_annotation(
    title = "Biofilm formation dynamics",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  )

figure1

ggsave(
  filename = "figure1.svg",
  plot = figure1,
  width = 20,
  height = 10,
  units = "in"
)

################################################
#C. FLUORIDE EXPERIMENT
################################################
#. Fluoride was added at 0 hours for the test group. 
#control groups included. 14 in each group, with paired samples
################################################
###################################################
# A) CLEANING INSPECTING AND FORMATTING THE DATA
#1. IMPORT AND INSPECT THE DATA
###################################################

library(tidyverse)
library(ggpubr)
library(mgcv)
library(patchwork)

fluoride_raw <- read_csv("fluoride_biofilm.csv", show_col_types = FALSE)

glimpse(fluoride_raw)
names(fluoride_raw)

#pairing fluoride groups in order
time_col <- names(fluoride_raw)[1]
sample_cols <- names(fluoride_raw)[-1]

# Number of paired samples
n_pairs <- length(sample_cols) / 2

control_cols <- sample_cols[1:n_pairs]
fluoride_cols <- sample_cols[(n_pairs + 1):length(sample_cols)]

pair_key <- tibble(
  sample_number = control_cols,
  control_col = control_cols,
  fluoride_col = fluoride_cols
)

pair_key
####################################################
#2. Converting to long format
####################################################
control_long <- fluoride_raw %>%
  dplyr::select(all_of(time_col), all_of(control_cols)) %>%
  pivot_longer(
    cols = all_of(control_cols),
    names_to = "sample_number",
    values_to = "cell_index"
  ) %>%
  rename(Time_Hour = all_of(time_col)) %>%
  mutate(condition = "Control")

fluoride_long <- fluoride_raw %>%
  dplyr::select(all_of(time_col), all_of(fluoride_cols)) %>%
  pivot_longer(
    cols = all_of(fluoride_cols),
    names_to = "fluoride_col",
    values_to = "cell_index"
  ) %>%
  rename(Time_Hour = all_of(time_col)) %>%
  left_join(pair_key, by = "fluoride_col") %>%
  dplyr::select(Time_Hour, sample_number, cell_index, condition = fluoride_col) %>%
  mutate(condition = "Fluoride")

fluoride_biofilm_long <- bind_rows(control_long, fluoride_long) %>%
  mutate(
    sample_number = factor(sample_number),
    condition = factor(condition, levels = c("Control", "Fluoride"))
  ) %>%
  arrange(sample_number, condition, Time_Hour)

head(fluoride_biofilm_long)
####################################################
#3. Plotting average fluoride versus control curves
####################################################
fluoride_summary <- fluoride_biofilm_long %>%
  group_by(condition, Time_Hour) %>%
  summarise(
    mean_cell_index = mean(cell_index, na.rm = TRUE),
    sd_cell_index = sd(cell_index, na.rm = TRUE),
    n = sum(!is.na(cell_index)),
    se_cell_index = sd_cell_index / sqrt(n),
    lower_ci = mean_cell_index - 1.96 * se_cell_index,
    upper_ci = mean_cell_index + 1.96 * se_cell_index,
    .groups = "drop"
  )

fluoride_curve <- ggplot(
  fluoride_summary,
  aes(x = Time_Hour, y = mean_cell_index, colour = condition, fill = condition)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(linewidth = 1.4) +
  scale_colour_manual(values = c(
    "Control" = "#1F78B4",
    "Fluoride" = "#C97C8A"
  )) +
  scale_fill_manual(values = c(
    "Control" = "#A6CEE3",
    "Fluoride" = "#F4C2C2"
  )) +
  labs(
    subtitle = "Mean cell index with 95% confidence interval",
    x = "Time (hours)",
    y = "Mean cell index",
    colour = "Condition",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

fluoride_curve #Figure 2A

####################################################
#4. Calculate the biofilm kinetics
#maximum biofilm capacity
#time to maximum biofilm
#AUC, total biofilm burden
#maximum growth rate
#time to half maximum
#final cell index
####################################################
trapz_auc <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

first_crossing_time <- function(t, y, threshold) {
  idx <- which(y >= threshold)
  if (length(idx) == 0) {
    return(NA_real_)
  } else {
    return(t[min(idx)])
  }
}

fluoride_features <- fluoride_biofilm_long %>%
  group_by(sample_number, condition) %>%
  group_modify(~ {
    
    d <- .x %>%
      arrange(Time_Hour)
    
    t <- d$Time_Hour
    y <- d$cell_index
    
    max_y <- max(y, na.rm = TRUE)
    max_index <- which.max(y)
    
    slopes <- diff(y) / diff(t)
    slope_time <- (head(t, -1) + tail(t, -1)) / 2
    
    tibble(
      maxcellindex = max_y,
      time_to_max = t[max_index],
      auc_0_12 = trapz_auc(t, y),
      max_growth_rate = max(slopes, na.rm = TRUE),
      time_of_max_growth_rate = slope_time[which.max(slopes)],
      time_to_half_max = first_crossing_time(t, y, 0.5 * max_y),
      end_cellindex = y[length(y)]
    )
  }) %>%
  ungroup()

fluoride_features

####################################################
#5. Creating paired fluoride effect variables
####################################################
fluoride_features_wide <- fluoride_features %>%
  pivot_wider(
    names_from = condition,
    values_from = c(
      maxcellindex,
      time_to_max,
      auc_0_12,
      max_growth_rate,
      time_of_max_growth_rate,
      time_to_half_max,
      end_cellindex
    ),
    names_sep = "_"
  )

fluoride_effects <- fluoride_features_wide %>%
  mutate(
    maxcellindex_diff = maxcellindex_Fluoride - maxcellindex_Control,
    maxcellindex_percent_change =
      100 * (maxcellindex_Fluoride - maxcellindex_Control) / maxcellindex_Control,
    maxcellindex_percent_inhibition =
      100 * (maxcellindex_Control - maxcellindex_Fluoride) / maxcellindex_Control,
    
    auc_diff = auc_0_12_Fluoride - auc_0_12_Control,
    auc_percent_change =
      100 * (auc_0_12_Fluoride - auc_0_12_Control) / auc_0_12_Control,
    auc_percent_inhibition =
      100 * (auc_0_12_Control - auc_0_12_Fluoride) / auc_0_12_Control,
    
    growth_rate_diff = max_growth_rate_Fluoride - max_growth_rate_Control,
    growth_rate_percent_inhibition =
      100 * (max_growth_rate_Control - max_growth_rate_Fluoride) / max_growth_rate_Control
  )

fluoride_effects

##############################################################################
#6. Paired statistical test
#using paired t test, because the same sample is measured with and without fluoride. 
##############################################################################
feature_vars <- c(
  "maxcellindex",
  "auc_0_12",
  "time_to_max",
  "max_growth_rate",
  "time_to_half_max",
  "end_cellindex"
)

fluoride_paired_tests <- map_dfr(feature_vars, function(v) {
  
  control_var <- paste0(v, "_Control")
  fluoride_var <- paste0(v, "_Fluoride")
  
  d <- fluoride_features_wide %>%
    dplyr::select(
      sample_number,
      control = all_of(control_var),
      fluoride = all_of(fluoride_var)
    ) %>%
    drop_na()
  
  test <- wilcox.test(
    d$fluoride,
    d$control,
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    feature = v,
    n = nrow(d),
    control_median = median(d$control, na.rm = TRUE),
    fluoride_median = median(d$fluoride, na.rm = TRUE),
    median_difference_fluoride_minus_control =
      median(d$fluoride - d$control, na.rm = TRUE),
    median_percent_change =
      median(100 * (d$fluoride - d$control) / d$control, na.rm = TRUE),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_value_adj_BH = p.adjust(p_value, method = "BH")
  )

fluoride_paired_tests

####################################################
# 7. Function for the Paired box plots 
####################################################
make_paired_boxplot <- function(feature, ylab, title) {
  
  control_var <- paste0(feature, "_Control")
  fluoride_var <- paste0(feature, "_Fluoride")
  
  p_value <- wilcox.test(
    fluoride_features_wide[[fluoride_var]],
    fluoride_features_wide[[control_var]],
    paired = TRUE,
    exact = FALSE
  )$p.value
  
  p_label <- paste0(
    "Paired Wilcoxon p = ",
    format.pval(p_value, digits = 3, eps = 0.001)
  )
  
  plot_data <- fluoride_features %>%
    dplyr::select(sample_number, condition, value = all_of(feature)) %>%
    mutate(condition = factor(condition, levels = c("Control", "Fluoride")))
  
  ggplot(plot_data, aes(x = condition, y = value, fill = condition)) +
    geom_boxplot(
      alpha = 0.85,
      width = 0.6,
      colour = "grey25",
      outlier.shape = NA
    ) +
    geom_line(
      aes(group = sample_number),
      colour = "grey60",
      alpha = 0.6,
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.08,
      size = 2.2,
      alpha = 0.8,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(plot_data$value, na.rm = TRUE) * 1.08,
      label = p_label,
      fontface = "bold",
      size = 4
    ) +
    scale_fill_manual(values = c(
      "Control" = "#A6CEE3",
      "Fluoride" = "#F4C2C2"
    )) +
    labs(
      x = "",
      y = ylab
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
}

# maximum biofilm capacity
plot_fluoride_max <- make_paired_boxplot(
  feature = "maxcellindex",
  ylab = "Maximum cell index"
)

plot_fluoride_max #Figure 2B

####################################################
#8. AUC, total biofilm burden
####################################################
plot_fluoride_auc <- make_paired_boxplot(
  feature = "auc_0_12",
  ylab = "AUC (0–12 hours)"
)

plot_fluoride_auc #Figure 2C

####################################################
#9. Time specific paired test
#This tests whether fluoride and control differ at each time point.
####################################################
time_specific_tests <- fluoride_biofilm_long %>%
  pivot_wider(
    id_cols = c(sample_number, Time_Hour),
    names_from = condition,
    values_from = cell_index
  ) %>%
  group_by(Time_Hour) %>%
  summarise(
    n = sum(!is.na(Control) & !is.na(Fluoride)),
    median_control = median(Control, na.rm = TRUE),
    median_fluoride = median(Fluoride, na.rm = TRUE),
    median_difference = median(Fluoride - Control, na.rm = TRUE),
    median_percent_inhibition =
      median(100 * (Control - Fluoride) / Control, na.rm = TRUE),
    p_value = wilcox.test(
      Fluoride,
      Control,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_value_adj_BH = p.adjust(p_value, method = "BH"),
    significant = p_value_adj_BH < 0.05
  )

time_specific_tests

fluoride_time<- ggplot(
  time_specific_tests,
  aes(x = Time_Hour, y = median_percent_inhibition)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(colour = "#C97C8A", linewidth = 1.3) +
  geom_point(
    aes(shape = significant),
    size = 2.3,
    colour = "black",
    fill = "#F4C2C2"
  ) +
  labs(
    x = "Time (hours)",
    y = "Median inhibition (%)",
    shape = "BH-adjusted p < 0.05"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )


fluoride_time #Figure 2D


##################################################################
# D. FLUORIDE 7 HOUR EXPIREMENT
###############################################################
################################################
#. Fluoride was added at 7 hours for the test group. 
#control groups included. 
################################################
###################################################
# CLEANING INSPECTING AND FORMATTING THE DATA
#1. IMPORT AND INSPECT THE DATA
###################################################
library(tidyverse)
library(ggpubr)
library(mgcv)
library(patchwork)

fluoride7_raw <- read_csv("fluoride_7h_biofilm.csv", show_col_types = FALSE)

# Rename first column as Time_Hour
names(fluoride7_raw)[1] <- "Time_Hour"

fluoride7_long <- fluoride7_raw %>%
  pivot_longer(
    cols = -Time_Hour,
    names_to = "sample_raw",
    values_to = "cell_index"
  ) %>%
  mutate(
    condition = if_else(str_detect(sample_raw, "F$"), "Fluoride at 7h", "Control"),
    sample_number = str_remove(sample_raw, "F$"),
    sample_number = factor(sample_number),
    condition = factor(condition, levels = c("Control", "Fluoride at 7h")),
    post7 = if_else(Time_Hour >= 7, "Post-7h", "Pre-7h"),
    post7 = factor(post7, levels = c("Pre-7h", "Post-7h")),
    time_after_7 = pmax(Time_Hour - 7, 0)
  )

head(fluoride7_long)

################################################################################################
#2. Plotting average biofilm curves with fluoride addition time
################################################################################################
fluoride7_summary <- fluoride7_long %>%
  group_by(condition, Time_Hour) %>%
  summarise(
    mean_cell_index = mean(cell_index, na.rm = TRUE),
    sd_cell_index = sd(cell_index, na.rm = TRUE),
    n = sum(!is.na(cell_index)),
    se_cell_index = sd_cell_index / sqrt(n),
    lower_ci = mean_cell_index - 1.96 * se_cell_index,
    upper_ci = mean_cell_index + 1.96 * se_cell_index,
    .groups = "drop"
  )

fluoride7_curve <- ggplot(
  fluoride7_summary,
  aes(x = Time_Hour, y = mean_cell_index, colour = condition, fill = condition)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(linewidth = 1.4) +
  geom_vline(
    xintercept = 7,
    linetype = "dashed",
    linewidth = 0.8,
    colour = "grey30"
  ) +
  annotate(
    "text",
    x = 7.1,
    y = max(fluoride7_summary$upper_ci, na.rm = TRUE),
    label = "Fluoride added",
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold"
  ) +
  scale_colour_manual(values = c(
    "Control" = "#1F78B4",
    "Fluoride at 7h" = "#C97C8A"
  )) +
  scale_fill_manual(values = c(
    "Control" = "#A6CEE3",
    "Fluoride at 7h" = "#F4C2C2"
  )) +
  labs(
    subtitle = "Fluoride was added at 7 hours",
    x = "Time (hours)",
    y = "Mean cell index",
    colour = "Condition",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

fluoride7_curve #Figure 2E

################################################
#3. Calculating pre and post 7 hours features
################################################
trapz_auc <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

auc_between <- function(t, y, lower, upper) {
  keep <- t >= lower & t <= upper
  
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  
  trapz_auc(t[keep], y[keep])
}

slope_between <- function(t, y, lower, upper) {
  keep <- t >= lower & t <= upper
  
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  
  coef(lm(y[keep] ~ t[keep]))[2]
}

fluoride7_features <- fluoride7_long %>%
  group_by(sample_number, condition) %>%
  group_modify(~ {
    
    d <- .x %>% arrange(Time_Hour)
    
    t <- d$Time_Hour
    y <- d$cell_index
    
    idx_7h <- which.min(abs(t - 7))
    y_7h <- y[idx_7h]
    
    tibble(
      maxcellindex_overall = max(y, na.rm = TRUE),
      time_to_max_overall = t[which.max(y)],
      
      auc_0_7 = auc_between(t, y, 0, 7),
      auc_7_12 = auc_between(t, y, 7, 12),
      auc_0_12 = auc_between(t, y, 0, 12),
      
      cellindex_7h = y_7h,
      end_cellindex = y[length(y)],
      change_7_to_12 = y[length(y)] - y_7h,
      percent_change_7_to_12 = 100 * (y[length(y)] - y_7h) / y_7h,
      
      slope_0_7 = slope_between(t, y, 0, 7),
      slope_7_12 = slope_between(t, y, 7, 12),
      slope_change_post_minus_pre = slope_between(t, y, 7, 12) - slope_between(t, y, 0, 7)
    )
  }) %>%
  ungroup()

fluoride7_features

################################################
#4. Creating paired dataset for comparison
################################################
fluoride7_features_wide <- fluoride7_features %>%
  pivot_wider(
    names_from = condition,
    values_from = c(
      maxcellindex_overall,
      time_to_max_overall,
      auc_0_7,
      auc_7_12,
      auc_0_12,
      cellindex_7h,
      end_cellindex,
      change_7_to_12,
      percent_change_7_to_12,
      slope_0_7,
      slope_7_12,
      slope_change_post_minus_pre
    ),
    names_sep = "_"
  )

fluoride7_effects <- fluoride7_features_wide %>%
  mutate(
    post_auc_diff = `auc_7_12_Fluoride at 7h` - auc_7_12_Control,
    post_auc_percent_change =
      100 * (`auc_7_12_Fluoride at 7h` - auc_7_12_Control) / auc_7_12_Control,
    
    post_slope_diff = `slope_7_12_Fluoride at 7h` - slope_7_12_Control,
    
    change_7_to_12_diff =
      `change_7_to_12_Fluoride at 7h` - change_7_to_12_Control,
    
    end_cellindex_diff =
      `end_cellindex_Fluoride at 7h` - end_cellindex_Control
  )

fluoride7_effects

################################################
#5. Paires tests focussed on post 7 hour outcomes
################################################
post7_vars <- c(
  "auc_7_12",
  "change_7_to_12",
  "percent_change_7_to_12",
  "slope_7_12",
  "end_cellindex"
)

fluoride7_post_tests <- map_dfr(post7_vars, function(v) {
  
  control_var <- paste0(v, "_Control")
  fluoride_var <- paste0(v, "_Fluoride at 7h")
  
  d <- fluoride7_features_wide %>%
    dplyr::select(
      sample_number,
      control = all_of(control_var),
      fluoride = all_of(fluoride_var)
    ) %>%
    drop_na()
  
  test <- wilcox.test(
    d$fluoride,
    d$control,
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    feature = v,
    n = nrow(d),
    control_median = median(d$control, na.rm = TRUE),
    fluoride_median = median(d$fluoride, na.rm = TRUE),
    median_difference_fluoride_minus_control =
      median(d$fluoride - d$control, na.rm = TRUE),
    median_percent_change =
      median(100 * (d$fluoride - d$control) / d$control, na.rm = TRUE),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_value_adj_BH = p.adjust(p_value, method = "BH")
  )

fluoride7_post_tests

#auc_7_12
#change_7_to_12
#slope_7_12
#end_cellindex

################################################
#6. Paired box plots
#Final cell index and AUC from 7 to 12 hours
################################################

make_fluoride7_boxplot <- function(feature, ylab, title) {
  
  control_var <- paste0(feature, "_Control")
  fluoride_var <- paste0(feature, "_Fluoride at 7h")
  
  d_wide <- fluoride7_features_wide %>%
    dplyr::select(
      sample_number,
      control = all_of(control_var),
      fluoride = all_of(fluoride_var)
    ) %>%
    drop_na()
  
  p_value <- wilcox.test(
    d_wide$fluoride,
    d_wide$control,
    paired = TRUE,
    exact = FALSE
  )$p.value
  
  p_label <- paste0(
    "Paired Wilcoxon p = ",
    format.pval(p_value, digits = 3, eps = 0.001)
  )
  
  plot_data <- d_wide %>%
    pivot_longer(
      cols = c(control, fluoride),
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(
      condition = recode(
        condition,
        control = "Control",
        fluoride = "Fluoride at 7h"
      ),
      condition = factor(condition, levels = c("Control", "Fluoride at 7h"))
    )
  
  ggplot(plot_data, aes(x = condition, y = value, fill = condition)) +
    geom_boxplot(
      alpha = 0.85,
      width = 0.6,
      colour = "grey25",
      outlier.shape = NA
    ) +
    geom_line(
      aes(group = sample_number),
      colour = "grey60",
      alpha = 0.6,
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.08,
      size = 2.2,
      alpha = 0.8,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(plot_data$value, na.rm = TRUE) * 1.08,
      label = p_label,
      fontface = "bold",
      size = 4
    ) +
    scale_fill_manual(values = c(
      "Control" = "#A6CEE3",
      "Fluoride at 7h" = "#F4C2C2"
    )) +
    labs(
      x = "",
      y = ylab
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
}

plot_end_cellindex <- make_fluoride7_boxplot(
  feature = "end_cellindex",
  ylab = "Final cell index",
)

plot_end_cellindex #Figure 2F

plot_post_auc <- make_fluoride7_boxplot(
  feature = "auc_7_12",
  ylab = "AUC from 7 to 12 hours",
)

plot_post_auc #Figure 2G

######################################################################################################
#7. Time specific paired test
#This tests whether fluoride (after being added at 7 hours) and control differ at each time point.
######################################################################################################
fluoride7_time_tests <- fluoride7_long %>%
  pivot_wider(
    id_cols = c(sample_number, Time_Hour),
    names_from = condition,
    values_from = cell_index
  ) %>%
  group_by(Time_Hour) %>%
  summarise(
    n = sum(!is.na(Control) & !is.na(`Fluoride at 7h`)),
    median_control = median(Control, na.rm = TRUE),
    median_fluoride = median(`Fluoride at 7h`, na.rm = TRUE),
    median_difference = median(`Fluoride at 7h` - Control, na.rm = TRUE),
    p_value = wilcox.test(
      `Fluoride at 7h`,
      Control,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_value_adj_BH = p.adjust(p_value, method = "BH"),
    significant = p_value_adj_BH < 0.05
  )

fluoride7_time_tests

fluoride7_difference_plot <- ggplot(
  fluoride7_time_tests,
  aes(x = Time_Hour, y = median_difference)
) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 7, linetype = "dashed", colour = "grey30") +
  geom_line(colour = "#C97C8A", linewidth = 1.3) +
  geom_point(
    aes(shape = significant),
    size = 2.5,
    colour = "black",
    fill = "#F4C2C2"
  ) +
  labs(
    x = "Time (hours)",
    y = "Median difference: Fluoride − Control",
    shape = "BH-adjusted p < 0.05"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

fluoride7_difference_plot # Figure 2H

######################################################################################################
#8. FIGURE 2
######################################################################################################

Figure_2 <- (
  fluoride_curve |plot_fluoride_max | plot_fluoride_auc |fluoride_time
) / (
  fluoride7_curve | plot_end_cellindex | plot_post_auc |fluoride7_difference_plot
) +
  plot_annotation(
    title = "Effect of Fluoride on Biofilm Formation",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

Figure_2


ggsave(
  filename = "Figure_2.svg",
  plot = Figure_2,
  width = 25,
  height = 10,
  units = "in"
)

##################################################################
# E. ARGININE EXPIREMENT
###############################################################
################################################
#. Arginie was added for the test group. 
#control groups included. 
################################################
###################################################
# CLEANING INSPECTING AND FORMATTING THE DATA
#1. IMPORT AND INSPECT THE DATA
###################################################
library(tidyverse)
library(ggpubr)
library(mgcv)
library(patchwork)

arginine_raw <- read_csv("arginine_biofilm.csv", show_col_types = FALSE)

# Rename first column as Time_Hour
names(arginine_raw)[1] <- "Time_Hour"

arginine_long <- arginine_raw %>%
  pivot_longer(
    cols = -Time_Hour,
    names_to = "sample_raw",
    values_to = "cell_index"
  ) %>%
  mutate(
    condition = if_else(str_detect(sample_raw, "A$"), "Arginine", "Control"),
    sample_number = str_remove(sample_raw, "A$"),
    sample_number = factor(sample_number),
    condition = factor(condition, levels = c("Control", "Arginine"))
  )

head(arginine_long)

###################################################
#2. Plot average control vs arginien biofim curves
###################################################
arginine_summary <- arginine_long %>%
  group_by(condition, Time_Hour) %>%
  summarise(
    mean_cell_index = mean(cell_index, na.rm = TRUE),
    sd_cell_index = sd(cell_index, na.rm = TRUE),
    n = sum(!is.na(cell_index)),
    se_cell_index = sd_cell_index / sqrt(n),
    lower_ci = mean_cell_index - 1.96 * se_cell_index,
    upper_ci = mean_cell_index + 1.96 * se_cell_index,
    .groups = "drop"
  )

arginine_curve <- ggplot(
  arginine_summary,
  aes(x = Time_Hour, y = mean_cell_index, colour = condition, fill = condition)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(linewidth = 1.4) +
  scale_colour_manual(values = c(
    "Control" = "#1F78B4",
    "Arginine" = "#7DAA72"
  )) +
  scale_fill_manual(values = c(
    "Control" = "#A6CEE3",
    "Arginine" = "#CDEAC0"
  )) +
  labs(
    subtitle = "Mean cell index with 95% confidence interval",
    x = "Time (hours)",
    y = "Mean cell index",
    colour = "Condition",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

arginine_curve #Figure 3A
###################################################
#3. Derive biofilm kinetic features
###################################################
trapz_auc <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

first_crossing_time <- function(t, y, threshold) {
  idx <- which(y >= threshold)
  if (length(idx) == 0) {
    return(NA_real_)
  } else {
    return(t[min(idx)])
  }
}

arginine_features <- arginine_long %>%
  group_by(sample_number, condition) %>%
  group_modify(~ {
    
    d <- .x %>%
      arrange(Time_Hour)
    
    t <- d$Time_Hour
    y <- d$cell_index
    
    max_y <- max(y, na.rm = TRUE)
    max_index <- which.max(y)
    
    slopes <- diff(y) / diff(t)
    slope_time <- (head(t, -1) + tail(t, -1)) / 2
    
    tibble(
      maxcellindex = max_y,
      time_to_max = t[max_index],
      auc_0_12 = trapz_auc(t, y),
      max_growth_rate = max(slopes, na.rm = TRUE),
      time_of_max_growth_rate = slope_time[which.max(slopes)],
      time_to_half_max = first_crossing_time(t, y, 0.5 * max_y),
      end_cellindex = y[length(y)]
    )
  }) %>%
  ungroup()

arginine_features

###################################################
#4. reate paired arginine effect variables
###################################################

arginine_features_wide <- arginine_features %>%
  pivot_wider(
    names_from = condition,
    values_from = c(
      maxcellindex,
      time_to_max,
      auc_0_12,
      max_growth_rate,
      time_of_max_growth_rate,
      time_to_half_max,
      end_cellindex
    ),
    names_sep = "_"
  )

arginine_effects <- arginine_features_wide %>%
  mutate(
    maxcellindex_diff = maxcellindex_Arginine - maxcellindex_Control,
    maxcellindex_percent_change =
      100 * (maxcellindex_Arginine - maxcellindex_Control) / maxcellindex_Control,
    maxcellindex_percent_reduction =
      100 * (maxcellindex_Control - maxcellindex_Arginine) / maxcellindex_Control,
    
    auc_diff = auc_0_12_Arginine - auc_0_12_Control,
    auc_percent_change =
      100 * (auc_0_12_Arginine - auc_0_12_Control) / auc_0_12_Control,
    auc_percent_reduction =
      100 * (auc_0_12_Control - auc_0_12_Arginine) / auc_0_12_Control,
    
    growth_rate_diff = max_growth_rate_Arginine - max_growth_rate_Control,
    growth_rate_percent_reduction =
      100 * (max_growth_rate_Control - max_growth_rate_Arginine) / max_growth_rate_Control
  )

arginine_effects

###################################################
#5. Paired statistical test
###################################################
feature_vars <- c(
  "maxcellindex",
  "auc_0_12",
  "time_to_max",
  "max_growth_rate",
  "time_to_half_max",
  "end_cellindex"
)

arginine_paired_tests <- map_dfr(feature_vars, function(v) {
  
  control_var <- paste0(v, "_Control")
  arginine_var <- paste0(v, "_Arginine")
  
  d <- arginine_features_wide %>%
    dplyr::select(
      sample_number,
      control = all_of(control_var),
      arginine = all_of(arginine_var)
    ) %>%
    drop_na()
  
  test <- wilcox.test(
    d$arginine,
    d$control,
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    feature = v,
    n = nrow(d),
    control_median = median(d$control, na.rm = TRUE),
    arginine_median = median(d$arginine, na.rm = TRUE),
    median_difference_arginine_minus_control =
      median(d$arginine - d$control, na.rm = TRUE),
    median_percent_change =
      median(100 * (d$arginine - d$control) / d$control, na.rm = TRUE),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_value_adj_BH = p.adjust(p_value, method = "BH")
  )

arginine_paired_tests

###################################################
#6.Paired box plots
###################################################
make_paired_boxplot <- function(feature, ylab, title) {
  
  control_var <- paste0(feature, "_Control")
  arginine_var <- paste0(feature, "_Arginine")
  
  p_value <- wilcox.test(
    arginine_features_wide[[arginine_var]],
    arginine_features_wide[[control_var]],
    paired = TRUE,
    exact = FALSE
  )$p.value
  
  p_label <- paste0(
    "Paired Wilcoxon p = ",
    format.pval(p_value, digits = 3, eps = 0.001)
  )
  
  plot_data <- arginine_features %>%
    dplyr::select(sample_number, condition, value = all_of(feature)) %>%
    mutate(condition = factor(condition, levels = c("Control", "Arginine")))
  
  ggplot(plot_data, aes(x = condition, y = value, fill = condition)) +
    geom_boxplot(
      alpha = 0.85,
      width = 0.6,
      colour = "grey25",
      outlier.shape = NA
    ) +
    geom_line(
      aes(group = sample_number),
      colour = "grey60",
      alpha = 0.6,
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.08,
      size = 2.2,
      alpha = 0.8,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(plot_data$value, na.rm = TRUE) * 1.08,
      label = p_label,
      fontface = "bold",
      size = 4
    ) +
    scale_fill_manual(values = c(
      "Control" = "#A6CEE3",
      "Arginine" = "#CDEAC0"
    )) +
    labs(
      x = "",
      y = ylab
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
}

plot_arginine_max <- make_paired_boxplot(
  feature = "maxcellindex",
  ylab = "Maximum cell index",
)

plot_arginine_max #Figure 3C

plot_arginine_auc <- make_paired_boxplot(
  feature = "auc_0_12",
  ylab = "AUC (0–12 hours)",
)

plot_arginine_auc #Figure 3D


plot_arginine_time <- make_paired_boxplot(
  feature = "time_to_max",
  ylab = "Time to maximum biofilm (hours)",
)

plot_arginine_time #Figure 3E

######################################################################################################
#7. Time specific paired test
#This tests whether Arginine  and control samples differ at each time point.
######################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

arginine_time_tests <- arginine_long %>%
  pivot_wider(
    id_cols = c(sample_number, Time_Hour),
    names_from = condition,
    values_from = cell_index
  ) %>%
  group_by(Time_Hour) %>%
  summarise(
    n = sum(!is.na(Control) & !is.na(Arginine)),
    
    median_control = median(Control, na.rm = TRUE),
    median_arginine = median(Arginine, na.rm = TRUE),
    
    median_difference = median(Arginine - Control, na.rm = TRUE),
    
    median_percent_reduction =
      median(100 * (Control - Arginine) / Control, na.rm = TRUE),
    
    p_value = wilcox.test(
      Arginine,
      Control,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  mutate(
    p_value_adj_BH = p.adjust(p_value, method = "BH"),
    significant = p_value_adj_BH < 0.05
  )

arginine_time <- ggplot(
  arginine_time_tests,
  aes(x = Time_Hour, y = median_percent_reduction)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_line(
    colour = "#7DAA72",
    linewidth = 1.3
  ) +
  geom_point(
    aes(shape = significant),
    size = 2.8,
    colour = "black",
    fill = "#CDEAC0",
    stroke = 0.7
  ) +
  scale_shape_manual(
    values = c(
      "FALSE" = 21,
      "TRUE" = 24
    )
  ) +
  labs(
    x = "Time (hours)",
    y = "Median reduction (%)",
    shape = "BH-adjusted p < 0.05"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "grey20"),
    panel.grid.major = element_line(
      colour = "#E6E6E6",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

arginine_time #Figure 3B


######################################################################################################
#8. FIGURE 3
######################################################################################################

arginine_figure <- (
  arginine_curve| arginine_time
) / (
  plot_arginine_max | plot_arginine_auc | plot_arginine_time
) +
  plot_annotation(
    title = "Effect of Arginine on Biofilm Formation",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

arginine_figure

ggsave(
  filename = "arginine_figure.svg",
  plot = arginine_figure,
  width = 20,
  height = 10,
  units = "in"
)

############################## END ########################################################################
