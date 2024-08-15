#Clear existing data and graphics
rm(list=ls())
graphics.off()

setwd("~/OMT_April2024/qiime2analysis_clean/")

#Library the packages
library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.

library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq

library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq

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

library(microbiomeSeq) # Data analysis and visualization

library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown

library(ranacapa) # Data analysis 

library(grid) # support data visualization

library(gridExtra)  # support data visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(png) # Figure download

library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'

library(ggpubr) # publication quality figures, based on ggplot2

library(RColorBrewer) # nice color options

library(microbiomeutilities) # some utility tools 

library(Maaslin2) ##for differential abundance

library(DT) #for interactive tables

library(reshape2)

library(scales)

library(data.table)

library(Biostrings)
####################################
####Creating a phyloseq object#####
###################################
physeq<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree.qza",
  "taxonomy.qza", 
  metadata= "metadata.tsv"
)
physeq

###Summarising the phyloseq object
summarize_phyloseq(physeq) #gives the brief summary of the samples and the data
print_ps(physeq) #prints the same summary
summary(sample_sums(physeq)) #fives the min, max, median, mean for the reads

#Assessing the phyloseq object
ntaxa(physeq) # gives the number of ASVs
nsamples(physeq) #number of samples
sample_names(physeq)[1:5] #overview of sample names 
rank_names(physeq)  #Rank names
sample_variables(physeq)  # the variables in the metadata

otu_table(physeq)[1:5, 1:5] ## the feature table 

tax_table(physeq)[1:5, 1:4]  # the taxonomy table

## Check the data

#```{r check, echo=TRUE, message=FALSE, warning=FALSE}
ps1<- physeq  #creating a copy of the original phyloseq object
any(taxa_sums(ps1) == 0) # check if any OTUs are not present in any samples
# The answer was FALSE, there are no OTU's that are not found in any samples. 
ps1a <- prune_taxa(taxa_sums(ps1) > 0, ps1)

# check again if any OTUs are not present in any samples
any(taxa_sums(ps1a) == 0)

# subtract the number of OTUs in original (ps1) with number of OTUs in new phyloseq (ps1a) object.
# no. of OTUs in original 
ntaxa(ps1) ##212
# no. of OTUs in new #210
ntaxa(ps1a)
ntaxa(ps1) - ntaxa(ps1a) ## 0 removed


###Use the interactive table to check under which taxonomic rank chloroplast and Mitrochondria are mentioned. Specify them Class or Order in the code below.    

ps1a <- subset_taxa(ps1,Class!="c__Chloroplast")

ps1b <- subset_taxa(ps1a,Order!="o__Mitochondria")

ps1 <- merge_phyloseq(ps1a, ps1b)
###In the next step we are sorting the samples
sort(sample_sums(ps1))

tax_table_df <- as.data.frame(tax_table(physeq)) ##converting to data frame
datatable(tax_table_df)# the table is interactive you can scroll and search through it for details.
tax_table_df <- as.data.frame(physeq@tax_table)
View(tax_table_df)
# to save the taxa table
## Save new metadata file
write.table(tax_table_df, file = "taxa_table.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

###using PS1 for further analysis########
##########################################
#####    ALPHA DIVERSITY #############
#########################################
richness <- estimate_richness(ps1)
head(richness)

write.table(richness,
            file = "richness.tsv", 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
            
###Check the normal distribution curve for alpha diversity
hist(richness$Shannon, main="Shannon index", xlab="")
hist(richness$Observed, main="Observed features", xlab="")
hist(richness$Chao1, main="Chao 1 index", xlab="")
# the distribution was not normal, using non-parametric test. For two groups comparison, Wilcoxon test and for 
# multiple groups using Kruskal Wallis rank sum  test.

###combining the richness file with the metadata.
alpha_diversity <- data.frame(sample_data(ps1), richness)
write.table(alpha_diversity,
            file = "alpha_diversity.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



#a_my_comparisons <- list( c("YY", "MM"), c("CA", "TT"), c("NN", "SS"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))



#####1. Age category NS
comps <- make_pairs(sample_data(ps1)$age_group)
print(comps)

#pairwise.wilcox.test(richness$Observed, sample_data(ps1)$age_group, p.adj = "bonf")
#kruskal.test(richness$Observed ~ sample_data(ps1)$age_group)
#pairwise.wilcox.test(richness$Shannon, sample_data(ps1)$age_group, p.adj = "bonf")
#kruskal.test(richness$Shannon ~ sample_data(ps1)$age_group)

wilcox.test(Observed ~ age_group, data = alpha_diversity)
wilcox.test(Shannon ~ age_group, data = alpha_diversity)
#,
#subset = (age_group %in% c("Adults", "Young")))

p1 <- plot_richness(ps1, x = "age_group", measures = c("Observed", "Shannon"), color = "age_group") +
  geom_boxplot(aes(fill = age_group), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Young adults" = "#D55E00", "Adults" = "#0072B2")) +  # Specify actual colors
  labs(x = "Age group") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = age_group), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


# Print the plot
print(p1)

######## 2. Sex NS
comps <- make_pairs(sample_data(ps1)$gender.factor.x)
print(comps)

wilcox.test(Observed ~ gender.factor.x, data = alpha_diversity)
wilcox.test(Shannon ~ gender.factor.x, data = alpha_diversity)

p2 <- plot_richness(ps1, x = "gender.factor.x)", measures = c("Observed", "Shannon"), color = "gender.factor.x") +
  geom_boxplot(aes(fill = gender.factor.x), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Male" = "#0072B2", "Female" = "#D55E00")) +  # Specify actual colors
  labs(x = "Gender") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = gender.factor.x), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6)) +  # Statistical comparison
  stat_compare_means(label.y = 6)  # Additional mean comparison

# Print the plot
print(p2)

#####3. Country_born S 
comps <- make_pairs(sample_data(ps1)$country_born.factor)
print(comps)
wilcox.test(Observed ~ country_born.factor, data = alpha_diversity)
wilcox.test(Shannon ~ country_born.factor, data = alpha_diversity)

# Assuming sample_data(ps1) has a 'country_born.factor' column
#sample_data(ps1)$country_born.factor <- factor(sample_data(ps1)$country_born.factor)

# Making pairwise comparisons for the 'country_born.factor'
#comps <- make_pairs(sample_data(ps1)$country_born.factor)
#print(comps)

p3 <- plot_richness(ps1, x = "country_born.factor", measures = c("Observed", "Shannon"), color = "country_born.factor") +
  geom_boxplot(aes(fill = country_born.factor), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Other" = "#D55E00", "Australia" = "#0072B2")) +  # Specify actual colors
  labs(x = "Country of Birth") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = country_born.factor), position = position_jitter(width = 0.2)) + 
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                       label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


# Print the plot
print(p3)

####. 4. Level of education NS
comps <- make_pairs(sample_data(ps1)$education_status)
print(comps)

wilcox.test(Observed ~ education_status, data = alpha_diversity)
wilcox.test(Shannon ~ education_status, data = alpha_diversity)

p4 <- plot_richness(ps1, x = "education_status", measures = c("Observed", "Shannon"), color = "education_status") +
  geom_boxplot(aes(fill = education_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Secondary education or less" = "#0072B2", "Tertiary education" = "#D55E00")) +  # Specify actual colors
  labs(x = "Level of Education") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = education_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

# Print the plot
print(p4)

###########5.  Current study status S
comps <- make_pairs(sample_data(ps1)$study_status)
print(comps)

wilcox.test(Observed ~ study_status, data = alpha_diversity)
wilcox.test(Shannon ~ study_status, data = alpha_diversity)


p5 <- plot_richness(ps1, x = "study_status", measures = c("Observed", "Shannon"), color = "study_status") +
  geom_boxplot(aes(fill = study_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Full time student" = "#D55E00", "Not studying" = "#0072B2")) +  # Specify actual colors
  labs(x = "Student status") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = study_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p5)
########### 6. Employment S
comps <- make_pairs(sample_data(ps1)$employment_status)
print(comps)

wilcox.test(Observed ~ employment_status, data = alpha_diversity)
wilcox.test(Shannon ~ employment_status, data = alpha_diversity)

p6 <- plot_richness(ps1, x = "employment_status", measures = c("Observed", "Shannon"), color = "employment_status") +
  geom_boxplot(aes(fill = employment_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Employed" = "#0072B2", "Not working" = "#D55E00")) +  # Specify actual colors
  labs(x = "Employment status") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = employment_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p6)

#############
###7. Physical activity NS
comps <- make_pairs(sample_data(ps1)$physical_activity_status)
print(comps)


#kruskal_result <- kruskal.test(Observed ~ physical_activity_status, data = alpha_diversity)
#print(kruskal_result)

kruskal.test(richness$Observed ~ sample_data(ps1)$physical_activity_status)
kruskal.test(richness$Shannon ~ sample_data(ps1)$physical_activity_status)

library(dunn.test)

# Perform Dunn's test for pairwise comparisons
dunn_pa_o <- dunn.test(alpha_diversity$Observed, alpha_diversity$physical_activity_status, 
                         method="bonferroni")
dunn_pa_s <- dunn.test(alpha_diversity$Shannon, alpha_diversity$physical_activity_status, 
                       method="bonferroni")
# Print the results
print(dunn_pa_o)
print(dunn_pa_s)


sample_data(ps1)$physical_activity_status <- factor(sample_data(ps1)$physical_activity_status,
                                                    levels = c("Minimal activity", "Moderate activity", "High activity"))




p7 <- plot_richness(ps1, x = "physical_activity_status", measures = c("Observed", "Shannon"), color = "physical_activity_status") +
  geom_boxplot(aes(fill = physical_activity_status), color = "black", outlier.color = NA) +
  scale_fill_manual(values = c("Minimal activity" = "#009E73", "Moderate activity" = "#0072B2", "High activity" = "#D55E00")) +
  labs(x = "Physical Activity Status") +
  theme_minimal() +
  theme(strip.background = element_blank(), legend.position = "none") +
  geom_jitter(aes(color = physical_activity_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

# Print the plot
print(p7)

###.8. Mental well being
comps <- make_pairs(sample_data(ps1)$mental_health_status)
print(comps)

wilcox.test(Observed ~ mental_health_status, data = alpha_diversity)
wilcox.test(Shannon ~ mental_health_status, data = alpha_diversity)

p8 <- plot_richness(ps1, x = "mental_health_status", measures = c("Observed", "Shannon"), color = "mental_health_status") +
  geom_boxplot(aes(fill = mental_health_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Low distress" = "#0072B2", "Severe distress" = "#D55E00")) +  # Specify actual colors
  labs(x = "Mental Health Status") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = mental_health_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p8)

######9. Self rated general health
comps <- make_pairs(sample_data(ps1)$generalhealth_status)
print(comps)
wilcox.test(Observed ~ generalhealth_status, data = alpha_diversity)
wilcox.test(Shannon ~ generalhealth_status, data = alpha_diversity)

p9 <- plot_richness(ps1, x = "generalhealth_status", measures = c("Observed", "Shannon"), color = "generalhealth_status") +
  geom_boxplot(aes(fill = generalhealth_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Excellent" = "#0072B2", "Good" = "#D55E00")) +  # Specify actual colors
  labs(x = "General Health Status") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = generalhealth_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p9)

#######10. Carbohydrate status
comps <- make_pairs(sample_data(ps1)$carbohydrate_status)
print(comps)

kruskal.test(richness$Observed ~ sample_data(ps1)$carbohydrate_status)
kruskal.test(richness$Shannon ~ sample_data(ps1)$carbohydrate_status)

dunn_carb_o <- dunn.test(alpha_diversity$Observed, alpha_diversity$carbohydrate_status, 
                       method="bonferroni")
dunn_carb_s <- dunn.test(alpha_diversity$Shannon, alpha_diversity$carbohydrate_status, 
                       method="bonferroni")
# Print the results
print(dunn_carb_o)
print(dunn_carb_s)


sample_data(ps1)$carbohydrate_status <- factor(sample_data(ps1)$carbohydrate_status,
                                                    levels = c("Low carb", "Average carb", "High carb"))

p10 <- plot_richness(ps1, x = "carbohydrate_status", measures = c("Observed", "Shannon"), color = "carbohydrate_status") +
  geom_boxplot(aes(fill = carbohydrate_status), color = "black", outlier.color = NA) +
  scale_fill_manual(values = c("Low carb" = "#009E73", "Average carb" = "#0072B2", "High carb" = "#D55E00")) +
  labs(x = "Daily Carbohydrate Intake") +
  theme_minimal() +
  theme(strip.background = element_blank(), legend.position = "none") +
  geom_jitter(aes(color = carbohydrate_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p10)

##########11.Daily protein intake
comps <- make_pairs(sample_data(ps1)$Protein_status)
print(comps)

kruskal.test(richness$Observed ~ sample_data(ps1)$Protein_status)
kruskal.test(richness$Shannon ~ sample_data(ps1)$Protein_status)

dunn_prot_o <- dunn.test(alpha_diversity$Observed, alpha_diversity$Protein_status, 
                         method="bonferroni")
dunn_prot_s <- dunn.test(alpha_diversity$Shannon, alpha_diversity$Protein_status, 
                         method="bonferroni")
# Print the results
print(dunn_prot_o)
print(dunn_prot_s)


sample_data(ps1)$Protein_status <- factor(sample_data(ps1)$Protein_status,
                                               levels = c("Low protein", "Average protein", "High protein"))

p11 <- plot_richness(ps1, x = "Protein_status", measures = c("Observed", "Shannon"), color = "Protein_status") +
  geom_boxplot(aes(fill = Protein_status), color = "black", outlier.color = NA) +
  scale_fill_manual(values = c("Low protein" = "#009E73", "Average protein" = "#0072B2", "High protein" = "#D55E00")) +
  labs(x = "Daily Protein Intake") +
  theme_minimal() +
  theme(strip.background = element_blank(), legend.position = "none") +
  geom_jitter(aes(color = Protein_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p11)


#######12. Daily Fat Intake
comps <- make_pairs(sample_data(ps1)$Fat_status)
print(comps)

kruskal.test(richness$Observed ~ sample_data(ps1)$Fat_status)
kruskal.test(richness$Shannon ~ sample_data(ps1)$Fat_status)

dunn_fat_o <- dunn.test(alpha_diversity$Observed, alpha_diversity$Fat_status, 
                         method="bonferroni")
dunn_fat_s <- dunn.test(alpha_diversity$Shannon, alpha_diversity$Fat_status, 
                         method="bonferroni")
# Print the results
print(dunn_fat_o)
print(dunn_fat_s)

sample_data(ps1)$Fat_status <- factor(sample_data(ps1)$Fat_status,
                                          levels = c("Low fat", "Average fat", "High fat"))

p12 <- plot_richness(ps1, x = "Fat_status", measures = c("Observed", "Shannon"), color = "Fat_status") +
  geom_boxplot(aes(fill = Fat_status), color = "black", outlier.color = NA) +
  scale_fill_manual(values = c("Low fat" = "#009E73", "Average fat" = "#0072B2", "High fat" = "#D55E00")) +
  labs(x = "Daily Fat Intake") +
  theme_minimal() +
  theme(strip.background = element_blank(), legend.position = "none") +
  geom_jitter(aes(color = Fat_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p12)

######## 13. Sugar intake
kruskal.test(richness$Observed ~ sample_data(ps1)$Sugars_status)
kruskal.test(richness$Shannon ~ sample_data(ps1)$Sugars_status)

dunn_sugar_o <- dunn.test(alpha_diversity$Observed, alpha_diversity$Sugars_status, 
                        method="bonferroni")
dunn_sugar_s <- dunn.test(alpha_diversity$Shannon, alpha_diversity$Sugars_status, 
                        method="bonferroni")
# Print the results
print(dunn_sugar_o)
print(dunn_sugar_s)

comps <- make_pairs(sample_data(ps1)$Sugars_status)
print(comps)
sample_data(ps1)$Sugars_status <- factor(sample_data(ps1)$Sugars_status,
                                      levels = c("Low sugar", "Average sugar", "High sugar"))

p13 <- plot_richness(ps1, x = "Sugars_status", measures = c("Observed", "Shannon"), color = "Sugars_status") +
  geom_boxplot(aes(fill = Sugars_status), color = "black", outlier.color = NA) +
  scale_fill_manual(values = c("Low sugar" = "#009E73", "Average sugar" = "#0072B2", "High sugar" = "#D55E00")) +
  labs(x = "Daily Sugar Intake") +
  theme_minimal() +
  theme(strip.background = element_blank(), legend.position = "none") +
  geom_jitter(aes(color = Sugars_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p13)

####### 14.Fiber intake
comps <- make_pairs(sample_data(ps1)$fiber_status)
print(comps)
wilcox.test(Observed ~ fiber_status, data = alpha_diversity)
wilcox.test(Shannon ~ fiber_status, data = alpha_diversity)

p14 <- plot_richness(ps1, x = "fiber_status", measures = c("Observed", "Shannon"), color = "fiber_status") +
  geom_boxplot(aes(fill = fiber_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Low Fiber Intake" = "#0072B2", "Adequate Fiber Intake" = "#D55E00")) +  # Specify actual colors
  labs(x = "General Health Status") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = fiber_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p14)
#############15. Alcohol consumption
comps <- make_pairs(sample_data(ps1)$alcohol_status)
print(comps)

wilcox.test(Observed ~ alcohol_status, data = alpha_diversity)
wilcox.test(Shannon ~ alcohol_status, data = alpha_diversity)

p15 <- plot_richness(ps1, x = "alcohol_status", measures = c("Observed", "Shannon"), color = "alcohol_status") +
  geom_boxplot(aes(fill = alcohol_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Alcohol consumers" = "#0072B2", "OccasionalNon drinkers" = "#D55E00")) +  # Specify actual colors
  labs(x = "Alcohol consumption per week") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = alcohol_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p15)

###### 16. Smoking 
comps <- make_pairs(sample_data(ps1)$smoking_status)
print(comps)

wilcox.test(Observed ~ smoking_status, data = alpha_diversity)
wilcox.test(Shannon ~ smoking_status, data = alpha_diversity)


p16 <- plot_richness(ps1, x = "smoking_status", measures = c("Observed", "Shannon"), color = "smoking_status") +
  geom_boxplot(aes(fill = smoking_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Non smoker" = "#0072B2", "Smoker" = "#D55E00")) +  # Specify actual colors
  labs(x = "Smoking status") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = smoking_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)
print(p16)


#####17. Tooth brushing frequency NS
wilcox.test(Observed ~ toothbrsh_freq_status, data = alpha_diversity)
wilcox.test(Shannon ~ toothbrsh_freq_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$toothbrsh_freq_status)
print(comps)

p17 <- plot_richness(ps1, x = "toothbrsh_freq_status", measures = c("Observed", "Shannon"), color = "toothbrsh_freq_status") +
  geom_boxplot(aes(fill = toothbrsh_freq_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Once daily" = "#0072B2", "Twice or more daily" = "#D55E00")) +  # Specify actual colors
  labs(x = "Tooth brushing frequency") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = toothbrsh_freq_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p17)

#####18. Flossing frequency NS
wilcox.test(Observed ~ flossing_status, data = alpha_diversity)
wilcox.test(Shannon ~ flossing_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$flossing_status)
print(comps)

p18 <- plot_richness(ps1, x = "flossing_status", measures = c("Observed", "Shannon"), color = "flossing_status") +
  geom_boxplot(aes(fill = flossing_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Flossing" = "#0072B2", "No Flossing" = "#D55E00")) +  # Specify actual colors
  labs(x = "Flossing frequency") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = flossing_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p18)

####### 19. Last dental visit
wilcox.test(Observed ~ last_visit_dent_status, data = alpha_diversity)
wilcox.test(Shannon ~ last_visit_dent_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$last_visit_dent_status)
print(comps)

p19 <- plot_richness(ps1, x = "last_visit_dent_status", measures = c("Observed", "Shannon"), color = "last_visit_dent_status") +
  geom_boxplot(aes(fill = last_visit_dent_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Less than 12 months" = "#0072B2", "More than 12 months" = "#D55E00")) +  # Specify actual colors
  labs(x = "Last dental visit") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = last_visit_dent_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p19)

########### 20. Fluoride application
wilcox.test(Observed ~ fluoride_prof_status, data = alpha_diversity)
wilcox.test(Shannon ~ fluoride_prof_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$fluoride_prof_status)
print(comps)

p20 <- plot_richness(ps1, x = "fluoride_prof_status", measures = c("Observed", "Shannon"), color = "fluoride_prof_status") +
  geom_boxplot(aes(fill = fluoride_prof_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Fluoride applied" = "#0072B2", "No fluoride" = "#D55E00")) +  # Specify actual colors
  labs(x = "Professional fluoride application in lifetime") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = fluoride_prof_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p20)

########### 21. Caries experience
wilcox.test(Observed ~ caries_experience_status, data = alpha_diversity)
wilcox.test(Shannon ~ caries_experience_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$caries_experience_status)
print(comps)

p21 <- plot_richness(ps1, x = "caries_experience_status", measures = c("Observed", "Shannon"), 
                     color = "caries_experience_status") +
  geom_boxplot(aes(fill = caries_experience_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Caries Experience" = "#0072B2", "No Caries Experience" = "#D55E00")) +  # Specify actual colors
  labs(x = "Caries experience in lifetime") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = caries_experience_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p21)

########### 22.Saliva flow rate
wilcox.test(Observed ~ saliva_flow_status, data = alpha_diversity)
wilcox.test(Shannon ~ saliva_flow_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$saliva_flow_status)
print(comps)

p22 <- plot_richness(ps1, x = "saliva_flow_status", measures = c("Observed", "Shannon"), 
                     color = "saliva_flow_status") +
  geom_boxplot(aes(fill = saliva_flow_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Low saliva flow rate" = "#0072B2", 
                               "Normal saliva flow rate" = "#D55E00")) +  # Specify actual colors
  labs(x = "Saliva flow rate") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = saliva_flow_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


print(p22)

########### 23.Saliva pH
comps <- make_pairs(sample_data(ps1)$ph_status)
print(comps)

kruskal.test(richness$Observed ~ sample_data(ps1)$ph_status)
kruskal.test(richness$Shannon ~ sample_data(ps1)$ph_status)

dunn_ph_o <- dunn.test(alpha_diversity$Observed, alpha_diversity$ph_status, 
                          method="bonferroni")
dunn_ph_s <- dunn.test(alpha_diversity$Shannon, alpha_diversity$ph_status, 
                          method="bonferroni")
# Print the results
print(dunn_ph_o)
print(dunn_ph_s)

sample_data(ps1)$ph_status<- factor(sample_data(ps1)$ph_status,
                                         levels = c("Low pH", "Average pH", "High pH"))

p23 <- plot_richness(ps1, x = "ph_status", measures = c("Observed", "Shannon"), color = "ph_status") +
  geom_boxplot(aes(fill = ph_status), color = "black", outlier.color = NA) +
  scale_fill_manual(values = c("Low pH" = "#009E73", "Average pH" = "#0072B2", "High pH" = "#D55E00")) +
  labs(x = "ph_status") +
  theme_minimal() +
  theme(strip.background = element_blank(), legend.position = "none") +
  geom_jitter(aes(color = ph_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)



print(p23)

########### 24.Saliva glucose test
wilcox.test(Observed ~ postglucose_status, data = alpha_diversity)
wilcox.test(Shannon ~ postglucose_status, data = alpha_diversity)

comps <- make_pairs(sample_data(ps1)$postglucose_status)
print(comps)

p24 <- plot_richness(ps1, x = "postglucose_status", measures = c("Observed", "Shannon"), 
                     color = "postglucose_status") +
  geom_boxplot(aes(fill = postglucose_status), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("Low pH" = "#0072B2", 
                               "High pH" = "#D55E00")) +  # Specify actual colors
  labs(x = "pH Post Glucose Challenge test") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = postglucose_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p24)

#########Combining the alpha diversity plot that were significant
library(patchwork)
p_alpha <- p3+ p5+p6+p10+p12+p13+p14+p19+p20 + plot_layout(ncol = 3)

# Save the plot as PNG
ggsave(filename = "alpha_diversity_sig.png", plot = p_alpha, width = 20, height = 10, dpi = 300)

# Save the plot as SVG
ggsave(filename = "alpha_diversity_sig.svg", plot = p_alpha, width = 20, height = 10)
########

##combining the alpha diversity plots for non-significant factors
#p1+p4+p5+p6+p7+p9+p12+p13+p15+p16+ plot_layout(ncol = 2)


##########################################
#############Beta diversity#########
###########################################
library(microbiome)
library(ggplot2)
library(dplyr)
library(vegan)

metadata <- data.frame(sample_data(ps1))
# Calculate both Bray-Curtis and UniFrac distances
bray_dist <- distance(ps1, method = "bray")
unifrac_dist <- distance(ps1, method = "unifrac", weighted = TRUE)

####### 1. Age#######
adonis2(bray_dist ~ age_group, data = metadata, permutations = 999)
adonis2(unifrac_dist ~ age_group,
                      data = metadata, permutations=999)

####### 2. Gender
adonis2(bray_dist ~ gender.factor.x,
                            data = metadata, permutations=999)
adonis2(unifrac_dist ~ gender.factor.x,
                           data = metadata, permutations=999)

#######3. Country born
adonis2(bray_dist ~ country_born.factor,
                    data = metadata, permutations=999)
adonis2(unifrac_dist ~ country_born.factor,
                             data = metadata, permutations=999)


#######4.Level of education
adonis2(bray_dist ~ education_status,
                             data = metadata, permutations=999)
adonis2(unifrac_dist ~ education_status,
                             data = metadata, permutations=999)

#######5. Study status
adonis2(bray_dist~ study_status,
                               data = metadata, permutations=999,)
adonis2(unifrac_dist ~ study_status,
                               data = metadata, permutations=999)


#######6. Employment
adonis2(bray_dist ~ employment_status,
                             data = metadata, permutations=999,)
adonis2(unifrac_dist ~ employment_status,
                             data = metadata, permutations=999)

#######7. Physical activity
adonis2(bray_dist ~ physical_activity_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ physical_activity_status,
        data = metadata, permutations=999)


#######8. Mental Well-being
adonis2(bray_dist ~ mental_health_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ mental_health_status,
        data = metadata, permutations=999)

######9. Self rated general health
adonis2(bray_dist ~ generalhealth_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ generalhealth_status,
        data = metadata, permutations=999)

###### 10. Carbohydrate intake
adonis2(bray_dist ~ carbohydrate_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ carbohydrate_status,
        data = metadata, permutations=999)

###### 11. Protein intake
adonis2(bray_dist ~ Protein_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ Protein_status,
        data = metadata, permutations=999)

###### 12. Fat intake
adonis2(bray_dist ~ Fat_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ Fat_status,
        data = metadata, permutations=999)

###### 13. Sugar intake
adonis2(bray_dist ~ Sugars_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ Sugars_status,
        data = metadata, permutations=999)


###### 14. Fiber intake
adonis2(bray_dist ~ fiber_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ fiber_status,
        data = metadata, permutations=999)


###### 15. Alcohol consumption
adonis2(bray_dist ~ alcohol_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ alcohol_status,
        data = metadata, permutations=999)

###### 16. Smoking
adonis2(bray_dist ~ smoking_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ smoking_status,
        data = metadata, permutations=999)


###### 17. Tooth brushing frequency
adonis2(bray_dist ~ toothbrsh_freq_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ smoking_status,
        data = metadata, permutations=999)


###### 18. Flossing frequency
adonis2(bray_dist ~ flossing_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ flossing_status,
        data = metadata, permutations=999)


###### 19. Last dental visit
adonis2(bray_dist ~ last_visit_dent_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ last_visit_dent_status,
        data = metadata, permutations=999)



###### 20. Fluoride application
adonis2(bray_dist ~ fluoride_prof_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ fluoride_prof_status,
        data = metadata, permutations=999)

###### 21. Caries experience
adonis2(bray_dist ~ caries_experience_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ caries_experience_status,
        data = metadata, permutations=999)


###### 22. Saliva flow rate
adonis2(bray_dist ~ saliva_flow_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ saliva_flow_status,
        data = metadata, permutations=999)


###### 23. Saliva pH
adonis2(bray_dist ~ ph_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ ph_status,
        data = metadata, permutations=999)


###### 24. Post glucose saliva test
adonis2(bray_dist~ postglucose_status,
        data = metadata, permutations=999)
adonis2(unifrac_dist ~ postglucose_status,
        data = metadata, permutations=999)


##########Bray Curtis######
adonis_combine_b <- adonis2(bray_dist ~ country_born.factor + study_status * employment_status + 
                            generalhealth_status + carbohydrate_status * Sugars_status + 
        fiber_status + flossing_status + last_visit_dent_status * fluoride_prof_status + ph_status *
          postglucose_status,
        data = metadata, permutations=999)
write.table(adonis_combine_b,
            file = "adonis_combined_b1.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

adonis_combine_b1 <- adonis2(bray_dist ~ country_born.factor + study_status * employment_status + 
                               carbohydrate_status  + 
                              last_visit_dent_status * fluoride_prof_status + ph_status *
                              postglucose_status,
                            data = metadata, permutations=999)
write.table(adonis_combine_b1,
            file = "adonis_combined_b2.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
adonis_combine_b2 <- adonis2(bray_dist ~ country_born.factor + study_status + employment_status + 
                                last_visit_dent_status * fluoride_prof_status + ph_status,
                             data = metadata, permutations=999)

write.table(adonis_combine_b2,
            file = "adonis_combined_b3.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#####################
##########Unifrac######
adonis_combine_u <-  adonis2(unifrac_dist ~ country_born.factor + study_status * employment_status + 
                               generalhealth_status + carbohydrate_status * Sugars_status + 
                               fiber_status + flossing_status + last_visit_dent_status * fluoride_prof_status 
                             + ph_status + postglucose_status,
                             data = metadata, permutations=999)
write.table(adonis_combine_u,
            file = "adonis_combined_u1.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

adonis_combine_u1 <-  adonis2(unifrac_dist ~ country_born.factor + study_status +  
                               generalhealth_status + carbohydrate_status +
                               flossing_status + last_visit_dent_status * fluoride_prof_status,
                             data = metadata, permutations=999)
write.table(adonis_combine_u1,
            file = "adonis_combined_u2.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

adonis_combine_u2 <-  adonis2(unifrac_dist ~ country_born.factor +
                                generalhealth_status + carbohydrate_status +
                                flossing_status + last_visit_dent_status,
                              data = metadata, permutations=999)
write.table(adonis_combine_u2,
            file = "adonis_combined_u3.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
########ANOVA########
anova_result_shannon <- aov(Shannon ~ country_born.factor + study_status +
                              employment_status + carbohydrate_status +
                              Fat_status + Sugars_status + fiber_status +
                              last_visit_dent_status + fluoride_prof_status,
                            data = alpha_diversity)
anova_result_shannon <- summary(anova_result_shannon)
anova_df_shannon <- as.data.frame(anova_result_observed[[1]])
write.table(anova_df_shannon,
            file = "anova_shannon.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

anova_result_observed <- aov(Observed ~ country_born.factor + study_status +
                               employment_status + carbohydrate_status +
                               Fat_status + Sugars_status + fiber_status +
                               last_visit_dent_status + fluoride_prof_status,
                             data = alpha_diversity)

anova_result_observed <- summary(anova_result_observed)
anova_df_observed <- as.data.frame(anova_result_observed[[1]])

write.table(anova_df_observed,
            file = "anova_observed.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Visualize with PCoA for one of the distances, e.g., Bray-Curtis
#ordi_bray <- ordinate(ps1, method = "PCoA", distance = bray_dist)
#plot_ordination(ps1, ordi_bray, color = "country_born.factor") + geom_point(aes(shape = country_born.factor)) + theme_minimal()
######################

##########Correlation analysis and RDA plots##############

#########################
set.seed(1) # for reproducible stochastic processes
library(phyloseq)
library(ggplot2)
library(patchwork) # for combining multiple plots
library(microViz)

##checking the data

sample_names(ps1) %>% head()
taxa_names(ps1) %>% head()
samdat_tbl(ps1)
otu_get(ps1, taxa = 1:3, samples = 1:5) 
rank_names(ps1)

##Fixing the tables
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
phyloseq_validate(ps3)

ps3 <- ps3 %>% tax_fix(unknowns = c("uncultured"))

##############

library(microViz)
library(corncob)
library(dplyr)
library(ggplot2)
phylo <- ps3 %>%
  tax_mutate(Species = NULL) %>%
  tax_fix()
phylo



#########correlation
corr_plot <- ps3 %>%
  tax_agg("Genus") %>%
  tax_sort(by = prev, at = "Genus") %>%
  cor_heatmap(
    seriation_method = "Identity",
    seriation_method_col = "OLO_ward",
    taxa = tax_top(ps3, 20, by = max, rank = "Genus"),
    vars = c( "age", "gender", "country_born", "emp_stat", "filled", "Sugars",
             "Protein","Fat", "Carbohydrate", "Dietary.Fibre", "k10score", "total_met_minutes", "baseline_pH",
             "post.glucose_2"),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:1),
      CLR = anno_tax_box(trans = "clr", zero_replace = "halfmin")
    )
  )


R1<- ps3 %>%
  tax_agg("Genus") %>%
  tax_sort(by = prev, at = "Genus") %>%
  cor_heatmap(
    seriation_method = "Identity",
    seriation_method_col = "OLO_ward",
    taxa = tax_top(ps3, 15, by = max, rank = "Genus"),
    vars = c("age", "gender", "country_born","emp_stat", "filled"),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:2),
      CLR = anno_tax_box(trans = "clr", zero_replace = "halfmin")
    )
  )
print(R1)

############RDA plot
library(stringr)
library(phyloseq)
library(ggplot2)
library(microViz)
renamer <- function(x) str_replace(x, pattern = "_", replacement = " ")# first we make a function that replaces any unwanted "_" in our taxa labels with spaceslibrary(stringr)

ps3 %>%
  ps_mutate(
    Older_adult = ifelse(age_group == "Adults", yes = 1, no = 0),
    Male = ifelse(gender.factor.x == "Male", yes = 1, no = 0),  
    Australian = ifelse(country_born.factor == "Australia", yes = 1, no = 0),
    Student = ifelse(study_status == "Full time student ", yes = 1, no = 0),
    Employed = ifelse(employment_status == "Employed", yes = 1, no = 0),
    High_carb = ifelse(carbohydrate_status == "High carb", yes = 1, no = 0),
    High_sugar = ifelse(Sugars_status == "High sugar", yes = 1, no = 0),
    High_fat= ifelse(Fat_status == "High fat",yes = 1, no = 0 ),
    Low_fibre = ifelse(fiber_status == "Low fiber intake", yes = 1, no = 0),
    Regular_dental_check = ifelse(last_visit_dent_status == "Less than 12 months", yes = 1, no = 0),
    Fluoride_application = ifelse(fluoride_prof_status == "Fluoride applied", yes = 1, no = 0),
    Twice_brushing = ifelse(toothbrsh_freq_status == "Twice or more daily", yes = 1, no = 0),
    Flossing = ifelse(flossing_status == "Flossing", yes = 1, no = 0),
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("Older_adult", "Male", "Australian","Student", "Employed", "High_carb","High_sugar", "High_fat", 
                    "Low_fibre", "Regular_dental_check", "Fluoride_application", "Twice_brushing", "Flossing"),
    method = "RDA",
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "age_group", size = 4, alpha = 0.5,shape = "age_group",
    auto_caption = NA, # remove the helpful automatic caption
    plot_taxa = 1:25, taxon_renamer = renamer, # renamer is the function we made earlier
    tax_vec_length = 5, # this value is a scalar multiplier for the biplot score vectors
    tax_lab_length = 6, # this multiplier moves the labels, independently of the arrowheads
    tax_lab_style = tax_lab_style(size = 1.8, alpha = 0.5), # create a list of options to tweak the taxa labels' default style
    constraint_vec_length = 4, # this adjusts the length of the constraint arrows, and the labels track these lengths by default
    constraint_vec_style = vec_constraint(2, alpha = 0.5), # this styles the constraint arrows
    constraint_lab_style = constraint_lab_style(size = 4) # this styles the constraint labels
  ) +
  scale_colour_brewer(palette = "Set1")


######2
# Load necessary libraries
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)

# Apply the transformations and perform RDA
ps_transformed <- ps3 %>%
  ps_mutate(
    Older_adult = ifelse(age_group == "Adults", yes = 1, no = 0),
    Male = ifelse(gender.factor.x == "Male", yes = 1, no = 0),  
    Australian = ifelse(country_born.factor == "Australia", yes = 1, no = 0),
    Tertiary_education = ifelse(education_status == "Tertiary education", yes = 1, no = 0),
    Low_physical_activity = ifelse(physical_activity_status == "Minimal activity", yes = 1, no = 0),
    Depression_anxiety = ifelse(mental_health_status == "Severe distress", yes = 1, no = 0),
    Student = ifelse(study_status == "Full time student", yes = 1, no = 0),
    Employed = ifelse(employment_status == "Employed", yes = 1, no = 0),
    High_carb = ifelse(carbohydrate_status == "High carb", yes = 1, no = 0),
    High_sugar = ifelse(Sugars_status == "High sugar", yes = 1, no = 0),
    Low_protein = ifelse(Protein_status == "Low protein", yes = 1, no = 0),
    High_fat = ifelse(Fat_status == "High fat", yes = 1, no = 0),
    Low_fibre = ifelse(fiber_status == "Low fiber intake", yes = 1, no = 0),
    Alcohol_consumption = ifelse(alcohol_status == "Alcohol consumers", yes = 1, no = 0),
    Smokers = ifelse(smoking_status == "Smoker", yes = 1, no = 0),
    Regular_dental_check = ifelse(last_visit_dent_status == "Less than 12 months", yes = 1, no = 0),
    Fluoride_application = ifelse(fluoride_prof_status == "Fluoride applied", yes = 1, no = 0),
    Twice_brushing = ifelse(toothbrsh_freq_status == "Twice or more daily", yes = 1, no = 0),
    Caries_experience = ifelse(caries_experience_status == "Caries experience", yes = 1, no = 0),
    Flossing = ifelse(flossing_status == "Flossing", yes = 1, no = 0),
    Low_saliva_flowrate = ifelse(saliva_flow_status == "Low saliva flow rate", yes = 1, no = 0),
    Low_saliva_pH = ifelse(ph_status == "Low pH", yes = 1, no = 0),
    Low_saliva_ph_glucose = ifelse(postglucose_status == "Low pH", yes = 1, no = 0)
  ) %>%
  tax_transform("clr", rank = "Genus")

# Print the column names to ensure the variables were created correctly
print(colnames(sample_data(ps_transformed)))

rda_result <- ps_transformed %>%
  ord_calc(
    constraints = c("Older_adult", "Male", "Australian", "Tertiary_education", "Low_physical_activity",
                    "Depression_anxiety", "Student", "Employed", "High_carb", "High_sugar", "Low_protein", "High_fat", 
                    "Low_fibre", "Alcohol_consumption", "Smokers", "Regular_dental_check", "Fluoride_application", 
                    "Twice_brushing", "Caries_experience", "Flossing",
                    "Low_saliva_flowrate", "Low_saliva_pH", "Low_saliva_ph_glucose"),
    method = "RDA",
    scale_cc = FALSE
  ) %>%
  ord_plot(
    colour = "age_group", size = 4, alpha = 0.5, shape = "age_group",
    auto_caption = NA, 
    plot_taxa = 1:30, taxon_renamer = renamer, 
    tax_vec_length = 6, 
    tax_lab_length = 8, 
    tax_lab_style = tax_lab_style(size = 1.8, alpha = 0.5),
    constraint_vec_length = 12, 
    constraint_vec_style = vec_constraint(2, alpha = 0.5),
    constraint_lab_style = constraint_lab_style(size = 4)
  ) +
  scale_colour_brewer(palette = "Set1")

# Save the plot as SVG
ggsave(filename = "rda_plot.svg", plot = rda_result, device = "svg", width = 10, height = 8)

# Print the RDA plot to verify
print(rda_result)
####################################
##################HeatMap ########
####################################
library(dplyr)
# Transform the data to compositional and then to CLR
ps_transformed <- ps3 %>%
  tax_transform("compositional", rank = "Genus") %>%
  tax_transform("clr", zero_replace = "halfmin", chain = TRUE)

# Calculate the total abundance of each genus and select the top 9
# Extract the abundance table and calculate total abundance for each genus
abundance_table <- otu_table(ps_transformed)
tax_table_df <- as.data.frame(tax_table(ps_transformed))

# Calculate total abundance per genus
total_abundance <- rowSums(abundance_table)
tax_table_df$total_abundance <- total_abundance

# Select the top 30 genera based on total abundance
top_genera <- tax_table_df %>%
  group_by(Genus) %>%
  summarize(total_abundance = sum(total_abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 50) %>%
  pull(Genus)

# Prune the phyloseq object to keep only the top 30 genera
ps_filtered <- prune_taxa(tax_table(ps_transformed)[, "Genus"] %in% top_genera, ps_transformed)

# Create the heatmap with the filtered data
heatmap_plot <- comp_heatmap(ps_filtered, colors = heat_palette(sym = TRUE), name = "Abundance CLR",
             tax_anno = taxAnnotation(
               Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
             ),
             heatmap_legend_param = list(
               at = log10(c(1, 0.1, 0.01, 0.001, 0.0001)),
               labels = c("100%", "10%", "1%", "0.1%", "0.01%")
             )
)
leg_plot <- comp_heatmap(ps_filtered, colors = heat_palette(sym = TRUE), name = "Abundance CLR",
                             tax_anno = taxAnnotation(
                               Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
                             )
)

ggsave(filename = "leg_plot.svg", plot = leg_plot, device = "svg", width = 10, height = 20)

# Save the plot as PNG
ggsave(filename = "heatmap_plot.png", plot = heatmap_plot, width = 10, height = 20, dpi = 300)

# Save the plot as SVG
ggsave(filename = "heatmap_plot.svg", plot = heatmap_plot, device = "svg", width = 10, height = 20)

##############ABDUNDNACE PLOT###########
# set up for alphabetical sorting
topTaxa <- ps3 %>%
  tax_top(n = 9, rank = "Genus") %>%
  sort() # this makes them alphabetical

# plot with alphabetical sorting
taxa_plot <- ps3 %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  comp_barplot(
    tax_order = topTaxa, # this brings the named taxa to the front
    tax_level = "Genus", n_taxa = 9, merge_other = FALSE, other_name = "Other",
    label = "pid",bar_outline_colour = "grey5"
  ) +
  coord_flip()

# Save the plot as PNG
ggsave(filename = "taxa_plot.png", plot = taxa_plot, width = 20, height = 20, dpi = 300)

# Save the plot as SVG
ggsave(filename = "taxa_plot.svg", plot = taxa_plot, device = "svg", width = 10, height = 20)

#############
library(rsvg)
library(ggplot2)
library(grid)

# Load the SVG file
svg_file <- "heatmap.svg"

# Render the SVG file as a raster image
svg_image <- rsvg::rsvg(svg_file)

# Convert to a grid object
svg_grob <- rasterGrob(svg_image, width = unit(1, "npc"), height = unit(1, "npc"))
svg_plot <- ggplot() +
  annotation_custom(svg_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()
##################
# Combine the plots with labels

library(patchwork)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)
library(ggplot2)

heatmap_grob <- grid::grid.grabExpr(draw(heatmap_plot))

# Convert the grob to a ggplot-compatible object using cowplot
heatmap_gg <- ggdraw() + draw_grob(heatmap_grob)

Figure1 <- (heatmap_gg | p_alpha) + 
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = 'a')
 

# Save the combined plot as an image file
ggsave(filename = "Figure1.png", plot = Figure1, width = 25, height = 10, dpi = 300)
ggsave(filename = "Figure1.svg", plot = Figure1, width = 25, height = 10)
