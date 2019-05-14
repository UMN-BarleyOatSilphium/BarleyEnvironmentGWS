## S2MET Prediction Models source script
## 
## A script that automatically loads the data relevant for the S2MET project

library(tidyverse)
library(readxl)
library(rrBLUP)
library(neyhart)
library(boot)
library(pbr)
library(lmerTest)

## Directories
proj_dir <- repo_dir
alt_proj_dir <- "C:/Users/jln54/GoogleDrive/BarleyLab/Projects/S2MET_Predictions//"

## Google drive directory
gdrive_dir <- "C:/Users/jln54//GoogleDrive"

# Geno, pheno, and enviro data
geno_dir <-  file.path(gdrive_dir, "BarleyLab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/")

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- pheno_dir <- file.path(alt_proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")



######
# MSI Source starts here
######


# Source the project functions
source(file.path(proj_dir, "source_functions.R"))


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))


# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Load the trial metadata
trial_info <- read_csv(file.path(data_dir, "trial_metadata.csv"))

# Vector of relevant traits
traits <- c("GrainYield", "HeadingDate", "PlantHeight")


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

vp <- entry_list %>% 
  filter(Class == "S2C1R") %>% 
  pull(Line)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  filter(Class != "Check") %>% 
  pull(Line)

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Calculate the K matrix
K <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)

# # Create a replacement vector for each trait that adds a space between words and
# # also included the unit
# trait_replace <- unique(S2_MET_BLUEs$trait) %>%
#   set_names(x = c("Grain Yield", "Heading Date", "Plant Height"), nm = .)
# 
# trait_replace_unit <- unique(S2_MET_BLUEs$trait) %>%
#   set_names(x = c("Grain Yield\n(kg ha^-1)", "Heading Date\n(days)", "Plant Height\n(cm)"), nm = .)
# 
# # Create a name replacement vector for the distance methods
# dist_method_replace <- c(
#   "ge_mean_D" = "Phenotypic\nDistance", 
#   "ge_PCA_dist" = "GxE BLUP PCA", 
#   "great_circle_dist" = "Great Circle\nDistance", 
#   # "ec_one_PCA_dist" = "1 yr Environmental\nCovariates", 
#   "ec_multi_PCA_dist" = "10 yr Environmental\nCovariates")


## Rank the environments according to heritability
env_herit_rank <- stage_one_data %>% 
  ungroup() %>% 
  filter(!str_detect(trial, "S2C1")) %>% 
  select(trait, environment, heritability) %>% 
  split(.$trait) %>%
  map(~arrange(., desc(heritability)) %>% 
        mutate(environment = factor(environment, levels = .$environment)))


# Pull out the trials with irrigation
irrig_env <- subset(trial_info, irrigated == "yes", environment)[[1]]

# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs_temp <- S2_MET_BLUEs %>% 
  filter(!environment %in% irrig_env) %>%
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "BCW16" & trait == "PlantHeight"))

## Remove environments with low heritability
high_herit_environments <- env_herit_rank %>%
  map(mutate, environment = as.character(environment)) %>% 
  map(~filter(., heritability >= 0.10))

S2_MET_BLUEs_temp1 <- S2_MET_BLUEs_temp %>% 
  split(.$trait) %>%
  map2_df(.x = ., .y = high_herit_environments, ~filter(.x, environment %in% .y$environment))


## Reassign BLUEs
S2_MET_BLUEs <- S2_MET_BLUEs_temp1


# Find environments in which just data on both the TP and VP is available
tp_vp_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()

# Find environments in which just data on the TP is available
tp_only_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) == 0) %>% 
  distinct(environment) %>% 
  pull()

# Find environments in which just data on the VP is available
vp_only_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) == 0, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()

## Split these vectors based on traits
tp_vp_env_trait <- S2_MET_BLUEs %>% 
  group_by(environment, trait) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) > 1) %>% 
  split(.$trait) %>% 
  map(~unique(.$environment))
  
tp_only_env_trait <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) == 0) %>%
  split(.$trait) %>% 
  map(~unique(.$environment))
  
vp_only_env_trait <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) == 0, sum(line_name %in% vp_geno) > 1) %>% 
  split(.$trait) %>% 
  map(~unique(.$environment))


## Create two list of environments - for complete or realistic training
complete_train_env <- tp_vp_env_trait

realistic_train_env <- tp_vp_env_trait %>%
  map(~str_subset(., "15|16"))

realistic_test_env <- tp_vp_env_trait %>%
  map(~str_subset(., "17"))


## Final filter of BLUEs
S2_MET_BLUEs <- filter(S2_MET_BLUEs, environment %in% tp_vp_env)


## Remove
rm(S2_MET_BLUEs_temp, S2_MET_BLUEs_temp1, s2_discrete_mat, s2_imputed_mat, s2_imputed_mat_use)


### Colors, abbreviations, etc.

## Significant level
alpha <- 0.05

## Create a factor for the distance methods
dist_method_replace <- c("pheno_dist" = "Phenotypic Distance", "pheno_loc_dist" = "Location Phenotypic Distance",  "great_circle_dist" = "Great Circle Distance", 
                         "OYEC_All" = "One Year All ECs", "OYEC_Mean" = "One Year Mean Cor EC", "OYEC_IPCA" = "One Year IPCA Cor EC", 
                         "MYEC_All" = "Multi Year All ECs", "MYEC_Mean" = "Multi Year Mean Cor EC", "MYEC_IPCA" = "Multi Year IPCA Cor EC",
                         "AMMI" = "AMMI", "sample" = "Random")
dist_method_abbr <- abbreviate(dist_method_replace)

## Alternative abbreviations
dist_method_abbr <- setNames(c("PD", "LocPD", "GCD", "1Yr-All-EC", "1Yr-Mean-EC", "1Yr-IPCA-EC", "All-EC", "Mean-EC", "IPCA-EC", "AMMI", "Random"),
                             names(dist_method_replace))

colors <- umn_palette(3)
colors2 <- umn_palette(2)
colors3 <- umn_palette(4)

colors_use <- c("#FFB71E", "#FFDE7A", colors[1], colors[5:7], colors[c(3, 8)], colors2[3], colors2[4], "grey75")

# colors <- wesanderson::wes_palette(name = "Darjeeling1", n = 9, type = "continuous")
# colors2 <- wesanderson::wes_palette(name = "Darjeeling2", n = 10, type = "continuous")
# colors_use <- c(colors[c(1:3, 5:7)], colors2[2:4], colors[4], "grey75")

dist_colors <- setNames(colors_use, dist_method_abbr)

## Subset for use
dist_method_abbr_use <- dist_method_abbr[c("AMMI", "pheno_dist", "pheno_loc_dist", "great_circle_dist", "MYEC_All", "MYEC_Mean", "MYEC_IPCA", "sample")]
dist_colors_use <- dist_colors[dist_method_abbr_use]



# ## Replacement vector for CV
cv_replace <- c("cv1", "pov1", "pocv1",  "cv2" , "pocv2", "cv0", "pov0", "pocv0", "cv00", "pov00", "pocv00") %>%
  setNames(object = toupper(.), .)

## Set replacement
set_replace <- c("complete" = "Leave-one-out", "realistic" = "Time-forward")

