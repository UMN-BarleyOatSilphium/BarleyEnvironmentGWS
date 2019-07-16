## S2MET Prediction Models source script
## 
## A script that automatically loads the data relevant for the S2MET project

library(sommer)
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
geno_dir <-  file.path(gdrive_dir, "BarleyLab/Breeding/GenotypicData/GBS_Genotype_Data/")

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
load(file.path(pheno_dir, "S2_MET_BLUEs_models.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))


# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Load the trial metadata
trial_info <- read_csv(file.path(data_dir, "trial_metadata.csv"))

# Vector of relevant traits
# traits <- c("GrainYield", "HeadingDate", "PlantHeight", "MaturityDate")
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




# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs_temp <- S2_MET_BLUEs %>% 
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



## Final filter of BLUEs
S2_MET_BLUEs <- filter(S2_MET_BLUEs, environment %in% tp_vp_env)

## Create a df with trait, environment, heritability
env_trait_herit <- stage_one_data %>% 
  filter(!str_detect(trial, "S2C1"),
         environment %in% tp_vp_env,
         trait %in% traits) %>% 
  distinct(trait, environment, heritability)

## Remove
rm(S2_MET_BLUEs_temp, S2_MET_BLUEs_temp1, s2_discrete_mat, s2_imputed_mat, s2_imputed_mat_use, env_herit_rank, stage_one_data)


### Colors, abbreviations, etc.


colors <- umn_palette(3)
colors2 <- umn_palette(2)
colors3 <- umn_palette(4)

colors_use <- c("#FFB71E", "#FFDE7A", colors[1], colors[5:7], colors[c(3, 8)], colors2[3], colors2[4], "grey75")


# ## Replacement vector for CV
cv_replace <- c("cv1", "pov1", "pocv1",  "cv2" , "pocv2", "cv0", "pov0", "pocv0", "cv00", "pov00", "pocv00") %>%
  setNames(object = toupper(.), .)


## Create codes for the scheme numbers (designate genotype/environment tested/untested)
scheme_number_replace <- c("00" = "G^u~E^u", "0" = "G^t~E^u", "1" = "G^u~E^t", "2" = "G^t~E^t")







## Sample environments common to all traits
## Sample data
set.seed(1227)
sample_env <- S2_MET_BLUEs %>%
  filter(trait %in% traits) %>%
  distinct(trait, environment) %>%
  group_by(environment) %>%
  filter(n() == length(traits)) %>%
  distinct(environment) %>%
  pull() %>%
  sample(5)




