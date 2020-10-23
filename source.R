## S2MET Prediction Models source script
## 
## A script that automatically loads the data relevant for the S2MET project


# List of packages
pkgs <- c("sommer", "tidyverse", "readxl", "rrBLUP", "neyhart", "pbr")
# Load these packages
invisible(lapply(X = pkgs, library, character.only = TRUE))


## Directories
proj_dir <- repo_dir

## Google drive directory
gdrive_dir <- "C:/GoogleDrive"

# Geno, pheno, and enviro data
geno_dir <-  file.path(gdrive_dir, "BarleyLab/Breeding/GenotypicData/GBS_Genotype_Data/")
pheno_dir <- file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/")
meta_dir <- file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Metadata/")
enviro_dir <- file.path(gdrive_dir, "BarleyLab/Breeding/EnvironmentalData/")

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


# Load variable csv
## Table of variables, nicknames, and units for HWSD
hwsd_variables <- read.csv(file = "C:/GoogleDrive/BarleyLab/Breeding/EnvironmentalData/RawData/SoilData/HWSD/hwsd_variable_reference.csv", stringsAsFactors = FALSE)





######
# MSI Source starts here
######


# Source the project functions
source(file.path(proj_dir, "source_functions.R"))


# Load the phenotypic data
load(file.path(pheno_dir, "S2_tidy_BLUE.RData"))

# Load the trial metadata
trial_info <- read_csv(file = file.path(data_dir, "trial_metadata.csv"))


# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))


# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Grab the entry names that are not checks
tp <- subset(entry_list, Class == "S2TP", Line, drop = T)
tp <- setdiff(tp, c("07MT-10")) # Remove hullness line
vp <- subset(entry_list, Class == "S2C1R", Line, drop = T)


# Vector of relevant traits
# traits <- c("GrainYield", "HeadingDate", "PlantHeight", "MaturityDate")
traits <- c("GrainYield", "HeadingDate", "PlantHeight", "TestWeight", "GrainProtein")
trials <- subset(trial_info, project2 == "S2MET" & str_detect(trial, "S2C1", negate = TRUE), trial, drop = TRUE)


# ## Traits by environment
# s2_tidy_BLUE %>% 
#   filter(trait %in% traits, environment %in% trial_info$environment, year > 2014) %>% 
#   mutate(nEnvAll = n_distinct(environment)) %>%
#   group_by(trait) %>% 
#   summarize(nEnv = n_distinct(environment), nEnvAll = mean(nEnvAll))


# Find the tp and vp that are genotyped with markers
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))


# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Calculate the K matrix
K <- Kgeno <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)

## Add missing entries as unrelated
K <- bdiag(K, diag(x = mean(diag(K)), nrow = length(c(tp, vp)) - length(c(tp_geno, vp_geno))))
K <- as.matrix(K); dimnames(K) <- replicate(2, c(row.names(s2_imputed_mat_use), setdiff(c(tp, vp), c(tp_geno, vp_geno))), simplify = FALSE)
K <- K[sort(row.names(K)), sort(row.names(K))]




## Rank the environments according to heritability
env_trait_herit <- s2_metadata %>% 
  select(trial, trait, heritability, varR) %>% 
  inner_join(., distinct(trial_info, trial, environment), by = "trial") %>%
  filter(trial %in% trials, trait %in% traits, heritability >= 0.10) %>%
  select(trait, environment, varR, heritability)

## 3 trait-trials were removed


## Filter the S2 tidy blues for S2MET
S2_MET_BLUEs <- s2_tidy_BLUE %>%
  filter(trial %in% subset(trial_info, project2 == "S2MET", trial, drop = T),
         line_name %in% c(tp, vp),
         trait %in% traits) %>%
  # Add full location names and rename environments according to trial_info
  select(-location, -environment) %>%
  inner_join(., distinct(trial_info, trial, environment, location), by = "trial") %>% # 150 env-traits here
  # Filter out environments with low heritability
  inner_join(., select(env_trait_herit, trait, environment), by = c("trait", "environment")) %>% # 147 env-traits here
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM")) %>% # 138 env-traits here
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "EON17" & trait == "HeadingDate"),
         !(environment == "KNY16" & trait == "TestWeight"),
         !(location %in% c("Charlottetown", "Alburgh", "Grande_rhonde_valley") & trait == "HeadingDate"),
         !(trait %in% c("PlantHeight", "TestWeight") & location %in% c("Grande_rhonde_valley")) ) %>%
  # Rename and reorder
  select(trial, environment, location, year, trait, line_name, value, std_error = std.error) # 130 env-traits here

# 17 more trait-trials removed




## Separate environments into those for training/testing and those for external validation
train_test_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull() %>%
  sort()

validation_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) == 0, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()

# Translate these to locations
train_test_loc <- sort(unique(subset(trial_info, environment %in% train_test_env, location, drop = TRUE)))
validation_loc <- sort(unique(subset(trial_info, environment %in% validation_env, location, drop = TRUE)))



## Final filter of BLUEs
S2_MET_BLUEs <- filter(S2_MET_BLUEs, environment %in% c(train_test_env, validation_env)) %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca")) %>%
  arrange(trait, environment)



## Trait units - renaming vector
trait_units <- setNames(object = c("kg~ha^-1", "days", "cm", "g~L^-1", "'%'"), nm = traits)
trait_units1 <- setNames(object = c("Mg~ha^-1", "days", "cm", "g~L^-1", "'%'"), nm = traits)

# The favorable sign for each trait
trait_sign <- tibble(trait = traits, sign = c(1, -1, -1, 1, -1))



## Functions that might be useful for plotting ##


## A vector to rename models
model_replace <- c("model1" = "g", "model2_id" = "g + e", "model2_cov" = "g + e", 
                   "model3_id" = "g + e + (ge)", "model3_cov" = "g + e + (ge)", "model3_cov1" = "g + e + (ge)",
                   "model4_id" = "g + l", "model4_cov" = "g + l", "model5_id" = "g + l + (gl)", 
                   "model5_cov" = "g + l + (gl)", "model5_cov1" = "g + l + (gl)")


f_model_replace <- function(x) model_replace[x]
# f_model_replace <- function(x) paste0("M", toupper(str_extract(x, "[0-9]{1}[a-z]{0,1}")))
# Vector to rename validation schemes
f_validation_replace <- function(x) str_replace_all(x, c("tp" = "Tested founders", "vp" = "Untested offspring"))
f_pop_replace <- function(x) str_replace_all(x, c("all" = "All", "tp" = "FP", "vp" = "OP"))
# Replace type
f_type_replace <- function(x) c("loeo" = "New environment", "lolo" = "New location", "loyo" = "New year",
                                "env_external" = "Holdout environment", "loc_external" = "Holdout location")[x]
# Replace ec selection
f_ec_selection_replace <- function(x)  c("stepwise_cv_adhoc" = "StepwiseCV", "stepAIC_adhoc" = "StepwiseAIC", 
                                         "lasso_cv_adhoc" = "LASSO",  "apriori" = "KnownStresses", 
                                         "all" = "All", "none" = "None")[x]
f_growth_stage_replace <- function(x) 
  c("early_vegetative" = "EV", "late_vegetative" = "LV", "heading" = "HD", "flowering" = "FL", "grain_fill" = "GF")[x]







## Remove
rm(s2_discrete_mat, s2_metadata, s2_imputed_mat, s2_tidy_BLUE)


# # Save relevant data in a single binary file
# phenotype_data <- S2_MET_BLUEs
# marker_data <- s2_imputed_mat_use
# marker_metadata <- snp_info
# phenotype_metadata <- trial_info
# training_environments <- train_test_env
# prediction_environments <- validation_env
# 
# save("phenotype_data", "phenotype_metadata", "marker_data", "marker_metadata",
#      "training_environments", "prediction_environments",
#      "trial_irrigation_data",
#      file = "../../GEPredictionPipeline/data/example_data.RData")
  
  
  