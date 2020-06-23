## S2MET Prediction Models source script
## 
## A script that automatically loads the data relevant for the S2MET project


# List of packages
pkgs <- c("sommer", "tidyverse", "readxl", "rrBLUP", "neyhart", "pbr")
# Load these packages
invisible(lapply(X = pkgs, library, character.only = TRUE))


## Directories
proj_dir <- repo_dir
alt_proj_dir <- "C:/GoogleDrive/BarleyLab/Projects/S2MET_Predictions//"

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



######
# MSI Source starts here
######


# Source the project functions
source(file.path(proj_dir, "source_functions.R"))


# Load the phenotypic data
load(file.path(pheno_dir, "S2_tidy_BLUE.RData"))

# Load the trial metadata
trial_info <- read_csv(file = file.path(meta_dir, "trial_metadata.csv")) %>%
  filter(population %in% c("s2tp", "s2c1", "s2c1r", "s2met"), type == "spy") %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca"))


# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))


# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Grab the entry names that are not checks
tp <- subset(entry_list, Class == "S2TP", Line, drop = T)
tp <- setdiff(tp, "07MT-10") # Remove hullness line
vp <- subset(entry_list, Class == "S2C1R", Line, drop = T)


# Vector of relevant traits
# traits <- c("GrainYield", "HeadingDate", "PlantHeight", "MaturityDate")
traits <- c("GrainYield", "HeadingDate", "PlantHeight", "TestWeight", "GrainProtein")
trials <- subset(trial_info, project2 == "S2MET", trial, drop = TRUE)


# Find the tp and vp that are genotyped with markers
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))


# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Calculate the K matrix
K <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)




## Rank the environments according to heritability
env_trait_herit <- s2_metadata %>% 
  select(trial, trait, heritability, varR) %>% 
  left_join(., distinct(s2_tidy_BLUE, trial, environment), by = "trial") %>%
  filter(!str_detect(trial, "S2C1"),
         trial %in% trials) %>%
  select(trait, environment, varR, heritability) %>%
  # filter for relevant traits
  filter(trait %in% traits) %>%
  filter(., heritability >= 0.10)


## Filter the S2 tidy blues for S2MET
S2_MET_BLUEs <- s2_tidy_BLUE %>%
  filter(trial %in% subset(trial_info, project2 == "S2MET", trial, drop = T),
         line_name %in% c(tp, vp),
         trait %in% traits) %>%
  # Filter out environments with low heritability
  left_join(select(env_trait_herit, -heritability), .) %>%
  # Add full location name
  select(-location) %>%
  left_join(., distinct(trial_info, environment, location)) %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM")) %>%
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "EON17" & trait == "HeadingDate"),
         !(environment == "KNY16" & trait == "TestWeight"),
         !(location == "Charlottetown" & trait == "HeadingDate")) %>%
  # Rename and reorder
  select(trial, environment, location, year, trait, line_name, value, std_error = std.error)


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
model_replace <- c("model1" = "g", "model2_id" = "g + E", "model2_cov" = "g + w", 
                   "model3_id" = "g + E + gE", "model3_cov" = "g + w + gu", "model3_cov1" = "g + w + gw")

## Models to present
# model_present <- model_replace[str_detect(names(model_replace), "_", negate = TRUE)]
model_present <- model_replace


f_model_replace <- function(x) model_replace[x]
f_model_replace2 <- function(x) model_present[x]
# f_model_replace <- function(x) paste0("M", toupper(str_extract(x, "[0-9]{1}[a-z]{0,1}")))
# Vector to rename validation schemes
f_pop_replace <- function(x) str_replace_all(x, c("tp" = "Cross-validation", "vp" = "Offspring validation"))
# Replace type
f_type_replace <- function(x) c("loeo" = "New environment", "lolo" = "New location", "loyo" = "New year")[x]
# Replace ec selection
f_ec_selection_replace <- function(x)  c("rfa_cv_adhoc" = "stepCV", "stepAIC_adhoc" = "stepAIC", "apriori" = "literature", 
                                         "all" = "All", "none" = "None")[x]
f_growth_stage_replace <- function(x) 
  c("early_vegetative" = "EV", "late_vegetative" = "LV", "heading" = "HD", "flowering" = "FL", "grain_fill" = "GF")[x]




## Remove
rm(s2_discrete_mat, s2_metadata, s2_imputed_mat, s2_tidy_BLUE)




