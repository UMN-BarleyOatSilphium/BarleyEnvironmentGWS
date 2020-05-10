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
  filter(population %in% c("s2tp", "s2c1", "s2met"), type == "spy") %>%
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


# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))


# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Calculate the K matrix
K <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)




## Rank the environments according to heritability
env_herit_rank <- s2_metadata %>% 
  select(trial, trait, heritability) %>% 
  left_join(., distinct(s2_tidy_BLUE, trial, environment)) %>%
  filter(!str_detect(trial, "S2C1"),
         trial %in% trials) %>%
  select(trait, environment, heritability) %>%
  # filter for relevant traits
  filter(trait %in% traits) %>%
  split(.$trait) %>%
  map(~arrange(., desc(heritability)) %>%
        mutate(environment = factor(environment, levels = .$environment)))

## Remove environments with low heritability
## This will be the df for correcting for heritability when calculating prediction accuracy
env_trait_herit <- env_herit_rank %>%
  map(mutate, environment = as.character(environment)) %>% 
  map_df(~filter(., heritability >= 0.10))


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
  # Remove irrigated trials
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "EON17" & trait == "HeadingDate"),
         !(environment == "KNY16" & trait == "TestWeight")) %>%
  # Rename and reorder
  select(trial, environment, location, year, trait, line_name, value, std_error = std.error)



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
S2_MET_BLUEs <- filter(S2_MET_BLUEs, environment %in% tp_vp_env) %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca")) %>%
  arrange(trait, environment)



# ######################
# ### Renaming and unit vectors
# ######################
# 
# ## Load the ec model file
# ec_model_filename <- file.path(result_dir, "ec_model_building.RData")
# 
# if (file.exists(ec_model_filename)) {
#   
#   load(ec_model_filename)
#   
#   ## Units for covariate variables
#   covariate_variable_rename <- c("mint" = "T[min]", "maxt" = "T[max]", "tmean" = "T[mean]",
#                                  "water_stress" = "drought", "radn" = "R")
#   
#   covariate_variable_unit <- covariate_variable_rename
#   names(covariate_variable_unit) <- c("degree*'C'", "degree*'C'", "degree*'C'", "mm", "MJ~m^-2")
#   
#   # Temporary renaming vector
#   covariate_variable_unit_temp <- str_replace_all(covariate_variable_unit, "\\[", "\\\\[") %>%
#     str_replace_all(., "\\]", "\\\\]") %>%
#     setNames(., names(covariate_variable_unit))
#   
#   ## Names for growth stages
#   growth_stage_rename <- c("flowering" = "FT", "vegetative" = "VG", "grain_fill" = "GF")
#   
#   ## Rename vector
#   covariate_rename <- str_replace_all(names(ec_tomodel)[-1], growth_stage_rename) %>% 
#     str_replace_all(., covariate_variable_rename) %>% 
#     str_replace_all(., "_", " ") %>% 
#     str_split(., " ") %>% 
#     map(rev) %>% 
#     map_chr(~paste0(.[1], "[(", .[2], ")]"))
#   
#   names(covariate_rename) <- names(ec_tomodel)[-1]
#   
#   ## Units
#   covariate_units <- names(ec_tomodel)[-1] %>%
#     setNames(nm = ., object = names(ec_tomodel)[-1] %>% 
#                str_replace_all(., covariate_rename) %>% 
#                map_chr(., ~names(covariate_variable_unit_temp)[str_which(string = .x, pattern = covariate_variable_unit_temp)]) )
#   
# 
# 
# }


## Trait units - renaming vector
trait_units <- setNames(object = c("kg~ha^-1", "days", "cm", "g~L^-1", "'%'"), nm = traits)
# The favorable sign for each trait
trait_sign <- tibble(trait = traits, sign = c(1, -1, -1, 1, -1))


## Remove
rm(s2_discrete_mat, s2_imputed_mat, env_herit_rank, s2_metadata, s2_tidy_BLUE)




