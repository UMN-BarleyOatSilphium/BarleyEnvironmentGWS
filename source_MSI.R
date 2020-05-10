## S2MET source for MSI
## 
## A script that automatically loads the data relevant for the S2MET project

# Load packages
packages <- c("sommer", "dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "rrBLUP", "purrr", "boot", "pbr", "lme4", "modelr", "neyhart")

invisible(lapply(packages, library, character.only = TRUE))

## Directories
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_TrainingDesign/"


geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
enviro_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/EnvironmentalData"
meta_dir <- pheno_dir

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(alt_proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


## Source the 'source.R' script from a define starting point
source_lines <- readLines(file.path(repo_dir, "source.R"))
source_lines_discard <- seq(which(grepl(pattern = "^# MSI", x = source_lines)))
source_lines_run <- source(textConnection(source_lines[-source_lines_discard]))

