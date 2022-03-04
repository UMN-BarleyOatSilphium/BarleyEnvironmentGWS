## S2MET source for SCINet
## 
## A script that automatically loads the data relevant for the S2MET project

# Load packages
packages <- c("sommer", "dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "rrBLUP", "purrr", "boot", "lme4", "modelr", "neyhart")

invisible(lapply(packages, library, character.only = TRUE))

## Directories
proj_dir <- "/project/gifvl_vaccinium/barley_work/BarleyEnvironmentGWS"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")

geno_dir <-  "/project/gifvl_vaccinium/barley_work/Data/GenotypicData"
pheno_dir <- "/project/gifvl_vaccinium/barley_work/Data/PhenotypicData"
meta_dir <- pheno_dir


## Source the 'source.R' script from a define starting point
source_lines <- readLines(file.path(repo_dir, "source.R"))
source_lines_discard <- seq(which(grepl(pattern = "^# MSI", x = source_lines))) # Keep the MSI line here, even for SCINet
source_lines_run <- source(textConnection(source_lines[-source_lines_discard]))

