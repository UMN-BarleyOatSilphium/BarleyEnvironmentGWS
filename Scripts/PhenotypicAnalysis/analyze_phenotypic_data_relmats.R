## S2MET Prediction Models
## 
## Fit full models to estimate variance components
## 
## Author: Jeff Neyhart
## Last modified: 8 January 2020
## 

# This script will analyze the phenotypic data using covariates and markers
# 

# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))

# Other packages
library(modelr)
library(broom)
# library(lme4qtl)

load(file = file.path(result_dir, "feature_selection_results.RData"))
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

# Set the number of cores
n_cores <- 8




## Data.frame to use to fit models
f_return <- function(x) x
data_to_model <- S2_MET_BLUEs %>% 
  filter(environment %in% train_test_env, line_name %in% c(tp_geno, vp_geno)) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter(!(population == "tp" & line_name %in% vp_geno), !(population == "vp" & line_name %in% tp_geno)) %>%
  # Copy columns for modeling covariates
  droplevels() %>%
  mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.))) %>%
  mutate(gxe = interaction(line_name, environment, drop = TRUE, sep = ":")) %>%
  mutate_at(vars(line_name, environment, gxe), list(cov = f_return))

## Assemble covariates to model
covariates_to_model <- concurrent_fact_reg_feature_selection %>% 
  select(trait, model, source, apriori) %>% 
  gather(feat_sel_type, features, apriori) %>%
  bind_rows(., rename(concurrent_feature_selection, features = adhoc)) %>%
  mutate(features = map(features, "optVariables") %>% map(~setdiff(x = ., y = "line_name"))) %>%
  filter(source == "daymet", model == "model3") %>%
  select(trait, feat_sel_type, features) %>%
  bind_rows(., tibble(trait = traits, feat_sel_type = "all", features = list(names(ec_tomodel_scaled$daymet)[-1:-2])))


# Matrix of covariates
covariate_mat <- ec_tomodel_scaled$daymet %>%
  as.data.frame() %>%
  select(-source) %>%
  column_to_rownames("environment") %>%
  as.matrix()
  

## Create a results df
data_to_model_split <- data_to_model %>% 
  group_by(trait, population) %>% 
  nest() %>%
  ungroup() %>%
  full_join(., covariates_to_model, by = "trait") %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE)

## Iterate over cores
pheno_variance_analysis_out <- coreApply(X = data_to_model_split, FUN = function(core_df) {
  
  # Output list
  out <- vector("list", length = nrow(core_df))
  
  # Iterate over rows
  for (i in seq_along(out)) {
  
    row <- core_df[i,]
    
    # Trait
    tr <- row$trait
    data <- droplevels(row$data[[1]])
    genotypes <- levels(data$line_name)
    environments <- levels(data$environment)
    covariates <- row$features[[1]]
    
    # Separate covariates into main effect or interaction
    covariate_list <- tibble(term = unique(covariates)) %>% 
      mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
             covariate = str_remove(term, "line_name:")) %>%
      split(.$class) %>% 
      map("covariate")
    
    
    # Create or subset relationship matrices
    K_use <- K[genotypes, genotypes]
    Zg <- model.matrix(~ -1 + line_name, data)
    Emain <- Eint2 <- Env_mat(x = covariate_mat[environments, covariate_list$main, drop = FALSE], method = "Jarq")
    Eint <- Env_mat(x = covariate_mat[environments, covariate_list$int, drop = FALSE], method = "Jarq")
    # Replace with diagonal matrix if NA
    if (all(is.na(Eint))) { 
      Eint[] <- 0
      diag(Eint) <- 1 
      no_gxe_cov <- TRUE
      
    } else {
      no_gxe_cov <- FALSE
      
    }
    
    Ze <- model.matrix(~ -1 + environment, data)
    KE <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% Eint, Ze)
    dimnames(KE) <- replicate(2, data$gxe, simplify = FALSE)
    
    # # Fit the model
    # fit <- relmatLmer(value ~ (1|line_name_cov) + (1|line_name) + (1|environment_cov) + (1|environment) + 
    #                     (1|gxe), 
    #                   data = data, relmat = list(line_name_cov = K_use, environment_cov = Emain))
    
    ## Sommer
    model_stdout <- capture.output({ 
      model_try <- try( {
        fit_mmer <- mmer(fixed = value ~ 1, 
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + environment + vs(environment_cov, Gu = Emain) +
                           gxe + vs(gxe_cov, Gu = KE),
                         data = data, verbose = TRUE) 
      }, silent = TRUE) })
    
    # Retry with fewer iterations if singular
    if (class(model_try) == "try-error") {
      
      # Parse the monitor for iterations
      model_monitor <- read_table(head(model_stdout, -1))
      
      # Select 1 iteration fewer
      iter_max <- max(model_monitor$iteration) - 1
      
      # refit the model
      fit_mmer <- mmer(fixed = value ~ 1, 
                       random = ~ line_name + vs(line_name_cov, Gu = K_use) + environment + vs(environment_cov, Gu = Emain) +
                         gxe + vs(gxe_cov, Gu = KE),
                       data = data, verbose = TRUE, iter = iter_max)
      
    }
    
    
    # Get the variance components; asemble into a tibble
    varcomp_df <- fit_mmer$sigma %>% 
      map_dbl(~.) %>% 
      tibble(term = names(.), variance = .) %>% 
      mutate(term = str_remove_all(term, "u:"), 
             source = str_remove_all(term, "_cov")) %>%
      group_by(source) %>%
      mutate(total_source_variance = sum(variance)) %>%
      ungroup() %>%
      mutate(variance_prop = variance / total_source_variance,
             note = list(NULL))
    
    # Return
    out[[i]] <- varcomp_df
    
  } # Close loop
  
  # Add out to core_df and return
  core_df %>%
    mutate(results = out) %>%
    select(-data, -core)
  
})


# Rename
pheno_variance_analysis <- pheno_variance_analysis_out
  

# Save the df
save("pheno_variance_analysis", file = file.path(result_dir, "full_model_variance_analysis.RData"))
