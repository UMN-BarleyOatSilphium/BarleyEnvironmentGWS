## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-environment-out prediction
## 
## Author: Jeff Neyhart
## Last modified: 8 January 2020
## 


# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))

## Number of cores
n_core <- detectCores()
n_core <- 16


# Other packages
library(modelr)
library(broom)
library(parallel)



## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))
# Rename
ec_model_building <- unified_ec_models

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))
# Load the full model variance components
load(file.path(result_dir, "full_models.RData"))


### Models

# Cross-validation will use the following models:
# M1 - random genotype (K)
# M2 - M1 + random environment (E)
# M3 - M2 + random genotype:environment (KE)

models <- paste0("model", 1:3)


### CV0 and POV0

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)),
         environment = as.factor(environment))


## 
## Leave-one-environment-out
## 


# Generate skeleton train/test sets for LOEO
loeo_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo2(data = group_by(., environment))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno))))
  

## Assign cores and split
data_train_test1 <- loeo_train_test %>% 
  assign_cores(df = ., n_core = n_core, split = TRUE)


# Iterate over rows
# for (i in seq(nrow(loeo_predictions_out))) {
# for (i in seq(i, nrow(loeo_predictions_out))) {


## Parallelize
loeo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      # Get training and test data
      row <- core_df[i,]
      tr <- unique(row$trait)
      
      # train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
      #   mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
      # test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
      #   select(environment, line_name, value)
      
      ## Different subsetting method
      train <- row$train[[1]]
      test <- row$test[[1]]
      
      ## Get the covariate model
      ec_model_i <- ec_model_building %>% 
        unnest(final_model) %>%
        subset(trait == row$trait & model == "model3_fwd", object, drop = TRUE)
      
      # Extract the fitted model object formula
      ec_model_form <- formula(ec_model_i[[1]])
      
      # main_environment_covariates
      main_environment_covariates <- all.vars(ec_model_form) %>% 
        subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
                                             pattern = .))))
      
      # interaction_environment_covariates
      interaction_environment_covariates <- all.vars(ec_model_form) %>% 
        subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
                                             pattern = .)))) %>% setdiff(., "line_name")
      
      # All covariates
      all_covariates <- union(main_environment_covariates, interaction_environment_covariates)
      
      ## Subset covariates and convert to a matrix
      ec_mat <- ec_tomodel_scaled %>% 
        select(environment, all_covariates) %>% 
        as.data.frame() %>%
        column_to_rownames("environment") %>%
        as.matrix()
      
      
      ## Create relationship matrices
      E <- Env_mat(x = ec_mat[,main_environment_covariates, drop = FALSE], method = "Rincent2019")
      KE <- Env_mat(x = ec_mat[,interaction_environment_covariates, drop = FALSE], method = "Rincent2019") %>%
        kronecker(X = K, Y = ., make.dimnames = TRUE)
      
      
      
      #################
      ## Fit models and extract predictions
      #################
      
      prediction_out <- setNames(object = vector("list", length(models)), models)
      
      # Iterate over models
      for (mod in models) {
        
        ## Get the variance components from the full model
        full_varcomp <- subset(full_models, trait == tr & model == mod)["sigma"] %>%
          map(1)
        
        # Fit the model, extract predictions
        model_out <- predict_gv(train = train, test = test, model = mod, relMat.list = list(G = K, E = E, GE = KE), 
                                object = full_varcomp, add.intercept = TRUE, verbose = TRUE)
        
        # Add predictions to the list
        prediction_out[[mod]] <- model_out
        
      }
      

      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- imap_dfr(prediction_out, ~tibble(model = .y, predictions = .x))
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in environment", as.character(row$environment), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()
