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


# Other packages
library(modelr)
library(broom)
library(parallel)


## Number of cores
# n_core <- detectCores()
n_core <- 8

# time_frame to use for the location relationship matrix
time_frame_use <- "time_frame5"


## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))
load(file = file.path(result_dir, "historical_ec_model_building.RData"))

## Load the factorial regression results
load(file.path(result_dir, "factorial_regression_results.RData"))


# For each environment, fit 3 models:
# model1: y = G + r
# model2: y = G + E + r
# model3: y = G + E + GE + r

# E will be represented by covariates

# For each model, use three sets of covariates:
# adhoc
# adhoc no soil
# apriori

# For each set of covariates:
# 1. Create an environmental relationship matrix
# 2. Use regular factorial regression


# Total models per environment:
# 3 models x 3 covariate sets x 2 ways of fitting = 18 models


### CV0 and POV0

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment, loc = location)


# Reorganize covariate df
covariates_tomodel <- fr_var_summary %>%
  unnest(var_prop_summary) %>%
  filter(timeframe == "concurrent",
         group %in% c("main", "interaction")) %>%
  select(trait, selection, term) %>%
  group_by(trait, selection) %>%
  nest(.key = "covariates")


## 
## Leave-one-environment-out
## 



# Generate skeleton train/test sets for LOEO
loeo_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo_grouped(data = droplevels(group_by(., env)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(env, loc), fct_contr_sum) %>%
                       mutate(env1 = env, loc1 = loc)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(env, loc), fct_contr_sum) %>%
                      mutate(env1 = env, loc1 = loc)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)
  

## Assign cores and split
data_train_test1 <- loeo_train_test %>% 
  assign_cores(df = ., n_core = n_core, split = TRUE)


## Parallelize
loeo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      # The genomewide prediction function is in the source_functions.R script
      prediction_out <- genomewide_prediction2(x = row)
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(prediction_out$prediction_out, train_n = list(prediction_out$train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in environment", as.character(row$env), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()


## Save the results
save("loeo_predictions_out", file = file.path(result_dir, "loeo_predictions_fact_reg.RData"))






