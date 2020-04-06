## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-location-out prediction
## 
## Author: Jeff Neyhart
## Last modified: 27 Feb 2020
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

# Load the full model variance components
load(file.path(result_dir, "full_models.RData"))






### Models

# Cross-validation will use the following models:
# model1: y = G + r
# model2: y = G + E + r
# model3: y = G + E + GE + r
# model4: y = G + L + r
# model5: y = G + L + GL + r
# 
# Models with the suffix 'a' use the AMMI covariance matrix, instead of the one
# derived from ECs
# 
# Models with the suffix 'b' use the same covariance matrix for E and GE
# 



### CV0 and POV00

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq_len(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment, loc = location)


## 
## Leave-one-location-out
## 



# Generate skeleton train/test sets for LOLO
lolo_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo2(data = droplevels(group_by(., location)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(env, loc), fct_contr_sum) %>%
                       mutate(env1 = env, loc1 = loc)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(env, loc), fct_contr_sum) %>%
                      mutate(env1 = env, loc1 = loc)) )
  

## Assign cores and split
data_train_test1 <- lolo_train_test %>% 
  assign_cores(df = ., n_core = n_core, split = TRUE)


# Iterate over rows
# for (i in seq(nrow(loeo_predictions_out))) {
# for (i in seq(i, nrow(loeo_predictions_out))) {


## Parallelize
lolo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      # The genomewide prediction function is in the source_functions.R script
      prediction_out <- genomewide_prediction(x = row)
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(prediction_out$prediction_out, train_n = list(prediction_out$train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in location", as.character(row$location), "complete.")
      
    } # Close loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()


## Save the results
save("lolo_predictions_out", file = file.path(result_dir, "lolo_predictions.RData"))







