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

# Load the full model variance components
load(file.path(result_dir, "full_models.RData"))



## In the location relationship matrices, convert Ithaca1 and Ithaca2 to Ithaca
location_relmat_df <- location_relmat_df %>%
  mutate_at(vars(contains("E_mat")), 
            # First convert all Ithaca[0-9] to Ithaca; then select the first one
            ~map(., ~`dimnames<-`(.x, lapply(dimnames(.x), str_replace, pattern = "Ithaca1|Ithaca2", replacement = "Ithaca")) %>%
                   subset(., !duplicated(row.names(.)), !duplicated(colnames(.))) ))
                   

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


## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value, 
  model1 = ~ 1,
  model2 = model1,
  # model2a = add_predictors(model2, ~ env),
  # model2b = reformulate(c("1", main_environment_covariates)),
  model3 = model2,
  # model3a = model2a,
  # model3b = model2b
  model4 = model1,
  model5 = model4,
  model2a = model2,
  model3a = model3,
  model4a = model4,
  model5a = model5
)


## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model1 = ~ vs(line_name, Gu = K),
  model2 = add_predictors(model1, ~ vs(env, Gu = E)),
  # model2a = model1,
  # model2b = model1,
  model3 = add_predictors(model2, ~ vs(line_name:env1, Gu = GE)),
  # model3a = add_predictors(model1, ~ vs(line_name:env1, Gu = GE)),
  # model3b = model3a
  model4 = add_predictors(model1, ~vs(loc, Gu = L)),
  model5 = add_predictors(model4, ~vs(line_name:loc, Gu = GL)),
  # Modified models
  model2a = add_predictors(model1, ~ vs(env, Gu = E_IPC)),
  model3a = add_predictors(model2, ~ vs(line_name:env1, Gu = GE_IPC)),
  model4a = add_predictors(model1, ~vs(loc, Gu = L_IPC)),
  model5a = add_predictors(model4, ~vs(line_name:loc, Gu = GL_IPC)),
) %>% map(~ formula(delete.response(terms(.)))) # Remove response

# Residual formula
resid_form <- ~ vs(units)


### CV0 and POV0

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca")) %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment, loc = location)


## 
## Leave-one-environment-out
## 



# Generate skeleton train/test sets for LOEO
loeo_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo2(data = droplevels(group_by(., env)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(env, loc), fct_contr_sum) %>%
                       mutate(env1 = env, loc1 = loc)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(env, loc), fct_contr_sum) %>%
                      mutate(env1 = env, loc1 = loc)) )
  

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
      prediction_out <- genomewide_prediction(x = row)
      
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
save("loeo_predictions_out", file = file.path(result_dir, "loeo_predictions.RData"))






