## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-location-out prediction
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
n_core <- 12

# time_frame to use for the location relationship matrix
time_frame_use <- "time_frame5_2010_2014"

## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables_nasapower.RData"))

## Load the factorial regression results
# load(file.path(result_dir, "factorial_regression_results.RData"))
load(file.path(result_dir, "feature_selection_results_nasapower.RData"))


# For each environment, fit 3 models:
# model1: y = G + r
# model2: y = G + L + r
# model3: y = G + L + GL + r

# L will be represented by covariates or by identity covariance matrices


# For each model, use three sets of covariates:
# adhoc
# adhoc no soil
# apriori

# For each set of covariates:
# 1. Create an environmental relationship matrix
# 2. Use regular factorial regression


# Total models per environment:
# 3 models x 3 covariate sets x 2 ways of fitting = 18 models


## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value,
  model1 = ~ 1,
  # model2_fr = reformulate(covariate_list$main),
  model4_cov = model1,
  model4_id = model4_cov,
  # model3_fr = model2_fr,
  model5_cov = model1,
  model5_id = model5_cov
)

## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model1 = ~ vs(line_name, Gu = K),
  # model2_fr = model1,
  model4_cov = add_predictors(model1, ~ vs(loc1, Gu = E)),
  model4_id = add_predictors(model1, ~ vs(loc1, Gu = I)),
  # model3_fr = model3_fr_rand,
  model5_cov = add_predictors(model4_cov, ~ vs(line_name:loc1, Gu = GE)),
  model5_id = add_predictors(model4_cov, ~ vs(line_name:loc1, Gu = GI))
) %>% map(~ formula(delete.response(terms(.)))) # Remove response

# Combine into list
model.list <- list(fixed = model_fixed_forms, random = model_rand_forms)

### CV0 and POV0

# Data to use - location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

data_to_model <- S2_MET_loc_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(location), as.factor) %>%
  mutate(loc = location)


# Reorganize covariate df
covariates_tomodel <- fr_var_summary %>%
  unnest(var_prop_summary) %>%
  filter(timeframe == "historical", # Historical for locations
         group %in% c("main", "interaction")) %>%
  select(trait, selection, term) %>%
  group_by(trait, selection) %>%
  nest(.key = "covariates")


## 
## Leave-one-location-out
## 



# Generate skeleton train/test sets for LOEO
lolo_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo_grouped(data = droplevels(group_by(., loc)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(loc), fct_contr_sum) %>%
                       mutate(loc1 = loc)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(loc), fct_contr_sum) %>%
                      mutate(loc1 = loc)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)
  

## Assign cores and split
data_train_test1 <- lolo_train_test %>% 
  assign_cores(df = ., n_core = n_core, split = TRUE)

## Only use models 2 and 3 with covariates; save the others for a separate run, 
## since they will be the same
model.list1 <- lapply(X = model.list, "[", c("model4_cov", "model5_cov"))


## Parallelize
lolo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nLoc = n_distinct(location), nObs = n())
      
      # List of covariates
      covariate_list <- row$covariates[[1]] %>% 
        mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
               covariate = str_remove(term, "line_name:")) %>%
        split(.$class) %>% 
        map("covariate")
      
      # Create a matrix of scaled and centered covariates
      covariate_mat <- historical_ec_tomodel_timeframe_scaled[[time_frame_use]] %>%
        filter(location %in% levels(row$train[[1]]$loc1)) %>%
        select(., location, unique(unlist(covariate_list))) %>%
        as.data.frame() %>%
        column_to_rownames("location") %>%
        as.matrix()
      
      ## Environmental relationship matrices
      Emain <- Env_mat(x = covariate_mat[,covariate_list$main, drop = FALSE], method = "Jarq")
      if (is.null(covariate_list$interaction)) {
        Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
      } else {
        Eint <- Env_mat(x = covariate_mat[,covariate_list$int, drop = FALSE], method = "Jarq")
      }
      
      # The genomewide prediction function is in the source_functions.R script
      prediction_out <- genomewide_prediction2(x = row, model.list = model.list1, K = K, E = Emain, KE = Eint)
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(prediction_out$prediction_out, train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in location", as.character(row$loc), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()


## Save the results
save("lolo_predictions_out", file = file.path(result_dir, "lolo_predictions_fact_reg.RData"))






