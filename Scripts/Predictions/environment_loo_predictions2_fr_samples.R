## S2MET Prediction Models
## 
## Environment-specific or location-specific predictions
## These predictions use the covariates selected using the LOO
## procedure; so environments and locations are totally unobserved, even
## in variable selection.
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

# Load the covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

# Load the environmental covariates selected through the sampling procedure
load(file.path(result_dir, "factorial_regression_results_sample.RData"))


## For each environment or location, fit models 2 and 3 or 4 and 5 using
## the covariates. No other models are necessary



# Leave-one-environment-out -----------------------------------------------

## List of models
## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value,
  model2_cov = ~ 1,
  model3_cov = model2_cov,
)

## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model2_cov = ~ vs(line_name, Gu = K) + vs(env1, Gu = E),
  model3_cov = add_predictors(model2_cov, ~ vs(line_name:env1, Gu = GE)),
) %>% map(~ formula(delete.response(terms(.)))) # Remove response

# Combine into list
model.list <- list(fixed = model_fixed_forms, random = model_rand_forms)




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
covariates_tomodel <- concurrent_fact_reg_sample %>%
  # Only subset model 3, models 2 and 3 are run; only use ad hoc covariates
  filter(model == "model3", selection != "apriori") %>%
  select(trait, env = dropped_group, selection, term = covariates) %>%
  unnest(term) %>%
  group_by(trait, env, selection) %>%
  nest(.key = "covariates")

## DF of covariate values to use
ec_data_tomodel <- ec_tomodel_scaled


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
loeo_predictions_sample_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nEnv = n_distinct(environment), nObs = n())
      
      # List of covariates
      covariate_list <- row$covariates[[1]] %>% 
        mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
               covariate = str_remove(term, "line_name:")) %>%
        split(.$class) %>% 
        map("covariate")
      
      # Create a matrix of scaled and centered covariates
      covariate_mat <- ec_data_tomodel %>%
        filter(environment %in% levels(row$train[[1]]$env1)) %>%
        select(., environment, unique(unlist(covariate_list))) %>%
        as.data.frame() %>%
        column_to_rownames("environment") %>%
        as.matrix()
      
      ## Environmental relationship matrices
      Emain <- Env_mat(x = covariate_mat[,covariate_list$main, drop = FALSE], method = "Jarq")
      if (is.null(covariate_list$interaction)) {
        Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
      } else {
        Eint <- Env_mat(x = covariate_mat[,covariate_list$int, drop = FALSE], method = "Jarq")
      }
      
      
      # The genomewide prediction function is in the source_functions.R script
      prediction_out <- genomewide_prediction2(x = row, model.list = model.list,  K = K, E = Emain, KE = Eint)
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(prediction_out$prediction_out, train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in environment", as.character(row$env), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  })




# Leave-one-year-out -----------------------------------------------

## List of models
## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value,
  model4_cov = ~ 1,
  model5_cov = model4_cov,
)

## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model4_cov = ~ vs(line_name, Gu = K) + vs(loc1, Gu = E),
  model5_cov = add_predictors(model4_cov, ~ vs(line_name:loc1, Gu = GE)),
) %>% map(~ formula(delete.response(terms(.)))) # Remove response

# Combine into list
model.list <- list(fixed = model_fixed_forms, random = model_rand_forms)




# Data to use - location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
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
covariates_tomodel <- historical_fact_reg_sample %>%
  mutate(model = ifelse(model == "model3", "model5", model)) %>%
  # Only subset model 5, models 4 and 5 are run; only use ad hoc covariates
  filter(model == "model5", selection != "apriori") %>%
  select(trait, loc = dropped_group, selection, term = covariates) %>%
  unnest(term) %>%
  group_by(trait, loc, selection) %>%
  nest(.key = "covariates")

## DF of covariate values to use
ec_data_tomodel <- historical_ec_tomodel_timeframe_scaled[[time_frame_use]]


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


## Parallelize
lolo_predictions_sample_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      # Record the number of location and observations used for training
      train_n <- summarize(row$train[[1]], nLoc = n_distinct(location), nObs = n())
      
      # List of covariates
      covariate_list <- row$covariates[[1]] %>% 
        mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
               covariate = str_remove(term, "line_name:")) %>%
        split(.$class) %>% 
        map("covariate")
      
      # Create a matrix of scaled and centered covariates
      covariate_mat <- ec_data_tomodel %>%
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
      prediction_out <- genomewide_prediction2(x = row, model.list = model.list,  K = K, E = Emain, KE = Eint)
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(prediction_out$prediction_out, train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in location", as.character(row$loc), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  })




## Save the results
save("loeo_predictions_sample_out", "lolo_predictions_sample_out", 
     file = file.path(result_dir, "loo_predictions_fact_reg_samples.RData"))






