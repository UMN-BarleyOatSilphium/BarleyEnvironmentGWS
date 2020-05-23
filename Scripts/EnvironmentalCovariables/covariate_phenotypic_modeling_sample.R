## S2MET phenotypic modeling with covariates
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## Model the impact of covariates on GxE and LxE
## 


# ## Local machine
# # Repository directory
# repo_dir <- getwd()
# # Source the main project script
# source(file.path(repo_dir, "source.R"))

## MSI
# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))


# Other packages
library(modelr)
library(broom)
library(parallel)
library(car)


## significance level
alpha <- 0.05

# Number of cores
n_cores <- 8



# Load data ---------------------------------------------------------------

# Load covariates for environments and historical covariates
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


## Filter BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  mutate_at(vars(line_name, environment), as.factor)


# Modeling ----------------------------------------------------------------


# Concurrent environmental covariates =====================================


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = names(ec_tomodel_centered)[-1]) %>%
  filter(
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 



# Use the following covariates for all traits (except where timing would make
# it inappropriate):
# 
# radiation during vegetative
# tmin during flowering
# water stress during flowering
# grain fill water stress
# grain fill tmax


## A priori covariates - quite minimal; each should have a citation 
apriori_covariate_df <- trait_covariate_df %>%
  filter(covariate %in% c("early_vegetative.gdd_sum", "late_vegetative.gdd_sum", "flowering.mint_mean",
                          "grain_fill.maxt_mean", "grain_fill.water_balance_sum"))



## Group by trait and nest
concurrent_fact_reg_sample_data <- S2_MET_BLUEs_tomodel %>%
  # Split by trait
  split(.$trait) %>%
  # LOO based on environment grouping
  imap_dfr(~group_by(.x, environment) %>% crossv_loo_grouped(.) %>% mutate(trait = .y)) %>%
  rename(dropped_group = environment)

concurrent_fact_reg_sample <- concurrent_fact_reg_sample_data %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
    
      row <- core_df[i,]
      df <- as_tibble(row$train[[1]])
      
      # Factorize
      df1 <- df %>%
        droplevels() %>%
        left_join(., ec_tomodel_centered, by = "environment") %>%
        mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
      
      test <- as_tibble(row$test[[1]]) %>%
        left_join(., ec_tomodel_centered, by = "environment")
      
      # Apriori
      ## Add covariates - filter
      covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
      apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "apriori")
      
      # Ad hoc
      ## Add covariates - filter
      covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
      adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step")
      
      # Ad hoc - without soil
      ## Add covariates - filter
      covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T), 
                               covariate, drop = TRUE)
      adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step")
      
      ## df of output
      results <- tibble(model = c("base", "base_alt", "model2", "model3"),
                        apriori = apriori_out,
                        adhoc = adhoc_out,
                        adhoc_nosoil = adhoc_nosoil_out)
      
      
      ## Predict the test set using each model
      ## return these results
      out[[i]] <- results %>%
        gather(selection, fit, -model) %>%
        filter(model != "base_alt") %>%
        mutate(test_predictions = map(fit, ~add_predictions(data = test, model = .) %>% select(line_name, value, pred)),
               covariates = map(fit, ~str_subset(string = attr(terms(formula(.)), "term.labels"), pattern = "line_name$", negate = TRUE)),
               rmse = map_dbl(fit, ~rmse(model = ., data = test)),
               acc = map_dbl(test_predictions, ~cor(.$value, .$pred))) %>%
        select(-fit)
      
    } # Close the loop
    
    # Add results to the core_df and return
    core_df %>% 
      mutate(results = out) %>% 
      select(dropped_group, trait, results) %>%
      unnest(results)
    
  }) %>% bind_rows()

  









# Historical covariates =====================================================


## Ideas
## 1. Calculate locations means and use that as the input (and validation)
## 2. Use all data and fit a model with random year, location:year, and g:year
## 3. Use all data and fit the g + location model


## Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp)


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = names(historical_ec_tomodel_centered$time_frame5)[-1:-2]) %>%
  filter(
    covariate != "awc_range",
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 

## Select the historical covariate data
historical_ec_tomodel_centered_use <- historical_ec_tomodel_centered$time_frame5
  


## Drop one location at a time and see if the same covariates
## are identified
## Group by trait and model
historical_fact_reg_sample_data <- S2_MET_loc_BLUEs_tomodel %>%
  # Split by trait
  split(.$trait) %>%
  # LOO based on environment grouping
  imap_dfr(~group_by(.x, location) %>% crossv_loo_grouped(.) %>% mutate(trait = .y)) %>%
  rename(dropped_group = location)
  
historical_fact_reg_sample <- historical_fact_reg_sample_data %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
  
      row <- core_df[i,]
      df <- as_tibble(row$train[[1]])
      
      # Factorize
      df1 <- df %>%
        droplevels() %>%
        left_join(., historical_ec_tomodel_centered_use, by = "location") %>%
        mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
      
      test <- as_tibble(row$test[[1]]) %>%
        left_join(., historical_ec_tomodel_centered_use, by = "location")
      
      # Apriori
      ## Add covariates - filter
      covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
      apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "apriori")
      
      # Ad hoc
      ## Add covariates - filter
      covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
      adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")
      
      # Ad hoc - without soil
      ## Add covariates - filter
      covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T), 
                               covariate, drop = TRUE)
      adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")
      
      ## df of output
      results <- tibble(model = c("base", "base_alt", "model2", "model3"),
                        apriori = apriori_out,
                        adhoc = adhoc_out,
                        adhoc_nosoil = adhoc_nosoil_out)
      
      
      ## Predict the test set using each model
      ## return these results
      out[[i]] <- results %>%
        gather(selection, fit, -model) %>%
        filter(model != "base_alt") %>%
        mutate(test_predictions = map(fit, ~add_predictions(data = test, model = .) %>% select(line_name, value, pred)),
               covariates = map(fit, ~str_subset(string = attr(terms(formula(.)), "term.labels"), pattern = "line_name$", negate = TRUE)),
               rmse = map_dbl(fit, ~rmse(model = ., data = test)),
               acc = map_dbl(test_predictions, ~cor(.$value, .$pred))) %>%
        select(-fit)
      
    } # Close the loop
    
    # Add results to the core_df and return
    core_df %>% 
      mutate(results = out) %>% 
      select(dropped_group, trait, results) %>%
      unnest(results)
    
  }) %>% bind_rows()



## Save results
save("historical_fact_reg_sample", "concurrent_fact_reg_sample",
     file = file.path(result_dir, "factorial_regression_results_sample.RData"))
