## S2MET Prediction Models
## 
## Fit full models to estimate variance components
## 
## Author: Jeff Neyhart
## Last modified: 8 January 2020
## 

# This script includes code to estimate variance components for each validation
# model using only the TP


# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(broom)


## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))
# Rename
ec_model_building <- unified_ec_models

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))




## Data.frame to use to fit models
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno)) %>%
  # add covariates
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno)),
         environment = as.factor(environment))

## Unnest the ec model
ec_model_touse <- ec_model_building %>% 
  unnest(final_model) %>%
  filter(model == "model3_ammi")
  


## List of models
## Lowercase = fixed, uppercase = random, eB = fixed environment using covariates
# model1: y = G + r
# model2: y = G + E + r
# model2a: y = G + e + r
# model2b: y = G + eB + r
# model3: y = G + E + GE + r
# model3a: y = G + e + GE + r
# model3b: y = G + eB + GE + r


## Character vector of models
model_list <- c("model1", "model2", "model2a", "model2b", "model3", "model3a", "model3b")




## Create a results df
full_model_df <- data_to_model %>% 
  group_by(trait) %>% 
  nest() %>%
  mutate(results = list(NULL))

## Loop over rows
for (i in seq(nrow(full_model_df))) {
# for (i in seq(i, nrow(full_model_df))) {
    
  
  # Trait
  tr <- full_model_df$trait[i]
  data <- full_model_df$data[[i]] %>%
    # Add covariates
    left_join(., ec_tomodel_centered) %>%
    # Need to convert environment back to factor
    mutate_at(vars(line_name, environment), as.factor) %>%
    rename(env = environment)
  
  # Grab the model
  ec_model_i <- subset(ec_model_touse, trait == tr, object, drop = TRUE)[[1]]
  # Extract the fitted model object formula
  ec_model_form <- formula(ec_model_i)
  
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
  
  
  
  ## Assign relationship matrices
  K <- K # Genomic
  E <- Env_mat(x = ec_mat[, main_environment_covariates, drop = FALSE], method = "Rincent2019")
  GE <- Env_mat(x = ec_mat[, interaction_environment_covariates, drop = FALSE], method = "Rincent2019") %>% 
    kronecker(X = K, Y = ., make.dimnames = TRUE)
  
  # # List of covariance matrices
  # relMat.list <- list(G = K, E = E_main, GE = GE)
  # 
  # # Covariate list
  # covariate.list <- list(E = main_environment_covariates)
  
  
  ###############
  ## Model fitting
  ###############
  
  ## Create a list of formulas
  model_fixed_forms <- formulas(
    .response = ~ value, 
    model1 = ~ 1,
    model2 = model1,
    model2a = add_predictors(model2, ~ env),
    model2b = reformulate(c("1", main_environment_covariates)),
    model3 = model2,
    model3a = model2a,
    model3b = model2b
  )
  
  model_rand_forms <- formulas(
    .response = ~ value,
    model1 = ~ vs(line_name, Gu = K),
    model2 = add_predictors(model1, ~ vs(env, Gu = E)),
    model2a = model1,
    model2b = model1,
    model3 = add_predictors(model2, ~ vs(line_name:env, Gu = GE)),
    model3a = add_predictors(model1, ~ vs(line_name:env, Gu = GE)),
    model3b = model3a
  ) %>% map(~ formula(delete.response(terms(.)))) # Remove response
    
  
  ## Map over pairs of formula and fit the models
  model_fits <- map2(.x = model_fixed_forms, .y = model_rand_forms, 
                     ~mmer(fixed = .x, random = .y, rcov = ~ units, data = data, date.warning = FALSE,
                           verbose = TRUE) )
  
  ## Convert to tibble, extract diagonistics and varcomp
  model_fits_diag <- model_fits %>%
    tibble(model = names(.), object = .) %>%
    mutate(logLik = map(object, "monitor") %>% map_dbl(~last(.[1,])), # LogLik
           predictions = map(object, predict),
           R2 = map_dbl(predictions, ~cor(.x$predictions$value, .x$predictions$predicted.value.value)^2),
           sigma = map(object, "sigma")) %>%
    select(-object, -predictions)
  
  # Add to list
  full_model_df$results[[i]] <- model_fits_diag
  
}


## Unnest, add logLik
full_models <- full_model_df %>%
  unnest(results)


## Plot R2 and logLik
full_models %>%
  ggplot(aes(x = model, y = logLik)) +
  geom_col() +
  facet_wrap(~ trait, scales = "free_y")

full_models %>%
  ggplot(aes(x = model, y = R2)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 1)) + 
  facet_wrap(~ trait, scales = "free_y")




# Save the df
save("full_models", file = file.path(result_dir, "full_models.RData"))



  
