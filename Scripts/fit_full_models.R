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
  filter(model == "model3_fwd")


## Create a results df
full_model_df <- data_to_model %>% 
  group_by(trait) %>% 
  nest() %>%
  mutate(results = list(NULL))

## Loop over rows
for (i in seq(nrow(full_model_df))) {
  
  # Trait
  tr <- full_model_df$trait[i]
  data <- full_model_df$data[[i]]
  
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
  E_main <- Env_mat(x = ec_mat[, main_environment_covariates, drop = FALSE], method = "Rincent2019")
  E_int <- Env_mat(x = ec_mat[, interaction_environment_covariates, drop = FALSE], method = "Rincent2019")
  GE <- kronecker(X = K, Y = E_int, make.dimnames = TRUE)
    
  
  # Fit each of models 1, 2, 3
  full_model_out <- tibble(model = paste0("model", 1:3), object = list(NULL)) 
  
  for (m in seq(nrow(full_model_out))) {
    full_model_out$object[[m]] <- predict_gv(train = data, model = full_model_out$model[m],
                                             relMat.list = list(G = K, E = E_main, GE = GE), verbose = T)$object

  }
  
  ## Add the model output to the results df
  full_model_df$results[[i]] <- full_model_out
  
}


## Unnest, add logLik
full_model_df1 <- full_model_df %>%
  unnest(results) %>%
  mutate(logLik = map(object, "monitor") %>% map_dbl(~last(.[1,])),
         sigma = map(object, "sigma"))


## Select only variance components and loglik
full_models <- select(full_model_df1, -object)

# Save the df
save("full_models", file = file.path(result_dir, "full_models.RData"))



  
