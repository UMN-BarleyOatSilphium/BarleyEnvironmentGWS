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
load(file = file.path(result_dir, "historical_ec_model_building.RData"))
# Rename
ec_model_building <- unified_ec_models



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
  


### Models

# Cross-validation will use the following models:
# model1: y = G + r
# model2: y = G + E + r
# model3: y = G + E + GE + r
# model4: y = G + L + r
# model5: y = G + L + GL + r


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
  model5 = model4
)


## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model1 = ~ vs(line_name, Gu = K),
  model2 = add_predictors(model1, ~ vs(env, Gu = E)),
  # model2a = model1,
  # model2b = model1,
  model3 = add_predictors(model2, ~ vs(line_name:env, Gu = GE)),
  # model3a = add_predictors(model1, ~ vs(line_name:env1, Gu = GE)),
  # model3b = model3a
  model4 = add_predictors(model1, ~vs(loc, Gu = L)),
  model5 = add_predictors(model4, ~vs(line_name:loc, Gu = GL))
) %>% map(~ formula(delete.response(terms(.)))) # Remove response

# Residual formula
resid_form <- ~ vs(units)



## Character vector of models
model_list <- names(model_rand_forms)




## Create a results df
full_model_df <- data_to_model %>% 
  mutate_at(vars(environment, location), list(fct = ~fct_contr_sum(as.factor(.)))) %>%
  rename(env = environment_fct, loc = location_fct) %>%
  group_by(trait) %>% 
  nest() %>%
  mutate(results = list(NULL))

## Loop over rows
for (i in seq_len(nrow(full_model_df))) {
# for (i in seq(i, nrow(full_model_df))) {
    
  
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
  E <- Env_mat(x = ec_mat[, main_environment_covariates, drop = FALSE], method = "Rincent2019")
  GE <- Env_mat(x = ec_mat[, interaction_environment_covariates, drop = FALSE], method = "Rincent2019") %>% 
    kronecker(X = K, Y = ., make.dimnames = TRUE)
  
  
  ## Create a new tibble with model formulas
  model_formula_df <- tibble(model = model_list, fixed = model_fixed_forms,
                             rand = model_rand_forms)
  model_formula_df <- bind_rows(
    filter(model_formula_df, ! model %in% c("model4", "model5")),
    crossing(filter(model_formula_df, model %in% c("model4", "model5")), time_frame = location_relmat_df$time_frame)
  )
  
  # Iterate over this df
  model_fits <- pmap(.l = model_formula_df, ~{
    fixed <- ..2
    rand <- ..3
    tf <- ..4
    
    if (!is.na(tf)) {
      L <- subset(location_relmat_df, trait == tr & time_frame == tf, E_mat_main, drop = TRUE)[[1]]
      GL <- subset(location_relmat_df, trait == tr & time_frame == tf, E_mat_int, drop = TRUE)[[1]] %>%
        kronecker(X = K, Y = ., make.dimnames = TRUE)
      
    }
    
    # Fit the model
    mmer(fixed = fixed, random = rand, rcov = resid_form, 
         data = data, date.warning = FALSE, verbose = TRUE)
    
  })
    
  
  ## Generate predictions
  model_predictions <- model_formula_df %>%
    mutate(object = model_fits) %>%
    mutate(R2 = map2_dbl(object, model, ~{
      
      mf <- model.frame(value ~ 1 + env + loc + line_name, data = .x$data)
      y <- model.response(mf)
      
      # Zg
      Zg <- model.matrix(~ -1 + line_name, mf)
      # Ze or Zl
      Ze <- if (.y %in% c("model2", "model3")) {
        model.matrix(~ -1 + env, mf) 
      } else if (.y %in% c("model4", "model5")) {
        model.matrix(~ -1 + loc, mf)
      } else {
        NULL
      }
      
      Zge <- if (.y %in% c("model3")) {
        model.matrix(~ -1 + line_name:env, mf)
      } else if (.y %in% c("model5")) {
        model.matrix(~ -1 + line_name:loc, mf)
      } else {
        NULL
      }
      
      
      # List of Zs
      Zs <- list(Zg, Ze, Zge) %>%
        subset(., map_lgl(., ~!is.null(.)))
      # List of ranefs
      ranefs <- map(.x$U, "value")
      
      ## Random predictions
      us <- map2(.x = Zs, .y = ranefs, ~.x %*% .y) %>%
        reduce(`+`)
      # Predict phenotypes
      y_hat <- .x$Beta$Estimate + us
      
      # Return R2
      cor(y, y_hat)^2
      
    }))
      

  ## Convert to tibble, extract diagonistics and varcomp
  model_fits_diag <- model_predictions %>%
    mutate(logLik = map(object, "monitor") %>% map_dbl(~last(.[1,])), # LogLik
           sigma = map(object, "sigma")) %>%
    select(-object)
  
  # Add to list
  full_model_df$results[[i]] <- model_fits_diag
  
}


## Unnest, add logLik
full_models <- full_model_df %>%
  unnest(results)



# Save the df
save("full_models", file = file.path(result_dir, "full_models.RData"))



  
