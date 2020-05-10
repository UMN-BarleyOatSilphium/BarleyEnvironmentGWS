## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-location-out prediction
## 
## Author: Jeff Neyhart
## Last modified: 27 Feb 2020
## 


# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# # Run the source script
# repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
# source(file.path(repo_dir, "source_MSI.R"))


# Other packages
library(modelr)
library(broom)
library(parallel)



# time_frame to use for the location relationship matrix
time_frame_use <- "time_frame5"



## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))
load(file = file.path(result_dir, "historical_ec_model_building.RData"))



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



### CV0 only

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq_len(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno))) %>%
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
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno)) %>% 
                      mutate_at(vars(env, loc), fct_contr_sum) %>%
                      mutate(env1 = env, loc1 = loc)) )


## Use the covariates identified in the model building to fit factorial
## regression models
lolo_factorial_regression_out <- lolo_train_test %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(lolo_factorial_regression_out))) {
  
  row <- lolo_factorial_regression_out[i,]
  
  tr <- unique(row$trait)
  
  # Get the training and testing data
  train <- row$train[[1]] %>%
    droplevels() %>%
    mutate_at(vars(env1, loc1), fct_contr_sum)
  test <- row$test[[1]] %>%
    droplevels() %>%
    mutate(env1 = as.factor(NA), loc1 = as.factor(NA))
  
  # Record the number of environment and observations used for training
  train_n <- summarize(train, nEnv = n_distinct(environment), nObs = n())
  
  
  ## Get a list of environmental and historical covariates
  env_covariates_use <- filter(environmental_relmat_df, trait == tr) %>% 
    mutate(covariates_use = `names<-`(final_covariates, term)) %>% 
    pull(covariates_use)
  hist_env_covariates_use <- filter(location_relmat_df, trait == tr, time_frame == time_frame_use) %>% 
    mutate(covariates_use = `names<-`(final_covariates, term)) %>% 
    pull(covariates_use)
  
  ## Create dfs of the covariates
  env_covariate_df <- map(env_covariates_use, ~select(ec_tomodel_centered, environment, .x)) %>%
    imap(~`names<-`(.x, c("environment", paste0(.y, ".", names(.x)[-1])))) %>%
    reduce(left_join, by = "environment")
  location_covariate_df <- map(hist_env_covariates_use, 
                               ~select(historical_ec_tomodel_centered[[time_frame_use]], location, .x)) %>%
    imap(~`names<-`(.x, c("location", paste0(.y, ".", names(.x)[-1])))) %>%
    reduce(left_join, by = "location")
  
  ############################
  # Model fitting
  ############################
  
  ## Models 1-3 (main g, e, gxe; no covariates)
  model_formulas_123 <- formulas(
    .response = ~ value,
    model1 = ~ 1 + line_name,
    model2 = add_predictors(model1, ~ env1),
    model3 = add_predictors(model2, ~ line_name:env1)
  )
  
  # Fit models
  model_fits_123 <- fit_with(data = train, lm, model_formulas_123)
  # Predict
  model_predictions_123 <- map(model_fits_123, ~{
    pred <- predict(object = .x, newdata = test, type = "terms", terms = "line_name")
    mutate(test, pred = as.numeric(pred) + attr(pred, "constant"))
  })
  
  
  ## Models 2a and 3a - factorial regression
  model_formulas_2a3a <- formulas(
    .response = ~ value,
    model2a = add_predictors(model_formulas_123$model1, reformulate(str_subset(names(env_covariate_df), "main"))),
    model3a = add_predictors(model2a, reformulate(paste0("line_name:", str_subset(names(env_covariate_df), "int"))))
  )
  
  ## Add covariates to train/test df
  train1 <- left_join(train, env_covariate_df, by = "environment")
  test1 <- left_join(test, env_covariate_df, by = "environment")
  
  # Fit models
  model_fits_2a3a <- fit_with(data = train1, lm, model_formulas_2a3a)
  # Predict
  model_predictions_2a3a <- map(model_fits_2a3a, ~add_predictions(test1, model = .))
  

  ## Models 4-5 (l, gxl; no covariates)
  model_formulas_45 <- formulas(
    .response = ~ value,
    model4 = add_predictors(model_formulas_123$model1, ~ loc1),
    model5 = add_predictors(model4, ~ line_name:loc1)
  )
  
  # Fit models
  model_fits_45 <- fit_with(data = train, lm, model_formulas_45)
  # Predict
  model_predictions_45 <- map(model_fits_45, ~{
    pred <- predict(object = .x, newdata = test, type = "terms", terms = "line_name")
    mutate(test, pred = as.numeric(pred) + attr(pred, "constant"))
  })
  
  ## Models 4a and 5a - factorial regression
  model_formulas_4a5a <- formulas(
    .response = ~ value,
    model4a = add_predictors(model_formulas_123$model1, reformulate(str_subset(names(location_covariate_df), "main"))),
    model5a = add_predictors(model4a, reformulate(paste0("line_name:", str_subset(names(location_covariate_df), "int"))))
  )
  
  ## Add covariates to train/test df
  train1 <- left_join(train, location_covariate_df, by = "location")
  test1 <- left_join(test, location_covariate_df, by = "location")
  
  # Fit models
  model_fits_4a5a <- fit_with(data = train1, lm, model_formulas_4a5a)
  # Predict
  model_predictions_4a5a <- map(model_fits_4a5a, ~add_predictions(test1, model = .))
  
  
  ## Combine the prediction dfs and return
  predictions_df <- c(model_predictions_123, model_predictions_2a3a, model_predictions_45, model_predictions_4a5a) %>% 
    imap_dfr(~mutate(.x, model = .y)) %>% 
    select(environment, location, year, trait, model, line_name, value, pred)
  
  lolo_factorial_regression_out$out[[i]] <- predictions_df
  
}


## Calculate overall and within-environment accuracy
lolo_fact_reg_overall_acc <- lolo_factorial_regression_out %>%
  unnest(out) %>%
  group_by(trait, model) %>%
  summarize(acc = cor(value, pred))

lolo_fact_reg_withinEnv_acc <- lolo_factorial_regression_out %>%
  unnest(out) %>%
  group_by(trait, model, environment) %>%
  summarize(acc = cor(value, pred))

## Summarize and combine
lolo_factorial_regression_acc <- lolo_fact_reg_overall_acc %>%
  left_join(., lolo_fact_reg_withinEnv_acc %>% summarize_at(vars(acc), list(~mean, ~min, ~max)))

as.data.frame(lolo_factorial_regression_acc)

lolo_factorial_regression_acc %>%
  ggplot(aes(x = trait, fill = model)) +
  geom_bar(aes(y = acc), stat = "identity", position = position_dodge(0.9)) +
  geom_point(aes(y = mean), position = position_dodge(0.9))



## Plot
lolo_fact_reg_plot_list <- lolo_factorial_regression_out %>%
  unnest(out) %>%
  group_by(trait) %>%
  do(plot = {
    ggplot(., aes(x = pred, y = value, color = environment)) +
      geom_point(size = 0.5) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(trait ~ model)
  }) %>% ungroup()

# Combine
cowplot::plot_grid(plotlist = lolo_fact_reg_plot_list$plot, ncol = 1)



## Conclusion - the factorial regression approach is not useful; maintain current course!






