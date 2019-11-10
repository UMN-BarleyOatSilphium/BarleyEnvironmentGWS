## S2MET Prediction Models
## 
## Leave-one-environment-out prediction
## CV and POV
## 
## Author: Jeff Neyhart
## Last modified: 14 October 2019
## 

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(broom)

## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))



### Models

# Cross-validation will use the following models:
# M1 - random genotype, random environment
# M2 - random genotype (K), random environments
# M3 - random genotype (K), random environment, predicted covariate slopes
# M4 - random genotype (K), covariate-predicted environment, predicted covariate slopes



### Parent-offspring validation

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  ## Sample environments
  # filter(environment %in% sample_env) %>%
  ##
  # add covariates
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)),
         environment = as.factor(environment))


## 
## Leave-one-environment-out
## 

# Generate skeleton train/test sets
data_train_test <- data_to_model %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    
    ## Get the integer ids for the train/test rows
    test_id <- subset(df, line_name %in% c(tp_geno, vp_geno), id)
    train_id <- subset(data_to_model, line_name %in% tp_geno & environment != unique(df$environment) & trait == unique(df$trait), id)
    
    tibble(train = list(train_id), test = list(test_id))
    
  }) %>% ungroup()


## Iterate over trait-environment combinations
predictions_out <- data_train_test %>%
  mutate(out = list(NULL))
  

# Iterate over rows
for (i in seq(nrow(predictions_out))) {
  
  # Get training and test data
  row <- predictions_out[i,]
  train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
    mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
    select(environment, line_name, value)
  
  ## Get the covariate model
  ec_model_i <- subset(ec_model_building, trait == row$trait, final_model, drop = T)[[1]]
  # # Convert the model to fixed
  # ec_model_form <- as.formula(paste0("value ~ 1 + ", paste0(str_remove(attr(terms(formula(ec_model_i)), "term.labels"), "1 \\| "), collapse = " + ")))
  # # Remove environment
  # ec_model_form_no_env <- formula(drop.terms(terms(ec_model_form), 2, keep.response = TRUE))

  # Keep same model
  ec_model_form <- formula(ec_model_i)
  
  
  ## Get the environment effect model
  env_effect_model <- subset(env_effect_models, trait == row$trait, model, drop = TRUE)[[1]]
  
  ## All covariates that will be used
  env_covariate_list <- intersect(names(ec_tomodel), union(all.vars(ec_model_form), all.vars(formula(env_effect_model)))) %>%
    str_subset(., "environment", negate = TRUE)
  
  ## Add covariates to the train df
  train1 <- ec_tomodel %>%
    select(., environment, env_covariate_list) %>%
    left_join(train, ., by = "environment") %>%
    mutate_at(vars(line_name, environment), as.factor) %>%
    mutate_at(vars(line_name, environment), droplevels) %>%
    mutate_at(vars(line_name, environment), ~`contrasts<-`(., value = `colnames<-`(contr.sum(levels(.)), head(levels(.), -1)))) %>%
    ungroup()
  
  ## Extract value of covariates
  test_covariates <- ec_tomodel %>% 
    subset(., environment == unique(test$environment), intersect(env_covariate_list, all.vars(ec_model_form))) %>%
    unlist()
  
  

  
  ###################
  ### Model 1 - random genotype and environment
  ###################
  
  # Create the formula
  model1_form <- terms(ec_model_form) %>% 
    attr(., "term.labels") %>% 
    str_subset(., "\\|") %>% 
    str_c("(", ., ")", collapse = " + ") %>%
    str_c("~ ", .) %>%
    as.formula() %>%
    update.formula(old = ec_model_form, new = .)
  
  
  ## fit the model
  model1_fit <- lmer(formula = model1_form, data = train1)
  
  ## Predict - this should simply be the mean + line BLUP
  model1_prediction <- test %>%
    mutate(prediction = predict(object = model1_fit, newdata = test, allow.new.levels = TRUE))
  
  
  
  ###################
  ### Model 2 - random genotype (with markers) and environment
  ###################
  
  ## Get matrices ready
  y <- train1$value
  Zg <- model.matrix(~ -1 + line_name, train1)
  Ze <- model.matrix(~ -1 + environment, train1)
  
  ## Fit with sommer
  model2_fit <- mmer(fixed = value ~ 1, random = ~vs(line_name, Gu = K))
  
  
  
 
  
  
  
   
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  

  ## Fit a model to estimate genotypic and environmental effects
  main_effect_model <- lm(formula = value ~ line_name + environment, data = train1)
  
  ## Get the grand mean
  grand_mean <- coef(main_effect_model)[[1]]
  
  # Get the effects
  gen_env_effects <- tidy(main_effect_model) %>%
    mutate(group = str_extract(term, "line_name|environment")) %>%
    filter(!is.na(group)) %>%
    # Edit the term
    mutate(term = str_remove(term, group)) %>%
    split(.$group) %>%
    map(~add_row(., estimate = -sum(.$estimate), term = tail(levels(train1[[unique(.$group)]]), 1),
                 group = unique(.$group))) %>%
    map2(.x = ., .y = names(.), ~select(.x, term, estimate) %>% `names<-`(., c(.y, "effect")))
  
  # The genotypic effects will be blues for model1
  model1_blues <- gen_env_effects$line_name %>%
    mutate(variable = "genotype_mean") %>%
    rename(estimate = effect)
  

  ## Fit the covariate model using the training data
  ec_model_train <- lm(formula = ec_model_form_no_env, data = train1)

  # Get coefficients
  model_coef <- coef(ec_model_train)
  # Calculate the missing coefficient
  missing_coef <- model_coef %>% 
    subset(., str_detect(names(.), "line_name[A-Za-z0-9-]*$")) %>%
    {-sum(.)}
  

  ## Tidy the coefficients
  model2_blues <- tibble(term = names(model_coef), estimate = model_coef) %>%
    ## add the missing coefficients
    add_row(term = paste0("line_name", tail(levels(train2$line_name), 1)), estimate = missing_coef) %>%
    filter(str_detect(term, "line_name")) %>%
    mutate(term = str_remove_all(term, "line_name")) %>%
    select(term, estimate) %>%
    ## Create variables
    separate(data = ., col = term, into = c("line_name", "variable"), sep = ":", fill = "right") %>%
    mutate(variable = ifelse(is.na(variable), "genotype_mean", variable),
           line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  
  ## M1
  
  
  
  
  
  # First predict genotype-specific effect
  model1_blups <- model1_blues %>%
    mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
    group_by(variable) %>%
    do({
      mf <- .
      rr_fit <- mixed.solve(y = mf$estimate, Z = model.matrix(~ -1 + line_name, mf), K = K)
      # Return predictions
      rr_fit$u %>% 
        tibble(line_name = names(.), prediction = .)
    }) %>% ungroup() %>%
    spread(variable, prediction) %>%
    as.data.frame() %>%
    column_to_rownames("line_name") %>%
    as.matrix()
  
  ## Now predict 
  model1_prediction <- model1_blups %>%
    as.data.frame() %>%
    rownames_to_column("line_name") %>%
    rename(prediction = genotype_mean) %>%
    left_join(test, ., by = "line_name") %>%
    mutate(model = "model1", scaled_prediction = prediction + grand_mean)
  
  ## M2
  # First predict genotype-specific slopes
  model2_blups <- model2_blues %>%
    mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
    group_by(variable) %>%
    do({
      mf <- .
      rr_fit <- mixed.solve(y = mf$estimate, Z = model.matrix(~ -1 + line_name, mf), K = K)
      # Return predictions
      rr_fit$u %>% 
        tibble(line_name = names(.), prediction = .)
    }) %>% ungroup() %>%
    spread(variable, prediction) %>%
    as.data.frame() %>%
    column_to_rownames("line_name") %>%
    as.matrix()
  
  ## Now predict 
  model2_prediction <- {model2_blups %*% c(1, test_covariates)} %>%
    as.data.frame() %>%
    rownames_to_column("line_name") %>%
    rename(prediction = V1) %>%
    left_join(test, ., by = "line_name") %>%
    mutate(model = "model2", scaled_prediction = prediction + grand_mean)
  
  
  ## Model 3
  # First train the environment effect model on the training data
  env_effect_retrain <- gen_env_effects$environment %>%
    left_join(., distinct(select(train1, trial, environment, env_covariate_list)), by = "environment") %>%
    update(env_effect_model, data = .)
  # Predict the environment effect
  env_effect_prediction <- predict(object = env_effect_retrain, 
                                   newdata = subset(ec_tomodel, environment == unique(test$environment), env_covariate_list))
  
  # Add predicted environmental effect
  model3_prediction <- model2_prediction %>%
    mutate(scaled_prediction = scaled_prediction + env_effect_prediction,
           model = "model3")
  
  
  ## Return these results
  predictions_out$out[[i]] <- bind_rows(model1_prediction, model2_prediction, model3_prediction)
  
}


### Calculate and summarize accuracy
predictions_out %>%
  unnest(out) %>%
  group_by(trait, environment, model) %>%
  summarize(accuracy = cor(prediction, value))


