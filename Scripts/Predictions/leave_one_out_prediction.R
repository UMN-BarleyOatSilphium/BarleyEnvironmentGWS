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


# # Run the source script
# repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
# source(file.path(repo_dir, "source_MSI.R"))
# 
# ## Number of cores
# n_core <- detectCores()
# n_core <- 8



# Other packages
library(modelr)
library(broom)
# For random regression
library(orthopolynom)



## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))



### Models

# Cross-validation will use the following models:
# M1 - random genotype (K), random environment
# M2 - random genotype (K), random environment (E)
# M3 - random genotype (K), random environment (E), random gxe (GE)
# M4 - random genotype (K), random environment, random regression slopes



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
  # filter(trait %in% str_subset(traits, "Grain|Weight")) %>%
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
  
  ## All covariates that will be used
  env_covariate_list <- intersect(names(ec_tomodel), all.vars(ec_model_form)) %>%
    str_subset(., "environment", negate = TRUE)
  
  ## Add covariates to the train df
  train1 <- ec_tomodel %>%
    select(., environment, env_covariate_list) %>%
    left_join(mutate(train, environment = as.character(environment)), ., by = "environment") %>%
    mutate(environment = as.factor(environment),
           line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  
  test1 <- ec_tomodel %>%
    select(., environment, env_covariate_list) %>%
    left_join(mutate(test, environment = as.character(environment)), ., by = "environment") %>%
    mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  
  ## Extract value of covariates
  test_covariates <- env_covariate_list
  

  
  
  
  
  ###################
  ### Model 1 - random genotype (covariate) and environment
  ###################
  
  ## Copy line_name as new factor / variable
  train2 <- bind_cols(train1, rename_all(as.data.frame(rerun(length(test_covariates), train1$line_name)), 
                                         ~paste0("line_name", seq_along(test_covariates)))) %>%
    mutate(environment = factor(environment, levels = c(as.character(row$environment), levels(train1$environment))),
           ge = interaction(line_name, environment, drop = FALSE, sep = ":"))
  
  
  # Random effect formula
  random_form1 <- ~vs(line_name, Gu = K) + environment
  
  ## Fit with sommer
  model1_fit <- mmer(fixed = value ~ 1, random = random_form1, data = train2, date.warning = FALSE)
  
  
  # Predict
  # Fixed effects
  model1_fixed_prediction <- model1_fit$Beta$Estimate
  
  # Random effects
  model1_prediction <- model1_fit$U %>% 
    subset(., str_detect(names(.), "u:")) %>% 
    map("value") %>%
    map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
    reduce(crossing) %>%
    mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
    select(line_name, prediction) %>% 
    mutate(prediction = prediction + model1_fixed_prediction) %>%
    left_join(test, ., by = c("line_name"))
  
  
  
  ###################
  ### Model 2 - random genotype (K) + random environment (E)
  ###################
  
  ## Copy line_name as new factor / variable
  train2 <- bind_cols(train1, rename_all(as.data.frame(rerun(length(test_covariates), train1$line_name)), 
                                         ~paste0("line_name", seq_along(test_covariates)))) %>%
    mutate(environment = factor(environment, levels = c(as.character(row$environment), levels(train1$environment))),
           ge = interaction(line_name, environment, drop = FALSE, sep = ":"))
  
  
  ## Create the environmental covariance matrix
  ec_df <- ec_tomodel %>% 
    filter(environment %in% levels(train2$environment)) %>% 
    select(environment, test_covariates) %>%
    as.data.frame() %>%
    column_to_rownames("environment")
  
  
  # Matrix of covariates for the test environment
  test_covariate_x <- t(as.matrix(ec_df[row$environment,, drop = FALSE]))

  
  ## Calculate standardized difference between environments for each covariate
  ec_dist_mat <- map(ec_df, ~{
    dist_mat <- as.matrix(dist(.x))
    dimnames(dist_mat) <- list(row.names(ec_df), row.names(ec_df))
    ## Standardize
    dist_mat / diff(range(.x))
  })
  
  ## Calculate the covariance matrix
  ## This line will sum each of the same coordinate element in the list of matrices
  E_mat <- 1 - reduce(ec_dist_mat, `+`)
  
  
  
  # Random effect formula
  random_form2 <- ~vs(line_name, Gu = K) + vs(environment, Gu = E_mat)
  
  ## Fit with sommer
  model2_fit <- mmer(fixed = value ~ 1, random = random_form2, data = train2, date.warning = FALSE)
  
  
  # Predict
  # Fixed effects
  model2_fixed_prediction <- model2_fit$Beta$Estimate
  
  # Random effects
  model2_prediction <- model2_fit$U %>% 
    subset(., str_detect(names(.), "u:")) %>% 
    map("value") %>%
    map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
    reduce(crossing) %>%
    mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
    select(line_name, environment, prediction) %>% 
    mutate(prediction = prediction + model2_fixed_prediction) %>%
    left_join(test, ., by = c("line_name", "environment"))
  
  
  
  ###################
  ### Model 3 - random genotype (K) + random environment (E) + random regression of K * E
  ###################


  
  ## Calculate GxE matrix
  GE_mat <- kronecker(K, E_mat, make.dimnames = TRUE)

  
  ## Add GxE term to formula
  random_form3 <- add_predictors(random_form2, ~vs(line_name:environment, Gu = GE_mat))

  ## Fit with sommer
  model3_fit <- mmer(fixed = value ~ 1, random = random_form3, data = train2, date.warning = FALSE)
  
  
  # Predict
  # Fixed effects
  model3_fixed_prediction <- model3_fit$Beta$Estimate[1]
  
  # Random effects
  model3_random_prediction <- model3_fit$U %>% 
    subset(., str_detect(names(.), "u:")) %>% 
    map("value") %>%
    map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), paste0(.y, "_estimate"))))
  
  model3_prediction <- model3_random_prediction %>%
    map(~{
      if (str_detect(names(.x)[1], ":")) {
        separate(.x, names(.x)[1], str_split(names(.x)[1], ":")[[1]], sep = ":")
      } else {
        .x
      }
      
    }) %>%
    rev(.) %>%
    reduce(., .f = left_join) %>%
    mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
    select(line_name, environment, prediction) %>% 
    mutate(prediction = prediction + model3_fixed_prediction) %>%
    left_join(test, ., by = c("line_name", "environment"))

  
  ###################
  ### Model 4 - random genotype (K) + random environment + random regression of covariate on K
  ###################
  
  # Random effect formula
  # vs(x, y) specifies an interaction
  random_form4 <- formula(paste0("~ vs(line_name, Gu = K) + environment + ", 
                                 paste0("vs(", test_covariates, ", line_name, Gu = K)", collapse = " + ")))
                                        
  ## Fit with sommer
  model4_fit <- mmer(fixed = value ~ 1, random = random_form4, data = train2, date.warning = FALSE)
  

  # Predict
  # Fixed effects
  model4_fixed_prediction <- model4_fit$Beta$Estimate[1]
  
  # Random effects
  model4_random_prediction <- model4_fit$U %>% 
    subset(., str_detect(names(.), ":")) %>% 
    map("value") %>%
    map(~tibble(term = names(.x), estimate = .x)) %>%
    imap(~`names<-`(.x, rev(str_split(.y, pattern = ":")[[1]]))) %>%
    reduce(left_join, by = "line_name") %>%
    rename_at(vars(u), ~"genotype_mean_estimate")
  
  ## Convert random regression coefficients to matrix
  ec_Beta <- model4_random_prediction %>%
    select(-genotype_mean_estimate) %>%
    as.data.frame(.) %>% 
    column_to_rownames("line_name") %>% 
    as.matrix()
  
  ## Predict the genotype-specific value for the target environment
  model4_random_prediction1 <- model4_random_prediction %>%
    mutate(response_estimate = c(ec_Beta %*% test_covariate_x))
  
  model4_prediction <- model4_random_prediction1 %>%
    mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
    select(line_name, prediction) %>% 
    mutate(prediction = prediction + model4_fixed_prediction) %>%
    left_join(test, ., by = c("line_name"))
  
  
  
  # ###################
  # ### Model 5 - random genotype (K) + covariate + random regression of covariate on K
  # ###################
  # 
  # # fixed effect formula
  # fixed_form5 <- formula(paste0("value ~ 1 + ", paste0(test_covariates, collapse = " + ")))
  # 
  # # Random effect formula
  # # vs(x, y) specifies an interaction
  # random_form5 <- formula(paste0("~ vs(line_name, Gu = K) + ", paste0("vs(", test_covariates, ", line_name, Gu = K)", collapse = " + ")))
  # 
  # ## Fit with sommer
  # model5_fit <- mmer(fixed = fixed_form5, random = random_form5, data = train2, date.warning = FALSE)
  # 
  # 
  # # Predict
  # # Fixed effects
  # model5_fixed_predictions_all <- model5_fit$Beta$Estimate
  # # Calculate grand mean and environmental main effect
  # model5_fixed_prediction <- c(model5_fixed_predictions_all %*% c(1, test_covariate_x))
  # 
  # 
  # # Random effects
  # model5_random_prediction <- model5_fit$U %>% 
  #   subset(., str_detect(names(.), ":")) %>% 
  #   map("value") %>%
  #   map(~tibble(term = names(.x), estimate = .x)) %>%
  #   imap(~`names<-`(.x, rev(str_split(.y, pattern = ":")[[1]]))) %>%
  #   reduce(left_join, by = "line_name") %>%
  #   rename_at(vars(u), ~"genotype_mean_estimate")
  # 
  # ## Convert random regression coefficients to matrix
  # ec_Beta <- model5_random_prediction %>%
  #   select(-genotype_mean_estimate) %>%
  #   as.data.frame(.) %>% 
  #   column_to_rownames("line_name") %>% 
  #   as.matrix()
  # 
  # ## Predict the genotype-specific value for the target environment
  # model5_random_prediction1 <- model5_random_prediction %>%
  #   mutate(response_estimate = c(ec_Beta %*% test_covariate_x))
  # 
  # model5_prediction <- model5_random_prediction1 %>%
  #   mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
  #   select(line_name, prediction) %>% 
  #   mutate(prediction = prediction + model5_fixed_prediction) %>%
  #   left_join(test, ., by = c("line_name"))

  
  ###################
  
  ###################

  ## Combine and return the predictions
  predictions_out$out[[i]] <- bind_rows(
    mutate(model1_prediction, model = "model1"),
    mutate(model2_prediction, model = "model2"),
    mutate(model3_prediction, model = "model3"),
    mutate(model4_prediction, model = "model4")
    )
  
  
  ## Notify user
  cat("\nPredictions for trait", row$trait, "in environment", as.character(row$environment), "complete.")
  
}


## Unnest
loeo_predictions_df <- unnest(predictions_out, out)



  




####################################
### Leave-one-year-out
####################################


# Generate skeleton train/test sets
data_train_test <- data_to_model %>%
  group_by(trait, year) %>%
  do({
    df <- .
    
    ## Get the integer ids for the train/test rows
    test_id <- subset(df, line_name %in% c(tp_geno, vp_geno), id)
    train_id <- subset(data_to_model, line_name %in% tp_geno & year != unique(df$year) & trait == unique(df$trait), id)
    
    tibble(train = list(train_id), test = list(test_id))
    
  }) %>% ungroup()



## Iterate over trait-environment combinations
predictions_out <- data_train_test %>%
  # filter(trait %in% str_subset(traits, "Grain|Weight")) %>%
  mutate(out = list(NULL))


# Iterate over rows
for (i in seq(nrow(predictions_out))) {
# for (i in seq(i, nrow(predictions_out))) {
  
  
  # Get training and test data
  row <- predictions_out[i,]
  test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
    select(environment, line_name, value) %>%
    mutate_at(vars(environment, line_name), as.character)
  train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
    mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  
  ## Get the covariate model
  ec_model_i <- subset(ec_model_building, trait == row$trait, final_model, drop = T)[[1]]
  # # Convert the model to fixed
  # ec_model_form <- as.formula(paste0("value ~ 1 + ", paste0(str_remove(attr(terms(formula(ec_model_i)), "term.labels"), "1 \\| "), collapse = " + ")))
  # # Remove environment
  # ec_model_form_no_env <- formula(drop.terms(terms(ec_model_form), 2, keep.response = TRUE))
  
  # Keep same model
  ec_model_form <- formula(ec_model_i)
  
  ## All covariates that will be used
  env_covariate_list <- intersect(names(ec_tomodel), all.vars(ec_model_form)) %>%
    str_subset(., "environment", negate = TRUE)
  
  ## Add covariates to the train df
  train1 <- ec_tomodel %>%
    select(., environment, env_covariate_list) %>%
    left_join(train, ., by = "environment") %>%
    mutate(environment = factor(environment, levels = c(unique(.$environment), unique(test$environment))),
           line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  
  test1 <- ec_tomodel %>%
    select(., environment, env_covariate_list) %>%
    left_join(test, ., by = "environment") %>%
    mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  
  ## Extract value of covariates
  test_covariates <- env_covariate_list
  
  
  
  
  ###################
  ### Model 1 - random genotype (covariate) and environment
  ###################
  
  ## Copy line_name as new factor / variable
  train2 <- bind_cols(train1, rename_all(as.data.frame(rerun(length(test_covariates), train1$line_name)), 
                                         ~paste0("line_name", seq_along(test_covariates)))) %>%
    mutate(environment = factor(environment, levels = c(as.character(row$environment), levels(train1$environment))),
           ge = interaction(line_name, environment, drop = FALSE, sep = ":"))
  
  
  # Random effect formula
  random_form1 <- ~vs(line_name, Gu = K) + environment
  
  ## Fit with sommer
  model1_fit <- mmer(fixed = value ~ 1, random = random_form1, data = train2,
                     verbose = FALSE, date.warning = FALSE)
  
  
  # Predict
  # Fixed effects
  model1_fixed_prediction <- model1_fit$Beta$Estimate
  
  # Random effects
  model1_prediction <- model1_fit$U %>% 
    subset(., str_detect(names(.), "u:")) %>% 
    map("value") %>%
    map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
    reduce(crossing) %>%
    mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
    select(line_name, prediction) %>% 
    mutate(prediction = prediction + model1_fixed_prediction) %>%
    left_join(test, ., by = c("line_name"))
  
  
  
  ###################
  ### Model 2 - random genotype (with markers) and environment (with covariates)
  ###################
  
  ## Copy line_name as new factor / variable
  train2 <- bind_cols(train1, rename_all(as.data.frame(rerun(length(test_covariates), train1$line_name)), 
                                         ~paste0("line_name", seq_along(test_covariates)))) %>%
    mutate(ge = interaction(line_name, environment, drop = FALSE, sep = ":"))
  
  ## Create the environmental covariance matrix
  ec_df <- ec_tomodel %>% 
    filter(environment %in% levels(train2$environment)) %>% 
    select(environment, test_covariates) %>%
    as.data.frame() %>%
    column_to_rownames("environment")
  
  ## Calculate standardized difference between environments for each covariate
  ec_dist_mat <- map(ec_df, ~{
    dist_mat <- as.matrix(dist(.x))
    dimnames(dist_mat) <- list(row.names(ec_df), row.names(ec_df))
    ## Standardize
    dist_mat / diff(range(.x))
  })
  
  ## Calculate the covariance matrix
  ## This line will sum each of the same coordinate element in the list of matrices
  E_mat <- 1 - reduce(ec_dist_mat, `+`)
  
  
  
  # Random effect formula
  random_form2 <- ~vs(line_name, Gu = K) + vs(environment, Gu = E_mat)
  
  ## Fit with sommer
  model2_fit <- mmer(fixed = value ~ 1, random = random_form2, data = train2,
                     verbose = FALSE, date.warning = FALSE)
  
  
  # Predict
  # Fixed effects
  model2_fixed_prediction <- model2_fit$Beta$Estimate
  
  # Random effects
  model2_prediction <- model2_fit$U %>% 
    subset(., str_detect(names(.), "u:")) %>% 
    map("value") %>%
    map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
    reduce(crossing) %>%
    mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
    select(line_name, environment, prediction) %>% 
    mutate(prediction = prediction + model2_fixed_prediction) %>%
    left_join(test, ., by = c("line_name", "environment"))
  
  ###################
  ### Model 3 - random genotype (with markers) and environment 
  ### with random GxE (with marker x covariate)
  ###################
  
  
  
  ## Calculate GxE matrix
  GE_mat <- kronecker(K, E_mat, make.dimnames = TRUE)
  
  
  ## Add GxE term to formula
  random_form3 <- add_predictors(random_form2, ~vs(line_name:environment, Gu = GE_mat))
  
  ## Fit with sommer
  model3_fit <- mmer(fixed = value ~ 1, random = random_form3, data = train2, verbose = FALSE,
                     date.warning = FALSE)
  
  
  # Predict
  # Fixed effects
  model3_fixed_prediction <- model3_fit$Beta$Estimate[1]
  
  # Random effects
  model3_random_prediction <- model3_fit$U %>% 
    subset(., str_detect(names(.), "u:")) %>% 
    map("value") %>%
    map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), paste0(.y, "_estimate"))))
  
  model3_prediction <- model3_random_prediction %>%
    map(~{
      if (str_detect(names(.x)[1], ":")) {
        separate(.x, names(.x)[1], str_split(names(.x)[1], ":")[[1]], sep = ":")
      } else {
        .x
      }
      
    }) %>%
    rev(.) %>%
    reduce(., .f = left_join) %>%
    mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
    select(line_name, environment, prediction) %>% 
    mutate(prediction = prediction + model3_fixed_prediction) %>%
    left_join(test, ., by = c("line_name", "environment"))
  
  
  
  ###################
  ### Model 4 - random genotype (with markers) and environment with individually predicted covariate slopes
  ###################
  
  ## Reset contrasts
  train2 <- train1 %>%
    mutate_at(vars(line_name, environment), as.factor) %>%
    mutate_at(vars(line_name, environment), droplevels) %>%
    mutate_at(vars(line_name, environment), ~`contrasts<-`(., value = `colnames<-`(contr.sum(levels(.)), head(levels(.), -1)))) %>%
    ungroup()
  
  
  ## Refit the mixed model
  ec_model_i_refit <- lmer(formula(ec_model_i), data = train2)
  
  ## Get the ranefs for line_name, add intercept, these are BLUPs for prediction
  genotype_mean_blups <- (ranef(ec_model_i_refit)$line_name + fixef(ec_model_i_refit)[1]) %>%
    rownames_to_column("line_name") %>%
    rename(genotype_mean = `(Intercept)`)
  
  ## Get the slope coefficients
  genotype_slope_blues <- fixef(ec_model_i_refit) %>%
    tibble(term = names(.), estimate = .) %>% 
    filter(str_detect(term, "line_name")) %>% 
    mutate(term = str_remove_all(term, "line_name")) %>% 
    separate(term, c("covariate", "line_name"), sep = ":") %>%
    spread(covariate, estimate) %>%
    add_row(line_name = tail(levels(train2$line_name), 1))
  
  genotype_slope_blues[nrow(genotype_slope_blues),-1] <- map_dbl(genotype_slope_blues[-1], ~-sum(., na.rm = TRUE))
  
  
  ## Create a model frame
  model4_mf <- inner_join(genotype_mean_blups, genotype_slope_blues, by = "line_name") %>%
    mutate(line_name = factor(line_name, levels = levels(train1$line_name)))
  
  Zg <- model.matrix(~ -1 + line_name, model4_mf)
  
  ## Run genomic prediction on both main blups and slopes
  model4_fit <- apply(X = model4_mf[,-1,drop = FALSE], MARGIN = 2, FUN = function(y) {
    fit <- mixed.solve(y = y, Z = Zg, K = K)$u
  })
  
  
  ## Predictions
  covariate_newdata <- test1 %>% 
    select(environment, test_covariates) %>% 
    distinct() %>% 
    gather(covariate, value, -environment) %>% 
    spread(environment, value) %>% 
    column_to_rownames("covariate") %>%
    as.matrix()
  
  # fixed effects
  model4_fixed_prediction <- fixef(ec_model_i_refit)[[1]]
  
  # Intercept (genotypic value)
  model4_prediction_genotype_effect <- cbind(estimate = model4_fit[,"genotype_mean"]) %>%
    as.data.frame() %>%
    rownames_to_column("line_name")
  
  ## Environmental effect
  model4_prediction_environment_effect <- (fixef(ec_model_i_refit)[test_covariates] %*% covariate_newdata) %>% 
    as.data.frame() %>% 
    gather(environment, estimate)
  
  # Reaction to covariate
  model4_prediction_slope <- model4_fit[,-1,drop = FALSE]
  
  ## Predictions of genotyoe-specific effects in each environment
  model4_prediction_random <- (model4_prediction_slope %*% covariate_newdata) %>% 
    as.data.frame() %>% 
    rownames_to_column("line_name") %>%
    gather(environment, estimate, -line_name)
  
  
  ## Combine to predict
  model4_prediction <- model4_prediction_random %>%
    left_join(., model4_prediction_genotype_effect, by = "line_name") %>% 
    left_join(., model4_prediction_environment_effect, by = "environment") %>%
    mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
    select(line_name, environment, prediction) %>% 
    mutate(prediction = prediction + model4_fixed_prediction) %>%
    left_join(test, ., by = c("line_name", "environment"))
  
  
  ## Combine and return the predictions
  predictions_out$out[[i]] <- bind_rows(
    mutate(model1_prediction, model = "model1"),
    mutate(model2_prediction, model = "model2"),
    mutate(model3_prediction, model = "model3"),
    mutate(model4_prediction, model = "model4")
  )
  
  
  ## Notify user
  cat("\nPredictions for trait", row$trait, "in year", as.character(row$year), "complete.")
  
}


## Unnest
loyo_predictions_df <- unnest(predictions_out, out)





## Save the results
save("loeo_predictions_df", "loyo_predictions_df", file = file.path(result_dir, "loo_predictions_out.RData"))











  
  
  
  
  
  
  
  
  
  
  
  
  
  

# # Iterate over rows
# for (i in seq(nrow(predictions_out))) {
#   
#   # Get training and test data
#   row <- predictions_out[i,]
#   train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
#     mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
#   test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
#     select(environment, line_name, value)
#   
#   ## Get the covariate model
#   ec_model_i <- subset(ec_model_building, trait == row$trait, final_model, drop = T)[[1]]
#   # Convert the model to fixed
#   ec_model_form <- as.formula(paste0("value ~ 1 + ", paste0(str_remove(attr(terms(formula(ec_model_i)), "term.labels"), "1 \\| "), collapse = " + ")))
#   # Remove environment
#   ec_model_form_no_env <- formula(drop.terms(terms(ec_model_form), 2, keep.response = TRUE))
#   
#   ## Get the environment effect model
#   env_effect_model <- subset(env_effect_models, trait == row$trait, model, drop = TRUE)[[1]]
#   
#   ## All covariates that will be used
#   env_covariate_list <- intersect(names(ec_tomodel), union(all.vars(ec_model_form), all.vars(formula(env_effect_model)))) %>%
#     str_subset(., "environment", negate = TRUE)
#   
#   ## Add covariates to the train df
#   train1 <- ec_tomodel %>%
#     select(., environment, env_covariate_list) %>%
#     left_join(train, ., by = "environment") %>%
#     mutate_at(vars(line_name, environment), as.factor) %>%
#     mutate_at(vars(line_name, environment), droplevels) %>%
#     mutate_at(vars(line_name, environment), ~`contrasts<-`(., value = `colnames<-`(contr.sum(levels(.)), head(levels(.), -1))))
#   
#   ## Extract value of covariates
#   test_covariates <- ec_tomodel %>% 
#     subset(., environment == unique(test$environment), intersect(env_covariate_list, all.vars(ec_model_form))) %>%
#     unlist()
# 
#   ## Fit a model to estimate genotypic and environmental effects
#   main_effect_model <- lm(formula = value ~ line_name + environment, data = train1)
#   
#   ## Get the grand mean
#   grand_mean <- coef(main_effect_model)[[1]]
#   
#   # Get the effects
#   gen_env_effects <- tidy(main_effect_model) %>%
#     mutate(group = str_extract(term, "line_name|environment")) %>%
#     filter(!is.na(group)) %>%
#     # Edit the term
#     mutate(term = str_remove(term, group)) %>%
#     split(.$group) %>%
#     map(~add_row(., estimate = -sum(.$estimate), term = tail(levels(train1[[unique(.$group)]]), 1),
#                  group = unique(.$group))) %>%
#     map2(.x = ., .y = names(.), ~select(.x, term, estimate) %>% `names<-`(., c(.y, "effect")))
#   
#   # The genotypic effects will be blues for model1
#   model1_blues <- gen_env_effects$line_name %>%
#     mutate(variable = "genotype_mean") %>%
#     rename(estimate = effect)
#   
# 
#   ## Fit the covariate model using the training data
#   ec_model_train <- lm(formula = ec_model_form_no_env, data = train1)
# 
#   # Get coefficients
#   model_coef <- coef(ec_model_train)
#   # Calculate the missing coefficient
#   missing_coef <- model_coef %>% 
#     subset(., str_detect(names(.), "line_name[A-Za-z0-9-]*$")) %>%
#     {-sum(.)}
#   
# 
#   ## Tidy the coefficients
#   model2_blues <- tibble(term = names(model_coef), estimate = model_coef) %>%
#     ## add the missing coefficients
#     add_row(term = paste0("line_name", tail(levels(train1$line_name), 1)), estimate = missing_coef) %>%
#     filter(str_detect(term, "line_name")) %>%
#     mutate(term = str_remove_all(term, "line_name")) %>%
#     select(term, estimate) %>%
#     ## Create variables
#     separate(data = ., col = term, into = c("line_name", "variable"), sep = ":", fill = "right") %>%
#     mutate(variable = ifelse(is.na(variable), "genotype_mean", variable),
#            line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
#   
#   ## M1
#   # First predict genotype-specific effect
#   model1_blups <- model1_blues %>%
#     mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
#     group_by(variable) %>%
#     do({
#       mf <- .
#       rr_fit <- mixed.solve(y = mf$estimate, Z = model.matrix(~ -1 + line_name, mf), K = K)
#       # Return predictions
#       rr_fit$u %>% 
#         tibble(line_name = names(.), prediction = .)
#     }) %>% ungroup() %>%
#     spread(variable, prediction) %>%
#     as.data.frame() %>%
#     column_to_rownames("line_name") %>%
#     as.matrix()
#   
#   ## Now predict 
#   model1_prediction <- model1_blups %>%
#     as.data.frame() %>%
#     rownames_to_column("line_name") %>%
#     rename(prediction = genotype_mean) %>%
#     left_join(test, ., by = "line_name") %>%
#     mutate(model = "model1", scaled_prediction = prediction + grand_mean)
#   
#   ## M2
#   # First predict genotype-specific slopes
#   model2_blups <- model2_blues %>%
#     mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
#     group_by(variable) %>%
#     do({
#       mf <- .
#       rr_fit <- mixed.solve(y = mf$estimate, Z = model.matrix(~ -1 + line_name, mf), K = K)
#       # Return predictions
#       rr_fit$u %>% 
#         tibble(line_name = names(.), prediction = .)
#     }) %>% ungroup() %>%
#     spread(variable, prediction) %>%
#     as.data.frame() %>%
#     column_to_rownames("line_name") %>%
#     as.matrix()
#   
#   ## Now predict 
#   model2_prediction <- {model2_blups %*% c(1, test_covariates)} %>%
#     as.data.frame() %>%
#     rownames_to_column("line_name") %>%
#     rename(prediction = V1) %>%
#     left_join(test, ., by = "line_name") %>%
#     mutate(model = "model2", scaled_prediction = prediction + grand_mean)
#   
#   
#   ## Model 3
#   # First train the environment effect model on the training data
#   env_effect_retrain <- gen_env_effects$environment %>%
#     left_join(., distinct(select(train1, trial, environment, env_covariate_list)), by = "environment") %>%
#     update(env_effect_model, data = .)
#   # Predict the environment effect
#   env_effect_prediction <- predict(object = env_effect_retrain, 
#                                    newdata = subset(ec_tomodel, environment == unique(test$environment), env_covariate_list))
#   
#   # Add predicted environmental effect
#   model3_prediction <- model2_prediction %>%
#     mutate(scaled_prediction = scaled_prediction + env_effect_prediction,
#            model = "model3")
#   
#   
#   ## Return these results
#   predictions_out$out[[i]] <- bind_rows(model1_prediction, model2_prediction, model3_prediction)
#   
# }





