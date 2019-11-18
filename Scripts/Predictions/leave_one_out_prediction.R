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


# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))

## Number of cores
n_core <- detectCores()
n_core <- 16


# Other packages
library(modelr)
library(broom)
library(parallel)



## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))



### Models

# Cross-validation will use the following models:
# M1 - random genotype (K), fixed environment
# M2 - random genotype (K), fixed environment by covariate
# M4 - random genotype (K), fixed environment by covariate, random regression slopes



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


# ## Iterate over trait-environment combinations
# loeo_predictions_out <- data_train_test %>%
#   # filter(trait %in% str_subset(traits, "Grain|Weight")) %>%
#   mutate(out = list(NULL))

## Assign cores and split
data_train_test1 <- data_train_test %>% 
  assign_cores(df = ., n_core = n_core) %>% 
  split(.$core)


# Iterate over rows
# for (i in seq(nrow(loeo_predictions_out))) {
# for (i in seq(i, nrow(loeo_predictions_out))) {
    

## Parallelize
loeo_predictions_out <- data_train_test1 %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      # Get training and test data
      row <- core_df[i,]
      
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
      ### Model 1 - random genotype (covariate) and fixed environment
      ###################
      
      ## Copy line_name as new factor / variable
      train2 <- bind_cols(train1, rename_all(as.data.frame(rerun(length(test_covariates), train1$line_name)), 
                                             ~paste0("line_name", seq_along(test_covariates)))) %>%
        mutate(environment = factor(environment, levels = c(as.character(row$environment), levels(train1$environment))),
               ge = interaction(line_name, environment, drop = FALSE, sep = ":"))
      
      
      # Random effect formula
      fixed_form1 <- value ~ 1 + environment
      random_form1 <- ~vs(line_name, Gu = K)
      
      ## Fit with sommer
      model1_fit <- mmer(fixed = fixed_form1, random = random_form1, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model1_fixed_prediction <- subset(model1_fit$Beta, Effect == "(Intercept)", Estimate, drop = TRUE)
      
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
      ### Model 2 - random genotype (K) + fixed environment based on covariate
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
      test_covariate_x <- t(as.matrix(ec_df[as.character(row$environment),, drop = FALSE]))
      
      
      # ## Calculate standardized difference between environments for each covariate
      # ec_dist_mat <- map(ec_df, ~{
      #   dist_mat <- as.matrix(dist(.x))
      #   dimnames(dist_mat) <- list(row.names(ec_df), row.names(ec_df))
      #   ## Standardize
      #   dist_mat / diff(range(.x))
      # })
      # 
      # ## Calculate the covariance matrix
      # ## This line will sum each of the same coordinate element in the list of matrices
      # E_mat <- 1 - reduce(ec_dist_mat, `+`)
      
      # Fixed formula
      fixed_form2 <- as.formula(paste0("value ~ 1 +", paste0(test_covariates, collapse = " + ")))
      # Random effect formula
      # random_form2 <- ~vs(line_name, Gu = K) + vs(environment, Gu = E_mat)
      random_form2 <- random_form1
      
      ## Fit with sommer
      model2_fit <- mmer(fixed = fixed_form2, random = random_form2, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model2_fixed_prediction <- c(model2_fit$Beta$Estimate %*% c(1, test_covariate_x))
      
      # Random effects
      model2_prediction <- model2_fit$U %>% 
        subset(., str_detect(names(.), "u:")) %>% 
        map("value") %>%
        map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
        reduce(crossing) %>%
        mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
        select(line_name, prediction) %>% 
        mutate(prediction = prediction + model2_fixed_prediction) %>%
        left_join(test, ., by = c("line_name"))
      
      
      
      # ###################
      # ### Model 3 - random genotype (K) + random environment (E) + random regression of K * E
      # ###################
      # 
      # 
      # 
      # ## Calculate GxE matrix
      # GE_mat <- kronecker(K, E_mat, make.dimnames = TRUE)
      # 
      # 
      # ## Add GxE term to formula
      # random_form3 <- add_predictors(random_form2, ~vs(line_name:environment, Gu = GE_mat))
      # 
      # ## Fit with sommer
      # model3_fit <- mmer(fixed = fixed_form1, random = random_form3, data = train2, date.warning = FALSE)
      # 
      # 
      # # Predict
      # # Fixed effects
      # model3_fixed_prediction <- model3_fit$Beta$Estimate[1]
      # 
      # # Random effects
      # model3_random_prediction <- model3_fit$U %>% 
      #   subset(., str_detect(names(.), "u:")) %>% 
      #   map("value") %>%
      #   map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), paste0(.y, "_estimate"))))
      # 
      # model3_prediction <- model3_random_prediction %>%
      #   map(~{
      #     if (str_detect(names(.x)[1], ":")) {
      #       separate(.x, names(.x)[1], str_split(names(.x)[1], ":")[[1]], sep = ":")
      #     } else {
      #       .x
      #     }
      #     
      #   }) %>%
      #   rev(.) %>%
      #   reduce(., .f = left_join) %>%
      #   mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
      #   select(line_name, environment, prediction) %>% 
      #   mutate(prediction = prediction + model3_fixed_prediction) %>%
      #   left_join(test, ., by = c("line_name", "environment"))
      # 
      # 
      # ###################
      # ### Model 4 - random genotype (K) + fixed environment covariate + random regression of covariate on K
      # ###################
      
      # Fixed formula
      fixed_form4 <- fixed_form2
      
      # Random effect formula
      # vs(x, y) specifies an interaction
      random_form4 <- formula(paste0("~ vs(line_name, Gu = K) + ", 
                                     paste0("vs(", test_covariates, ", line_name, Gu = K)", collapse = " + ")))
      
      ## Fit with sommer
      model4_fit <- mmer(fixed = fixed_form4, random = random_form4, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model4_fixed_prediction <- c(model4_fit$Beta$Estimate %*% c(1, test_covariate_x))
      
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
      
      
      
      
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- bind_rows(
        mutate(model1_prediction, model = "model1"),
        mutate(model2_prediction, model = "model2"),
        # mutate(model3_prediction, model = "model3"),
        mutate(model4_prediction, model = "model4")
      )
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in environment", as.character(row$environment), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()







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



# ## Iterate over trait-environment combinations
# loyo_predictions_out <- data_train_test %>%
#   # filter(trait %in% str_subset(traits, "Grain|Weight")) %>%
#   mutate(out = list(NULL))
# 
# # Iterate over rows
# for (i in seq(nrow(loyo_predictions_out))) {


## Assign cores and split
data_train_test1 <- data_train_test %>% 
  assign_cores(df = ., n_core = n_core) %>% 
  split(.$core)


# Iterate over rows
# for (i in seq(nrow(loeo_predictions_out))) {
# for (i in seq(i, nrow(loeo_predictions_out))) {


## Parallelize
loyo_predictions_out <- data_train_test1 %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      # Get training and test data
      row <- core_df[i,]

  
      # Get training and test data
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
      fixed_form1 <- value ~ 1 + environment
      random_form1 <- ~vs(line_name, Gu = K)
      
      ## Fit with sommer
      model1_fit <- mmer(fixed = fixed_form1, random = random_form1, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model1_fixed_prediction <- subset(model1_fit$Beta, Effect == "(Intercept)", Estimate, drop = TRUE)
      
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
      
      
      ## Create the environmental covariance matrix
      ec_df <- ec_tomodel %>% 
        filter(environment %in% levels(train2$environment)) %>% 
        select(environment, test_covariates) %>%
        as.data.frame() %>%
        column_to_rownames("environment")
      
      # Matrix of covariates for the test environment
      test_covariate_x <- t(as.matrix(ec_df[unique(test$environment),, drop = FALSE]))
    
      
      # ## Calculate standardized difference between environments for each covariate
      # ec_dist_mat <- map(ec_df, ~{
      #   dist_mat <- as.matrix(dist(.x))
      #   dimnames(dist_mat) <- list(row.names(ec_df), row.names(ec_df))
      #   ## Standardize
      #   dist_mat / diff(range(.x))
      # })
      # 
      # ## Calculate the covariance matrix
      # ## This line will sum each of the same coordinate element in the list of matrices
      # E_mat <- 1 - reduce(ec_dist_mat, `+`)
      
      
      # Fixed formula
      fixed_form2 <- as.formula(paste0("value ~ 1 +", paste0(test_covariates, collapse = " + ")))
      # Random effect formula
      # random_form2 <- ~vs(line_name, Gu = K) + vs(environment, Gu = E_mat)
      random_form2 <- random_form1
      
      ## Fit with sommer
      model2_fit <- mmer(fixed = fixed_form2, random = random_form2, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model2_fixed_prediction <- (cbind(mu = 1, t(test_covariate_x)) %*% model2_fit$Beta$Estimate) %>%
        as.data.frame() %>%
        rename(fixed_estimate = V1) %>%
        rownames_to_column("environment")
      
      # Random effects
      model2_prediction <- model2_fit$U %>% 
        subset(., str_detect(names(.), "u:")) %>% 
        map("value") %>%
        map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
        reduce(crossing) %>%
        crossing(., model2_fixed_prediction) %>%
        mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
        select(line_name, environment, prediction) %>% 
        left_join(test, ., by = c("line_name", "environment"))
      
      
      
      # ###################
      # ### Model 3 - random genotype (K) + random environment (E) + random regression of K * E
      # ###################
      # 
      # 
      # 
      # ## Calculate GxE matrix
      # GE_mat <- kronecker(K, E_mat, make.dimnames = TRUE)
      # 
      # 
      # ## Add GxE term to formula
      # random_form3 <- add_predictors(random_form2, ~vs(line_name:environment, Gu = GE_mat))
      # 
      # ## Fit with sommer
      # model3_fit <- mmer(fixed = value ~ 1, random = random_form3, data = train2, date.warning = FALSE)
      # 
      # 
      # # Predict
      # # Fixed effects
      # model3_fixed_prediction <- model3_fit$Beta$Estimate[1]
      # 
      # # Random effects
      # model3_random_prediction <- model3_fit$U %>% 
      #   subset(., str_detect(names(.), "u:")) %>% 
      #   map("value") %>%
      #   map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), paste0(.y, "_estimate"))))
      # 
      # model3_prediction <- model3_random_prediction %>%
      #   map(~{
      #     if (str_detect(names(.x)[1], ":")) {
      #       separate(.x, names(.x)[1], str_split(names(.x)[1], ":")[[1]], sep = ":")
      #     } else {
      #       .x
      #     }
      #     
      #   }) %>%
      #   rev(.) %>%
      #   reduce(., .f = left_join) %>%
      #   mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
      #   select(line_name, environment, prediction) %>% 
      #   mutate(prediction = prediction + model3_fixed_prediction) %>%
      #   left_join(test, ., by = c("line_name", "environment"))
      # 
      # 
      # ###################
      # ### Model 4 - random genotype (K) + random environment + random regression of covariate on K
      # ###################
      

      # Fixed formula
      fixed_form4 <- fixed_form2
      
      # Random effect formula
      # vs(x, y) specifies an interaction
      random_form4 <- formula(paste0("~ vs(line_name, Gu = K) + ", 
                                     paste0("vs(", test_covariates, ", line_name, Gu = K)", collapse = " + ")))
      
      ## Fit with sommer
      model4_fit <- mmer(fixed = fixed_form4, random = random_form4, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model4_fixed_prediction <- (cbind(mu = 1, t(test_covariate_x)) %*% model4_fit$Beta$Estimate) %>%
        as.data.frame() %>%
        rename(fixed_estimate = V1) %>%
        rownames_to_column("environment")
      
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
      model4_random_response_prediction <- (ec_Beta %*% test_covariate_x) %>%
        as.data.frame() %>%
        rownames_to_column("line_name") %>%
        gather(environment, response_estimate, -line_name)
      
      model4_prediction <- left_join(model4_random_prediction, model4_random_response_prediction, by = "line_name") %>%
        crossing(., model4_fixed_prediction) %>%
        mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
        select(line_name, environment, prediction) %>% 
        left_join(test, ., by = c("line_name", "environment"))
      
      
      
      
      ###################
      
      ###################
      
      
      
      
      ## Combine and return the predictions
      out[[i]] <- bind_rows(
        mutate(model1_prediction, model = "model1"),
        mutate(model2_prediction, model = "model2"),
        # mutate(model3_prediction, model = "model3"),
        mutate(model4_prediction, model = "model4")
      )
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in year", as.character(row$year), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()



## Save the results
save("loeo_predictions_out", "loyo_predictions_out", file = file.path(result_dir, "loo_predictions_out.RData"))
