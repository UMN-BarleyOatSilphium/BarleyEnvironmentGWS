## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-environment-out prediction
## 
## Author: Jeff Neyhart
## Last modified: 8 January 2020
## 


# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# # Run the source script
# repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
# source(file.path(repo_dir, "source_MSI.R"))

## Number of cores
# n_core <- detectCores()
n_core <- 4


# Other packages
library(modelr)
library(broom)
library(parallel)



## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))
# Rename
ec_model_building <- unified_ec_models

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))
# Load the full model variance components
load(file.path(result_dir, "full_models.RData"))


### Models

# Cross-validation will use the following models:
# model1: y = G + r
# model2: y = G + E + r
# model2a: y = G + e + r
# model2b: y = G + eB + r
# model3: y = G + E + GE + r
# model3a: y = G + e + GE + r
# model3b: y = G + eB + GE + r

model_list <- unique(full_models$model)



### CV0 and POV0

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)),
         environment = as.factor(environment)) %>%
  ## Rename environment
  rename(env = environment)


## 
## Leave-one-environment-out
## 


# Generate skeleton train/test sets for LOEO
loeo_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo2(data = droplevels(group_by(., env)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno))))
  

## Assign cores and split
data_train_test1 <- loeo_train_test %>% 
  assign_cores(df = ., n_core = n_core, split = TRUE)


# Iterate over rows
# for (i in seq(nrow(loeo_predictions_out))) {
# for (i in seq(i, nrow(loeo_predictions_out))) {


## Parallelize
loeo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      # Get training and test data
      row <- core_df[i,]
      tr <- unique(row$trait)
      
      # train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
      #   mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
      # test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
      #   select(environment, line_name, value)
      
      ## Different subsetting method
      train <- row$train[[1]]
      test <- row$test[[1]]
      
      
      ## Get the covariate model
      ec_model_i <- ec_model_building %>% 
        unnest(final_model) %>%
        subset(trait == row$trait & model == "model3_ammi", object, drop = TRUE)
      
      # Extract the fitted model object formula
      ec_model_form <- formula(ec_model_i[[1]])
      
      # main_environment_covariates
      main_environment_covariates <- all.vars(ec_model_form) %>% 
        subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
                                             pattern = .))))
      
      # interaction_environment_covariates
      interaction_environment_covariates <- all.vars(ec_model_form) %>% 
        subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
                                             pattern = .)))) %>% setdiff(., "line_name")
      
      
      ## Add covariates to train and test
      train <- train %>% 
        left_join(., filter(ec_tomodel_scaled, environment %in% levels(train$env)) %>% 
                    select(env = environment, union(main_environment_covariates, interaction_environment_covariates)) %>% 
                    mutate(env = factor(env, levels = levels(train$env))))
      test <- test %>%
        left_join(., filter(ec_tomodel_scaled, environment %in% levels(test$env)) %>% 
                    select(env = environment, union(main_environment_covariates, interaction_environment_covariates)) %>% 
                    mutate(env = factor(env, levels = levels(test$env))))
      
      ## Subset covariates and convert to a matrix
      ec_mat <- ec_tomodel_scaled %>% 
        filter(environment %in% levels(train$env)) %>%
        select(environment, union(main_environment_covariates, interaction_environment_covariates)) %>% 
        as.data.frame() %>%
        column_to_rownames("environment") %>%
        as.matrix()
      
      
      
      ## Create relationship matrices
      K <- K # Genomic
      E <- Env_mat(x = ec_mat[,main_environment_covariates, drop = FALSE], method = "Rincent2019")
      GE <- Env_mat(x = ec_mat[,interaction_environment_covariates, drop = FALSE], method = "Rincent2019") %>%
        kronecker(X = K, Y = ., make.dimnames = TRUE)
      
      
      ## Create a list of model formulas
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
        model1 = ~ vs(line_name, Gu = K, Gt = varG, Gtc = fixm(1)),
        model2 = add_predictors(model1, ~ vs(env, Gu = E, Gt = varE, Gtc = fixm(1))),
        model2a = model1,
        model2b = model1,
        model3 = add_predictors(model2, ~ vs(line_name:env, Gu = GE, Gt = varGE, Gtc = fixm(1))),
        model3a = add_predictors(model1, ~ vs(line_name:env, Gu = GE, Gt = varGE, Gtc = fixm(1))),
        model3b = model3a
      ) %>% map(~ formula(delete.response(terms(.)))) # Remove response
      
      # Residual formula
      resid_form <- ~ vs(units, Gt = varR, Gtc = fixm(1))
      
      
      
      #################
      ## Fit models and extract predictions
      #################
      
      prediction_out <- full_models %>%
        filter(trait == tr) %>%
        mutate(prediction = list(NULL))
      
      # Iterate over models
      for (m in seq(nrow(prediction_out))) {
        
        # Model name and formulas
        mod <- prediction_out$model[m]
        fixed_form <- model_fixed_forms[[mod]]
        rand_form <- model_rand_forms[[mod]]
        
        ## Get the variance components from the full model
        full_varcomp <- subset(prediction_out, model == mod, sigma, drop = T)[[1]]
        # Assign values to separate objects
        varG <- full_varcomp$`u:line_name`
        varE <- full_varcomp$`u:env`
        varGE <- full_varcomp$`u:;line_name:env`
        varR <- full_varcomp$units
        
        ## Change contrasts, if necessary
        if ("env" %in% attr(terms(fixed_form), "term.labels")) {
          train1 <- train %>%
            mutate(env = droplevels(env))
          contrasts(train1$env) <- contr.sum(levels(train1$env)) %>% `colnames<-`(., head(levels(train1$env), -1))
          
          test1 <- test %>%
            mutate(env = droplevels(env))
          contrasts(test1$env) <- contr.sum(levels(test1$env)) %>% `colnames<-`(., head(levels(test1$env), -1))
          
          
          ## Create an X matrix for test
          Xtest <- model.matrix(fixed_form, test)
          
        } else {
          train1 <- train
          
          ## Create an X matrix for test
          Xtest <- model.matrix(fixed_form, test)
        
        }
        
        
        ## Fit the model
        model_fit <- mmer(fixed = fixed_form, random = rand_form, rcov = resid_form, 
                          data = train1, date.warning = FALSE, verbose = TRUE)
        
        ## Fixed effects
        fixed_eff <- coef(model_fit) %>%
          select(term = Effect, estimate = Estimate) %>%
          column_to_rownames("term") %>%
          as.matrix()
        
        fixed_pred <- Xtest %*% fixed_eff
        
        
        
        ## Random effects
        rand_eff <- randef(model_fit) %>%
          map("value")
        
        Ztest_list <- attr(terms(rand_form), "term.labels") %>% 
          # Extract the terms
          str_extract("\\(.*\\)") %>% 
          str_remove_all("\\(|\\)") %>% 
          str_split(., ", ") %>%
          map_chr(first) %>%
          # Put into a formula
          map(~reformulate(c("-1", .))) %>%
          map(~model.matrix(object = ., data = test))
        
        rand_pred <- reduce(map2(Ztest_list, rand_eff, `%*%`), `+`)
        
        # Add predictions to the list
        prediction_out$prediction[[m]] <- mutate(test, pred_complete = as.numeric(fixed_pred + rand_pred), 
                                                 pred_incomplete = as.numeric(rand_pred))
        
      }
      

      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- imap_dfr(prediction_out, ~tibble(model = .y, predictions = .x))
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in environment", as.character(row$environment), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()


## Save the results
save("loeo_predictions_out", file = file.path(result_dir, "loeo_predictions.RData"))


## Simple analysis
loeo_predictions_out %>% 
  unnest(predictions) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  group_by(trait, environment, pop, model) %>%
  summarize(acc = cor(prediction, value)) %>%
  ggplot(aes(x = pop, y = acc, color = model)) +
  geom_boxplot() +
  facet_grid(~ trait)






