## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-location-out prediction
## 
## Author: Jeff Neyhart
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

# Source of covariates
source_use <- "daymet"

## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

## Load the factorial regression results
load(file.path(result_dir, "feature_selection_results.RData"))

## Data.frame of timeframes to use per trait
time_frame_use_df <- bind_rows(historical_feature_selection, historical_feature_importance) %>%
  mutate(feat_sel_type = str_remove(feat_sel_type, "_nosoil")) %>%
  distinct(trait, time_frame, feat_sel_type)

## For either environments or locations, fit the models:
## 
## model1: G + r
## model2: G + E (L) + r
## model3: G + E (L) + GE (GL) + r
## 


## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value,
  model1 = ~ 1,
  model2_cov = model1,
  model2_id = model1,
  model3_cov = model1,
  model3_id = model1
)

## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model1 = ~ vs(line_name, Gu = K),
  model2_cov = add_predictors(model1, ~ vs(site1, Gu = E)),
  model2_id = add_predictors(model1, ~ vs(site1, Gu = I)),
  model3_cov = add_predictors(model2_cov, ~ vs(line_name:site1, Gu = GE)),
  model3_id = add_predictors(model2_id, ~ vs(line_name:site1, Gu = GI))
) %>% map(~ formula(delete.response(terms(.)))) # Remove response


# Combine into list
model.list <- list(fixed = model_fixed_forms, random = model_rand_forms)






# Locations ================================================


# Data to use - location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

## Pull out location covariates to use
loc_ec_tomodel_scaled <- bind_rows(historical_ec_tomodel_timeframe_scaled, historical_ec_tomodel_window_scaled) %>%
  filter(source == source_use, time_frame %in% unique(time_frame_use_df$time_frame))



# Leave-one-out -----------------------------------------------
# Data to use
data_to_model <- S2_MET_loc_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         location %in% train_test_loc) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(location), as.factor) %>%
  mutate(loc = location)


# Remove soil variable from all variables to create a "no soil" group
historical_all_features <- historical_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(historical_all_features, .)

# Combine feature selection df
feature_selection_df <- bind_rows(
  historical_all_features,
  historical_apriori_feature_selection,
  mutate(mutate(historical_feature_selection, source = "daymet"), feat_sel_type = str_replace(feat_sel_type, "rfa", "stepwise")),
  mutate(historical_feature_importance, source = "daymet", covariates = map(covariates, "importance"),
         covariates = map(covariates, ~subset(rownames_to_column(as.data.frame(.x), "covariate"), importance != 0)))
) %>% select(-contains("time_frame"))


# Give equal-weight importance to covariates except those from the LASSO
feature_selection_df <- feature_selection_df %>%
  mutate(covariates = pmap(select(., covariates, feat_sel_type, model), ~{
    if (grepl("lasso", ..2)) {
      ..1
    } else {
      # Subset for interaction covariates, if model is model3
      if (..3 == "model5") {
        vars <- str_subset(..1$optVariables, "line_name:") %>% str_remove_all(., "line_name:")
      } else {
        vars <- ..1$optVariables
      }
      
      data.frame(covariate = setdiff(vars, "line_name"), stringsAsFactors = FALSE) %>%
        mutate(importance = 1 / nrow(.))
    }
  }))


# Reorganize covariate df
covariates_tomodel <- feature_selection_df %>%
  # Filter the source to use
  filter(source %in% source_use) %>%
  rename(feature_selection = feat_sel_type) %>%
  unnest(covariates) %>%
  filter(covariate != "line_name") %>%
  group_by(source, trait, feature_selection, model) %>%
  nest(.key = "covariates") %>%
  ungroup() %>%
  group_by(trait) %>% 
  nest(.key = "model_covariates") %>%
  ungroup()



# Generate skeleton train/test sets for LOEO
lolo_train_test <- data_to_model %>%
  group_by(trait) %>%
  mutate(site = loc) %>%
  do({crossv_loo_grouped(data = droplevels(group_by(., site)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel, by = "trait")

## Assign cores and split
data_train_test1 <- lolo_train_test %>%
  assign_cores(df = ., n_core = n_core, split = TRUE)


## Parallelize
lolo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      # Get covariates
      covariates_row <- row$model_covariates[[1]]

      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nSite = n_distinct(site), nObs = n())
      
      ##############
      ##############
      
      ## First fit models with identity covariance for environments
      models_run <- lapply(X = model.list, "[", c("model1", "model2_id", "model3_id"))
      
      Emain <- Eint <- diag(nlevels(row$train[[1]]$site))
      dimnames(Emain) <- replicate(2, levels(row$train[[1]]$site), simplify = FALSE)
      dimnames(Eint) <- dimnames(Emain)
      
      # The genomewide prediction function is in the source_functions.R script
      prediction_out_id <- genomewide_prediction2(x = row, model.list = models_run, K = Kgeno, E = Emain, KE = Eint)
      
      ##############
      ##############
      
      ## Relationships for E
      # Subset models
      models_run <- lapply(X = model.list, "[", c("model2_cov", "model3_cov"))
      # Get covariates
      covariates_use <- covariates_row %>%
        spread(model, covariates) %>%
        mutate(out = list(NULL))
      
      # Iterate over rows
      for (r in seq_len(nrow(covariates_use))) {
        
        # What is the feature selection method?
        fsm <- covariates_use$feature_selection[r]
        # Get the timeframe for this trait
        time_frame_use <- subset(time_frame_use_df, 
                                 trait == row$trait & feat_sel_type == ifelse(grepl("lasso", fsm), "lasso_cv_adhoc", "stepwise_cv_adhoc"),
                                 time_frame, drop = TRUE)
        
        # List of covariates
        covariate_list <- covariates_use[r,] %>%
          select(main = model4, interaction = model5) %>%
          as.list() %>%
          map(1)
        
        # Create a matrix of scaled and centered covariates
        covariate_mat <- loc_ec_tomodel_scaled %>%
          filter(time_frame == time_frame_use,
                 location %in% levels(row$train[[1]]$site1)) %>%
          select(-source, -time_frame) %>%
          as.data.frame() %>%
          column_to_rownames("location") %>%
          as.matrix()
        
        ## Environmental relationship matrices
        Emain <- Env_mat(x = covariate_mat, terms = covariate_list$main$covariate,
                         weights = covariate_list$main$importance, method = "weightedJarq")
        if (is.null(covariate_list$interaction)) {
          Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
        } else {
          Eint <- Env_mat(x = covariate_mat, terms = covariate_list$interaction$covariate,
                          weights = covariate_list$interaction$importance, method = "weightedJarq")
        }
        
        # run predictions
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = Kgeno, E = Emain, KE = Eint)
        
        # Add to df
        covariates_use$out[[r]] <- prediction_out_cov$prediction_out
        
      }
      
      prediction_out_cov <- unnest(covariates_use, out)
      
      
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(bind_rows(prediction_out_id$prediction_out, prediction_out_cov), train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in site", as.character(row$site), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()


# External validation -----------------------------------------------


# Data to use
data_to_model <- S2_MET_loc_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         location %in% c(train_test_loc, validation_loc)) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(location), as.factor) %>%
  mutate(loc = location)


# Generate skeleton train/test sets for validation
loc_external_train_val <- data_to_model %>%
  group_by(trait) %>%
  mutate(site = loc) %>%
  do({ tibble(train = list(resample(data = droplevels(.), idx = which(.$site %in% train_test_loc))),
              test = list(resample(data = droplevels(.), idx = which(.$site %in% validation_loc)))) }) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% vp_geno) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel, by = "trait")

## Assign cores and split
loc_external_train_val1 <- loc_external_train_val %>%
  assign_cores(df = ., n_core = 1, split = TRUE)



## Don't parallelize
loc_external_predictions_out <- loc_external_train_val1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of sites and observations used for training
      train_n <- summarize(row$train[[1]], nSite = n_distinct(site), nObs = n())
      
      ##############
      ##############
      
      ## First fit models with identity covariance for environments
      models_run <- lapply(X = model.list, "[", c("model1", "model2_id", "model3_id"))
      
      Emain <- Eint <- diag(nlevels(row$train[[1]]$site))
      dimnames(Emain) <- replicate(2, levels(row$train[[1]]$site), simplify = FALSE)
      dimnames(Eint) <- dimnames(Emain)
      
      # The genomewide prediction function is in the source_functions.R script
      prediction_out_id <- genomewide_prediction2(x = row, model.list = models_run, K = Kgeno, E = Emain, KE = Eint)
      
      ##############
      ##############
      
      ## Relationships for E and GxE
      # Subset models
      models_run <- lapply(X = model.list, "[", c("model2_cov", "model3_cov"))
      # Get covariates
      covariates_use <- covariates_row %>%
        spread(model, covariates) %>%
        mutate(out = list(NULL))
      
      # Iterate over rows
      for (r in seq_len(nrow(covariates_use))) {
        
        # What is the feature selection method?
        fsm <- covariates_use$feature_selection[r]
        # Get the timeframe for this trait
        time_frame_use <- subset(time_frame_use_df, 
                                 trait == row$trait & feat_sel_type == ifelse(grepl("lasso", fsm), "lasso_cv_adhoc", "stepwise_cv_adhoc"),
                                 time_frame, drop = TRUE)
        
        # List of covariates
        covariate_list <- covariates_use[r,] %>%
          select(main = model4, interaction = model5) %>%
          as.list() %>%
          map(1)
        
        # Create a matrix of scaled and centered covariates
        covariate_mat <- loc_ec_tomodel_scaled %>%
          filter(time_frame == time_frame_use,
                 location %in% levels(row$train[[1]]$site1)) %>%
          select(-source, -time_frame) %>%
          as.data.frame() %>%
          column_to_rownames("location") %>%
          as.matrix()
        
        ## Environmental relationship matrices
        Emain <- Env_mat(x = covariate_mat, terms = covariate_list$main$covariate,
                         weights = covariate_list$main$importance, method = "weightedJarq")
        if (is.null(covariate_list$interaction)) {
          Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
        } else {
          Eint <- Env_mat(x = covariate_mat, terms = covariate_list$interaction$covariate,
                          weights = covariate_list$interaction$importance, method = "weightedJarq")
        }
        
        # run predictions
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = Kgeno, E = Emain, KE = Eint)
        
        # Add to df
        covariates_use$out[[r]] <- prediction_out_cov$prediction_out
        
      }
      
      prediction_out_cov <- unnest(covariates_use, out)
      
      
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(bind_rows(prediction_out_id$prediction_out, prediction_out_cov), train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()

## Save the results
save("lolo_predictions_out", "loc_external_predictions_out",
     file = file.path(result_dir, "lolo_predictions_fact_reg.RData"))






