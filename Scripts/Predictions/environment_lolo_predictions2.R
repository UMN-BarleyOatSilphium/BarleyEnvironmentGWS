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

# time_frame to use for the location relationship matrix
time_frame_use <- "time_frame5_2010_2014"

## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

## Load the factorial regression results
load(file.path(result_dir, "feature_selection_results.RData"))

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
  model3_cov1 = model3_cov,
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
  model3_cov1 = add_predictors(model2_cov, ~ vs(line_name:site1, Gu = GxE)),
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
loc_ec_tomodel_scaled <-  historical_ec_tomodel_timeframe_scaled %>%
  subset(., grepl(pattern = time_frame_use, x = names(.)))


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



## Adjust feature selection df
historical_fact_reg_feature_selection <- historical_fact_reg_feature_selection %>%
  gather(selection_type, covariates, apriori, adhoc, adhoc_nosoil) %>% 
  unite(feat_sel_type, feat_sel_type, selection_type, sep = "_") %>% 
  mutate(feat_sel_type = ifelse(str_detect(feat_sel_type, "apriori"), "apriori", feat_sel_type))

historical_feature_selection <- historical_feature_selection %>%
  gather(selection_type, covariates, adhoc, adhoc_nosoil) %>% 
  unite(feat_sel_type, feat_sel_type, selection_type, sep = "_") %>%
  mutate(model = case_when(model == "model2" ~ "model4", model == "model3" ~ "model5"))

# Use the features identified in the environment-specific analysis
concurrent_feature_selection <- concurrent_feature_selection %>%
  gather(selection_type, covariates, adhoc, adhoc_nosoil) %>% 
  unite(feat_sel_type, feat_sel_type, selection_type, sep = "_") %>%
  mutate(model = case_when(model == "model2" ~ "model4", model == "model3" ~ "model5"),
         feat_sel_type = paste0("concurrent_", feat_sel_type))


# Combine feature selection df
feature_selection_df <- bind_rows(historical_fact_reg_feature_selection, historical_feature_selection)
  
# A vector of all covariates
covariates_by_source_df <- tibble(source = map_chr(map(loc_ec_tomodel_scaled, "source"), unique),
                                  covariates = map(loc_ec_tomodel_scaled, ~names(.)[-1:-3]))
# Data.frame with all covariates
all_covariates_df <- crossing(covariates_by_source_df, trait = traits, feature_selection = "all", 
                              model = c("model4", "model5")) %>% 
  mutate(covariates = map2(covariates, model, ~ifelse(.y == "model4", list(.x), list(c(.x, paste0("line_name:", .x))))) %>%
           map(1)) %>%
  filter(source %in% source_use)



# Reorganize covariate df
covariates_tomodel <- feature_selection_df %>%
  bind_rows(., concurrent_feature_selection) %>% # Add covariates selected in the concurrent analysis
  # Filter the source to use
  filter(source %in% source_use) %>%
  rename(feature_selection = feat_sel_type) %>%
  mutate(covariates = map(covariates, "optVariables"),
         covariates = modify_if(covariates, is.null, ~character(0)),
         direction = ifelse(is.na(direction), "forward", direction)) %>%
  # Add all covariates
  bind_rows(., all_covariates_df) %>%
  unnest() %>%
  filter(covariates != "line_name") %>%
  group_by(source, trait, feature_selection, model, direction) %>%
  nest(.key = "covariates") %>%
  ungroup() %>%
  # Do not use the no soil selections or pls
  filter(str_detect(feature_selection, "nosoil|pls", negate = TRUE)) %>%
  # Collapse covariates
  mutate(covariates = map(covariates, "covariates")) %>%
  # nest
  group_by(trait) %>% 
  nest(.key = "model_covariates")



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
      prediction_out_id <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
      
      ##############
      ##############
      
      ## Relationships for E
      # Subset models
      models_run <- lapply(X = model.list, "[", c("model2_cov", "model3_cov", "model3_cov1"))
      # Get covariates
      covariates_use <- covariates_row %>%
        select(-direction) %>%
        spread(model, covariates) %>%
        mutate(out = list(NULL))
      
      # Iterate over rows
      for (r in seq_len(nrow(covariates_use))) {
        
        # List of covariates
        covariate_list <- tibble(term = unique(unlist(covariates_use[r, c("model4", "model5")]))) %>% 
          mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
                 covariate = str_remove(term, "line_name:")) %>%
          split(.$class) %>% 
          map("covariate")
        
        # Create a matrix of scaled and centered covariates
        covariate_mat <- loc_ec_tomodel_scaled %>%
          subset(., grepl(pattern = covariates_use$source[r], x = names(.))) %>% .[[1]] %>%
          filter(location %in% levels(row$train[[1]]$site1)) %>%
          select(., location, unique(unlist(covariate_list))) %>%
          as.data.frame() %>%
          column_to_rownames("location") %>%
          as.matrix()
        
        ## Environmental relationship matrices
        Emain <- Env_mat(x = covariate_mat[,covariate_list$main, drop = FALSE], method = "Jarq")
        if (is.null(covariate_list$interaction)) {
          Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
        } else {
          Eint <- Env_mat(x = covariate_mat[,covariate_list$int, drop = FALSE], method = "Jarq")
        }
        
        # run predictions
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
        
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
      prediction_out_id <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
      
      ##############
      ##############
      
      ## Relationships for E and GxE
      # Subset models
      models_run <- lapply(X = model.list, "[", c("model2_cov", "model3_cov", "model3_cov1"))
      # Get covariates
      covariates_use <- covariates_row %>%
        select(-direction) %>%
        spread(model, covariates) %>%
        mutate(out = list(NULL))
      
      # Iterate over rows
      for (r in seq_len(nrow(covariates_use))) {
        
        # List of covariates
        covariate_list <- tibble(term = unique(unlist(covariates_use[r, c("model4", "model5")]))) %>% 
          mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
                 covariate = str_remove(term, "line_name:")) %>%
          split(.$class) %>% 
          map("covariate")
        
        # Create a matrix of scaled and centered covariates
        covariate_mat <- loc_ec_tomodel_scaled %>%
          subset(., grepl(pattern = covariates_use$source[r], x = names(.))) %>% .[[1]] %>%
          filter(location %in% levels(row$train[[1]]$site1)) %>%
          select(., location, unique(unlist(covariate_list))) %>%
          as.data.frame() %>%
          column_to_rownames("location") %>%
          as.matrix()
        
        ## Environmental relationship matrices
        Emain <- Env_mat(x = covariate_mat[,covariate_list$main, drop = FALSE], method = "Jarq")
        if (is.null(covariate_list$interaction)) {
          Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
        } else {
          Eint <- Env_mat(x = covariate_mat[,covariate_list$int, drop = FALSE], method = "Jarq")
        }
        
        # run predictions
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
        
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






