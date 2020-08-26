## S2MET Prediction Models
## 
## 75-25 environmental or location cross-validation
## 
## Author: Jeff Neyhart
## 


# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Run the source script for MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))


# Other packages
library(modelr)
library(broom)
library(parallel)


## Number of cores
# n_core <- detectCores()
n_core <- 12


# Proportion of data to use as test
pTest <- 0.25
pTrain <- 1 - pTest
# Number of train/test partitions
nCV <- 25

# Source of covariates
source_use <- "daymet"

## Data.frame of timeframes to use per trait
time_frame_use_df <- historical_feature_selection %>%
  distinct(trait, time_frame)

# If re-running predictions, should all be re-run?
rerun_all <- TRUE



## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

## Load the factorial regression results
load(file.path(result_dir, "feature_selection_results.RData"))




## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value,
  model1 = ~ 1,
  model2_cov = model1,
  model2_id = model1,
  model3_id = model1,
  model3_cov = model1,
)

## Models for de novo fitting 
## Create a list of model formulas
model_rand_forms <- formulas(
  .response = ~ value,
  model1 = ~ vs(line_name, Gu = K),
  model2_cov = add_predictors(model1, ~ vs(site1, Gu = E)),
  model2_id = add_predictors(model1, ~ vs(site1, Gu = I)),
  model3_id = add_predictors(model2_id, ~ vs(line_name:site1, Gu = GI)),
  model3_cov = add_predictors(model2_cov, ~ vs(line_name:site1, Gu = GE))
) %>% map(~ formula(delete.response(terms(.)))) # Remove response

# Combine into list
model.list <- list(fixed = model_fixed_forms, random = model_rand_forms)





# Environments -----------------------------------------------

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% train_test_env) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment, loc = location)

# Rename the covariate list
concurrent_all_features <- concurrent_all_features %>%
  mutate(covariates = map(covariates, ~setNames(.x, "optVariables")))

# Remove soil variable from all variables
concurrent_all_features <- concurrent_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(concurrent_all_features, .)

# Combine feature selection df
feature_selection_df <- bind_rows(concurrent_fact_reg_feature_selection, concurrent_feature_selection, 
                                  concurrent_all_features)

# Reorganize covariate df
covariates_tomodel <- feature_selection_df %>%
  select(-direction) %>%
  # Filter the source to use
  filter(source %in% source_use,
         # Do not use stepwiseAIC covariates
         str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  rename(feature_selection = feat_sel_type) %>%
  mutate(covariates = map(covariates, "optVariables")) %>%
  unnest(covariates) %>%
  filter(covariates != "line_name") %>%
  group_by(source, trait, feature_selection, model) %>%
  nest(.key = "covariates") %>%
  ungroup() %>%
  # Collapse covariates
  mutate(covariates = map(covariates, "covariates")) %>%
  # nest
  group_by(trait) %>% 
  nest(.key = "model_covariates") %>%
  ungroup()


# Generate skeleton train/test sets for LOEO
cv_train_test <- data_to_model %>%
  group_by(trait) %>%
  mutate(site = env) %>%
  do({crossv_mc_grouped(data = droplevels(group_by(., site)), n = nCV, test = pTest)}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)


## Pull the scenarios to run, assign cores and split
data_train_test1 <- cv_train_test %>%
  # inner_join(scenarios_torun, .) %>%
  assign_cores(df = ., n_core = n_core, split = TRUE)


## Parallelize
cv_env_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nSite = n_distinct(environment), nObs = n())
      test_n <- summarize(row$test[[1]], nSite = n_distinct(environment), nObs = n())
      
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
        
        # Covariate data source
        src <- covariates_use$source[r]
        
        # List of covariates
        covariate_list <- tibble(term = unique(unlist(covariates_use[r, c("model2", "model3")]))) %>% 
          mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
                 covariate = str_remove(term, "line_name:")) %>%
          split(.$class) %>% 
          map("covariate")
        
        # Create a matrix of scaled and centered covariates
        covariate_mat <- ec_tomodel_scaled[[src]] %>%
          filter(environment %in% levels(row$train[[1]]$site1)) %>%
          select(., environment, unique(unlist(covariate_list))) %>%
          as.data.frame() %>%
          column_to_rownames("environment") %>%
          as.matrix()
        
        ## Environmental relationship matrices
        Emain <- Env_mat(x = covariate_mat[,covariate_list$main, drop = FALSE], method = "Jarq")
        if (is.null(covariate_list$interaction)) {
          Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
        } else {
          Eint <- Env_mat(x = covariate_mat[,covariate_list$int, drop = FALSE], method = "Jarq")
        }
        
        # run predictions
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = Kgeno, E = Emain, KE = Eint)
        
        # Add to df
        covariates_use$out[[r]] <- prediction_out_cov$prediction_out
        
      }
      
      prediction_out_features <- unnest(covariates_use, out)
      
      
      ##############
      ##############
      
      
      ## Combine and return the predictions
      out[[i]] <- mutate(bind_rows(prediction_out_id$prediction_out, prediction_out_features),
                         train_n = list(train_n)) %>% unnest(train_n)
      # out[[i]] <- mutate(bind_rows(prediction_out_features, prediction_out_model3_cov_all), 
      #                    train_n = list(train_n)) %>% unnest(train_n)
      
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()




# Locations -----------------------------------------------


# Data to use - location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

## Pull out location covariates to use
loc_ec_tomodel_scaled <-  historical_ec_tomodel_timeframe_scaled %>%
  bind_rows() %>%
  filter(source == source_use, time_frame %in% unique(time_frame_use_df$time_frame))


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

# Rename the covariate list
historical_all_features <- historical_all_features %>%
  mutate(covariates = map(covariates, ~setNames(.x, "optVariables")))

# Remove soil variable from all variables
historical_all_features <- historical_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(historical_all_features, .)

# Combine feature selection df
feature_selection_df <- bind_rows(historical_fact_reg_feature_selection, mutate(historical_feature_selection, source = "daymet"),
                                  historical_all_features) %>%
  select(-contains("time_frame"))

# Reorganize covariate df
covariates_tomodel <- feature_selection_df %>%
  # bind_rows(., concurrent_feature_selection) %>% # Add covariates selected in the concurrent analysis
  # Filter the source to use
  filter(source %in% source_use) %>%
  # Do not use the AIC stepwise covariates
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  rename(feature_selection = feat_sel_type) %>%
  select(-direction) %>%
  mutate(covariates = map(covariates, "optVariables"),
         covariates = modify_if(covariates, is.null, ~character(0))) %>%
  unnest() %>%
  filter(covariates != "line_name") %>%
  group_by(source, trait, feature_selection, model) %>%
  nest(.key = "covariates") %>%
  ungroup() %>%
  # Collapse covariates
  mutate(covariates = map(covariates, "covariates")) %>%
  # nest
  group_by(trait) %>% 
  nest(.key = "model_covariates")


# Generate skeleton train/test sets for LOEO
cv_train_test <- data_to_model %>%
  group_by(trait) %>%
  mutate(site = loc) %>%
  do({crossv_mc_grouped(data = droplevels(group_by(., site)), n = nCV, test = pTest)}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)


## Pull the scenarios to run, assign cores and split
data_train_test1 <- cv_train_test %>%
  # inner_join(scenarios_torun, .) %>%
  assign_cores(df = ., n_core = n_core, split = TRUE)


## Parallelize
cv_loc_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      # Get the timeframe for this trait
      time_frame_use <- subset(time_frame_use_df, trait == row$trait, time_frame, drop = TRUE)
      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nSite = n_distinct(site), nObs = n())
      test_n <- summarize(row$test[[1]], nSite = n_distinct(site), nObs = n())
      
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
        
        # Covariate data source
        src <- covariates_use$source[r]
        
        # List of covariates
        covariate_list <- tibble(term = unique(unlist(covariates_use[r, c("model4", "model5")]))) %>% 
          mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
                 covariate = str_remove(term, "line_name:")) %>%
          split(.$class) %>% 
          map("covariate")
        
        # Create a matrix of scaled and centered covariates
        covariate_mat <- loc_ec_tomodel_scaled %>%
          filter(time_frame == time_frame_use,
                 location %in% levels(row$train[[1]]$site1)) %>%
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
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = Kgeno, E = Emain, KE = Eint)
        
        # Add to df
        covariates_use$out[[r]] <- prediction_out_cov$prediction_out
        
      }
      
      prediction_out_features <- unnest(covariates_use, out)
      
      
      ##############
      ##############
      
      
      ## Combine and return the predictions
      out[[i]] <- mutate(bind_rows(prediction_out_id$prediction_out, prediction_out_features),
                         train_n = list(train_n)) %>% unnest(train_n)
      # out[[i]] <- mutate(bind_rows(prediction_out_features, prediction_out_model3_cov_all), 
      #                    train_n = list(train_n)) %>% unnest(train_n)
      
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()








## Save the results
save("cv_env_predictions_out", "cv_loc_predictions_out", 
     file = file.path(result_dir, "cv7525_predictions_fact_reg.RData"))
