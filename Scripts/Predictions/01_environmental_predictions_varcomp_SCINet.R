## S2MET Prediction Models
## 
## All predictions using previously estimated variance components
## 
## Leave-one-environment-out prediction
## Leave-one-location-out prediction
## External environment/location prediction
## 
## Author: Jeff Neyhart
## 

if (Sys.info()["sysname"] == "Windows") {

  # Run on a local machine
  repo_dir <- getwd()
  source(file.path(repo_dir, "source.R"))
  
} else {
  
  # Run on SCINet
  repo_dir <- "/project/gifvl_vaccinium/barley_work/BarleyEnvironmentGWS"
  source(file.path(repo_dir, "source_SCINet.R"))
  
}



# Other packages
library(modelr)
library(broom)
library(parallel)


# Source of covariates
source_use <- "daymet"

# Number of cores
n_core <- 16



# Load data ---------------------------------------------------------------


## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

## Load the factorial regression results
load(file.path(result_dir, "concurrent_feature_selection_results.RData"))

## Load the interval covariate results
load(file.path(result_dir, "concurrent_intervalCovariates_feature_selection_results.RData"))

## Load the variance component estimates
load(file.path(result_dir, "stage_two_phenotypic_analysis.RData"))


# REmove the marker matrix
rm(s2_imputed_mat_use)



# Environment predictions -------------------------------------------------


#######
# Set model formula
#######

## Create a list of model formulas
model_fixed_forms <- formulas(
  .response = ~ value,
  model1 = ~ 1,
  model2_cov = model1,
  model2_id = model1,
  model3_id = model1,
  model3_cov = model1,
)

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


#####
# LOEO - Set up data sets
#####

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

# Subset stage-two variance components
env_varcomp_df <- stage_two_varcomp_env %>%
  unnest(var_comp) %>%
  filter(population == "tp", scaled == "scaled") %>%
  unnest(varcomp) %>%
  select(trait, model, term = grp, variance = vcov) %>%
  mutate(term = str_replace_all(term, "environment", "site1"), 
         term = ifelse(term == "Residual", "units", term))



# Rename the covariate list
concurrent_all_features <- concurrent_all_features %>%
  mutate(covariates = map(covariates, ~setNames(.x, "optVariables")))

# Remove soil variable from all variables
concurrent_all_features <- concurrent_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(concurrent_all_features, .)

# Combine feature selection df
feature_selection_df <- bind_rows(
  concurrent_all_features,
  concurrent_apriori_feature_selection,
  mutate(concurrent_stepwise_feature_selection, feat_sel_type = str_replace(feat_sel_type, "rfa", "stepwise")),
) %>%
  # Create a data.frame of covariates
  mutate(covariates = map(covariates, "optVariables"),
         covariates = map(covariates, ~tibble(covariate = .x)))
  
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

# Assign the list of matrices of environmental covariates
ec_mat_tomodel <- ec_tomodel_scaled



#####
# LOEO - run predictions
#####

# Generate skeleton train/test sets for LOEO
loeo_train_test <- data_to_model %>%
  mutate(site = env) %>%
  group_by(trait) %>%
  do({crossv_loo_grouped(data = droplevels(group_by(., site)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)


# Assign cores, split, and parallelize

## Pull the scenarios to run, assign cores and split
loeo_predictions_out <- loeo_train_test %>%
  assign_cores(df = ., n_core = n_core, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    # Output list
    predictions_out <- vector("list", length = nrow(core_df))
    
    # Iterate over rows in core_df
    for (i in seq_along(predictions_out)) {
      
      row <- core_df[i,]
      data <- row$train[[1]] %>%
        # Modify the weights
        mutate(weights = std_error^2, site1 = site) %>%
        as.data.frame()
      
      # Test data
      test <- row$test[[1]]
        
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(data, nSite = n_distinct(environment), nObs = n())
      
      ##############
      ##############
      
      # Create a df of model and covariate combinations
      models_df <- tibble(source = "daymet", feature_selection = "none", model = paste0("model", 1:3),
                          covariates = list(tibble(covariate = character()))) %>%
        bind_rows(., covariates_row) %>%
        mutate(submodel = ifelse(feature_selection == "none", paste0(model, "_id"), paste0(model, "_cov")),
               submodel = ifelse(submodel == "model1_id", "model1", submodel),
               out = list(NULL))
      
      # Iterate over rows of models
      for (r in seq_len(nrow(models_df))) {
        
        # Get the submodel, declare relationship matrices
        submodel <- models_df$submodel[[r]]
        model_r <- models_df$model[[r]]
        
        # Extract the variance components
        varcomp_r <- env_varcomp_df %>% 
          filter(trait == row$trait, model == model_r) %>%
          select(term, variance) %>%
          spread(term, variance) %>%
          mutate_all(~list(varcomp = matrix(.x))) %>%
          as.list()
        
        
        if (submodel == "model1" | endsWith(submodel, "id")) {
          Emain <- Eint <- diag(nlevels(data$site))
          dimnames(Emain) <- replicate(2, levels(data$site), simplify = FALSE)
          dimnames(Eint) <- dimnames(Emain)
          
        } else {
          # Covariate data source
          src <- models_df$source[r]
          
          # Grab the covariates (main effect and interaction covariates are listed together)
          covariate_list <- models_df$covariates[[r]] %>%
            mutate(group = ifelse(str_detect(covariate, ":"), "interaction", "main")) %>%
            split(.$group) %>%
            map(~mutate(.x, covariate = str_remove(covariate, "line_name:")))
          
          # Create main matrix of scaled and centered covariate values
          covariate_mat <- ec_mat_tomodel[[src]] %>%
            select(-source) %>%
            filter(environment %in% levels(data$site1)) %>%
            as.data.frame() %>%
            column_to_rownames("environment") %>%
            as.matrix()
          
          ## Environmental relationship matrices
          Emain <- Env_mat(x = covariate_mat, terms = covariate_list$main$covariate, method = "Jarq")
          if (is.null(covariate_list$interaction)) {
            Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
          } else {
            Eint <- Env_mat(x = covariate_mat, terms = covariate_list$interaction$covariate, method = "Jarq")
          }
          
        }
        
        # GxE matrix
        GE <- kronecker(K, Eint, make.dimnames = TRUE)
        
        # Add genetic covariance matrix for line_names
        random <- varcomp_r %>%
          subset(., names(.) != "units") %>%
          modify_at(., "line_name", ~c(.x, list(K = K)))
        
        # Add environmental information if called for
        if (model_r == "model2") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain)))
          
        } else if (model_r == "model3") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain))) %>%
            modify_at("line_name:site1", ~c(.x, list(K = GE)))
          
        }
        
        # Fixed and random formula
        model_run <- map(model.list, submodel)
        fixed <- model_run$fixed
        
        # Modify the training data
        train <- data %>%
          select(line_name, site1, value, weights) %>%
          mutate(line_name = factor(line_name, levels = row.names(K)),
                 site1 = factor(site1, levels = row.names(Emain)))
    
        # Fit the model
        fit <- fit_lmm(fixed = fixed, random_cov = random, r_cov = varcomp_r["units"], data = train,
                       weights = "weights", verbose = TRUE)
        
        # Generate predictions
        predictions <- predict_lmm(object = fit, terms = names(random)) %>%
          # Join with the test set
          inner_join(test, ., by = intersect(names(random), names(test))) %>%
          select(trait, line_name, site, value, std_error, pred, pev, reliability, predicted_value)
        
        # Store the predictions in the list
        models_df$out[[r]] <- predictions
        
      }
      
      # Store the models_df data frame in the predictions output list
      predictions_out[[i]] <- models_df
      
    } # End of core_df loop
    
    # Add the predictions out list to the core_df data.frame
    core_df %>%
      mutate(predictions = predictions_out) %>%
      select(-train, -test, -model_covariates, -core)
    
  }) %>% bind_rows() # End of parallelization

# Save as tibble
loeo_varcomp_predictions_out <- as_tibble(loeo_predictions_out)




#####
# External validation - Set up data sets
#####

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% c(train_test_env, validation_env)) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment, loc = location)

#####
# External validation - run predictions
#####


# Generate skeleton train/test sets for validation
env_external_train_val <- data_to_model %>%
  mutate(site = env) %>%
  group_by(trait) %>%
  do({ tibble(train = list(resample(data = droplevels(.), idx = which(.$site %in% train_test_env))),
              test = list(resample(data = droplevels(.), idx = which(.$site %in% validation_env)))) }) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% vp_geno) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)


## Pull the scenarios to run, assign cores and split
env_external_predictions_out <- env_external_train_val %>%
  assign_cores(df = ., n_core = 5, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    # Output list
    predictions_out <- vector("list", length = nrow(core_df))
    
    # Iterate over rows in core_df
    for (i in seq_along(predictions_out)) {
      
      row <- core_df[i,]
      data <- row$train[[1]] %>%
        # Modify the weights
        mutate(weights = std_error^2, site1 = site) %>%
        as.data.frame()
      
      # Test data
      test <- row$test[[1]]
      
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(data, nSite = n_distinct(environment), nObs = n())
      
      ##############
      ##############
      
      # Create a df of model and covariate combinations
      models_df <- tibble(source = "daymet", feature_selection = "none", model = paste0("model", 1:3),
                          covariates = list(tibble(covariate = character()))) %>%
        bind_rows(., covariates_row) %>%
        mutate(submodel = ifelse(feature_selection == "none", paste0(model, "_id"), paste0(model, "_cov")),
               submodel = ifelse(submodel == "model1_id", "model1", submodel),
               out = list(NULL))
      
      # Iterate over rows of models
      for (r in seq_len(nrow(models_df))) {
        
        # Get the submodel, declare relationship matrices
        submodel <- models_df$submodel[[r]]
        model_r <- models_df$model[[r]]
        
        # Extract the variance components
        varcomp_r <- env_varcomp_df %>% 
          filter(trait == row$trait, model == model_r) %>%
          select(term, variance) %>%
          spread(term, variance) %>%
          mutate_all(~list(varcomp = matrix(.x))) %>%
          as.list()
        
        
        if (submodel == "model1" | endsWith(submodel, "id")) {
          Emain <- Eint <- diag(nlevels(data$site))
          dimnames(Emain) <- replicate(2, levels(data$site), simplify = FALSE)
          dimnames(Eint) <- dimnames(Emain)
          
        } else {
          # Covariate data source
          src <- models_df$source[r]
          
          # Grab the covariates (main effect and interaction covariates are listed together)
          covariate_list <- models_df$covariates[[r]] %>%
            mutate(group = ifelse(str_detect(covariate, ":"), "interaction", "main")) %>%
            split(.$group) %>%
            map(~mutate(.x, covariate = str_remove(covariate, "line_name:")))
          
          # Create main matrix of scaled and centered covariate values
          covariate_mat <- ec_mat_tomodel[[src]] %>%
            select(-source) %>%
            filter(environment %in% levels(data$site1)) %>%
            as.data.frame() %>%
            column_to_rownames("environment") %>%
            as.matrix()
          
          ## Environmental relationship matrices
          Emain <- Env_mat(x = covariate_mat, terms = covariate_list$main$covariate, method = "Jarq")
          if (is.null(covariate_list$interaction)) {
            Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
          } else {
            Eint <- Env_mat(x = covariate_mat, terms = covariate_list$interaction$covariate, method = "Jarq")
          }
          
        }
        
        # GxE matrix
        GE <- kronecker(K, Eint, make.dimnames = TRUE)
        
        # Add genetic covariance matrix for line_names
        random <- varcomp_r %>%
          subset(., names(.) != "units") %>%
          modify_at(., "line_name", ~c(.x, list(K = K)))
        
        # Add environmental information if called for
        if (model_r == "model2") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain)))
          
        } else if (model_r == "model3") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain))) %>%
            modify_at("line_name:site1", ~c(.x, list(K = GE)))
          
        }
        
        # Fixed and random formula
        model_run <- map(model.list, submodel)
        fixed <- model_run$fixed
        
        # Modify the training data
        train <- data %>%
          select(line_name, site1, value, weights) %>%
          mutate(line_name = factor(line_name, levels = row.names(K)),
                 site1 = factor(site1, levels = row.names(Emain)))
        
        # Fit the model
        fit <- fit_lmm(fixed = fixed, random_cov = random, r_cov = varcomp_r["units"], data = train,
                       weights = "weights", verbose = TRUE)
        
        # Generate predictions
        predictions <- predict_lmm(object = fit, terms = names(random)) %>%
          # Join with the test set
          inner_join(test, ., by = intersect(names(random), names(test))) %>%
          select(trait, line_name, site, value, std_error, pred, pev, reliability, predicted_value)
        
        # Store the predictions in the list
        models_df$out[[r]] <- predictions
        
      }
      
      # Store the models_df data frame in the predictions output list
      predictions_out[[i]] <- models_df
      
    } # End of core_df loop
    
    # Add the predictions out list to the core_df data.frame
    core_df %>%
      mutate(predictions = predictions_out) %>%
      select(-train, -test, -model_covariates, -core)
    
  }) %>% bind_rows() # End of parallelization


# Save as tibble
env_external_varcomp_predictions_out <- as_tibble(env_external_predictions_out)




# Environment predictions - interval covariates ---------------------------


#####
# LOEO - Set up data sets
#####

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
concurrent_interval_all_features <- concurrent_interval_all_features %>%
  mutate(covariates = map(covariates, ~setNames(.x, "optVariables")))

# Remove soil variable from all variables
concurrent_interval_all_features <- concurrent_interval_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(concurrent_interval_all_features, .)

# Combine feature selection df
feature_selection_df <- bind_rows(
  concurrent_interval_all_features,
  mutate(concurrent_stepwise_interval_feature_selection, feat_sel_type = str_replace(feat_sel_type, "rfa", "stepwise")),
) %>%
  # Create a data.frame of covariates
  mutate(covariates = map(covariates, "optVariables"),
         covariates = map(covariates, ~tibble(covariate = .x)))

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


# Assign the list of matrices of environmental covariates
ec_mat_tomodel <- ec_tomodel_interval_scaled



#####
# LOEO - run predictions
#####

# Generate skeleton train/test sets for LOEO
loeo_train_test <- data_to_model %>%
  mutate(site = env) %>%
  group_by(trait) %>%
  do({crossv_loo_grouped(data = droplevels(group_by(., site)))}) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% c(tp_geno, vp_geno)) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)


# Assign cores, split, and parallelize

## Pull the scenarios to run, assign cores and split
loeo_predictions_out <- loeo_train_test %>%
  assign_cores(df = ., n_core = n_core, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    # Output list
    predictions_out <- vector("list", length = nrow(core_df))
    
    # Iterate over rows in core_df
    for (i in seq_along(predictions_out)) {
      
      row <- core_df[i,]
      data <- row$train[[1]] %>%
        # Modify the weights
        mutate(weights = std_error^2, site1 = site) %>%
        as.data.frame()
      
      # Test data
      test <- row$test[[1]]
      
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(data, nSite = n_distinct(environment), nObs = n())
      
      ##############
      ##############
      
      # Create a df of model and covariate combinations
      models_df <- tibble(source = "daymet", feature_selection = "none", model = paste0("model", 1:3),
                          covariates = list(tibble(covariate = character()))) %>%
        bind_rows(., covariates_row) %>%
        mutate(submodel = ifelse(feature_selection == "none", paste0(model, "_id"), paste0(model, "_cov")),
               submodel = ifelse(submodel == "model1_id", "model1", submodel),
               out = list(NULL))
      
      # Iterate over rows of models
      for (r in seq_len(nrow(models_df))) {
        
        # Get the submodel, declare relationship matrices
        submodel <- models_df$submodel[[r]]
        model_r <- models_df$model[[r]]
        
        # Extract the variance components
        varcomp_r <- env_varcomp_df %>% 
          filter(trait == row$trait, model == model_r) %>%
          select(term, variance) %>%
          spread(term, variance) %>%
          mutate_all(~list(varcomp = matrix(.x))) %>%
          as.list()
        
        
        if (submodel == "model1" | endsWith(submodel, "id")) {
          Emain <- Eint <- diag(nlevels(data$site))
          dimnames(Emain) <- replicate(2, levels(data$site), simplify = FALSE)
          dimnames(Eint) <- dimnames(Emain)
          
        } else {
          # Covariate data source
          src <- models_df$source[r]
          
          # Grab the covariates (main effect and interaction covariates are listed together)
          covariate_list <- models_df$covariates[[r]] %>%
            mutate(group = ifelse(str_detect(covariate, ":"), "interaction", "main")) %>%
            split(.$group) %>%
            map(~mutate(.x, covariate = str_remove(covariate, "line_name:")))
          
          # Create main matrix of scaled and centered covariate values
          covariate_mat <- ec_mat_tomodel[[src]] %>%
            select(-source) %>%
            filter(environment %in% levels(data$site1)) %>%
            as.data.frame() %>%
            column_to_rownames("environment") %>%
            as.matrix()
          
          ## Environmental relationship matrices
          Emain <- Env_mat(x = covariate_mat, terms = covariate_list$main$covariate, method = "Jarq")
          if (is.null(covariate_list$interaction)) {
            Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
          } else {
            Eint <- Env_mat(x = covariate_mat, terms = covariate_list$interaction$covariate, method = "Jarq")
          }
          
        }
        
        # GxE matrix
        GE <- kronecker(K, Eint, make.dimnames = TRUE)
        
        # Add genetic covariance matrix for line_names
        random <- varcomp_r %>%
          subset(., names(.) != "units") %>%
          modify_at(., "line_name", ~c(.x, list(K = K)))
        
        # Add environmental information if called for
        if (model_r == "model2") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain)))
          
        } else if (model_r == "model3") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain))) %>%
            modify_at("line_name:site1", ~c(.x, list(K = GE)))
          
        }
        
        # Fixed and random formula
        model_run <- map(model.list, submodel)
        fixed <- model_run$fixed
        
        # Modify the training data
        train <- data %>%
          select(line_name, site1, value, weights) %>%
          mutate(line_name = factor(line_name, levels = row.names(K)),
                 site1 = factor(site1, levels = row.names(Emain)))
        
        # Fit the model
        fit <- fit_lmm(fixed = fixed, random_cov = random, r_cov = varcomp_r["units"], data = train,
                       weights = "weights", verbose = TRUE)
        
        # Generate predictions
        predictions <- predict_lmm(object = fit, terms = names(random)) %>%
          # Join with the test set
          inner_join(test, ., by = intersect(names(random), names(test))) %>%
          select(trait, line_name, site, value, std_error, pred, pev, reliability, predicted_value)
        
        # Store the predictions in the list
        models_df$out[[r]] <- predictions
        
      }
      
      # Store the models_df data frame in the predictions output list
      predictions_out[[i]] <- models_df
      
    } # End of core_df loop
    
    # Add the predictions out list to the core_df data.frame
    core_df %>%
      mutate(predictions = predictions_out) %>%
      select(-train, -test, -model_covariates, -core)
    
  }) %>% bind_rows() # End of parallelization

# Save as tibble
loeo_varcomp_interval_predictions_out <- as_tibble(loeo_predictions_out)




#####
# External validation - Set up data sets
#####

# Data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% c(train_test_env, validation_env)) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment, loc = location)

#####
# External validation - run predictions
#####


# Generate skeleton train/test sets for validation
env_external_train_val <- data_to_model %>%
  mutate(site = env) %>%
  group_by(trait) %>%
  do({ tibble(train = list(resample(data = droplevels(.), idx = which(.$site %in% train_test_env))),
              test = list(resample(data = droplevels(.), idx = which(.$site %in% validation_env)))) }) %>%
  ungroup() %>%
  mutate(train = map(train, ~filter(as.data.frame(.), line_name %in% tp_geno) %>% 
                       mutate_at(vars(site), fct_contr_sum) %>%
                       mutate(site1 = site)),
         test = map(test, ~filter(as.data.frame(.), line_name %in% vp_geno) %>% 
                      mutate_at(vars(site), fct_contr_sum) %>%
                      mutate(site1 = site)) ) %>%
  # Combine with the different covariate sets
  left_join(., covariates_tomodel)


## Pull the scenarios to run, assign cores and split
env_external_predictions_out <- env_external_train_val %>%
  assign_cores(df = ., n_core = 5, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    # Output list
    predictions_out <- vector("list", length = nrow(core_df))
    
    # Iterate over rows in core_df
    for (i in seq_along(predictions_out)) {
      
      row <- core_df[i,]
      data <- row$train[[1]] %>%
        # Modify the weights
        mutate(weights = std_error^2, site1 = site) %>%
        as.data.frame()
      
      # Test data
      test <- row$test[[1]]
      
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(data, nSite = n_distinct(environment), nObs = n())
      
      ##############
      ##############
      
      # Create a df of model and covariate combinations
      models_df <- tibble(source = "daymet", feature_selection = "none", model = paste0("model", 1:3),
                          covariates = list(tibble(covariate = character()))) %>%
        bind_rows(., covariates_row) %>%
        mutate(submodel = ifelse(feature_selection == "none", paste0(model, "_id"), paste0(model, "_cov")),
               submodel = ifelse(submodel == "model1_id", "model1", submodel),
               out = list(NULL))
      
      # Iterate over rows of models
      for (r in seq_len(nrow(models_df))) {
        
        # Get the submodel, declare relationship matrices
        submodel <- models_df$submodel[[r]]
        model_r <- models_df$model[[r]]
        
        # Extract the variance components
        varcomp_r <- env_varcomp_df %>% 
          filter(trait == row$trait, model == model_r) %>%
          select(term, variance) %>%
          spread(term, variance) %>%
          mutate_all(~list(varcomp = matrix(.x))) %>%
          as.list()
        
        
        if (submodel == "model1" | endsWith(submodel, "id")) {
          Emain <- Eint <- diag(nlevels(data$site))
          dimnames(Emain) <- replicate(2, levels(data$site), simplify = FALSE)
          dimnames(Eint) <- dimnames(Emain)
          
        } else {
          # Covariate data source
          src <- models_df$source[r]
          
          # Grab the covariates (main effect and interaction covariates are listed together)
          covariate_list <- models_df$covariates[[r]] %>%
            mutate(group = ifelse(str_detect(covariate, ":"), "interaction", "main")) %>%
            split(.$group) %>%
            map(~mutate(.x, covariate = str_remove(covariate, "line_name:")))
          
          # Create main matrix of scaled and centered covariate values
          covariate_mat <- ec_mat_tomodel[[src]] %>%
            select(-source) %>%
            filter(environment %in% levels(data$site1)) %>%
            as.data.frame() %>%
            column_to_rownames("environment") %>%
            as.matrix()
          
          ## Environmental relationship matrices
          Emain <- Env_mat(x = covariate_mat, terms = covariate_list$main$covariate, method = "Jarq")
          if (is.null(covariate_list$interaction)) {
            Eint <- diag(ncol(Emain)); dimnames(Eint) <- dimnames(Emain)
          } else {
            Eint <- Env_mat(x = covariate_mat, terms = covariate_list$interaction$covariate, method = "Jarq")
          }
          
        }
        
        # GxE matrix
        GE <- kronecker(K, Eint, make.dimnames = TRUE)
        
        # Add genetic covariance matrix for line_names
        random <- varcomp_r %>%
          subset(., names(.) != "units") %>%
          modify_at(., "line_name", ~c(.x, list(K = K)))
        
        # Add environmental information if called for
        if (model_r == "model2") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain)))
          
        } else if (model_r == "model3") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain))) %>%
            modify_at("line_name:site1", ~c(.x, list(K = GE)))
          
        }
        
        # Fixed and random formula
        model_run <- map(model.list, submodel)
        fixed <- model_run$fixed
        
        # Modify the training data
        train <- data %>%
          select(line_name, site1, value, weights) %>%
          mutate(line_name = factor(line_name, levels = row.names(K)),
                 site1 = factor(site1, levels = row.names(Emain)))
        
        # Fit the model
        fit <- fit_lmm(fixed = fixed, random_cov = random, r_cov = varcomp_r["units"], data = train,
                       weights = "weights", verbose = TRUE)
        
        # Generate predictions
        predictions <- predict_lmm(object = fit, terms = names(random)) %>%
          # Join with the test set
          inner_join(test, ., by = intersect(names(random), names(test))) %>%
          select(trait, line_name, site, value, std_error, pred, pev, reliability, predicted_value)
        
        # Store the predictions in the list
        models_df$out[[r]] <- predictions
        
      }
      
      # Store the models_df data frame in the predictions output list
      predictions_out[[i]] <- models_df
      
    } # End of core_df loop
    
    # Add the predictions out list to the core_df data.frame
    core_df %>%
      mutate(predictions = predictions_out) %>%
      select(-train, -test, -model_covariates, -core)
    
  }) %>% bind_rows() # End of parallelization


# Save as tibble
env_external_interval_varcomp_predictions_out <- as_tibble(env_external_predictions_out)



save("loeo_varcomp_predictions_out", "env_external_varcomp_predictions_out",
     "loeo_varcomp_interval_predictions_out", "env_external_interval_varcomp_predictions_out",
     file = file.path(result_dir, "environment_varcomp_predictions.RData"))

















