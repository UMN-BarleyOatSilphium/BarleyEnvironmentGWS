## S2MET Prediction Models
## 
## All predictions using previously estimated variance components
## 
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
load(file.path(result_dir, "historical_feature_selection_results.RData"))

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

## Data.frame of timeframes to use per trait
time_frame_use_df <- historical_feature_selection %>%
  distinct(trait, time_frame) %>%
  arrange(trait)

# Data to use - location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value), .groups = "drop")



#####
# LOLO - Set up data sets
#####

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

# Subset stage-two variance components
loc_varcomp_df <- stage_two_varcomp_loc %>%
  unnest(var_comp) %>%
  filter(population == "tp", scaled == "scaled") %>%
  unnest(varcomp) %>%
  select(trait, model, term = grp, variance = vcov) %>%
  mutate(term = str_replace_all(term, "location", "site1"), 
         term = ifelse(term == "Residual", "units", term))



# Rename the covariate list
historical_all_features <- historical_all_features %>%
  mutate(covariates = map(covariates, ~setNames(.x, "optVariables")))

# Remove soil variable from all variables
historical_all_features <- historical_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(historical_all_features, .)

# Combine feature selection df
feature_selection_df <- bind_rows(
  historical_all_features,
  historical_apriori_feature_selection,
  historical_feature_selection,
) %>%
  select(-contains("time_frame")) %>%
  # Create a data.frame of covariates
  mutate(covariates = map(covariates, "optVariables"),
         covariates = map(covariates, ~tibble(covariate = .x))) %>%
  # Add time_frame information
  left_join(., time_frame_use_df)
  
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
ec_mat_tomodel <- c(historical_ec_tomodel_timeframe_scaled, historical_ec_tomodel_window_scaled)
ec_mat_tomodel <- bind_rows(ec_mat_tomodel[paste0(source_use, ".", unique(time_frame_use_df$time_frame))])



#####
# LOLO and External validation  - run predictions
#####

# Generate skeleton train/test sets for LOLO
lolo_train_test <- data_to_model %>%
  mutate(site = loc) %>%
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
  left_join(., covariates_tomodel) %>%
  mutate(type = "lolo")

# Generate skeleton train/test sets for external predictions
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
  left_join(., covariates_tomodel, by = "trait") %>%
  mutate(site = "holdout", type = "loc_external")


# Combine train/test data frames
# Assign cores, split, and parallelize

## Pull the scenarios to run, assign cores and split
location_predictions_out <- bind_rows(lolo_train_test, loc_external_train_val) %>%
  assign_cores(df = ., n_core = n_core, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    # Output list
    predictions_out <- vector("list", length = nrow(core_df))
    
    # Iterate over rows in core_df
    for (i in seq_along(predictions_out)) {
      
      row <- core_df[i,]
      data <- row$train[[1]] %>%
        # Modify the weights
        mutate(site1 = site) %>%
        as.data.frame()
      
      # Test data
      test <- row$test[[1]]
        
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(data, nSite = n_distinct(site1), nObs = n())
      
      ##############
      ##############
      
      # Create a df of model and covariate combinations
      models_df <- tibble(source = "daymet", feature_selection = "none", model = paste0("model", c(1, 4, 5)),
                          covariates = list(tibble(covariate = character()))) %>%
        bind_rows(., covariates_row) %>%
        mutate(submodel = str_replace_all(model, c("4" = "2", "5" = "3")),
               submodel = ifelse(feature_selection == "none", paste0(submodel, "_id"), paste0(submodel, "_cov")),
               submodel = ifelse(submodel == "model1_id", "model1", submodel),
               out = list(NULL))
      
      # Iterate over rows of models
      for (r in seq_len(nrow(models_df))) {
        
        # Get the submodel, declare relationship matrices
        submodel <- models_df$submodel[[r]]
        model_r <- str_remove(submodel, "_id|_cov")
        
        # Extract the variance components
        varcomp_r <- loc_varcomp_df %>% 
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
          
          # Covariate time frame
          time_frame_use <- subset(time_frame_use_df, trait == row$trait, time_frame, drop = TRUE)
          
          
          # Grab the covariates (main effect and interaction covariates are listed together)
          covariate_list <- models_df$covariates[[r]] %>%
            mutate(group = ifelse(str_detect(covariate, ":"), "interaction", "main")) %>%
            split(.$group) %>%
            map(~mutate(.x, covariate = str_remove(covariate, "line_name:")))
          
          # Create main matrix of scaled and centered covariate values
          covariate_mat <- ec_mat_tomodel %>%
            filter(time_frame == time_frame_use, source == src) %>%
            select(-source, -time_frame) %>%
            filter(location %in% levels(data$site1)) %>%
            as.data.frame() %>%
            column_to_rownames("location") %>%
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
        if (model_r == "model4") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain)))
          
        } else if (model_r == "model5") {
          random <- random %>%
            modify_at("site1", ~c(.x, list(K = Emain))) %>%
            modify_at("line_name:site1", ~c(.x, list(K = GE)))
          
        }
        
        # Fixed and random formula
        model_run <- map(model.list, submodel)
        fixed <- model_run$fixed
        
        # Modify the training data
        train <- data %>%
          select(line_name, site1, value) %>%
          mutate(line_name = factor(line_name, levels = row.names(K)),
                 site1 = factor(site1, levels = row.names(Emain)))
    
        # Fit the model
        fit <- fit_lmm(fixed = fixed, random_cov = random, r_cov = varcomp_r["units"], data = train,
                       verbose = TRUE)
        
        # Generate predictions
        predictions <- predict_lmm(object = fit, terms = names(random)) %>%
          # Join with the test set
          inner_join(test, ., by = intersect(names(random), names(test))) %>%
          select(trait, line_name, site, value, pred, pev, reliability, predicted_value)
        
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
location_varcomp_predictions_out <- as_tibble(location_predictions_out)



save("location_varcomp_predictions_out", file = file.path(result_dir, "location_varcomp_predictions.RData"))

















