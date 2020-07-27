## S2MET Prediction Models
## 
## Environment-specific predictions
## 
## Leave-one-environment-out prediction
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

# If re-running predictions, should all be re-run?
rerun_all <- FALSE



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



# Environments ================================================

# Leave-one-out -----------------------------------------------

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
loeo_train_test <- data_to_model %>%
  group_by(trait) %>%
  mutate(site = env) %>%
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



## Use previous runs to determine which scenarios to predict ##

if (!rerun_all) {
  # Read in the data
  load(file.path(result_dir, "loeo_predictions_fact_reg.RData"))
  
  # Find those scenarios that failed to run
  if (!is.data.frame(loeo_predictions_out)) {
    loeo_predictions_out <- loeo_predictions_out %>%
      subset(., !sapply(., is.null)) %>%
      bind_rows()
  }
  
  # Which scenarios were already run?
  scenarios_ran <- select(loeo_predictions_out, trait, site, .id)
  # Which need to be run?
  scenarios_torun <- anti_join(select(loeo_train_test, trait, site, .id), scenarios_ran)
  
} else {
  scenarios_torun <- select(loeo_train_test, trait, site, .id)
    
}


## Pull the scenarios to run, assign cores and split
data_train_test1 <- loeo_train_test %>%
  inner_join(scenarios_torun, .) %>%
  assign_cores(df = ., n_core = n_core, split = TRUE)


## Parallelize
loeo_predictions_out <- data_train_test1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nSite = n_distinct(environment), nObs = n())
      
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
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
        
        # Add to df
        covariates_use$out[[r]] <- prediction_out_cov$prediction_out
        
      }
      
      prediction_out_features <- unnest(covariates_use, out)
      
      
      # ##############
      # ##############
      # 
      # # Use stepwise covariates for E and all covariates for GxE
      # # Subset models
      # models_run <- lapply(X = model.list, "[", c("model3_cov_all"))
      # # Get covariates
      # covariates_use <- covariates_row %>%
      #   select(-direction) %>%
      #   unite(feature_selection_model, feature_selection, model) %>%
      #   filter(feature_selection_model %in% c("rfa_cv_adhoc_model2", "all_model3")) %>%
      #   spread(feature_selection_model, covariates) %>%
      #   mutate(out = list(NULL))
      # 
      # # Iterate over rows
      # for (r in seq_len(nrow(covariates_use))) {
      #   
      #   # Covariate data source
      #   src <- covariates_use$source[r]
      #   
      #   # List of covariates
      #   covariate_list <- list(
      #     main = str_subset(string = covariates_use$rfa_cv_adhoc_model2[[1]], pattern = "line_name:", negate = TRUE),
      #     interaction = str_remove(str_subset(string = covariates_use$all_model3[[1]], pattern = "line_name:"), "line_name:")
      #   )
      #   
      #   # Create a matrix of scaled and centered covariates
      #   covariate_mat <- ec_tomodel_scaled[[src]] %>%
      #     filter(environment %in% levels(row$train[[1]]$site1)) %>%
      #     select(., environment, unique(unlist(covariate_list))) %>%
      #     as.data.frame() %>%
      #     column_to_rownames("environment") %>%
      #     as.matrix()
      #   
      #   ## Environmental relationship matrices
      #   Emain <- Env_mat(x = covariate_mat[,covariate_list$main, drop = FALSE], method = "Jarq")
      #   Eint <- Env_mat(x = covariate_mat[,covariate_list$int, drop = FALSE], method = "Jarq")
      # 
      #   # run predictions
      #   prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
      #   
      #   # Add to df
      #   covariates_use$out[[r]] <- prediction_out_cov$prediction_out
      #   
      # }
      # 
      # prediction_out_model3_cov_all <- unnest(covariates_use, out) %>%
      #   mutate(feature_selection = "rfa_cv_adhoc") %>% 
      #   select(source, feature_selection, trait, model, prediction)
      # 
      # ###################
      # 
      # ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(bind_rows(prediction_out_id$prediction_out, prediction_out_features),
                         train_n = list(train_n)) %>% unnest(train_n)
      # out[[i]] <- mutate(bind_rows(prediction_out_features, prediction_out_model3_cov_all), 
      #                    train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in site", as.character(row$site), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) # %>% bind_rows()




# External environment validation -----------------------------------------------

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


# Generate skeleton train/test sets for validation
env_external_train_val <- data_to_model %>%
  group_by(trait) %>%
  mutate(site = env) %>%
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

## Assign cores and split
env_external_train_val1 <- env_external_train_val %>%
  assign_cores(df = ., n_core = 1, split = TRUE)


## Don't parallelize
env_external_predictions_out <- env_external_train_val1 %>%
  coreApply(X = ., FUN = function(core_df) {
    
    ## Output list
    out <- vector("list", nrow(core_df))
    
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      # Get covariates
      covariates_row <- row$model_covariates[[1]]
      
      # Record the number of environment and observations used for training
      train_n <- summarize(row$train[[1]], nSite = n_distinct(environment), nObs = n())
      
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
        prediction_out_cov <- genomewide_prediction2(x = row, model.list = models_run, K = K, E = Emain, KE = Eint)
        
        # Add to df
        covariates_use$out[[r]] <- prediction_out_cov$prediction_out
        
      }
      
      prediction_out_cov <- unnest(covariates_use, out)
      
      
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- mutate(bind_rows(prediction_out_id$prediction_out, prediction_out_cov), train_n = list(train_n)) %>% unnest(train_n)
      # out[[i]] <- mutate(prediction_out_cov, train_n = list(train_n)) %>% unnest(train_n)
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    mutate(core_df, out = out)
    
  }) %>% bind_rows()




## Save the results
save("loeo_predictions_out", "env_external_predictions_out", 
     file = file.path(result_dir, "loeo_predictions_fact_reg.RData"))





# ## Combine results into one file
# if (file.exists(file.path(result_dir, "loeo_predictions_fact_reg_supp.RData"))) {
#   
#   load(file.path(result_dir, "loeo_predictions_fact_reg.RData"))
#   loeo_predictions_out_reg <- loeo_predictions_out
#   env_external_predictions_out_reg <- env_external_predictions_out
# 
#   # Load the supplemental
#   load(file.path(result_dir, "loeo_predictions_fact_reg_supp.RData"))
#   loeo_predictions_out_supp <- loeo_predictions_out
# 
#   ## Combine
#   loeo_predictions_out <- loeo_predictions_out_reg %>%
#     subset(., !sapply(., is.null)) %>%
#     bind_rows() %>%
#     bind_rows(., bind_rows(loeo_predictions_out_supp)) %>%
#     arrange(trait, site)
#   # Rename
#   env_external_predictions_out <- env_external_predictions_out_reg
# 
#   # Resave
#   save("loeo_predictions_out", "env_external_predictions_out",
#        file = file.path(result_dir, "loeo_predictions_fact_reg.RData"))
# }





