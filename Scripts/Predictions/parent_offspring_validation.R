## S2MET Prediction Models
## 
## Parent-offspring validation
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: 17 July 2019
## 

# # Parent-offspring cross-validation



# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# 
# # Other packages
# library(modelr)



# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))
library(modelr)

## Number of cores
n_core <- detectCores()
n_core <- 16





### Parent-offspring validation

# Data to use
pov_data <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  ## Sample environments
  # filter(environment %in% sample_env) %>%
  ##
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))


## POV00 - prediction of the validation population in untested environments
##
## Leave-one-environment-out

# Generate skeleton train/test sets
pov00_train_test <- pov_data %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    
    ## Get the integer ids for the train/test rows
    test_id <- subset(df, line_name %in% vp_geno, id)
    train_id <- subset(pov_data, line_name %in% tp_geno & environment != unique(df$environment) & trait == unique(df$trait), id)

    tibble(train = list(train_id), test = list(test_id))
    
  }) %>% ungroup()


## Iterate over trait-environment combinations
pov00_predictions <- pov00_train_test %>%
  assign_cores(n_core = n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    # list of predictions
    out <- vector("list", nrow(core_df))
    # Iterate over the rows of core_df
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      train <- left_join(row$train[[1]], pov_data, by = c("id"))
      test <- left_join(row$test[[1]], pov_data, by = c("id"))
      
      ## Model 1 - fixed environment, random genotypic main effect
      m1_out <- model1(train = train, test = test, Kg = K)
      m2_out <- model2(train = train, test = test, Kg = K)
      
      ## Add accuracy to out
      out[[i]] <- data.frame(model = c("M1", "M2"), scheme = "pov00", accuracy = c(m1_out$accuracy, m2_out$accuracy))
      
    }
    
    # Add out to core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()
    

## POV1 - predict the untested VP in tested environments
# Generate skeleton train/test sets
pov1_train_test <- pov_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Get the integer ids for the train/test rows
    test_id <- subset(df, line_name %in% vp_geno, id)
    train_id <- subset(pov_data, line_name %in% tp_geno & trait == unique(df$trait), id)

    tibble(train = list(train_id), test = list(test_id))
    
  }) %>% ungroup()

    

## Iterate over trait-environment combinations
pov1_predictions <- pov1_train_test %>%
  group_by(trait) %>%
  do({
    
    row <- .
    
    train <- left_join(row$train[[1]], pov_data, by = "id")
    test <- left_join(row$test[[1]], pov_data, by = "id")
    
    ## Model 1 - fixed environment, random genotypic main effect
    m1_out <- model1(train = train, test = test, Kg = K)
    m2_out <- model2(train = train, test = test, Kg = K)
      
    ## Calculate accuracy per environment
    m1_acc <- m1_out$pgv %>% 
      group_by(environment) %>% 
      summarize(M1 = cor(value, pred_value))
    
    m2_acc <- m2_out$pgv %>% 
      group_by(environment) %>% 
      summarize(M2 = cor(value, pred_value))
    
    full_join(m1_acc, m2_acc, by = "environment") %>% 
      mutate(scheme = "pov1") %>%
      gather(model, accuracy, -environment, -scheme)
    
  })




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "pov_results.RData"))

