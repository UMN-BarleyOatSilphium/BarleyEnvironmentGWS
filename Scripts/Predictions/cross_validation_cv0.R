## S2MET Prediction Models
## 
## Cross validation 0
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: 15 July 2019
## 


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

## Number of folds
k <- 5
# Number of reps
nCV <- 25



### Parent-offspring validation

# Data to use
cv_data <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  ## Sample environments
  # filter(environment %in% sample_env) %>%
  ##
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))




## CV0 - prediction of tested TP (and VP) lines in untested environments

# Generate skeleton train/test sets
cv0_train_test <- cv_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## List of unique environments for the trait
    trait_env <- unique(df$environment)

    ## Generate LOEO randomization
    tibble(environment = trait_env) %>%
      mutate(train = map(trait_env, ~subset(df, environment != ., id)),
             test = map(trait_env, ~subset(df, environment == ., id)))
    
  }) %>% ungroup()
  
  
  


## Iterate over trait-environment combinations
cv0_predictions <- cv0_train_test %>%
  assign_cores(n_core = n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    ## List for accuracy
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      row <- core_df[i,]

      train <- left_join(row$train[[1]], cv_data, by = c("id"))
      test <- left_join(row$test[[1]], cv_data, by = c("id"))
      
      ## Model 1 - fixed environment, random genotypic main effect
      m1_out <- model1(train = train, test = test, Kg = K)
      ## Model 2 - fixed environment, random genotypic main effect, random GxE
      m2_out <- model2(train = train, test = test, Kg = K)
      
      
      ## Calculate accuracy
      m1_acc <- m1_out$pgv %>%
        mutate(scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")) %>%
        group_by(scheme) %>%
        summarize(accuracy = cor(pred_value, value)) %>%
        mutate(model = "M1")
      
      m2_acc <- m2_out$pgv %>%
        mutate(scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")) %>%
        group_by(scheme) %>%
        summarize(accuracy = cor(pred_value, value)) %>%
        mutate(model = "M2")
      
      ## Combine and return
      out[[i]] <- bind_rows(m1_acc, m2_acc)
      
    }
    
    ## Add out to the core_df and return
    core_df %>%
      select(-core, -train, -test) %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "cv0_results.RData"))

