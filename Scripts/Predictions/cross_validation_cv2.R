## S2MET Prediction Models
## 
## Cross validation 2
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



### Cross-validation

# Data to use
cv_data <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  ## Sample environments
  filter(environment %in% sample_env) %>%
  ##
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))


## Data.frame of TP lines for generating CV folds
cv_tp_df <- data.frame(line_name = factor(tp_geno, levels = c(tp_geno, vp_geno)))
vp_df <- tibble(line_name = vp_geno)






## CV2 - prediction of tested training population in tested environments

# Generate skeleton train/test sets
cv2_train_test <- cv_data %>%
  filter(line_name %in% tp_geno) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## 
    cv_data_train <- distinct(df, environment, id) %>%
      group_by(environment)
    
    ## Generate train/test folds per environment
    replicate(n = nCV, expr = {
      do(cv_data_train, crossv_kfold(data = ., k = k))  %>% 
        ungroup() %>% 
        mutate_at(vars(train, test), ~map(., as.data.frame)) %>% 
        split(.$.id) %>% 
        map_df(~do(., tibble(train = list(bind_rows(.$train)), test = list(bind_rows(.$test)), .id = .$.id[1])))
    }, simplify = FALSE) %>%
      map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y))
                         
    
  }) %>% ungroup()



## Iterate over trait-environment combinations

cv2_predictions <- cv2_train_test %>%
  group_by(trait, rep) %>% nest() %>%
  assign_cores(n_core = n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    ## List for accuracy
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
    
      df <- core_df[[3]][[i]]
      
      ## List of training sets
      train_list <- map(df$train, ~left_join(., cv_data, by = c("environment", "id")))
      test_list <- map(df$test, ~left_join(., cv_data, by = c("environment", "id")))
      
      ## Model 1 - fixed environment, random genotypic main effect
      m1_out_list <- map2(.x = train_list, .y = test_list, ~model1(train = .x, test = .y, Kg = K))
      ## Model 2 - fixed environment, random genotypic main effect, random GxE
      m2_out_list <- map2(.x = train_list, .y = test_list, ~model2(train = .x, test = .y, Kg = K))
      
      ## Combine, calculate accuracy across reps
      m1_cv_acc <- m1_out_list %>%
        map("pgv") %>%
        map_df(~subset(., line_name %in% tp_geno)) %>%
        group_by(environment) %>% 
        summarize(accuracy = cor(value, pred_value)) %>%
        mutate(scheme = "cv2", model = "M1")
      
      ## Combine, calculate accuracy across reps
      m2_cv_acc <- m2_out_list %>%
        map("pgv") %>%
        map_df(~subset(., line_name %in% tp_geno)) %>%
        group_by(environment) %>% 
        summarize(accuracy = cor(value, pred_value)) %>%
        mutate(scheme = "cv2", model = "M2")
      
      ## Combine and return
      out[[i]] <- bind_rows(m1_cv_acc, m2_cv_acc)
      
    }
    
    ## Add out to the core_df and return
    core_df %>%
      select(-core) %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "cv2_results.RData"))

