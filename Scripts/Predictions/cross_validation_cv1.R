## S2MET Prediction Models
## 
## Cross validation 1
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: May 13, 2019
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
  filter(environment %in% tp_vp_env) %>%
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.)))

# ## Sample data
# set.seed(1512)
# sample_env <- cv_data %>% 
#   distinct(trait, environment) %>% 
#   group_by(environment) %>% 
#   filter(n() == 3) %>% 
#   distinct(environment) %>% 
#   pull() %>% 
#   sample(5)

cv_data <- cv_data %>% 
  # filter(environment %in% sample_env) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))

## Data.frame of TP lines for generating CV folds
cv_tp_df <- data.frame(line_name = factor(tp_geno, levels = c(tp_geno, vp_geno)))
vp_df <- data_frame(line_name = vp_geno)


## CV1 - prediction of the training (and validation) population in tested environments

# Generate skeleton train/test sets
cv1_train_test <- cv_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## The test set is just the VP
    cv_data_test <- distinct(cv_data, trait, line_name, environment) %>% 
      filter(trait == unique(df$trait))
    cv_data_train <- filter(cv_data_test, line_name %in% tp_geno)
    
    ## Generate train/test folds
    replicate(n = nCV, expr = crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>% 
      map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
      mutate_at(vars(train, test), ~map(., as.data.frame)) %>%
      mutate(test = map(test, ~bind_rows(., vp_df))) %>%
      mutate(train = map(train, ~left_join(., cv_data_train, by = "line_name")),
             test = map(test, ~left_join(., cv_data_test, by = "line_name")))
    
  }) %>% ungroup()



## Iterate over trait-environment combinations
cv1_predictions <- cv1_train_test %>%
  group_by(trait, rep) %>%
  do({
    
    df <- .
    
    ## List of training sets
    train_list <- map(df$train, ~left_join(., cv_data, by = c("environment", "line_name", "trait")))
    test_list <- map(df$test, ~left_join(., cv_data, by = c("environment", "line_name", "trait")))
    
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
      mutate(scheme = "cv1", model = "M1")
    
    m1_pocv_acc <- m1_out_list %>%
      map("pgv") %>%
      map(~subset(., line_name %in% vp_geno)) %>%
      map_df(~group_by(., environment) %>% summarize(M1 = cor(value, pred_value))) %>%
      group_by(., environment) %>%
      summarize(accuracy = mean(M1)) %>%
      mutate(scheme = "pocv1", model = "M1")
    
    ## Combine, calculate accuracy across reps
    m2_cv_acc <- m2_out_list %>%
      map("pgv") %>%
      map_df(~subset(., line_name %in% tp_geno)) %>%
      group_by(environment) %>% 
      summarize(accuracy = cor(value, pred_value)) %>%
      mutate(scheme = "cv1", model = "M2")
    
    m2_pocv_acc <- m2_out_list %>%
      map("pgv") %>%
      map(~subset(., line_name %in% vp_geno)) %>%
      map_df(~group_by(., environment) %>% summarize(M2 = cor(value, pred_value))) %>%
      group_by(., environment) %>%
      summarize(accuracy = mean(M2)) %>%
      mutate(scheme = "pocv1", model = "M2")
    
    ## Combine and return
    bind_rows(m1_cv_acc, m1_pocv_acc, m2_cv_acc, m2_pocv_acc)
    
  })




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "cv1_results.RData"))

