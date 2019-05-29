## S2MET Prediction Models
## 
## Cross validation 0 and 00
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: May 21, 2019
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


## CV00 - prediction of untested TP (and VP) lines in untested environments

# Generate skeleton train/test sets
cv00_train_test <- cv_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## List of unique environments for the trait
    trait_env <- unique(df$environment)
    df1 <- distinct(df, line_name, trait, environment)
    
    ## Generate 25 randomizations per environment
    map(trait_env, ~{
      env <- .
      df_train <- subset(df1, environment != env)
      df_test <- subset(df1, environment == env)
      
      replicate(n = nCV, expr = crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>% 
        map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
        mutate_at(vars(train, test), ~map(., as.data.frame)) %>%
        mutate(test = map(test, ~bind_rows(., vp_df))) %>%
        mutate(train = map(train, ~left_join(., df_train, by = "line_name")),
               test = map(test, ~left_join(., df_test, by = "line_name")))
    }) %>%
      map2_df(.x = ., .y = trait_env, ~mutate(.x, environment = .y))
    
  }) %>% ungroup()
  
  
  


## Iterate over trait-environment combinations
cv00_predictions <- cv00_train_test %>%
  group_by(trait, rep, environment) %>% nest() %>%
  assign_cores(n_core = n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    ## List for accuracy
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      df <- core_df[[4]][[i]]
      
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
        filter(!is.na(value)) %>%
        summarize(accuracy = cor(value, pred_value)) %>%
        mutate(scheme = "cv00", model = "M1")
      
      m1_pocv_acc <- m1_out_list %>%
        map("pgv") %>%
        map(~subset(., line_name %in% vp_geno)) %>%
        map_df(~group_by(., environment) %>% summarize(M1 = cor(value, pred_value))) %>%
        summarize(accuracy = mean(M1)) %>%
        mutate(scheme = "pocv00", model = "M1")
      
      ## Combine, calculate accuracy across reps
      m2_cv_acc <- m2_out_list %>%
        map("pgv") %>%
        map_df(~subset(., line_name %in% tp_geno)) %>%
        filter(!is.na(value)) %>%
        summarize(accuracy = cor(value, pred_value)) %>%
        mutate(scheme = "cv00", model = "M2")
      
      m2_pocv_acc <- m2_out_list %>%
        map("pgv") %>%
        map(~subset(., line_name %in% vp_geno)) %>%
        map_df(~summarize(., M2 = cor(value, pred_value))) %>%
        summarize(accuracy = mean(M2)) %>%
        mutate(scheme = "pocv00", model = "M2")
      
      ## Combine and return
      out[[i]] <- bind_rows(m1_cv_acc, m1_pocv_acc, m2_cv_acc, m2_pocv_acc)
      
    }
    
    ## Add out to the core_df and return
    core_df %>%
      select(-core, -data) %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "cv00_results.RData"))

