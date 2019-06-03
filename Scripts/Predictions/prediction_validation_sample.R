## S2MET Prediction Models
## 
## Validation samples
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: May 13, 2019
## 

# # Parent-offspring validation and cross-validation
# 
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
n_core <- 8

## Number of CV folds and reps
nCV <- 25
k <- 5


## Load the distance methods
load(file.path(alt_proj_dir, "Results/distance_method_results.RData"))


## Number of sample environments
n_env_sample <- 8



### Parent-offspring validation

# Data to use
pov_data <- S2_MET_BLUEs %>% 
  filter(environment %in% tp_vp_env) %>%
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.)))

## Sample data
set.seed(1512)
sample_env <- pov_data %>%
  distinct(trait, environment) %>%
  group_by(environment) %>%
  filter(n() == 3) %>%
  distinct(environment) %>%
  pull() %>%
  sample(n_env_sample)

pov_data_use <- pov_data %>% 
  filter(environment %in% sample_env) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))


# Data to use
cv_data_use <- S2_MET_BLUEs %>% 
  filter(environment %in% tp_vp_env) %>%
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  filter(environment %in% sample_env) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))


## Subset the correlation matrix for environments
## Also use IPCA-EC
E_mats_use <- env_rank_df %>% 
  filter(trait %in% traits, model %in% c("MYEC_IPCA", "pheno_dist"), set == "complete") %>%
  # filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>%
  filter(!mat_set %in% c("Malosetti", "MalosettiStand")) %>%
  select(trait, model, cov)


## POV00 - prediction of the validation population in untested environments
## Leave-one-environment-out

# Generate skeleton train/test sets
pov00_train_test <- pov_data_use %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    
    ## The training set is just the TP and not the environment in df
    train <- distinct(pov_data_use, trait, line_name, environment) %>% 
      filter(line_name %in% tp_geno, environment != unique(df$environment), trait == unique(df$trait))
    ## The test set is just the VP and ONLY the environment in df
    test <- distinct(pov_data_use, trait, line_name, environment) %>% 
      filter(line_name %in% vp_geno, environment == unique(df$environment), trait == unique(df$trait))
    
    data_frame(val_environment = unique(df$environment), train = list(train), test = list(test))
    
  }) %>% ungroup() %>%
  select(-environment)


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
      train <- left_join(row$train[[1]], pov_data_use, by = c("environment", "line_name", "trait"))
      test <- left_join(row$test[[1]], pov_data_use, by = c("environment", "line_name", "trait"))
      
      ## Model 1 - fixed environment, random genotypic main effect
      m1_out <- model1(train = train, test = test, Kg = K)
      ## Model 2 - fixed environment, random genotypic main effect, and random GxE (identity)
      m2_out <- model2(train = train, test = test, Kg = K)
      
      ## Model 3a - fixed environment, random genotypic main effect, and random GxE (correlation)
      E_mat1 <- subset(E_mats_use, trait == row$trait & model == "pheno_dist", cov, drop = T)[[1]]
      Ke <- E_mat1[unique(c(train$environment, test$environment)), unique(c(train$environment, test$environment))]
      
      ## Environments to exnclude in the relationship matrix (unobserved)
      env_rm <- setdiff(unique(test$environment), unique(train$environment))
      env_keep <- setdiff(unique(train$environment), unique(test$environment))
      Ke_use <- as.matrix(.bdiag(lst = list(Ke[env_keep, env_keep], diag(length(env_rm)))))
      dimnames(Ke_use) <- dimnames(Ke)
      
      m3a_out <- model3(train = train, test = test, Kg = K, Ke = Ke_use)
      
      ## Model 3b - fixed environment, random genotypic main effect, and random GxE (ECs)
      E_mat1 <- subset(E_mats_use, trait == row$trait & model == "MYEC_IPCA", cov, drop = T)[[1]]
      Ke_use <- E_mat1[unique(c(train$environment, test$environment)), unique(c(train$environment, test$environment))]
      
      m3b_out <- model3(train = train, test = test, Kg = K, Ke = Ke_use)
      
      
      ## Add accuracy to out
      out[[i]] <- data.frame(model = c("M1", "M2", "M3A", "M3B"), scheme = "pov00", 
                             accuracy = c(m1_out$accuracy, m2_out$accuracy, m3a_out$accuracy, m3b_out$accuracy))
      
    }
    
    # Add out to core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()


## POV1 - predict the untested VP in tested environments
# Generate skeleton train/test sets
pov1_train_test <- pov_data_use %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## The training set is just the TP
    train <- distinct(pov_data_use, trait, line_name, environment) %>% 
      filter(line_name %in% tp_geno, trait == unique(df$trait))
    ## The test set is just the VP
    test <- distinct(pov_data_use, trait, line_name, environment) %>% 
      filter(line_name %in% vp_geno, trait == unique(df$trait))
    
    data_frame(train = list(train), test = list(test))
    
  }) %>% ungroup()



## Iterate over trait-environment combinations
pov1_predictions <- pov1_train_test %>%
  group_by(trait) %>%
  do({
    
    row <- .
    
    train <- left_join(row$train[[1]], pov_data_use, by = c("environment", "line_name", "trait"))
    test <- left_join(row$test[[1]], pov_data_use, by = c("environment", "line_name", "trait"))
    
    ## Model 1 - fixed environment, random genotypic main effect
    m1_out <- model1(train = train, test = test, Kg = K)
    m2_out <- model2(train = train, test = test, Kg = K)
    
    ## Model 3a - fixed environment, random genotypic main effect, and random GxE (correlation)
    E_mat1 <- subset(E_mats_use, trait == row$trait & model == "pheno_dist", cov, drop = T)[[1]]
    Ke_use <- E_mat1[unique(c(train$environment, test$environment)), unique(c(train$environment, test$environment))]
    
    m3a_out <- model3(train = train, test = test, Kg = K, Ke = Ke_use)
    
    ## Model 3b - fixed environment, random genotypic main effect, and random GxE (ECs)
    E_mat1 <- subset(E_mats_use, trait == row$trait & model == "MYEC_IPCA", cov, drop = T)[[1]]
    Ke_use <- E_mat1[unique(c(train$environment, test$environment)), unique(c(train$environment, test$environment))]
    
    m3b_out <- model3(train = train, test = test, Kg = K, Ke = Ke_use)
    
    ## Combine and calculate accuracy per model and environment
    list(M1 = m1_out, M2 = m2_out, M3A = m3a_out, M3B = m3b_out) %>%
      map("pgv") %>%
      map2_df(.x = ., .y = names(.), ~mutate(.x, model = .y)) %>%
      group_by(model, environment) %>%
      summarize(accuracy = cor(value, pred_value)) %>%
      ungroup() %>%
      mutate(scheme = "pov1")
    
  })



## Data.frame of TP lines for generating CV folds
cv_tp_df <- data.frame(line_name = factor(tp_geno, levels = c(tp_geno, vp_geno)))
vp_df <- data_frame(line_name = vp_geno)


## CV00 - predict the untested TP (and VP) in untest tested environments
# Generate skeleton train/test sets
cv00_train_test <- cv_data_use %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## All environments
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
  group_by(trait, environment, rep) %>%
  do({
    
    df <- .
    
    ## List of training sets
    train_list <- map(df$train, ~left_join(., cv_data_use, by = c("environment", "line_name", "trait")))
    test_list <- map(df$test, ~left_join(., cv_data_use, by = c("environment", "line_name", "trait")))
    
    ## Model 1 - fixed environment, random genotypic main effect
    m1_out_list <- map2(.x = train_list, .y = test_list, ~model1(train = .x, test = .y, Kg = K))
    ## Model 2 - fixed environment, random genotypic main effect, random GxE
    m2_out_list <- map2(.x = train_list, .y = test_list, ~model2(train = .x, test = .y, Kg = K))
    
    ## Model 3a - fixed environment, random genotypic main effect, and random GxE (correlation)
    E_mat1 <- subset(E_mats_use, trait == df$trait[[1]] & model == "pheno_dist", cov, drop = T)[[1]]
    # List of environments
    env_use <- na.omit(c(unique(df$train[[1]]$environment), unique(df$test[[1]]$environment)))
    Ke_use <- E_mat1[env_use, env_use]
    
    m3a_out_list <- map2(.x = train_list, .y = test_list, ~model3(train = .x, test = .y, Kg = K, Ke = Ke_use))
    
    ## Model 3b - fixed environment, random genotypic main effect, and random GxE (ECs)
    E_mat1 <- subset(E_mats_use, trait == df$trait[[1]] & model == "MYEC_IPCA", cov, drop = T)[[1]]
    Ke_use <- E_mat1[env_use, env_use]
    
    m3b_out_list <- map2(.x = train_list, .y = test_list, ~model3(train = .x, test = .y, Kg = K, Ke = Ke_use))
    
    ## CV accuracy
    cv_acc <- list(M1 = m1_out_list, M2 = m2_out_list, M3A = m3a_out_list, M3B = m3b_out_list) %>%
      map(~map_df(., "pgv")) %>%
      map(~filter(., line_name %in% tp_geno, !is.na(value))) %>%
      map2_df(.x = ., .y = names(.), ~mutate(.x, model = .y)) %>%
      group_by(model, environment) %>%
      summarize(accuracy = cor(value, pred_value)) %>%
      ungroup() %>%
      mutate(scheme = "cv00")
    
    ## POCV accuracy
    pocv_acc <- list(M1 = m1_out_list, M2 = m2_out_list, M3A = m3a_out_list, M3B = m3b_out_list) %>%
      map(~map(., "pgv") %>% map_df(., ~filter(., line_name %in% vp_geno, !is.na(value)) %>% summarize(accuracy = cor(value, pred_value)))) %>%
      map2_df(.x = ., .y = names(.), ~mutate(.x, model = .y)) %>%
      group_by(model) %>%
      summarize(accuracy = mean(accuracy)) %>%
      ungroup() %>%
      mutate(scheme = "pocv00")
    
    bind_rows(cv_acc, pocv_acc) %>%
      mutate(environment = df$environment[1])
    
  }) %>% bind_rows()




## CV00 - predict the test TP (and VP) in test environments
# Generate skeleton train/test sets
cv0_train_test <- cv_data_use %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## All environments
    trait_env <- unique(df$environment)
    df1 <- distinct(df, line_name, trait, environment)
    
    tibble(environment = trait_env) %>%
      mutate(train = map(environment, ~subset(df1, environment != .)),
             test = map(environment, ~subset(df1, environment == .)))
    
  }) %>% ungroup()



## Iterate over trait-environment combinations
cv0_predictions <- cv0_train_test %>%
  group_by(trait, environment) %>%
  do({
    
    row <- .
    
    train <- left_join(row$train[[1]], pov_data_use, by = c("environment", "line_name", "trait"))
    test <- left_join(row$test[[1]], pov_data_use, by = c("environment", "line_name", "trait"))
    
    ## Model 1 - fixed environment, random genotypic main effect
    m1_out <- model1(train = train, test = test, Kg = K)
    ## Model 2 - fixed environment, random genotypic main effect, random GxE
    m2_out <- model2(train = train, test = test, Kg = K)
    
    ## Model 3a - fixed environment, random genotypic main effect, and random GxE (correlation)
    E_mat1 <- subset(E_mats_use, trait == row$trait & model == "pheno_dist", cov, drop = T)[[1]]
    # List of environments
    env_use <- na.omit(c(unique(row$train[[1]]$environment), unique(row$test[[1]]$environment)))
    Ke_use <- E_mat1[env_use, env_use]
    
    m3a_out <- model3(train = train, test = test, Kg = K, Ke = Ke_use)
    
    ## Model 3b - fixed environment, random genotypic main effect, and random GxE (ECs)
    E_mat1 <- subset(E_mats_use, trait == row$trait & model == "MYEC_IPCA", cov, drop = T)[[1]]
    Ke_use <- E_mat1[env_use, env_use]
    
    m3b_out <- model3(train = train, test = test, Kg = K, Ke = Ke_use)
    
    ## Combine and calculate accuracy per model and environment
    list(M1 = m1_out, M2 = m2_out, M3A = m3a_out, M3B = m3b_out) %>%
      map("pgv") %>%
      map2_df(.x = ., .y = names(.), ~mutate(.x, model = .y)) %>%
      mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00")) %>%
      group_by(model, environment, scheme) %>%
      summarize(accuracy = cor(value, pred_value)) %>%
      ungroup()
    
  }) %>% bind_rows()




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "prediction_results_sample.RData"))

