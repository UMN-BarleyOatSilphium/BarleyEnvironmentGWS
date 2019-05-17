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



## Load the distance methods
load(file.path(alt_proj_dir, "Results/distance_method_results.RData"))


## Number of sample environments
n_env_sample <- 10



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




save_list <- ls(pattern = "predictions$")

## Save
save(list = save_list, file = file.path(result_dir, "pov_results_sample.RData"))

