## Predictions using different environmental covariance matrices
## 
## Author: Jeff Neyhart
## Last Updated: November 20, 2018
## 
## This script will look at prediction accuracies from reaction norm models of GxE using different
## environmental covariance matrices. We will use LOEO and 60-40 environment predictions.
## 


### Run on MSI
# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))


# # Run the source script
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(modelr)


# Number of cores
n_core <- 16
n_core <- detectCores()


# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)


## Construct different covariance matrices using distance, covariates, or phenotypic correlation
env_cov_mats <- env_rank_df %>% 
  select(set, trait, model, env_cov_mat = cov) %>%
  arrange(trait, model)

## Subset the matrices for environments in which both the TP and VP were phenotyped
env_cov_mats_use <- env_cov_mats %>%
  filter(! model %in% c("great_circle_dist", "pheno_loc_dist")) %>%
  mutate(env_cov_mat = map(env_cov_mat, ~.[row.names(.) %in% tp_vp_env, colnames(.) %in% tp_vp_env]),
         # Create the GxE covariance matrices
         env_cov_mat = map(env_cov_mat, ~kronecker(., K, make.dimnames = T)))



## Run different proportions of training environments
## This will use all data
prop_train_env <- c(0.25, 0.5, 0.75)
n_iter <- 10

# Generate training and test sets
environment_mc_samples <- S2_MET_BLUEs_use %>%
  filter(environment %in% tp_vp_env) %>%
  group_by(trait) %>%
  do({
    df <- .
    # Number the rows
    df1 <- mutate(df, row = seq(nrow(df))) %>%
      droplevels()
    envs <- distinct(df1, environment)
    
    # Generate environment samples
    samples <- data_frame(pTrainEnv = prop_train_env) %>%
      mutate(envSample = map(pTrainEnv, ~rerun(.n = n_iter, sample_frac(tbl = envs, size = .)))) %>% 
      unnest() %>% 
      group_by(pTrainEnv) %>% 
      mutate(iter = seq(n_iter)) %>%
      ungroup()
    
    # Generate the samples
    samples %>% 
      mutate(train = map(envSample, ~left_join(., df1, by = "environment") %>% filter(line_name %in% tp_geno)), 
             test = map(envSample, ~anti_join(df1, ., by = "environment") %>% filter(line_name %in% vp_geno))) %>% 
      mutate_at(vars(train, test), ~map(., ~pull(., row) %>% resample(data = df1, .)))
    
  }) %>% ungroup()
    
  



## Run predictions
prediction_model_split <- environment_mc_samples %>%
  assign_cores(n_core) %>%
  split(.$core)



# Use mclapply to parallelize
environment_mc_predictions  <- mclapply(X = prediction_model_split, FUN = function(core_df) {

  # # For local machine
  # i <- 7
  # core_df <- prediction_model_split[[i]]
  # #

  results_out <- vector("list", nrow(core_df))


  # Iterate over rows
  for (i in seq_along(results_out)) {

    # Create model matrices
    mf <- model.frame(value ~ line_name + environment, core_df$train[[i]])
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    X <- model.matrix(~ 1 + environment, droplevels(mf))

    pgvs <- subset(env_cov_mats_use, trait == core_df$trait[i] & set == "complete", env_cov_mat, drop = T) %>%
      map(~mixed.solve(y = y, Z = Z, X = X, K = .)) %>%
      map("u") %>%
      # Combine the PGVs with the observations
      map(~data.frame(term = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
            separate(term, c("environment", "line_name"), sep = ":") %>%
            left_join(as.data.frame(core_df$test[[i]]), ., by = c("environment", "line_name")) %>%
            select(trait, environment, line_name, value, pgv))

    K_mat_pgv <- env_cov_mats_use %>%
      filter(trait == core_df$trait[i], set == "complete") %>%
      mutate(predictions = pgvs) %>%
      select(-env_cov_mat)



    ## Predict just using the mean
    Zg <- model.matrix(~ -1 + line_name, mf)
    pgv_mean <- mixed.solve(y = y, Z = Zg, X = X, K = K)$u %>%
      data.frame(line_name = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
      left_join(as.data.frame(core_df$test[[i]]), ., by = c("line_name")) %>%
      select(trait, environment, line_name, value, pgv)

    # Add to the K_mat df

    # Add to the K list
    results_out[[i]] <- K_mat_pgv %>%
      add_row(set = "complete", trait = unique(.$trait), model = "mean", predictions = list(pgv_mean))

  }

  core_df %>%
    mutate(out = results_out) %>%
    select(trait, pTrainEnv, iter, out)

})

environment_mc_predictions <- bind_rows(environment_mc_predictions)

# Save
save_file <- file.path(result_dir, "env_cov_mat_mc_predictions.RData")
save("environment_mc_predictions", file = save_file)






### Leave-one-environment-out
# Generate training and test sets
environment_loeo_samples <- S2_MET_BLUEs_use %>%
  filter(environment %in% tp_vp_env) %>%
  group_by(trait) %>%
  do({
    df <- .
    # Number the rows
    df1 <- mutate(df, row = seq(nrow(df))) %>%
      droplevels()
    envs <- as.character(unique(df1$environment))

    # Generate environment samples
    samples <- data_frame(testEnv = envs) %>%
      mutate(train = map(testEnv, ~filter(df1, environment != ., line_name %in% tp_geno)),
             test = map(testEnv, ~filter(df1, environment == ., line_name %in% vp_geno))) %>%
      mutate_at(vars(train, test), ~map(., ~pull(., row) %>% resample(data = df1, .)))


  }) %>% ungroup()




## Run predictions
prediction_model_split <- environment_loeo_samples %>%
  assign_cores(n_core) %>%
  split(.$core)

# Use mclapply to parallelize
environment_loeo_predictions  <- mclapply(X = prediction_model_split, FUN = function(core_df) {

  # # For local machine
  # i <- 7
  # core_df <- prediction_model_split[[i]]
  # #

  results_out <- vector("list", nrow(core_df))


  # Iterate over rows
  for (i in seq_along(results_out)) {

    # Create model matrices
    mf <- model.frame(value ~ line_name + environment, core_df$train[[i]])
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    X <- model.matrix(~ 1 + environment, droplevels(mf))

    # Fit the model
    pgvs <- subset(env_cov_mats_use, trait == core_df$trait[i] & set == "complete", env_cov_mat, drop = T) %>%
      map(~mixed.solve(y = y, Z = Z, X = X, K = .)) %>%
      map("u") %>%
      # Combine the PGVs with the observations
      map(~data.frame(term = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
            separate(term, c("environment", "line_name"), sep = ":") %>%
            left_join(as.data.frame(core_df$test[[i]]), ., by = c("environment", "line_name")) %>%
            select(trait, environment, line_name, value, pgv))

    K_mat_pgv <- env_cov_mats_use %>%
      filter(trait == core_df$trait[i], set == "complete") %>%
      mutate(predictions = pgvs) %>%
      select(-env_cov_mat)
    
    ## Predict just using the mean
    Zg <- model.matrix(~ -1 + line_name, mf)
    pgv_mean <- mixed.solve(y = y, Z = Zg, X = X, K = K)$u %>%
      data.frame(line_name = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
      left_join(as.data.frame(core_df$test[[i]]), ., by = c("line_name")) %>%
      select(trait, environment, line_name, value, pgv)

    
    # Add to the K list
    results_out[[i]] <- K_mat_pgv %>%
      add_row(set = "complete", trait = unique(.$trait), model = "mean", predictions = list(pgv_mean))
  }


  core_df %>%
    mutate(out = results_out) %>%
    select(trait, testEnv, out)

})

environment_covmat_loeo_predictions <- bind_rows(environment_loeo_predictions)









### Realistic predictions
# Generate training and test sets
environment_loeo_samples_realistic <- S2_MET_BLUEs_use %>%
  filter(environment %in% tp_vp_env) %>%
  group_by(trait) %>%
  do({
    df <- .
    # Number the rows
    df1 <- mutate(df, row = seq(nrow(df))) %>%
      droplevels()
    
    ## Get the available environments from the covariance matrices
    envs <- env_cov_mats %>% 
      filter(set == "realistic", trait == unique(df$trait), ! model %in% c("great_circle_dist", "pheno_loc_dist")) %>% 
      pull() %>% 
      map(rownames) %>% 
      reduce(intersect)
    
    testEnvs <- df1 %>% 
      filter(year == 2017) %>% 
      distinct(environment) %>% 
      pull() %>% 
      as.character() %>%
      intersect(., envs)
    
    trainEnvs <- df1 %>% 
      filter(year != 2017) %>% 
      distinct(environment) %>% 
      pull() %>% 
      as.character() %>%
      intersect(., envs)
    
    ## Filter the dataset again
    df2 <- df1 %>%
      filter(environment %in% c(testEnvs, trainEnvs)) %>%
      droplevels()
    
    
    ## Pull df rows for resampling
    testRows <- df2 %>% 
      filter(environment %in% testEnvs, line_name %in% vp_geno) %>% 
      pull(row)
    
    trainRows <- df2 %>% 
      filter(environment %in% trainEnvs, line_name %in% tp_geno) %>% 
      pull(row)
    
    # Create samples
    data_frame(test = list(resample(df2, testRows)), train = list(resample(df2, trainRows)))
    
  }) %>% ungroup()





# Run locally
environment_loeo_realistic_predictions  <- environment_loeo_samples_realistic %>%
  mutate(predictions = list(NULL))
  
for (i in seq(nrow(environment_loeo_realistic_predictions))) {

  # Create model matrices
  mf <- model.frame(value ~ line_name + environment, environment_loeo_realistic_predictions$train[[i]])
  y <- model.response(mf)
  Z <- model.matrix(~ -1 + line_name:environment, mf)
  X <- model.matrix(~ 1 + environment, droplevels(mf))
  
  # Fit the model
  pgvs <- subset(env_cov_mats_use, trait == environment_loeo_realistic_predictions$trait[i] & set == "realistic", env_cov_mat, drop = T) %>%
    map(~mixed.solve(y = y, Z = Z, X = X, K = .)) %>%
    map("u") %>%
    # Combine the PGVs with the observations
    map(~data.frame(term = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
          separate(term, c("environment", "line_name"), sep = ":") %>%
          left_join(as.data.frame(environment_loeo_realistic_predictions$test[[i]]), ., by = c("environment", "line_name")) %>%
          select(trait, environment, line_name, value, pgv))
  
  ## Predict just using the mean
  Zg <- model.matrix(~ -1 + line_name, mf)
  pgv_mean <- mixed.solve(y = y, Z = Zg, X = X, K = K)$u %>%
    data.frame(line_name = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
    left_join(as.data.frame(environment_loeo_realistic_predictions$test[[i]]), ., by = c("line_name")) %>%
    select(trait, environment, line_name, value, pgv)
  
  environment_loeo_realistic_predictions$predictions[[i]] <- env_cov_mats_use %>%
    filter(trait == environment_loeo_realistic_predictions$trait[i] & set == "realistic") %>%
    mutate(predictions = pgvs) %>%
    select(-env_cov_mat) %>%
    add_row(set = unique(.$set), trait = unique(.$trait), model = "mean", predictions = list(pgv_mean))
  
}


environment_loeo_realistic_predictions <- environment_loeo_realistic_predictions %>% 
  select(-test, -train)


# Save
save_file <- file.path(result_dir, "env_cov_mat_realistic_predictions.RData")
save("environment_loeo_realistic_predictions", file = save_file)

















#### Testing #####


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
library(modelr)
library(broom)


# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))


# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits) %>%
  mutate_at(vars(environment:line_name), as.factor)


## Create training and testing sets for each trait
train_test <- S2_MET_BLUEs_use %>% 
  filter(environment %in% c(tp_vp_env, tp_only_env)) %>% 
  group_by(trait, location) %>% 
  filter(n_distinct(year) > 1) %>%
  group_by(trait) %>%
  do(resample_prediction(data = ., train.exp = "line_name %in% tp_geno & year != 2017", test.exp = "line_name %in% vp_geno & year == 2017")) %>%
  ungroup()
  

## Line mean predictions (ignores GxE)
mean_predictions <- train_test %>%
  group_by(trait) %>%
  do({
    df <- .

    # model.frame
    mf <- model.frame(value ~ line_name + environment + location, data = df$train[[1]])
    y <- model.response(mf)

    # Fixed effect of environment
    X <- model.matrix(~ environment, data = droplevels(mf))
    Zg <- model.matrix(~ -1 + line_name, data = mf)
    colnames(Zg) <- colnames(K)

    ## Fit a mean model
    fit1 <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = K)))
 
    # Extract predictions and calculate accuracy
    predictions <- fit1$u.hat %>% 
      map(fix_data_frame) %>% 
      reduce(left_join, by = "term") %>% 
      rename(line_name = term) %>%
      mutate(pred_value = rowSums(select(., contains("T1")))) %>% 
      select(line_name, pred_value)

    left_join(as_data_frame(df$test[[1]]), predictions) %>% 
      group_by(environment) %>% 
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup()


## Genotype and environment interaction by location distance
ggl_predictions <- train_test %>%
  group_by(trait) %>%
  do({
    df <- droplevels(.) # Drop unused environments for this trait
    
    # model.frame
    mf <- model.frame(value ~ line_name + environment + location, data = df$train[[1]])
    y <- model.response(mf)
    
    # Fixed effect of environment
    X <- model.matrix(~ environment, data = droplevels(mf))
    # Random effect of genotype and gxe
    Zg <- model.matrix(~ -1 + line_name, data = mf)
    colnames(Zg) <- colnames(K)
    
    Zge <- model.matrix(~ -1 + line_name:environment, data = mf)
    ## Create the GxE relationship matrix
    # First the relationship matrix of locations
    Lcov <- subset(env_rank_df, model == "pheno_loc_dist" & set == "realistic" & trait == df$trait, cov, drop = T)[[1]]
    Kge <- kronecker(X = Lcov[levels(mf$environment), levels(mf$environment)], Y = K, make.dimnames = TRUE)
    colnames(Zge) <- colnames(Kge)
    
    
    ## Fit a model
    fit1 <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = K), ge = list(Z = Zge, K = Kge)))
    # fit2 <- sommer::mmer(Y = y, X = X, Z = list(ge = list(Z = Zge, K = Kge)))
    
    
    # Extract predictions and calculate accuracy
    predictions <- fit1$u.hat %>% 
      map(fix_data_frame) %>% 
      map(~separate(., term, c("environment", "line_name"), sep = ":", fill = "left")) %>%
      map(~select_if(., ~!all(is.na(.)))) %>%
      reduce(full_join, by = reduce(map(., names), intersect) %>% .[!str_detect(., "T1")]) %>%
      mutate(pred_value = rowSums(select(., contains("T1")))) %>% 
      select(line_name, environment, pred_value)
    
    left_join(as_data_frame(df$test[[1]]), predictions, by = c("environment", "line_name")) %>% 
      group_by(environment) %>% 
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup()


## Genotype and environment interaction by location distance and ECs
gge_predictions <- train_test %>%
  left_join(., filter(env_rank_df, set == "realistic", str_detect(model, "EC"))) %>%
  group_by(trait, model, mat_set) %>%
  do({
    df <- droplevels(.) # Drop unused environments for this trait
    
    # model.frame
    mf <- model.frame(value ~ line_name + environment, data = df$train[[1]])
    y <- model.response(mf)
    
    # Fixed effect of environment
    X <- model.matrix(~ environment, data = droplevels(mf))
    # Random effect of genotype and gxe
    Zg <- model.matrix(~ -1 + line_name, data = mf)
    colnames(Zg) <- colnames(K)
    
    Zge <- model.matrix(~ -1 + line_name:environment, data = mf)
    ## Create the GxE relationship matrix
    # First the relationship matrix of locations
    Lcov <- subset(env_rank_df, model == "pheno_loc_dist" & set == "realistic" & trait == df$trait, cov, drop = T)[[1]]
    ECcov <- df$cov[[1]]
    # Add the elements
    Ecov <- Lcov + ECcov
    
    Kge <- kronecker(X = Ecov[levels(mf$environment), levels(mf$environment)], Y = K, make.dimnames = TRUE)
    colnames(Zge) <- colnames(Kge)
    
    
    ## Fit a model
    fit1 <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = K), ge = list(Z = Zge, K = Kge)))

    # Extract predictions and calculate accuracy
    predictions <- fit1$u.hat %>% 
      map(fix_data_frame) %>% 
      map(~separate(., term, c("environment", "line_name"), sep = ":", fill = "left")) %>%
      map(~select_if(., ~!all(is.na(.)))) %>%
      reduce(full_join, by = reduce(map(., names), intersect) %>% .[!str_detect(., "T1")]) %>%
      mutate(pred_value = rowSums(select(., contains("T1")))) %>% 
      select(line_name, environment, pred_value)
    
    left_join(as_data_frame(df$test[[1]]), predictions, by = c("environment", "line_name")) %>% 
      group_by(environment) %>% 
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup()



## Model
accuracy_data <- bind_rows(mutate(mean_predictions, model = "G"), mutate(ggl_predictions, model = "GGL"), 
                      unite(gge_predictions, model, model, mat_set) %>% 
                        mutate(model = str_remove_all(model, "_NA"))) %>%
  mutate_at(vars(environment, model), as.factor)

save("accuracy_data", file = file.path(result_dir, "env_cov_mat_testing_results.RData"))

model_fit <- accuracy_data %>%
  filter(model %in% c("G", "GGL") | str_detect(model, "MYEC_All")) %>%
  group_by(trait) %>%
  do(fit = lm(accuracy ~ environment + model, data = .))

# Effects
effects::Effect("model", model_fit$fit[[2]]) %>% plot




### Test disconnected CV versus VP predictions
## Sample environments
sample_envs <- sample(x = reduce(complete_train_env, intersect), size = 10)


# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait == "GrainYield",
         environment %in% sample_envs) %>%
  mutate_at(vars(environment:line_name), as.factor)

## Subset the correlation matrix
env_cov_use <- filter(env_rank_df, model == "pheno_dist", trait %in% unique(S2_MET_BLUEs_use$trait)) %>% 
  mutate(cov = map(cov, ~.[sample_envs, sample_envs])) %>% 
  select(trait, cov)

# Number of cv iterations
cv_iter <- 10
# Testing set size
cv_test_size <- length(vp_geno)


## Create training/testing cv sets
cv_test <- replicate(n = cv_iter, sample(tp_geno, size = cv_test_size), simplify = FALSE)
cv_train <- map(cv_test, ~setdiff(tp_geno, .))


## Create training and testing sets for each trait
train_test <- distinct(S2_MET_BLUEs_use, trait, environment) %>% 
  crossing(., data_frame(iter = seq(cv_iter), cv_train_set = cv_train, cv_test_set = cv_test)) %>%
  mutate(cv_train = map2(environment, cv_train_set, ~filter(S2_MET_BLUEs_use, line_name %in% .y, environment != .x)),
         cv_test = map2(environment, cv_test_set, ~filter(S2_MET_BLUEs_use, line_name %in% .y, environment == .x)),
         vp_test = map(environment, ~filter(S2_MET_BLUEs_use, line_name %in% vp_geno, environment == .))) %>%
  select(-contains("set")) %>%
  gather(test_set, data, cv_test, vp_test)
  
  
## Line mean predictions (ignores GxE)
mean_predictions <- train_test %>%
  group_by(trait, environment, iter, test_set) %>%
  do(out = gblup(K = K, train = .$cv_train[[1]], test = .$data[[1]], fit.env = TRUE)) %>%
  ungroup()

## Line x environment predictions using correlations
ge_predictions <- train_test %>%
  group_by(trait, environment, iter, test_set) %>%
  do({
    df <- .
    
    # model.frame
    mf <- model.frame(value ~ line_name + environment, data = df$cv_train[[1]])
    y <- model.response(mf)
    X <- model.matrix(~ 1 + environment, droplevels(mf))
    
    # Kronecker of K and E
    E <- subset(env_cov_use, trait == df$trait, cov, drop = T)[[1]]
    KE <- kronecker(E, K, make.dimnames = TRUE)
    
    # Random effect of g
    Zg <- model.matrix(~ -1 + line_name, mf)
    colnames(Zg) <- colnames(K)
    # Random effect of gxe
    Zge <- model.matrix(~ -1 + line_name:environment, data = mf)
    colnames(Zge) <- colnames(KE)
    # 
    
    ## fit the model
    fit <- mixed.solve(y = y, X = X, Z = Zge, K = KE)
    
    predictions <- data.frame(term = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE) %>% 
      separate(term, c("environment", "line_name"), sep = ":")
    
    # fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = K), ge = list(Z = Zge, K = KE)), silent = TRUE)
    # 
    # predictions <- full_join(x = fix_data_frame(fit$u.hat$g, newcol = "line_name"), y = separate(fix_data_frame(fit$u.hat$ge), term, c("environment", "line_name"), sep = ":"), by = c("line_name")) %>% 
    #   mutate(pred_value = {rowSums(select(., contains("T1")))}) %>% 
    #   select(-contains("T1"))
    
    
    left_join(df$data[[1]], predictions, by = c("environment", "line_name")) %>% 
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup()



## Combine results
predictions_out <- bind_rows(
  mutate(mean_predictions, accuracy = map_dbl(out, "accuracy"), model = "G"),
  mutate(ge_predictions, model = "GE")
) %>% select(-out)


## model
## 
predictions_out1 <- predictions_out %>%
  mutate_at(vars(model, test_set, iter), as.factor)



fit <- lmer(accuracy ~ model + test_set + model:test_set + (1|environment) + (1|environment:test_set) + (1|iter), data = predictions_out1)

plot(effects::allEffects(fit))



