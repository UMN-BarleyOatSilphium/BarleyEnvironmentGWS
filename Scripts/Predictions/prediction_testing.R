## S2MET Prediction Models
## 
## Prediction testing script
## 
## Author: Jeff Neyhart
## Last modified: 3 January 2020
## 

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(broom)


## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))
# Rename
ec_model_building <- unified_ec_models

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))




### Tests ####
 
## 1. Test which covariance matrix is best for modeling the main effect of environment
## 2. Test whether estimating variance components using the full dataset is appropriate (instead
## of fitting a new model with each iteration of CV)


## Create a sample dataset to use
## Trait: Grain yield
## Number of environments: 7

set.seed(1250)

data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait == "GrainYield") %>%
  # Sample environments
  filter(environment %in% sample(unique(environment), 10)) %>%
  mutate(.id = seq(nrow(.))) %>%
  # add covariates
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)),
         environment = as.factor(environment))


# Generate skeleton train/test sets
data_train_test <- data_to_model %>%
  group_by(trait) %>%
  do({crossv_loo2(data = group_by(., environment))}) %>%
  ungroup()




# Test 1 - which covariance matrix is best for modeling the main effect of environment?

## Fit full models
data_to_model1 <- data_to_model %>%
  filter(line_name %in% tp_geno)




## Use the model building data.frame to determine covariates etc.
ec_model_building_use <- ec_model_building %>% 
  unnest(final_model) %>% 
  filter(model == "model3_fwd",
         trait %in% unique(data_to_model1$trait))

# Extract the fitted model object formula
ec_model_form <- formula(ec_model_building_use$object[[1]])

# Fixed covariates
fixed_covariates <- all.vars(ec_model_form) %>% 
  subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
                                       pattern = .))))

# Random regression covariates
random_covariates <- all.vars(ec_model_form) %>% 
  subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
                                       pattern = .)))) %>% setdiff(., "line_name")

# All covariates
all_covariates <- union(fixed_covariates, random_covariates)

## Subset covariates and convert to a matrix
ec_mat <- ec_tomodel_scaled %>% 
  select(environment, all_covariates) %>% 
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()


## Create broad relationship matrices
# Use a DF to store

rel_mat_df <- tibble(method = eval(formals(Env_mat)$method), Gmat = list(K)) %>%
  mutate(Emat_main = map(method, ~Env_mat(x = ec_mat[, fixed_covariates, drop = FALSE], method = .)),
         Emat_gxe = map(method, ~Env_mat(x = ec_mat[, random_covariates, drop = FALSE], method = .)),
         GxEmat = map2(.x = Gmat, .y = Emat_gxe, kronecker, make.dimnames = TRUE),
         model = list(NULL))

## Create incidence matrices
mf <- model.frame(value ~ line_name + environment, data = data_to_model1, drop.unused.levels = TRUE)


## Loop over rows
for (i in seq(nrow(rel_mat_df))) {
  
  # Extract the matrices
  G <- rel_mat_df$Gmat[[i]]
  E <- rel_mat_df$Emat_main[[i]]
  GE <- rel_mat_df$GxEmat[[i]]
  
  ## 
  
  # Fit a model
  fit <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = G) + vs(environment, Gu = E) +
                vs(line_name:environment, Gu = GE), data = mf)
  
  ## Return the model
  rel_mat_df$model[[i]] <- fit
  
}

# Grab the loglikelihood
rel_mat_df1 <- rel_mat_df %>% 
  mutate(logLik = map(model, "monitor") %>% map_dbl(., ~max(.[1,])))

## The Rincent2019 method of computing relationshp matrices is best





# Test 2 - Test whether estimating variance components using the full dataset is appropriate (instead
# of fitting a new model with each iteration of CV)

# Use the full model from the previous test
full_model_df <- subset(rel_mat_df1, method == "Rincent2019")

# Get the relationship matrices
G <- full_model_df$Gmat[[1]]
E <- full_model_df$Emat_main[[1]]
GE <- full_model_df$GxEmat[[1]]

# Extract the model
full_mmer_model <- full_model_df$model[[1]]

# Get the variance components
full_mmer_varcomp <- full_mmer_model$sigma
varG <- full_mmer_varcomp$`u:line_name`
varE <- full_mmer_varcomp$`u:environment`
varGE <- full_mmer_varcomp$`u:int`
varR <- full_mmer_varcomp$units

# Get a vector of variance components
full_mmer_varcom_vect <- set_names(x = full_mmer_model$sigmaVector, nm = names(full_mmer_varcomp))


## CV0 and POV0 randomizations
## Remember that the previous variance components were estimated using only the TP
cv0_pov0_rand <- data_train_test %>%
  mutate_at(vars(train, test), ~map(., as.data.frame)) %>%
  mutate(results = list(NULL))


## For each randomization, fit a new model or use the variance components from
## the full model fit
for (i in seq(nrow(cv0_pov0_rand))) {
  
  # Train the test data
  train <- cv0_pov0_rand$train[[i]]
  test <- cv0_pov0_rand$test[[i]]
  
  # Model frame and design matrices
  mf <- model.frame(value ~ line_name + environment, train)
  y <- model.response(mf)
  X <- model.matrix(~ 1, mf)
  Zg <- model.matrix(~ -1 + line_name, mf)
  Ze <- model.matrix(~ -1 + environment, mf)
  Zge <- model.matrix(~ -1 + line_name:environment, mf)
  
  ## Fit a new model
  fit_new <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = G) + vs(environment, Gu = E) +
                    vs(line_name:environment, Gu = GE), data = train)
  
  ## Extract predictions
  fit_new_pred <- randef(fit_new) %>% 
    map("value") %>%
    map(~tibble(term = names(.x), pred = .x), .vars = vars(term)) %>% 
    imap(~`names<-`(.x, c(.y, "pred"))) %>% 
    unname() %>%
    imap(~`names<-`(.x, c(names(.x)[1], paste0("pred", .y)))) %>%
    map(~rename_at(.x, vars(1), ~str_remove(., "u:"))) %>% 
    modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":")) %>%
    # Order by number of cols
    .[order(map_dbl(., ncol), decreasing = TRUE)] %>%
    reduce(left_join) %>%
    ## Sum predictions
    mutate(new_pred = rowMeans(select(., contains("pred")))) %>%
    select(line_name, environment, new_pred)

  
  ## Use the old model
  # Fit a model using previous variance components
  fit_old <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = G, Gt = varG, Gtc = fixm(1)) +
                    vs(environment, Gu = E, Gt = varE, Gtc = fixm(1)) +
                    vs(line_name:environment, Gu = GE, Gt = varGE, Gtc = fixm(1)),
                  data = train, rcov = ~vs(units, Gt = varR, Gtc = fixm(1)))
  
  ## Extract predictions
  fit_old_pred <- randef(fit_old) %>% 
    map("value") %>%
    map(~tibble(term = names(.x), pred = .x), .vars = vars(term)) %>% 
    imap(~`names<-`(.x, c(.y, "pred"))) %>% 
    unname() %>%
    imap(~`names<-`(.x, c(names(.x)[1], paste0("pred", .y)))) %>%
    map(~rename_at(.x, vars(1), ~str_remove(., "u:"))) %>% 
    modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":")) %>%
    # Order by number of cols
    .[order(map_dbl(., ncol), decreasing = TRUE)] %>%
    reduce(left_join) %>%
    ## Sum predictions
    mutate(old_pred = rowMeans(select(., contains("pred")))) %>%
    select(line_name, environment, old_pred)
  
  
  ## Combine with the test df
  test_pred <- left_join(test, fit_new_pred) %>% 
    left_join(., fit_old_pred)
  
  # Add to df
  cv0_pov0_rand$results[[i]] <- test_pred
  
}


## Calculate prediction accuracy
cv0_pov0_rand_results <- cv0_pov0_rand %>%
  unnest(results) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  group_by(trait, environment, pop) %>%
  summarize_at(vars(contains("pred")), ~cor(., value))

# Plot
cv0_pov0_rand_results %>%
  gather(type, acc, contains("pred")) %>%
  ggplot(aes(x = pop, y = acc, color = type)) +
  geom_boxplot()

## Using the variance components from the original model is clearly better


## POV00
pov00_rand <- data_train_test %>%
  mutate_at(vars(train, test), ~map(., as.data.frame)) %>%
  mutate(train = map(train, ~filter(., line_name %in% tp)),
         test = map(test, ~filter(., line_name %in% vp))) %>%
  mutate(results = list(NULL))


## For each randomization, fit a new model or use the variance components from
## the full model fit
for (i in seq(nrow(pov00_rand))) {
  
  # Train the test data
  train <- pov00_rand$train[[i]]
  test <- pov00_rand$test[[i]]
  
  # Model frame and design matrices
  mf <- model.frame(value ~ line_name + environment, train)
  y <- model.response(mf)
  X <- model.matrix(~ 1, mf)
  Zg <- model.matrix(~ -1 + line_name, mf)
  Ze <- model.matrix(~ -1 + environment, mf)
  Zge <- model.matrix(~ -1 + line_name:environment, mf)
  
  ## Fit a new model
  fit_new <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = G) + vs(environment, Gu = E) +
                    vs(line_name:environment, Gu = GE), data = train)
  
  ## Extract predictions
  fit_new_pred <- randef(fit_new) %>% 
    map("value") %>%
    map(~tibble(term = names(.x), pred = .x), .vars = vars(term)) %>% 
    imap(~`names<-`(.x, c(.y, "pred"))) %>% 
    unname() %>%
    imap(~`names<-`(.x, c(names(.x)[1], paste0("pred", .y)))) %>%
    map(~rename_at(.x, vars(1), ~str_remove(., "u:"))) %>% 
    modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":")) %>%
    # Order by number of cols
    .[order(map_dbl(., ncol), decreasing = TRUE)] %>%
    reduce(left_join) %>%
    ## Sum predictions
    mutate(new_pred = rowMeans(select(., contains("pred")))) %>%
    select(line_name, environment, new_pred)
  
  
  ## Use the old model
  # Fit a model using previous variance components
  fit_old <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = G, Gt = varG, Gtc = fixm(1)) +
                    vs(environment, Gu = E, Gt = varE, Gtc = fixm(1)) +
                    vs(line_name:environment, Gu = GE, Gt = varGE, Gtc = fixm(1)),
                  data = train, rcov = ~vs(units, Gt = varR, Gtc = fixm(1)))
  
  ## Extract predictions
  fit_old_pred <- randef(fit_old) %>% 
    map("value") %>%
    map(~tibble(term = names(.x), pred = .x), .vars = vars(term)) %>% 
    imap(~`names<-`(.x, c(.y, "pred"))) %>% 
    unname() %>%
    imap(~`names<-`(.x, c(names(.x)[1], paste0("pred", .y)))) %>%
    map(~rename_at(.x, vars(1), ~str_remove(., "u:"))) %>% 
    modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":")) %>%
    # Order by number of cols
    .[order(map_dbl(., ncol), decreasing = TRUE)] %>%
    reduce(left_join) %>%
    ## Sum predictions
    mutate(old_pred = rowMeans(select(., contains("pred")))) %>%
    select(line_name, environment, old_pred)
  
  
  ## Combine with the test df
  test_pred <- left_join(test, fit_new_pred) %>% 
    left_join(., fit_old_pred)
  
  # Add to df
  pov00_rand$results[[i]] <- test_pred
  
}


## Calculate prediction accuracy
pov00_rand_results <- pov00_rand %>%
  unnest(results) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  group_by(trait, environment, pop) %>%
  summarize_at(vars(contains("pred")), ~cor(., value))


#
pov00_rand_results %>%
  gather(type, acc, contains("pred")) %>%
  ggplot(aes(x = pop, y = acc, color = type)) +
  geom_boxplot()
  







## Test the AMMI correlation methods of selecting covariates


## Use the model building data.frame to determine covariates etc.
ec_model_building_use <- ec_model_building %>% 
  unnest(final_model) %>% 
  filter(model == "model3_ammi",
         trait %in% unique(data_to_model1$trait))

# Extract the fitted model object formula
ec_model_form <- formula(ec_model_building_use$object[[1]])

# Fixed covariates
fixed_covariates <- all.vars(ec_model_form) %>% 
  subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
                                       pattern = .))))

# Random regression covariates
random_covariates <- all.vars(ec_model_form) %>% 
  subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
                                       pattern = .)))) %>% setdiff(., "line_name")

# All covariates
all_covariates <- union(fixed_covariates, random_covariates)

## Subset covariates and convert to a matrix
ec_mat <- ec_tomodel_scaled %>% 
  select(environment, all_covariates) %>% 
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()


## Create broad relationship matrices
# Use a DF to store

rel_mat_df_1 <- tibble(method = eval(formals(Env_mat)$method), Gmat = list(K)) %>%
  mutate(Emat_main = map(method, ~Env_mat(x = ec_mat[, fixed_covariates, drop = FALSE], method = .)),
         Emat_gxe = map(method, ~Env_mat(x = ec_mat[, random_covariates, drop = FALSE], method = .)),
         GxEmat = map2(.x = Gmat, .y = Emat_gxe, kronecker, make.dimnames = TRUE),
         model = list(NULL))

## Create incidence matrices
mf <- model.frame(value ~ line_name + environment, data = data_to_model1, drop.unused.levels = TRUE)


## Loop over rows
for (i in seq(nrow(rel_mat_df_1))) {
  
  # Extract the matrices
  G <- rel_mat_df_1$Gmat[[i]]
  E <- rel_mat_df_1$Emat_main[[i]]
  GE <- rel_mat_df_1$GxEmat[[i]]
  
  ## 
  
  # Fit a model
  fit <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = G) + vs(environment, Gu = E) +
                vs(line_name:environment, Gu = GE), data = mf)
  
  ## Return the model
  rel_mat_df_1$model[[i]] <- fit
  
}

# Grab the loglikelihood
rel_mat_df_11 <- rel_mat_df_1 %>% 
  mutate(logLik = map(model, "monitor") %>% map_dbl(., ~max(.[1,])))

## The Rincent2019 method of computing relationshp matrices is best



## Compare LL
rel_mat_df1
rel_mat_df_11

rel_mat_df1$model[[3]]$sigma
rel_mat_df_11$model[[3]]$sigma









  
  
# # R matrix - inverse
# Rinv <- solve(diag(length(y)) * full_mmer_varcom_vect["units"])
# 
# # Thetas
# theta_g <- ( full_mmer_varcom_vect["units"] / full_mmer_varcom_vect["u:line_name"] ) * ginv(G)
# theta_e <- ( full_mmer_varcom_vect["units"] / full_mmer_varcom_vect["u:environment"] ) * ginv(E)
# theta_ge <- ( full_mmer_varcom_vect["units"] / full_mmer_varcom_vect["u:int"] ) * ginv(GE)
# 
# 
# ## Right hand side of the equation
# rhs <- rbind(
#   t(X) %*% Rinv %*% y,
#   t(Zg) %*% Rinv %*% y,
#   t(Ze) %*% Rinv %*% y,
#   t(Zge) %*% Rinv %*% y
# )
# 
# # Left hand side
# lhs <- rbind(
#   cbind( t(X) %*% Rinv %*% X , t(X) %*% Rinv %*% Zg, t(X) %*% Rinv %*% Ze, t(X) %*% Rinv %*% Zge ),
#   cbind( t(Zg) %*% Rinv %*% X , t(X) %*% Rinv %*% Zg, t(Zg) %*% Rinv %*% Ze, t(Zg) %*% Rinv %*% Zge ),
#   cbind( t(Ze) %*% Rinv %*% X , t(Ze) %*% Rinv %*% Zg, t(Ze) %*% Rinv %*% Ze, t(Ze) %*% Rinv %*% Zge ),
#   cbind( t(Zge) %*% Rinv %*% X , t(Zge) %*% Rinv %*% Zg, t(Zge) %*% Rinv %*% Ze, t(Zge) %*% Rinv %*% Zge )
#   
#   
# )
  
  
  
  
  
  
  
  
  
  
  




      
      train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
        mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
      test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
        select(environment, line_name, value)
      
      ## Get the covariate model
      ec_model_i <- ec_model_building %>% 
        unnest(final_model) %>%
        subset(trait == row$trait & model == "model3", object, drop = TRUE)
      
      
      
      # # Convert the model to fixed
      # ec_model_form <- as.formula(paste0("value ~ 1 + ", paste0(str_remove(attr(terms(formula(ec_model_i)), "term.labels"), "1 \\| "), collapse = " + ")))
      # # Remove environment
      # ec_model_form_no_env <- formula(drop.terms(terms(ec_model_form), 2, keep.response = TRUE))
      
      # Keep same model
      ec_model_form <- formula(ec_model_i[[1]])
      
      ## Fixed covariates
      fixed_covariates <- all.vars(ec_model_form) %>% 
        subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
                                             pattern = .))))
      
      ## Random regression covariates
      random_covariates <- all.vars(ec_model_form) %>% 
        subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
                                             pattern = .)))) %>%
        setdiff(., "line_name")
      
      ## all covariates
      all_covariates <- union(fixed_covariates, random_covariates)
      
      ## Add covariates to the train df
      train1 <- ec_tomodel_centered %>%
        select(., environment, all_covariates) %>%
        left_join(mutate(train, environment = as.character(environment)), ., by = "environment") %>%
        mutate(environment = as.factor(environment),
               line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
      
      test1 <- ec_tomodel_centered %>%
        select(., environment, all_covariates) %>%
        left_join(mutate(test, environment = as.character(environment)), ., by = "environment") %>%
        mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
      
      ## Extract value of covariates
      test_covariates <- setdiff(all_covariates, "line_name")
      
      
      
      
      
      ###################
      ### Model 1 - random genotype (covariate) and fixed environment
      ###################
      
      ## Copy line_name as new factor / variable
      train2 <- bind_cols(train1, rename_all(as.data.frame(rerun(length(test_covariates), train1$line_name)), 
                                             ~paste0("line_name", seq_along(test_covariates)))) %>%
        mutate(environment = factor(environment, levels = c(as.character(row$environment), levels(train1$environment))),
               ge = interaction(line_name, environment, drop = FALSE, sep = ":"))
      
      
      # Random effect formula
      fixed_form1 <- value ~ 1 + environment
      random_form1 <- ~vs(line_name, Gu = K)
      
      ## Fit with sommer
      model1_fit <- mmer(fixed = fixed_form1, random = random_form1, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model1_fixed_prediction <- subset(model1_fit$Beta, Effect == "(Intercept)", Estimate, drop = TRUE)
      
      # Random effects
      model1_prediction <- model1_fit$U %>% 
        subset(., str_detect(names(.), "u:")) %>% 
        map("value") %>%
        map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
        reduce(crossing) %>%
        mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
        select(line_name, prediction) %>% 
        mutate(prediction = prediction + model1_fixed_prediction) %>%
        left_join(test, ., by = c("line_name"))
      
      
      
      ###################
      ### Model 2 - random genotype (K) + fixed environment based on covariate
      ###################
      
      
      ## Create the environmental covariance matrix
      ec_df <- ec_tomodel_centered %>% 
        filter(environment %in% levels(train2$environment)) %>% 
        select(environment, test_covariates) %>%
        as.data.frame() %>%
        column_to_rownames("environment")
      
      
      # Matrix of covariates for the test environment
      fixed_covariate_x <- t(as.matrix(ec_df[as.character(row$environment), fixed_covariates, drop = FALSE]))
      random_covariate_x <- t(as.matrix(ec_df[as.character(row$environment), random_covariates, drop = FALSE]))
      
      # ## Calculate standardized difference between environments for each covariate
      # ec_dist_mat <- map(ec_df, ~{
      #   dist_mat <- as.matrix(dist(.x))
      #   dimnames(dist_mat) <- list(row.names(ec_df), row.names(ec_df))
      #   ## Standardize
      #   dist_mat / diff(range(.x))
      # })
      # 
      # ## Calculate the covariance matrix
      # ## This line will sum each of the same coordinate element in the list of matrices
      # E_mat <- 1 - reduce(ec_dist_mat, `+`)
      
      # Fixed formula
      fixed_form2 <- as.formula(paste0("value ~ 1 +", paste0(fixed_covariates, collapse = " + ")))
      # Random effect formula
      # random_form2 <- ~vs(line_name, Gu = K) + vs(environment, Gu = E_mat)
      random_form2 <- random_form1
      
      ## Fit with sommer
      model2_fit <- mmer(fixed = fixed_form2, random = random_form2, data = train2, date.warning = FALSE)
      
      
      # Predict
      # Fixed effects
      model2_fixed_prediction <- c(model2_fit$Beta$Estimate %*% c(1, fixed_covariate_x))
      
      # Random effects
      model2_prediction <- model2_fit$U %>% 
        subset(., str_detect(names(.), "u:")) %>% 
        map("value") %>%
        map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), "estimate"))) %>%
        reduce(crossing) %>%
        mutate(prediction = rowSums(select(., contains("estimate")))) %>% 
        select(line_name, prediction) %>% 
        mutate(prediction = prediction + model2_fixed_prediction) %>%
        left_join(test, ., by = c("line_name"))
      
      
      
      # ###################
      # ### Model 3 - random genotype (K) + random environment (E) + random regression of K * E
      # ###################
      # 
      # 
      # 
      # ## Calculate GxE matrix
      # GE_mat <- kronecker(K, E_mat, make.dimnames = TRUE)
      # 
      # 
      # ## Add GxE term to formula
      # random_form3 <- add_predictors(random_form2, ~vs(line_name:environment, Gu = GE_mat))
      # 
      # ## Fit with sommer
      # model3_fit <- mmer(fixed = fixed_form1, random = random_form3, data = train2, date.warning = FALSE)
      # 
      # 
      # # Predict
      # # Fixed effects
      # model3_fixed_prediction <- model3_fit$Beta$Estimate[1]
      # 
      # # Random effects
      # model3_random_prediction <- model3_fit$U %>% 
      #   subset(., str_detect(names(.), "u:")) %>% 
      #   map("value") %>%
      #   map2(., names(.), ~tibble(term = names(.x), estimate = .x) %>% `names<-`(., c(str_remove(.y, "^u:"), paste0(.y, "_estimate"))))
      # 
      # model3_prediction <- model3_random_prediction %>%
      #   map(~{
      #     if (str_detect(names(.x)[1], ":")) {
      #       separate(.x, names(.x)[1], str_split(names(.x)[1], ":")[[1]], sep = ":")
      #     } else {
      #       .x
      #     }
      #     
      #   }) %>%
      #   rev(.) %>%
      #   reduce(., .f = left_join) %>%
      #   mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
      #   select(line_name, environment, prediction) %>% 
      #   mutate(prediction = prediction + model3_fixed_prediction) %>%
      #   left_join(test, ., by = c("line_name", "environment"))
      # 
      
      
      
      ###################
      ### Model 4 - random genotype (K) + fixed environment covariate + random regression of covariate on K
      ###################
      
      # Fixed formula
      fixed_form4 <- fixed_form2
      
      # Random effect formula
      # vs(x, y) specifies an interaction
      random_form4 <- formula(paste0("~ vs(line_name, Gu = K) + ", 
                                     paste0("vs(", random_covariates, ", line_name, Gu = K)", collapse = " + ")))
      
      ## Fit with sommer
      model4_fit <- mmer(fixed = fixed_form4, random = random_form4, data = train2, date.warning = FALSE)
      
      ## If the model fails, do not return predictions for this model
      if (length(model4_fit) == 0) {
        
        model4_prediction <- mutate(test, prediction = as.numeric(NA))
        
      } else {
        
        # Predict
        # Fixed effects
        model4_fixed_prediction <- c(model4_fit$Beta$Estimate %*% c(1, fixed_covariate_x))
        
        # Random effects
        model4_random_prediction <- model4_fit$U %>% 
          subset(., str_detect(names(.), ":")) %>% 
          map("value") %>%
          map(~tibble(term = names(.x), estimate = .x)) %>%
          imap(~`names<-`(.x, rev(str_split(.y, pattern = ":")[[1]]))) %>%
          reduce(left_join, by = "line_name") %>%
          rename_at(vars(u), ~"genotype_mean_estimate")
        
        ## Convert random regression coefficients to matrix
        ec_Beta <- model4_random_prediction %>%
          select(-genotype_mean_estimate) %>%
          as.data.frame(.) %>% 
          column_to_rownames("line_name") %>% 
          as.matrix()
        
        ## Predict the genotype-specific value for the target environment
        model4_random_prediction1 <- model4_random_prediction %>%
          mutate(response_estimate = c(ec_Beta %*% random_covariate_x))
        
        model4_prediction <- model4_random_prediction1 %>%
          mutate(prediction = rowSums(select(., contains("_estimate")))) %>% 
          select(line_name, prediction) %>% 
          mutate(prediction = prediction + model4_fixed_prediction) %>%
          left_join(test, ., by = c("line_name"))
        
      }
      
      
      
      
      
      ###################
      
      ###################
      
      ## Combine and return the predictions
      out[[i]] <- bind_rows(
        mutate(model1_prediction, model = "model1"),
        mutate(model2_prediction, model = "model2"),
        # mutate(model3_prediction, model = "model3"),
        mutate(model4_prediction, model = "model4")
      )
      
      
      ## Notify user
      cat("\nPredictions for trait", row$trait, "in environment", as.character(row$environment), "complete.")
      
    } # CLose loop
    
    ## Add results to the core_df
    core_df %>%
      mutate(out = out) %>%
      unnest(out)
    
  }) %>% bind_rows()
