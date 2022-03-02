## S2MET Genomewide Environment Predictions
## 
## Phenotypic modeling with covariates
## 
## 

# Repository directory
repo_dir <- getwd()
# Source the main project script
source(file.path(repo_dir, "source.R"))

library(paletteer)
library(cowplot)


## Add new packages to load
pkgs <- union(pkgs, c("modelr", "broom", "lme4", "car", "patchwork", "caret", "glmnet"))
# Load these packages
invisible(lapply(X = pkgs, library, character.only = TRUE))

## significance level
alpha <- 0.05



# Load data ---------------------------------------------------------------

# Load covariates for environments and historical covariates
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


## Filter BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp,
         environment %in% train_test_env) %>% # only model train/test environments
  mutate_at(vars(line_name, environment), as.factor)

# Vector of covariate names
covariate_names <- names(ec_tomodel_interval_centered$daymet)[-1:-2]


# Modeling ----------------------------------------------------------------


# Concurrent environmental covariates =====================================


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = covariate_names) %>%
  filter(!(str_detect(covariate, "bulk_density"))) 



# Stepwise feature selection ----------------------------------------------

## Use some feature selection procedures to identify covariates
## 
## Wrap the feature selection within cross-validation to avoid selection
## bias
## 
## Use lm with CV - implemented through recursive feature addition
## 

concurrent_stepwise_feature_selection_list <- S2_MET_BLUEs_tomodel %>%
  crossing(., source = names(ec_tomodel_interval_centered)) %>%
  group_by(trait, source) %>%
  nest() %>%
  mutate(out = list(NULL))

first_null <- min(which(sapply(concurrent_stepwise_feature_selection_list$out, is.null)))

for (i in seq(first_null, nrow(concurrent_stepwise_feature_selection_list))) {
    
    df <- concurrent_stepwise_feature_selection_list$data[[i]] %>%
      mutate(trait = concurrent_stepwise_feature_selection_list$trait[i])
    src <- concurrent_stepwise_feature_selection_list$source[i]
    
    # Factorize
    df1 <- df %>%
      # filter(location != "Aberdeen") %>%
      droplevels() %>%
      left_join(., ec_tomodel_interval_centered[[src]], by = "environment") %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    

    loo_indices <- df1 %>%
      group_by(environment) %>%
      crossv_loo_grouped() %>%
      pull(train) %>%
      map("idx")
    
    # Vector of covariates for this trait
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    
    ## Recursive feature addition
    rfa_out_df <- select_features_met(data = df1, covariates.use = covariates_use, env.col = "environment", search.method = "stepwise")
    
    concurrent_stepwise_feature_selection_list$out[[i]] <- bind_rows(rfa_out_df)
    
}

concurrent_stepwise_feature_selection <- unnest(concurrent_stepwise_feature_selection_list, out) %>%
  gather(selection_type, covariates, adhoc, adhoc_nosoil) %>% 
  unite(feat_sel_type, feat_sel_type, selection_type, sep = "_")





# Determine feature importance using LASSO --------------------------------

# Calculate environmental means
env_means <- S2_MET_BLUEs_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Factorize
    df1 <- df %>%
      droplevels() %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.))) %>%
      mutate(weight = std_error^2)
    
    # Fit the model
    fit <- lmer(value ~ (1|line_name) + environment, data = df1, weights = weight)
    
    ## Return a df of environmental effects
    fixef(fit) %>% 
      tibble(environment = names(.), effect = .) %>% 
      filter(environment != "(Intercept)") %>% 
      mutate(environment = str_remove(environment, "environment"),
             mu = fixef(fit)[1]) %>% 
      add_row(environment = last(levels(df1$environment)), effect = -sum(.$effect), 
              mu = .$mu[1])
    
  }) %>% ungroup()





concurrent_feature_importance_list <- env_means %>%
  crossing(., source = names(ec_tomodel_centered)) %>%
  group_by(trait, source) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(concurrent_feature_importance_list))) {
  
  df <- concurrent_feature_importance_list$data[[i]] %>%
    mutate(trait = concurrent_feature_importance_list$trait[i])
  src <- concurrent_feature_importance_list$source[i]
  
  # Factorize
  df1 <- df %>%
    # filter(location != "Aberdeen") %>%
    droplevels() %>%
    left_join(., ec_tomodel_scaled[[src]], by = "environment") %>%
    mutate_at(vars(environment), ~fct_contr_sum(as.factor(.))) %>%
    dplyr::rename(value = effect)
  
  loo_indices <- df1 %>%
    group_by(environment) %>%
    crossv_loo_grouped() %>%
    pull(train) %>%
    map("idx")
  
  # Vector of covariates for this trait
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
  
  ## Include interactions between rainfall and soil variables
  interaction_covariates <- cross(list(str_subset(covariates_use, "water_balance"), str_subset(covariates_use, "soil"))) %>%
    map_chr(~map_chr(., 1) %>% paste0(., collapse = ":"))

  covariates_use <- c(covariates_use, interaction_covariates)
  
  ## Feature importance using the LASSO
  lasso_importance_out <- select_features_met(data = df1, covariates.use = covariates_use, env.col = "environment", search.method = "lasso")
  
  concurrent_feature_importance_list$out[[i]] <- lasso_importance_out
  
}


concurrent_feature_importance <- unnest(concurrent_feature_importance_list, out) %>%
  select(-contains("RMSE")) %>%
  gather(selection_type, covariates, adhoc, adhoc_nosoil) %>% 
  unite(feat_sel_type, feat_sel_type, selection_type, sep = "_")



# ## Test relationship matrix construction
# impt <- concurrent_feature_importance_list$out[[1]]$adhoc[[1]]
# 
# x <- ec_tomodel_scaled$daymet %>%
#   select(-source) %>%
#   as.data.frame() %>%
#   column_to_rownames("environment") %>%
#   as.matrix()
# 
# # Multiply covariate data by importance
# wt <- matrix(data = c(impt), nrow = nrow(x), ncol = nrow(impt), byrow = TRUE, dimnames = list(row.names(x), row.names(impt)))
# 
# # Subset covariates with non-zero wts - multiply weights
# nonzero <- colMeans(wt == 0) != 1
# x1 <- x[,colnames(wt),drop = FALSE][,nonzero, drop = FALSE] * wt[,nonzero, drop = FALSE]
# 
# Emat <- tcrossprod(x1)





## Create a data.frame of all covariates
concurrent_all_features <- trait_covariate_df %>%
  group_by(trait) %>% 
  nest(.key = "covariates") %>% 
  mutate(covariates = map(covariates, "covariate")) %>%
  crossing(., source = names(ec_tomodel_centered), feat_sel_type = "all", model = c("model2", "model3")) %>%
  mutate(covariates = modify_if(covariates, model == "model3", ~c(., paste0("line_name:", .))),
         covariates = map(covariates, ~list(optVariables = .)))
           

## Save

save("concurrent_stepwise_feature_selection", "concurrent_apriori_feature_selection", 
     "concurrent_feature_importance", "concurrent_all_features", 
     file = file.path(result_dir, "concurrent_feature_selection_results.RData"))

