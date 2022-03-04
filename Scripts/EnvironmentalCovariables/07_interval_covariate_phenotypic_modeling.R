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

concurrent_stepwise_interval_feature_selection <- concurrent_stepwise_feature_selection_list %>%
  ungroup() %>%
  select(-data) %>%
  unnest(cols = out) %>%
  gather(selection_type, covariates, adhoc, adhoc_nosoil) %>% 
  unite(feat_sel_type, feat_sel_type, selection_type, sep = "_")




## Create a data.frame of all covariates
concurrent_interval_all_features <- trait_covariate_df %>%
  group_by(trait) %>% 
  nest(.key = "covariates") %>% 
  mutate(covariates = map(covariates, "covariate")) %>%
  crossing(., source = names(ec_tomodel_interval_centered), feat_sel_type = "all", model = c("model2", "model3")) %>%
  mutate(covariates = modify_if(covariates, model == "model3", ~c(., paste0("line_name:", .))),
         covariates = map(covariates, ~list(optVariables = .)))
           

## Save

save("concurrent_stepwise_interval_feature_selection", "concurrent_interval_all_features", 
     file = file.path(result_dir, "concurrent_interval_feature_selection_results.RData"))

