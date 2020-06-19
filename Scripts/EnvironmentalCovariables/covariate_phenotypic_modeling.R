## S2MET phenotypic modeling with covariates
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## Model the impact of covariates on GxE and LxE
## 
## This script will look at variable selection to determine optimal models that use
## covariates
## 

# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))

library(paletteer)
library(cowplot)


## Add new packages to load
pkgs <- union(pkgs, c("modelr", "broom", "lme4", "car", "patchwork", "caret"))
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


# Modeling ----------------------------------------------------------------


# Concurrent environmental covariates =====================================


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = reduce(map(ec_tomodel_centered, ~names(.)[-1:-2]), union)) %>%
  filter(
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 



# Use the following covariates for all traits (except where timing would make
# it inappropriate):
# 
# radiation during vegetative
# tmin during flowering
# water stress during flowering
# grain fill water stress
# grain fill tmax
# 


## A priori covariates - quite minimal; each should have a citation 
## 
## HeadingDate - everything before flowering
## 
## PlantHeight - everything before grain fill
##  
## Grain yield, grain protein, and test weight
## - Everything for plant height plus...
## - Elevated temperature during grain fill (Passarella et al 2005)
## - Drought during grain fill - (Savin and Nicholas  1996)
## 
## 
## 

apriori_covariate_df <- trait_covariate_df %>%
  filter(str_detect(covariate, "soil", negate = TRUE)) %>%
  {bind_rows(
    filter(., trait %in% c("HeadingDate", "PlantHeight")),
    filter(., trait %in% c("GrainProtein", "GrainYield", "TestWeight"),
           covariate %in% c(subset(., trait == "PlantHeight", covariate, drop = TRUE), 
                            "grain_fill.tmean_mean", "grain_fill.water_balance_sum"))
  )}


## Factorial regression with AIC stepwise selection ##

## Group by trait and model
concurrent_fact_reg <- S2_MET_BLUEs_tomodel %>%
  crossing(., source = names(ec_tomodel_centered)) %>%
  group_by(trait, source) %>%
  do({

    df <- .
    # Factorize
    df1 <- df %>%
      # filter(location != "Aberdeen") %>%
      droplevels() %>%
      left_join(., bind_rows(ec_tomodel_centered), by = c("source", "environment")) %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    
    # Apriori
    ## Add covariates - filter
    covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    apriori_out <- fact_reg(data = df1, covariates = covariates_use, method = "apriori")
    


    # Ad hoc
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step", criterion = "BIC")


    # Ad hoc - without soil
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T),
                             covariate, drop = TRUE)
    adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step", criterion = "BIC")


    ## Return results
    tibble(model = c("base", "base_alt", "model2", "model3"),
           apriori = apriori_out,
           adhoc = adhoc_out,
           adhoc_nosoil = adhoc_nosoil_out)

  }) %>% ungroup()


# ## Assess via LOO cv
# concurrent_fact_reg_loo_test <- S2_MET_BLUEs_tomodel %>%
#   left_join(., bind_rows(ec_tomodel_centered)) %>%
#   group_by(source, trait) %>%
#   do(crossv_loo_grouped(group_by(., environment))) %>%
#   ungroup() %>%
#   left_join(., concurrent_fact_reg) %>%
#   filter(str_detect(model, "base", negate = TRUE)) %>%
#   filter(trait == "GrainProtein") %>%
#   mutate(train_fit = map2(.x = train, .y = adhoc, ~update(.y, data = .x)),
#          test_pred = map2(.x = test, .y = train_fit, ~add_predictions(as.data.frame(.x), .y)))
# 
# concurrent_fact_reg_loo_test %>% 
#   unnest(test_pred) %>%
#   qplot(x = pred, y = value, color = environment, data = .) +
#   facet_wrap(~ source + model, scales = "free")


## Prepare results for saving
concurrent_fact_reg_feature_selection <- concurrent_fact_reg %>%
  mutate(feat_sel_type = "stepAIC", direction = "both") %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  mutate_at(vars(apriori, adhoc, adhoc_nosoil), ~map(., ~list(optVariables = attr(terms(.x), "term.labels"))))
  
  




# Feature selection ############################################################

## Use some feature selection procedures to identify covariates
## 
## Wrap the feature selection within cross-validation to avoid selection
## bias
## 
## Use lm with CV - implemented through recursive feature elimination
## 

concurrent_feature_selection_list <- S2_MET_BLUEs_tomodel %>%
  crossing(., source = names(ec_tomodel_centered)) %>%
  group_by(trait, source) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(concurrent_feature_selection_list))) {
    
    df <- concurrent_feature_selection_list$data[[i]] %>%
      mutate(trait = concurrent_feature_selection_list$trait[i])
    src <- concurrent_feature_selection_list$source[i]
    
    # Factorize
    df1 <- df %>%
      # filter(location != "Aberdeen") %>%
      droplevels() %>%
      left_join(., ec_tomodel_centered[[src]], by = "environment") %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    

    loo_indices <- df1 %>%
      group_by(environment) %>%
      crossv_loo_grouped() %>%
      pull(train) %>%
      map("idx")
    
    ## Recursive feature addition
    rfa_out_df <- select_features_met(data = df1, env.col = "environment", search.method = "hill.climb")
    
    concurrent_feature_selection_list$out[[i]] <- bind_rows(rfa_out_df)
    
}

concurrent_feature_selection <- unnest(concurrent_feature_selection_list, out)



# ## Use a genetic algorithm to find ideal models
# 
# concurrent_genalg_features_list <- S2_MET_BLUEs_tomodel %>%
#   crossing(., source = names(ec_tomodel_centered)) %>%
#   group_by(trait, source) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# for (i in seq_len(nrow(concurrent_genalg_features_list))) {
#   
#   df <- concurrent_genalg_features_list$data[[i]] %>%
#     mutate(trait = concurrent_genalg_features_list$trait[i])
#   src <- concurrent_genalg_features_list$source[i]
#   
#   # Factorize
#   df1 <- df %>%
#     # filter(location != "Aberdeen") %>%
#     droplevels() %>%
#     left_join(., ec_tomodel_centered[[src]], by = "environment") %>%
#     mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
#   
#   
#   loo_indices <- df1 %>%
#     group_by(environment) %>%
#     crossv_loo_grouped() %>%
#     pull(train) %>%
#     map("idx")
#   
#   ## Recursive feature addition
#   rfa_out_df <- select_features_met(data = df1, env.col = "environment", search.method = "hill.climb")
#   
#   concurrent_genalg_features_list$out[[i]] <- bind_rows(rfa_out_df)
#   
# }
# 
# concurrent_genalg_features <- unnest(concurrent_genalg_features_list, out)
# 


## Save

save("concurrent_feature_selection", 
     "concurrent_fact_reg_feature_selection",
     file = file.path(result_dir, "concurrent_feature_selection_results.RData"))




# Historical covariates =====================================================


## Ideas
## 1. Calculate locations means and use that as the input (and validation)
## 2. Use all data and fit a model with random year, location:year, and g:year
## 3. Use all data and fit the g + location model


## Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp)


## Select the historical covariate data
historical_ec_tomodel_centered_use <- historical_ec_tomodel_timeframe_centered %>%
  subset(., map_lgl(names(.), ~str_detect(., "time_frame5_2010_2014")))
names(historical_ec_tomodel_centered_use) <- str_remove(string = names(historical_ec_tomodel_centered_use), pattern = ".time_frame5_2010_2014")


# Create a data.frame of covariates per trait
trait_covariate_df <- historical_ec_tomodel_centered_use %>%
  map(~names(.)[-1:-3]) %>%
  reduce(union) %>%
  crossing(trait = traits, covariate = .) %>%
  filter(
    covariate != "awc_range",
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 



# Use the following covariates for all traits (except where timing would make
# it inappropriate):
# 
# radiation during vegetative
# tmin during flowering
# water stress during flowering
# grain fill water stress
# grain fill tmax


# # Load packages
# invisible(lapply(X = pkgs, library, character.only = TRUE))

historical_fact_reg <- S2_MET_loc_BLUEs_tomodel %>%
  crossing(., source = map_chr(historical_ec_tomodel_centered_use, ~unique(pull(., "source")))) %>%
  group_by(trait, source) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(historical_fact_reg))) {

    df <- historical_fact_reg$data[[i]]
    df$trait <- historical_fact_reg$trait[i]
    
    src <- historical_fact_reg$source[i]

    # Factorize
    df1 <- df %>%
      droplevels() %>%
      left_join(., historical_ec_tomodel_centered_use[[src]], by = "location") %>%
      mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))

    # Apriori
    ## Add covariates - filter
    covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "apriori")

    # Ad hoc
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")

    # Ad hoc - without soil
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T),
                             covariate, drop = TRUE)
    adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")


    ## Return results
    historical_fact_reg$out[[i]] <- tibble(model = c("base", "base_alt", "model4", "model5"),
           apriori = apriori_out,
           adhoc = adhoc_out,
           adhoc_nosoil = adhoc_nosoil_out)

}


historical_fact_reg1 <- unnest(historical_fact_reg, out)


## Assess via LOO cv
historical_fact_reg_loo_test <- historical_fact_reg %>% 
  mutate(loo_indices = map(data, ~crossv_loo_grouped(group_by(., location))) %>%
           map("train") %>% map(., ~map(., "idx")),
         out = map2(out, loo_indices, ~mutate(.x, loo_indices = list(.y)))) %>%
  unnest(out) %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  mutate(adhoc_cv = map2(adhoc, loo_indices, ~rapid_cv(object = .x, index = .y, return.predictions = TRUE)))
  
historical_fact_reg_loo_test %>%
  mutate(out = map(adhoc_cv, "cv_metrics") %>% map(~as.data.frame(as.list(.)))) %>%
  unnest(out)
  


## Prepare results for saving
historical_fact_reg_feature_selection <- historical_fact_reg1 %>%
  mutate(feat_sel_type = "stepAIC", direction = "forward") %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  mutate_at(vars(apriori, adhoc, adhoc_nosoil), ~map(., ~list(optVariables = attr(terms(.x), "term.labels"))))






# Feature selection ############################################################

historical_feature_selection_list <- historical_fact_reg %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(historical_feature_selection_list))) {

    df <- historical_feature_selection_list$data[[i]] %>% 
      mutate(trait = historical_feature_selection_list$trait[i])
    
    src <- historical_fact_reg$source[i]
    
    
    # Factorize
    df1 <- df %>%
      droplevels() %>%
      left_join(., historical_ec_tomodel_centered_use[[src]], by = "location") %>%
      mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
    
    loo_indices <- df1 %>%
      group_by(location) %>%
      crossv_loo_grouped() %>%
      pull(train) %>%
      map("idx")
    
    ## Recursive feature addition
    rfa_out_df <- select_features_met(data = df1, env.col = "location", search.method = "hill.climb")
    
    
    historical_feature_selection_list$out[[i]] <- bind_rows(rfa_out_df)
    
}

historical_feature_selection <- unnest(historical_feature_selection_list, out)


## Save historical feature selection
save("historical_feature_selection", 
     "historical_fact_reg_feature_selection",
     file = file.path(result_dir, "historical_feature_selection_results.RData"))


## Reload the concurrent data and save everything together
load(file.path(result_dir, "concurrent_feature_selection_results.RData"))

save("concurrent_feature_selection", "historical_feature_selection", 
     "concurrent_fact_reg_feature_selection", "historical_fact_reg_feature_selection",
     file = file.path(result_dir, "feature_selection_results.RData"))


