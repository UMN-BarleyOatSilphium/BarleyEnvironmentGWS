## S2MET phenotypic modeling with covariates
## 
## Author: Jeff Neyhart
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

# data source
source_use <- "daymet"

# Number of cores
n_cores <- 8



# Load data ---------------------------------------------------------------

# Load covariates for environments and historical covariates
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


## Filter BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp,
         environment %in% train_test_env) %>% # only model train/test environments
  mutate_at(vars(line_name, environment), as.factor)


# Modeling ----------------------------------------------------------------

# Historical covariates ===================================================

## Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

# Use just the TP for modeling
S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp)



# Identify the most suitable time frame for summarization -----------------


# Implement the stepwise regression algorithm per timeframe to identify the
# timeframe with the lowest PRESS

historical_ec_tomodel_timeframe_centered_use <- historical_ec_tomodel_timeframe_centered %>%
  subset(., grepl(pattern = source_use, x = names(.)))

# Subset different windows of years
historical_ec_tomodel_window_centered_use <- historical_ec_tomodel_window_centered %>%
  subset(., grepl(pattern = source_use, x = names(.)))

# Combine the list
historical_ec_tomodel_centered_use <- c(historical_ec_tomodel_timeframe_centered_use, historical_ec_tomodel_window_centered_use)

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


# CReate a setup df
historical_timeframe_selection <- S2_MET_loc_BLUEs_tomodel %>%
  left_join(., bind_rows(historical_ec_tomodel_centered_use), by = "location") %>%
  group_by(trait, time_frame) %>%
  nest() %>%
  ungroup()


# Split and parallelize
historical_timeframe_selection_out <- historical_timeframe_selection %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE) %>%
  coreApply(X = ., FUN = function(core_df) {
    
    # Vector to store output
    out <- vector("list", length = nrow(core_df))
    
    # Iterate over rows of core_df
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      # Factorize
      df1 <- mutate_at(row$data[[1]], vars(line_name, location), ~fct_contr_sum(as.factor(.)))
      
      loo_indices <- df1 %>%
        group_by(location) %>%
        crossv_loo_grouped() %>%
        pull(train) %>%
        map("idx")
      
      # 1. Fit a base model
      base_fit <- lm(value ~ 1 + line_name, data = df1)
      covariates_use <- subset(trait_covariate_df, trait == row$trait, covariate, drop = TRUE)
      
      ## Recursive feature addition ##
      ## Main effect
      
      # 2. Define the scope
      scope <- list(lower = formula(base_fit), upper = reformulate(c("line_name", covariates_use), response = "value"))
      # Run rfa
      rfa_out <- rfa_loo(object = base_fit, data = df1, scope = scope, metric = "RMSE", index = loo_indices, env.col = "location")
      
      ## Interactions
      # 1. Fit a base model
      base_fit_int <- update(base_fit, formula = reformulate(rfa_out$optVariables, response = "value"))
      # 2. Define the scope
      scope <- list(lower = formula(base_fit_int), 
                    upper = reformulate(c(rfa_out$optVariables, paste0("line_name:", covariates_use)), response = "value"))
      # Run rfa
      rfa_out_int <- rfa_loo(object = base_fit_int, data = df1, scope = scope, metric = "RMSE", 
                             index = loo_indices, env.col = "location")
      
      ## Create a tibble
      out[[i]] <- tibble(
        feat_sel_type = "rfa_cv",
        model = c("model2", "model3"),
        adhoc = list(rfa_out, rfa_out_int),
      )
      
    }
    
    # Add out to core_df and return
    core_df %>% 
      mutate(out = out) %>% 
      unnest(out)
    
  }) %>% bind_rows()
    

# Save these results
save("historical_timeframe_selection_out", file = file.path(result_dir, "historical_covariate_timeframe_selection.RData"))






# ## Ideas
# ## 1. Calculate locations means and use that as the input (and validation)
# ## 2. Use all data and fit a model with random year, location:year, and g:year
# ## 3. Use all data and fit the g + location model
# 
# 
# # The timeframe to use for prediction, as suggested in the location analysis script
# time_frame_use <- "time_frame5_2010_2014"
# 
# 
# 
# 
# 
# ## Select the historical covariate data
# historical_ec_tomodel_centered_use <- historical_ec_tomodel_timeframe_centered %>%
#   subset(., map_lgl(names(.), ~str_detect(., time_frame_use)))
# names(historical_ec_tomodel_centered_use) <- str_remove(string = names(historical_ec_tomodel_centered_use), pattern = ".time_frame5_2010_2014")
# 
# 
# # Create a data.frame of covariates per trait
# trait_covariate_df <- historical_ec_tomodel_centered_use %>%
#   map(~names(.)[-1:-3]) %>%
#   reduce(union) %>%
#   crossing(trait = traits, covariate = .) %>%
#   filter(
#     covariate != "awc_range",
#     !(str_detect(covariate, "bulk_density")),
#     !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
#     !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
#   ) 
# 
# 
# 
# # Use the following covariates for all traits (except where timing would make
# # it inappropriate):
# # 
# # radiation during vegetative
# # tmin during flowering
# # water stress during flowering
# # grain fill water stress
# # grain fill tmax
# 
# 
# # # Load packages
# # invisible(lapply(X = pkgs, library, character.only = TRUE))
# 
# historical_fact_reg <- S2_MET_loc_BLUEs_tomodel %>%
#   crossing(., source = map_chr(historical_ec_tomodel_centered_use, ~unique(pull(., "source")))) %>%
#   group_by(trait, source) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# for (i in seq_len(nrow(historical_fact_reg))) {
# 
#     df <- historical_fact_reg$data[[i]]
#     df$trait <- historical_fact_reg$trait[i]
#     
#     src <- historical_fact_reg$source[i]
# 
#     # Factorize
#     df1 <- df %>%
#       droplevels() %>%
#       left_join(., historical_ec_tomodel_centered_use[[src]], by = "location") %>%
#       mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
# 
#     # Apriori
#     ## Add covariates - filter
#     covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
#     apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "apriori")
# 
#     # Ad hoc
#     ## Add covariates - filter
#     covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
#     adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")
# 
#     # Ad hoc - without soil
#     ## Add covariates - filter
#     covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T),
#                              covariate, drop = TRUE)
#     adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")
# 
# 
#     ## Return results
#     historical_fact_reg$out[[i]] <- tibble(model = c("base", "base_alt", "model4", "model5"),
#            apriori = apriori_out,
#            adhoc = adhoc_out,
#            adhoc_nosoil = adhoc_nosoil_out)
# 
# }
# 
# 
# historical_fact_reg1 <- unnest(historical_fact_reg, out)
# 
# 
# ## Assess via LOO cv
# historical_fact_reg_loo_test <- historical_fact_reg %>% 
#   mutate(loo_indices = map(data, ~crossv_loo_grouped(group_by(., location))) %>%
#            map("train") %>% map(., ~map(., "idx")),
#          out = map2(out, loo_indices, ~mutate(.x, loo_indices = list(.y)))) %>%
#   unnest(out) %>%
#   filter(str_detect(model, "base", negate = TRUE)) %>%
#   mutate(adhoc_cv = map2(adhoc, loo_indices, ~rapid_cv(object = .x, index = .y, return.predictions = TRUE)))
#   
# historical_fact_reg_loo_test %>%
#   mutate(out = map(adhoc_cv, "cv_metrics") %>% map(~as.data.frame(as.list(.)))) %>%
#   unnest(out)
#   
# 
# 
# ## Prepare results for saving
# historical_fact_reg_feature_selection <- historical_fact_reg1 %>%
#   mutate(feat_sel_type = "stepAIC", direction = "forward") %>%
#   filter(str_detect(model, "base", negate = TRUE)) %>%
#   mutate_at(vars(apriori, adhoc, adhoc_nosoil), ~map(., ~list(optVariables = attr(terms(.x), "term.labels"))))
# 
# 
# 
# 
# 
# 
# # Feature selection ############################################################
# 
# historical_feature_selection_list <- historical_fact_reg %>%
#   mutate(out = list(NULL))
# 
# for (i in seq_len(nrow(historical_feature_selection_list))) {
# 
#     df <- historical_feature_selection_list$data[[i]] %>% 
#       mutate(trait = historical_feature_selection_list$trait[i])
#     
#     src <- historical_fact_reg$source[i]
#     
#     
#     # Factorize
#     df1 <- df %>%
#       droplevels() %>%
#       left_join(., historical_ec_tomodel_centered_use[[src]], by = "location") %>%
#       mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
#     
#     loo_indices <- df1 %>%
#       group_by(location) %>%
#       crossv_loo_grouped() %>%
#       pull(train) %>%
#       map("idx")
#     
#     ## Recursive feature addition
#     rfa_out_df <- select_features_met(data = df1, env.col = "location", search.method = "hill.climb")
#     
#     
#     historical_feature_selection_list$out[[i]] <- bind_rows(rfa_out_df)
#     
# }
# 
# historical_feature_selection <- unnest(historical_feature_selection_list, out)
# 
# 
# ## Save historical feature selection
# save("historical_feature_selection", 
#      "historical_fact_reg_feature_selection",
#      file = file.path(result_dir, "historical_feature_selection_results.RData"))
# 
# 
# ## Reload the concurrent data and save everything together
# load(file.path(result_dir, "concurrent_feature_selection_results.RData"))
# 
# save("concurrent_feature_selection", "historical_feature_selection", 
#      "concurrent_fact_reg_feature_selection", "historical_fact_reg_feature_selection",
#      file = file.path(result_dir, "feature_selection_results.RData"))
# 
# 
