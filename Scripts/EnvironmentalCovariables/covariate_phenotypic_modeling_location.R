## S2MET phenotypic modeling with covariates
## 
## Author: Jeff Neyhart
## 
## This script will look at variable selection to determine optimal models that use
## covariates
## 

# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))

library(modelr)
library(broom)
library(lme4)
library(car)

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


## Ideas
## 1. Calculate locations means and use that as the input (and validation)
## 2. Use all data and fit a model with random year, location:year, and g:year
## 3. Use all data and fit the g + location model


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

# Vector of covariate names
covariate_names <- names(historical_ec_tomodel_timeframe_centered[[1]])[-1:-3]




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


# Create a setup df
historical_timeframe_selection <- S2_MET_loc_BLUEs_tomodel %>%
  left_join(., bind_rows(historical_ec_tomodel_centered_use), by = "location") %>%
  group_by(trait, time_frame) %>%
  nest() %>%
  ungroup() %>%
  # Output list
  mutate(out = list(NULL)) %>%
  # Adhoc and adhoc no soil
  crossing(., feat_sel_type = c("adhoc", "adhoc_nosoil"))


# Iterate over rows
for (i in seq_len(nrow(historical_timeframe_selection))) {
  
  row <- historical_timeframe_selection[i,]
  
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
  
  # Remove soil variables if called for
  if (row$feat_sel_type == "adhoc_nosoil") {
    covariates_use <- str_subset(string = covariates_use, pattern = "soil", negate = TRUE)
  }
  
  
  ## Recursive feature addition ##
  ## Main effect
  
  # 2. Define the scope
  scope <- list(lower = formula(base_fit), upper = reformulate(c("line_name", covariates_use), response = "value"))
  # Run rfa
  rfa_out <- try(rfa_loo(object = base_fit, data = df1, scope = scope, metric = "RMSE", 
                         index = loo_indices, env.col = "location"))
  
  # If try is error, set rfa_out and rfa_out_int both as null
  if (class(rfa_out) == "try-error") {
    rfa_out <- rfa_out_int <- NULL
    
    # else proceed
  } else {
    
    ## Interactions
    # 1. Fit a base model
    base_fit_int <- update(base_fit, formula = reformulate(rfa_out$optVariables, response = "value"))
    # 2. Define the scope
    scope <- list(lower = formula(base_fit_int), 
                  upper = reformulate(c(rfa_out$optVariables, paste0("line_name:", covariates_use)), response = "value"))
    # Run rfa
    rfa_out_int <- try(rfa_loo(object = base_fit_int, data = df1, scope = scope, metric = "RMSE", 
                               index = loo_indices, env.col = "location"))
    
    # Set rfa_out_int to null if error
    if (class(rfa_out_int) == "try-error") rfa_out_int <- NULL
    
  }
  
  ## Create a tibble
  historical_timeframe_selection$out[[i]] <- tibble(
    model = c("model4", "model5"),
    covariates = list(rfa_out, rfa_out_int),
  )
  
  # Notify
  cat("Stepwise selection for trait", row$trait, "with timeframe", row$time_frame, "complete.\n")
  
}

historical_timeframe_selection_out <- unnest(historical_timeframe_selection, out)


# Save these results
save("historical_timeframe_selection_out", file = file.path(result_dir, "historical_covariate_timeframe_selection.RData"))


# 
# 
# ## Analyze the results ##
# 
# # What timeframe returns the lowest RMSE per trait?
# historical_timeframe_selection_out1 <- historical_timeframe_selection_out %>%
#   filter(model == "model3") %>%
#   mutate(results = map(adhoc, "finalResults") %>% map(as.list) %>% map(as_tibble)) %>%
#   unnest(results) %>%
#   # Parse timeframe
#   mutate(time_frame_type = str_extract(time_frame, "time_frame|window"),
#          time_frame1 = str_remove(time_frame, "time_frame|window")) %>%
#   separate(time_frame1, c("length", "start_year", "end_year"), sep = "_") %>%
#   mutate_at(vars(length, contains("year")), parse_guess) %>%
#   # Round the RMSE
#   mutate_at(vars(RMSE, R2), ~round(., 5))
# 
# 
# ## First, display the results of the naive selection of time_frame5
# historical_timeframe_selection_out1 %>%
#   filter(time_frame == "time_frame5_2010_2014")
# 
# 
# ## Find the the best, shortest, most recent timeframe based on RMSE
# ## Timeframe
# historical_timeframe_sorted_timeframe <- historical_timeframe_selection_out1 %>%
#   filter(time_frame_type == "time_frame") %>%
#   # Round the 
#   split(.$trait) %>%
#   map(~arrange(., RMSE, end_year, length))
# 
# ## Window
# historical_timeframe_sorted_window <- historical_timeframe_selection_out1 %>%
#   filter(time_frame_type == "window") %>%
#   # Round the 
#   split(.$trait) %>%
#   map(~arrange(., RMSE, end_year, length))
# 
# 
# ## Time frames appear to be best (or at least comparable to windows); use these 
# historical_timeframe_best_timeframe <- historical_timeframe_sorted_timeframe %>%
#   map(~filter(., length < 15)) %>% # USe a shorter length
#   map_df(~slice(., 1))
# 
# 
# ## These will be the feature selection covariates
# historical_feature_selection <- historical_timeframe_selection_out %>%
#   inner_join(., select(historical_timeframe_best_timeframe, trait, time_frame)) %>%
#   select(-data) %>%
#   gather(selection_type, covariates, adhoc, adhoc_nosoil) %>% 
#   unite(feat_sel_type, feat_sel_type, selection_type, sep = "_") %>%
#   mutate(model = case_when(model == "model2" ~ "model4", model == "model3" ~ "model5"))
#   
# 
# ## Use this timeframe for the apriori and stepwiseAIC procedures ## 
# 
# 
# ## A priori covariates - quite minimal; each should have a citation 
# ## 
# ## HeadingDate - everything before flowering
# ## 
# ## PlantHeight - everything before grain fill
# ##  
# ## Grain yield, grain protein, and test weight
# ## - Everything for plant height plus...
# ## - Elevated temperature during grain fill (Passarella et al 2005)
# ## - Drought during grain fill - (Savin and Nicholas  1996)
# ## 
# ## 
# ## 
# 
# apriori_covariate_df <- trait_covariate_df %>%
#   filter(str_detect(covariate, "soil", negate = TRUE)) %>%
#   {bind_rows(
#     filter(., trait %in% c("HeadingDate", "PlantHeight")),
#     filter(., trait %in% c("GrainProtein", "GrainYield", "TestWeight"),
#            covariate %in% c(subset(., trait == "PlantHeight", covariate, drop = TRUE), 
#                             "grain_fill.tmean_mean", "grain_fill.water_balance_sum"))
#   )}
# 
# 
# historical_fact_reg <- S2_MET_loc_BLUEs_tomodel %>%
#   crossing(., source = map_chr(historical_ec_tomodel_timeframe_centered, ~unique(pull(., "source")))) %>%
#   left_join(., distinct(historical_feature_selection, trait, time_frame)) %>%
#   group_by(trait, source, time_frame) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# for (i in seq_len(nrow(historical_fact_reg))) {
# 
#     df <- historical_fact_reg$data[[i]]
#     df$trait <- historical_fact_reg$trait[i]
# 
#     # data source
#     src <- historical_fact_reg$source[i]
#     # time frame
#     tf <- historical_fact_reg$time_frame[i]
#     
#     # Extract the covariates data to use
#     historical_ec_tomodel_centered_use <- bind_rows(historical_ec_tomodel_timeframe_centered) %>%
#       filter(source == src, time_frame == tf)
# 
#     # Factorize
#     df1 <- df %>%
#       droplevels() %>%
#       left_join(., historical_ec_tomodel_centered_use, by = "location") %>%
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
# ## Prepare results for saving
# historical_fact_reg_feature_selection <- historical_fact_reg1 %>%
#   mutate(feat_sel_type = "stepAIC", direction = "forward") %>%
#   filter(str_detect(model, "base", negate = TRUE)) %>%
#   mutate_at(vars(apriori, adhoc, adhoc_nosoil), ~map(., ~list(optVariables = attr(terms(.x), "term.labels")))) %>%
#   gather(selection_type, covariates, apriori, adhoc, adhoc_nosoil) %>% 
#   unite(feat_sel_type, feat_sel_type, selection_type, sep = "_") %>% 
#   mutate(feat_sel_type = ifelse(str_detect(feat_sel_type, "apriori"), "apriori", feat_sel_type))
# 
# 
# 
# 
# ## Create a data.frame of all covariates
# historical_all_features <- trait_covariate_df %>%
#   group_by(trait) %>% 
#   nest(.key = "covariates") %>% 
#   mutate(covariates = map(covariates, "covariate")) %>%
#   crossing(., source = names(ec_tomodel_centered), feat_sel_type = "all", model = c("model4", "model5")) %>%
#   mutate(covariates = modify_if(covariates, model == "model5", ~c(., paste0("line_name:", .))),
#          covariates = map(covariates, ~list(optVariables = .)))
# 
# 
# 
# 
# ## Save historical feature selection
# save("historical_feature_selection", "historical_fact_reg_feature_selection", "historical_all_features",
#      file = file.path(result_dir, "historical_feature_selection_results.RData"))
# 
# 
# ## Reload the concurrent data and save everything together
# load(file.path(result_dir, "concurrent_feature_selection_results.RData"))
# 
# save("concurrent_feature_selection", "historical_feature_selection",
#      "concurrent_fact_reg_feature_selection", "historical_fact_reg_feature_selection",
#      "concurrent_all_features", "historical_all_features",
#      file = file.path(result_dir, "feature_selection_results.RData"))
# 
# 
