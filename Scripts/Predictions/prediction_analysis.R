## S2MET Crop Model Prediction
## 
## Script for analyzing cross-validation results
## 
## Author: Jeff Neyhart
## Last modified: 5 March 2020
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# # Load other packages
library(broom)
library(car)
library(lmerTest)
library(patchwork)
library(cowplot)
library(modelr)
library(paletteer)


## Load the validation results
file_list <- list.files(result_dir, pattern = "reml_predictions.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))

# Load feature selections
load(file.path(result_dir, "feature_selection_results.RData"))



## Significant level
alpha <- 0.05


# 
# Modify prediction objects
# 

loeo_predictions_out <- bind_rows(loeo_predictions_out) %>%
  select(trait, site, .id, out) %>%
  unnest(out, names_repair = tidyr_legacy) %>%
  select(trait, model, feature_selection, prediction) %>%
  unnest(prediction) %>%
  rename(predicted_value = pred_complete)
  
env_external_predictions_out <- env_external_predictions_out %>%
  select(trait, out) %>%
  unnest(out, names_repair = tidyr_legacy) %>%
  select(trait, model, feature_selection, prediction) %>%
  unnest(prediction) %>%
  rename(predicted_value = pred_complete)


lolo_predictions_out <- bind_rows(lolo_predictions_out) %>%
  select(trait, site, .id, out) %>%
  unnest(out, names_repair = tidyr_legacy) %>%
  select(trait, model, time_frame, feature_selection, prediction) %>%
  unnest(prediction) %>%
  rename(predicted_value = pred_complete)

loc_external_predictions_out <- loc_external_predictions_out %>%
  select(trait, out) %>%
  unnest(out, names_repair = tidyr_legacy) %>%
  select(trait, model, time_frame, feature_selection, prediction) %>%
  unnest(prediction) %>%
  rename(predicted_value = pred_complete)


## Grab the prediction outputs and combine
prediction_list <- object_list %>%
  set_names(., .) %>%
  map(get) %>%
  subset(., map_lgl(., ~nrow(.) > 1))




# Accuracy across environments --------------------------------------------


## Combine data.frames and mutate columns
predictions_df_temp <- prediction_list %>%
  setNames(., gsub(pattern = "_predictions.*", replacement = "", x = names(.))) %>%
  imap(~mutate(.x, type = .y) ) %>%
  modify_if(., str_detect(names(.), "lo[a-z]o"), ~mutate(.x, .id = as.character("01"))) %>%
  # Combine time_frame and feature selection for LOLO longterm
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~mutate_if(., is.character, parse_guess) %>%
           mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  rename(selection = feature_selection, cv_rep = .id) %>%
  mutate(selection = ifelse(is.na(selection), "none", selection),
         pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  select(-which(names(.) %in% c(".id", "core", "trait1"))) %>%
  # Coalesce columns
  mutate(test_group = site) %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1", "nSite", "site1", "source"))) %>%
  # Subset trait-site combinations that are relevant
  inner_join(distinct(select(gather(distinct(S2_MET_BLUEs, trait, environment, location), var, site, -trait), -var)), .) %>%
  # Sort
  arrange(type, model, selection, cv_rep, site)


# Add information about time frame selection
time_frame_selection_info <- historical_feature_selection %>%
  select(trait, time_frame, time_frame_selection = selection, feat_sel_type) %>%
  distinct() %>%
  full_join(., map_df(list(historical_all_features, historical_apriori_feature_selection, historical_feature_selection), 
                      ~distinct_at(., vars(trait, time_frame, feat_sel_type, feat_sel_type)))) %>%
  bind_rows(., distinct_at(historical_feature_importance, vars(trait, time_frame, time_frame_selection = selection, feat_sel_type))) %>%
  rename(selection = feat_sel_type)

predictions_df <- predictions_df_temp %>%
  left_join(., time_frame_selection_info) %>%
  # Remove columns that are lists
  select_if(~!is.list(.)) %>%
  select(trait, site, model, type, selection, time_frame, time_frame_selection, test_group, cv_rep,
         pop, line_name, value, predicted_value)




## Calculate accuracy and bias per train group, model, and population
predictive_ability <- predictions_df %>%
  group_by(trait, model, pop, type, cv_rep, selection, time_frame, time_frame_selection, test_group) %>%
  # First calculate accuracy per environment within a cv rep
  mutate(ability = cor(predicted_value, value), 
         bias = bias(obs = value, pred = predicted_value),
         rmse = sqrt(mean((value - predicted_value)^2))) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, cv_rep, selection, time_frame, time_frame_selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(predicted_value, value), 
         bias_all = bias(obs = value, pred = predicted_value),
         rmse_all = sqrt(mean((value - predicted_value)^2))) %>%
  ungroup() %>%
  # Reduce observations
  distinct_at(vars(-line_name, -value, -contains("pred")))




## Adjust ability using heritability
within_environment_prediction_accuracy <- predictive_ability %>%
  filter(str_detect(type, "loeo|env_external")) %>%
  select(-contains("_all")) %>%
  left_join(., env_trait_herit, by = c("trait", "test_group" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability))





# Calculate predictive ability over all observations; also calculate RMSE
accuracy_bias_all <- predictive_ability %>% 
  distinct(trait, cv_rep, model, pop, type, selection, time_frame, time_frame_selection, ability_all, bias_all)

# First create annotation df
across_site_prediction_accuracy_annotation <- accuracy_bias_all %>%
  mutate(bias_all = bias_all * 100) %>%
  mutate_at(vars(ability_all, bias_all), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_ann = paste0("r[MP]==", ability_all),
         bias_ann = paste0("Bias==", bias_all, "*'%'"),
         annotation = paste0(ability_ann, "*','~", bias_ann))




## Save everything

output_filename <- paste0("prediction_accuracy_compiled_reml_", format(Sys.Date(), "%Y%m%d"), ".RData")

save("predictions_df", "predictive_ability", "across_site_prediction_accuracy_annotation",
     "within_environment_prediction_accuracy", 
     file = file.path(result_dir, output_filename))
