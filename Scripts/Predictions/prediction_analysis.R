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
file_list <- list.files(result_dir, pattern = "predictions.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))

# Load feature selections
load(file.path(result_dir, "feature_selection_results.RData"))



## Significant level
alpha <- 0.05


# 
# Modify prediction objects
# 

loeo_varcomp_predictions_out <- loeo_varcomp_predictions_out %>%
  unnest(predictions, names_repair = tidyr_legacy) %>%
  select(model, feature_selection, out) %>%
  unnest(out)

loeo_varcomp_interval_predictions_out <- loeo_varcomp_interval_predictions_out %>%
  unnest(predictions, names_repair = tidyr_legacy) %>%
  select(model, feature_selection, out) %>%
  unnest(out)

env_external_varcomp_predictions_out <- env_external_varcomp_predictions_out %>%
  unnest(predictions, names_repair = tidyr_legacy) %>%
  select(model, feature_selection, out) %>%
  unnest(out)

env_external_interval_varcomp_predictions_out <- env_external_interval_varcomp_predictions_out %>%
  unnest(predictions, names_repair = tidyr_legacy) %>%
  select(model, feature_selection, out) %>%
  unnest(out)



lolo_predictions_out <- bind_rows(lolo_predictions_out) %>%
  select(trait, site, .id, out) %>%
  unnest(out, names_repair = tidyr_legacy) %>%
  select(trait, model, time_frame, feature_selection, prediction) %>%
  unnest(prediction)

loc_external_predictions_out <- loc_external_predictions_out %>%
  select(trait, out) %>%
  unnest(out, names_repair = tidyr_legacy) %>%
  select(trait, model, time_frame, feature_selection, prediction) %>%
  unnest(prediction)


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
  select_if(~!is.list(.))





## Calculate accuracy and bias per train group, model, and population
predictive_ability <- predictions_df %>%
  group_by(trait, model, pop, type, cv_rep, selection, time_frame, time_frame_selection, test_group) %>%
  # First calculate accuracy per environment within a cv rep
  mutate(ability = cor(pred_complete, value), 
         bias = bias(obs = value, pred = pred_complete),
         rmse = sqrt(mean((value - pred_complete)^2))) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, cv_rep, selection, time_frame, time_frame_selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = bias(obs = value, pred = pred_complete),
         rmse_all = sqrt(mean((value - pred_complete)^2))) %>%
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
save("predictions_df", "predictive_ability", "across_site_prediction_accuracy_annotation",
     "within_environment_prediction_accuracy", 
     file = file.path(result_dir, "prediction_accuracy_compiled.RData"))











## Quick plot of accuracy across all environments/locations
## 
## The plot below is specifically designed for LOLO
## 
# predictive_ability %>%
#   filter(type == "loc_external", selection != "none", model == "model3_cov") %>%
#   distinct(trait, model, pop, type, selection, time_frame, time_frame_selection, ability_all, rmse_all) %>%
#   unite(selection1, selection, time_frame, sep = ":", remove = FALSE) %>%
#   ggplot(aes(x = time_frame_selection, y = ability_all, fill = selection, group = selection1)) +
#   geom_col(position = position_dodge(0.9)) +
#   facet_grid(pop + type ~ trait, labeller = labeller(.multi_line = FALSE)) +
#   theme(legend.position = "bottom")

# Distinct predictive abilities per trait, pop, selection, type, model
predictive_ability_distinct <- predictive_ability %>%
  distinct(trait, model, selection, type, pop, ability_all, rmse_all) %>%
  filter(selection != "none", str_detect(selection, "nosoil|lasso|apriori", negate = TRUE))


predictive_ability_distinct_loeo <- predictive_ability_distinct %>%
  filter(str_detect(type, "loeo")) %>%
  mutate(annotation = paste0("rMP=", format_numbers(ability_all)))

# Create a big plot of predicted vs observed values
g_big_predictions <- predictions_df %>%
  inner_join(., predictive_ability_distinct_loeo) %>%
  ggplot(aes(x = pred_complete, y = value, color = test_group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 0.5) +
  geom_text(data = predictive_ability_distinct_loeo, aes(x = Inf, y = -Inf, label = annotation), 
            vjust = -1, hjust = 1, size = 2, inherit.aes = FALSE) +
  scale_color_discrete(guide = FALSE) +
  facet_wrap(~ trait + pop + type + selection + model, scales = "free", nrow = 10, labeller = labeller(.multi_line = FALSE)) +
  theme_presentation2(6)

# Save
ggsave(filename = "loeo_compare_covariate_predictions.png", plot = g_big_predictions, path = fig_dir,
       height = 20, width = 20)



predictive_ability_distinct_env_external <- predictive_ability_distinct %>%
  filter(str_detect(type, "env_external")) %>%
  mutate(annotation = paste0("rMP=", format_numbers(ability_all)))

# Create a big plot of predicted vs observed values
g_big_predictions <- predictions_df %>%
  inner_join(., predictive_ability_distinct_env_external) %>%
  ggplot(aes(x = pred_complete, y = value, color = test_group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 0.5) +
  geom_text(data = predictive_ability_distinct_env_external, aes(x = Inf, y = -Inf, label = annotation), 
            vjust = -1, hjust = 1, size = 2, inherit.aes = FALSE) +
  scale_color_discrete(guide = FALSE) +
  facet_wrap(~ trait + pop + type + selection + model, scales = "free", nrow = 10, labeller = labeller(.multi_line = FALSE)) +
  theme_presentation2(6)

# Save
ggsave(filename = "env_external_compare_covariate_predictions.png", plot = g_big_predictions, path = fig_dir,
       height = 20, width = 8)






# Compare prediction accuracies within environments


# # Plot predicted versus observed values for a subset
# predictions_df %>%
#   filter(str_detect(type, "external"), model == "model3_cov", str_detect(selection, "lasso|stepwise")) %>%
#   ggplot(aes(x = pred_complete, y = value, color = test_group)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point(size = 0.5) +
#   scale_color_discrete(guide = FALSE) +
#   facet_wrap(~ type + trait + selection + model + pop, scales = "free", ncol = 4, labeller = labeller(.multi_line = FALSE)) +
#   theme_presentation2(10)
# 
# 
# predictive_ability %>%
#   filter(model == "model3_cov", pop == "tp", selection == "rfa_cv_adhoc") %>%
#   distinct(trait, model, type, selection, ability_all)




# ## Quick plot of LOLO long-term covariates
# 
# lolo_longterm_plotlist <- predictions_df %>%
#   filter(type == "lolo_longterm_covariates") %>%
#   split(.$trait) %>%
#   map(~{
#     ann_df <- left_join(.x, across_site_prediction_accuracy_annotation) %>%
#       distinct(trait, model, selection, pop, ability_ann)
#     
#     ggplot(data = .x, aes(x = pred_complete, y = value, color = test_group)) +
#       geom_abline(slope = 1, intercept = 0) +
#       geom_point(size = 0.5) +
#       geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = ability_ann), size = 2, inherit.aes = FALSE, 
#                 vjust = -1, parse = TRUE, hjust = 1.1) +
#       scale_color_discrete(guide = FALSE) +
#       facet_grid(selection ~ model + pop) +
#       theme_presentation2(10)
#   })














