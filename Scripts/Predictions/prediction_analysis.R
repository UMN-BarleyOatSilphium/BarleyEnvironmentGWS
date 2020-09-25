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
# file_list <- list.files(result_dir, pattern = "predictions.RData$", full.names = TRUE)
file_list <- list.files(result_dir, pattern = "fact_reg.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))


## Significant level
alpha <- 0.05




## Grab the prediction outputs and combine
prediction_list <- object_list %>%
  set_names(., .) %>%
  map(get) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_remove(names(.), "_predictions_out"))



# Accuracy across environments --------------------------------------------


## Combine data.frames and mutate columns
predictions_df <- prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  modify_if(., str_detect(names(.), "lo[a-z]o"), ~mutate(.x, .id = as.character("01"))) %>%
  # Combine time_frame and feature selection for LOLO longterm
  modify_if(.x = ., .p = map_lgl(., ~!all(is.na(.$time_frame))), ~unite(.x, feature_selection, feature_selection, time_frame, sep = ":")) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~mutate_if(., is.character, parse_guess) %>%
           mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  rename(selection = feature_selection, cv_rep = .id) %>%
  mutate(selection = ifelse(is.na(selection), "none", selection),
         pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  select(-which(names(.) %in% c(".id", "core", "trait1"))) %>%
  # Coalesce columns
  mutate(test_group = site, nGroup = nSite) %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1", "nSite", "site1", "source"))) %>%
  # Subset trait-site combinations that are relevant
  inner_join(distinct(select(gather(distinct(S2_MET_BLUEs, trait, environment, location), var, site, -trait), -var)), .) %>%
  # Sort
  arrange(type, model, selection, cv_rep, site)



## Calculate accuracy and bias per train group, model, and population
predictive_ability <- predictions_df %>%
  group_by(trait, model, pop, type, cv_rep, test_group, selection) %>%
  # First calculate accuracy per environment within a cv rep
  mutate(ability = cor(pred_complete, value), 
         bias = bias(obs = value, pred = pred_complete),
         rmse = sqrt(mean((value - pred_complete)^2))) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection, cv_rep) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = bias(obs = value, pred = pred_complete),
         rmse_all = sqrt(mean((value - pred_complete)^2))) %>%
  ungroup() %>%
  # Reduce observations
  distinct_at(vars(-line_name, -value, -contains("pred")))




## Adjust ability using heritability
within_environment_prediction_accuracy <- predictive_ability %>%
  filter(type %in% c("loeo", "env_external", "cv_env")) %>%
  select(-contains("_all")) %>%
  left_join(., env_trait_herit, by = c("trait", "test_group" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability))





## Quick plot of accuracy across all environments/locations
predictive_ability %>%
  filter(selection != "none") %>%
  distinct(trait, model, pop, type, selection, ability_all, rmse_all) %>%
  ggplot(aes(x = model, y = ability_all, fill = selection)) +
  geom_col(position = position_dodge(0.9)) +
  facet_grid(trait ~ type + pop)


# Plot predicted versus observed values for a subset
predictions_df %>%
  filter(str_detect(type, "external"), model == "model3_cov", str_detect(selection, "lasso|stepwise")) %>%
  ggplot(aes(x = pred_complete, y = value, color = test_group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 0.5) +
  scale_color_discrete(guide = FALSE) +
  facet_wrap(~ type + trait + selection + model + pop, scales = "free", ncol = 4, labeller = labeller(.multi_line = FALSE)) +
  theme_presentation2(10)


predictive_ability %>%
  filter(model == "model3_cov", pop == "tp", selection == "rfa_cv_adhoc") %>%
  distinct(trait, model, type, selection, ability_all)




# Calculate predictive ability over all observations; also calculate RMSE
accuracy_bias_all <- predictive_ability %>% 
  distinct(trait, cv_rep, model, pop, type, selection, ability_all, bias_all)

# First create annotation df
across_site_prediction_accuracy_annotation <- accuracy_bias_all %>%
  mutate(bias_all = bias_all * 100) %>%
  mutate_at(vars(ability_all, bias_all), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_ann = paste0("r[MP]==", ability_all),
         bias_ann = paste0("Bias==", bias_all, "*'%'"),
         annotation = paste0(ability_ann, "*','~", bias_ann))



## Quick plot of LOLO long-term covariates

lolo_longterm_plotlist <- predictions_df %>%
  filter(type == "lolo_longterm_covariates") %>%
  split(.$trait) %>%
  map(~{
    ann_df <- left_join(.x, across_site_prediction_accuracy_annotation) %>%
      distinct(trait, model, selection, pop, ability_ann)
    
    ggplot(data = .x, aes(x = pred_complete, y = value, color = test_group)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = ability_ann), size = 2, inherit.aes = FALSE, 
                vjust = -1, parse = TRUE, hjust = 1.1) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(selection ~ model + pop) +
      theme_presentation2(10)
  })







## Save everything
save("predictions_df", "predictive_ability", "across_site_prediction_accuracy_annotation",
     "within_environment_prediction_accuracy", 
     file = file.path(result_dir, "prediction_accuracy_compiled.RData"))

