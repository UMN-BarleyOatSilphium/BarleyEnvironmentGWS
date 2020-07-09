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


# model1: y = G + r
# model2: y = G + E + r
# model3: y = G + E + GE + r
# model4: y = G + L + r
# model5: y = G + L + GL + r




## Grab the prediction outputs and combine
prediction_list <- map(set_names(object_list, object_list), get) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_remove(names(.), "_predictions_out"))



# Accuracy across environments --------------------------------------------


## Combine data.frames and mutate columns
predictions_df <- prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~mutate_if(., is.character, parse_guess) %>%
           mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  rename(selection = feature_selection) %>%
  mutate(selection = ifelse(is.na(selection), "none", selection),
         pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  select(-which(names(.) %in% c(".id", "core", "trait1"))) %>%
  # Coalesce columns
  mutate(leave_one_group = site, nGroup = nSite) %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1", "nSite", "site1", "source"))) %>%
  # Subset trait-site combinations that are relevant
  inner_join(distinct(select(gather(distinct(S2_MET_BLUEs, trait, environment, location), var, site, -trait), -var)), .)



## Calculate accuracy and bias per train group, model, and population
predictive_ability <- predictions_df %>%
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = bias(obs = value, pred = pred_complete),
         rmse = sqrt(mean((value - pred_complete)^2))) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = bias(obs = value, pred = pred_complete),
         rmse_all = sqrt(mean((value - pred_complete)^2))) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  summarize_at(vars(ability, bias, rmse, ability_all, bias_all, rmse_all, nObs, nGroup), mean) %>%
  ungroup()


## Drop crookston as a location and recalculate location prediction accuracy for grain yield
predictions_df %>%
  filter(type == "lolo", leave_one_group != "Crookston", trait == "GrainYield") %>%
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  summarize(ability_all = cor(pred_complete, value), 
            bias_all = bias(obs = value, pred = pred_complete),
            rmse_all = sqrt(mean((value - pred_complete)^2))) %>%
  as.data.frame()
  


## Adjust ability using heritability
within_environment_prediction_accuracy <- predictive_ability %>%
  filter(type %in% c("loeo", "env_external")) %>%
  select(-contains("_all")) %>%
  left_join(., env_trait_herit, by = c("trait","leave_one_group" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability))





## Quick plot of accuracy across all data points
predictive_ability %>%
  filter(type %in% c("lolo", "loeo"), selection != "stepAIC_adhoc") %>%
  distinct(trait, model, pop, type, selection, ability_all, rmse_all) %>%
  ggplot(aes(x = model, y = ability_all, fill = selection)) +
  geom_col(position = position_dodge(0.9)) +
  facet_grid(trait ~ type + pop)

## Print some predictive abilities
predictive_ability %>%
  filter(type == "lolo", str_detect(model, "model3"), pop == "tp") %>%
  select(trait, model, selection, contains("_all")) %>%
  distinct() %>%
  as.data.frame()


predictions_df %>%
  filter(type == "lolo", trait == "GrainProtein", !selection %in% c("none", "concurrent_rfa_cv_adhoc")) %>%
  ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 0.5) +
  scale_color_discrete(guide = FALSE) +
  facet_grid(type + selection ~ model + pop) +
  theme_presentation2(10)



predictive_ability %>%
  filter(model == "model3_cov", pop == "tp", selection == "rfa_cv_adhoc") %>%
  distinct(trait, model, type, selection, ability_all)




# Calculate predictive ability over all observations; also calculate RMSE
accuracy_bias_all <- predictive_ability %>% 
  distinct(trait, model, pop, type, selection, ability_all, bias_all)

# First create annotation df
across_site_prediction_accuracy_annotation <- accuracy_bias_all %>%
  mutate(bias_all = bias_all * 100) %>%
  mutate_at(vars(ability_all, bias_all), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(annotation = paste0("r[MP]==", ability_all, "*','~Bias==", bias_all, "%"))





## Save everything
save("predictions_df", "predictive_ability", "across_site_prediction_accuracy_annotation",
     "within_environment_prediction_accuracy", 
     file = file.path(result_dir, "prediction_accuracy_compiled.RData"))












# #### Appendix ####
# 
# 
# 
# # Plot predicted versus observed value
# g_loo_prediction_list <- loo_predictions_df %>%
#   filter(model %in% names(model_replace)) %>%
#   filter(selection != "concurrent_rfa_cv_adhoc") %>%
#   # Max character length of units
#   mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
#   group_by(trait, type) %>%
#   do({
#     df <- .
#     
#     # Convert yield to t ha^-1
#     if (unique(df$trait) == "GrainYield") {
#       df1 <- mutate_at(df, vars(value, pred_complete), ~. / 1000)
#     } else {
#       df1 <- df
#     }
#     
#     breaks <- map(subset(df1, pred_complete > 0, c(pred_complete, value)), pretty) 
#     
#     
#     ## Extract the appropriate accuracy annotation
#     r_mp_annotation <- left_join(distinct(df1, trait, type, model, pop, selection), loo_prediction_accuracy_annotation,
#                                  by = c("trait", "type", "model", "pop", "selection")) %>%
#       mutate(annotation = ability_all_annotation,
#              x = breaks$pred_complete[1], y = c(last(breaks$value)))
#     
#     df1_none <- filter(df1, selection == "none")
#     
#     ## First plot selection == none (identity covariance for environments and gxe)
#     g_plot_none <- df1_none %>%
#       ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
#       geom_abline(slope = 1, intercept = 0) +
#       geom_point(size = 0.5, alpha = 0.5) +
#       geom_text(data = left_join(distinct(df1_none, trait, type, model, selection), r_mp_annotation), 
#                 aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE, hjust = 0, size = 2) +
#       scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
#       scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
#       scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
#       facet_grid(selection + pop ~ model, switch = "y", 
#                  labeller = labeller(selection = f_ec_selection_replace, pop = f_pop_replace, model = f_model_replace)) +
#       labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
#       theme_presentation2(10)
#     
#     # Filter the df
#     df1_rest <- filter(df1, selection != "none")
#     
#     ## Create the plot of each feature selection procedure
#     g_plot <- df1_rest %>%
#       ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
#       geom_abline(slope = 1, intercept = 0) +
#       geom_point(size = 0.5, alpha = 0.5) +
#       geom_text(data = left_join(distinct(df1_rest, trait, type, model, selection), r_mp_annotation), 
#                 aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE, hjust = 0, size = 2) +
#       scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
#       scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
#       scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
#       facet_grid(selection + pop ~ model, switch = "y", 
#                  labeller = labeller(selection = f_ec_selection_replace, pop = f_pop_replace, model = f_model_replace)) +
#       labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
#       theme_presentation2(10)
#     
#     ## Output a list of plots
#     tibble(base_model_plot = list(g_plot_none),
#            full_covariate_selection_plot = list(g_plot))
#     
#   }) %>% ungroup()
# 
# 
# 
# ## Combine base models and all covariate selection
# for (i in seq_len(nrow(g_loo_prediction_list))) {
#   
#   # Combine the plots
#   plot_comb <- plot_grid(g_loo_prediction_list$base_model_plot[[i]], g_loo_prediction_list$full_covariate_selection_plot[[i]],
#                          rel_heights = c(1/3, 1), ncol = 1)
#   
#   # Save
#   filename <- paste0("loo_model_predictions_observations_", g_loo_prediction_list$trait[i], 
#                      "_", g_loo_prediction_list$type[i], ".jpg")
#   
#   ggsave(filename = filename, plot = plot_comb, path = fig_dir, width = 6, height = 16, dpi = 1000)
#   
# }
# 
# 
# ## Just plot the reduced set of feature selection
# for (i in seq_len(nrow(g_loo_prediction_list))) {
#   
#   # Combine the plots
#   plot_comb <- g_loo_prediction_list$full_covariate_selection_plot[[i]] %>%
#     modify_at("data", ~filter(., selection != "stepAIC_adhoc"))
#   plot_comb$layers[[3]]$data <- filter(plot_comb$layers[[3]]$data, selection != "stepAIC_adhoc")
#   
#   # Save
#   filename <- paste0("loo_covariate_models_predictions_observations_", g_loo_prediction_list$trait[i], 
#                      "_", g_loo_prediction_list$type[i], ".jpg")
#   
#   ggsave(filename = filename, plot = plot_comb, path = fig_dir, width = 6, height = 10, dpi = 1000)
#   
# }
# 
# 
# 
# 
# ## Visualize
# # Plot mean and range of predictions
# g_loo_predictions_summ <- loo_prediction_accuracy_summ1 %>%
#   filter(model %in% names(model_replace)) %>%
#   mutate(selection = f_ec_selection_replace(selection),
#          selection = fct_relevel(selection, "None")) %>% # Put "None" first
#   ggplot(aes(x = selection, group = model, color = model)) +
#   # geom_linerange(aes(ymin = ability_lower, ymax = ability_upper, group = model), 
#   #                position = position_dodge(0.9), color = "grey85") +
#   geom_point(aes(y = ability), position = position_dodge(0.9)) +
#   geom_text(data = annotation_df, aes(x = 1, y = 1, label = annotation), size = 2, hjust = 0,
#             inherit.aes = FALSE) +
#   geom_segment(data = annotation_df, aes(x = 3, xend = 3, y = 1, yend = 1 - LSD), inherit.aes = FALSE) +
#   scale_color_paletteer_d(package = "dutchmasters", palette = "milkmaid",
#                           name = "Model", labels = f_model_replace) +
#   scale_y_continuous(name = "Predictive ability", breaks = pretty) +
#   # scale_x_discrete(name = "Validation scheme", labels = f_pop_replace) +
#   scale_x_discrete(name = "Covariate selection", labels = function(x) parse(text = x)) +
#   facet_grid(type ~ trait + pop, labeller = labeller(trait = str_add_space, type = toupper, pop = f_pop_replace),
#              switch = "y") +
#   theme_presentation2(10) +
#   theme(strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Save
# ggsave(filename = "loo_model_predictions_summary.jpg", plot = g_loo_predictions_summ,
#        path = fig_dir, width = 12, height = 6, dpi = 1000)