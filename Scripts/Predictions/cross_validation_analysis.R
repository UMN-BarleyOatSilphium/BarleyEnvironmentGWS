## S2MET Predictions
## 
## Script for analyzing cross-validation results
## 
## Author: Jeff Neyhart
## Last modified: January 28, 2019
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# # Load some packages
library(lubridate)
library(effects) # For ls means / marginal means
library(ggforce)
library(ggridges)
library(gridExtra)
library(broom)
library(car)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))

# Load results
cv_files <- list.files(result_dir, full.names = TRUE, pattern = "cross_validation_")
lapply(cv_files, load, envir = .GlobalEnv)


## Significant level
alpha <- 0.05


## Create a factor for the distance methods
dist_method_replace <- c("pheno_dist" = "Phenotypic Distance", "pheno_loc_dist" = "Location Phenotypic Distance",  "great_circle_dist" = "Great Circle Distance", 
                         "OYEC_All" = "One Year All ECs", "OYEC_Mean" = "One Year Mean Cor EC", "OYEC_IPCA" = "One Year IPCA Cor EC", 
                         "MYEC_All" = "Multi Year All ECs", "MYEC_Mean" = "Multi Year Mean Cor EC", "MYEC_IPCA" = "Multi Year IPCA Cor EC",
                         "sample" = "Random")
dist_method_abbr <- abbreviate(dist_method_replace)

## Alternative abbreviations
dist_method_abbr <- setNames(c("PD", "LocPD", "GCD", "1Yr-All-EC", "1Yr-Mean-EC", "1Yr-IPCA-EC", "All-EC", "Mean-EC", "IPCA-EC", "Random"),
                             names(dist_method_replace))

colors <- umn_palette(3)
colors_use <- c(colors[2], "#FFDE7A", colors[c(1, 3, 8, 4, 9, 5, 10)], "grey75")
dist_colors <- setNames(colors_use, dist_method_abbr)

## A vector to rename models
model_replace <- c("M2" = "M1 (Main effect)", "M3" = "M2 (GxE)", "M4" = "M1 (Main effect)", "M5" = "M2 (GxE)")












###
### Cross-validation results
### 
### 

# Regular cross-validation



## CV1 and CV2 - predicting untested genotypes or partially tested genotypes
cv12_accuracy <- cv12_prediction %>% 
  unnest(prediction) %>% 
  separate(.id, c(".id", "fold"), sep = "_") %>%
  group_by(cv, trait, model, .id, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>% 
  mutate(zscore = ztrans(accuracy)) %>% 
  mutate_at(vars(model, environment), as.factor)


## Fit a model for each CV
cv12_tomodel <- cv12_accuracy %>%
  group_by(cv, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(cv12_tomodel))) {
  
  df <- unnest(cv12_tomodel[i,], data)
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment) + (1|model:environment), data = df)
  
  # Significance testing
  anova_df <- tidy(anova(fit))
  ranova_df <- tidy(ranova(fit))
  effects_df <- as.data.frame(Effect("model", fit))
  
  cv12_tomodel$out[[i]] <- data_frame(anova = list(anova_df), ranova = list(ranova_df), effects = list(effects_df))
  
}

cv12_summ <- cv12_tomodel %>% 
  unnest(out) %>% 
  select(-data) %>%
  unnest(effects) %>%
  mutate_at(vars(fit, lower, upper), zexp)


## Plot
cv12_summ %>%
  mutate(cv = toupper(cv)) %>%
  ggplot(aes(x = cv, y = fit, ymin = lower, ymax = upper, shape = model, color = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(~ trait)




## CV0 and CV00 - predicting future years
cv_zero_future_prediction %>% split(.$cv) %>% map(~group_by(., trait, model) %>% summarize(n = n())) 


# Calculate accuracy
cv_zero_future_acc <- cv_zero_future_prediction %>% 
  unnest(prediction) %>% 
  separate(.id, c(".id", "fold"), sep = "_") %>%
  group_by(trait, cv, model, .id, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(zscore = ztrans(accuracy)) %>%
  mutate_at(vars(model, environment), as.factor)

## Summarize for CV0
cv0_future_summ <- cv_zero_future_acc %>% 
  filter(cv == "cv0") %>% 
  group_by(cv, trait, model) %>% 
  summarize(fit = mean(accuracy), sd = sd(accuracy), n = n()) %>%
  ungroup()

# cv    trait       model   fit
# 1 cv0   GrainYield  M4    0.366
# 2 cv0   GrainYield  M5_PD 0.383
# 3 cv0   HeadingDate M4    0.834
# 4 cv0   HeadingDate M5_PD 0.798
# 5 cv0   PlantHeight M4    0.568
# 6 cv0   PlantHeight M5_PD 0.513


## Fit a model for CV00
cv00_future_tomodel <- cv_zero_future_acc %>%
  filter(cv == "cv00") %>%
  group_by(cv, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(cv00_future_tomodel))) {

  df <- unnest(cv00_future_tomodel[i,], data)
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment) + (1|model:environment), data = df)
  
  # Significance testing
  anova_df <- tidy(anova(fit))
  ranova_df <- tidy(ranova(fit))
  effects_df <- as.data.frame(Effect("model", fit))
  
  cv00_future_tomodel$out[[i]] <- data_frame(anova = list(anova_df), ranova = list(ranova_df), effects = list(effects_df))
  
}

cv00_future_summ <- cv00_future_tomodel %>% 
  unnest(out) %>% 
  select(-data) %>%
  unnest(effects) %>%
  mutate_at(vars(fit, lower, upper), zexp)


## Combine and plot
cv_future_combine <- bind_rows(cv0_future_summ, cv00_future_summ) %>%
  mutate(cv = paste0(cv, "_future"))

cv_future_combine %>% 
  mutate(cv = toupper(cv)) %>%
  ggplot(aes(x = cv, y = fit, ymin = lower, ymax = upper, shape = model, color = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(~ trait)





## CV0 and CV00 - leave one environment out
cases <- cv_zero_loeo_prediction %>% 
  unnest(prediction) %>%
  distinct(cv, trait, .id, model, environment) %>%
  split(.$cv) %>%
  map(~group_by(., trait, model, environment) %>% summarize(n = n())) 


# Calculate accuracy
cv_zero_loeo_acc <- cv_zero_loeo_prediction %>% 
  unnest(prediction) %>% 
  group_by(trait, cv, model, .id, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(zscore = ztrans(accuracy)) %>%
  mutate_at(vars(model, environment), as.factor)

## Summarize for CV0
cv0_loeo_summ <- cv_zero_loeo_acc %>% 
  filter(cv == "cv0") %>% 
  group_by(cv, trait, model) %>% 
  summarize(fit = mean(accuracy)) %>%
  ungroup()

# cv    trait       model   fit
# 1 cv0   GrainYield  M4    0.406
# 2 cv0   GrainYield  M5_PD 0.353
# 3 cv0   HeadingDate M4    0.832
# 4 cv0   HeadingDate M5_PD 0.803
# 5 cv0   PlantHeight M4    0.544
# 6 cv0   PlantHeight M5_PD 0.514


## Fit a model for CV00
cv00_loeo_tomodel <- cv_zero_loeo_acc %>%
  filter(cv == "cv00") %>%
  group_by(cv, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(cv00_loeo_tomodel))) {
  
  df <- unnest(cv00_loeo_tomodel[i,], data)
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment) + (1|model:environment), data = df)
  
  # Significance testing
  anova_df <- tidy(anova(fit))
  ranova_df <- tidy(ranova(fit))
  effects_df <- as.data.frame(Effect("model", fit))
  
  cv00_loeo_tomodel$out[[i]] <- data_frame(anova = list(anova_df), ranova = list(ranova_df), effects = list(effects_df))
  
}

cv00_loeo_summ <- cv00_loeo_tomodel %>% 
  unnest(out) %>% 
  select(-data) %>%
  unnest(effects) %>%
  mutate_at(vars(fit, lower, upper), zexp)


## Combine and plot
cv_loeo_combine <- bind_rows(cv0_loeo_summ, cv00_loeo_summ) %>%
  mutate(cv = paste0(cv, "_loeo"))

cv_loeo_combine %>%
  mutate(cv = toupper(cv)) %>%
  ggplot(aes(x = cv, y = fit, ymin = lower, ymax = upper, shape = model, color = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(~ trait)



cv_combine <- bind_rows(cv12_summ, cv_loeo_combine, cv_future_combine) %>% mutate(cv_class = "cv")













#####
##### Parent-offspring validation
##### 




pov1_summ <- pov1_prediction %>% 
  unnest() %>%
  group_by(trait, model, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize(fit = mean(accuracy)) %>%
  ungroup() %>%
  mutate(cv = "pov1")


## Summarize for POV0
pov0_future_summ <- pov0_future_prediction %>%
  bind_rows() %>%
  unnest(out) %>% unnest() %>%
  group_by(trait, model, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize(fit = mean(accuracy)) %>%
  ungroup() %>%
  mutate(cv = "pov0_future")



cases <- pov0_loeo_prediction %>%
  bind_rows() %>%
  unnest(out) %>% 
  group_by(., trait, model) %>% 
  summarize(n = n())



pov0_loeo_summ <- pov0_loeo_prediction %>%
  bind_rows() %>%
  unnest(out) %>% unnest(prediction) %>%
  group_by(trait, model, environment) %>%
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize(fit = mean(accuracy)) %>%
  ungroup() %>%
  mutate(cv = "pov0_loeo")


pov00_future_tomodel <- pov00_future_prediction %>%
  bind_rows() %>% 
  unnest(out) %>% unnest(prediction)

pov00_loeo_tomodel <- pov00_loeo_prediction %>%
  bind_rows() %>% 
  unnest(out) %>% unnest(prediction)


## Fit models to POV00
pov00_tomodel <- bind_rows(mutate(pov00_future_tomodel, cv = "pov00_future"), mutate(pov00_loeo_tomodel, cv = "pov00_loeo")) %>%
  group_by(cv, trait, model, environment) %>%
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(model = as.factor(model),
         zscore = ztrans(accuracy)) %>%
  group_by(cv, trait) %>%
  nest() %>%
  mutate(out = list(NULL))


for (i in seq(nrow(pov00_tomodel))) {
  
  df <- unnest(pov00_tomodel[i,], data)
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment), data = df)
  
  # Significance testing
  anova_df <- tidy(anova(fit))
  ranova_df <- tidy(ranova(fit))
  effects_df <- as.data.frame(Effect("model", fit))
  
  pov00_tomodel$out[[i]] <- data_frame(anova = list(anova_df), ranova = list(ranova_df), effects = list(effects_df))
  
}

pov00_summ <- pov00_tomodel %>% 
  unnest(out) %>% 
  select(-data) %>%
  unnest(effects) %>%
  mutate_at(vars(fit, lower, upper), zexp)


pov_combine <- bind_rows(pov1_summ, pov0_future_summ, pov0_loeo_summ, pov00_summ) %>%
  mutate(cv_class = "pov")



## Combine and plot
pov_combine %>%
  separate(cv, c("cv", "type"), sep = "_") %>%
  mutate(cv_ann = ifelse(!is.na(type), paste0(toupper(cv), "\n(", str_to_title(type), ")"), toupper(cv))) %>%
  ggplot(aes(x = cv_ann, y = fit, ymin = lower, ymax = upper, shape = model, color = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  
  
  



## 
## Combine all CV and POV results
## 

all_val_combine <- bind_rows(cv_combine, pov_combine) %>%
  separate(cv, c("cv", "type"), sep = "_") %>%
  mutate(cv_ann = ifelse(!is.na(type), paste0(toupper(cv), "\n(", str_to_title(type), ")"), toupper(cv)),
         cv_class = toupper(cv_class))


g_all_cv <- all_val_combine %>%
  ggplot(aes(x = cv_ann, y = fit, ymin = lower, ymax = upper, color = model, shape = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(trait ~ cv_class, scales = "free_x", space = "free_x") +
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  scale_color_discrete(name = "Model") +
  scale_shape_discrete(name = "Model") +
  theme_presentation2(base_size = 12) +
  theme(axis.title.x = element_blank(), legend.position = "bottom")


ggsave(filename = "all_cross_validation_accuracy.jpg", plot = g_all_cv, path = fig_dir, width = 6, height = 5, dpi = 1000)



#























################
################
## Appendix ####
################
################



# # Rename the distance methods
# cumulative_pred_results <- cluster_pred_out %>% 
#   unnest(out) %>%
#   rename(dist_method = model) %>%
#   mutate(iter = parse_number(dist_method),
#          iter = ifelse(str_detect(dist_method, "sample"), iter, 1),
#          dist_method = ifelse(str_detect(dist_method, "sample"), "sample", dist_method),
#          dist_method = str_replace_all(dist_method, dist_method_replace),
#          dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
#   left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait")) %>%
#   group_by(environment, trait, dist_method, iter) %>%
#   mutate(max_accuracy_adding = accuracy[n_e == max(n_e)]) %>%
#   ungroup()
# 
# 
# ## Summarize the accuracies and create CIs for the random method
# cumulative_pred_results_summ <- cumulative_pred_results %>%
#   group_by(environment, trait, dist_method, dist_method_abbr, n_e, max_accuracy, max_accuracy_adding) %>% 
#   summarize_at(vars(accuracy), funs(mean, sd, n())) %>% 
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
#   ungroup()
# 
# 
# # Split by trait
# cumulative_pred_results_split <- cumulative_pred_results_summ %>% 
#   split(.$trait) %>% 
#   purrr::map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy_adding, decreasing = TRUE)]))))
# 
# # Plot
# g_plotlist <- cumulative_pred_results_split %>%
#   purrr::map(~ggplot(data = ., aes(x = n_e, y = mean, color = dist_method_abbr)) + 
#                geom_hline(aes(yintercept = max_accuracy_adding, lty = "Accuracy\nUsing\nAll Data")) +
#                # geom_point() + 
#                geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#                geom_line() + 
#                scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#                scale_x_continuous(breaks = pretty) +
#                scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#                facet_wrap(~ environment, ncol = 5) + 
#                xlab("Number of training environments") +
#                ylab("Prediction accuracy") +
#                theme_minimal() )
# 
# # Save
# for (i in seq_along(g_plotlist)) {
#   ggsave(filename = paste0("cumulative_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
#          height = 10, width = 8, dpi = 1000, path = fig_dir)
# }
# 
# 
# ## Highlight STP16
# stp_subset <- cumulative_pred_results_summ %>%
#   filter(environment == "STP16", trait != "HeadingDateAGDD") %>%
#   filter(dist_method_abbr %in% c("GrCD", "PhnD", "TYEC", "Rndm"))
# 
# g_cumulative_pred_example <- stp_subset %>%
#   # filter(dist_method_abbr %in% c("GrCD", "Rndm")) %>%
#   ggplot(data = ., aes(x = n_e, y = mean, color = dist_method_abbr)) +
#   geom_hline(aes(yintercept = max_accuracy_adding, lty = "Accuracy using\nall environments")) +
#   geom_point() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   geom_line() +
#   geom_blank(data = stp_subset, aes(x = n_e, y = mean), inherit.aes = FALSE) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure", guide = guide_legend(order = 1)) +
#   scale_x_continuous(breaks = pretty) +
#   scale_linetype_manual(values = 2, name = NULL, guide = guide_legend(order = 2)) +
#   facet_wrap(~ trait, scales = "free_y") +
#   xlab("Number of training environments") +
#   ylab("Prediction accuracy") +
#   theme_presentation2() + theme(legend.position = "bottom", legend.direction = "horizontal")
# 
# # Save
# ggsave(filename = "cumulative_environment_prediction_STP16_1.jpg", plot = g_cumulative_pred_example,
#        height = 5, width = 10, dpi = 1000, path = fig_dir)
# 
# 
# 
# 
# ## For each trait and environment, calculate the difference between the accuracy
# ## using the nth training environment and the max accuracy
# ## Then summarize for each trait
# cumulative_pred_diff <- cumulative_pred_results %>% 
#   mutate(diff_accuracy = accuracy - max_accuracy_adding) %>%
#   group_by(trait, dist_method_abbr, n_e, iter) %>% 
#   summarize(diff_accuracy = mean(diff_accuracy), max_accuracy = mean(max_accuracy)) %>% ## Take the mean over all environments for a method/iteration
#   summarize_at(vars(diff_accuracy, max_accuracy), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   select(trait:n_e, max_accuracy = max_accuracy_mean, mean = diff_accuracy_mean, sd = diff_accuracy_sd, n = diff_accuracy_n) %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot - Random
# g_diff_accuracy <- cumulative_pred_diff %>% 
#   ggplot(aes(x = n_e, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = 0, lty = "Accuracy using\nall environments")) +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = 2, name = NULL) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
#   xlab("Number of training environments") +
#   ylab("Prediction accuracy (relative to using all data)") +
#   theme_presentation() +
#   theme(strip.placement = "outside")
# 
# # Save
# ggsave(filename = "cumulative_environment_prediction_relative.jpg", plot = g_diff_accuracy,
#        height = 8, width = 7, dpi = 1000, path = fig_dir)
# 
# 
# 
# # Plot - examples
# g_diff_accuracy <- cumulative_pred_diff %>% 
#   mutate_at(vars(mean, lower, upper), funs(. + max_accuracy)) %>%
#   filter(trait != "HeadingDateAGDD", dist_method_abbr %in% c("GrCD", "PhnD", "TYEC", "Rndm")) %>% 
#   ggplot(aes(x = n_e, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy using\nall environments")) +
#   geom_point() +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure", guide = guide_legend(order = 1)) +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = 2, name = NULL, guide = guide_legend(order = 3)) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_wrap(~ trait, scales = "free_y", nrow = 1) + 
#   xlab("Number of training environments") +
#   ylab("Prediction accuracy") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# # Save
# ggsave(filename = "cumulative_environment_prediction_relative_example.jpg", plot = g_diff_accuracy,
#        height = 5, width = 10, dpi = 1000, path = fig_dir)
# 
# 
# ## Choose an arbitrary number of environments and show the average accuracy advantage
# cumulative_pred_diff_example <- cumulative_pred_diff %>% 
#   filter(n_e %in% c(3, 5, 10)) %>% 
#   select(trait:n_e, mean) %>% 
#   mutate(mean = round(mean, 3)) %>% 
#   spread(dist_method_abbr, mean) %>% 
#   arrange(n_e)
# 
# # Save as an image
# t_cumulative_pred_diff <- grid.arrange(tableGrob(cumulative_pred_diff_example, rows = NULL, theme = ttheme_minimal()))
# 
# ggsave(filename = "cumulative_environment_prediction_example.jpg", plot = t_cumulative_pred_diff, path = fig_dir,
#        height = 5, width = 7, dpi = 1000)
# 
# 
# ## For each distance method, find the average number of environments in which the accuracy is maximized
# cumulative_pred_nE <- cumulative_pred_results %>% 
#   group_by(trait, dist_method_abbr, iter, environment) %>%
#   top_n(x = ., n = 1, wt = accuracy) %>% summarize(n_e = mean(n_e))
# 
# cumulative_pred_nE_summ <- cumulative_pred_nE %>%
#   summarize(n_e = mean(n_e)) %>%
#   summarize_at(vars(n_e), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot
# g_n_train <- cumulative_pred_nE_summ %>% 
#   ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
#   geom_point(position = position_dodge(0.6), size = 3) +
#   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_y_continuous(breaks = pretty) +
#   ylab("Peak accuracy training environments") +
#   theme_presentation() +
#   theme(axis.title.x = element_blank())
# 
# ggsave(filename = "cumulative_max_accuracy_nE.jpg", plot = g_n_train, path = fig_dir,
#        height = 5, width = 8, dpi = 1000)
# 
# 
# ## Using the phenotypic distance measure, calculate the average distance at which
# ## the maximum prediction accuracy is achieved
# 
# # First calculate the mean distance for each training set
# average_pheno_distance <- pred_env_dist_rank$tp %>% 
#   filter(model == "pheno_dist") %>% 
#   mutate(env_rank = map(env_rank, ~unlist(.) %>% data_frame(pred_environment = names(.), distance = .) %>% mutate(n_e = seq(nrow(.))))) %>% 
#   unnest() %>%
#   group_by(trait, environment) %>% 
#   mutate(distance = cummean(distance)) %>%
#   ungroup() %>%
#   select(-model)
# 
# 
# cumulative_pred_nE_summ_avg <- cumulative_pred_nE %>% 
#   left_join(., average_pheno_distance) %>%
#   summarize(distance = mean(distance)) %>%
#   summarize_at(vars(distance), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# 
# # Plot
# g_dist<- cumulative_pred_nE_summ_avg %>% 
#   ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
#   geom_point(position = position_dodge(0.6), size = 3) +
#   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_y_continuous(breaks = function(x) pretty(x, n = 3)) +
#   ylab("Peak accuracy training environments") +
#   facet_wrap(~trait, scales = "free", nrow = 1) +
#   theme_presentation() +
#   theme(axis.title.x = element_blank(), strip.text = element_blank())
# 
# ggsave(filename = "cumulative_max_accuracy_pheno_dist.jpg", plot = g_dist, path = fig_dir,
#        height = 5, width = 10, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Window predictions
# 
# 
# # Rename the distance methods
# window_pred_results <- cluster_pred_out_window %>% 
#   mutate(out = map(out, ~mutate(., window = seq(nrow(.))))) %>%
#   unnest(out) %>%
#   rename(dist_method = model) %>%
#   mutate(iter = parse_number(dist_method),
#          iter = ifelse(str_detect(dist_method, "sample"), iter, 1),
#          dist_method = ifelse(str_detect(dist_method, "sample"), "sample", dist_method),
#          dist_method = str_replace_all(dist_method, dist_method_replace),
#          dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
#   left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait")) 
# 
# 
# ## Summarize the accuracies and create CIs for the random method
# window_pred_results_summ <- window_pred_results %>%
#   group_by(environment, trait, dist_method, dist_method_abbr, window, max_accuracy) %>% 
#   summarize_at(vars(accuracy), funs(mean, sd, n())) %>% 
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
#   ungroup()
# 
# 
# # Split by trait
# window_pred_results_split <- window_pred_results_summ %>% 
#   split(.$trait) %>% 
#   map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy, decreasing = TRUE)]))))
# 
# # Plot
# g_plotlist <- window_pred_results_split %>%
#   map(~ggplot(data = ., aes(x = window, y = mean, color = dist_method_abbr)) + 
#         geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy\nUsing\nAll Data")) +
#         # geom_point() + 
#         geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#         geom_line() + 
#         scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#         scale_x_continuous(breaks = pretty) +
#         scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#         facet_wrap(~ environment, ncol = 5) + 
#         xlab("Window number") +
#         ylab("Prediction accuracy") +
#         theme_minimal() )
# 
# # Save
# for (i in seq_along(g_plotlist)) {
#   ggsave(filename = paste0("window_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
#          height = 10, width = 8, dpi = 1000, path = fig_dir)
# }
# 
# 
# 
# 
# 
# 
# # For each trait and environment, calculate the difference between the accuracy
# ## using the nth training environment and the max accuracy
# ## Then summarize for each trait
# window_pred_diff <- window_pred_results %>% 
#   mutate(diff_accuracy = accuracy - max_accuracy) %>%
#   group_by(trait, dist_method_abbr, window, iter) %>% 
#   summarize(diff_accuracy = mean(diff_accuracy)) %>% ## Take the mean over all environments for a method/iteration
#   summarize_at(vars(diff_accuracy), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot - Random
# g_diff_accuracy <- window_pred_diff %>% 
#   ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
#   xlab("Window number") +
#   ylab("Prediction accuracy (relative to using all data)") +
#   theme_presentation() +
#   theme(strip.placement = "outside")
# 
# # Save
# ggsave(filename = "window_environment_prediction_relative.jpg", plot = g_diff_accuracy,
#        height = 8, width = 7, dpi = 1000, path = fig_dir)
# 
# 
# ## Fit smoothing lines
# g_diff_accuracy_smooth <- window_pred_diff %>% 
#   ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) + 
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
#   xlab("Window number") +
#   ylab("Prediction accuracy (relative to using all data)") +
#   theme_presentation() +
#   theme(strip.placement = "outside")
# 
# # Save
# ggsave(filename = "window_environment_prediction_relative_smooth.jpg", plot = g_diff_accuracy_smooth,
#        height = 8, width = 7, dpi = 1000, path = fig_dir)
# 
# 
# 
# ## For each distance method, find the average window number in which the accuracy is maximized
# window_pred_wn <- window_pred_results %>% 
#   group_by(trait, dist_method_abbr, iter, environment) %>%
#   top_n(x = ., n = 1, wt = accuracy) %>% summarize(window = mean(window))
# 
# window_pred_wn_summ <- window_pred_wn %>%
#   summarize(window = mean(window)) %>%
#   summarize_at(vars(window), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot
# g_window_train <- window_pred_wn_summ %>% 
#   ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
#   geom_point(position = position_dodge(0.6), size = 3) +
#   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_y_continuous(breaks = pretty) +
#   ylab("Peak accuracy training environments") +
#   theme_presentation() +
#   theme(axis.title.x = element_blank())
# 
# ggsave(filename = "window_max_accuracy_wn.jpg", plot = g_window_train, path = fig_dir,
#        height = 5, width = 8, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### Environment covariance matrix prediction ####
# 
# # Load the results
# load(file.path(result_dir, "env_cov_mat_predictions.RData"))
# 
# 
# 
# ## MC predictions
# mc_pred_tidy <- environment_mc_predictions %>%
#   unnest() %>% unnest() %>% 
#   select(-trait1, -trait2) %>%
#   mutate_at(vars(trait, model, environment), as.factor)
# 
# ## Summarize the correlations across all environments for each iteration
# mc_pred_summ <- mc_pred_tidy %>% 
#   group_by(trait, pTrainEnv, model, environment, iter) %>% 
#   summarize(accuracy = cor(value, pgv))
# 
# # Now take the mean over iterations
# mc_pred_summ1 <- mc_pred_summ %>% 
#   summarize(accuracy = mean(accuracy))
# 
# # Plot for each model and trait
# g_model_acc <- mc_pred_summ1 %>%
#   ggplot(aes(x = trait, y = accuracy, fill = model)) +
#   geom_boxplot(position = "dodge", alpha = 0.5) +
#   xlab("Trait") + 
#   ylab("Prediction accuracy") +
#   scale_fill_discrete(name = "Model") +
#   facet_grid(~ pTrainEnv) +
#   theme_presentation2()
# 
# ggsave(filename = "environmental_cov_model_accuracy.jpg", plot = g_model_acc, path = fig_dir, width = 10, height = 4, dpi = 1000)
# 








