## S2MET Predictions
## 
## Script for analyzing cross-validation results
## 
## Author: Jeff Neyhart
## Last modified: May 14, 2019
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# # Load other packages
library(effects) # For ls means / marginal means
library(broom)
library(car)
library(lmerTest)


## Load the validation results
file_list <- list.files(result_dir, pattern = "predictions_out.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))


## Significant level
alpha <- 0.05


## A vector to rename models
model_replace <- c("model1" = "G + e", "model2" = "G + E", "model3" = "G + E + GE", "model4" = "G + E(b + 1)")


###################################
### Leave-one-environment-out
###################################

#

## Calculate accuracy per environment, model, and population
loeo_predictions_ability <- loeo_predictions_df %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  group_by(trait, environment, model, pop) %>%
  summarize(ability = cor(prediction, value), bias = mean(prediction - value)) %>%
  ungroup() %>%
  mutate(model = str_replace_all(model, model_replace))

## Adjust ability using heritability
loeo_predictions_accuracy <- loeo_predictions_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))


### Plot accuracy
loeo_predictions_accuracy %>%
  ggplot(aes(x = trait, y = accuracy, color = model)) +
  geom_boxplot(position = position_dodge(0.9)) +
  facet_grid(~ pop) +
  theme_presentation(16)


## Calculate average
loeo_predictions_accuracy %>%
  group_by(trait, pop, model) %>%
  summarize(accuracy = mean(accuracy)) %>%
  spread(model, accuracy)




###################################
### Leave-one-year-out
###################################

#

## Calculate accuracy per environment, model, and population
loyo_predictions_ability <- loyo_predictions_df %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  group_by(trait, year, environment, model, pop) %>%
  summarize(ability = cor(prediction, value), bias = mean(prediction - value)) %>%
  ungroup() %>%
  mutate(model = as.factor(str_replace_all(model, model_replace)))

## Adjust ability using heritability
loyo_predictions_accuracy <- loyo_predictions_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))


### Plot accuracy
loyo_predictions_accuracy %>%
  ggplot(aes(x = pop, y = accuracy, color = model)) +
  geom_boxplot() +
  scale_color_discrete(name = "Model", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  facet_grid(trait ~ ., switch = "y") +
  theme_presentation2(16)



##### Compare LOEO and LOYO #####

loeo_predictions_accuracy %>%
  group_by(trait, model, pop) %>%
  summarize(accuracy = mean(accuracy)) %>%
  spread(model, accuracy)



## Calculate average - weight by number of environments in a year
loyo_predictions_accuracy %>%
  group_by(trait, pop, model, year) %>%
  mutate(weight = n()) %>%
  summarize(accuracy = mean(accuracy), weight = mean(weight)) %>%
  summarize(accuracy = weighted.mean(x = accuracy, w = weight)) %>%
  spread(model, accuracy)



































# # ## combine all data and assign the scheme number
# # pov00_predictions <- rename(pov00_predictions, environment = val_environment)
# # 
# # results_df <- map(object_list, get) %>%
# #   map_df(ungroup) %>%
# #   select(-core) %>%
# #   mutate(scheme_number = str_extract(scheme, "[0-9]{1,2}")) %>%
# #   mutate_at(vars(model, scheme), as.factor)
# # 
# # 
# # ## 
# # ## Fit a model
# # ## 
# # 
# # model_results <- results_df %>% 
# #   group_by(trait, scheme_number) %>%
# #   nest() %>%
# #   mutate(out = list(NULL))
# # 
# # for (i in seq(nrow(model_results))) {
# #   
# #   df <- model_results$data[[i]]
# #   
# #   if (model_results$scheme_number[i] == "00") {
# #     next
# #     
# #   } else {
# #     
# #     fit <- lmer(accuracy ~ model + scheme + model:scheme + (1|environment) + (1|environment:model) + 
# #                   (1|environment:scheme)+ (1|environment:model:scheme), 
# #                 data = df)
# #   }
# #   
# #   out <- data_frame(model = list(fit),
# #                     effects = list(as.data.frame(Effect(focal.predictors = c("model", "scheme"), fit))),
# #                     ranova = list(tidy(ranova(fit))),
# #                     anova = list(tidy(anova(fit))))
# #   
# #   model_results$out[[i]] <- out
# #   
# # }
# #                     
# # 
# # ## Unnest
# # model_results1 <- model_results %>%
# #   filter(!map_lgl(out, is.null)) %>%
# #   unnest(out)
# # 
# # 
# # ## Plot
# # model_results1 %>% 
# #   unnest(effects) %>%
# #   # filter(scheme_number == "2") %>%
# #   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
# #   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
# #   geom_point(position = position_dodge(0.7)) +
# #   facet_grid(trait ~ scheme_number, scales = "free_x") +
# #   theme_presentation2()
# 
# 
# 
# 
# ## Model within scheme
# ## 
# ## POV00
# ## 
# pov00_analysis <- pov00_predictions %>% 
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# for (i in seq(nrow(pov00_analysis))) {
#   
#   df <- pov00_analysis$data[[i]]
#   
#   # Fit a model
#   # fit <- lmer(zscore ~ model + (1|environment), data = df)
#   fit <- lmer(accuracy ~ model + (1|environment), data = df)
#   
#   
#   ## Rescale accuracy and return
#   pov00_analysis$out[[i]] <- tibble(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# pov00_analysis1 <- unnest(pov00_analysis, out)
# 
# ## Plot
# g_pov00 <- pov00_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "pov00_accuracy_analysis.jpg", plot = g_pov00, path = fig_dir, width = 7, height = 4, dpi = 1000)
# 
# 
# 
# 
# 
# ## Model within scheme
# ## 
# ## POV1
# ## 
# 
# # Histogram
# qplot(x = accuracy, data = pov1_predictions, geom = "density", facets = ~trait)
# 
# 
# pov1_analysis <- pov1_predictions %>% 
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# for (i in seq(nrow(pov1_analysis))) {
#   
#   df <- pov1_analysis$data[[i]]
#   
#   # Fit a model
#   # fit <- lmer(zscore ~ model + (1|environment), data = df, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
#   fit <- lmer(accuracy ~ model + (1|environment), data = df)
#   
#   
#   ## Rescale accuracy and return
#   pov1_analysis$out[[i]] <- data_frame(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# pov1_analysis1 <- unnest(pov1_analysis, out)
# 
# ## Plot
# g_pov1 <- pov1_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "pov1_accuracy_analysis.jpg", plot = g_pov1, path = fig_dir, width = 7, height = 4, dpi = 1000)
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
# ## Model within scheme
# ## 
# ## CV1
# ## 
# 
# # Histogram
# qplot(x = accuracy, data = cv1_predictions, geom = "density", facets = scheme~trait)
# 
# 
# ## Variance heterogeneity between schemes within models?
# cv1_predictions %>%
#   mutate(zscore = ztrans(accuracy),
#          model = as.factor(model)) %>%
#   group_by(trait, scheme) %>%
#   do(test = leveneTest(zscore ~ model, data = .)) %>%
#   pull(test)
# 
# ## Model scheme separately
# cv1_analysis <- cv1_predictions %>% 
#   mutate_at(vars(model, rep), as.factor) %>%
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# 
# for (i in seq(nrow(cv1_analysis))) {
#   
#   df <- cv1_analysis$data[[i]]
#   
#   # Fit a model
#   # fit <- lmer(zscore ~ model + (1|environment) + (1|environment:model) + (1|rep:model), data = df)
#   fit <- lmer(accuracy ~ model + (1|environment:model) + (1|rep:model), data = df)
#   
#   
#   ## Rescale accuracy and return
#   cv1_analysis$out[[i]] <- data_frame(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# cv1_analysis1 <- unnest(cv1_analysis, out)
# 
# ## Plot
# g_cv1 <- cv1_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "cv1_accuracy_analysis.jpg", plot = g_cv2, path = fig_dir, width = 7, height = 4, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# ## Model within scheme
# ## 
# ## CV2
# ## 
# 
# # Histogram
# qplot(x = accuracy, data = cv2_predictions, geom = "density", facets = ~trait)
# 
# 
# ## Model scheme separately
# cv2_analysis <- cv2_predictions %>% 
#   mutate_at(vars(model, rep), as.factor) %>%
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# 
# for (i in seq(nrow(cv2_analysis))) {
#   
#   df <- cv2_analysis$data[[i]]
#   
#   # Fit a model
#   fit <- lmer(zscore ~ model + (1|environment) + (1|environment:model) + (1|model:rep), data = df)
#   
#   ## Rescale accuracy and return
#   cv2_analysis$out[[i]] <- data_frame(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# cv2_analysis1 <- unnest(cv2_analysis, out)
# 
# ## Plot
# g_cv2 <- cv2_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "cv2_accuracy_analysis.jpg", plot = g_cv2, path = fig_dir, width = 7, height = 4, dpi = 1000)
# 
# 
# 
# 
# ## Model within scheme
# ## 
# ## POCV2
# ## 
# 
# # Histogram
# qplot(x = accuracy, data = pocv2_predictions, geom = "density", facets = scheme~trait)
# 
# ## Model scheme separately
# pocv2_analysis <- pocv2_predictions %>% 
#   mutate_at(vars(model, rep), as.factor) %>%
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# 
# 
# for (i in seq(nrow(pocv2_analysis))) {
#   
#   df <- pocv2_analysis$data[[i]]
#   
#   # Fit a model
#   fit <- lmer(zscore ~ model + (1|environment) + (1|environment:model) + (1|model:rep), data = df)
#   
#   ## Rescale accuracy and return
#   pocv2_analysis$out[[i]] <- data_frame(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# pocv2_analysis1 <- unnest(pocv2_analysis, out)
# 
# ## Plot
# g_pocv2 <- pocv2_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "pocv2_accuracy_analysis.jpg", plot = g_pocv2, path = fig_dir, width = 7, height = 4, dpi = 1000)
# 
# 
# 
# 
# 
# 
# ## Model within scheme
# ## 
# ## CV0 and POCV0
# ## 
# 
# # Histogram
# qplot(x = accuracy, data = cv0_predictions, geom = "density", facets = scheme~trait)
# 
# ## Model scheme separately
# cv0_analysis <- cv0_predictions %>% 
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# 
# 
# for (i in seq(nrow(cv0_analysis))) {
#   
#   df <- cv0_analysis$data[[i]]
#   
#   # Fit a model
#   fit <- lmer(zscore ~ model + (1|environment), data = df)
#   
#   ## Rescale accuracy and return
#   cv0_analysis$out[[i]] <- data_frame(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# cv0_analysis1 <- unnest(cv0_analysis, out)
# 
# ## Plot
# g_cv0 <- cv0_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "cv0_accuracy_analysis.jpg", plot = g_cv0, path = fig_dir, width = 7, height = 4, dpi = 1000)
# 
# 
# 
# 
# ## Model within scheme
# ## 
# ## CV00 and POCV00
# ## 
# 
# # Histogram
# qplot(x = accuracy, data = cv00_predictions, geom = "density", facets = scheme~trait)
# 
# ## Model scheme separately
# cv00_analysis <- cv00_predictions %>% 
#   mutate_at(vars(model, rep), as.factor) %>%
#   group_by(trait, scheme) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# 
# 
# for (i in seq(nrow(cv00_analysis))) {
#   
#   df <- cv00_analysis$data[[i]]
#   
#   # Fit a model
#   fit <- lmer(zscore ~ model + (1|environment) + (1|environment:model) + (1|rep:model), data = df)
#   
#   ## Rescale accuracy and return
#   cv00_analysis$out[[i]] <- data_frame(
#     model = list(fit),
#     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
#     ranova = list(tidy(ranova(fit))),
#     anova = list(tidy(anova(fit)))
#   )
#   
# }
# 
# cv00_analysis1 <- unnest(cv00_analysis, out)
# 
# ## Plot
# g_cv00 <- cv00_analysis1 %>% 
#   unnest(effects) %>%
#   mutate(model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 3) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
#   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = "cv00_accuracy_analysis.jpg", plot = g_cv00, path = fig_dir, width = 7, height = 4, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# ####
# # All analyses
# ####
# 
# 
# ## Plot everything
# all_analyses <- ls(pattern = "analysis1") %>%
#   subset(. != "cv2_analysis1") %>%
#   map_df(get) %>%
#   mutate(scheme_number = str_extract(scheme, "[0-9]{1,2}"),
#          scheme_code = str_replace_all(scheme_number, scheme_number_replace),
#          scheme_code = factor(scheme_code, levels = rev(scheme_number_replace)))
# 
# 
# g_all_validation <- all_analyses %>% 
#   unnest(effects) %>%
#   # Remove some schema
#   filter(!scheme %in% c("pocv00", "pocv1")) %>%
#   # Rename some schema
#   mutate(scheme = str_replace(scheme, "pocv", "pov"),
#          model = str_replace_all(model, model_replace),
#          model = factor(model, levels = model_replace)) %>%
#   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
#   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
#   geom_point(position = position_dodge(0.7), size = 2) +
#   # facet_grid(trait ~ scheme_number, scales = "free_x") +
#   facet_grid(trait ~ scheme_code, scales = "free_x", space = "free_x", 
#              labeller = labeller(trait = str_add_space, scheme_code = label_parsed)) +
#   scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
#   scale_x_discrete(name = "Validation scheme", labels = toupper) +
#   scale_color_discrete(name = "Model") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# 
# ggsave(filename = "all_accuracy_analysis.jpg", plot = g_all_validation, path = fig_dir, width = 8, height = 6, dpi = 1000)
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
# # #### Sample results
# # ## Load the validation results
# # file_list <- list.files(result_dir, "results_sample.RData$", full.names = TRUE)
# # object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))
# # # Rename
# # invisible(map(object_list, ~assign(x = paste0(., "_sample"), value = get(.), envir = .GlobalEnv)))
# # object_list <- ls(pattern = "predictions_sample$")
# # 
# # ## Fix cv0
# # cv0_predictions_sample <- mutate(cv0_predictions_sample, scheme = str_replace(scheme, "00", "0"))
# # 
# # 
# # results_df <- map(object_list, get) %>%
# #   map_df(ungroup) %>%
# #   select(-core) %>%
# #   rename(rep = .id) %>%
# #   mutate(scheme_number = str_extract(scheme, "[0-9]{1,2}")) %>%
# #   mutate_at(vars(model, scheme, rep), as.factor)
# # 
# # 
# # ## 
# # ## Fit models
# # ## 
# # ## 
# # 
# # 
# # ## POV00
# # ## 
# # ## 
# # 
# # pov00_analysis_sample <- results_df %>% 
# #   filter(scheme == "pov00") %>%
# #   mutate(zscore = ztrans(accuracy)) %>%
# #   group_by(trait, scheme) %>%
# #   nest() %>%
# #   mutate(out = list(NULL))
# # 
# # for (i in seq(nrow(pov00_analysis_sample))) {
# #   
# #   df <- pov00_analysis_sample$data[[i]]
# #   
# #   # Fit a model
# #   fit <- lmer(zscore ~ model + (1|environment:rep), data = df)
# #   
# #   ## Rescale accuracy and return
# #   pov00_analysis_sample$out[[i]] <- data_frame(
# #     model = list(fit),
# #     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
# #     ranova = list(tidy(ranova(fit))),
# #     anova = list(tidy(anova(fit)))
# #   )
# #   
# # }
# # 
# # pov00_analysis_sample1 <- unnest(pov00_analysis_sample, out)
# # 
# # ## Plot
# # g_pov <- pov00_analysis_sample1 %>% 
# #   unnest(effects) %>%
# #   # mutate(model = str_replace_all(model, model_replace),
# #   #        model = factor(model, levels = model_replace)) %>%
# #   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
# #   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
# #   geom_point(position = position_dodge(0.7), size = 3) +
# #   # facet_grid(trait ~ scheme_number, scales = "free_x") +
# #   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
# #   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
# #   scale_x_discrete(name = "Validation scheme", labels = toupper) +
# #   scale_color_discrete(name = "Model") +
# #   theme_presentation2() +
# #   theme(legend.position = "bottom")
# # 
# # 
# # 
# # 
# # ## POV1
# # ## 
# # ## 
# # 
# # pov1_analysis_sample <- results_df %>% 
# #   filter(scheme == "pov1") %>%
# #   mutate(zscore = ztrans(accuracy)) %>%
# #   group_by(trait, scheme) %>%
# #   nest() %>%
# #   mutate(out = list(NULL))
# # 
# # for (i in seq(nrow(pov1_analysis_sample))) {
# #   
# #   df <- pov1_analysis_sample$data[[i]]
# #   
# #   # Fit a model
# #   fit <- lmer(zscore ~ 1 + model + (1|environment), data = df)
# #   
# #   ## Rescale accuracy and return
# #   pov1_analysis_sample$out[[i]] <- data_frame(
# #     model = list(fit),
# #     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
# #     ranova = list(tidy(ranova(fit))),
# #     anova = list(tidy(anova(fit)))
# #   )
# #   
# # }
# # 
# # pov1_analysis_sample1 <- unnest(pov1_analysis_sample, out)
# # 
# # ## Plot
# # g_pov1 <- pov1_analysis_sample1 %>% 
# #   unnest(effects) %>%
# #   # mutate(model = str_replace_all(model, model_replace),
# #   #        model = factor(model, levels = model_replace)) %>%
# #   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
# #   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
# #   geom_point(position = position_dodge(0.7), size = 3) +
# #   # facet_grid(trait ~ scheme_number, scales = "free_x") +
# #   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
# #   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
# #   scale_x_discrete(name = "Validation scheme", labels = toupper) +
# #   scale_color_discrete(name = "Model") +
# #   theme_presentation2() +
# #   theme(legend.position = "bottom")
# # 
# # 
# # 
# # ## CV00
# # ## 
# # ## 
# # 
# # cv00_analysis_sample <- results_df %>% 
# #   filter(scheme %in% c("cv00", "pocv00")) %>%
# #   mutate(zscore = ztrans(accuracy)) %>%
# #   group_by(trait, scheme) %>%
# #   nest() %>%
# #   mutate(out = list(NULL))
# # 
# # for (i in seq(nrow(cv00_analysis_sample))) {
# #   
# #   df <- cv00_analysis_sample$data[[i]]
# #   
# #   # Fit a model
# #   fit <- lmer(zscore ~ 1 + model + (1|environment:rep), data = df)
# #   
# #   ## Rescale accuracy and return
# #   cv00_analysis_sample$out[[i]] <- data_frame(
# #     model = list(fit),
# #     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
# #     ranova = list(tidy(ranova(fit))),
# #     anova = list(tidy(anova(fit)))
# #   )
# #   
# # }
# # 
# # cv00_analysis_sample1 <- unnest(cv00_analysis_sample, out)
# # 
# # ## Plot
# # g_cv00 <- cv00_analysis_sample1 %>% 
# #   unnest(effects) %>%
# #   # mutate(model = str_replace_all(model, model_replace),
# #   #        model = factor(model, levels = model_replace)) %>%
# #   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
# #   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
# #   geom_point(position = position_dodge(0.7), size = 3) +
# #   # facet_grid(trait ~ scheme_number, scales = "free_x") +
# #   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
# #   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
# #   scale_x_discrete(name = "Validation scheme", labels = toupper) +
# #   scale_color_discrete(name = "Model") +
# #   theme_presentation2() +
# #   theme(legend.position = "bottom")
# # 
# # 
# # 
# # 
# # 
# # ## CV0
# # ## 
# # ## 
# # 
# # cv0_analysis_sample <- results_df %>% 
# #   filter(scheme %in% c("cv0", "pocv0")) %>%
# #   mutate(zscore = ztrans(accuracy)) %>%
# #   group_by(trait, scheme) %>%
# #   nest() %>%
# #   mutate(out = list(NULL))
# # 
# # for (i in seq(nrow(cv0_analysis_sample))) {
# #   
# #   df <- cv0_analysis_sample$data[[i]]
# #   
# #   # Fit a model
# #   fit <- lmer(zscore ~ 1 + model + (1|environment:rep), data = df)
# #   
# #   ## Rescale accuracy and return
# #   cv0_analysis_sample$out[[i]] <- data_frame(
# #     model = list(fit),
# #     effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
# #     ranova = list(tidy(ranova(fit))),
# #     anova = list(tidy(anova(fit)))
# #   )
# #   
# # }
# # 
# # cv0_analysis_sample1 <- unnest(cv0_analysis_sample, out)
# # 
# # ## Plot
# # g_cv0 <- cv0_analysis_sample1 %>% 
# #   unnest(effects) %>%
# #   # mutate(model = str_replace_all(model, model_replace),
# #   #        model = factor(model, levels = model_replace)) %>%
# #   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
# #   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
# #   geom_point(position = position_dodge(0.7), size = 3) +
# #   # facet_grid(trait ~ scheme_number, scales = "free_x") +
# #   facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
# #   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
# #   scale_x_discrete(name = "Validation scheme", labels = toupper) +
# #   scale_color_discrete(name = "Model") +
# #   theme_presentation2() +
# #   theme(legend.position = "bottom")
# # 
# # 
# # 
# # ## Plot all
# # all_analyses_sample <- ls(pattern = "analysis_sample1") %>%
# #   subset(. != "cv2_analysis1") %>%
# #   map_df(get) %>%
# #   mutate(scheme_number = str_extract(scheme, "[0-9]{1,2}"))
# # 
# # 
# # g_all_validation <- all_analyses_sample %>% 
# #   unnest(effects) %>%
# #   # mutate(model = str_replace_all(model, model_replace),
# #   #        model = factor(model, levels = model_replace)) %>%
# #   ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
# #   geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
# #   geom_point(position = position_dodge(0.7), size = 2) +
# #   # facet_grid(trait ~ scheme_number, scales = "free_x") +
# #   facet_grid(trait ~ scheme_number, scales = "free_x", space = "free_x", labeller = labeller(trait = str_add_space)) +
# #   scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
# #   scale_x_discrete(name = "Validation scheme", labels = toupper) +
# #   scale_color_discrete(name = "Model") +
# #   theme_presentation2() +
# #   theme(legend.position = "bottom")
# # 
# # 
# # ggsave(filename = "all_accuracy_analysis_sample.jpg", plot = g_all_validation, path = fig_dir, width = 8, height = 6, dpi = 1000)



