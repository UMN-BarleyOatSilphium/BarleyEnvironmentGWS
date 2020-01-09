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
library(patchwork)
library(cowplot)


## Load the validation results
file_list <- list.files(result_dir, pattern = "predictions_out.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))


## Significant level
alpha <- 0.05


## A vector to rename models
model_replace <- c("model1" = "G + e", "model2" = "G + E", "model3" = "G + E + GE", "model4" = "G + E(b + 1)")
f_model_replace <- function(x) str_replace_all(x, model_replace)

###################################
### Leave-one-environment-out
###################################

## Add pop information
loeo_predictions_df <- loeo_predictions_out %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp"))
  


## Calculate accuracy and bias per environment, model, and population
## Also calculate accuracy and bias across all environments for a model/population
loeo_predictions_ability <- loeo_predictions_df %>%
  group_by(trait, environment, model, pop) %>%
  summarize(ability = cor(prediction, value), bias = mean(prediction - value, na.rm = T)) %>%
  bind_rows(., mutate(loeo_predictions_df, environment = "All") %>% group_by(trait, environment, model, pop) %>%
              summarize(ability = cor(prediction, value, use = "complete.obs"), bias = mean(prediction - value, na.rm = T))) %>%
  ungroup()

## Adjust ability using heritability
loeo_predictions_accuracy <- loeo_predictions_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))


## Create a table with trait, model, population, total accuracy and range of accuracy within environment
loeo_predictions_accuracy_summ <- loeo_predictions_accuracy %>%
  group_by(trait, model, pop) %>%
  summarize_at(vars(ability, bias), list(all = ~.[environment == "All"], mean = ~mean(.[environment != "All"], na.rm = T),
                                         min = ~min(.[environment != "All"], na.rm = T), max = ~max(.[environment != "All"], na.rm = T))) %>%
  ungroup()



## Plot predicted and observed value, with accuracy for each environment
## and overall

# First create annotation for accuracy
loeo_predictions_df1 <- loeo_predictions_accuracy %>%
  mutate(annotation = paste0(environment, "~(r[MP]==", formatC(x = ability, 2, format = "g"), ")")) %>%
  left_join(loeo_predictions_df, .)


# Plot predicted versus observed value
g_loeo_prediction_list <- loeo_predictions_df1 %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Add accuracy summaries
    df1 <- left_join(df, loeo_predictions_accuracy_summ, by = c("trait", "model", "pop"))
    
    
    max_wid <- unique(df1$max_nchar)
    # Scale truncation function
    scale_trunc <- function(x) str_pad(string = x, width = max_wid, side = "left", pad = " ")
    
    # Create a separate legend df
    x1 <- mutate(df1, x = ifelse(model == "model1", 0.20, 0.80), y = 0.10)
    
    ## Split x by model and population
    x_split <- split(x1, list(x1$pop, x1$model))
    
    # Map over this split
    g_list <- vector("list", length = length(x_split))
    
    for (i in seq_along(g_list)) {
      x2 <- x_split[[i]]
    
      g_x <- ggplot(data = x2, aes(x = prediction, y = value, color = annotation)) +
        geom_point(size = 0.5) +
        scale_y_continuous(name = "Observation", breaks = pretty, labels = scale_trunc) +
        scale_x_continuous(name = "Prediction", breaks = pretty) +
        scale_color_discrete(name = NULL, labels = function(x) parse(text = x), guide = guide_legend(ncol = 1)) +
        theme_presentation2(10) +
        theme(legend.position = c(x2$x[1], x2$y[1]), legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"),
              legend.text = element_text(size = 6),
              # No axis titles
              axis.title = element_blank())
      
      ## Add faceting depending on position
      if (i == 1) {
        g_list[[i]] <- g_x + 
          facet_grid(trait ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace))
      } else {
        g_list[[i]] <- g_x + 
          facet_grid(. ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace))
        
      }
      
    } # Close loop
    
    ## Combine the list of plots and return - this is plot version 1
    combined_plot1 <- plot_grid(plotlist = g_list, nrow = 1)
    
    
    
    ## Plot version 2 is simply a grid without the complicated legend
    ## Instead, add 3 summaries: total accuracy, mean accuracy, accuracy range.
    # First create a annotation df
    ann_df <- select(df1, trait, model, pop, contains("ability_"), bias_all) %>% 
      distinct() %>%
      mutate(
        # prefix = ifelse(statistic == "ability_all", "r[MP]==", "Bias[M]=="),
        # annotation = paste0(prefix, formatC(x = value, digits = 2, format = "fg")),
        annotation = paste0("atop(All~r[MP]==", formatC(x = ability_all, digits = 2, format = "fg"), ", ",
                            "Mean~r[MP]==", formatC(x = ability_mean, digits = 2, format = "fg"), ", ",
                            "Range~r[MP]*':'~", formatC(x = ability_min, digits = 2, format = "fg"), "*','~", 
                            formatC(x = ability_max, digits = 2, format = "fg"), ")"))
    
    # Next, plot
    combined_plot2 <- df1 %>%
      ggplot(aes(x = prediction, y = value, color = annotation, group = 1)) +
      geom_point(size = 0.5) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), 
                color = "black", parse = TRUE, size = 1.5, vjust = -0.1, hjust = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "black", lwd = 0.5) +
      scale_y_continuous(name = "Observation", breaks = pretty, labels = scale_trunc) +
      scale_x_continuous(name = "Prediction", breaks = pretty) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(trait ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace)) +
      theme_presentation2(10) +
      theme(axis.title = element_blank())
    
    ## Return both plots
    tibble(plot_name = c("plot1", "plot2"), plot = list(combined_plot1, combined_plot2))
    
    }) %>% ungroup()


### Combine plots (type 1) ###



### Combine plots (type 2) ###
g_loeo_predictions_list2 <- g_loeo_prediction_list %>%
  filter(plot_name == "plot2") %>%
  pull(plot) %>%
  # Remove strip names from all but first element
  modify_at(.x = ., .at = -1, .f = ~. + theme(strip.text.x = element_blank()))

# Use cowplot
g_loeo_predictions2 <- plot_grid(plotlist = g_loeo_predictions_list2, align = "v", ncol = 1, 
                                 rel_heights = c(1, rep(0.8, length(g_loeo_predictions_list2) - 1)))

# Save
ggsave(filename = "loeo_model_predictions_observations.jpg", plot = g_loeo_predictions2,
       path = fig_dir, width = 9, height = 9, dpi = 1000)







## Visualize
g_loeo_predictions_summ <- loeo_predictions_accuracy_summ %>%
  ggplot(aes(x = pop, color = model)) +
  geom_linerange(aes(ymin = ability_min, ymax = ability_max, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  geom_point(aes(y = ability_all, shape = "all"), position = position_dodge(0.9)) +
  geom_point(aes(y = ability_mean, shape = "mean"), position = position_dodge(0.9)) +
  scale_shape_discrete(name = "Statistic", labels = str_to_title) + 
  scale_color_discrete(name = "Model", labels = f_model_replace) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Population", labels = str_to_upper) +
  facet_grid(~ trait) +
  theme_presentation2(10)

# Save
ggsave(filename = "loeo_model_predictions_summary.jpg", plot = g_loeo_predictions_summ,
       path = fig_dir, width = 6, height = 3, dpi = 1000)




### Plot accuracy
loeo_predictions_accuracy %>%
  ggplot(aes(x = trait, y = accuracy, color = model)) +
  geom_boxplot(position = position_dodge(0.9)) +
  facet_grid(~ pop) +
  theme_presentation(16)


## Calculate average and 1 sd
loeo_predictions_accuracy_summ <- loeo_predictions_accuracy %>%
  group_by(trait, pop, model) %>%
  summarize_at(vars(accuracy, ability, bias), list(~mean, ~sd, lower = ~mean(.) - sd(.), upper = ~mean(.) + sd(.))) %>%
  ungroup()


# Function to add truncation to axis text
axis_mod <- function(x, w) str_pad(string = x, width = w, side = "left", pad = " ")

## Plot accuracy versus bias
g_loeo_predictions_accuracy_list <- loeo_predictions_accuracy_summ %>%
  # Determine the max width of axis variables
  mutate(x_width = max(nchar(pretty(accuracy_lower))),
         y_width = max(nchar(pretty(bias_lower)))) %>%
  split(.$trait) %>%
  map(~{
    ggplot(data = .x, aes(x = accuracy_mean, y = bias_mean, xmin = accuracy_lower, xmax = accuracy_upper, 
               ymin = bias_lower, ymax = bias_upper)) +
      geom_errorbar(color = "grey85") + 
      geom_errorbarh(color = "grey85") +
      geom_point(aes(color = model, shape = pop), size = 2) +
      scale_x_continuous(breaks = pretty) + 
      scale_y_continuous(breaks = pretty, labels = function(x) str_pad(string = x, width = unique(.x$y_width), 
                                                                       side = "left", pad = " ")) +
      facet_grid(~ trait, scales = "free") +
      theme_presentation2(base_size = 10) +
      theme(axis.title = element_blank(), panel.grid.minor = element_blank())
  })


## Combine plots
g_loeo_predictions_accuracy_summ <- plot_grid(plotlist = map(g_loeo_predictions_accuracy_list, ~.+theme(legend.position = "none")),
                                              nrow = 1)
# Save
ggsave(filename = "loeo_accuracy_bias_combined.jpg", plot = g_loeo_predictions_accuracy_summ, path = fig_dir,
       height = 2, width = 8.5, dpi = 1000)
                                              




###################################
### Leave-one-year-out
###################################

## Add pop information
loyo_predictions_df <- loyo_predictions_out %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp"))



## Calculate accuracy and bias per environment, model, and population
## Also calculate accuracy and bias across all environments for a model/population
loyo_predictions_ability <- loyo_predictions_df %>%
  group_by(trait, environment, model, pop) %>%
  summarize(ability = cor(prediction, value), bias = mean(prediction - value, na.rm = T)) %>%
  bind_rows(., mutate(loyo_predictions_df, environment = "All") %>% group_by(trait, environment, model, pop) %>%
              summarize(ability = cor(prediction, value, use = "complete.obs"), bias = mean(prediction - value, na.rm = T))) %>%
  ungroup()

## Adjust ability using heritability
loyo_predictions_accuracy <- loyo_predictions_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))

## Create a table with trait, model, population, total accuracy and range of accuracy within environment
loyo_predictions_accuracy_summ <- loyo_predictions_accuracy %>%
  group_by(trait, model, pop) %>%
  summarize_at(vars(ability, bias), list(all = ~.[environment == "All"], mean = ~mean(.[environment != "All"], na.rm = T),
                                         min = ~min(.[environment != "All"], na.rm = T), max = ~max(.[environment != "All"], na.rm = T))) %>%
  ungroup()



## Plot predicted and observed value, with accuracy for each environment
## and overall

# First create annotation for accuracy
loyo_predictions_df1 <- loyo_predictions_accuracy %>%
  mutate(annotation = paste0(environment, "~(r[MP]==", formatC(x = ability, 2, format = "g"), ")")) %>%
  left_join(loyo_predictions_df, .)


# Plot predicted versus observed value
g_loyo_prediction_list <- loyo_predictions_df1 %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait) %>%
  do({
    df <- .
    max_wid <- unique(df$max_nchar)
    # Scale truncation function
    scale_trunc <- function(x) str_pad(string = x, width = max_wid, side = "left", pad = " ")
    
    # Create a separate legend df
    x1 <- mutate(df, x = ifelse(model == "model1", 0.20, 0.80), y = 0.10)
    
    ## Split x by model and population
    x_split <- split(x1, list(x1$pop, x1$model))
    
    # Map over this split
    g_list <- vector("list", length = length(x_split))
    
    for (i in seq_along(g_list)) {
      x2 <- x_split[[i]]
      
      g_x <- ggplot(data = x2, aes(x = prediction, y = value, color = annotation)) +
        geom_point(size = 0.5) +
        scale_y_continuous(name = "Observation", breaks = pretty, labels = scale_trunc) +
        scale_x_continuous(name = "Prediction", breaks = pretty) +
        scale_color_discrete(name = NULL, labels = function(x) parse(text = x), guide = guide_legend(ncol = 1)) +
        theme_presentation2(10) +
        theme(legend.position = c(x2$x[1], x2$y[1]), legend.key.height = unit(0.5, "line"), legend.key.width = unit(0.5, "line"),
              legend.text = element_text(size = 6),
              # No axis titles
              axis.title = element_blank())
      
      ## Add faceting depending on position
      if (i == 1) {
        g_list[[i]] <- g_x + 
          facet_grid(trait ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace))
      } else {
        g_list[[i]] <- g_x + 
          facet_grid(. ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace))
        
      }
      
    } # Close loop
    
    ## Combine the list of plots and return - this is plot version 1
    combined_plot1 <- plot_grid(plotlist = g_list, nrow = 1)
    
    ## Plot version 2 is simply a grid without the complicated legend
    ## Instead, add 3 summaries: total accuracy, mean accuracy, accuracy range.
    # First create a annotation df
    ann_df <- select(df1, trait, model, pop, contains("ability_"), bias_all) %>% 
      distinct() %>%
      mutate(
        # prefix = ifelse(statistic == "ability_all", "r[MP]==", "Bias[M]=="),
        # annotation = paste0(prefix, formatC(x = value, digits = 2, format = "fg")),
        annotation = paste0("atop(All~r[MP]==", formatC(x = ability_all, digits = 2, format = "fg"), ", ",
                            "Mean~r[MP]==", formatC(x = ability_mean, digits = 2, format = "fg"), ", ",
                            "Range~r[MP]*':'~", formatC(x = ability_min, digits = 2, format = "fg"), "*','~", 
                            formatC(x = ability_max, digits = 2, format = "fg"), ")"))
    
    
    combined_plot2 <- df %>%
      ggplot(aes(x = prediction, y = value, color = annotation, group = 1)) +
      geom_point(size = 0.5) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), 
                color = "black", parse = TRUE, size = 3, vjust = -0.5, hjust = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "black", lwd = 0.5) +
      scale_y_continuous(name = "Observation", breaks = pretty, labels = scale_trunc) +
      scale_x_continuous(name = "Prediction", breaks = pretty) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(trait ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace)) +
      theme_presentation2(10) +
      theme(axis.title = element_blank())
    
    ## Return both plots
    tibble(plot_name = c("plot1", "plot2"), plot = list(combined_plot1, combined_plot2))
    
  }) %>% ungroup()




### Combine plots (type 1) ###



### Combine plots (type 2) ###
g_loyo_predictions_list2 <- g_loyo_prediction_list %>%
  filter(plot_name == "plot2") %>%
  pull(plot) %>%
  # Remove strip names from all but first element
  modify_at(.x = ., .at = -1, .f = ~. + theme(strip.text.x = element_blank()))

# Use cowplot
g_loyo_predictions2 <- plot_grid(plotlist = g_loyo_predictions_list2, align = "v", ncol = 1, 
                                 rel_heights = c(1, rep(0.8, length(g_loyo_predictions_list2) - 1)))

# Save
ggsave(filename = "loyo_model_predictions_observations.jpg", plot = g_loyo_predictions2,
       path = fig_dir, width = 9, height = 9, dpi = 1000)







## Visualize
g_loyo_predictions_summ <- loyo_predictions_accuracy_summ %>%
  ggplot(aes(x = pop, color = model)) +
  geom_linerange(aes(ymin = ability_min, ymax = ability_max, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  geom_point(aes(y = ability_all, shape = "all"), position = position_dodge(0.9)) +
  geom_point(aes(y = ability_mean, shape = "mean"), position = position_dodge(0.9)) +
  scale_shape_discrete(name = "Statistic", labels = str_to_title) + 
  scale_color_discrete(name = "Model", labels = f_model_replace) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Population", labels = str_to_upper) +
  facet_grid(~ trait) +
  theme_presentation2(10)

# Save
ggsave(filename = "loyo_model_predictions_summary.jpg", plot = g_loyo_predictions_summ,
       path = fig_dir, width = 6, height = 3, dpi = 1000)










##### Compare LOEO and LOYO #####

loo_predictions_accuracy_summ <- bind_rows(
  mutate(loeo_predictions_accuracy_summ, scheme = "LOEO"), 
  mutate(loyo_predictions_accuracy_summ, scheme = "LOYO")
)

## Visualize
g_loo_predictions_summ <- loo_predictions_accuracy_summ %>%
  ggplot(aes(x = scheme, color = model)) +
  geom_linerange(aes(ymin = ability_min, ymax = ability_max, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  geom_point(aes(y = ability_all, shape = "all"), position = position_dodge(0.9)) +
  geom_point(aes(y = ability_mean, shape = "mean"), position = position_dodge(0.9)) +
  scale_shape_discrete(name = "Statistic", labels = str_to_title) + 
  scale_color_discrete(name = "Model", labels = f_model_replace) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Population", labels = str_to_upper) +
  facet_grid(trait ~ pop, switch = "y", labeller = labeller(trait = str_add_space, pop = toupper)) +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary.jpg", plot = g_loo_predictions_summ,
       path = fig_dir, width = 4, height = 5, dpi = 1000)




#






























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






