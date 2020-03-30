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


## Significant level
alpha <- 0.05


# model1: y = G + r
# model2: y = G + E + r
# model3: y = G + E + GE + r
# model4: y = G + L + r
# model5: y = G + L + GL + r


## A vector to rename models
model_replace <- c("model1" = "G", "model2" = "G + E", "model2a" = "G + E (AMMI)", 
                   "model3" = "G + E + GE", "model3a" = "G + E + GE (AMMI)",
                   "model4" = "G + L", "model4a" = "G + L (AMMI)",
                   "model5" = "G + L + GL", "model5a" = "G + L + GL (AMMI)")


f_model_replace <- function(x) model_replace[x]
# f_model_replace <- function(x) paste0("M", toupper(str_extract(x, "[0-9]{1}[a-z]{0,1}")))
# Vector to rename validation schemes
f_pop_replace <- function(x) str_replace_all(x, c("tp" = "CV0", "vp" = "POV00"))

# Color scheme for models
model_colors <- c(neyhart_palette("umn1")[1], neyhart_palette("umn3")[3], neyhart_palette("umn1")[3],
                  neyhart_palette("umn2")[3], neyhart_palette("umn3")[4], neyhart_palette("umn1")[4],
                  neyhart_palette("umn2")[4])


# Phenotye data to use
data_to_model <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.))) %>%
  droplevels() %>%
  mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)),
         # Collapse Ithacas
         location = str_remove_all(location, "[0-9]{1}")) %>%
  mutate_at(vars(environment, location, year), as.factor) %>%
  mutate(env = environment)



## Grab the prediction outputs
loo_prediction_list <- map(set_names(object_list, object_list), get) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  # Select relevant columns
  map(~select(., c("trait", "model", "prediction", "nEnv", "nObs"))) %>%
  set_names(x = ., nm = str_extract(names(.), "^[a-z]{4}"))


loo_predictions_df <- loo_prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~select(., trait, environment, location, year, type, model, line_name, value, 
                 contains("nEnv"), contains("nObs"), contains("pred")) %>%
        mutate_if(is.character, parse_guess) %>%
        mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp"))


## Calculate accuracy and bias per environment, model, and population
loo_predictive_ability <- loo_predictions_df %>%
  group_by(trait, model, pop, type, environment) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, environment) %>%
  summarize_at(vars(ability, bias, ability_all, bias_all, nObs, nEnv), mean) %>%
  ungroup()


## Adjust ability using heritability
loo_prediction_accuracy <- loo_predictive_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))



## Plot predicted and observed value, with mean and range of accuracy per
## environment and accuracy overall

# Calculate prediction accuracy over all observations; bootstrap to get a confidence
# interval
loo_predictive_ability_all <- loo_predictions_df %>% 
  group_by(trait, type, model, pop) %>% 
  do(neyhart::bootstrap(x = .$value, y = .$pred_complete, fun = "cor", boot.reps = 1000)) %>%
  ungroup()


# First create annotation df
loo_prediction_accuracy_annotation <- loo_predictive_ability_all %>%
  rename(ability_all = base, ability_lower = ci_lower, ability_upper = ci_upper) %>%
  mutate_at(vars(contains("ability")), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_all_annotation = paste0("r[MP]==", ability_all, "~(", ability_lower, "*','~", ability_upper, ")"))


# Plot predicted versus observed value
g_loo_prediction_list <- loo_predictions_df %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait) %>%
  do(plot = {
    df <- .
    
    df1 <- df
    
    breaks <- map(select(df1, pred_complete, value), pretty) 
    
    
    ## Extract the appropriate accuracy annotation
    r_mp_annotation <- left_join(distinct(df1, trait, type, model, pop), loo_prediction_accuracy_annotation,
                                 by = c("trait", "type", "model", "pop")) %>%
      select(trait:pop, contains("annotation")) %>% 
      mutate(annotation = ability_all_annotation) %>%
      mutate(x = breaks$pred_complete[1], y = c(last(breaks$value)))
    
    ## Create the plot
    g_plot <- df1 %>%
      ggplot(aes(x = pred_complete, y = value, color = environment)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = r_mp_annotation, aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE,
                hjust = 0, size = 2) +
      scale_y_continuous(name = "Observation", breaks = pretty, labels = scale_trunc) +
      scale_x_continuous(name = "Prediction", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(type + pop ~ model, switch = "y", 
                 labeller = labeller(type = toupper, pop = f_pop_replace, model = f_model_replace)) +
      labs(subtitle = str_add_space(unique(df1$trait))) +
      theme_presentation2(10)
    
    ## Return the plot
    g_plot
    
    }) %>% ungroup()


for (i in seq_len(nrow(g_loo_prediction_list))) {
  
  # Save
  ggsave(filename = paste0("loo_model_predictions_observations_", g_loo_prediction_list$trait[[i]], ".jpg"), 
         plot = g_loo_prediction_list$plot[[i]], path = fig_dir, width = 15, height = 10, dpi = 1000)
  
}




## Barplot of overall prediction accuracy
g_loo_predictions_all_summ <- loo_predictive_ability_all %>%
  ggplot(aes(x = pop, group = model)) +
  geom_col(aes(y = base, fill = model), position = "dodge") + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(0.9),
                width = 0.5, color = "black") + 
  # scale_fill_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_fill_paletteer_d(package = "dutchmasters", palette = "milkmaid",
                         name = "Model", labels = f_model_replace) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Validation scheme", labels = f_pop_replace) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_all_summary.jpg", plot = g_loo_predictions_all_summ,
       path = fig_dir, width = 10, height = 6, dpi = 1000)





## Create a function list for calculating quantiles
q <- c(alpha / 2, 0.5, 1 - (alpha / 2))
# Function vector
q_funs <- c(map(q, ~partial(quantile, probs = .x, na.rm = TRUE)), partial(mean, na.rm = TRUE)) %>%
  set_names(., c("lower", "median", "upper", "mean"))


## Summarize mean and range
loo_prediction_accuracy_summ <- loo_prediction_accuracy %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, model, pop) %>%
  summarize_at(vars(ability, bias, accuracy), funs(!!!q_funs))
                 



## Visualize
# Plot mean and range of predictions
g_loo_predictions_summ <- loo_prediction_accuracy_summ %>%
  ggplot(aes(x = pop, color = model)) +
  geom_linerange(aes(ymin = ability_lower, ymax = ability_upper, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  # geom_linerange(aes(ymin = ability_q25, ymax = ability_q75, group = model), 
  #                position = position_dodge(0.9), color = "grey85", lwd = 1) +
  geom_point(aes(y = ability_mean), position = position_dodge(0.9)) +
  # geom_point(aes(y = ability_wmean), position = position_dodge(0.9)) +
  # scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_color_paletteer_d(package = "dutchmasters", palette = "milkmaid",
                          name = "Model", labels = f_model_replace) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Validation scheme", labels = f_pop_replace) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary.jpg", plot = g_loo_predictions_summ,
       path = fig_dir, width = 10, height = 6, dpi = 1000)



## Create a table

# Report the mean (and range) in predictive ability across environments for each
# validation scheme, population, and type
loo_prediction_accuracy_table <- loo_prediction_accuracy_summ %>% 
  mutate_at(vars(contains("ability")), ~formatC(., digits = 2, width = 2, format = "g")) %>%
  # mutate(annotation = paste0(ability_mean, " (", ability_min, ", ", ability_max, ")")) %>%
  mutate(annotation = paste0(ability_wmean, " (", ability_min, ", ", ability_max, ")")) %>%
  # Rename
  mutate(model = f_model_replace(model),
         pop = f_pop_replace(pop),
         type = toupper(type)) %>%
  select(trait, type, model, pop, annotation) %>%
  spread(model, annotation) %>%
  arrange(trait, pop, type)



# ## Calculate weighted average of accuracy using the number of environments or observations
# loo_prediction_accuracy_table <- loo_prediction_accuracy %>%
#   group_by(type, trait, model, pop) %>%
#   ## Calculate weighted average of prediction accuracy and bias using number of environments
#   summarize_at(vars(accuracy, bias), ~ weighted.mean(., w = nEnv)) %>%
#   ungroup()

write_csv(x = loo_prediction_accuracy_table, path = file.path(fig_dir, "loo_prediction_accuracy_table.csv"))



# Remove range - just report mean
loo_prediction_accuracy_table1 <- loo_prediction_accuracy_summ %>% 
  mutate_at(vars(contains("ability")), ~formatC(., digits = 2, width = 2, format = "g")) %>%
  mutate(annotation = ability_mean) %>%
  # Rename
  mutate(model = f_model_replace(model),
         pop = f_pop_replace(pop),
         type = toupper(type)) %>%
  select(trait, type, model, pop, annotation) %>%
  spread(model, annotation) %>%
  arrange(trait, pop, type)

write_csv(x = loo_prediction_accuracy_table1, path = file.path(fig_dir, "loo_prediction_accuracy_table1.csv"))





## Analyze bias

## Plot bias
g_loo_bias_summ <- loo_prediction_accuracy_summ %>%
  ggplot(aes(x = pop, color = model)) +
  geom_linerange(aes(ymin = bias_min, ymax = bias_max, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  geom_linerange(aes(ymin = bias_q25, ymax = bias_q75, group = model), 
                 position = position_dodge(0.9), color = "grey85", lwd = 1) +
  geom_point(aes(y = bias_mean), position = position_dodge(0.9)) +
  scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Validation scheme", labels = f_pop_replace) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_bias_summary.jpg", plot = g_loo_bias_summ,
       path = fig_dir, width = 8, height = 6, dpi = 1000)


## Plot accuracy versus bias
g_loo_ability_bias_each <- loo_prediction_accuracy_summ %>%
  ggplot(aes(x = bias_mean, y = ability_mean, shape = pop, color = model)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = bias_q25, xend = bias_q75, y = ability_mean, yend = ability_mean), color = "grey85", lwd = 0.5) +
  geom_segment(aes(x = bias_mean, xend = bias_mean, y = ability_q25, yend = ability_q75), color = "grey85", lwd = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_x_continuous(name = "Bias", breaks = pretty) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y", scales = "free_x") +
  theme_light(base_size = 8)

# Save
ggsave(filename = "loo_model_accuracy_bias_each.jpg", plot = g_loo_ability_bias_each,
       path = fig_dir, width = 8, height = 4.5, dpi = 1000)



# Plot summaries of accuracy and bias
g_loo_ability_bias_all <- loo_prediction_accuracy %>% 
  select(., trait, model, pop, type, contains("all")) %>% 
  distinct() %>%
  ggplot(aes(x = bias_all, y = ability_all, shape = pop, color = model)) +
  geom_vline(xintercept = 0) +
  geom_point(size = 2) +
  scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_shape_discrete(name = "Validation\nscheme", labels = f_pop_replace) +
  scale_x_continuous(name = "Bias", breaks = pretty) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y", scales = "free_x") +
  theme_light(base_size = 8)

# Save
ggsave(filename = "loo_model_accuracy_bias_all.jpg", plot = g_loo_ability_bias_all,
       path = fig_dir, width = 8, height = 4.5, dpi = 1000)



# Convert to plotly
g_loo_ability_bias_all_plotly <- plotly::ggplotly(p = g_loo_ability_bias_all)
# Save as HTML widget
htmlwidgets::saveWidget(widget = plotly::as_widget(g_loo_ability_bias_all_plotly),
                        file = file.path(fig_dir, "loo_model_accuracy_bias_all.html"))







