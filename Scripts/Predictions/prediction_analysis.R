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


## Filter loeo results
loeo_predictions_out <- loeo_predictions_out %>% 
  subset(., map_lgl(., is_tibble)) %>%
  bind_rows()


## Significant level
alpha <- 0.05


# model1: y = G + r
# model2: y = G + E + r
# model3: y = G + E + GE + r
# model4: y = G + L + r
# model5: y = G + L + GL + r


# Color scheme for models
model_colors <- paletteer_d(package = "ggsci", palette = "default_nejm", n = length(model_replace)) %>%
  setNames(., names(model_replace))





# Cross-validation results ------------------------------------------------




## Grab the prediction outputs - cross-validation
loo_prediction_list <- map(set_names(object_list, object_list), get) %>%
  subset(., str_detect(names(.), "^lo[a-z]o")) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_extract(names(.), "^[a-z]{4}"))


## Combine data.frames and mutate columns
loo_predictions_df <- loo_prediction_list %>%
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
  mutate(leave_one_group = site1, nGroup = nSite) %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1", "nSite", "site1", "source")))



## Calculate accuracy and bias per train group, model, and population
loo_predictive_ability <- loo_predictions_df %>%
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  summarize_at(vars(ability, bias, ability_all, bias_all, nObs, nGroup), mean) %>%
  ungroup()


# ## Adjust ability using heritability
# loo_prediction_accuracy <- loo_predictive_ability %>%
#   left_join(., env_trait_herit, by = c("leave_one_group" = "environment")) %>%
#   mutate(accuracy = ability / sqrt(heritability))



## Quick plot of accuracy across all data points
loo_predictive_ability %>%
  distinct(trait, model, pop, type, selection, ability_all, bias_all) %>%
  ggplot(aes(x = model, y = ability_all, fill = selection)) +
  geom_col(position = position_dodge(0.9)) +
  facet_grid(trait ~ type + pop)

loo_predictions_df %>%
  filter(trait == "GrainYield", type == "lolo", selection != "none") %>%
  filter(selection != "concurrent_rfa_cv_adhoc") %>%
  ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 2) +
  scale_color_discrete(guide = FALSE) +
  facet_grid(type + selection ~ model + pop) +
  theme_presentation2(10)






## Plot predicted and observed value, with mean and range of accuracy per
## environment and accuracy overall

# Calculate predictive ability over all observations; bootstrap to get a confidence
# interval; also calculate RMSE
loo_accuracy_rmse_all <- loo_predictions_df %>% 
  group_by(trait, type, model, pop, selection) %>% 
  do(tibble(cor_test = list(cor.test(x = .$value, y = .$pred_complete)),
            rmse = sqrt(mean((.$value - .$pred_complete)^2)))) %>%
  ungroup() %>%
  mutate(ability = map_dbl(cor_test, "estimate"),
         ci = map(cor_test, "conf.int"),
         lower = map_dbl(ci, 1), upper = map_dbl(ci, 2)) %>%
  select(-ci)

# First create annotation df
loo_prediction_accuracy_annotation <- loo_accuracy_rmse_all %>%
  rename(ability_all = ability, ability_lower = lower, ability_upper = upper) %>%
  mutate_at(vars(contains("ability")), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_all_annotation = paste0("r[MP]==", ability_all, "~(", ability_lower, "*','~", ability_upper, ")"))


# Plot predicted versus observed value
g_loo_prediction_list <- loo_predictions_df %>%
  filter(model %in% names(model_replace)) %>%
  filter(selection != "concurrent_rfa_cv_adhoc") %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait, type) %>%
  do({
    df <- .
    
    # Convert yield to t ha^-1
    if (unique(df$trait) == "GrainYield") {
      df1 <- mutate_at(df, vars(value, pred_complete), ~. / 1000)
    } else {
      df1 <- df
    }
    
    breaks <- map(subset(df1, pred_complete > 0, c(pred_complete, value)), pretty) 
    
    
    ## Extract the appropriate accuracy annotation
    r_mp_annotation <- left_join(distinct(df1, trait, type, model, pop, selection), loo_prediction_accuracy_annotation,
                                 by = c("trait", "type", "model", "pop", "selection")) %>%
      mutate(annotation = ability_all_annotation,
             x = breaks$pred_complete[1], y = c(last(breaks$value)))
    
    df1_none <- filter(df1, selection == "none")
    
    ## First plot selection == none (identity covariance for environments and gxe)
    g_plot_none <- df1_none %>%
      ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = left_join(distinct(df1_none, trait, type, model, selection), r_mp_annotation), 
                aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE, hjust = 0, size = 2) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(selection + pop ~ model, switch = "y", 
                 labeller = labeller(selection = f_ec_selection_replace, pop = f_pop_replace, model = f_model_replace)) +
      labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
      theme_presentation2(10)
    
    # Filter the df
    df1_rest <- filter(df1, selection != "none")
    
    ## Create the plot of each feature selection procedure
    g_plot <- df1_rest %>%
      ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = left_join(distinct(df1_rest, trait, type, model, selection), r_mp_annotation), 
                aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE, hjust = 0, size = 2) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(selection + pop ~ model, switch = "y", 
                 labeller = labeller(selection = f_ec_selection_replace, pop = f_pop_replace, model = f_model_replace)) +
      labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
      theme_presentation2(10)
    
    ## Output a list of plots
    tibble(base_model_plot = list(g_plot_none),
           full_covariate_selection_plot = list(g_plot))
    
    }) %>% ungroup()



## Combine base models and all covariate selection
for (i in seq_len(nrow(g_loo_prediction_list))) {
  
  # Combine the plots
  plot_comb <- plot_grid(g_loo_prediction_list$base_model_plot[[i]], g_loo_prediction_list$full_covariate_selection_plot[[i]],
                         rel_heights = c(1/3, 1), ncol = 1)
  
  # Save
  filename <- paste0("loo_model_predictions_observations_", g_loo_prediction_list$trait[i], 
                     "_", g_loo_prediction_list$type[i], ".jpg")
  
  ggsave(filename = filename, plot = plot_comb, path = fig_dir, width = 6, height = 16, dpi = 1000)
  
}


## Just plot the reduced set of feature selection
for (i in seq_len(nrow(g_loo_prediction_list))) {
  
  # Combine the plots
  plot_comb <- g_loo_prediction_list$full_covariate_selection_plot[[i]] %>%
    modify_at("data", ~filter(., selection != "stepAIC_adhoc"))
  plot_comb$layers[[3]]$data <- filter(plot_comb$layers[[3]]$data, selection != "stepAIC_adhoc")
  
  # Save
  filename <- paste0("loo_covariate_models_predictions_observations_", g_loo_prediction_list$trait[i], 
                     "_", g_loo_prediction_list$type[i], ".jpg")
  
  ggsave(filename = filename, plot = plot_comb, path = fig_dir, width = 6, height = 10, dpi = 1000)
  
}




## Barplot of overall prediction accuracy
g_loo_predictions_all_summ <- loo_accuracy_bias_all %>%
  mutate(selection = f_ec_selection_replace(selection),
         selection = fct_relevel(selection, "None")) %>% # Put "None" first
  filter(measure == "ability") %>%
  ggplot(aes(x = selection, group = model)) +
  geom_hline(yintercept = 0, color = "grey85") +
  geom_col(aes(y = base, fill = model), position = "dodge", color = "black", lwd = 0.1) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(0.9),
                width = 0.5, color = "black") + 
  # scale_fill_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_fill_paletteer_d(package = "dutchmasters", palette = "milkmaid",
                         name = "Model", labels = f_model_replace) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  # scale_x_discrete(name = "Validation scheme", labels = f_pop_replace) +
  scale_x_discrete(name = "Covariate selection", labels = function(x) parse(text = x)) +
  facet_grid(type ~ trait + pop, labeller = labeller(trait = str_add_space, type = toupper, pop = f_pop_replace),
             switch = "y") +
  theme_presentation2(10) +
  theme(strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave(filename = "loo_model_predictions_all_accuracy.jpg", plot = g_loo_predictions_all_summ,
       path = fig_dir, width = 14, height = 6, dpi = 1000)






### Summarize accuracy within each environment ###


## Create a function list for calculating quantiles
q <- c(alpha / 2, 0.5, 1 - (alpha / 2))
# Function vector
q_funs <- c(map(q, ~partial(quantile, probs = .x, na.rm = TRUE)), partial(mean, na.rm = TRUE)) %>%
  set_names(., c("lower", "median", "upper", "mean"))


## Summarize mean and range
loo_prediction_accuracy_summ <- loo_predictive_ability %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, model, pop, selection) %>%
  summarize_at(vars(ability, bias), funs(!!!q_funs)) %>%
  ungroup()

## Determine the LSD between model-selection groups
loo_prediction_accuracy_summ1 <- loo_predictive_ability %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, pop) %>%
  do({
    dat <- .
    fit <- lm(ability ~ model + selection + model:selection, data = dat)
    n <- max(as.numeric(xtabs(~ model + selection, dat)))
    MS_within <- sigma(fit)^2 # Within-group error (residuals)
    df <- df.residual(fit)
    LSD <- qt(p = 1 - (0.05 / 2), df = df) * sqrt(MS_within * (2 / n))
    
    # Return mean and LSD of prediction accuracy
    aggregate(ability ~ model + selection, dat, mean) %>% 
      mutate(LSD = LSD)
  }) %>% ungroup()

annotation_df <- loo_prediction_accuracy_summ1 %>%
  distinct(type, trait, pop, LSD) %>%
  mutate(annotation = paste0("LSD: ", round(LSD, 2)))

## Visualize
# Plot mean and range of predictions
g_loo_predictions_summ <- loo_prediction_accuracy_summ1 %>%
  filter(model %in% names(model_replace)) %>%
  mutate(selection = f_ec_selection_replace(selection),
         selection = fct_relevel(selection, "None")) %>% # Put "None" first
  ggplot(aes(x = selection, group = model, color = model)) +
  # geom_linerange(aes(ymin = ability_lower, ymax = ability_upper, group = model), 
  #                position = position_dodge(0.9), color = "grey85") +
  geom_point(aes(y = ability), position = position_dodge(0.9)) +
  geom_text(data = annotation_df, aes(x = 1, y = 1, label = annotation), size = 2, hjust = 0,
            inherit.aes = FALSE) +
  geom_segment(data = annotation_df, aes(x = 3, xend = 3, y = 1, yend = 1 - LSD), inherit.aes = FALSE) +
  scale_color_paletteer_d(package = "dutchmasters", palette = "milkmaid",
                          name = "Model", labels = f_model_replace) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  # scale_x_discrete(name = "Validation scheme", labels = f_pop_replace) +
  scale_x_discrete(name = "Covariate selection", labels = function(x) parse(text = x)) +
  facet_grid(type ~ trait + pop, labeller = labeller(trait = str_add_space, type = toupper, pop = f_pop_replace),
             switch = "y") +
  theme_presentation2(10) +
  theme(strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave(filename = "loo_model_predictions_summary.jpg", plot = g_loo_predictions_summ,
       path = fig_dir, width = 12, height = 6, dpi = 1000)


## Create a table

# Report the mean (and range) in predictive ability across environments for each
# validation scheme, population, and type
loo_prediction_accuracy_table <- loo_prediction_accuracy_summ %>% 
  mutate_at(vars(contains("ability")), ~formatC(., digits = 2, width = 2, format = "g")) %>%
  mutate(annotation = paste0(ability_mean, " (", ability_lower, ", ", ability_upper, ")")) %>%
  # mutate(annotation = paste0(ability_wmean, " (", ability_min, ", ", ability_max, ")")) %>%
  # Rename
  mutate(model = f_model_replace(model),
         pop = f_pop_replace(pop),
         type = toupper(type)) %>%
  select(trait, type, model, pop, annotation) %>%
  spread(model, annotation) %>%
  arrange(trait, pop, type)

write_csv(x = loo_prediction_accuracy_table, path = file.path(fig_dir, "loo_prediction_accuracy_table.csv"))

## Report accuracy for POV
pov_loo_prediction_accuracy_table <- loo_prediction_accuracy_table %>%
  filter(pop == "POV00") %>%
  select(trait, type, unname(model_present)) %>%
  mutate(type = f_type_replace(tolower(type)))
write_csv(x = pov_loo_prediction_accuracy_table, path = file.path(fig_dir, "pov_loo_prediction_accuracy_table.csv"))










# External validation -----------------------------------------------------




## Grab the prediction outputs - cross-validation
external_prediction_list <- map(set_names(object_list, object_list), get) %>%
  subset(., str_detect(names(.), "external")) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_extract(names(.), "env|loc"))


## Combine data.frames and mutate columns
external_predictions_df <- external_prediction_list %>%
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
  mutate(leave_one_group = site,
         nGroup = nSite) %>%
  # Remove charlottetown
  # filter(leave_one_group != "Charlottetown") %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1", "nSite", "site1")))



## Calculate accuracy and bias per train group, model, and population
external_predictive_ability <- external_predictions_df %>%
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  summarize_at(vars(ability, bias, ability_all, bias_all, nObs, nGroup), mean) %>%
  ungroup()


## Quick plot of accuracy across all data points
external_predictive_ability %>%
  distinct(trait, model, pop, type, selection, ability_all, bias_all) %>%
  ggplot(aes(x = model, y = ability_all, fill = selection)) +
  geom_col(position = position_dodge(0.9)) +
  facet_grid(trait ~ type + pop, scales = "free_y") +
  theme_genetics()

external_predictions_df %>%
  filter(trait == "PlantHeight", selection != "none", type == "loc") %>%
  ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  # scale_color_discrete(guide = FALSE) +
  facet_grid(type + selection ~ model + pop) +
  theme_presentation2(10)



## Plot predicted and observed value, with mean and range of accuracy per
## environment and accuracy overall

# Calculate predictive ability over all observations; bootstrap to get a confidence
# interval; also calculate bias
external_accuracy_bias_all <- external_predictions_df %>% 
  group_by(trait, type, model, pop, selection) %>% 
  do({
    tibble(measure = c("ability", "bias"),
           out = list(neyhart::bootstrap(x = .$value, y = .$pred_complete, fun = "cor", boot.reps = 1000),
                      neyhart::bootstrap(x = .$pred_complete, y = .$value, fun = "bias", boot.reps = 1000)))
  }) %>% ungroup() %>%
  unnest() %>%
  select(-statistic)


# First create annotation df
external_prediction_accuracy_annotation <- external_accuracy_bias_all %>%
  filter(measure == "ability") %>%
  rename(ability_all = base, ability_lower = ci_lower, ability_upper = ci_upper) %>%
  mutate_at(vars(contains("ability")), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_all_annotation = paste0("r[MP]==", ability_all, "~(", ability_lower, "*','~", ability_upper, ")"))


# Plot predicted versus observed value
g_external_prediction_list <- external_predictions_df %>%
  filter(model %in% names(model_replace)) %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait, type) %>%
  do(plot = {
    df <- .
    
    # Convert yield to t ha^-1
    if (unique(df$trait) == "GrainYield") {
      df1 <- mutate_at(df, vars(value, pred_complete), ~. / 1000)
    } else {
      df1 <- df
    }
    
    breaks <- map(subset(df1, pred_complete > 0, c(pred_complete, value)), pretty) 
    
    
    ## Extract the appropriate accuracy annotation
    r_mp_annotation <- left_join(distinct(df1, trait, type, model, pop, selection), loo_prediction_accuracy_annotation,
                                 by = c("trait", "type", "model", "pop", "selection")) %>%
      mutate(annotation = ability_all_annotation,
             x = breaks$pred_complete[1], y = c(last(breaks$value)),
             # Edit the covariate selection variable
             selection = paste0("Covariates:~", f_ec_selection_replace(selection)))
    
    ## First plot selection == none (identity covariance for environments and gxe)
    g_plot_none <- df1 %>%
      filter(selection == "none") %>%
      # Edit the covariate selection variable
      mutate(selection = paste0("Covariates:~", f_ec_selection_replace(selection))) %>%
      ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = filter(r_mp_annotation, str_detect(selection, "None")), 
                aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE, hjust = 0, size = 2) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(selection + pop ~ model, switch = "y", 
                 labeller = labeller(selection = label_parsed, pop = f_pop_replace, model = f_model_replace)) +
      labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
      theme_presentation2(10)
    
    ## Create the plot
    g_plot <- df1 %>%
      filter(selection != "none") %>%
      # Edit the covariate selection variable
      mutate(selection = paste0("Covariates:~", f_ec_selection_replace(selection))) %>%
      ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = filter(r_mp_annotation, str_detect(selection, "None", negate = TRUE)), 
                aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE, hjust = 0, size = 2) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(selection + pop ~ model, switch = "y", 
                 labeller = labeller(selection = label_parsed, pop = f_pop_replace, model = f_model_replace)) +
      labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
      theme_presentation2(10)
    
    ## Combine plots and return
    plot_grid(g_plot_none, g_plot, ncol = 1, rel_heights = c(1/3, 1))
    
  }) %>% ungroup()


for (i in seq_len(nrow(g_external_prediction_list))) {
  
  # Save
  filename <- paste0("external_model_predictions_observations_", 
                     paste0(map_chr(select(g_external_prediction_list, -plot), i), collapse = "_"), 
                     ".jpg")
  ggsave(filename = filename, plot = g_external_prediction_list$plot[[i]], 
         path = fig_dir, width = 6, height = 8, dpi = 1000)
  
}



#





































# Analyze output of factorial regression samples --------------------------

## Predictions from this exercise were genuinely leave-one-out.
## A location or environment was removed, covariates were selected, and
## predictions were made. These results should be compared to the previous
## leave-one-out strategy

# Load the data
load(file.path(result_dir, "loo_predictions_fact_reg_samples.RData"))

# Combine the data into list
## Grab the prediction outputs
loo_samples_prediction_list <- ls(pattern = "sample_out") %>%
  setNames(., .) %>%
  map(get) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_extract(names(.), "^[a-z]{4}"))


## Combine data.frames and mutate columns
loo_samples_prediction_df <- loo_samples_prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~mutate_if(., is.character, parse_guess) %>%
           mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  select(-which(names(.) %in% c(".id", "core", "trait1"))) %>%
  # Coalesce columns
  mutate(leave_one_group = coalesce(environment, location),
         nGroup = coalesce(nEnv, nLoc)) %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1")))


## Calculate accuracy and bias per train group, model, and population
loo_samples_predictive_ability <- loo_samples_prediction_df %>%
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  summarize_at(vars(ability, bias, ability_all, bias_all, nObs, nGroup), mean) %>%
  ungroup()

## Annotation df
loo_samples_predictive_ability_annotation <- loo_samples_predictive_ability %>%
  distinct(trait, model, pop, type, selection, ability_all) %>%
  mutate(annotation = paste0("r[MP]==", round(ability_all,  2)))


## Plot predicted versus observed phenotypes
g_loo_samples_prediction <- loo_samples_prediction_df %>%
  # Split by type and pop
  group_by(type, pop) %>%
  do(plot = {
    df <- .
    
    # Group by trait
    trait_plots <- df %>%
      group_by(trait) %>%
      do(plot = {
        df1 <- .
        
        ## Subset the annotation df
        df1_annotation <- left_join(x = distinct(df1, trait, selection, model, type, pop), loo_samples_predictive_ability_annotation)
        
        # Plot
        ggplot(data = df1, aes(x = pred_complete, y = value, color = leave_one_group)) +
          geom_point(size = 0.5) +
          geom_text(data = df1_annotation, aes(x = Inf, y = -Inf, label = annotation), 
                    parse = TRUE, inherit.aes = FALSE, vjust = -1.2, hjust = 1, size = 2) +
          facet_grid(trait ~ model + selection, switch = "y") +
          scale_color_discrete(guide = FALSE) +
          scale_x_continuous(breaks = pretty) +
          scale_y_continuous(breaks = pretty) +
          theme_genetics(10)
        
      }) %>% ungroup()
    
    # Combine the plots
    plot_grid(plotlist = trait_plots$plot, ncol = 1, align = "hv")
    
  })


## Visualize
# LOEO - TP
g_loo_samples_prediction$plot[[1]]

# LOEO - VP
g_loo_samples_prediction$plot[[2]]

# LOLO - TP
g_loo_samples_prediction$plot[[3]]

# LOLO - VP
g_loo_samples_prediction$plot[[4]]













