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


# ## A vector to rename models
# model_replace <- c("model1" = "G", "model2" = "G + E", "model2_P" = "G + E (AMMI)", 
#                    "model3" = "G + E + GE", "model3_P" = "G + E + GE (AMMI)", "model3_GxE" = "G + E + GxE",
#                    "model4" = "G + L", "model4_P" = "G + L (AMMI)",
#                    "model5" = "G + L + GL", "model5_P" = "G + L + GL (AMMI)", "model5_GxL" = "G + L + GxL")
                   
## A vector to rename models
model_replace <- c("model1" = "g", "model2_id" = "g + E", "model2_cov" = "g + w", 
                   "model3_id" = "g + w + gE", "model3_cov" = "g + w + gw")

## Models to present
# model_present <- model_replace[str_detect(names(model_replace), "_", negate = TRUE)]
model_present <- model_replace


f_model_replace <- function(x) model_replace[x]
f_model_replace2 <- function(x) model_present[x]
# f_model_replace <- function(x) paste0("M", toupper(str_extract(x, "[0-9]{1}[a-z]{0,1}")))
# Vector to rename validation schemes
f_pop_replace <- function(x) str_replace_all(x, c("tp" = "CV0", "vp" = "POV00"))
# Replace type
f_type_replace <- function(x) c("loeo" = "New environment", "lolo" = "New location", "loyo" = "New year")[x]
# Replace ec selection
f_ec_selection_replace <- function(x) 
  c("adhoc" = "italic(ad~hoc)", "adhoc_nosoil" = "italic(ad~hoc)~(no~soil)", "apriori" = "italic(a~priori)")[x]

# Color scheme for models
model_colors <- c(neyhart_palette("umn1")[1], neyhart_palette("umn3")[3], neyhart_palette("umn1")[3],
                  neyhart_palette("umn3")[4], neyhart_palette("umn1")[4])
# names(model_colors) <- grep(pattern = "_", x = names(model_replace), value = TRUE, invert = TRUE)
names(model_colors) <- names(model_replace)



## Grab the prediction outputs
loo_prediction_list <- map(set_names(object_list, object_list), get) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_extract(names(.), "^[a-z]{4}"))


loo_predictions_df <- loo_prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~mutate_if(., is.character, parse_guess) %>%
           mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  select(-which(names(.) %in% c(".id", "core", "trait1")))


## Calculate accuracy and bias per environment, model, and population
loo_predictive_ability <- loo_predictions_df %>%
  group_by(trait, model, pop, type, environment, selection) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, environment, selection) %>%
  summarize_at(vars(ability, bias, ability_all, bias_all, nObs, nEnv), mean) %>%
  ungroup()


## Adjust ability using heritability
loo_prediction_accuracy <- loo_predictive_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))



# ## Quick plot of accuracy for each trait, model, pop, type, and selection
# loo_prediction_accuracy %>%
#   group_by(trait, model, pop, type, selection) %>%
#   summarize_at(vars(accuracy, bias), mean) %>%
#   ggplot(aes(x = model, y = accuracy, fill = selection)) +
#   geom_col(position = position_dodge(0.9)) +
#   facet_grid(trait ~ type + pop)
# 
# ## Quick plot of accuracy across all data points
# loo_prediction_accuracy %>%
#   distinct(trait, model, pop, type, selection, ability_all, bias_all) %>%
#   ggplot(aes(x = model, y = ability_all, fill = selection)) +
#   geom_col(position = position_dodge(0.9)) +
#   facet_grid(trait ~ type + pop)
# 
# loo_predictions_df %>%
#   filter(trait == "TestWeight") %>%
#   ggplot(aes(x = pred_complete, y = value, color = environment)) +
#   geom_point() +
#   scale_color_discrete(guide = FALSE) +
#   facet_grid(type + selection ~ model + pop) +
#   theme_presentation2(10)






## Plot predicted and observed value, with mean and range of accuracy per
## environment and accuracy overall

# Calculate predictive ability over all observations; bootstrap to get a confidence
# interval; also calculate bias
loo_accuracy_bias_all <- loo_predictions_df %>% 
  group_by(trait, type, model, pop, selection) %>% 
  do({
    tibble(measure = c("ability", "bias"),
           out = list(neyhart::bootstrap(x = .$value, y = .$pred_complete, fun = "cor", boot.reps = 1000),
                      neyhart::bootstrap(x = .$pred_complete, y = .$value, fun = "bias", boot.reps = 1000)))
    }) %>% ungroup() %>%
  unnest() %>%
  select(-statistic)



# First create annotation df
loo_prediction_accuracy_annotation <- loo_accuracy_bias_all %>%
  filter(measure == "ability") %>%
  rename(ability_all = base, ability_lower = ci_lower, ability_upper = ci_upper) %>%
  mutate_at(vars(contains("ability")), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_all_annotation = paste0("r[MP]==", ability_all, "~(", ability_lower, "*','~", ability_upper, ")"))


# Plot predicted versus observed value
g_loo_prediction_list <- loo_predictions_df %>%
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
    
    breaks <- map(select(df1, pred_complete, value), pretty) 
    
    
    ## Extract the appropriate accuracy annotation
    r_mp_annotation <- left_join(distinct(df1, trait, type, model, pop, selection), loo_prediction_accuracy_annotation,
                                 by = c("trait", "type", "model", "pop", "selection")) %>%
      mutate(annotation = ability_all_annotation,
             x = breaks$pred_complete[1], y = c(last(breaks$value)),
             # Edit the covariate selection variable
             selection = paste0("Covariates:~", f_ec_selection_replace(selection)))
    
    ## Create the plot
    g_plot <- df1 %>%
      # Edit the covariate selection variable
      mutate(selection = paste0("Covariates:~", f_ec_selection_replace(selection))) %>%
      ggplot(aes(x = pred_complete, y = value, color = environment)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = r_mp_annotation, aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE,
                hjust = 0, size = 2) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(selection + pop ~ model, switch = "y", 
                 labeller = labeller(selection = label_parsed, pop = f_pop_replace, model = f_model_replace)) +
      labs(subtitle = paste0(toupper(unique(df1$type)), ": ", str_add_space(unique(df1$trait)))) +
      theme_presentation2(10)
    
    ## Return the plot
    g_plot
    
    }) %>% ungroup()


for (i in seq_len(nrow(g_loo_prediction_list))) {
  
  # Save
  filename <- paste0("loo_model_predictions_observations_", 
                     paste0(map_chr(select(g_loo_prediction_list, -plot), i), collapse = "_"), 
                     ".jpg")
  ggsave(filename = filename, plot = g_loo_prediction_list$plot[[i]], 
         path = fig_dir, width = 12, height = 10, dpi = 1000)
  
}



## Save realistic data
for (i in seq_len(nrow(g_loo_prediction_list))) {
  
  # Edit the plot
  g_new <- g_loo_prediction_list$plot[[i]] + 
    facet_grid(type ~ model, switch = "y", labeller = labeller(type = f_type_replace, model = f_model_replace))
  g_new <- g_new %>%
    modify_at(.x = ., .at = "data", ~
                filter(., model %in% names(model_present), pop == "tp"))
  g_new$layers[[3]]$data <- filter(g_new$layers[[3]]$data, model %in% names(model_present), pop == "tp") 
  
  # Save
  ggsave(filename = paste0("loo_model_predictions_observations_", g_loo_prediction_list$trait[[i]], "_realistic.jpg"), 
         plot = g_new, path = fig_dir, width = 8, height = 5, dpi = 1000)
  
}




## Barplot of overall prediction accuracy
g_loo_predictions_all_summ <- loo_accuracy_bias_all %>%
  mutate(selection = f_ec_selection_replace(selection)) %>%
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
       path = fig_dir, width = 12, height = 6, dpi = 1000)

## Add points for bias
g_loo_predictions_all_summ_alt <- g_loo_predictions_all_summ +
  geom_point(data = subset(loo_accuracy_bias_all, measure == "bias"),
             aes(y = base), position = position_dodge(0.9)) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty,
                     sec.axis = sec_axis(trans = ~ ., name = "Bias"))
  


# Remove unrealistic levels
g_loo_predictions_all_summ1 <- g_loo_predictions_all_summ$data %>%
  filter(., model %in% names(model_present), pop == "tp") %>%
  mutate(trait = str_add_space(trait) %>% str_replace(., " ", "\n")) %>%
  ggplot(aes(x = trait, group = model)) +
  geom_hline(yintercept = 0, color = "grey85") +
  geom_col(aes(y = base, fill = model), position = position_dodge(0.75), color = "black", width = 0.75) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(0.75),
                width = 0.5, color = "black") + 
  scale_fill_manual(name = "Model", labels = f_model_replace, values = model_colors) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Trait", labels = f_pop_replace) +
  facet_grid(type ~ ., labeller = labeller(type = f_type_replace), switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_all_accuracy_realistic.jpg", plot = g_loo_predictions_all_summ1,
       path = fig_dir, width = 5, height = 6, dpi = 1000)





## Plot bias versus accuracy
# Transform
loo_accuracy_bias_all1 <- loo_accuracy_bias_all %>% 
  select(trait, type, model, pop, measure, base, ci_upper, ci_lower) %>%
  gather(x, y, base, ci_upper, ci_lower) %>% 
  unite(x, measure, x) %>% 
  spread(x, y)


## Plot bar plot of accuracy with points for bias
g_loo_predictions_all_summ1_alt <- loo_accuracy_bias_all1 %>%
  filter(., model %in% names(model_present), pop == "tp") %>%
  mutate(trait = str_add_space(trait) %>% str_replace(., " ", "\n")) %>%
  ggplot(aes(x = trait, group = model)) +
  geom_hline(yintercept = 0, color = "grey85") +
  geom_col(aes(y = ability_base, fill = model), position = position_dodge(0.75), color = "black", width = 0.75) + 
  geom_errorbar(aes(ymin = ability_ci_lower, ymax = ability_ci_upper), position = position_dodge(0.75),
                width = 0.5, color = "black") + 
  geom_point(aes(y = bias_base * 10, shape = "Bias"), position = position_dodge(0.75)) +
  scale_fill_manual(name = "Model", labels = f_model_replace, values = model_colors) +
  scale_y_continuous(name = expression("Predictive ability"~(italic(r[MP]))), breaks = pretty,
                     sec.axis = sec_axis(trans = ~ . * 10, name = "Bias (%)", breaks = pretty)) +
  scale_x_discrete(name = "Trait", labels = f_pop_replace) +
  scale_shape_discrete(name = NULL) +
  facet_grid(type ~ ., labeller = labeller(type = f_type_replace), switch = "y") +
  theme_presentation2(10) +
  theme(strip.placement = "outside")
  
# Save
ggsave(filename = "loo_model_predictions_all_accuracy_bias_realistic.jpg", plot = g_loo_predictions_all_summ1_alt,
       path = fig_dir, width = 5, height = 6, dpi = 1000)




g_loo_ability_bias_all <- loo_accuracy_bias_all1 %>%
  ggplot(aes(x = bias_base, y = ability_base, shape = pop, color = model)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = bias_ci_lower, xend = bias_ci_upper, y = ability_base, yend = ability_base), 
               color = "grey85", lwd = 0.5) +
  geom_segment(aes(x = bias_base, xend = bias_base, y = ability_ci_lower, yend = ability_ci_upper), 
               color = "grey85", lwd = 0.5) +
  geom_point(size = 2) +
  scale_color_paletteer_d(package = "dutchmasters", palette = "milkmaid",
                         name = "Model", labels = f_model_replace) +
  scale_shape_discrete(name = "Validation\nscheme", labels = f_pop_replace) +
  scale_x_continuous(name = "Bias (%)", breaks = pretty, labels = scales::percent_format(accuracy = 1, suffix = NULL)) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y", scales = "free_x") +
  theme_presentation2(base_size = 8)

# Save
ggsave(filename = "loo_model_predictions_all_summary.jpg", plot = g_loo_ability_bias_all,
       path = fig_dir, width = 10, height = 6, dpi = 1000)



# Remove unrealistic levels
g_loo_ability_bias_all1 <- g_loo_ability_bias_all %>%
  modify_at(.x = ., .at = "data", ~
              filter(., model %in% names(model_present), pop == "tp"))
g_loo_ability_bias_all1 <- g_loo_ability_bias_all1 +
  facet_grid(type ~ trait, scales = "free_x", labeller = labeller(type = f_type_replace), switch = "y") +
  scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) +
  scale_shape_discrete(guide = FALSE) +
  theme_presentation2(base_size = 10)

# Save
ggsave(filename = "loo_model_predictions_all_summary_realistic.jpg", plot = g_loo_ability_bias_all1,
       path = fig_dir, width = 8, height = 4, dpi = 1000)

# Save as HTML widget
htmlwidgets::saveWidget(widget = plotly::as_widget(plotly::ggplotly(g_loo_ability_bias_all)),
                        file = file.path(fig_dir, "loo_model_predictions_all_summary.html"))









### Summarize accuracy within each environment ###




## Create a function list for calculating quantiles
q <- c(alpha / 2, 0.5, 1 - (alpha / 2))
# Function vector
q_funs <- c(map(q, ~partial(quantile, probs = .x, na.rm = TRUE)), partial(mean, na.rm = TRUE)) %>%
  set_names(., c("lower", "median", "upper", "mean"))


## Summarize mean and range
loo_prediction_accuracy_summ <- loo_prediction_accuracy %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, model, pop, selection) %>%
  summarize_at(vars(ability, bias, accuracy), funs(!!!q_funs)) %>%
  ungroup()

## Determine the LSD between model-selection groups
loo_prediction_accuracy_summ1 <- loo_prediction_accuracy %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, pop) %>%
  do({
    dat <- .
    fit <- lm(accuracy ~ model + selection + model:selection, data = dat)
    n <- unique(as.numeric(xtabs(~ model + selection, dat)))
    MS_within <- sigma(fit)^2 # Within-group error (residuals)
    df <- df.residual(fit)
    LSD <- qt(p = 1 - (0.05 / 2), df = df) * sqrt(MS_within * (2 / n))
    
    # Return mean and LSD of prediction accuracy
    aggregate(accuracy ~ model + selection, dat, mean) %>% 
      mutate(LSD = LSD)
  }) %>% ungroup()

annotation_df <- loo_prediction_accuracy_summ1 %>%
  distinct(type, trait, pop, LSD) %>%
  mutate(annotation = paste0("LSD: ", round(LSD, 2)))

## Visualize
# Plot mean and range of predictions
g_loo_predictions_summ <- loo_prediction_accuracy_summ1 %>%
  filter(model %in% names(model_replace)) %>%
  mutate(selection = f_ec_selection_replace(selection)) %>%
  ggplot(aes(x = selection, group = model, color = model)) +
  # geom_linerange(aes(ymin = ability_lower, ymax = ability_upper, group = model), 
  #                position = position_dodge(0.9), color = "grey85") +
  geom_point(aes(y = accuracy), position = position_dodge(0.9)) +
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




# Remove unrealistic levels
g_loo_predictions_summ1 <- loo_prediction_accuracy_summ %>%
  filter(model %in% names(model_present),
         pop == "vp") %>%
  # rename(mean = ability_mean, upper = ability_upper, lower = ability_lower) %>% # Ability
  rename(mean = accuracy_mean, upper = accuracy_upper, lower = accuracy_lower) %>% # Accuracy
  mutate(trait = str_add_space(trait) %>% str_replace_all(" ", "\n")) %>%
  ggplot(aes(x = trait, color = model)) +
  geom_linerange(aes(ymin = lower, ymax = upper, group = model), position = position_dodge(0.65), color = "grey85") +
  geom_point(aes(y = mean), position = position_dodge(0.65)) +
  scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) +
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  scale_x_discrete(name = "Trait", labels = f_pop_replace) +
  facet_grid(type ~ ., labeller = labeller(type = f_type_replace),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary_realistc.jpg", plot = g_loo_predictions_summ1,
       path = fig_dir, width = 5, height = 4, dpi = 1000)


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







## Practical exercises ##

## Load the fitted ammi models
load(file.path(result_dir, "ammi_model_fit.RData"))
n_select <- 20

# Make predictions in each location using the AMMI results
ammiN_fit_location1 <- ammiN_fit_location %>%
  mutate(y_hat = map2(g_scores, loc_scores, ~{
    # Effects
    g_eff <- select(.x, line_name, effect)
    l_eff <- select(.y, location, effect)
    gl_eff_mat <- outer(X = .x[["score"]], Y = .y[["score"]], FUN = "*")
    dimnames(gl_eff_mat) <- list(g_eff$line_name, l_eff$location)
    gl_eff <- gl_eff_mat %>%
      as.data.frame() %>%
      rownames_to_column("line_name") %>%
      gather(location, effect, -line_name)
    # Calculate y_hat
    full_join(gl_eff, g_eff, by = "line_name") %>% 
      full_join(., l_eff, by = "location") %>% 
      mutate(y_hat = effect.x + effect.y + effect) %>%
      select(line_name, location, y_hat)
  })) %>%
  mutate(y_hat = map2(mu, y_hat, ~mutate(.y, y_hat = y_hat + .x)))


# For each location, find the best 10 varieties, according to the AMMI model
tp_ammi_location_winners <- ammiN_fit_location1 %>% 
  unnest(y_hat) %>%
  filter(line_name %in% tp_geno) %>%
  # Add trait sign and location information
  left_join(., trait_sign, by = "trait") %>%
  group_by(trait, location) %>%
  top_n(x = ., n = n_select, wt = (y_hat * sign)) %>%
  ungroup()


# Find the ideal variety-location combinations in the TP
tp_prediction_location_winners <- loo_predictions_df %>%
  filter(pop == "tp", type == "lolo") %>%
  group_by(trait, location, model, line_name) %>%
  summarize_at(vars(value, pred_complete), mean) %>%
  # Assign trait signs
  left_join(., trait_sign, by = "trait") %>%
  top_n(x = ., n = n_select, wt = (pred_complete * sign)) %>%
  ungroup()

## For each trait, location, and model, determine the proportion of recovered genotypes
tp_location_winners_recovery <- full_join(
  x = nest(group_by(tp_prediction_location_winners, trait, location, model), line_name, .key = "prediction"),
  y = nest(group_by(tp_ammi_location_winners, trait, location), line_name, .key = "ammi"),
  by = c("trait", "location")
) %>%
  mutate(pRecover = map2_dbl(.x = prediction, .y = ammi, ~sum(.x$line_name %in% .y$line_name)/nrow(.y)))

## Calculate mean and range for each trait and model
(tp_location_winners_recovery_summ <- tp_location_winners_recovery %>% 
  filter(str_detect(model, "a$", negate = TRUE)) %>%
  group_by(trait, model) %>% 
  summarize_at(vars(pRecover), list(~mean, ~min, ~max)) %>%
  as.data.frame())



## Plot
g_tp_location_winners <- tp_location_winners_recovery_summ %>%
  ggplot(aes(x = trait, y = mean, fill = model)) +
  geom_col(position = position_dodge(0.8), color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.5, position = position_dodge(0.8)) +
  scale_fill_paletteer_d(package = "dutchmasters", palette = "milkmaid",
                         name = "Model", labels = f_model_replace) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = paste0("Proportion of best lines recovered (n/", n_select, ")"), breaks = pretty, limits = c(0,1)) +
  theme_acs(10)

# Save
ggsave(filename = "tp_location_winner_recovery.jpg", plot = g_tp_location_winners,
       path = fig_dir, width = 8, height = 4, dpi = 1000)


## Plot realistic levels
g_tp_location_winners1 <- g_tp_location_winners %>%
  modify_at(., "data", ~filter(., model %in% names(model_present))) +
  scale_fill_manual(name = "Model", labels = f_model_replace, values = model_colors) +

# Save
ggsave(filename = "tp_location_winner_recovery_realistic.jpg", plot = g_tp_location_winners1,
       path = fig_dir, width = 6, height = 4, dpi = 1000)

# Remove range bars
g_tp_location_winners2 <- g_tp_location_winners1
g_tp_location_winners2$layers <- g_tp_location_winners2$layers[1]
ggsave(filename = "tp_location_winner_recovery_realistic2.jpg", plot = g_tp_location_winners2,
       path = fig_dir, width = 6, height = 4, dpi = 1000)
 









