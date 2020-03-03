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
library(broom)
library(car)
library(lmerTest)
library(patchwork)
library(cowplot)
library(modelr)


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
model_replace <- c("model1" = "G", "model2" = "G + E", "model3" = "G + E + GE",
                   "model4" = "G + L", "model5" = "G + L + GL")
f_model_replace <- function(x) str_replace_all(x, model_replace)
# f_model_replace <- function(x) paste0("M", toupper(str_extract(x, "[0-9]{1}[a-z]{0,1}")))
# Vector to rename validation schemes
f_type_replace <- function(x) str_replace_all(x, c("tp" = "CV0", "vp" = "POV00"))

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
loo_prediction_list <- map(set_names(object_list, object_list), get)



loyo_predictions_out1 <- loyo_predictions_out %>%
  group_by(trait, year) %>% 
  nest(out) %>% 
  mutate(out = map(data, ~mutate(.[[1]][[1]], train_n = list(.[[1]][[2]])))) %>% 
  unnest(out) %>% 
  unnest(train_n)

loeo_predictions_out1 <- loeo_predictions_out %>%
  group_by(trait, env) %>% 
  nest(out) %>%
  mutate(out = map(data, ~mutate(.[[1]][[1]], train_n = list(.[[1]][[2]])))) %>% 
  unnest(out) %>% 
  unnest(train_n)


## Create a list manually
loo_prediction_list <- list(
  lolo = lolo_predictions_out,
  loeo = loeo_predictions_out1,
  loyo = loyo_predictions_out1
)


loo_prediction_out <- loo_prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~select(., trait, environment, location, year, type, model, line_name, value, contains("pred")) %>%
        mutate_if(is.character, parse_guess) %>%
        mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  # Filter for models listed above
  filter(model %in% names(model_replace))

loo_predictions_df <- loo_prediction_out


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
  # Now summarize
  group_by(trait, model, pop, type, environment) %>%
  summarize_at(vars(ability, bias, ability_all, bias), mean) %>%
  ungroup()


## Adjust ability using heritability
loo_prediction_accuracy <- loo_predictive_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))



## Plot predicted and observed value, with mean and range of accuracy per
## environment and accuracy overall

# First create annotation df
loo_prediction_accuracy_annotation <- loo_prediction_accuracy %>%
  group_by(trait, model, pop, type) %>%
  ## Calculate min, max, mean environment-specific ability
  summarize(ability_min = min(ability), ability_max = max(ability),
            ability_mean = mean(ability), ability_all = mean(ability_all)) %>%
  ungroup() %>%
  mutate_at(vars(contains("ability")), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_annotation = paste0("bar(r)[MP]==", ability_mean, "~(", ability_min, "*','~", ability_max, ")"),
         ability_all_annotation = paste0("r[MP]==", ability_all))

# Plot predicted versus observed value
g_loo_prediction_list <- loo_predictions_df %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait, type) %>%
  do(plot = {
    df <- .
    
    df1 <- df

    max_wid <- unique(df1$max_nchar)
    # Scale truncation function
    scale_trunc <- function(x) str_pad(string = x, width = max_wid, side = "left", pad = " ")
    # Calculate breaks for this trait
    breaks <- map(select(df1, pred_complete, value), pretty) 
    
    # Create a separate legend df
    x1 <- mutate(df1, x = ifelse(model == "model1", 0.20, 0.80), y = 0.10)
    ## Split x by model and population
    x_split <- split(x1, list(x1$pop, x1$model))
    
    # Map over this split
    g_list <- vector("list", length = length(x_split))
    
    for (i in seq_along(g_list)) {
      x2 <- x_split[[i]]
      
      ## Extract the appropriate accuracy annotation
      r_mp_annotation <- left_join(distinct(x2, trait, type, model, pop), loo_prediction_accuracy_annotation,
                                   by = c("trait", "type", "model", "pop")) %>%
        select(trait:pop, contains("annotation")) %>% 
        mutate(annotation = paste0("atop(", ability_annotation, ", ", ability_all_annotation, ")")) %>%
        mutate(x = breaks$pred_complete[1], y = c(last(breaks$value)))
    
      g_x <- ggplot(data = x2, aes(x = pred_complete, y = value, color = environment)) +
        geom_abline(slope = 1, intercept = 0) +
        geom_point(size = 0.5) +
        geom_text(data = r_mp_annotation, aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE,
                  hjust = 0, size = 2) +
        scale_y_continuous(name = "Observation", breaks = breaks$value, labels = scale_trunc, limits = range(breaks$value)) +
        scale_x_continuous(name = "Prediction", breaks = breaks$pred_complete, limits = range(breaks$pred_complete)) +
        scale_color_discrete(name = NULL, labels = function(x) parse(text = x), guide = FALSE) +
        theme_presentation2(10) +
        theme(axis.title = element_blank())
      
      ## Add faceting depending on position
      if (i == 1) {
        g_list[[i]] <- g_x + 
          facet_grid(trait ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace))
      } else {
        g_list[[i]] <- g_x + 
          facet_grid(. ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace)) +
          theme(axis.text.y = element_blank())
        
      }
      
    } # Close loop
    
    ## Combine the list of plots and return - this is plot version 1
    combined_plot1 <- plot_grid(plotlist = g_list, nrow = 1, 
                                rel_widths = c(1, rep(0.75, length(g_list) - 1)))
    
    ## Extract the appropriate accuracy annotation
    r_mp_annotation <- left_join(distinct(df1, trait, type), loo_prediction_accuracy_annotation,
                                 by = c("trait", "type")) %>%
      select(trait:pop, contains("annotation")) %>% 
      # Code below create two lines that are left justified using atop
      mutate(annotation = paste0("atop(", ability_annotation, ", ", 
                                 paste0(str_pad(ability_all_annotation, width = nchar(ability_annotation)-2, 
                                                side = "right", pad = "~")), "' ')")) %>%
      mutate(x = breaks$pred_complete[1], y = c(last(breaks$value)))


    ## Plot version 2 is simply a grid without the complicated legend
    # Next, plot
    combined_plot2 <- df1 %>%
      ggplot(aes(x = pred_complete, y = value, color = environment, group = 1)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point(size = 0.5) +
      geom_text(data = r_mp_annotation, aes(x = x, y = y, label = annotation), parse = TRUE, inherit.aes = FALSE,
                hjust = 0, size = 2) +
      scale_y_continuous(name = "Observation", breaks = pretty, labels = scale_trunc) +
      scale_x_continuous(name = "Prediction", breaks = pretty) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(trait ~ model + pop, switch = "y", labeller = labeller(trait = str_add_space, model = f_model_replace)) +
      theme_presentation2(10) +
      theme(axis.title = element_blank())

    ## Return both plots
    tibble(plot_name = c("plot1", "plot2"), plot = list(combined_plot1, combined_plot2))
    
    }) %>% ungroup()


# ### Combine plots (type 1) ###
# g_loo_predictions2 <- g_loo_prediction_list %>%
#   split(.$type) %>%
#   map(~mutate(.x, plot = modify_at(.x = plot, .at = -1, ~. + theme(strip.text.x = element_blank())))) %>%
#   map(~plot_grid(plotlist = .x$plot, align = "v", ncol = 1, rel_heights = c(1, rep(0.8, length(.) - 1))))


### Combine plots (type 2) ###
g_loo_prediction_list1 <- g_loo_prediction_list %>%
  unnest(plot) %>%
  filter(plot_name == "plot2") %>%
  split(.$type) %>%
  # Remove strip names from all but first element
  map(~mutate(., plot = modify_at(.x = plot, .at = -1, .f = ~. + theme(strip.text.x = element_blank()))))

# Use cowplot
g_loo_predictions2 <- g_loo_prediction_list1 %>%
  map(., ~pull(., plot) %>%
        plot_grid(plotlist = ., align = "v", ncol = 1, rel_heights = c(1, rep(0.8, length(.) - 1))) )


for (type in names(g_loo_predictions2)) {
  
  # Save
  ggsave(filename = paste0(type, "_model_predictions_observations.jpg"), plot = g_loo_predictions2[[type]],
         path = fig_dir, width = 15, height = 10, dpi = 1000)
  
}







## Summarize mean and range
loo_prediction_accuracy_summ <- loo_prediction_accuracy %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, model, pop) %>%
  summarize_at(vars(ability, bias, accuracy), 
               list(~mean, ~min, q25 = ~quantile(., 0.25), q75 = ~quantile(., 0.75), ~max)) %>%
  ungroup()



## Visualize
# Plot mean and range of predictions
g_loo_predictions_summ <- loo_prediction_accuracy_summ %>%
  ggplot(aes(x = pop, color = model)) +
  geom_linerange(aes(ymin = ability_min, ymax = ability_max, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  geom_linerange(aes(ymin = ability_q25, ymax = ability_q75, group = model), 
                 position = position_dodge(0.9), color = "grey85", lwd = 1) +
  geom_point(aes(y = ability_mean), position = position_dodge(0.9)) +
  scale_color_manual(name = "Model", labels = f_model_replace, values = model_colors) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Validation scheme", labels = f_type_replace) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary.jpg", plot = g_loo_predictions_summ,
       path = fig_dir, width = 8, height = 6, dpi = 1000)




## Create a table

# Report the mean (and range) in predictive ability across environments for each
# validation scheme, population, and type
loo_prediction_accuracy_table <- loo_prediction_accuracy_summ %>% 
  mutate_at(vars(contains("ability")), ~formatC(., digits = 2, width = 2, format = "g")) %>%
  mutate(annotation = paste0(ability_mean, " (", ability_min, ", ", ability_max, ")")) %>%
  # Rename
  mutate(model = f_model_replace(model),
         pop = f_type_replace(pop),
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
         pop = f_type_replace(pop),
         type = toupper(type)) %>%
  select(trait, type, model, pop, annotation) %>%
  spread(model, annotation) %>%
  arrange(trait, pop, type)

write_csv(x = loo_prediction_accuracy_table1, path = file.path(fig_dir, "loo_prediction_accuracy_table1.csv"))



