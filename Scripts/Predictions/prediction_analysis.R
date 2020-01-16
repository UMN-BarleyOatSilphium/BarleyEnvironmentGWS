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
library(modelr)


## Load the validation results
file_list <- list.files(result_dir, pattern = "predictions.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))


## Significant level
alpha <- 0.05


## A vector to rename models
model_replace <- c("model1" = "G", "model2" = "G + E", "model3" = "G + E + GE")
f_model_replace <- function(x) str_replace_all(x, model_replace)


## For LOYO and LOEO, determine the number of size of the training population (nObs and nEnv)
## for each environment
loeo_prediction_info <- loeo_predictions_out %>% 
  unnest(predictions) %>% 
  distinct(trait, environment, line_name) %>%
  filter(line_name %in% tp) %>%
  split(.$trait) %>%
  map_df(~group_by(., trait, environment) %>% crossv_loo2()) %>%
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         nObs = map_dbl(train, nrow),
         nEnv = map_dbl(train, ~n_distinct(.$environment))) %>%
  unnest(test) %>%
  distinct(trait, environment, nObs, nEnv)

loyo_prediction_info <- loyo_predictions_out %>% 
  unnest(predictions) %>% 
  distinct(trait, year, environment, line_name) %>%
  filter(line_name %in% tp) %>%
  split(.$trait) %>%
  map_df(~group_by(., trait, year) %>% crossv_loo2()) %>%
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         nObs = map_dbl(train, nrow),
         nEnv = map_dbl(train, ~n_distinct(.$environment))) %>%
  unnest(test) %>%
  distinct(trait, environment, nObs, nEnv)

# Combine
loo_prediction_info <- bind_rows(
  mutate(loeo_prediction_info, type = "loeo"),
  mutate(loyo_prediction_info, type = "loyo")
)



## Combine LOEO and LOYO results
loo_predictions_df <- bind_rows(
  mutate(loeo_predictions_out, type = "loeo") %>% select(-environment, -core),
  mutate(loyo_predictions_out, type = "loyo") %>% select(-year, -core)
) %>%
  unnest(predictions) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp"),
         environment = as.character(environment))
  
## Calculate accuracy and bias per environment, model, and population
loo_predictive_ability <- loo_predictions_df %>%
  group_by(trait, environment, model, pop, type) %>%
  summarize(ability = cor(prediction, value), 
            bias = mean(prediction - value)) %>%
  ungroup() %>%
  left_join(., loo_prediction_info)


## Adjust ability using heritability
loo_prediction_accuracy <- loo_predictive_ability %>%
  left_join(., env_trait_herit) %>%
  mutate(accuracy = ability / sqrt(heritability))




## Plot predicted and observed value, with accuracy for each environment
## and overall

# First create annotation for accuracy
loo_predictions_df1 <- loo_prediction_accuracy %>%
  mutate(annotation = paste0(environment, "~(r[MP]==", formatC(x = ability, 2, format = "g"), ")")) %>%
  left_join(loo_predictions_df, .)


# Plot predicted versus observed value
g_loo_prediction_list <- loo_predictions_df1 %>%
  # Max character length of units
  mutate(max_nchar = max(nchar(last(pretty(value))))) %>%
  group_by(trait, type) %>%
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
g_loo_prediction_list1 <- g_loo_prediction_list %>%
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
         path = fig_dir, width = 9, height = 9, dpi = 1000)
  
}







## Summarize mean and range


loo_prediction_accuracy_summ <- loo_prediction_accuracy %>%
  group_by(type, trait, model, pop) %>%
  summarize_at(vars(ability, bias, accuracy), list(~mean, ~min, ~max)) %>%
  ungroup()



## Visualize
# Plot mean and range of predictions
g_loo_predictions_summ <- loo_prediction_accuracy_summ %>%
  ggplot(aes(x = pop, color = model)) +
  geom_linerange(aes(ymin = ability_min, ymax = ability_max, group = model), 
                 position = position_dodge(0.9), color = "grey85") +
  geom_point(aes(y = ability_mean), position = position_dodge(0.9)) +
  scale_color_discrete(name = "Model", labels = f_model_replace) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Population", labels = str_to_upper) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary.jpg", plot = g_loo_predictions_summ,
       path = fig_dir, width = 6, height = 4, dpi = 1000)


# alternate
g_loo_predictions_summ1 <- loo_prediction_accuracy %>%
  left_join(., loo_prediction_accuracy_summ) %>%
  ggplot(aes(x = pop, fill = model)) +
  geom_boxplot(aes(y = ability), position = position_dodge(0.9), color = "grey85", lwd = 0.5) +
  geom_point(aes(y = ability_mean), position = position_dodge(0.9), size = 1.5) +
  scale_shape_discrete(name = "Statistic", labels = str_to_title) + 
  scale_fill_discrete(name = "Model", labels = f_model_replace) + 
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = "Population", labels = str_to_upper) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary1.jpg", plot = g_loo_predictions_summ1,
       path = fig_dir, width = 6, height = 4, dpi = 1000)


# Use accuracy instead of ability
g_loo_predictions_summ2 <- loo_prediction_accuracy %>%
  left_join(., loo_prediction_accuracy_summ) %>%
  ggplot(aes(x = pop, fill = model)) +
  geom_boxplot(aes(y = accuracy), position = position_dodge(0.9), color = "grey85", lwd = 0.5) +
  geom_point(aes(y = accuracy_mean), position = position_dodge(0.9), size = 1.5) +
  scale_shape_discrete(name = "Statistic", labels = str_to_title) + 
  scale_fill_discrete(name = "Model", labels = f_model_replace) + 
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  scale_x_discrete(name = "Population", labels = str_to_upper) +
  facet_grid(type ~ trait, labeller = labeller(trait = str_add_space, type = toupper),
             switch = "y") +
  theme_presentation2(10)

# Save
ggsave(filename = "loo_model_predictions_summary2.jpg", plot = g_loo_predictions_summ2,
       path = fig_dir, width = 6, height = 4, dpi = 1000)



## Create a table
loo_prediction_accuracy_table <- loo_prediction_accuracy %>%
  group_by(type, trait, model, pop) %>%
  summarize(accuracy1 = weighted.mean(x = accuracy, w = nObs),
            accuracy2 = weighted.mean(x = accuracy, w = nEnv),
            accuracy = mean(accuracy)) %>%
  ungroup()
  
loo_prediction_accuracy_table1 <- loo_prediction_accuracy_table %>%
  select(type, trait, model, pop, accuracy2) %>%
  spread(pop, accuracy2)

write_csv(x = loo_prediction_accuracy_table1, path = file.path("loo_prediction_accuracy_table.csv"))

