## S2MET genotype-environment ecophysiological analysis
## 
## Determine environmental indices by examining the effect of
## environmental covariates during critical growth stages
## 
## Author: Jeff Neyhart
## Last modified: 4 October  2019
## 


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load other libraries
library(lubridate)
library(modelr)
library(broom)
library(pbr)
library(cowplot)
library(pls)
library(car)

# Significance level
alpha <- 0.05

## Environments to use with corresponding trials
env_trials <- S2_MET_BLUEs %>% 
  distinct(trial, environment) %>%
  filter(environment %in% tp_vp_env,
         str_detect(trial, "S2C1", negate = TRUE))



############################
### Load data
############################


## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))

## Load the fitted ammi models
load(file.path(result_dir, "ammi_model_fit.RData"))



############################
### Prepare ECs for modeling
############################


## Create the ECs and select the relevant ones for modeling
ec_select <- growth_stage_covariates %>%
  ## Remove heading as a growth stage
  filter(stage != "heading") %>%
  ## Only use TP environments
  inner_join(., env_trials, by = "trial") %>% 
  select(-trial) %>%
  gather(covariate, value, -environment, -stage) %>%
  unite(covariate, stage, covariate, sep = ".") %>%
  ## Remove trange, relhum, and rain (water stress will cover this)
  filter(str_detect(covariate, "radn_mean", negate = TRUE)) %>%
  spread(covariate, value)

## Summarize min/max/var for each covariate
ec_select_summ <- ec_select %>%
  gather(covariate, value, -environment) %>%
  group_by(covariate) %>%
  summarize_at(vars(value), list(min = min, max = max, var = var, cv = ~sd(.) / mean(.)))


##### 
# Look at correlations between covariates
##### 


## Calculate pairwise correlations
ec_pairwise_cor <- cor(ec_select[,-1]) %>% 
  as.dist() %>%
  tidy() %>%
  rename(covariate1 = item1, covariate2 = item2, correlation = distance) %>%
  # Sort by descending absolution cor coef
  arrange(desc(abs(correlation)))


## Are there covariates that tend to be highly correlated with others?
ec_pairwise_cor %>% 
  group_by(covariate1) %>% 
  summarize(correlation = mean(correlation)) %>% 
  arrange(desc(abs(correlation)))



## Separate covariates into stage and measurement
ec_pairwise_cor1 <- ec_pairwise_cor %>%
  separate(covariate1, c("stage1", "measurement1"), sep = "\\.") %>%
  separate(covariate2, c("stage2", "measurement2"), sep = "\\.")


## Are there measurements that tend to be highly correlated with others?
ec_pairwise_cor1 %>% 
  group_by(measurement1, measurement2) %>% 
  summarize(correlation = mean(correlation)) %>%
  arrange(desc(abs(correlation)))



## Which covariates are enriched beyond some threshold?
cor_threshold <- 0.6
ec_pairwise_cor %>% 
  filter(correlation >= cor_threshold) %>%
  group_by(covariate2) %>%
  summarize(n = n()) %>%
  arrange(desc(n))



## Filter out some covariates
to_remove <- c("radn_mean", "tmean_mean")





## Prepare the covariates for modeling
## Center, but do not scale. Save the mean for later
ec_tomodel_temp <- ec_select %>%
  mutate_at(vars(-environment), scale, scale = FALSE, center = TRUE)

ec_tomodel_centers <- ec_tomodel_temp %>%
  summarize_at(vars(-environment), ~attr(., "scaled:center")) %>%
  gather(covariate, center)

## Convert scaled to numeric
ec_tomodel_centered <- ec_tomodel_temp %>%
  mutate_at(vars(-environment), as.numeric)

## Center and scale
ec_tomodel_scaled <- ec_select %>%
  mutate_at(vars(-environment), scale, scale = TRUE, center = TRUE) %>%
  mutate_at(vars(-environment), as.numeric)


# Vector of covariates
environmental_covariates <- names(ec_tomodel_centered)[-1]







## Test for normality using ks test
ec_tomodel_normality <- ec_select %>%
  gather(covariate, value, -environment) %>%
  group_by(covariate) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value), sd = sd(.$value))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value"))

subset(ec_tomodel_normality, p_value < 0.1)





## Prepare data for modelling:
## 1. only use the tp
## 2. add ECs to the df

s2_met_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  left_join(select(ec_tomodel_centered, c("environment", subset(ec_tomodel_normality, p_value >= 0.10, covariate, drop = TRUE))), by = "environment") %>%
  mutate_at(vars(line_name, environment), as.factor)


## What is the level of coverage in this dataset?
s2_met_tomodel %>% 
  group_by(trait, line_name) %>%
  summarize(n = n_distinct(environment)) %>% 
  mutate(mean_env = n / max(n)) %>%
  top_n(x = ., n = 5, wt = -mean_env)



# New vector of covariates
environmental_covariates <- s2_met_tomodel %>% 
  select(matches("vegetative|heading|flowering|grain_fill")) %>%
  names()

## Scatterplot matrix
scatterplotMatrix(distinct(select(s2_met_tomodel, environmental_covariates)), smooth = FALSE)


## Determine if there is sufficient variation for a covariate
distinct(select(s2_met_tomodel, environmental_covariates)) %>%
  map_dbl(var)

## Assign covariates for each trait
ec_by_trait_apriori <- list(
  GrainProtein = c("grain_fill.maxt_mean", "grain_fill.water_balance_sum"),
  GrainYield = c("flowering.mint_mean", "grain_fill.maxt_mean", "grain_fill.water_balance_sum"),
  TestWeight = c("flowering.mint_mean", "grain_fill.maxt_mean", "grain_fill.water_balance_sum"),
  PlantHeight = str_subset(environmental_covariates, pattern = "grain_fill", negate = TRUE),
  HeadingDate = str_subset(environmental_covariates, pattern = "grain_fill|flowering|heading", negate = TRUE)
) %>%
  tibble(trait = names(.), covariates = .)

# # New vector of covariates
# environmental_covariates_apriori <- reduce(ec_by_trait_apriori$covariates, union)



## Model all covariates without apriori information
ec_by_trait_all <- tibble(
  trait = sort(traits),
  covariates = list(
    environmental_covariates, # Grain protein
    environmental_covariates, # Grain yield
    str_subset(environmental_covariates, "vegetative"), # Heading date
    str_subset(environmental_covariates, "grain", negate = TRUE), # Plant height
    environmental_covariates # Test weight
  )
)

# New vector of covariates
environmental_covariates_all <- reduce(ec_by_trait_all$covariates, union)




############################
### Model covariates that are predictive of the environmental mean
############################

## Get the environmental effects from the AMMI model
base_model_env_effect <- ammiN_fit_location %>%
  unnest(e_scores) %>% 
  select(trait, environment, effect)



## Use stepwise regression to find the best model
env_effect_to_model <- base_model_env_effect %>%
  left_join(s2_met_tomodel)

# Group by trait
env_effect_models <- env_effect_to_model %>%
  group_by(trait) %>%
  nest() %>%
  mutate(model = list(NULL))



for (i in seq(nrow(env_effect_models))) {

  df <- env_effect_models$data[[i]]
  tr <- env_effect_models$trait[i]

  ## Remove some covariates depending on when the trait is measured
  df1 <- select(df, environment, effect, subset(ec_by_trait_all, trait == tr, covariates, drop = T)[[1]]) %>%
    distinct()

  # Minimal model
  min_model <- effect ~ 1
  # Max model
  max_model <- add_predictors(min_model, as.formula(paste0("~", paste0(names(df1)[-1:-2], collapse = " + "))))

  # Stepwise regression
  fit_base <- lm(formula = min_model, data = df1, contrasts = "contr.sum")
  fit_step <- step(object = fit_base, scope = max_model, direction = "forward")
  
  
  # Did we run out of dfs? If so, reduce terms
  while (df.residual(fit_step) <= 1) {
    
    ## Drop1
    fit_step_drop1 <- drop1(fit_step) %>%
      as.data.frame() %>%
      rownames_to_column("term")
    
    ## Remove the two covariates that gives the lowest AIC, besides "none"
    to_remove <- fit_step_drop1 %>%
      filter(str_detect(term, "none", negate = T)) %>% 
      top_n(x = ., n = 2, wt = -AIC) %>%
      pull(term) %>% 
      sort() %>%
      first()
    
    ## Reduce the number of parameters by 1
    fit_step_new_formula <- terms(formula(fit_step)) %>% 
      drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% to_remove), keep.response = TRUE) %>%
      formula()
    
    ## Refit the model
    fit_step <- update(object = fit_step, formula = fit_step_new_formula)

  }
  
  fit_step1 <- fit_step

  
  ## Remove covariates based on VIF
  vif_i <- vif(fit_step1)
  fit_stepi <- fit_step1
  vif_cutoff <- 4
  
  while (any(vif_i > vif_cutoff)) {
    form_i <- terms(formula(fit_stepi)) %>%
      drop.terms(., dropx = which( attr(., "term.labels") %in% names(which.max(vif_i))), keep.response = TRUE ) %>%
      formula()
    
    fit_stepi <- update(object = fit_stepi, formula = form_i)
    vif_i <- vif(fit_stepi)
  }
  
  # Backward elimination from here
  fit_step_final <- step(object = fit_stepi, direction = "backward")

  ## Return the model
  env_effect_models$model[[i]] <- fit_step_final

}


## Examine coefficients of covariates
env_effect_models$model %>% 
  map(summary)


### This seems to have helped.




## Plot fitted values versus observed environmental mean
env_effect_models %>%
  mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>%
  unnest(model_data) %>%
  ggplot(aes(x = pred, y = effect)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ trait, scales = "free") +
  theme_acs()
  
## Perform leave-one-out cross-validation
env_effect_models_looCV <- env_effect_models %>%
  mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>% 
  group_by(trait) %>%
  do({
    row <- .
    cv_df <- unnest(row, model_data) %>%
      crossv_loo() %>%
      mutate_at(vars(test, train), ~map(., as.data.frame))
    model_fits <- map(cv_df$train, ~update(object = row$model[[1]], data = .))
    
    # Predictions
    predictions <- map2_df(cv_df$test, model_fits, add_predictions)
    
    # RMSE
    rmse_df <- mean(map2_dbl(model_fits, cv_df$test, rmse))
    
    # Accuracy
    accuracy <- cor(predictions$effect, predictions$pred)
    
    ## Return summaries
    tibble(predictions = list(predictions), accuracy = accuracy, rmse = rmse_df)
    
  }) %>% ungroup()
  

# Plot
env_effect_models_looCV %>%
  unnest(predictions) %>%
  ggplot(aes(x = pred, y = effect)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ trait, scales = "free") +
  theme_acs()



# trait        predictions       accuracy     rmse
# 1 GrainProtein <tibble [10 x 7]>    0.727    0.970
# 2 GrainYield   <tibble [22 x 8]>    0.715 1083.   
# 3 HeadingDate  <tibble [25 x 6]>    0.751    3.37 
# 4 PlantHeight  <tibble [27 x 7]>    0.491    8.78 
# 5 TestWeight   <tibble [12 x 9]>    0.621   36.8 





############################
### Identify covariates that are correlated with the AMMI distance matrix
############################


## Using the AMMI output and phi, calculate the distance matrix between environments
## This will be denoted as W_ammi, consistent with Rincent 2019
ammi_dist <- ammiN_fit_location %>%
  mutate(W = map(phi, ~as.matrix(dist(t(.)))),
         W_ammi = map(W, ~1 - (. / max(.)))) %>%
  select(trait, W_ammi)

## Add the coviarates assigned to each trait
ec_tomodel_ammi <- left_join(ammi_dist, ec_by_trait_all, by = "trait")

# Create a matrix of scaled and centered covariates
ec_tomodel_scaled_mat <- ec_tomodel_scaled %>% 
  as.data.frame() %>%
  column_to_rownames("environment") %>% 
  as.matrix()


# Tolerance of difference of correlations
cor_tol <- 0.01


## For each trait, use an algorithm to determine the covariates to use
## This is based on the approach of Rincent 2019
ec_ammi_dist <- ec_tomodel_ammi %>%
  group_by(trait) %>%
  do({
    
    row <- .
    
    # AMMI dist matrix
    W_ammi <- row$W_ammi[[1]]
    # Vector of covariates
    covariates <- row$covariates[[1]]
    
    ## Subset the EC scaled matrix
    ec_tomodel_scaled_mat_use <- ec_tomodel_scaled_mat[rownames(W_ammi), covariates]
    
    ## Choose the starting covariate as the one whose distance matrix has
    ## the highest correlation with the W_ammi distance matrix
    initial_W_cov <- map(covariates, ~make_dist_mat(ec_tomodel_scaled_mat_use[, ., drop = FALSE])) %>%
      map_dbl(~cor(as.numeric(.), as.numeric(W_ammi)))
    
    starting_covariate <- covariates[which.max(initial_W_cov)]
    # Create two vectors of available and used covariates
    used_covariates <- tested_covariates <- starting_covariate
    available_covariates <- setdiff(covariates, used_covariates)
    
    old_cor <- max(initial_W_cov)
    
    ## Vector of correlations
    i <- 1
    ec_ammi_cor_list <- vector("list", length(covariates))
    ec_ammi_cor_list[[i]] <- tibble(added_covariate = starting_covariate, cor_with_ammi = max(initial_W_cov))
    
    ## For each remaining covariate, find the next one that gives the highest gain (or lowest loss)
    while (length(available_covariates) > 0) {
      
      # Advance i
      i <- i + 1
      
      # Scan the available covariates, add to the matrix, calculate distance, compare to W_ammi
      new_W_cov <- map(available_covariates, ~c(used_covariates, .)) %>%
        map(~make_dist_mat(ec_tomodel_scaled_mat_use[, ., drop = FALSE])) %>%
        map_dbl(~cor(as.numeric(.), as.numeric(W_ammi)))
      
      old_cor <- max(new_W_cov)
      
      # Find the max difference
      which_max_diff <- which.max(new_W_cov - old_cor)
      
      # Get the covariate
      selected_covariate <- available_covariates[which_max_diff]
      
      ## Add to the list
      ec_ammi_cor_list[[i]] <- tibble(added_covariate = selected_covariate, cor_with_ammi = new_W_cov[which_max_diff])
      
      # Determine remaining covariates
      used_covariates <- union(used_covariates, selected_covariate)
      # Determine new available covariate vector
      available_covariates <- setdiff(covariates, used_covariates)
      
    }
    
    ## Collapse the list; calculate difference from previous
    ec_ammi_cor_df <- bind_rows(ec_ammi_cor_list) %>%
      mutate(difference = c(diff(cor_with_ammi), 0),
             stop = which.min(difference >= cor_tol),
             number = seq(nrow(.)))
    
    
    ## List of final covariates
    final_covariates <- subset(ec_ammi_cor_df, number <= stop, added_covariate, drop = TRUE)
    # Final distance matrix
    ec_dist_mat_final <- make_dist_mat(ec_tomodel_scaled_mat_use[, final_covariates, drop = FALSE])
    # final correlation
    final_cor <- subset(ec_ammi_cor_df, number == stop, cor_with_ammi, drop = TRUE)
    
    ## Return these results
    tibble(final_cor = final_cor, final_covariates = list(final_covariates), ec_dist_mat_final = list(ec_dist_mat_final),
           test_results = list(ec_ammi_cor_df))
    
  }) %>% ungroup()





### Create relationship matrices for environments ####

environmental_relmat_df <- ec_ammi_dist %>%
  # Add location main effect models
  left_join(., select(env_effect_models, -data)) %>%
  # Get a list of main effect models
  mutate(main_environment_covariates = map(model, ~attr(terms(formula(.)), "term.labels"))) %>%
  mutate_at(vars(main_environment_covariates, final_covariates), ~map(., ~ec_tomodel_scaled_mat[, .x, drop = FALSE] )) %>%
  # Calculate E mat for the main environmental covariates
  mutate(E_mat_main = map(main_environment_covariates, ~{
    if (ncol(.x) == 0) {
      d1 <- diag(x = nrow(.x))
      `dimnames<-`(x = d1, value = replicate(2, row.names(.x), simplify = FALSE))
    } else {
      Env_mat(x = .x, method = "Rincent2")
    } }) ) %>%
  select(trait, interaction_covariate_mat = final_covariates, main_covariate_mat = main_environment_covariates, 
         E_mat_main, E_mat_int = ec_dist_mat_final)



## Save these results
save("environmental_relmat_df", "ec_ammi_dist", "ec_tomodel_centered", "ec_tomodel_scaled", "ec_tomodel_centers", 
     file = file.path(result_dir, "ec_model_building.RData"))












































##################################
### Appendix
##################################



############################
### Visualization of daily weather data
############################


# ## Look at daily stats
# daily_ec_select <- growth_stage_weather %>%
#   inner_join(., env_trials) %>%
#   select(environment, dap, stage, mint:water_balance)
# 
# 
# ## Summarize max temperatures during grain fill
# grain_fill_maxt_summary <- daily_ec_select %>% 
#   filter(growth_stage == "grain_fill") %>%
#   ## Count number of days with 15-18 max temp,
#   ## < 15, 18-25 (moderate high), 25 - 30 (high), > 30 (very high)
#   mutate(grain_fill_condition = case_when(
#     maxt < 15 ~ "suboptimal",
#     between(maxt, 15, 18) ~ "optimal",
#     between(maxt, 18, 30) ~ "moderate_high",
#     between(maxt, 30, 35) ~ "high",
#     maxt > 35 ~ "very_high"
#   )) %>%
#   group_by(environment, grain_fill_condition) %>%
#   summarize(nDays = n())



# #### Example daily temperature for a drought, non-drought location-year
# 
# ## color vector for growth stage
# growth_stage_color <- setNames(neyhart_palette("umn2")[c(4, 5, 2)], names(growth_stage_rename))
# growth_stage_color["vegetative"] <- neyhart_palette()[5]
# 
# # what covariate to plot
# ec_to_plot <- "maxt"
# 
# daily_ec_example_plots <- daily_ec_select %>%
#   gather(covariate, value, -environment, -dap, -growth_stage) %>%
#   # Re-order growth stages
#   mutate(growth_stage = factor(growth_stage, levels = c("vegetative", "flowering", "grain_fill"))) %>%
#   filter(environment %in% c("EON16", "CRM16"), covariate == ec_to_plot) %>%
#   mutate(min_dap = min(dap), max_dap = max(dap)) %>%
#   ## calculate ranges for each covariate
#   group_by(covariate) %>%
#   mutate_at(vars(value), list(~min, ~max)) %>%
#   group_by(covariate, environment) %>%
#   do(plot = {
#     df <- .
#     
#     ## Create segments for growth stages
#     gs_seg <- df %>% 
#       group_by(growth_stage) %>% 
#       summarize(start = min(dap) - 1, end = max(dap))
#     
#     ec_name <- covariate_variable_rename[unique(df$covariate)]
#     ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
#     
#     # x axis limits
#     x_limit <- c(df$min_dap[1] - 1, df$max_dap[1])
#     
#     # y axis limits
#     y_limit <- c(df$min[1], df$max[1])
#     y_end <- quantile(y_limit, 0.75)
#     
#     ## Plot
#     ggplot(data = df, aes(x = dap, y = value)) +
#       geom_segment(data = gs_seg, mapping = aes(x = start, xend = end, y = y_end, yend = y_end, color = growth_stage), lwd = 10) +
#       geom_line() +
#       scale_y_continuous(name = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), breaks = pretty, limits = y_limit) +
#       scale_x_continuous(breaks = pretty, name = "Days after planting", limits = x_limit) +
#       scale_color_manual(values = growth_stage_color, name = "Predicted\ngrowth stage",
#                          labels = function(x) str_to_title(str_remove_all(x, "_"))) +
#       labs(subtitle = unique(df$environment)) +
#       theme_presentation2() +
#       theme(legend.position = "bottom")
#       
#   })
# 
#   
# # Extract sample plots
# plot_list <- daily_ec_example_plots %>%
#   pull(plot) %>%
#   map(~. + theme(legend.position = "none", axis.title = element_blank()))
# plot_list[[1]] <- plot_list[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# # Merge
# select_ec_plot_merge <- plot_grid(plotlist = plot_list, ncol = 1)
# 
# # y_axis label
# ec_name <- covariate_variable_rename[ec_to_plot]
# ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
# y_label <- grid::textGrob(label = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), rot = 90) 
# 
# x_label <- grid::textGrob(label = "Days after planting")
# 
# select_ec_plot_merge1 <- plot_grid(y_label, select_ec_plot_merge, nrow = 1, rel_widths = c(0.05, 1))  %>%
#   plot_grid(., x_label, get_legend(daily_ec_example_plots$plot[[1]]),ncol = 1, rel_heights = c(1, 0.05, 0.2)) 
# 
# 
# ggsave(filename = "sample_maxt_growth_stage.jpg", plot = select_ec_plot_merge1, path = fig_dir,
#        height = 6, width = 7, dpi = 1000)




# # what covariate to plot
# ec_to_plot <- "water_stress"
# 
# daily_ec_example_plots <- daily_ec_select %>%
#   gather(covariate, value, -environment, -dap, -growth_stage) %>%
#   # Re-order growth stages
#   mutate(growth_stage = factor(growth_stage, levels = c("vegetative", "flowering", "grain_fill"))) %>%
#   filter(environment %in% c("EON16", "CRM16"), covariate == ec_to_plot) %>%
#   mutate(min_dap = min(dap), max_dap = max(dap)) %>%
#   ## calculate ranges for each covariate
#   group_by(covariate) %>%
#   mutate_at(vars(value), list(~min, ~max)) %>%
#   group_by(covariate, environment) %>%
#   do(plot = {
#     df <- .
#     
#     ## Create segments for growth stages
#     gs_seg <- df %>% 
#       group_by(growth_stage) %>% 
#       summarize(start = min(dap) - 1, end = max(dap))
#     
#     ec_name <- covariate_variable_rename[unique(df$covariate)]
#     ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
#     
#     # x axis limits
#     x_limit <- c(df$min_dap[1] - 1, df$max_dap[1])
#     
#     # y axis limits
#     y_limit <- c(df$min[1], df$max[1])
#     y_end <- quantile(y_limit, 0.75)
#     
#     ## Plot
#     ggplot(data = df, aes(x = dap, y = value)) +
#       geom_segment(data = gs_seg, mapping = aes(x = start, xend = end, y = y_end, yend = y_end, color = growth_stage), lwd = 10) +
#       geom_line() +
#       scale_y_continuous(name = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), breaks = pretty, limits = y_limit) +
#       scale_x_continuous(breaks = pretty, name = "Days after planting", limits = x_limit) +
#       scale_color_manual(values = growth_stage_color, name = "Predicted\ngrowth stage",
#                          labels = function(x) str_to_title(str_remove_all(x, "_"))) +
#       labs(subtitle = unique(df$environment)) +
#       theme_presentation2() +
#       theme(legend.position = "bottom")
#     
#   })
# 
# 
# # Extract sample plots
# plot_list <- daily_ec_example_plots %>%
#   pull(plot) %>%
#   map(~. + theme(legend.position = "none", axis.title = element_blank()))
# plot_list[[1]] <- plot_list[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# # Merge
# select_ec_plot_merge <- plot_grid(plotlist = plot_list, ncol = 1)
# 
# # y_axis label
# ec_name <- covariate_variable_rename[ec_to_plot]
# ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
# y_label <- grid::textGrob(label = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), rot = 90) 
# 
# x_label <- grid::textGrob(label = "Days after planting")
# 
# select_ec_plot_merge1 <- plot_grid(y_label, select_ec_plot_merge, nrow = 1, rel_widths = c(0.05, 1))  %>%
#   plot_grid(., x_label, get_legend(daily_ec_example_plots$plot[[1]]),ncol = 1, rel_heights = c(1, 0.05, 0.2)) 
# 
# 
# ggsave(filename = "sample_maxt_growth_stage.jpg", plot = select_ec_plot_merge1, path = fig_dir,
#        height = 6, width = 7, dpi = 1000)



##########################
## Investiate time-dependent stressors
##########################

# ## Look at "number of days above some threshold temperature x during grainfill" with x = 15, 16, ..., max(T)
# # Determine the threshold that explains the most GxE variance
# temp_thresh <- seq(20, 40)
# 
# grain_fill_temperature_stress <- map(temp_thresh, ~{
#   growth_stage_weather %>%
#     filter(growth_stage == "grain_fill") %>%
#     group_by(environment) %>%
#     summarize(stress_days = sum(maxt > .x), threshold = .x) 
# })
# 
# ## format s2met data for modeling
# s2_met_tomodel <- S2_MET_BLUEs %>%
#   filter(line_name %in% tp) %>%
#   filter(trait == "GrainYield") %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(data = map(data, ~mutate(., line_name = as.factor(line_name),
#                                   line_name = `contrasts<-`(line_name, value = `colnames<-`(contr.sum(levels(line_name)), head(levels(line_name), -1))))
#   ))
# 
# 
# ## fit models
# s2_met_tomodel_stress_days <- s2_met_tomodel %>%
#   crossing(., stress = grain_fill_temperature_stress) %>%
#   mutate(data = map2(data, stress, ~left_join(.x, .y, by = "environment")),
#          threshold = map_dbl(stress, ~unique(.$threshold)))
# 
# 
# ## Extract p-values from anova
# stress_days_fit1 <- s2_met_tomodel_stress_days %>%
#   group_by(trait, threshold) %>%
#   do(fit = lm(value ~ line_name + line_name:stress_days, data = .$data[[1]])) %>%
#   ungroup() %>%
#   mutate(pvalue = map_dbl(fit, ~as.data.frame(anova(.))[2,5]))
# 
# 
# plot(-log10(pvalue) ~ threshold, stress_days_fit1)
# 
# ## Look at regression from lowest pvalue
# best_fit <- stress_days_fit1 %>%
#   filter(pvalue == min(pvalue, na.rm = TRUE)) %>%
#   pull(fit)
# 
# plot(best_fit[[1]])
# 


# ## Plot environment and maxt during grain fill
# daily_ec_select %>% 
#   filter(growth_stage == "grain_fill") %>%
#   mutate(environment = factor(environment, levels = ec_select$environment[order(ec_select$grain_fill_maxt, decreasing = TRUE)])) %>%
#   ggplot(aes(x = environment, y = maxt)) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 35, ymax = 40, fill = "heat stress"), alpha = 0.5) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 35, fill = "high temperature"), alpha = 0.5) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30, fill = "moderate high temperature"), alpha = 0.5) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 15, ymax = 20, fill = "optimal temperature"), alpha = 0.5) +
#   geom_violin(fill = alpha("white", 0.5)) +
#   # geom_boxplot(fill = NA) +
#   scale_y_continuous(breaks = pretty) +
#   scale_fill_manual(values = rev(viridis::inferno(direction = -1, n = 10)[1:4]), name = "Temperature stress level",
#                     guide = guide_legend(title.position = "top")) +
#   labs(subtitle = "Maximum temperature stress during grain fill") +
#   theme_presentation2(base_size = 12) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
# ggsave(filename = "grain_fill_max_temp_stress.jpg", path = fig_dir, width = 8, height = 5, dpi = 300)
# 

 







############################
### Model covariates that explain environment effect
############################

# ## Fit a base, random effect model for all traits
# base_model_fit <- s2_met_tomodel %>%
#   group_by(trait) %>%
#   do(fit = lmer(value ~ 1 + (1|line_name) + (1|environment), data = .)) %>%
#   ungroup()
# 
# 
# 
# 
# ## Get estimates of environmental effect
# base_model_env_effect <- base_model_fit %>%
#   mutate(effect = map(fit, ~ranef(.)$environment %>% rownames_to_column("environment") %>% rename(effect = 2))) %>%
#   unnest(effect)
#   
# ## Histogram
# par(mfrow = c(2,3))
# for (tr in unique(base_model_env_effect$trait)) {
#   hist(subset(base_model_env_effect, trait == tr, effect, drop = T), main = tr)
# }
# par(mfrow = c(1,1))
# 
# 
# # # environments as fixed
# # base_model_env_effect_alt <- s2_met_tomodel %>%
# #   group_by(trait) %>%
# #   do(fit = lmer(value ~ 1 + (1|line_name) + environment, data = ., contrasts = list(environment = "contr.sum"))) %>%
# #   ungroup() %>%
# #   mutate(mf = map(fit, model.frame),
# #          effect = map2(fit, mf, ~tibble(environment = levels(.y$environment), effect = c(fixef(.x)[-1], -sum(fixef(.x)[-1]))))) %>%
# #   unnest(effect)
# # 
# # ## Compare
# # left_join(base_model_env_effect, base_model_env_effect_alt, by = c("trait", "environment")) %>%
# #   qplot(x = effect.x, y = effect.y, data = .) +
# #   facet_wrap(~ trait, scales = "free") +
# #   geom_abline(slope = 1, intercept = 0)
# # 
# # # Good! - proceed with random effect
# 
# 
# 
# 
# ## Use stepwise regression to find the best model
# env_effect_to_model <- base_model_env_effect %>%
#   left_join(s2_met_tomodel)
# 
# # Group by trait
# env_effect_models <- env_effect_to_model %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(model = list(NULL))
# 
# for (i in seq(nrow(env_effect_models))) {
#   
#   df <- env_effect_models$data[[i]]
#   tr <- env_effect_models$trait[i]
#   
#   ## Remove some covariates depending on when the trait is measured
#   df1 <- select(df, environment, effect, subset(ec_by_trait, trait == tr, covariates, drop = T)[[1]])
#   
#   # Minimal model
#   min_model <- effect ~ 1
#   # Max model
#   max_model <- add_predictors(min_model, as.formula(paste0("~", paste0(names(df1)[-1:-2], collapse = " + "))))
#   
#   # Stepwise regression
#   fit_base <- lm(formula = min_model, data = df1)
#   fit_step <- step(object = fit_base, scope = max_model, direction = "both")
#   
#   ## Return the model
#   env_effect_models$model[[i]] <- fit_step
#   
# }
# 
# 
# 
# ## Plot
# env_effect_models %>% 
#   mutate(data = map2(data, model, ~mutate(.x, predicted_effect = predict(.y)))) %>%
#   unnest(data) %>%
#   ggplot(aes(x = predicted_effect, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   facet_wrap(~ trait, scales = "free") +
#   theme_acs()






# ############################
# ### Model covariates that explain GxE
# ############################
# 
# 
# 
# #####
# # Model building
# #####
# 
# 
# ## Data.frame to store final models and covariates
# ec_model_building <- s2_met_tomodel %>% 
#   distinct(trait) %>%
#   mutate(all_model = list(NULL), apriori_model = list(NULL), final_model = list(NULL))
# 
# 
# 
# 
# 
# 
# ## 
# ## Grain yield, grain protein, and test weight
# ## 
# ## Start with the following a prior covariates:
# ## 
# ## tmin during flowering
# ## radiation during vegetative
# ## water stress during flowering
# ## grain fill water stress
# ## grain fill tmax
# ## 
# ## 
# ## For heading date and plant height, just use all covariates and backwards elimination
# ## 
# ## 
# 
# 
# ## Control
# control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
# 
# 
# ## Loop over traits
# for (i in seq(nrow(ec_model_building))) {
#   
#   # Get the trait
#   tr <- ec_model_building$trait[i]
#   
#   # Subset the data
#   tr_data <- filter(s2_met_tomodel, trait == tr)
#   
#   # Get the base model
#   tr_base_lmer <- lmer(value ~ 1 + (1|line_name) + environment, data = tr_data)
#   tr_base_lm <- lm(value ~ 1 + line_name + environment, data = tr_data)
#   
#   
#   
#   # ### Apriori model ###
#   # 
#   # ## Assign a priori covariates
#   # apriori_tr_covariates <- subset(ec_by_trait_apriori, trait == tr, covariates, drop = T)[[1]]
#   # 
#   # ## Create a new formula
#   # # tr_covariate_form1 <- as.formula(paste0("~", paste0(apriori_tr_covariates, collapse = " + "), " + ", 
#   # #                                        paste0("line_name:", apriori_tr_covariates, collapse = " + ")))
#   # 
#   # # Diagonal covariance structure
#   # tr_covariate_form2 <- as.formula(paste0("~ ", paste0("(0 +", apriori_tr_covariates, "|line_name)", collapse = " + ")))
#   # # # Unstructured covariance
#   # # tr_covariate_form2 <- as.formula(paste0("~ (0 +", paste0(apriori_tr_covariates, collapse = " + "), "|line_name)"))
#   # 
#   # # gxe_formula1 <- add_predictors(formula(tr_base_model), tr_covariate_form1)
#   # gxe_formula2 <- add_predictors(formula(tr_base_model), tr_covariate_form2)
#   # 
#   # 
#   # # Refit the model
#   # # tr_ec_model_mixed1 <- lmer(formula = gxe_formula1, data = tr_data)
#   # tr_ec_model_mixed2_apriori <- lmer(formula = gxe_formula2, data = tr_data, control = control)
#   # 
#   # 
#   # # Stepwise elimination of extract covariates
#   # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed2_apriori)
#   # # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed, reduce.random = FALSE, alpha.fixed = alpha)
#   # 
#   # # Get the model
#   # tr_ec_model_mixed_apriori_final <- get_model(tr_ec_model_mixed_step)
# 
# 
#   
#   
#   ### All EC model ###
#   # assign covariates
#   all_tr_covariates <- subset(ec_by_trait_all, trait == tr, covariates, drop = T)[[1]]
#   
#   ## Create a new formula
#   # LMER
#   # tr_covariate_form <- as.formula(paste0("~ (0 +", paste0(all_tr_covariates, collapse = " + "), "|line_name)"))
#   tr_covariate_form_lmer <- as.formula(paste0("~ ", paste0("(0 +", all_tr_covariates, "|line_name)", collapse = " + ")))
#   gxe_formula_lmer <- add_predictors(formula(tr_base_lmer), tr_covariate_form_lmer)
#   
#   # Refit the model
#   tr_ec_model_mixed_all <- lmer(formula = gxe_formula_lmer, data = tr_data, control = control)
#   
#   # Stepwise elimination of extra covariates
#   tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed_all)
#   # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed_all, reduce.random = FALSE, alpha.fixed = alpha)
#   
#   # Get the model
#   tr_ec_model_mixed_all_final <- get_model(tr_ec_model_mixed_step)
#   
#   
#   
#   
#   # LM
#   tr_covariate_form_lm <- as.formula(paste0("~ ", paste0(all_tr_covariates, ":line_name", collapse = " + ")))
#   gxe_formula_lm <- add_predictors(formula(tr_base_lm), tr_covariate_form_lm)
#   
#   ## Stepwise addition of covariates
#   tr_ec_model_fixed_step <- step(object = tr_base_lm, scope = gxe_formula_lm, direction = "both")
#   
#   test <- lm(value ~ 1 + line_name + environment + line_name:grain_fill.radn_sum, tr_data)
#   
#   
#   # Refit the model
#   tr_ec_model_fixed_all <- lm(formula = gxe_formula_lm, data = tr_data)
#   
#   # Stepwise elimination of extra covariates
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
#   ## Add results to the tibble
#   ec_model_building$apriori_model[i] <- list(tr_ec_model_mixed2_apriori)
#   ec_model_building$final_model[i] <- list(tr_ec_model_mixed_apriori_final)
#   ec_model_building$all_model[i] <- list(tr_ec_model_mixed_all_final)
#   
# }
# 
# 
# ## Tidy the output
# ec_model_building_tidy <- ec_model_building %>%
#   gather(model, object, -trait) %>%
#   mutate(model = str_remove(model, "_model"),
#          # Calculate varprop
#          varprop = map(object, ~as.data.frame(VarCorr(.)) %>% 
#                          mutate(term = case_when(is.na(var1) ~ grp, var1 == "(Intercept)" ~ grp, TRUE ~ var1), 
#                                 variance = vcov, varprop = variance / sum(variance)) ),
#          # Calculate statistics
#          AIC = map_dbl(object, AIC), BIC = map_dbl(object, BIC), logLik = map_dbl(object, logLik))




# 
# 
# ############################
# ### Heritability
# ############################
# 
# 
# ## Fit genomic heritability models for slopes
# # Extract the coefficients for each model
# # Keep final and apriori models
# ec_interaction_coef <- covariate_reg_coefs %>%
#   rename_all(~str_remove(., "_model")) %>%
#   gather(model, out, -trait) %>% 
#   unnest(out)
# 
# 
# # Subset the K matrix
# K_use <- K[tp_geno, tp_geno]
# 
# ## Fit models
# ec_interaction_coef_herit <- ec_interaction_coef %>%
#   group_by(trait, model, covariate) %>%
#   do({
#     df <- .
#     
#     df1 <- subset(df, line_name %in% tp_geno)
#     
#     ## Calculate heritability
#     invisible(capture.output(herit_fit <- marker_h2(data.vector = df1$estimate, geno.vector = df1$line_name, K = K_use, alpha = alpha)))
#     
#     ## Return df
#     tibble(heritability = herit_fit$h2, lower = herit_fit$conf.int1[1], upper = herit_fit$conf.int1[2])
#     
#   })
# 
# ## Write table
# write_csv(x = ec_interaction_coef_herit, path = file.path(fig_dir, "covariate_slope_heritability.csv"))
# 


