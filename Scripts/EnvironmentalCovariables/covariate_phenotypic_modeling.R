## S2MET phenotypic modeling with covariates
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## Model the impact of covariates on GxE and LxE
## 
## This script will look at variable selection to determine optimal models that use
## covariates
## 

# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))

library(paletteer)
library(cowplot)


## Add new packages to load
pkgs <- union(pkgs, c("modelr", "broom", "lme4", "car", "patchwork", "caret"))
# Load these packages
invisible(lapply(X = pkgs, library, character.only = TRUE))

## significance level
alpha <- 0.05



# Load data ---------------------------------------------------------------

# Load covariates for environments and historical covariates
load(file.path(result_dir, "concurrent_historical_covariables_nasapower.RData"))
# load(file.path(result_dir, "concurrent_historical_covariables_daymet.RData"))


## Filter BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp,
         environment %in% train_test_env) %>% # only model train/test environments
  mutate_at(vars(line_name, environment), as.factor)


# Modeling ----------------------------------------------------------------


# Concurrent environmental covariates =====================================


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = names(ec_tomodel_centered)[-1]) %>%
  filter(
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 



# Use the following covariates for all traits (except where timing would make
# it inappropriate):
# 
# radiation during vegetative
# tmin during flowering
# water stress during flowering
# grain fill water stress
# grain fill tmax
# 


## A priori covariates - quite minimal; each should have a citation 
## 
## HeadingDate - everything before flowering
## 
## PlantHeight - everything before grain fill
##  
## Grain yield, grain protein, and test weight
## - Everything for plant height plus...
## - Elevated temperature during grain fill (Passarella et al 2005)
## - Drought during grain fill - (Savin and Nicholas  1996)
## 
## 
## 

apriori_covariate_df <- trait_covariate_df %>%
  filter(str_detect(covariate, "soil", negate = TRUE)) %>%
  {bind_rows(
    filter(., trait %in% c("HeadingDate", "PlantHeight")),
    filter(., trait %in% c("GrainProtein", "GrainYield", "TestWeight"),
           covariate %in% c(subset(., trait == "PlantHeight", covariate, drop = TRUE), 
                            "grain_fill.tmean_mean", "grain_fill.water_balance_sum"))
  )}


## Factorial regression with AIC stepwise selection ##

## Group by trait and model
concurrent_fact_reg <- S2_MET_BLUEs_tomodel %>%
  group_by(trait) %>%
  do({

    df <- .
    # Factorize
    df1 <- df %>%
      # filter(location != "Aberdeen") %>%
      droplevels() %>%
      left_join(., ec_tomodel_centered, by = "environment") %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    
    # Apriori
    ## Add covariates - filter
    covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    apriori_out <- fact_reg(data = df1, covariates = covariates_use, method = "apriori")
    


    # Ad hoc
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step", criterion = "BIC")


    # Ad hoc - without soil
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T),
                             covariate, drop = TRUE)
    adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step")


    ## Return results
    tibble(model = c("base", "base_alt", "model2", "model3"),
           apriori = apriori_out,
           adhoc = adhoc_out,
           adhoc_nosoil = adhoc_nosoil_out)

  }) %>% ungroup()


## Assess via LOO cv
concurrent_fact_reg_loo_test <- S2_MET_BLUEs_tomodel %>%
  left_join(., ec_tomodel_centered) %>%
  group_by(trait) %>%
  do(crossv_loo_grouped(group_by(., environment))) %>%
  ungroup() %>%
  left_join(., concurrent_fact_reg) %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  filter(trait == "TestWeight") %>%
  mutate(train_fit = map2(.x = train, .y = adhoc, ~update(.y, data = .x)),
         test_pred = map2(.x = test, .y = train_fit, ~add_predictions(as.data.frame(.x), .y)))

concurrent_fact_reg_loo_test %>% 
  unnest(test_pred) %>%
  qplot(x = pred, y = value, color = environment, data = ., facets = ~model)


## Prepare results for saving
concurrent_fact_reg_feature_selection <- concurrent_fact_reg %>%
  mutate(feat_sel_type = "stepAIC", direction = "forward") %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  mutate_at(vars(apriori, adhoc, adhoc_nosoil), ~map(., ~list(optVariables = attr(terms(.x), "term.labels"))))
  
  




# Feature selection ############################################################

## Use some feature selection procedures to identify covariates
## 
## Wrap the feature selection within cross-validation to avoid selection
## bias
## 
## Use lm with CV - implemented through recursive feature elimination
## 

concurrent_feature_selection_list <- S2_MET_BLUEs_tomodel %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(concurrent_feature_selection_list))) {
    
    df <- concurrent_feature_selection_list$data[[i]] %>%
      mutate(trait = concurrent_feature_selection_list$trait[i])
    # Factorize
    df1 <- df %>%
      # filter(location != "Aberdeen") %>%
      droplevels() %>%
      left_join(., ec_tomodel_centered, by = "environment") %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    

    loo_indices <- df1 %>%
      group_by(environment) %>%
      crossv_loo_grouped() %>%
      pull(train) %>%
      map("idx")
    
    ## Recursive feature addition
    rfa_out_df <- rfa_proc()

    
    ## Recursive feature elimination using PLS
    
    # 1. Estimate environmental means
    # main_fit <- lm(value ~ line_name + environment, data = df1, weights = 1 / (std_error^2))
    main_fit <- lm(value ~ line_name + environment, data = df1)
    rfs_out_df <- rfe_proc(main_fit = main_fit, env.col = "environment")
    
    concurrent_feature_selection_list$out[[i]] <- bind_rows(rfa_out_df, rfs_out_df)
    
}

concurrent_feature_selection <- unnest(concurrent_feature_selection_list, out)




# Historical covariates =====================================================


## Ideas
## 1. Calculate locations means and use that as the input (and validation)
## 2. Use all data and fit a model with random year, location:year, and g:year
## 3. Use all data and fit the g + location model


## Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp)


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = names(historical_ec_tomodel_timeframe_centered$time_frame5_2010_2014)[-1:-2]) %>%
  filter(
    covariate != "awc_range",
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 

## Select the historical covariate data
historical_ec_tomodel_centered_use <- historical_ec_tomodel_timeframe_centered$time_frame5_2010_2014
  

# Use the following covariates for all traits (except where timing would make
# it inappropriate):
# 
# radiation during vegetative
# tmin during flowering
# water stress during flowering
# grain fill water stress
# grain fill tmax


# # Load packages
# invisible(lapply(X = pkgs, library, character.only = TRUE))

historical_fact_reg <- S2_MET_loc_BLUEs_tomodel %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(historical_fact_reg))) {

    df <- historical_fact_reg$data[[i]]
    df$trait <- historical_fact_reg$trait[i]
    # Factorize
    df1 <- df %>%
      droplevels() %>%
      left_join(., historical_ec_tomodel_centered_use, by = "location") %>%
      mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))

    # Apriori
    ## Add covariates - filter
    covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "apriori")

    # Ad hoc
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")

    # Ad hoc - without soil
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T),
                             covariate, drop = TRUE)
    adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "location", method = "step")


    ## Return results
    historical_fact_reg$out[[i]] <- tibble(model = c("base", "base_alt", "model4", "model5"),
           apriori = apriori_out,
           adhoc = adhoc_out,
           adhoc_nosoil = adhoc_nosoil_out)

}


historical_fact_reg <- unnest(historical_fact_reg, out)


## Assess via LOO cv
historical_fact_reg_loo_test <- S2_MET_loc_BLUEs_tomodel %>%
  left_join(., historical_ec_tomodel_centered_use) %>%
  group_by(trait) %>%
  do(crossv_loo_grouped(group_by(., location))) %>%
  ungroup() %>%
  left_join(., historical_fact_reg) %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  filter(trait == "TestWeight") %>%
  mutate(train_fit = map2(.x = train, .y = adhoc_nosoil, ~update(.y, data = .x)),
         test_pred = map2(.x = test, .y = train_fit, ~add_predictions(as.data.frame(.x), .y)))

historical_fact_reg_loo_test %>% 
  unnest(test_pred) %>%
  qplot(x = pred, y = value, color = location, data = .) +
  facet_wrap(~ model, scales = "free")


## Prepare results for saving
historical_fact_reg_feature_selection <- historical_fact_reg %>%
  mutate(feat_sel_type = "stepAIC", direction = "forward") %>%
  filter(str_detect(model, "base", negate = TRUE)) %>%
  mutate_at(vars(apriori, adhoc, adhoc_nosoil), ~map(., ~list(optVariables = attr(terms(.x), "term.labels"))))






# Feature selection ############################################################

historical_feature_selection_list <- S2_MET_loc_BLUEs_tomodel %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(historical_feature_selection_list))) {

    df <- historical_feature_selection_list$data[[i]] %>% 
      mutate(trait = historical_feature_selection_list$trait[i])
    # Factorize
    df1 <- df %>%
      arrange(location, line_name) %>%
      # filter(location != "Aberdeen") %>%
      droplevels() %>%
      left_join(., historical_ec_tomodel_centered_use, by = "location") %>%
      mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
    
    loo_indices <- df1 %>%
      group_by(location) %>%
      crossv_loo_grouped() %>%
      pull(train) %>%
      map("idx")
    
    ## Recursive feature addition
    rfa_out_df <- rfa_proc(env.col = "location")
    
    
    ## Recursive feature elimination using PLS
    
    # 1. Estimate environmental means
    main_fit <- lm(value ~ line_name + location, data = df1)
    rfs_out_df <- rfe_proc(main_fit = main_fit, env.col = "location")
    
    historical_feature_selection_list$out[[i]] <- bind_rows(rfa_out_df, rfs_out_df)
    
}

historical_feature_selection <- unnest(historical_feature_selection_list, out)


save("concurrent_feature_selection", "historical_feature_selection", 
     "concurrent_fact_reg_feature_selection", "historical_fact_reg_feature_selection",
     file = file.path(result_dir, "feature_selection_results_nasapower.RData"))











# Analyze results ---------------------------------------------------------

# Load the data
load(file = file.path(result_dir, "feature_selection_results_nasapower.RData"))

ec_tomodel_scaled_mat <- ec_tomodel_scaled %>%
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()

## Combine concurrent feature selection df
concurrent_features <- bind_rows(
  concurrent_fact_reg_feature_selection %>% select(trait, model, apriori, stepAIC = adhoc) %>% 
    gather(feat_sel_type, features, apriori, stepAIC),
  select(concurrent_feature_selection, trait, model, feat_sel_type, features = adhoc)
)

## combine model2 and model3 covariates
concurrent_features1 <- concurrent_features %>%
  filter(feat_sel_type != "rfs_pls") %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model2, model3, union)) %>%
  ## Add all covariates
  add_row(trait = unique(.$trait), feat_sel_type = "all", 
          features = list(c(names(ec_tomodel_centered)[-1], paste0("line_name:", names(ec_tomodel_centered)[-1])))) %>%
  mutate(features = map(features, ~setdiff(., "line_name"))) %>%
  select(-contains("model")) %>%
  mutate(interaction_features = map(features, ~str_subset(., ":")),
         main_features = map2(features, interaction_features, setdiff))
  


## How many covariates per trait and variable selection procedure?
(concurrent_features2 <- concurrent_features1 %>%
  mutate_at(vars(contains("features")), ~map_dbl(., length)) %>%
  arrange(trait, feat_sel_type))


# trait        feat_sel_type features interaction_features main_features
# 1 GrainProtein all                 82                   41            41
# 2 GrainProtein apriori             40                   20            20
# 3 GrainProtein rfa_cv               9                    0             9
# 4 GrainProtein stepAIC             10                    2             8
# 5 GrainYield   all                 82                   41            41
# 6 GrainYield   apriori             40                   20            20
# 7 GrainYield   rfa_cv              12                    0            12
# 8 GrainYield   stepAIC             28                    2            26
# 9 HeadingDate  all                 82                   41            41
# 10 HeadingDate  apriori             24                   12            12
# 11 HeadingDate  rfa_cv               4                    0             4
# 12 HeadingDate  stepAIC             19                    0            19
# 13 PlantHeight  all                 82                   41            41
# 14 PlantHeight  apriori             36                   18            18
# 15 PlantHeight  rfa_cv               5                    0             5
# 16 PlantHeight  stepAIC             24                    1            23
# 17 TestWeight   all                 82                   41            41
# 18 TestWeight   apriori             40                   20            20
# 19 TestWeight   rfa_cv               5                    0             5
# 20 TestWeight   stepAIC             16                    3            13

# Save as csv
write_csv(x = concurrent_features2, path = file.path(fig_dir, "concurrent_environmental_covariate.csv"))


## Determine overlap between select pairs of covariates
## 
## apriori vs rfa_cv
## apriori vs stepAIC
## rfa_cv vs stepAIC
## 
covariate_feature_type_comparison <- crossing(concurrent_features1, concurrent_features1) %>% 
  filter(trait == trait1, 
         feat_sel_type != feat_sel_type1,
         feat_sel_type %in% c("apriori", "rfa_cv"), 
         feat_sel_type1 %in% c("rfa_cv", "stepAIC")) %>%
  select(-trait1, -starts_with("features"))


# Create a table and export
covariate_feature_type_comparison_table <- covariate_feature_type_comparison %>%
  mutate(compare_interaction = map2(interaction_features, interaction_features1, intersect),
         compare_main = map2(main_features, main_features1, intersect)) %>%
  select(trait, starts_with("feat_sel_type"), starts_with("compare")) %>%
  gather(term, comparison, contains("compare")) %>%
  unnest(comparison) %>%
  mutate(term = str_remove(term, "compare_"))
write_csv(x = covariate_feature_type_comparison_table, path = file.path("concurrent_feature_selection_comparison.csv"))




## Calculate environmental relationship matrices based on these features
concurrent_features_env_relmat <- concurrent_features1 %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate(interaction_features_ec_mat = map(interaction_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE]),
         main_features_ec_mat = map(main_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE])) %>%
  mutate_at(vars(ends_with("ec_mat")), ~map(., ~Env_mat(x = .x, method = "Jarq")))

concurrent_features_env_relmat1 <- concurrent_features_env_relmat %>%
  select(trait, feat_sel_type, interaction = interaction_features_ec_mat, main = main_features_ec_mat) %>%
  gather(covariate_type, Emat, main, interaction) %>%
  ## If a matrix is all NA, convert to diagonal
  mutate(Emat = modify_if(Emat, ~all(is.na(.)), ~`diag<-`(ifelse(is.na(.), 0, 1), 1)))


## Compare the relationship between training/test and external environments for 
## each trait and feature selection type
concurrent_features_env_relmat1 %>%
  mutate(train_val_relat = map_dbl(Emat, ~mean(.[train_test_env, validation_env])),
         train_test_relat = ) %>% 
  arrange(covariate_type, trait, feat_sel_type) %>% 
  select(-Emat) %>% 
  View



# Define heat colors
heat_colors <- wesanderson::wes_palette("Zissou1", n = 5)[c(1,3,5)]


## Plot heatmaps of relationship matrices
# Separate plots by main/int covariate types
concurrent_features_heatmap_plots <- concurrent_features_env_relmat1 %>%
  group_by(trait, feat_sel_type, covariate_type) %>%
  do(plot = {
    row <- .
    
    # Create the heatmap to order the environments
    row_heat <- heatmap(x = row$Emat[[1]])
    
    # Factor order of environments
    env_order <- fct_inorder(row.names(row$Emat[[1]])[row_heat$rowInd])
    
    # Create the plotting data.frame
    dat <- row$Emat[[1]] %>%
      as.data.frame(.) %>% 
      rownames_to_column(., "environment") %>%
      gather(environment2, relationship, -environment) %>%
      # Refactor the environments
      mutate_at(vars(contains("environment")), ~factor(., levels = levels(env_order))) %>%
      mutate_at(vars(contains("environment")), list(group = ~ifelse(. %in% train_test_env, "training", "external"))) %>%
      mutate(environment_group = factor(environment_group, levels = c("training", "external")))
    
    # Plot
    g_heat <- dat %>%
      ggplot(aes(x = environment, y = environment2, fill = relationship)) +
      geom_tile() +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3]) +
      facet_grid(environment2_group ~ environment_group, scales = "free", space = "free") +
      labs(subtitle = paste(str_add_space(row$trait), row$feat_sel_type, sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = "none")
      
    # Return the plot
    g_heat
    
  }) %>% ungroup()

## Create plots by trait
for (plotList in split(concurrent_features_heatmap_plots, concurrent_features_heatmap_plots$trait)) {
  # Create plot
  plot_to_save <- plot_grid(plotlist = plotList$plot, nrow = n_distinct(plotList$feat_sel_type),
                            labels = c("Int.", "Main"), label_x = 0)
  
  # File name
  filename <- paste0("environment_covariate_relationship_", unique(plotList$trait), ".jpg")
  ggsave(filename = filename, plot = plot_to_save, path = fig_dir, 
         width = 8, height = 14, dpi = 1000)
  
}







## Analyze historical covariates ##

# Subset the data
historical_fact_reg_sample_rmse <- historical_fact_reg_sample %>%
  filter(model == "model3") %>% # Skip the base model
  select(-test_predictions, -covariates)

# Plot RMSE
g_historical_rmse <- historical_fact_reg_sample_rmse %>%
  # filter(!(trait == "GrainProtein" & rmse > 300),
  #        !(trait == "TestWeight" & rmse > 2000)) %>%
  ggplot(aes(x = selection, y = rmse, color = selection)) +
  geom_boxplot() +
  facet_wrap(~ trait, scales = "free_y", nrow = 1) +
  theme_genetics(base_size = 10)

# Remove outliers
g_historical_rmse1 <- g_historical_rmse %>%
  modify_at("data", ~filter(., !(trait == "GrainProtein" & rmse > 30), !(trait == "GrainYield" & rmse > 20000),
                            !(trait == "HeadingDate" & rmse > 30),
                            !(trait == "PlantHeight" & rmse > 10000), !(trait == "TestWeight" & rmse > 50000)))

## Combine plots
g_historical_rmse + g_historical_rmse1 + plot_layout(ncol = 1)


# Plot accuracy
historical_fact_reg_sample_rmse %>%
  ggplot(aes(x = selection, y = acc, color = selection)) +
  geom_boxplot() +
  facet_wrap(~ trait, scales = "free_y") +
  theme_genetics(base_size = 10)

## Summarize both
historical_fact_reg_sample_summary <- historical_fact_reg_sample_rmse %>% 
  gather(stat, value, rmse, acc) %>% 
  group_by(trait, model, selection, stat) %>% 
  do({
    x <- .$value
    stats <- boxplot.stats(x)
    tibble(median = median(x), lower = stats$conf[1], upper = stats$conf[2])
  }) %>% ungroup()

historical_fact_reg_sample_summary %>%
  arrange(stat) %>%
  as.data.frame()



## Compare covariates selected in each model with those selected using the 
## full dataset

# Prepare the full data results
full_data_fact_reg <- fr_var_summary %>%
  filter(timeframe == "historical", selection != "apriori") %>%
  mutate(full_covariates = map(var_prop_summary, ~subset(., ! term %in% c("line_name", "Residuals"), term, drop = TRUE))) %>%
  select(-var_prop_summary)


# How often do we get the same exact model?
# How often do we get the same full-data covariate
# How often do we get a new covariate?
historical_fact_reg_sample_covariate_check <- historical_fact_reg_sample %>%
  filter(model == "model3") %>%
  select(trait, dropped_group, selection, covariates) %>%
  inner_join(., full_data_fact_reg) %>%
  # Calculate the proportion of full-data covariates that were recovered
  mutate(full_cov_recovered = map2(full_covariates, covariates, intersect),
         prop_full_cov_recovered = map2_dbl(full_covariates, covariates, ~mean(.x %in% .y)),
         new_covariates = map2(covariates, full_covariates, setdiff))

## Calculate some averages, min, and max
historical_fact_reg_sample_covariate_check %>%
  group_by(trait, selection) %>%
  summarize(prop_full_cov_recovered_mean = mean(prop_full_cov_recovered),
            prop_full_cov_recovered_min = min(prop_full_cov_recovered),
            prop_full_cov_recovered_max = max(prop_full_cov_recovered))

## Calculate the frequency of recovery for full-data covariates
historical_fact_reg_sample_covariate_prop <- historical_fact_reg_sample_covariate_check %>%
  group_by(trait, selection) %>%
  do({
    df <- .
    ## Calculate contingency table for recovered covariates
    prop_full_cov <- table(unlist(df$full_cov_recovered)) %>% 
      as.data.frame() %>% 
      arrange(desc(Freq)) %>%
      mutate(Freq = Freq / nrow(df))
    
    ## Calculate contingency table for new covariates
    prop_new_cov <- table(unlist(df$new_covariates)) %>% 
      as.data.frame() %>% 
      arrange(desc(Freq)) %>%
      mutate(Freq = Freq / nrow(df))
    
    # Return tibble
    tibble(prop_full_cov = list(prop_full_cov), prop_new_cov = list(prop_new_cov))
    
  }) %>% ungroup()


## Compare the frequency of recovering a full-data covariate with the proportion
## of variance that it explains
historical_fact_reg_sample_covariate_prop %>%
  unnest(prop_full_cov) %>%
  rename(term = Var1) %>%
  inner_join(., unnest(subset(fr_var_summary, timeframe == "concurrent"))) %>%
  ggplot(aes(x = Freq, y = prop_var_exp)) +
  geom_point() +
  facet_grid(selection ~ trait)










### Test different timeframes of historical covariates ###

# Only use historial models
historical_fr_results <- ec_fr_results_df %>%
  filter(timeframe == "historical") %>%
  unnest(calls) %>%
  filter(model %in% c("model4", "model5"))


## Test the historical model above with the different timeframe
historical_ec_tomodel_timeframe_fact_reg_test <- historical_ec_tomodel_timeframe_centered %>%
  map(~NULL)

# Iterate over this list
for (i in seq_along(historical_ec_tomodel_timeframe_fact_reg_test)) {
  
  # Subset the covariate data
  covariate_data <- historical_ec_tomodel_timeframe_centered[[i]]
  
  # Nest the location blues with covariate data
  model_data <- left_join(S2_MET_loc_BLUEs_tomodel, covariate_data, by = "location") %>%
    group_by(trait) %>%
    nest()
  
  # Map over the model calls and refit models
  historical_fr_refit_i <- historical_fr_results %>%
    # Add data
    left_join(., model_data, by = "trait") %>%
    # Fit models
    mutate(fit = map2(call, data, ~lm(formula = as.formula(str_remove_all(.x, "lm\\(formula = |\\, data = data\\)")), data = .y)))
  
  # Return R2 and RMSE
  historical_ec_tomodel_timeframe_fact_reg_test[[i]] <- historical_fr_refit_i %>%
    mutate(rsquare = map_dbl(fit, ~summary(.)$r.squared),
           adj_rsquare = map_dbl(fit, ~summary(.)$adj.r.squared),
           AIC = map_dbl(fit, AIC),
           rmse = map_dbl(fit, ~rmse(model = ., data = model.frame(.)))) %>%
    select(trait, selection, model, rsquare:rmse)
  
}



## Create a tibble to plot
historical_timeframe_fact_reg_test_df <- historical_ec_tomodel_timeframe_fact_reg_test %>%
  tibble(time_frame = names(.), results = .) %>%
  unnest(results) %>%
  # Parse the timeframe
  mutate(time_frame_length = str_remove(time_frame, "time_frame")) %>%
  separate(time_frame_length, c("length", "start_year", "end_year"), sep = "_") %>%
  mutate_at(vars(length, start_year, end_year), parse_guess)


## Plot - only model 5
historical_timeframe_fact_reg_test_df %>%
  # Set negative adj_rsquare to 0
  mutate(adj_rsquare = ifelse(adj_rsquare < 0, 0, adj_rsquare)) %>%
  filter(model == "model5") %>%
  ggplot(aes(x = length, y = adj_rsquare, color = selection)) +
  # geom_point() +
  geom_line() +
  # facet_wrap(~ trait, scales = "free_y") +
  facet_grid( ~ trait) +
  theme_presentation2(8)




### Test different windows of historical covariates ###


## Test the historical model above with the different timeframe
historical_ec_tomodel_window_fact_reg_test <- historical_ec_tomodel_window_centered %>%
  map(~NULL)

# Iterate over this list
for (i in seq_along(historical_ec_tomodel_window_fact_reg_test)) {
  
  # Subset the covariate data
  covariate_data <- historical_ec_tomodel_window_centered[[i]]
  
  # Nest the location blues with covariate data
  model_data <- left_join(S2_MET_loc_BLUEs_tomodel, covariate_data, by = "location") %>%
    group_by(trait) %>%
    nest()
  
  # Map over the model calls and refit models
  historical_fr_refit_i <- historical_fr_results %>%
    # Add data
    left_join(., model_data, by = "trait") %>%
    # Fit models
    mutate(fit = map2(call, data, ~lm(formula = as.formula(str_remove_all(.x, "lm\\(formula = |\\, data = data\\)")), data = .y)))
  
  # Return R2 and RMSE
  historical_ec_tomodel_window_fact_reg_test[[i]] <- historical_fr_refit_i %>%
    mutate(rsquare = map_dbl(fit, ~summary(.)$r.squared),
           adj_rsquare = map_dbl(fit, ~summary(.)$adj.r.squared),
           AIC = map_dbl(fit, AIC),
           rmse = map_dbl(fit, ~rmse(model = ., data = model.frame(.)))) %>%
    select(trait, selection, model, rsquare:rmse)
  
}


## Create a tibble to plot
historical_window_fact_reg_test_df <- historical_ec_tomodel_window_fact_reg_test %>%
  tibble(time_frame = names(.), results = .) %>%
  unnest(results) %>%
  # Parse the timeframe
  mutate(time_frame_length = str_remove(time_frame, "window")) %>%
  separate(time_frame_length, c("size", "start_year", "end_year"), sep = "_") %>%
  mutate_at(vars(size, start_year, end_year), parse_guess)


## Plot - only model 5
historical_window_fact_reg_test_df %>%
  # Set negative adj_rsquare to 0
  mutate(adj_rsquare = ifelse(adj_rsquare < 0, 0, adj_rsquare)) %>%
  filter(model == "model5") %>%
  ggplot(aes(x = end_year, y = adj_rsquare, color = selection)) +
  # geom_point() +
  geom_line() +
  # facet_wrap(~ trait, scales = "free_y") +
  facet_grid(size ~ trait) +
  theme_presentation2(8)



historical_timeframe_fact_reg_test_df %>%
  group_by(trait, selection, model) %>%
  top_n(x = ., n = 1, wt = adj_rsquare) %>%
  ungroup() %>%
  arrange(trait, selection, model) %>%
  select(trait, selection, model, time_frame, adj_rsquare) %>%
  as.data.frame()

## Find the timeline or window with the highest R square
bind_rows(historical_timeframe_fact_reg_test_df, historical_window_fact_reg_test_df) %>%
  group_by(trait, selection, model) %>%
  top_n(x = ., n = 1, wt = adj_rsquare) %>%
  ungroup() %>%
  arrange(trait, selection, model) %>%
  select(trait, selection, model, time_frame, adj_rsquare) %>%
  as.data.frame()


















