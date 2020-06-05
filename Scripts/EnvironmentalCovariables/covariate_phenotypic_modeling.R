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


# Plot prediction/observation
tr <- "PlantHeight"
# form <- concurrent_feature_selection %>%
form <- historical_feature_selection %>%
  filter(trait == tr, feat_sel_type == "rfa_cv", model == "model3") %>%
  pull(adhoc) %>%
  map(1) %>%
  unlist() %>%
  reformulate(termlabels = ., response = "value")
# dat <- subset(S2_MET_BLUEs_tomodel, trait == tr) %>%
dat <- subset(S2_MET_loc_BLUEs_tomodel, trait == tr) %>%
  droplevels() %>%
  # left_join(., ec_tomodel_centered, by = "environment") %>%
  left_join(., historical_ec_tomodel_centered_use, by = "location") %>%
  mutate_at(vars(matches("line_name|environment|location")), ~fct_contr_sum(as.factor(.)))
loo_indices <- dat %>%
  group_by_at(vars(matches("environment|location"))) %>%
  crossv_loo_grouped() %>%
  pull(train) %>%
  map("idx")
train_list <- map(loo_indices, ~lm(form, dat[.x,]))
test_pred <- map2_df(loo_indices, train_list, ~add_predictions(data = dat[-.x,], model = .y))
plot(value ~ pred, test_pred); with(test_pred, cor(value, pred))




## Combine everything
ec_fr_results <- bind_rows(
  mutate(concurrent_fact_reg, timeframe = "concurrent"),
  mutate(historical_fact_reg, timeframe = "historical")
) %>% gather(selection, fit, -trait, -model, -timeframe)
  

## Use the base model to determine the proportion of E and GE variance explained by
## covariates
fr_var_summary <- ec_fr_results %>%
  group_by(trait, timeframe, selection) %>%
  do(var_prop_summary = {
    df <- .
    
    # Get the base model ANOVA
    base_anova <- tidy(anova(subset(df, model == "base_alt", fit, drop = TRUE)[[1]])) %>%
      # Remove line name term
      filter(term != "line_name") %>%
      # Rename residuals to GxE or GxL
      mutate(term = ifelse(term == "Residuals", paste0("line_name:", term[1]), term),
             group = ifelse(str_detect(term, ":"), "interaction", "main")) %>%
      select(term, group, sumsq)
      
    # Get model 3 or 5
    full_anova <- tidy(Anova(subset(df, model %in% c("model3", "model5"), fit, drop = TRUE)[[1]], type = "II"))
    
    # Use the full model to determine the proportion of variance explained by main 
    # effect or interaction covariates
    full_anova %>%
      # Assign main effect or interaction term
      mutate(group = case_when(
        term %in% c("line_name", "Residuals") ~ term,
        str_detect(term, ":") ~ "interaction",
        TRUE ~ "main")) %>%
      # Add group
      left_join(., base_anova, by = "group") %>%
      # Calculate proportion of variance explained
      mutate(prop_var_exp = sumsq.x / sumsq.y) %>%
      select(term = term.x, group, sumsq = sumsq.x, group_sumsq = sumsq.y, prop_var_exp)
    
  }) %>% ungroup()

## Number of covariates per group
fr_var_summary %>%
  unnest() %>%
  filter(group %in% c("main", "interaction")) %>%
  group_by(trait, timeframe, selection, group) %>%
  summarize(nEC = n_distinct(term)) %>%
  spread(group, nEC) %>%
  as.data.frame()


# trait  timeframe    selection interaction main
# 1  GrainProtein concurrent        adhoc           2    6
# 2  GrainProtein concurrent adhoc_nosoil           2    8
# 3  GrainProtein concurrent      apriori           5    5
# 4  GrainProtein historical        adhoc           1    5
# 5  GrainProtein historical adhoc_nosoil          NA    5
# 6  GrainProtein historical      apriori           5    5
# 7    GrainYield concurrent        adhoc           2   23
# 8    GrainYield concurrent adhoc_nosoil          NA   20
# 9    GrainYield concurrent      apriori           5    5
# 10   GrainYield historical        adhoc          NA   11
# 11   GrainYield historical adhoc_nosoil          NA   12
# 12   GrainYield historical      apriori           5    5
# 13  HeadingDate concurrent        adhoc          NA   18
# 14  HeadingDate concurrent adhoc_nosoil          NA    9
# 15  HeadingDate concurrent      apriori           2    2
# 16  HeadingDate historical        adhoc          NA    9
# 17  HeadingDate historical adhoc_nosoil          NA   10
# 18  HeadingDate historical      apriori           2    2
# 19  PlantHeight concurrent        adhoc           2   23
# 20  PlantHeight concurrent adhoc_nosoil          NA   14
# 21  PlantHeight concurrent      apriori           3    3
# 22  PlantHeight historical        adhoc           2   11
# 23  PlantHeight historical adhoc_nosoil           1   11
# 24  PlantHeight historical      apriori           3    3
# 25   TestWeight concurrent        adhoc           9   10
# 26   TestWeight concurrent adhoc_nosoil           7   10
# 27   TestWeight concurrent      apriori           5    5
# 28   TestWeight historical        adhoc           5    5
# 29   TestWeight historical adhoc_nosoil           2    6
# 30   TestWeight historical      apriori           5    5







## Print a summary for each trait/time_frame/selection
fr_var_summary_print <- fr_var_summary %>% 
  unnest() %>% 
  mutate(total_SS = ifelse(is.na(group_sumsq), sumsq, group_sumsq)) %>% 
  group_by(trait, timeframe, selection, group) %>% 
  summarize(total_SS = mean(total_SS), prop_SS = sum(prop_var_exp)) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("line_name", "main", "interaction", "Residuals"))) %>%
  arrange(trait, timeframe, selection, group)

# Save as table
write_csv(x = fr_var_summary_print, path = file.path(fig_dir, "covariate_fact_reg_variance_explained.csv"))


# Just interaction and main
fr_var_summary_print %>%
  filter(group %in% c("main", "interaction")) %>%
  select(-total_SS) %>%
  spread(selection, prop_SS) %>%
  as.data.frame()

## Proportion of total main effect or interaction variance explained

# trait        timeframe  group         adhoc apriori
# 1 GrainProtein concurrent main         0.0264   1.23 
# 2 GrainProtein concurrent interaction  0.430    0.680
# 3 GrainProtein historical main         0.734    1.99 
# 4 GrainProtein historical interaction  0.364    0.912
# 5 GrainYield   concurrent main         0.323    0.208
# 6 GrainYield   concurrent interaction  0.151    0.234
# 7 GrainYield   historical main         0.271    0.719
# 8 GrainYield   historical interaction NA        0.512
# 9 HeadingDate  concurrent main         0.127    0.545
# 10 HeadingDate  concurrent interaction NA        0.105
# 11 HeadingDate  historical main         0.0293   0.933
# 12 HeadingDate  historical interaction NA        0.350
# 13 PlantHeight  concurrent main         0.491    0.142
# 14 PlantHeight  concurrent interaction  0.266    0.133
# 15 PlantHeight  historical main        NA        0.173
# 16 PlantHeight  historical interaction  0.258    0.221
# 17 TestWeight   concurrent main        NA        0.343
# 18 TestWeight   concurrent interaction  0.664    0.482
# 19 TestWeight   historical main        NA        0.527
# 20 TestWeight   historical interaction  0.814    1.21



## Copy formulas and a single model frame
ec_fr_results_df <- ec_fr_results %>%
  group_by(trait, timeframe, selection) %>%
  do({
    df <- .
    
    # Create a tibble of models and calls
    formula_df <- select(df, model, fit) %>%
      mutate(call = map(fit, "call"),
             call = map(call, ~deparse(.) %>% str_trim(.) %>% paste0(., collapse = " "))) %>%
      select(-fit)
    
    # Get the model.frame for the largest model
    mf <- model.frame(subset(df, model %in% c("model3", "model5"), fit, drop = TRUE)[[1]]) %>%
      # Strip attributes to save disk space
      merTools:::stripAttributes()
      
    
    # Return a tibble
    tibble(calls = list(formula_df), model_frame = list(mf))
    
  }) %>% ungroup()






## Save results
save("ec_fr_results_df", "fr_var_summary", file = file.path(result_dir, "factorial_regression_results_nasapower.RData"))
# save("ec_fr_results_df", "fr_var_summary", file = file.path(result_dir, "factorial_regression_results_daymet.RData"))




## Double-check models for overfitting

# Load the results
load(file.path(result_dir, "factorial_regression_results_keep.RData"))


## Refit models
ec_fr_refit <- ec_fr_results_df %>%
  mutate(calls = map(calls, "call") %>% map(last)) %>%
  mutate(fit = list(NULL),
         overparam_main = NA,
         overparam_int = NA,
         R2 = as.numeric(NA),
         R2_adj = as.numeric(NA),
         LL = as.numeric(NA))

# Iterate over rows
for (i in seq_len(nrow(ec_fr_refit))) {
  row <- ec_fr_refit[i, ]
  
  data <- row$model_frame[[1]]
  fit <- eval(parse(text = row$calls[[1]]))
  
  ## Compare number of covariates versus number of environmnents/locations
  covariates <- formula(fit) %>%
    terms() %>% 
    attr(., "term.labels") %>% 
    subset(., . != "line_name") %>%
    split(., ifelse(str_detect(., ":"), "interaction", "main"))
  
  nCov <- map_dbl(covariates, length)
  
  ## The number of total df for covariates is 2 * (number of E or L - 1)
  J <- max(map_dbl(data[covariates$main], n_distinct)) - 1
  H <- 2 * J
  
  ec_fr_refit$fit[[i]] <- fit
  ec_fr_refit$overparam_main[i] <- nCov["main"] > H
  ec_fr_refit$overparam_int[i] <- nCov["interaction"] > H
  ec_fr_refit$R2[i] <- summary(fit)$r.squared
  ec_fr_refit$R2_adj[i] <- summary(fit)$adj.r.squared
  ec_fr_refit$LL[i] <- as.numeric(logLik(fit))
  
}


## Cross-validation
# Demonstrate with grain yield
# Leave one location out
fit_example <- subset(ec_fr_refit, trait == "GrainYield" & timeframe == "concurrent" & selection == "adhoc", fit, drop = TRUE)[[1]]


## Get the data, assign locations
mf <- model.frame(fit_example)
cov_use <- which.max(sapply(mf[,-1:-2], n_distinct))
mf1 <- mf %>%
  as_tibble() %>%
  mutate_at(vars(2 + cov_use), list(location = ~paste0("loc", as.numeric(as.factor(.)))))

## Leave-one-location-out
crossv_data <- mf1 %>%
  group_by(location) %>%
  crossv_loo_grouped()


## The original
crossv_fit <- map(crossv_data$train, ~update(fit_example, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

# RMSE ~ 500


## Testing ##
# Remove some soil variables
# 1. base saturation
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% "topsoil_base_saturation"), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

# RMSE ~ 1000


# 2. subsoil_teb
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% "subsoil_teb"), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

# RMSE ~ 730


# 3. subsoil_calcium_carbonate
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% "subsoil_calcium_carbonate"), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

# RMSE ~ 1016

# 4. pH
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which( str_detect(string = attr(., "term.labels"), pattern = "pH")), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

# 4. pH
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which( str_detect(string = attr(., "term.labels"), pattern = "topsoil_pH_h2o")), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))



## Weather
# 4. mean temperatures
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which( str_detect(string = attr(., "term.labels"), pattern = "tmean")), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))


# 5. Choose min, max, or mean temperature for a growth stage
# Choose the min, max, or mean with highest SS
to_remove <- c("late_vegetative.mint_mean", "late_vegetative.mint_mean", "grain_fill.maxt_mean")
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% to_remove), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))


# Test a reduced model using the tests from above
## Choose the min, max, or mean with highest SS
to_remove <- c("early_vegetative.maxt_mean", "late_vegetative.maxt_mean", "subsoil_teb")
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% to_remove), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))



## Effect plots
all_effects <- allEffects(fit_test1)
plot(all_effects[str_detect(names(all_effects), "line_name", negate = TRUE)])




## Try another example

## Cross-validation
# Demonstrate with grain yield
# Leave one location out
fit_example <- subset(ec_fr_refit, trait == "HeadingDate" & timeframe == "concurrent" & selection == "adhoc_nosoil", fit, drop = TRUE)[[1]]


## Get the data, assign locations
mf <- model.frame(fit_example)
cov_use <- which.max(sapply(mf[,-1:-2], n_distinct))
mf1 <- mf %>%
  as_tibble() %>%
  mutate_at(vars(2 + cov_use), list(location = ~paste0("loc", as.numeric(as.factor(.)))))

## Leave-one-location-out
crossv_data <- mf1 %>%
  group_by(location) %>%
  crossv_loo_grouped()


## The original
crossv_fit <- map(crossv_data$train, ~update(fit_example, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

# RMSE ~ 500


## Testing ##
## Choose the min, max, or mean with highest SS
to_remove <- c("early_vegetative.maxt_mean", "late_vegetative.tmean_mean")
fit_test1 <- terms(fit_example) %>% 
  drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% to_remove), keep.response = TRUE) %>% 
  formula() %>%
  update(object = fit_example, formula = ., data = mf1)

## Fit models
crossv_fit <- map(crossv_data$train, ~update(fit_test1, data = .))
crossv_predict <- map2_df(crossv_fit, crossv_data$test, ~add_predictions(as.data.frame(.y), .x))
# Plot
plot(value ~ pred, crossv_predict); cor(crossv_predict$value, crossv_predict$pred)
mean(map2_dbl(crossv_fit, crossv_data$test, rmse))

## Effect plots
all_effects <- allEffects(fit_test1)
plot(all_effects[str_detect(names(all_effects), "line_name", negate = TRUE)])










## Run for all traits
ec_fr_refit_crossv <- ec_fr_refit %>%
  mutate(full_rmse = as.numeric(NA),
         crossv_rmse = list(NULL),
         predictions = list(NULL))

for (i in seq_len(nrow(ec_fr_refit_crossv))) {
  
  row <- ec_fr_refit_crossv[i, ]
  fit_example <- row$fit[[1]]
  
  ## Get the data, assign locations
  mf <- model.frame(fit_example)
  cov_use <- which.max(sapply(mf[,-1:-2], n_distinct))
  mf1 <- mf %>%
    as_tibble() %>%
    mutate_at(vars(2 + cov_use), list(location = ~paste0("loc", as.numeric(as.factor(.)))))
  
  ## Leave-one-location-out
  crossv_data <- mf1 %>%
    group_by(location) %>%
    crossv_loo_grouped()
  
  ## Fit models
  crossv_fit <- map(crossv_data$train, ~update(fit_example, data = .))
  
  # Predict, measure RMSE
  ec_fr_refit_crossv$crossv_rmse[[i]] <- map2_dbl(crossv_fit, crossv_data$test, rmse)
  ec_fr_refit_crossv$full_rmse[i] <- rmse(fit_example, mf1)
  ec_fr_refit_crossv$predictions[[i]] <- crossv_data %>% 
    mutate(prediction = map2(test, crossv_fit, ~add_predictions(as.data.frame(.x), .y))) %>% 
    unnest(prediction)
  
  
}


## Summarize rmse from crossv
ec_fr_refit_crossv %>%
  mutate(crossv_rmse = map(crossv_rmse, unlist)) %>%
  mutate(crossv_rmse_mean = map_dbl(crossv_rmse, mean)) %>%
  select(trait, timeframe, selection, full_rmse, crossv_rmse_mean) %>%
  as.data.frame()

# trait  timeframe    selection    full_rmse crossv_rmse_mean
# 1  GrainProtein concurrent        adhoc    0.5964271         2.730672
# 2  GrainProtein concurrent adhoc_nosoil    0.5957886         1.630300
# 3  GrainProtein concurrent      apriori    0.5130104         1.543583
# 4  GrainProtein historical        adhoc    0.5611728         2.664219
# 5  GrainProtein historical adhoc_nosoil    0.7035249         1.791760
# 6  GrainProtein historical      apriori    0.3603326        22.134577
# 7    GrainYield concurrent        adhoc  506.1560506       572.140996
# 8    GrainYield concurrent adhoc_nosoil  859.4216257      3036.594374
# 9    GrainYield concurrent      apriori 1605.4889060      1744.402187
# 10   GrainYield historical        adhoc  410.7683455       435.729908
# 11   GrainYield historical adhoc_nosoil  410.7682400      3021.077893
# 12   GrainYield historical      apriori 1313.1826189      2104.068339
# 13  HeadingDate concurrent        adhoc    3.1507774         9.791189
# 14  HeadingDate concurrent adhoc_nosoil    3.4593926         4.785399
# 15  HeadingDate concurrent      apriori    4.8175084         4.778512
# 16  HeadingDate historical        adhoc    1.8426994         2.676330
# 17  HeadingDate historical adhoc_nosoil    1.8424856       465.127939
# 18  HeadingDate historical      apriori    3.8517715         4.689140
# 19  PlantHeight concurrent        adhoc    4.3562909        27.700869
# 20  PlantHeight concurrent adhoc_nosoil    8.4167638        14.777691
# 21  PlantHeight concurrent      apriori   12.1293827        12.186621
# 22  PlantHeight historical        adhoc    2.6918425       214.137583
# 23  PlantHeight historical adhoc_nosoil    2.9277105       117.649132
# 24  PlantHeight historical      apriori    9.2823767        13.235534
# 25   TestWeight concurrent        adhoc    6.5344199       339.627867
# 26   TestWeight concurrent adhoc_nosoil   10.0019010       221.750772
# 27   TestWeight concurrent      apriori   36.9320143        63.332485
# 28   TestWeight historical        adhoc    6.0143878      5296.860374
# 29   TestWeight historical adhoc_nosoil   14.6656554        66.636482
# 30   TestWeight historical      apriori   29.8547520     30152.414863





## Analyze results of factorial regression sampling
## 
## In this procedure, one environment or location was dropped, and
## the same covariate selection procedure was used.
## 

# Load the results
load(file.path(result_dir, "factorial_regression_results_sample.RData"))


## Analyze concurrent covariates ##

# Subset the data
concurrent_fact_reg_sample_rmse <- concurrent_fact_reg_sample %>%
  filter(model == "model3") %>% # Skip the base model
  select(-test_predictions, -covariates)
  
# Plot RMSE
concurrent_fact_reg_sample_rmse %>%
  filter(!(trait == "GrainProtein" & rmse > 300),
         !(trait == "TestWeight" & rmse > 2000)) %>%
  ggplot(aes(x = selection, y = rmse, color = selection)) +
  geom_boxplot() +
  facet_wrap(~ trait, scales = "free_y") +
  theme_genetics(base_size = 10)

# Plot accuracy
concurrent_fact_reg_sample_rmse %>%
  # filter(!(trait == "GrainProtein" & rmse > 300),
  #        !(trait == "TestWeight" & rmse > 2000)) %>%
  ggplot(aes(x = selection, y = acc, color = selection)) +
  geom_boxplot() +
  facet_wrap(~ trait, scales = "free_y") +
  theme_genetics(base_size = 10)

## Summarize both
concurrent_fact_reg_sample_summary <- concurrent_fact_reg_sample_rmse %>% 
  gather(stat, value, rmse, acc) %>% 
  group_by(trait, model, selection, stat) %>% 
  do({
    x <- .$value
    stats <- boxplot.stats(x)
    tibble(median = median(x), lower = stats$conf[1], upper = stats$conf[2])
  }) %>% ungroup()

concurrent_fact_reg_sample_summary %>%
  arrange(stat) %>%
  as.data.frame()
   


## Compare covariates selected in each model with those selected using the 
## full dataset

# Prepare the full data results
full_data_fact_reg <- fr_var_summary %>%
  filter(timeframe == "concurrent", selection != "apriori") %>%
  mutate(full_covariates = map(var_prop_summary, ~subset(., ! term %in% c("line_name", "Residuals"), term, drop = TRUE))) %>%
  select(-var_prop_summary)


# How often do we get the same exact model?
# How often do we get the same full-data covariate
# How often do we get a new covariate?
concurrent_fact_reg_sample_covariate_check <- concurrent_fact_reg_sample %>%
  filter(model == "model3") %>%
  select(trait, dropped_group, selection, covariates) %>%
  inner_join(., full_data_fact_reg) %>%
  # Calculate the proportion of full-data covariates that were recovered
  mutate(full_cov_recovered = map2(full_covariates, covariates, intersect),
         prop_full_cov_recovered = map2_dbl(full_covariates, covariates, ~mean(.x %in% .y)),
         new_covariates = map2(covariates, full_covariates, setdiff))

## Calculate some averages, min, and max
concurrent_fact_reg_sample_covariate_check %>%
  group_by(trait, selection) %>%
  summarize(prop_full_cov_recovered_mean = mean(prop_full_cov_recovered),
            prop_full_cov_recovered_min = min(prop_full_cov_recovered),
            prop_full_cov_recovered_max = max(prop_full_cov_recovered))

## The proportion of times when the same covariates are recovered seems to be correlated
## with the number of environments available for a trait

## Calculate the frequency of recovery for full-data covariates
concurrent_fact_reg_sample_covariate_prop <- concurrent_fact_reg_sample_covariate_check %>%
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
concurrent_fact_reg_sample_covariate_prop %>%
  unnest(prop_full_cov) %>%
  rename(term = Var1) %>%
  inner_join(., unnest(subset(fr_var_summary, timeframe == "concurrent"))) %>%
  ggplot(aes(x = Freq, y = prop_var_exp)) +
  geom_point() +
  facet_grid(selection ~ trait)





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


















