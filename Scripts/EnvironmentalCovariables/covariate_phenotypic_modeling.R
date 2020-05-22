## S2MET phenotypic modeling with covariates
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## Model the impact of covariates on GxE and LxE
## 

# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))


## Add new packages to load
pkgs <- union(pkgs, c("modelr", "broom", "lme4", "car"))
# Load these packages
invisible(lapply(X = pkgs, library, character.only = TRUE))

## significance level
alpha <- 0.05



# Load data ---------------------------------------------------------------

# Load covariates for environments and historical covariates
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


## Filter BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
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


## A priori covariates - quite minimal; each should have a citation 
apriori_covariate_df <- trait_covariate_df %>%
  filter(covariate %in% c("early_vegetative.gdd_sum", "late_vegetative.gdd_sum", "flowering.mint_mean",
                          "grain_fill.maxt_mean", "grain_fill.water_balance_sum"))




# # Load packages for local job
# invisible(lapply(X = pkgs, library, character.only = TRUE))
 

## Group by trait and model
concurrent_fact_reg <- S2_MET_BLUEs_tomodel %>%
  group_by(trait) %>%
  do({
    
    df <- .
    # Factorize
    df1 <- df %>%
      droplevels() %>%
      left_join(., ec_tomodel_centered, by = "environment") %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    
    # Apriori
    ## Add covariates - filter
    covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "apriori")
    
    # Ad hoc
    ## Add covariates - filter
    covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
    adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step")
    
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



## Load packages for running local job
invisible(lapply(X = pkgs, library, character.only = TRUE))

## Drop one environment at a time and see if the same covariates
## are identified


## Group by trait and nest
concurrent_fact_reg_sample <- S2_MET_BLUEs_tomodel %>%
  # Split by trait
  split(.$trait) %>%
  # LOO based on environment grouping
  imap_dfr(~group_by(.x, environment) %>% crossv_loo_grouped(.) %>% mutate(trait = .y)) %>%
  rename(dropped_group = environment) %>%
  mutate(out = list(NULL))

# The first null element
first_null <- min(which(sapply(concurrent_fact_reg_sample$out, is.null)))

for (i in seq(first_null, nrow(concurrent_fact_reg_sample))) {
    
  row <- concurrent_fact_reg_sample[i,]
  df <- as_tibble(row$train[[1]])
  
  # Factorize
  df1 <- df %>%
    droplevels() %>%
    left_join(., ec_tomodel_centered, by = "environment") %>%
    mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
  
  test <- as_tibble(row$test[[1]]) %>%
    left_join(., ec_tomodel_centered, by = "environment")
  
  # Apriori
  ## Add covariates - filter
  covariates_use <- subset(apriori_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
  apriori_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "apriori")
  
  # Ad hoc
  ## Add covariates - filter
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
  adhoc_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step")
  
  # Ad hoc - without soil
  ## Add covariates - filter
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T), 
                           covariate, drop = TRUE)
  adhoc_nosoil_out <- fact_reg(data = df1, covariates = covariates_use, env = "environment", method = "step")
  
  ## df of output
  results <- tibble(model = c("base", "base_alt", "model2", "model3"),
                    apriori = apriori_out,
                    adhoc = adhoc_out,
                    adhoc_nosoil = adhoc_nosoil_out)
  
  
  ## Predict the test set using each model
  ## return these results
  concurrent_fact_reg_sample$out[[i]] <- results %>%
    gather(selection, fit, -model) %>%
    filter(model != "base_alt") %>%
    mutate(test_predictions = map(fit, ~add_predictions(data = test, model = .) %>% select(line_name, value, pred)),
           covariates = map(fit, ~str_subset(string = attr(terms(formula(.)), "term.labels"), pattern = "line_name$", negate = TRUE)),
           rmse = map_dbl(fit, ~rmse(model = ., data = test)),
           acc = map_dbl(test_predictions, ~cor(.$value, .$pred))) %>%
    select(-fit)
  
}

  









# Historical covariates =====================================================


## Ideas
## 1. Calculate locations means and use that as the input (and validation)
## 2. Use all data and fit a model with random year, location:year, and g:year
## 3. Use all data and fit the g + location model


## Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp)


# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = names(historical_ec_tomodel_centered$time_frame5)[-1:-2]) %>%
  filter(
    covariate != "awc_range",
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 

## Select the historical covariate data
historical_ec_tomodel_centered_use <- historical_ec_tomodel_centered$time_frame5
  



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

## Group by trait and model
historical_fact_reg <- S2_MET_loc_BLUEs_tomodel %>%
  group_by(trait) %>%
  do({
    
    df <- .
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
    tibble(model = c("base", "base_alt", "model4", "model5"),
           apriori = apriori_out,
           adhoc = adhoc_out,
           adhoc_nosoil = adhoc_nosoil_out)
    
  }) %>% ungroup()





## Load packages for running local job
invisible(lapply(X = pkgs, library, character.only = TRUE))

## Drop one location at a time and see if the same covariates
## are identified
## Group by trait and model
historical_fact_reg_sample <- S2_MET_loc_BLUEs_tomodel %>%
  # Split by trait
  split(.$trait) %>%
  # LOO based on location grouping
  imap_dfr(~group_by(.x, location) %>% crossv_loo_grouped(.) %>% mutate(trait = .y)) %>%
  rename(dropped_group = location) %>%
  group_by(trait, dropped_group) %>%
  do({
    
    row <- .
    df <- as_tibble(row$train[[1]])
    
    
    # Factorize
    df1 <- df %>%
      droplevels() %>%
      left_join(., historical_ec_tomodel_centered_use, by = "location") %>%
      mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
    
    test <- as_tibble(row$test[[1]]) %>%
      left_join(., historical_ec_tomodel_centered_use, by = "location")
    
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
    
    ## df of output
    results <- tibble(model = c("base", "base_alt", "model2", "model3"),
                      apriori = apriori_out,
                      adhoc = adhoc_out,
                      adhoc_nosoil = adhoc_nosoil_out)
    
    
    ## Predict the test set using each model
    results %>%
      gather(selection, fit, -model) %>%
      filter(model != "base_alt") %>%
      mutate(test_predictions = map(fit, ~add_predictions(data = test, model = .) %>% select(line_name, value, pred)),
             rmse = map_dbl(fit, ~rmse(model = ., data = test)),
             acc = map_dbl(test_predictions, ~cor(.$value, .$pred)))
    
  }) %>% ungroup()










# Analyze results ---------------------------------------------------------




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
save("ec_fr_results_df", "fr_var_summary", file = file.path(result_dir, "factorial_regression_results.RData"))




## Double-check models for overfitting

# Load the results
load(file.path(result_dir, "factorial_regression_results.RData"))


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
fit_example <- subset(ec_fr_refit, trait == "GrainYield" & timeframe == "historical" & selection == "adhoc", fit, drop = TRUE)[[1]]

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
crossv_predict <- map2(crossv_fit, crossv_data$test, predict)
# Predict, measure RMSE
crossv_rmse <- map2_dbl(crossv_fit, crossv_data$test, rmse)
rmse(fit_example, mf1)

## Average error in CV is not different that error from full model!!



## Run for all traits
ec_fr_refit_crossv <- ec_fr_refit %>%
  mutate(full_rmse = as.numeric(NA),
         crossv_rmse = list(NULL))

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
