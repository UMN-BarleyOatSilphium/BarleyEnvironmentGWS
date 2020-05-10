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
load(file.path(result_dir, "ec_model_building.RData"))
load(file.path(result_dir, "historical_ec_model_building.RData"))


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
      left_join(., historical_ec_tomodel_centered$time_frame5, by = "location") %>%
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










