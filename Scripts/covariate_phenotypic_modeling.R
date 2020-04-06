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

## Load additional packages
library(modelr)
library(broom)
library(lme4)
library(car)



## significance level
alpha <- 0.05


####################
## Load data
####################

# Load the AMMI output
load(file.path(result_dir, "ammi_model_fit.RData"))

# Load covariates for environments and historical covariates
load(file.path(result_dir, "ec_model_building.RData"))
load(file.path(result_dir, "historical_ec_model_building.RData"))

# Unnest the fitted mixed model
s2_met_blups_tomodel <- mixed_model_fit %>%
  unnest(y_hat) %>%
  select(trait, line_name, environment, value = y_hat) %>%
  left_join(., distinct(S2_MET_BLUEs, environment, location, year))


####################
## Model GxE
####################

# Group by trait and model
ec_phenotypic_modeling <- s2_met_blups_tomodel %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))

# Iterate over rows
for (i in seq_len(nrow(ec_phenotypic_modeling))) {
  
  tr <- ec_phenotypic_modeling$trait[i]
  dat <- ec_phenotypic_modeling$data[[i]]
  # Add covariates
  covariates_by_term <- environmental_relmat_df %>% 
    filter(trait == tr) %>%
    mutate_at(vars(contains("covariate")), ~map(., colnames)) %>%
    select(contains("covariate")) %>%
    as.list() %>%
    map(1)
  
  # Add covariates to the dat
  dat1 <- left_join(dat, ec_tomodel_centered)
  
  ## Base model
  base_fit <- lm(value ~ 1 + line_name + environment, data = dat1)
  base_fit2 <- lm(value ~ 1 + line_name, data = dat1)
  fit1_anova <- tidy(Anova(base_fit, type = "II"))
  
  ## Model with covariates for environment
  form2 <- add_predictors(formula(base_fit2), reformulate(termlabels = covariates_by_term$main_covariate_mat))
  # form2 <- add_predictors(value ~ 1 + line_name, 
  #                         reformulate(termlabels = c(rev(covariates_by_term$main_covariate_mat), "environment")))
  
  env_ec_fit <- update(base_fit, formula = form2)
  fit2_anova <- tidy(Anova(env_ec_fit, type = "II"))
  
  ## Forward stepwise to identify covariates for environmental mean
  fit2_step <- step(object = base_fit2, direction = "forward",
                    scope = add_predictors(formula(base_fit2), reformulate(names(ec_tomodel_centered)[-1])))
  fit2_step_anova <- tidy(Anova(fit2_step, type = "II"))
  
  
  
  # Model with covariates for GxE
  form3 <- add_predictors(form2, reformulate(termlabels = paste0("line_name:", covariates_by_term$interaction_covariate_mat)))
  gxe_ec_fit <- update(base_fit, formula = form3)
  fit3_anova <- tidy(Anova(gxe_ec_fit, type = "II"))
  
  ## Forward stepwise to identify covariates for GxE
  step_out <- capture.output( {fit3_step <- step(object = fit2_step, direction = "forward",
                    scope = add_predictors(formula(formula(fit2_step)), 
                                           reformulate(paste0("line_name:", names(ec_tomodel_centered)[-1])))) })
  # If no DF, drop the last term
  if (df.residual(fit3_step) == 0) {
    # Drop the last two terms
    
    # new formula
    form3_new <- terms(formula(fit3_step)) %>%
      drop.terms(termobj = ., dropx = tail(seq_along(attr(., "term.labels")), 2), keep.response = TRUE) %>%
      formula()
    # refit
    fit3_step1 <- update(fit3_step, formula = form3_new)
    fit3_step <- fit3_step1
    
  }
  
  fit3_step_anova <- tidy(Anova(fit3_step, type = "II"))
  
  ## Get anova tables for each model
  # Find the proportion of variance explained by covariates and leftover
  main_cov_prop <- subset(fit3_anova, term %in% covariates_by_term$main_covariate_mat) %>%
    add_row(term = "environment", sumsq = subset(fit1_anova, term == "environment", sumsq, drop = TRUE)) %>%
    mutate(propSS = sumsq / last(sumsq))
  int_cov_prop <- subset(fit3_anova, term %in% paste0("line_name:", covariates_by_term$interaction_covariate_mat)) %>%
    add_row(term = "line_name:environment", sumsq = subset(fit1_anova, term == "Residuals", sumsq, drop = TRUE)) %>%
    mutate(propSS = sumsq / last(sumsq))
  # Combine
  var_prop <- bind_rows(
    mutate(main_cov_prop, term1 = "environment"), 
    mutate(int_cov_prop, term1 = "line_name:environment")
  )
  
  ## Now for stepwise
  main_cov_prop_step <- subset(fit2_step_anova, ! term %in% c("line_name", "Residuals")) %>%
    add_row(term = "environment", sumsq = subset(fit1_anova, term == "environment", sumsq, drop = TRUE)) %>%
    mutate(propSS = sumsq / last(sumsq))
  int_cov_prop_step <- subset(fit3_step_anova, str_detect(term, "line_name:")) %>%
    add_row(term = "line_name:environment", sumsq = subset(fit1_anova, term == "Residuals", sumsq, drop = TRUE)) %>%
    mutate(propSS = sumsq / last(sumsq))
  # Combine
  var_prop_step <- bind_rows(
    mutate(main_cov_prop_step, term1 = "environment"), 
    mutate(int_cov_prop_step, term1 = "line_name:environment")
  )
  
  # Return things
  ec_phenotypic_modeling$out[[i]] <- list(anova_tables = list(fit1 = fit1_anova, fit2 = fit2_anova, fit3 = fit3_anova),
                                          var_prop = list(main_ec = main_cov_prop, int_ec = int_cov_prop, var_prop = var_prop),
                                          var_prop_step = list(main_ec = main_cov_prop_step, int_ec = int_cov_prop_step, 
                                                               var_prop = var_prop_step))
  
}


## Output varprop tables
ec_var_prop <- ec_phenotypic_modeling %>% 
  mutate(var_prop = map(out, ~.$var_prop$var_prop)) %>% 
  unnest(var_prop)
write_csv(x = ec_var_prop, path = file.path(fig_dir, "environment_ec_varprop.csv"))

ec_var_prop %>% 
  group_by(trait, term1) %>% 
  summarize(cumulative_propSS = sum(propSS)) %>% 
  mutate(cumulative_propSS = cumulative_propSS - 1)







####################
## Model GxL
####################

# Group by trait and model
historical_ec_phenotypic_modeling <- s2_met_blups_tomodel %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))

# Iterate over rows
for (i in seq_len(nrow(historical_ec_phenotypic_modeling))) {
  
  tr <- historical_ec_phenotypic_modeling$trait[i]
  dat <- historical_ec_phenotypic_modeling$data[[i]]
  # Add covariates
  covariates_by_term <- location_relmat_df %>% 
    filter(trait == tr, time_frame == "time_frame5") %>%
    mutate_at(vars(contains("covariate")), ~map(., colnames)) %>%
    select(contains("covariate")) %>%
    as.list() %>%
    map(1)
  
  # Add covariates to the dat
  dat1 <- left_join(dat, select(historical_ec_tomodel_centered$time_frame5, location, unique(unlist(covariates_by_term)))) %>%
    mutate_at(vars(line_name, environment, location, year), ~fct_contr_sum(as.factor(.)))
  
  ## Base model
  base_fit <- lm(value ~ 1 + line_name + location + year + line_name:year, data = dat1)
  fit1_anova <- tidy(Anova(base_fit, type = "II"))
  
  ## Model with covariates for environment
  form2 <- add_predictors(value ~ 1 + line_name + year + line_name:year, 
                          reformulate(termlabels = covariates_by_term$main_covariate_mat))
  # form2 <- add_predictors(value ~ 1 + line_name, 
  #                         reformulate(termlabels = c(rev(covariates_by_term$main_covariate_mat), "environment")))
  
  env_ec_fit <- update(base_fit, formula = form2)
  fit2_anova <- tidy(Anova(env_ec_fit, type = "II"))
  
  
  # Model with covariates for GxE
  form3 <- add_predictors(form2, reformulate(termlabels = paste0("line_name:", covariates_by_term$interaction_covariate_mat)))
  gxe_ec_fit <- update(base_fit, formula = form3)
  fit3_anova <- tidy(Anova(gxe_ec_fit, type = "II"))
  
  ## Get anova tables for each model
  # Find the proportion of variance explained by covariates and leftover
  main_cov_prop <- subset(fit3_anova, term %in% covariates_by_term$main_covariate_mat) %>%
    add_row(term = "location", sumsq = subset(fit1_anova, term == "location", sumsq, drop = TRUE)) %>%
    mutate(propSS = sumsq / last(sumsq))
  int_cov_prop <- subset(fit3_anova, term %in% paste0("line_name:", covariates_by_term$interaction_covariate_mat)) %>%
    add_row(term = "line_name:location", sumsq = subset(fit1_anova, term == "Residuals", sumsq, drop = TRUE)) %>%
    mutate(propSS = sumsq / last(sumsq))
  # Combine
  var_prop <- bind_rows(
    mutate(main_cov_prop, term1 = "environment"), 
    mutate(int_cov_prop, term1 = "line_name:environment")
  )
  
  # Return things
  historical_ec_phenotypic_modeling$out[[i]] <- list(
    anova_tables = list(fit1 = fit1_anova, fit2 = fit2_anova, fit3 = fit3_anova),
    var_prop = list(main_ec = main_cov_prop, int_ec = int_cov_prop, var_prop = var_prop))
  
}


## Output varprop tables
historical_ec_var_prop <- historical_ec_phenotypic_modeling %>% 
  mutate(var_prop = map(out, ~.$var_prop$var_prop)) %>% 
  unnest(var_prop)
write_csv(x = historical_ec_var_prop, path = file.path(fig_dir, "location_ec_varprop.csv"))

ec_var_prop %>% 
  group_by(trait, term1) %>% 
  summarize(cumulative_propSS = sum(propSS)) %>% 
  mutate(cumulative_propSS = cumulative_propSS - 1)























