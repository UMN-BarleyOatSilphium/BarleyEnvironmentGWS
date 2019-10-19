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
library(lme4)
library(modelr)
library(broom)
library(pbr)
library(cowplot)
library(pls)
library(car)
library(heritability)

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


############################
### Visualization of daily weather data
############################


## Look at daily stats
daily_ec_select <- growth_stage_weather %>%
  inner_join(., env_trials) %>%
  select(environment, dap, growth_stage, mint:water_stress)


## Summarize max temperatures during grain fill
grain_fill_maxt_summary <- daily_ec_select %>% 
  filter(growth_stage == "grain_fill") %>%
  ## Count number of days with 15-18 max temp,
  ## < 15, 18-25 (moderate high), 25 - 30 (high), > 30 (very high)
  mutate(grain_fill_condition = case_when(
    maxt < 15 ~ "suboptimal",
    between(maxt, 15, 18) ~ "optimal",
    between(maxt, 18, 30) ~ "moderate_high",
    between(maxt, 30, 35) ~ "high",
    maxt > 35 ~ "very_high"
  )) %>%
  group_by(environment, grain_fill_condition) %>%
  summarize(nDays = n())

  




############################
### Prepare ECs for modeling
############################


## Create the ECs and select the relevant ones for modeling
ec_select <- growth_stage_covariates %>%
  ## Only use TP environments
  inner_join(., env_trials, by = "trial") %>% 
  select(-trial) %>%
  gather(covariate, value, -environment, -growth_stage) %>%
  unite(covariate, growth_stage, covariate, sep = "_") %>%
  ## Remove trange, relhum, and rain (water stress will cover this)
  filter(str_detect(covariate, "_trange|_rain|_rh2m", negate = TRUE)) %>%
  spread(covariate, value)



## Plot environment and maxt during grain fill
daily_ec_select %>% 
  filter(growth_stage == "grain_fill") %>%
  mutate(environment = factor(environment, levels = ec_select$environment[order(ec_select$grain_fill_maxt, decreasing = TRUE)])) %>%
  ggplot(aes(x = environment, y = maxt)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 35, ymax = 40, fill = "heat stress"), alpha = 0.5) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 35, fill = "high temperature"), alpha = 0.5) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30, fill = "moderate high temperature"), alpha = 0.5) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 15, ymax = 20, fill = "optimal temperature"), alpha = 0.5) +
  geom_violin(fill = alpha("white", 0.5)) +
  # geom_boxplot(fill = NA) +
  scale_y_continuous(breaks = pretty) +
  scale_fill_manual(values = rev(viridis::inferno(direction = -1, n = 10)[1:4]), name = "Temperature stress level",
                    guide = guide_legend(title.position = "top")) +
  labs(subtitle = "Maximum temperature stress during grain fill") +
  theme_presentation2(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(filename = "grain_fill_max_temp_stress.jpg", path = fig_dir, width = 8, height = 5, dpi = 300)





## Prepare the covariates for modeling
## Center, but do not scale. Save the mean for later
ec_tomodel_temp <- ec_select %>%
  mutate_at(vars(-environment), scale, scale = FALSE, center = TRUE)

ec_tomodel_centers <- ec_tomodel_temp %>%
  summarize_at(vars(-environment), ~attr(., "scaled:center")) %>%
  gather(covariate, center)

## Convert scaled to numeric
ec_tomodel <- ec_tomodel_temp %>%
  mutate_at(vars(-environment), as.numeric)


# Vector of covariates
environmental_covariates <- names(ec_tomodel)[-1]

## Scatterplot matrix
scatterplotMatrix(ec_tomodel[,-1], smooth = FALSE)



## Prepare data for modelling:
## 1. only use the tp
## 2. add ECs

s2_met_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  left_join(ec_tomodel, by = "environment") %>%
  mutate_at(vars(line_name, environment), as.factor) %>%
  # Grain fill tmin is not variable enough
  select(-grain_fill_mint) %>%
  # vegetative water stress is not variable enough
  select(-vegetative_water_stress) %>%
  # Vegetative tmin is not variable enough
  select(-vegetative_mint) %>%
  ## Remove CRM15 and CRM16 - outliers
  ## Remove plant height and CRM17 - outlier
  filter(!( trait == "GrainYield" & environment %in% c("CRM15", "CRM16") ),
         !( trait == "PlantHeight" & environment %in% c("CRM17") ))

## What is the level of coverage in this dataset?
s2_met_tomodel %>% 
  group_by(trait, line_name) %>%
  summarize(n = n_distinct(environment)) %>% 
  mutate(mean_env = n / max(n)) %>%
  top_n(x = ., n = 5, wt = -mean_env)






# New vector of covariates
environmental_covariates <- s2_met_tomodel %>% 
  select(matches("vegetative|flowering|grain_fill")) %>%
  names()

## Scatterplot matrix
scatterplotMatrix(distinct(select(s2_met_tomodel, environmental_covariates)), smooth = FALSE)


## Determine if there is sufficient variation for a covariate
distinct(select(s2_met_tomodel, environmental_covariates)) %>%
  map_dbl(var)

## Assign covariates for each trait
ec_by_trait <- list(
  GrainProtein = c("grain_fill_maxt", "grain_fill_water_stress"),
  GrainYield = c("flowering_mint", "grain_fill_maxt", "grain_fill_water_stress"),
  TestWeight = c("flowering_mint", "grain_fill_maxt", "grain_fill_water_stress"),
  PlantHeight = str_subset(environmental_covariates, pattern = "grain_fill", negate = TRUE) %>% 
    str_subset(string = ., pattern = "tmean|radn"),
  HeadingDate = str_subset(environmental_covariates, pattern = "grain_fill|flowering", negate = TRUE) %>% 
    str_subset(string = ., pattern = "tmean|radn")
) %>%
  tibble(trait = names(.), covariates = .)

# New vector of covariates
environmental_covariates <- reduce(ec_by_trait$covariates, union)




############################
### Model covariates that explain environment effect
############################

## Fit a base, random effect model for all traits
base_model_fit <- s2_met_tomodel %>%
  group_by(trait) %>%
  do(fit = lmer(value ~ 1 + (1|line_name) + (1|environment), data = .)) %>%
  ungroup()




## Get estimates of environmental effect
base_model_env_effect <- base_model_fit %>%
  mutate(effect = map(fit, ~ranef(.)$environment %>% rownames_to_column("environment") %>% rename(effect = 2))) %>%
  unnest(effect)
  
## Histogram
par(mfrow = c(2,3))
for (tr in unique(base_model_env_effect$trait)) {
  hist(subset(base_model_env_effect, trait == tr, effect, drop = T), main = tr)
}
par(mfrow = c(1,1))


# # environments as fixed
# base_model_env_effect_alt <- s2_met_tomodel %>%
#   group_by(trait) %>%
#   do(fit = lmer(value ~ 1 + (1|line_name) + environment, data = ., contrasts = list(environment = "contr.sum"))) %>%
#   ungroup() %>%
#   mutate(mf = map(fit, model.frame),
#          effect = map2(fit, mf, ~tibble(environment = levels(.y$environment), effect = c(fixef(.x)[-1], -sum(fixef(.x)[-1]))))) %>%
#   unnest(effect)
# 
# ## Compare
# left_join(base_model_env_effect, base_model_env_effect_alt, by = c("trait", "environment")) %>%
#   qplot(x = effect.x, y = effect.y, data = .) +
#   facet_wrap(~ trait, scales = "free") +
#   geom_abline(slope = 1, intercept = 0)
# 
# # Good! - proceed with random effect




## Use stepwise regression to find the best model
env_effect_to_model <- base_model_env_effect %>%
  left_join(s2_met_tomodel) %>%
  select(trait, environment, effect, environmental_covariates) %>%
  distinct()

# Group by trait
env_effect_models <- env_effect_to_model %>%
  group_by(trait) %>%
  nest() %>%
  mutate(model = list(NULL))

for (i in seq(nrow(env_effect_models))) {
  
  df <- env_effect_models$data[[i]]
  tr <- env_effect_models$trait[i]
  
  ## Remove some covariates depending on when the trait is measured
  df1 <- select(df, environment, effect, subset(ec_by_trait, trait == tr, covariates, drop = T)[[1]])
  
  # Minimal model
  min_model <- effect ~ 1
  # Max model
  max_model <- add_predictors(min_model, as.formula(paste0("~", paste0(names(df1)[-1:-2], collapse = " + "))))
  
  # Stepwise regression
  fit_base <- lm(formula = min_model, data = df1)
  fit_step <- step(object = fit_base, scope = max_model, direction = "both")
  
  ## Return the model
  env_effect_models$model[[i]] <- fit_step
  
}



## Plot
env_effect_models %>% 
  mutate(data = map2(data, model, ~mutate(.x, predicted_effect = predict(.y)))) %>%
  unnest(data) %>%
  ggplot(aes(x = predicted_effect, y = effect)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ trait, scales = "free") +
  theme_acs()






############################
### Model covariates that explain GxE
############################



#####
# Model building
#####





## Data.frame to store final models and covariates
ec_model_building <- s2_met_tomodel %>% 
  distinct(trait) %>%
  mutate(apriori_model = list(NULL), final_model = list(NULL))
  

### 
### Ad-hoc, a priori model approach
### 
### Use knowledge of traits to determine the covariates likely to explain the most variation
### in GxE
### 


## 
## Grain yield, grain protein, and test weight
## 
## Start with the following a prior covariates:
## 
## tmin during flowering
## radiation during vegetative
## water stress during flowering
## grain fill water stress
## grain fill tmax
## 
## 
## For heading date and plant height, just use all covariates and backwards elimination
## 
## 


## Loop over traits
for (i in seq(nrow(ec_model_building))) {
  
  # Get the trait
  tr <- ec_model_building$trait[i]
  
  # Subset the data
  tr_data <- filter(s2_met_tomodel, trait == tr)
  
  # Get the base model
  tr_base_model <- subset(base_model_fit, trait == tr, fit, drop = T)[[1]]
  
  ## Assign a priori covariates
  apriori_tr_covariates <- subset(ec_by_trait, trait == tr, covariates, drop = T)[[1]]
  
  ## Create a new formula
  tr_covariate_form <- as.formula(paste0("~", paste0(apriori_tr_covariates, collapse = " + "), " + ", 
                                         paste0("line_name:", apriori_tr_covariates, collapse = " + ")))
  # tr_covariate_form <- as.formula(paste0("~", paste0("(0 +", apriori_tr_covariates, "|line_name)", collapse = " + ")))

  gxe_formula <- add_predictors(formula(tr_base_model), tr_covariate_form)

  # Refit the model
  tr_ec_model_mixed <- lmer(formula = gxe_formula, data = tr_data)

  # Stepwise elimination of extract covariates
  # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed)
  tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed, reduce.random = FALSE, alpha.fixed = alpha)
  
  # Get the model
  tr_ec_model_mixed_step1 <- get_model(tr_ec_model_mixed_step)
  
  
  # # Adjust
  # model_keep <- tr_ec_model_mixed_step1
  # # model_keep <- lmer(formula = value ~ (1 | line_name) + (1 | environment) + line_name:flowering_mint + line_name:grain_fill_maxt + 
  # #                      line_name:vegetative_radn , data = tr_data)

  
  ## Add results to the tibble
  ec_model_building$apriori_model[i] <- list(tr_ec_model_mixed)
  ec_model_building$final_model[i] <- list(tr_ec_model_mixed_step1)
  
}



## Get the sensitivity coefficients for each covariate
covariate_reg_coefs <- ec_model_building %>%
  mutate_at(vars(contains("model")), 
            ~map(., ~fixef(.) %>%
                   tibble(term = names(.), estimate = .) %>% 
                   filter(str_detect(term, "line_name")) %>% 
                   mutate(term = str_remove_all(term, "line_name")) %>% 
                   separate(term, c("covariate", "line_name"), sep = ":") ))
                   

## Plot
g_coef_hist <- covariate_reg_coefs %>% 
  unnest(apriori_model) %>%
  ggplot(aes(x = estimate)) +
  geom_histogram() +
  facet_wrap(~ trait + covariate, scale = "free") +
  theme_presentation2(base_size = 8)
ggsave(filename = "apriori_covariate_coef.jpg", plot = g_coef_hist, path = fig_dir, height = 6, width = 6, dpi = 1000)

g_coef_hist <- covariate_reg_coefs %>% 
  unnest(final_model) %>%
  ggplot(aes(x = estimate)) +
  geom_histogram() +
  facet_wrap(~ trait + covariate, scale = "free") +
  theme_presentation2(base_size = 8)
ggsave(filename = "final_covariate_coef.jpg", plot = g_coef_hist, path = fig_dir, height = 6, width = 6, dpi = 1000)


## Plot phenotypic value against covariate
s2_met_tomodel %>%
  select(trait, environment, line_name, value, environmental_covariates) %>%
  gather(covariate, env_value, environmental_covariates) %>%
  ## Add top and bottom coefficients from regression
  left_join(group_by(unnest(covariate_reg_coefs, apriori_model), trait, covariate) %>% mutate(rank = min_rank(estimate)) %>% 
              filter(rank %in% c(1:5, c(max(rank) - 4, max(rank)))), ., by = c("trait", "line_name", "covariate")) %>%
  ggplot(aes(x = env_value, y = value, group = line_name)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ trait + covariate, scales = "free")+
  theme_presentation2(base_size = 8)
ggsave(filename = "geno_coef_response_apriori.jpg", path = fig_dir, height = 12, width = 8, dpi = 1000)

## Plot phenotypic value against covariate
s2_met_tomodel %>%
  select(trait, environment, line_name, value, environmental_covariates) %>%
  gather(covariate, env_value, environmental_covariates) %>%
  ## Add top and bottom coefficients from regression
  left_join(group_by(unnest(covariate_reg_coefs, final_model), trait, covariate) %>% mutate(rank = min_rank(estimate)) %>% 
              filter(rank %in% c(1:5, c(max(rank) - 4, max(rank)))), ., by = c("trait", "line_name", "covariate")) %>%
  ggplot(aes(x = env_value, y = value, group = line_name)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ trait + covariate, scales = "free")+
  theme_presentation2(base_size = 8)
ggsave(filename = "geno_coef_response_final.jpg", path = fig_dir, height = 12, width = 8, dpi = 1000)



## Save these results
save("ec_model_building", "env_effect_models", "covariate_reg_coefs", "s2_met_tomodel", "ec_tomodel", "ec_tomodel_centers", 
     file = file.path(result_dir, "ec_model_building.RData"))


  


############################
### Cross-validation
############################



## Leave one environment out
loeo_model_cv <- s2_met_tomodel %>%
  distinct(trait, environment) %>%
  mutate(out = list(NULL))

## Iterate over rows
for (i in seq(nrow(loeo_model_cv))) {
  
  # Get the environment to EXCLUDE
  env_drop <- as.character(loeo_model_cv$environment[i])
  # Trait
  tr <- loeo_model_cv$trait[i]
  
  # Mask phenotypic data for that environment
  train_data <- s2_met_tomodel %>%
    filter(trait == tr) %>%
    mutate(value = ifelse(environment == env_drop, NA, value)) %>%
    droplevels()
  
  ## Extract the base model
  base_model_tr <- subset(base_model_fit, trait == tr, fit, drop = T)[[1]]
  # Extract the environment effect model
  env_effect_model_tr <- subset(env_effect_models, trait == tr, model, drop = T)[[1]]
  # Extract the gxe ec model
  ec_model_tr <- subset(ec_model_building, trait == tr, final_model, drop = T)[[1]]
  
  ## Create new formulae - all fixed
  base_model_form <- as.formula(paste0("value ~ ", paste0(all.vars(formula(base_model_tr)[[3]]), collapse = " + ")))
  ec_model_form <- as.formula(paste0("value ~ ", paste0(str_remove(attr(terms(formula(ec_model_tr)), "term.labels"), "1 \\| "), collapse = " + ")))
    
  
  # model 1 is only main effects
  model1 <- base_model_form
  # Model 2 is main effects with environment represented by covariates
  model2 <- add_predictors(f = as.formula(paste0("value ~ 1 + ", paste0(all.vars(model1[[3]])[-2] , collapse = " + "))), 
                           as.formula(paste0("~ ", paste0(all.vars(formula(env_effect_model_tr)[[3]]), collapse = " + ")))) 
  # Model 3 is main effects with ECs
  model3 <- ec_model_form
  # Model 4 is main effects + environment by covariates + ECs
  model4 <- add_predictors(model2, formula(drop.terms(terms(model3), 1:2)))
  
  # Create a list
  models <- ls(pattern = "model[0-9]{1,}")
  model_list <- map(setNames(models, models), get)
  
  # Fit all
  fit_list <- fit_with(data = train_data, .f = lm, .formulas = model_list)
  ## Adjust levels of environments
  fit_list1 <- map(fit_list, ~{
    if ("environment" %in% names(.$xlevels)) {
      .$xlevels$environment <- levels(train_data$environment)
    }
    return(.)
  })
  
  ## Prediction environment - ECs for prediction
  pred_data <- subset(s2_met_tomodel, trait == tr & environment == env_drop)
  # Predict
  prediction_list <- map(fit_list1, ~{
    # Get the coefficients
    mod <- .
    cf <- coef(mod)
    
    # Create mm formula
    if ("environment" %in% all.vars(formula(mod))) {
      mm_form <- terms(formula(mod)) %>% 
        drop.terms(termobj = ., dropx = which(attr(., "term.labels") == "environment")) %>% 
        formula()
    } else {
      mm_form <- formula(mod)
      
    }
    X <- model.matrix(mm_form, pred_data, )
    
    cf_use <- which(names(cf) %in% colnames(X))
    # convert NA to 0
    cf[is.na(cf)] <- 0
    # Return y_hat
    c(X %*% cf[cf_use])
    
  })
                           
  ## Return predictions
  pred_data1 <- pred_data %>% 
    select(line_name:value) %>%
    bind_cols(., as_tibble(prediction_list))
  
  loeo_model_cv$out[[i]] <- pred_data1
  
}

  
## Summarize accuracy
loeo_model_pred <- loeo_model_cv %>% 
  unnest(out) %>%
  gather(model, prediction, starts_with("model"))

correlation_accuracy <- loeo_model_pred %>% 
  group_by(trait, environment, model) %>% 
  summarize(accuracy = cor(value, prediction),
            bias = mean(prediction - value)) %>%
  ungroup()

## Calculate average
correlation_accuracy %>%
  group_by(trait, model) %>%
  summarize_at(vars(accuracy, bias), mean)




############################
### Heritability
############################


## Fit genomic heritability models for slopes
# Extract the coefficients for each model
# Keep final and apriori models
ec_interaction_coef <- covariate_reg_coefs %>%
  rename_all(~str_remove(., "_model")) %>%
  gather(model, out, -trait) %>% 
  unnest(out)


# Subset the K matrix
K_use <- K[tp_geno, tp_geno]

## Fit models
ec_interaction_coef_herit <- ec_interaction_coef %>%
  group_by(trait, model, covariate) %>%
  do({
    df <- .
    
    df1 <- subset(df, line_name %in% tp_geno)
    
    ## Calculate heritability
    invisible(capture.output(herit_fit <- marker_h2(data.vector = df1$estimate, geno.vector = df1$line_name, K = K_use, alpha = alpha)))

    ## Return df
    tibble(heritability = herit_fit$h2, lower = herit_fit$conf.int1[1], upper = herit_fit$conf.int1[2])
    
  })

## Write table
write_csv(x = ec_interaction_coef_herit, path = file.path(fig_dir, "covariate_slope_heritability.csv"))











### Calculation of daily temperature and precipitation stats

# Load both climate data
load(file.path(data_dir, "Climate_Data/NOAA_Data/noaa_stations_trial_data.RData"))

### One-year environmental covariates
oneyear_env_data_unnest <- noaa_trial_data_oneyear_complete %>%
  unnest(data)

## Add trial data
oneyear_trial_env_data <- trial_info %>% 
  distinct(environment, location, year, planting_date) %>%
  mutate_all(parse_guess) %>%
  left_join(., oneyear_env_data_unnest) %>%
  mutate_at(vars(planting_date, date), ymd) %>%
  # Convert tenth of value to value
  mutate(value = value / 10)


## Subset data only after planting, then add 30-day intervals
oneyear_trial_env_data_use <- oneyear_trial_env_data %>% 
  select(environment, planting_date, date, month, datatype, value) %>% 
  filter(date >= planting_date) %>% 
  mutate(dap = as.numeric(date - planting_date))

## Calculate some transformed stats
## Range = max - min
## Avg = mean
one_year_daily_stats <- oneyear_trial_env_data_use %>% 
  # Filter for environments used in this analysis
  filter(environment %in% unique(S2_MET_BLUEs$environment)) %>%
  select(environment, dap, datatype, value) %>% 
  spread(datatype, value) %>% 
  mutate(TAVG = (TMAX + TMIN) / 2,
         TRANGE = TMAX - TMIN)

## Remove daps after which not all environments have data
max_dap <- one_year_daily_stats %>% 
  group_by(dap) %>% 
  summarize(nE = n_distinct(environment)) %>%
  filter(nE != n_distinct(one_year_daily_stats$environment)) %>% 
  pull(dap) %>% 
  min()



## For each environment, create intervals of 0 dap to max(dap) in steps of 1 day
## Calculate the average of all daily ECs within each interval
intervals <- crossing(start = seq(0, max(one_year_daily_stats$dap)), end = start) %>%
  filter(end >= start, end <= max_dap)

one_year_interval_ecs <- intervals %>%
  mutate(data = map2(start, end, ~subset(one_year_daily_stats, between(dap, .x, .y)) %>% 
                       group_by(environment) %>% summarize_at(vars(PRCP, TMAX, TMIN, TAVG, TRANGE), mean, na.rm = T))) %>%
  unnest(data)


## For each trait, calculate environmental means
## Fit a general linear model with weights
## Only use TP data
main_effect_model <- S2_MET_BLUEs %>%
  filter(line_name %in% tp, trait %in% traits) %>%
  group_by(trait) %>%
  do({
    
    df <- .
    # Contrasts
    df1 <- mutate_at(df, vars(line_name, environment), as.factor)
    
    contrasts(df1$line_name) <- contr.sum(levels(df1$line_name)) %>%
      `colnames<-`(., head(row.names(.), -1))
    contrasts(df1$environment) <- contr.sum(levels(df1$environment)) %>%
      `colnames<-`(., head(row.names(.), -1))
    
    # Fit the linear model
    fit <- lm(value ~ 1 + line_name + environment, data = df1, weights = 1 / df1$std_error)
    
    # Get environmenal effects (and means)
    tidy(fit) %>%
      filter(str_detect(term, "environment")) %>%
      mutate(term = str_remove(term, "environment")) %>%
      select(environment = term, effect = estimate) %>%
      add_row(environment = tail(levels(df1$environment), 1), effect = -sum(.$effect)) %>%
      # Add intercept
      mutate(mean = effect + coef(fit)[1])
    
  }) %>% ungroup()



## For each interval, correlate the observed EC in each environment with the environmental effect
environment_effect_and_ec <- full_join(main_effect_model, one_year_interval_ecs)

## Calculate correlation coefficients
interval_correlations <- environment_effect_and_ec %>%
  group_by(trait, start, end) %>%
  summarize_at(vars(PRCP, TMAX, TMIN, TAVG, TRANGE), funs(cor(., effect))) %>%
  ungroup() %>%
  gather(covariate, cor, -trait, -start, -end)


## Visualize via heatmap
interval_correlations %>%
  ggplot(aes(x = start, y = end, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_grid(trait ~ covariate)


## Grab the interval with the greatest abs(correlation) for each trait






## Fit a bunch of models
interval_models <- environment_effect_and_ec %>%
  mutate(HEAT = TAVG - TMIN) %>%
  group_by(trait, start, end) %>%
  do(fit = lm(effect ~ PRCP * HEAT, data = .))

## Extract R2
interval_models_R2 <- interval_models %>%
  ungroup() %>%
  mutate(R2 = map_dbl(fit, ~summary(.)$r.square))

## Visualize via heatmap
interval_models_R2 %>%
  ggplot(aes(x = start, y = end, fill = R2)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_grid(~ trait)

## Find the model with the largest R2
best_models <- interval_models_R2 %>%
  group_by(trait) %>%
  filter(R2 == max(R2))










