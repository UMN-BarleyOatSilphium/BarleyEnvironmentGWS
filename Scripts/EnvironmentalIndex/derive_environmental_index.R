## S2MET genotype-environment ecophysiological analysis
## 
## Determine environmental indices
## 
## Author: Jeff Neyhart
## Last modified: 9 June 2019
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

# Significance level
alpha <- 0.05


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










