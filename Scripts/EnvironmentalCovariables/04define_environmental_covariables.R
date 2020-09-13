## S2MET Genomewide Environment Predictions
## 
## Calculate environmental covariables from a crop model output
## 


# Repository directory
repo_dir <- getwd()
# Source the main project script
source(file.path(repo_dir, "source.R"))


## Add new packages to load
pkgs <- union(pkgs, c("modelr", "broom", "lubridate", "car", "car", "ggrepel", "cowplot"))
# Load these packages
invisible(lapply(X = pkgs, library, character.only = TRUE))

# Significance level
alpha <- 0.05

# Correlations cutoff
max_cor <- 0.90



## Environments to use with corresponding trials
env_trials <- S2_MET_BLUEs %>% 
  distinct(trial, environment) %>%
  filter(environment %in% c(train_test_env, validation_env),
         str_detect(trial, "S2C1", negate = TRUE))



## Fit models to calculate environmental means
env_means <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp, vp)) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Factorize
    df1 <- mutate_at(df, vars(line_name, environment), ~fct_contr_sum(as.factor(.))) %>%
      mutate(weight = std_error^2)
    
    # Fit the model
    fit <- lmer(value ~ (1|line_name) + environment, data = df1, weights = weight)
    
    ## Return a df of environmental effects
    fixef(fit) %>% 
      tibble(environment = names(.), effect = .) %>% 
      filter(environment != "(Intercept)") %>% 
      mutate(environment = str_remove(environment, "environment"),
             mu = fixef(fit)[1]) %>% 
      add_row(environment = last(levels(df1$environment)), effect = -sum(.$effect), 
              mu = .$mu[1])
    
  }) %>% ungroup()

# Just heading date
env_means_heading <- subset(env_means, trait == "HeadingDate") %>%
  mutate(obs_HD = mu + effect)



# Prepare environment-concurrent covariates -------------------------------


## Load the environmental covariates
load(file.path(result_dir, "environmental_covariates.RData"))



## Create the ECs and select the relevant ones for modeling
concurrent_growth_stage_covariates1 <- concurrent_growth_stage_covariates %>%
  inner_join(., env_trials, by = "trial") %>% 
  select(-trial) %>%
  gather(covariate, value, -environment, -stage, -source) %>%
  unite(covariate, stage, covariate, sep = ".") %>%
  ## Remove trange, relhum, and rain (water stress will cover this)
  filter(str_detect(covariate, "radn_mean", negate = TRUE)) %>%
  spread(covariate, value)

soil_covariates1 <- soil_covariates %>%
  inner_join(., env_trials, by = "trial") %>%
  select(-trial, -location)

# Add soil covariates
ec_select <- full_join(concurrent_growth_stage_covariates1, soil_covariates1, by = "environment")


## Summarize min/max/var for each covariate
ec_select_summ <- ec_select %>%
  gather(covariate, value, -environment, -source) %>%
  group_by(covariate, source) %>%
  summarize_at(vars(value), list(min = min, max = max, sd = sd, cv = cv), 
               na.rm = T) %>%
  ungroup()

## Remove covariates with low CV
ec_select1 <- ec_select_summ %>% 
  filter(!is.na(cv)) %>%
  select(source, covariate) %>%
  inner_join(., gather(ec_select, covariate, value, -environment, -source)) %>%
  spread(covariate, value)


### Covariate diagnostics ###

ec_select2 <- ec_select1

## Test for normality using ks test
ec_tomodel_normality <- ec_select2 %>%
  gather(covariate, value, -environment, -source) %>%
  group_by(source, covariate) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value, na.rm = T), sd = sd(.$value, na.rm = T))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value"),
         p_adj = p.adjust(p = p_value, method = "bonf"))

## Which covariates fail?
(ec_sign_nonnormal <- subset(ec_tomodel_normality, p_value < 0.10 * (0.1 / n_distinct(ec_tomodel_normality$covariate))))

# Remove the ECs that show significant departures from normality
ec_to_keep <- ec_tomodel_normality %>%
  anti_join(., select(ec_sign_nonnormal, source, covariate)) %>%
  select(source, covariate)



# Impute with the mean
ec_select3 <- ec_select2 %>%
  gather(covariate, value, -environment, -source) %>%
  inner_join(., ec_to_keep) %>%
  spread(covariate, value) %>%
  split(.$source) %>%
  # Inpute using the mean
  map_df(., ~mutate_at(., vars(-environment, -source), impute))


## Calculate pairwise correlations
ec_tidy_cor <- ec_select3 %>%
  split(.$source) %>%
  imap_dfr(~select(., -environment, -source) %>%
        cor() %>%
        as.dist() %>%
        tidy() %>%
        rename_all(~c("covariate1", "covariate2", "correlation")) %>%
        mutate_at(vars(-correlation), as.character) %>%
        mutate(source = .y)
  )

## Find those pairs of covariates with inflated correlations, above 0.9
(ec_inflated_cor <- ec_tidy_cor %>% 
  filter(abs(correlation) > max_cor))

# covariate1                  covariate2                 correlation source   
# 1 early_vegetative.tmean_mean early_vegetative.maxt_mean       0.919 daymet   
# 2 early_vegetative.radn_sum   early_vegetative.mint_mean      -0.909 daymet   
# 3 early_vegetative.tmean_mean early_vegetative.mint_mean       0.937 daymet   
# 4 subsoil_teb                 subsoil_pH_h2o                   0.924 daymet   
# 5 topsoil_teb                 subsoil_teb                      0.908 daymet   
# 6 subsoil_teb                 subsoil_pH_h2o                   0.924 nasapower
# 7 topsoil_teb                 subsoil_teb                      0.908 nasapower

# Remove nothing 

ec_select4 <- ec_select3 # %>% select(-subsoil_teb, -subsoil_ref_bulk_density)



## Prepare the covariates for modeling
## Center, but do not scale. Save the mean for later
ec_tomodel_temp <- ec_select4 %>%
  split(.$source) %>%
  map(~mutate_at(., vars(-environment, -source), scale, scale = FALSE, center = TRUE))

ec_tomodel_centers <- ec_tomodel_temp %>%
  map(~summarize_at(., vars(-environment, -source), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
ec_tomodel_centered <- ec_tomodel_temp %>%
  map(~mutate_at(., vars(-environment, -source), as.numeric))

## Center and scale
ec_tomodel_scaled <- ec_select4 %>%
  split(.$source) %>%
  map(~mutate_at(., vars(-environment, -source), ~as.numeric(scale(., scale = TRUE, center = TRUE))))







# Prepare historical covariates -------------------------------------------

## Environments to use with corresponding trials
loc_trials <- S2_MET_BLUEs %>% 
  # Remove irrigated trials
  filter(!str_detect(environment, "AID|HTM")) %>%
  filter(environment %in% c(train_test_env, validation_env),
         str_detect(trial, "S2C1", negate = TRUE)) %>%
  distinct(trial, location)



## Create the ECs and select the relevant ones for modeling
historical_location_growth_stage_covariates1 <- historical_location_growth_stage_covariates %>%
  ## Separate trial into location and year
  mutate(year = as.numeric(str_extract(trial, "[0-9]{4}")),
         location = str_remove(trial, "_[0-9]{4}")) %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca"),
         ## Replace Columbus with Wooster
         location = ifelse(location == "Columbus", "Wooster", location)) %>%
  select(-trial) %>%
  inner_join(., loc_trials, by = "location") %>% 
  # Filter for years before those observed in the S2MET project
  filter(year < min(S2_MET_BLUEs$year)) %>%
  gather(covariate, value, -trial, -stage, -location, -year, -source) %>%
  unite(covariate, stage, covariate, sep = ".")

soil_covariates1 <- soil_covariates %>%
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca"),
         ## Replace Columbus with Wooster
         location = ifelse(location == "Columbus", "Wooster", location)) %>%
  # Distinct locations
  distinct_at(vars(-trial)) %>%
  inner_join(., loc_trials, by = "location") %>%
  gather(covariate, value, -trial, -location) %>%
  left_join(distinct(historical_location_growth_stage_covariates1, year, trial, location), .) %>%
  # Cross with sources of weather data
  crossing(., source = unique(historical_location_growth_stage_covariates1$source))

## Combine
hist_ec_select <- bind_rows(historical_location_growth_stage_covariates1, soil_covariates1) %>%
  # Select only those covariates used for individual environments
  filter(covariate %in% names(ec_select4))

## Create timeframes ##
#
# Create two sets of time intervals:
# 1. Increasing number of years
# 2. Same number of years (5, 10, and 15) in a sliding window 
# 

# Create the vector of years
all_historical_years <- seq(max(hist_ec_select$year) - 30 + 1, max(hist_ec_select$year))

# Create accumulating years
year_time_frame <- accumulate(rev(all_historical_years), c)

# Create the sliding windows
year_sliding_window <- c(window5 = 5, window10 = 10, window15 = 15) %>%
  map(~slide::slide(.x = all_historical_years, ~ ., .before = floor(.x / 2), .after = floor(.x / 2), .step = 1, .complete = TRUE)) %>%
  map(~subset(., !sapply(., is.null)))

## Use the time frames to summarize covariates
ec_select_timeframe <- tibble(
  time_frame = paste0("time_frame", map_chr(year_time_frame, ~paste(length(.), min(.), max(.), sep = "_"))),
  years = year_time_frame) %>%
  # Subset the historical EC
  mutate(ec_data = map(years, ~subset(hist_ec_select, year %in% .x)))

# Summarize - calculate the mean covariate for each location across years
ec_select_timeframe_summary <- ec_select_timeframe %>%
  unnest(ec_data) %>%
  group_by(source, time_frame, covariate, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

## Use sliding windows to summarize covariates
ec_select_window <- tibble(
  time_frame = unlist(imap(year_sliding_window, function(l, name) map_chr(l, ~paste0(name, "_", min(.x), "_", max(.x))))),
  years = unlist(year_sliding_window, recursive = FALSE)) %>%
  # Subset the historical EC
  mutate(ec_data = map(years, ~subset(hist_ec_select, year %in% .x)))

# Summarize
ec_select_window_summary <- ec_select_window %>%
  unnest(ec_data) %>%
  group_by(source, time_frame, covariate, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()



## Spread the data
ec_select_timeframe_summary_wide <- spread(ec_select_timeframe_summary, covariate, value)
ec_select_window_summary_wide <- spread(ec_select_window_summary, covariate, value)




### Covariate diagnostics ###


# Number of distinct covariates (number of independent tests)
nTests <- n_distinct(ec_select_timeframe_summary$covariate)

## Test for normality using ks test
ec_timeframe_normality <- ec_select_timeframe_summary %>%
  group_by(source, covariate, time_frame) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value, na.rm = T), sd = sd(.$value, na.rm = T))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value"),
         p_adj = p.adjust(p = p_value, method = "bonf"))

# Find covariates that failed
subset(ec_timeframe_normality, p_value < (0.1 / nTests))

## None to remove


# Repeat for window
## Test for normality using ks test
ec_window_normality <- ec_select_window_summary %>%
  group_by(source, covariate, time_frame) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value, na.rm = T), sd = sd(.$value, na.rm = T))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value"),
         p_adj = p.adjust(p = p_value, method = "bonf"))

# Find covariates that failed
subset(ec_window_normality, p_value < (0.1 / nTests))

# None to remove


# Impute with the mean
ec_select_timeframe_summary_wide1 <- ec_select_timeframe_summary_wide %>%
  split(.$source) %>%
  # Inpute using the mean
  map_df(., ~mutate_at(., vars(-time_frame, -source, -location), impute))


ec_select_window_summary_wide1 <- ec_select_window_summary_wide %>%
  split(.$source) %>%
  # Inpute using the mean
  map_df(., ~mutate_at(., vars(-time_frame, -source, -location), impute))


## Calculate pairwise correlations ##
# Timeframe
ec_timeframe_tidy_cor <- ec_select_timeframe_summary_wide1 %>%
  select(-location) %>%
  split(list(.$source, .$time_frame)) %>%
  imap_dfr(~cor(.x[,-1:-2]) %>% as.dist() %>% tidy() %>% 
         rename_all(~c("covariate1", "covariate2", "correlation")) %>%
         mutate_at(vars(-correlation), as.character) %>%
         mutate(time_frame = .y))

## Find those pairs of covariates with inflated correlations, above 0.9
ec_timeframe_inflated_cor <- ec_timeframe_tidy_cor %>% 
  filter(abs(correlation) > max_cor)



## Many of the correlated covariates are the same covariate between growth stages,
## which suggests that the boundaries between growth stages are being blurred
## 

# Window
ec_window_tidy_cor <- ec_select_window_summary_wide1 %>%
  select(-location) %>%
  split(list(.$source, .$time_frame)) %>%
  imap_dfr(~cor(.x[,-1:-2]) %>% as.dist() %>% tidy() %>% 
             rename_all(~c("covariate1", "covariate2", "correlation")) %>%
             mutate_at(vars(-correlation), as.character) %>%
             mutate(time_frame = .y))

## Find those pairs of covariates with inflated correlations, above 0.9
ec_window_inflated_cor<- ec_window_tidy_cor %>% 
  filter(abs(correlation) > max_cor)


# Same problem




## Prepare the covariates for modeling

## Time frame ##
## 
## Center, but do not scale. Save the mean for later
historical_ec_tomodel_timeframe_temp <- ec_select_timeframe_summary_wide1 %>%
  split(list(.$source, .$time_frame)) %>%
  map(~mutate_at(.x, vars(-source, -time_frame, -location), scale, scale = FALSE, center = TRUE))

historical_ec_tomodel_timeframe_centers <- historical_ec_tomodel_timeframe_temp %>%
  map(~summarize_at(.x, vars(-source, -time_frame, -location), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
historical_ec_tomodel_timeframe_centered <- historical_ec_tomodel_timeframe_temp %>%
  map(~mutate_at(.x, vars(-time_frame, -location, -source), as.numeric) )

## Center and scale
historical_ec_tomodel_timeframe_scaled <- ec_select_timeframe_summary_wide1 %>%
  split(list(.$source, .$time_frame)) %>%
  map(~mutate_at(.x, vars(-time_frame, -location, -source), scale, scale = TRUE, center = TRUE) %>%
        mutate_at(vars(-time_frame, -location, -source), as.numeric) )


## Window ##
## 
## Center, but do not scale. Save the mean for later
historical_ec_tomodel_window_temp <- ec_select_window_summary_wide1 %>%
  split(list(.$source, .$time_frame)) %>%
  map(~mutate_at(.x, vars(-time_frame, -location, -source), scale, scale = FALSE, center = TRUE))

historical_ec_tomodel_window_centers <- historical_ec_tomodel_window_temp %>%
  map(~summarize_at(.x, vars(-time_frame, -location, -source), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
historical_ec_tomodel_window_centered <- historical_ec_tomodel_window_temp %>%
  map(~mutate_at(.x, vars(-time_frame, -location, -source), as.numeric) )

## Center and scale
historical_ec_tomodel_window_scaled <- ec_select_window_summary_wide1 %>%
  split(list(.$source, .$time_frame)) %>%
  map(~mutate_at(.x, vars(-time_frame, -location, -source), scale, scale = TRUE, center = TRUE) %>%
        mutate_at(vars(-time_frame, -location, -source), as.numeric) )




## Save these results
save("ec_tomodel_centered", "ec_tomodel_scaled", "ec_tomodel_centers", 
     "historical_ec_tomodel_timeframe_centered", "historical_ec_tomodel_timeframe_scaled", "historical_ec_tomodel_timeframe_centers",
     "historical_ec_tomodel_window_centered", "historical_ec_tomodel_window_scaled", "historical_ec_tomodel_window_centers",
     file = file.path(result_dir, "concurrent_historical_covariables.RData"))


