## S2MET genotype-environment ecophysiological analysis
## 
## Calculate environmental covariables from a crop model output
## 
## Author: Jeff Neyhart
## Last modified: 4 October  2019
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





# Assess accuracy of heading date predictions from multiple cultivar  --------



load(file.path(enviro_dir, "GrowthStaging/apsim_s2met_model_cultivar_test_results.RData"))

## Reorganize output
apsim_s2met_cultivar_out1 <- apsim_s2met_cultivar_out %>% 
  unnest(apsim_out) %>%
  # Assign cultivar
  mutate(cultivar = str_extract(simulation_name, "cultivar1\\=.*$") %>% str_remove(., "cultivar1=")) %>%
  mutate(stage = case_when(
    between(zadok_stage, 10, 30) ~ "early_vegetative",
    between(zadok_stage, 30, 50) ~ "late_vegetative",
    between(zadok_stage, 50, 70) ~ "flowering",
    between(zadok_stage, 70, 91) ~ "grain_fill")) %>%
  filter(sowing_das == 1) %>% 
  # Sort by trial, cultivar, dap
  arrange(trial, cultivar, day) %>%
  group_by(trial, cultivar) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(., dap = seq(nrow(.))))) %>%
  unnest() %>%
  select(trial, environment, cultivar, date = Date, day, dap, stage)


# Calculate avereage DAP of heading date predictions from CGM
cgm_predicted_heading <- apsim_s2met_cultivar_out1 %>% 
  filter(stage %in% c("heading", "flowering")) %>% 
  group_by(trial, environment, cultivar, stage) %>% 
  summarize(pred_HD = mean(dap)) %>%
  # summarize(pred_HD = min(dap)) %>%
  ungroup() %>%
  # Remove S2C1 trials
  filter(str_detect(trial, "S2C1", negate = TRUE))

## Combine
pred_obs_heading <- inner_join(env_means_heading, cgm_predicted_heading)

## Summarize
cgm_pred_HD_summary <- pred_obs_heading %>%
  group_by(cultivar, stage) %>%
  summarize(cor = cor(obs_HD, pred_HD), 
            mae = mean(abs(pred_HD - obs_HD)),
            rmse = sqrt(mean((pred_HD - obs_HD)^2)))

## Find the cultivar that results in the most accuracy prediction
cgm_pred_HD_summary %>%
  filter(stage == "flowering") %>%
  # arrange(desc(cor))
  arrange(rmse) %>%
  as.data.frame()

## Create an annotation df
cgm_pred_HD_summary_annotate <- cgm_pred_HD_summary %>%
  filter(stage == "flowering") %>%
  mutate(cor = paste0("R^2==", round(cor^2, 3)),
         rmse = paste0("RMSE==", round(rmse, 3)))

## Plot each cultivar
(g_cgm_cultivar_summary <- pred_obs_heading %>%
  filter(stage == "flowering") %>%
  ggplot(aes(y = obs_HD, x = pred_HD)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ cultivar) +
  scale_y_continuous(name = "Observed mean heading (dap)", breaks = pretty) +
  scale_x_continuous(name = "CGM predicted flowering (dap)", breaks = pretty)  +
  geom_text(data = cgm_pred_HD_summary_annotate, 
            aes(x = 89, y = 54, label = cor), parse = T, hjust = 1, size = 2) +
  geom_text(data = cgm_pred_HD_summary_annotate, 
            aes(x = 89, y = 52, label = rmse), parse = T, hjust = 1, size = 2) +
  theme_presentation2(10))

ggsave(filename = "figure_SXX_cgm_cultivar_testing.jpg", plot = g_cgm_cultivar_summary, path = fig_dir,
       height = 8, width = 10, dpi = 500)




# Assess accuracy of heading date predictions from the chosen crop --------


## Load the environmental covariates
## This contains output from both nasapower and dayment
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))




# Calculate avereage DAP of heading date predictions from CGM
cgm_predicted_heading <- growth_stage_weather %>% 
  filter(stage %in% c("heading", "flowering")) %>% 
  group_by(source, trial, environment, stage) %>% 
  rename(pred_HD = dap) %>%
  summarize_at(vars(pred_HD, day), list(mean = ~mean, min = ~min)) %>%
  ungroup() %>%
  # Remove S2C1 trials
  filter(str_detect(trial, "S2C1", negate = TRUE))


## Combine
pred_obs_heading <- inner_join(env_means_heading, cgm_predicted_heading)

## Summarize
(cgm_pred_HD_summary <- pred_obs_heading %>%
  group_by(source, stage) %>%
  summarize_at(vars(pred_HD_mean, pred_HD_min), list(
    cor = ~cor(obs_HD, .), mae = ~mean(abs(. - obs_HD)), rmse = ~sqrt(mean((. - obs_HD)^2))
  )))


## Plot by stage
g_pred_obs_HD <- pred_obs_heading %>%
  ggplot(aes(x = pred_HD_mean, y = obs_HD, label = environment)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ source + stage) +
  theme_presentation2(12)

## Replot and save
cgm_pred_HD_summary_annotate <- cgm_pred_HD_summary %>%
  mutate(pred_HD_mean_cor = paste0("R^2==", round(pred_HD_mean_cor^2, 3)),
         pred_HD_mean_rmse = paste0("RMSE==", round(pred_HD_mean_rmse, 3)))

g_cgm_summary <- pred_obs_heading %>%
  filter(stage == "flowering", source == "daymet") %>%
  ggplot(aes(y = obs_HD, x = pred_HD_mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  scale_y_continuous(name = "Observed mean heading (dap)", breaks = pretty) +
  scale_x_continuous(name = "CGM predicted flowering (dap)", breaks = pretty) +
  geom_text(data = filter(cgm_pred_HD_summary_annotate, stage == "flowering", source == "daymet"), 
            aes(x = 51, y = 72, label = pred_HD_mean_cor), parse = T, hjust = 0) +
  geom_text(data = filter(cgm_pred_HD_summary_annotate, stage == "flowering", source == "daymet"), 
            aes(x = 51, y = 70, label = pred_HD_mean_rmse), parse = T, hjust = 0) +
  theme_genetics()

ggsave(filename = "figure_SXX_cgm_pred_obs_HD.jpg", plot = g_cgm_summary, path = fig_dir,
       height = 3, width = 3, dpi = 1000)
  
  




# Prepare environment-concurrent covariates -------------------------------



## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))


## Compare sources
growth_stage_covariates %>%
  gather(covariate, value, -trial, -stage, -source) %>%
  spread(source, value) %>%
  ggplot(aes(x = daymet, y = nasapower)) +
  geom_point() +
  facet_wrap(~ covariate + stage, scales = "free", ncol = n_distinct(growth_stage_covariates$stage)) +
  theme_genetics()





## Create the ECs and select the relevant ones for modeling
growth_stage_covariates1 <- growth_stage_covariates %>%
  ## Remove heading as a growth stage
  filter(stage != "heading") %>%
  inner_join(., env_trials, by = "trial") %>% 
  select(-trial) %>%
  gather(covariate, value, -environment, -stage, -source) %>%
  unite(covariate, stage, covariate, sep = ".") %>%
  ## Remove trange, relhum, and rain (water stress will cover this)
  filter(str_detect(covariate, "radn_mean", negate = TRUE)) %>%
  spread(covariate, value)

soil_covariates1 <- soil_covariates %>%
  inner_join(., env_trials, by = "trial") %>%
  select(-trial)

# Add soil covariates
ec_select <- full_join(growth_stage_covariates1, soil_covariates1, by = "environment")


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
subset(ec_tomodel_normality, p_value < 0.10 * (0.1 / n_distinct(ec_tomodel_normality$covariate)))

## Which ECs should be kept
ec_to_keep <- ec_tomodel_normality %>%
  filter(p_value >= 0.10 * (0.1 / n_distinct(ec_tomodel_normality$covariate))) %>%
  select(source, covariate)



# Remove some covariate and impute with the mean
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



# Load the covariate data
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_historical_environmental_covariates.RData"))

## Environments to use with corresponding trials
loc_trials <- S2_MET_BLUEs %>% 
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "AID|HTM")) %>%
  filter(environment %in% c(train_test_env, validation_env),
         str_detect(trial, "S2C1", negate = TRUE)) %>%
  distinct(trial, location)



## Create the ECs and select the relevant ones for modeling
growth_stage_covariates1 <- growth_stage_covariates %>%
  ## Remove heading as a growth stage
  filter(stage != "heading") %>%
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
  inner_join(., loc_trials, by = "location") %>%
  gather(covariate, value, -trial, -location) %>%
  left_join(distinct(growth_stage_covariates1, year, trial, location), .) %>%
  # Cross with sources of weather data
  crossing(., source = unique(growth_stage_covariates1$source))

## Combine
hist_ec_select <- bind_rows(growth_stage_covariates1, soil_covariates1) %>%
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


