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

## Environments to use with corresponding trials
env_trials <- S2_MET_BLUEs %>% 
  distinct(trial, environment) %>%
  filter(environment %in% tp_vp_env,
         str_detect(trial, "S2C1", negate = TRUE))




############################
### Assess accuracy of heading date predictions from multi-cultivar crop models
############################

load(file.path(enviro_dir, "GrowthStaging/apsim_s2met_model_cultivar_test_results.RData"))

## Reorganize output
apsim_s2met_cultivar_out1 <- apsim_s2met_cultivar_out %>% 
  unnest(apsim_out) %>%
  # Assign cultivar
  mutate(cultivar = str_extract(out_name, "cultivar1\\=.*$") %>% str_remove(., "cultivar1=")) %>%
  mutate(stage = case_when(
    between(zadok_stage, 10, 30) ~ "early_vegetative",
    between(zadok_stage, 30, 50) ~ "late_vegetative",
    between(zadok_stage, 50, 60) ~ "heading",
    between(zadok_stage, 60, 70) ~ "flowering",
    between(zadok_stage, 70, 91) ~ "grain_fill")) %>%
  filter(sowing_das == 1) %>% 
  # Sort by trial, cultivar, dap
  arrange(trial, cultivar, day) %>%
  group_by(trial, cultivar) %>%
  nest() %>%
  mutate(data = map(data, ~mutate(., dap = seq(nrow(.))))) %>%
  unnest() %>%
  select(trial, environment, cultivar, date, day, dap, stage)


# Calculate avereage DAP of heading date predictions from CGM
cgm_predicted_heading <- apsim_s2met_cultivar_out1 %>% 
  filter(stage %in% c("heading", "flowering")) %>% 
  group_by(trial, environment, cultivar, stage) %>% 
  summarize(pred_HD = mean(dap)) %>%
  # summarize(pred_HD = min(dap)) %>%
  ungroup() %>%
  # Remove S2C1 trials
  filter(str_detect(trial, "S2C1", negate = TRUE))

# Get the environmental means of heading date from the AMMI model
env_mean_heading <- ammi_fit %>% 
  filter(trait == "HeadingDate") %>% 
  mutate(env_mean = map(fit_ammi, ~.x$mu + .x$Eeffect),
         env_mean = map(env_mean, ~tibble(environment = names(.), 
                                          obs_HD = .))) %>% 
  unnest(env_mean) %>%
  filter(environment != "EON17")


## Combine
pred_obs_heading <- inner_join(env_mean_heading, cgm_predicted_heading)

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
  mutate(cor = paste0("r==", round(cor, 3)),
         rmse = paste0("RMSE==", round(rmse, 3)))

## Plot each cultivar
(g_cgm_cultivar_summary <- pred_obs_heading %>%
  filter(stage == "flowering") %>%
  ggplot(aes(y = obs_HD, x = pred_HD)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ cultivar) +
  scale_y_continuous(name = "Observed mean heading (dap)") +
  scale_x_continuous(name = "CGM predicted flowering (dap)")  +
  geom_text(data = cgm_pred_HD_summary_annotate, 
            aes(x = 89, y = 54, label = cor), parse = T, hjust = 1, size = 2) +
  geom_text(data = cgm_pred_HD_summary_annotate, 
            aes(x = 89, y = 52, label = rmse), parse = T, hjust = 1, size = 2) +
  theme_presentation2(10))

ggsave(filename = "cgm_pred_obs_HD_cultivars.jpg", plot = g_cgm_cultivar_summary, path = fig_dir,
       height = 8, width = 10, dpi = 500)



############################
### Assess accuracy of heading date predictions from the chosen crop model
############################


# Calculate avereage DAP of heading date predictions from CGM
cgm_predicted_heading <- growth_stage_weather %>% 
  filter(stage %in% c("heading", "flowering")) %>% 
  group_by(trial, environment, stage) %>% 
  rename(pred_HD = dap) %>%
  summarize_at(vars(pred_HD, day), list(mean = ~mean, min = ~min)) %>%
  ungroup() %>%
  # Remove S2C1 trials
  filter(str_detect(trial, "S2C1", negate = TRUE))

# Get the environmental means of heading date from the AMMI model
env_mean_heading <- ammi_fit %>% 
  filter(trait == "HeadingDate") %>% 
  mutate(env_mean = map(fit_ammi, ~.x$mu + .x$Eeffect),
         env_mean = map(env_mean, ~tibble(environment = names(.), 
                                          obs_HD = .))) %>% 
  unnest(env_mean) %>%
  filter(environment != "EON17")


## Combine
pred_obs_heading <- inner_join(env_mean_heading, cgm_predicted_heading)

## Summarize
(cgm_pred_HD_summary <- pred_obs_heading %>%
  group_by(stage) %>%
  summarize_at(vars(pred_HD_mean, pred_HD_min), list(
    cor = ~cor(obs_HD, .), mae = ~mean(abs(. - obs_HD)), rmse = ~sqrt(mean((. - obs_HD)^2))
  )))


## Plot by stage
pred_obs_heading %>%
  ggplot(aes(x = pred_HD_mean, y = obs_HD)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ stage)



## Plot
(g_hd <- qplot(x = obs_HD, y = pred_HD_mean, data = pred_obs_heading, 
               geom = "point", facets = "stage", group = environment) +
    geom_abline(slope = 1, intercept = 0))
plotly::ggplotly(g_hd)

## Replot and save
cgm_pred_HD_summary_annotate <- cgm_pred_HD_summary %>%
  mutate(pred_HD_mean_cor = paste0("r==", round(pred_HD_mean_cor, 3)),
         pred_HD_mean_rmse = paste0("RMSE==", round(pred_HD_mean_rmse, 3)))
g_cgm_summary <- pred_obs_heading %>%
  filter(stage == "flowering") %>%
  ggplot(aes(x = obs_HD, y = pred_HD_mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  scale_x_continuous(name = "Observed mean heading (dap)") +
  scale_y_continuous(name = "CGM predicted flowering (dap)") +
  geom_text(data = filter(cgm_pred_HD_summary_annotate, stage == "flowering"), 
            aes(x = 51, y = 72, label = pred_HD_mean_cor), parse = T, hjust = 0) +
  geom_text(data = filter(cgm_pred_HD_summary_annotate, stage == "flowering"), 
            aes(x = 51, y = 70, label = pred_HD_mean_rmse), parse = T, hjust = 0)

ggsave(filename = "cgm_pred_obs_HD.jpg", plot = g_cgm_summary, path = fig_dir,
       height = 5, width = 5, dpi = 500)
  
  




############################
### Prepare environment-concurrent covariates
############################

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))


## Create the ECs and select the relevant ones for modeling
growth_stage_covariates1 <- growth_stage_covariates %>%
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

soil_covariates1 <- soil_covariates %>%
  ## Only use TP environments
  inner_join(., env_trials, by = "trial") %>%
  select(-trial)
  

# Add soil covariates
ec_select <- full_join(growth_stage_covariates1, soil_covariates1, by = "environment")


## Summarize min/max/var for each covariate
ec_select_summ <- ec_select %>%
  gather(covariate, value, -environment) %>%
  group_by(covariate) %>%
  summarize_at(vars(value), list(min = min, max = max, sd = sd, cv = cv), 
               na.rm = T)

## Remove covariates with low CV
ec_select <- select(ec_select, environment, subset(ec_select_summ, !is.na(cv), covariate, drop = TRUE))


##### 
# Look at correlations between covariates
##### 


## Calculate pairwise correlations
ec_pairwise_cor <- cor(ec_select[,-1], use = "pairwise.complete.obs") %>% 
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



## Test for normality using ks test
ec_tomodel_normality <- ec_select %>%
  gather(covariate, value, -environment) %>%
  group_by(covariate) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value, na.rm = T), sd = sd(.$value, na.rm = T))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value"),
         p_adj = p.adjust(p = p_value, method = "bonf"))

subset(ec_tomodel_normality, p_adj < 0.1)

## Which ECs should be kept
ec_to_keep <- subset(ec_tomodel_normality, p_adj >= 0.10, covariate, drop = TRUE)



# Remove some covariate and impute with the mean
ec_select1 <- ec_select %>%
  select(c("environment", ec_to_keep)) %>%
  # Inpute using the mean
  mutate_at(vars(-environment), impute)


## Prepare the covariates for modeling
## Center, but do not scale. Save the mean for later
ec_tomodel_temp <- ec_select1 %>%
  mutate_at(vars(-environment), scale, scale = FALSE, center = TRUE)

ec_tomodel_centers <- ec_tomodel_temp %>%
  summarize_at(vars(-environment), ~attr(., "scaled:center")) %>%
  gather(covariate, center)

## Convert scaled to numeric
ec_tomodel_centered <- ec_tomodel_temp %>%
  mutate_at(vars(-environment), as.numeric)

## Center and scale
ec_tomodel_scaled <- ec_select1 %>%
  mutate_at(vars(-environment), scale, scale = TRUE, center = TRUE) %>%
  mutate_at(vars(-environment), as.numeric)


# Vector of covariates
environmental_covariates <- names(ec_tomodel_centered)[-1]


## Visualize histogram
ec_select %>%
  gather(covariate, value, -environment) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ covariate, scales = "free_x") +
  theme_presentation2(10)



############################
### Prepare historical covariates
############################

# Load the covariate data
load(file.path(enviro_dir, "EnvironmentalCovariates/historical_environmental_covariates.RData"))

## Environments to use with corresponding trials
loc_trials <- S2_MET_BLUEs %>% 
  filter(environment %in% tp_vp_env,
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
  ## Only use TP environments
  inner_join(., loc_trials, by = "location") %>% 
  # Filter for years before those observed in the S2MET project
  filter(year < min(S2_MET_BLUEs$year)) %>%
  gather(covariate, value, -trial, -stage, -location, -year) %>%
  unite(covariate, stage, covariate, sep = ".") %>%
  ## Remove "radn_mean" 
  filter(str_detect(covariate, "radn_mean", negate = TRUE))

soil_covariates1 <- soil_covariates %>%
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca"),
         ## Replace Columbus with Wooster
         location = ifelse(location == "Columbus", "Wooster", location)) %>%
  ## Only use TP environments
  inner_join(., loc_trials, by = "location") %>%
  gather(covariate, value, -trial, -location) %>%
  left_join(distinct(growth_stage_covariates1, year, trial, location), .)

## Combine
ec_select <- bind_rows(growth_stage_covariates1, soil_covariates1)


## Summarize min/max/var for each covariate
ec_select_summ <- ec_select %>%
  group_by(covariate, year) %>%
  summarize_at(vars(value), list(min = min, max = max, sd = sd, cv = cv), 
               na.rm = T) %>%
  ungroup()

ec_year_tokeep <- ec_select_summ %>% 
  filter(!is.na(cv)) %>%
  select(covariate, year)

## Remove covariates with low CV
ec_select <- ec_select %>%
  inner_join(ec_year_tokeep, .)



# Vector of covariate names
ec_names <- unique(ec_select$covariate)


## Create timeframes ##
#
# Create two sets of time intervals:
# 1. Increasing number of years
# 2. Same number of years (5, 10, and 15) in a sliding window 
# 

# Create the vector of years
all_historical_years <- seq(max(ec_select$year) - 30 + 1, max(ec_select$year))

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
  mutate(ec_data = map(years, ~subset(ec_select, year %in% .x)))

# Summarize
ec_select_timeframe_summary <- ec_select_timeframe %>%
  unnest(ec_data) %>%
  group_by(time_frame, covariate, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

## Use sliding windows to summarize covariates
ec_select_window <- tibble(
  time_frame = unlist(imap(year_sliding_window, function(l, name) map_chr(l, ~paste0(name, "_", min(.x), "_", max(.x))))),
  years = unlist(year_sliding_window, recursive = FALSE)) %>%
  # Subset the historical EC
  mutate(ec_data = map(years, ~subset(ec_select, year %in% .x)))

# Summarize
ec_select_window_summary <- ec_select_window %>%
  unnest(ec_data) %>%
  group_by(time_frame, covariate, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()




## Measure the relationship between covariates from different time frames
ec_select_timeframe_corr <- ec_select_timeframe_summary %>% 
  group_by(covariate) %>%
  do(cor_mat = {
    df <- .
    select(df, -covariate) %>% 
      spread(time_frame, value) %>% 
      as.data.frame() %>% 
      column_to_rownames("location") %>% 
      cor()
  })

setNames(object = ec_select_timeframe_corr$cor_mat, ec_select_timeframe_corr$covariate)

## Spread the data
ec_select_timeframe_summary_wide <- spread(ec_select_timeframe_summary, covariate, value)
ec_select_window_summary_wide <- spread(ec_select_window_summary, covariate, value)


##### 
# Look at correlations between covariates
##### 


## Calculate pairwise correlations
ec_pairwise_cor <- ec_select_timeframe_summary_wide %>% 
  split(.$time_frame) %>% 
  map(~cor(.[,-1:-2]) %>%
        as.dist() %>%
        tidy() %>%
        rename(covariate1 = item1, covariate2 = item2, correlation = distance) %>%
        # Sort by descending absolution cor coef
        arrange(desc(abs(correlation))) )

## Are there covariates that tend to be highly correlated with others?
ec_pairwise_cor %>% 
  map(~group_by(., covariate1) %>% 
        summarize(correlation = mean(correlation)) %>% 
        arrange(desc(abs(correlation))) )

## Which covariates are enriched in correlations beyond some threshold?
cor_threshold <- 0.6
ec_pairwise_cor %>% 
  map(~filter(., correlation >= cor_threshold) %>%
        group_by(covariate2) %>%
        summarize(n = n()) %>%
        arrange(desc(n)) )


# Vector of covariates
environmental_covariates <- ec_names


## Test for normality using ks test
ec_tomodel_normality_timeframe <- ec_select_timeframe_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  group_by(time_frame, covariate) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value), sd = sd(.$value))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value")) %>%
  split(.$time_frame) %>%
  map_df(~mutate(., p_adj = p.adjust(p_value, method = "bonf")))

(ec_sig_nonnormal_timeframe <- subset(ec_tomodel_normality_timeframe, p_adj < 0.10))

# Check covariate by time_frame
# Determine which covariates to exclude
ec_to_remove_timeframe <- group_by(ec_sig_nonnormal_timeframe, covariate) %>% 
  summarize(n = n_distinct(time_frame)) %>% 
  filter(n > 20) %>% 
  pull(covariate)


# Repeat with window
ec_tomodel_normality_window <- ec_select_window_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  group_by(time_frame, covariate) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value), sd = sd(.$value))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value")) %>%
  split(.$time_frame) %>%
  map_df(~mutate(., p_adj = p.adjust(p_value, method = "bonf")))

(ec_sig_nonnormal_window <- subset(ec_tomodel_normality_window, p_adj < 0.10))

# Check covariate by time_frame
# Determine which covariates to exclude
ec_to_remove_window <- group_by(ec_sig_nonnormal_window, covariate) %>% 
  summarize(n = n_distinct(time_frame)) %>% 
  filter(n > 20) %>% 
  pull(covariate)



## Visualize
ec_select_timeframe_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  filter(covariate %in% unique(ec_sig_nonnormal_timeframe$covariate)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(time_frame ~ covariate)

## Remove this covariate from the df
ec_select_timeframe_summary_wide1 <- ec_select_timeframe_summary_wide %>%
  select(which(! names(.) %in% ec_to_remove_timeframe) ) %>%
  group_by(time_frame) %>%
  # Inpute using the mean
  mutate_at(vars(-location, -time_frame), impute) %>%
  ungroup()

## Visualize
ec_select_window_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  filter(covariate %in% unique(ec_sig_nonnormal_window$covariate)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(time_frame ~ covariate)

## Remove this covariate from the df
ec_select_window_summary_wide1 <- ec_select_window_summary_wide %>%
  select(which(! names(.) %in% ec_to_remove_window) ) %>%
  group_by(time_frame) %>%
  # Inpute using the mean
  mutate_at(vars(-location, -time_frame), impute) %>%
  ungroup()



## Prepare the covariates for modeling

## Time frame ##
## 
## Center, but do not scale. Save the mean for later
historical_ec_tomodel_timeframe_temp <- ec_select_timeframe_summary_wide1 %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = FALSE, center = TRUE))

historical_ec_tomodel_timeframe_centers <- historical_ec_tomodel_timeframe_temp %>%
  map(~summarize_at(.x, vars(-time_frame, -location), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
historical_ec_tomodel_timeframe_centered <- historical_ec_tomodel_timeframe_temp %>%
  map(~mutate_at(.x, vars(-time_frame, -location), as.numeric) )

## Center and scale
historical_ec_tomodel_timeframe_scaled <- ec_select_timeframe_summary_wide1 %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = TRUE, center = TRUE) %>%
        mutate_at(vars(-time_frame, -location), as.numeric) )


## Window ##
## 
## Center, but do not scale. Save the mean for later
historical_ec_tomodel_window_temp <- ec_select_window_summary_wide1 %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = FALSE, center = TRUE))

historical_ec_tomodel_window_centers <- historical_ec_tomodel_window_temp %>%
  map(~summarize_at(.x, vars(-time_frame, -location), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
historical_ec_tomodel_window_centered <- historical_ec_tomodel_window_temp %>%
  map(~mutate_at(.x, vars(-time_frame, -location), as.numeric) )

## Center and scale
historical_ec_tomodel_window_scaled <- ec_select_window_summary_wide1 %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = TRUE, center = TRUE) %>%
        mutate_at(vars(-time_frame, -location), as.numeric) )




## Save these results
save("ec_tomodel_centered", "ec_tomodel_scaled", "ec_tomodel_centers", 
     "historical_ec_tomodel_timeframe_centered", "historical_ec_tomodel_timeframe_scaled", "historical_ec_tomodel_timeframe_centers",
     "historical_ec_tomodel_window_centered", "historical_ec_tomodel_window_scaled", "historical_ec_tomodel_window_centers",
     file = file.path(result_dir, "concurrent_historical_covariables.RData"))

