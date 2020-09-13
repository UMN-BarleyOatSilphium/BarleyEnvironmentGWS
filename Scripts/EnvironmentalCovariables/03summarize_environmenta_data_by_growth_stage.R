## S2MET Genomewide Environment Predictions
## 
## Summarize environmental covariates
## 
##  

# Repository directory
repo_dir <- getwd()
# Source the main project script
source(file.path(repo_dir, "source.R"))

# Other packages
library(lubridate)


## Filter the trials for those without missing lat/lon coordinates
trial_meta1 <- trial_info %>%
  filter_at(vars(latitude, longitude), any_vars(!is.na(.))) %>%
  # Convert any dates to ymd
  mutate_at(vars(contains("date")), ymd)


# mint, maxt, rh2m are summarized by mean
# rain radn, and water_stress are summarized by sum
mean_vars <- c("maxt", "mint", "tmean", "radn")
sum_vars <- c("radn",  "water_balance", "gdd")


# Load data ---------------------------------------------------------------

# Load crop model output
load(file.path(result_dir, "apsim_growth_model_results.RData"))
# Load raw environmental data
load(file.path(result_dir, "raw_weather_soil_data.RData"))


# Add source
concurrent_growth_staging_daymet <- mutate(concurrent_growth_staging_daymet, source = "daymet")
location_historical_growth_staging_daymet <- mutate(location_historical_growth_staging_daymet, source = "daymet")

# Use the results of the APSIM models to determine growth stages
# Z10 - Z30: early vegetative
# Z30 - Z49: late vegetative
# Z50 - Z59: heading
# Z60 - Z69: flowering
# Z70 - Z91: grain fill



# Concurrent data ---------------------------------------------------------

concurrent_trial_growth_stages <- concurrent_growth_staging_daymet %>%
  mutate(growth_stages = map(apsim_out, ~{
    
    df <- .
    
    ## Determine stages
    stages <- mutate(df, stage = case_when(
      between(zadok_stage, 10, 30) ~ "early_vegetative",
      between(zadok_stage, 30, 50) ~ "late_vegetative",
      # between(zadok_stage, 50, 60) ~ "heading",
      between(zadok_stage, 50, 70) ~ "flowering",
      between(zadok_stage, 70, 91) ~ "grain_fill"))
    
    
    ## Return a df with date, dap, growth stage
    stages %>% 
      filter(sowing_das == 1) %>% 
      mutate(dap = seq(nrow(.))) %>% 
      select(date = Date, day, dap, stage)
    
  }))


## Add daily weather observations during for each day during a growth stage
concurrent_growth_stage_weather <- concurrent_trial_growth_stages %>%
  unnest(growth_stages) %>%
  left_join(., unnest(concurrent_trial_growth_stages, data)) %>%
  # Add GDD info
  left_join(., unnest(concurrent_trial_growth_stages, apsim_out) %>% select(source, trial, date = Date, tt = TT)) %>%
  ## Calculate water stress as the difference between pet and rain
  mutate(water_balance = rain - pet) %>%
  rename(gdd = tt)


## Summarize covariates in each growth stage
concurrent_growth_stage_covariates <- concurrent_growth_stage_weather %>%
  filter(!is.na(stage)) %>%
  # Calculate diurnal range and mean
  mutate(tmean = (maxt + mint) / 2,
         trange = maxt - mint) %>%
  select(source, trial, stage, mean_vars, sum_vars) %>%
  group_by(source, trial, stage) %>%
  ## Summarize mean covariates
  { full_join(
    x = summarize_at(., vars(mean_vars), list(mean = ~mean(.)), na.rm = TRUE),
    y = summarize_at(., vars(sum_vars), list(sum = ~sum(.)), na.rm = TRUE)
  ) } %>%
  ungroup()




# Historical location data ------------------------------------------------

# Add trial column 
location_historical_growth_staging_daymet <- location_historical_growth_staging_daymet %>%
  mutate(year = year(ymd(predicted_planting_date)),
         trial = paste(location, year, sep = "_"),
         year = year(ymd(predicted_planting_date)))



# Determine discrete growth stages based on zadok stage

historical_location_predicted_growth_stages <- location_historical_growth_staging_daymet %>%
  mutate(growth_stages = map(growth_model, ~{
    
    df <- .
    
    ## Determine stages
    stages <- mutate(df, stage = case_when(
      between(zadok_stage, 10, 30) ~ "early_vegetative",
      between(zadok_stage, 30, 50) ~ "late_vegetative",
      # between(zadok_stage, 50, 60) ~ "heading",
      between(zadok_stage, 50, 70) ~ "flowering",
      between(zadok_stage, 70, 91) ~ "grain_fill"))
    
    
    ## Return a df with date, dap, growth stage
    stages %>% 
      filter(sowing_das == 1) %>% 
      mutate(dap = seq(nrow(.))) %>% 
      select(date = Date, day, dap, stage)
    
  }))


## Add daily weather observations during for each day during a growth stage
historical_location_growth_stage_weather <- historical_location_predicted_growth_stages %>%
  unnest(growth_stages) %>%
  left_join(., unnest(historical_location_predicted_growth_stages, daily_data)) %>%
  # Add GDD info
  left_join(., unnest(historical_location_predicted_growth_stages, growth_model) %>% select(source, trial, date = Date, tt = TT, pet = eo)) %>%
  ## Calculate water stress as the difference between pet and rain
  mutate(water_balance = rain - pet) %>%
  rename(gdd = tt)


## Summarize covariates in each growth stage
historical_location_growth_stage_covariates <- historical_location_growth_stage_weather %>%
  filter(!is.na(stage)) %>%
  # Calculate diurnal range and mean
  mutate(tmean = (maxt + mint) / 2,
         trange = maxt - mint) %>%
  select(source, trial, location, stage, mean_vars, sum_vars) %>%
  group_by(source, trial, location, stage) %>%
  ## Summarize mean covariates
  { full_join(
    x = summarize_at(., vars(mean_vars), list(mean = ~mean(.)), na.rm = TRUE),
    y = summarize_at(., vars(sum_vars), list(sum = ~sum(.)), na.rm = TRUE)
  ) } %>%
  ungroup()




# Soil covariates ------------------------------------


## Calculated the weighted mean of numeric variables; these are weighted by the percent
## sharing of a soil type
soil_covariates <- soil_data %>%
  select(trial, location, soil_data) %>% 
  unnest() %>%
  filter(trial %in% concurrent_growth_stage_covariates$trial) %>%
  mutate_if(is.numeric, as.numeric) %>%
  select(trial, location, share_in_mapping_unit, matches("^topsoil|^subsoil")) %>%
  mutate(share_in_mapping_unit = share_in_mapping_unit / 100) %>%
  group_by(trial, location) %>%
  select_if(is.numeric) %>% # Remove character columns
  summarize_at(vars(-share_in_mapping_unit), ~weighted.mean(x = ., w = share_in_mapping_unit)) %>%
  select(-contains("texture")) %>%
  ungroup()



## Save
save("concurrent_growth_stage_covariates", "concurrent_growth_stage_weather", 
     "historical_location_growth_stage_covariates", "historical_location_growth_stage_weather",
     "soil_covariates",
     file = file.path(result_dir, "environmental_covariates.RData"))

