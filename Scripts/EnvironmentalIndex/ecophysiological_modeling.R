## S2MET ecophysiology
## 
## See if heading and maturity can be accurately predicted using GDD information
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

