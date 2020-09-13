## S2MET Genomewide Environment Predictions
## 
## Query data sources for external environmental data
## 
## 


# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))

# Extra libraries
library(lubridate)
library(rvest)
library(daymetr) # Also for meterological data
library(geosphere) # For GCD calculation and daylength
library(RSQLite)
library(dbplyr)


## Filter the trials for those without missing lat/lon coordinates
trial_meta1 <- trial_info %>%
  filter_at(vars(latitude, longitude), any_vars(!is.na(.))) %>%
  # Convert any dates to ymd
  mutate_at(vars(contains("date")), ymd)


# Change the enviro_dir
envir_dir <- "Z:/BARLEY_LAB/Climate Data/"


# Number of years of historical data
historical_length <- 30

## Parameters to query (others will be included by default)
pars <- c("tmin" = "TS_MIN", "tmax" = "TS_MAX", "prcp" = "PRECTOT", "relhum" = "RH2M")



# Daylength ---------------------------------------------------------------


## Iterate over trials
daily_daylength_data <- trial_meta1 %>%
  mutate(data = list(NULL))


## Loop over trials
for (i in seq(nrow(daily_daylength_data))) {
  
  ## Vector of dates
  dates <- seq(ymd(paste0(year(daily_daylength_data$planting_date[i]), "0101")), 
               ymd(paste0(year(daily_daylength_data$planting_date[i]), "1231")), by = "day")
  
  # query daylength data
  daylen <- daylength(lat = daily_daylength_data$latitude[i], doy = yday(dates))
  
  ## Return a df
  daily_daylength_data$data[[i]] <- data.frame(date = dates, daylength = daylen, stringsAsFactors = FALSE)
}



# Concurrent weather data -------------------------------------------------


## Use a function that iterates over trials
daily_weather_data <- get_weather_data(trial.info = trial_meta1, site = "daymet")


## Remove some variables and save

# Daymet does not count leap years (all years are 365 days)
# For those years, add a dummy dec. 31 by coping data from dec. 30

## Combine and save
concurrent_daily_weather_daymet <- daily_weather_data %>%
  mutate(data = modify_at(data, which(leap_year(year)), ~bind_rows(.x, slice(.x, nrow(.x))) %>% mutate(yday = seq_len(nrow(.))))) %>%
  unnest(data) %>% 
  dplyr::select(-type:-project3, -t3_trial_name:-plot_dim) %>%
  group_by(trial, location, year, environment, planting_date, harvest_date, latitude, longitude, elevation) %>% 
  nest() %>%
  mutate(data = map(data, ~rename(., year = year1)))





# Historical daylength and weather data -----------------------------------

## Find distinct location/lat/long coordinates from the trial metadata
location_meta <- trial_meta1 %>%
  filter_at(vars(latitude, longitude), any_vars(!is.na(.))) %>%
  group_by(location) %>% 
  slice(1) %>% 
  distinct(location, latitude, longitude) %>%
  ungroup() %>%
  mutate(data = list(NULL))

# Two ending years:
# 1. The last year of data in the metadata
# 2. The year before the first year of data (this is where the 30-year history will begin)
end_years <- range(trial_meta1$year)
first_year <- end_years[1] - historical_length

# Vector of years
seq_years <- seq(first_year, end_years[2])



## Daylength

## Iterate over trials
daylength_data <- location_meta %>%
  mutate(daylength_data = list(NULL))


## Loop over trials
for (i in seq(nrow(daylength_data))) {
  
  ## Vector of dates
  dates <- seq(ymd(paste0(first_year, "0101")), ymd(paste0(end_years[2], "1231")), by = "day")
  
  # query daylength data
  daylen <- daylength(lat = daylength_data$latitude[i], doy = yday(dates))
  
  ## Return a df
  daylength_data$daylength_data[[i]] <- tibble(date = dates, daylength = daylen)
}


## Weather data
historical_weather_data_daymet <- location_meta


## Loop over locations and query data
for (i in seq(nrow(historical_weather_data_daymet))) {
  
  meta_i <- historical_weather_data_daymet[i,]
  
  ##
  ## Create dates using the end year and historical timeframe
  dates <- paste0(c(first_year, end_years[2]), c("-01-01", "-12-31"))
  
  # Get lat/long
  lonlat <- unlist(meta_i[c("longitude", "latitude")])
  
  ## Pull data
  data_out <- download_daymet(lat = lonlat["latitude"], lon = lonlat["longitude"], 
                              start = min(year(dates)), end = max(year(dates)), silent = TRUE)
  
  ## Get relevant data and rename
  data_out1 <- data_out$data
  names(data_out1) <- c("year", "yday", "daylength_sec", "prcp", "srad", "swe", "tmax", "tmin", "vp")
  
  # Convert units
  # srad to daily total radiation
  data_out1$radn <- (data_out1$srad * data_out1$daylength_sec) / 1000000
  # Daylength to hours
  data_out1$daylength <- data_out1$daylength_sec / 3600
  
  # Remove some variables
  data_out2 <- data_out1[,-which(names(data_out1) %in% c("daylength_sec", "srad"))]
  
  # Extract elevation
  elevation <- data_out$altitude
  
  ## Output data
  historical_weather_data_daymet$data[[i]] <- list(elevation = elevation, data = data_out2)
  
  ## Notify user
  cat("\nData collected for location", meta_i$location, "going back", historical_length, "years.")
  
}


## Ungroup
historical_weather_data_daymet1 <- historical_weather_data_daymet %>%
  mutate(daily_weather_data = map(data, "data") %>% map(as_tibble),
         elevation = map_dbl(data, "elevation")) %>%
  select(-data)

## Merge daylength with daily weather data
location_historical_daily_weather_daymet <- historical_weather_data_daymet1 %>%
  select(location, latitude, longitude, elevation, daily_data = daily_weather_data)





# Soil data ---------------------------------------------------------------


## Location of HWSD .bil file
hwsd_bil <- file.path(envir_dir, "RawData/SoilData/HWSD/SpatialFiles/hwsd.bil")

## Get the soil data
soil_data <- get_soil_data(trial.info = trial_meta1, hwsd.bil = hwsd_bil)



# Save data ---------------------------------------------------------------

save("concurrent_daily_weather_daymet", "location_historical_daily_weather_daymet","soil_data", 
     file = file.path(result_dir, "raw_weather_soil_data.RData"))

