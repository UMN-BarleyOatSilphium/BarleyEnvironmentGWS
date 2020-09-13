## S2MET Genomewide Environment Predictions
## 
## Run APSIM crop model using the raw weather data
## 
## 

# Repository directory
repo_dir <- getwd()
# Source the main project script
source(file.path(repo_dir, "source.R"))

# Other packages
library(lubridate)

# Path for the apsim.exe
apsim_exe <- "C:/Program Files (x86)/Apsim79-r4044/Model/Apsim.exe"

## Filter the trials for those without missing lat/lon coordinates
trial_meta1 <- trial_info %>%
  filter_at(vars(latitude, longitude), any_vars(!is.na(.))) %>%
  # Convert any dates to ymd
  mutate_at(vars(contains("date")), ymd)


# Read in irrigation information
trial_irrigation_data <- map_df(list.files(meta_dir, recursive = TRUE, pattern = "irrigation", full.names = TRUE), read_csv) %>%
  mutate(Date = ymd(Date), yday = yday(Date)) %>%
  rename_all(tolower)


# Load weather and soil data
load(file.path(result_dir, "raw_weather_soil_data.RData"))



# Concurrent crop growth models ------------------------------------------------------

## Adjust precipitation using irrigation information
concurrent_daily_weather_daymet1 <- concurrent_daily_weather_daymet %>%
  left_join(., nest(group_by(trial_irrigation_data, trial), .key = "irrigation"), by = "trial") %>%
  # Edit the weather data
  mutate(data = map2(data, irrigation, ~{
    # If .y is null, pass
    if (is.null(.y)) {
      .x
    } else {
      # Left join
      left_join(.x, .y, by = "yday") %>%
        # edit prcp
        mutate(irrigation_mm = ifelse(is.na(irrigation_mm), 0, irrigation_mm),
               prcp = prcp + irrigation_mm) %>%
        # remove excess columns
        select(-date:-irrigation_mm)
    } })) %>%
  # Remove the irrigation column
  select(-irrigation)

# Add date to data
concurrent_daily_weather_daymet2 <- concurrent_daily_weather_daymet1 %>%
  mutate(data = map(data, ~mutate(., date = today(), date = `year<-`(date, year), date = `yday<-`(date, yday))),
         # Edit columns
         data = map(data, ~rename(., maxt = tmax, mint = tmin, rain = prcp)))



# Create MET files --------------------------------------------------------

# Create a temporary directory to store the files
met_dir <- file.path(data_dir, "METfiles")
dir.create(path = met_dir)

# Use the create_MET function to output MET files
create_MET(trial.info = concurrent_daily_weather_daymet2, data.col = "data", dir = met_dir)







# Run APSIM crop growth model with multiple cultivars ---------------------

# This will help determine the representative cultivar for future models

# File path of the base apsim simulation file
apsim_base <- file.path(data_dir, "barley_base_crop_simulation_factorial_cultivars.apsim")

# Create a temporary directory to store the output files
apsim_dir <- file.path(data_dir, "APSIMfiles")
dir.create(apsim_dir)

# Copy the apsim base file to this directory; change the path of the base apsim file
file.copy(from = apsim_base, to = file.path(apsim_dir, basename(apsim_base)))
apsim_base_use <- file.path(apsim_dir, basename(apsim_base))

# Run for each trial
concurrent_growth_staging_cultivars_daymet <- concurrent_daily_weather_daymet2 %>%
  filter(trial %in% subset(trial_meta1, project2 == "S2MET" | project3 == "S2MET", trial, drop = T)) %>%
  # filter(trial == "S2_MET_CRM16") %>%
  run_apsim(trial.info = ., met.dir = met_dir, apsim.base = apsim_base_use, apsim.exe = apsim_exe) %>%
  ## Add the evapotranspiration to the daily "data"
  mutate(apsim_out = map(apsim_out, "barley1"),
         apsim_out = map(apsim_out, as_tibble), 
         data = map2(data, apsim_out, ~left_join(.x, select(.y, date = Date, pet = eo), by = "date")))



# Assess accuracy of heading date predictions from multiple cultivar  --------


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


## Reorganize output
concurrent_growth_staging_cultivars_daymet1 <- concurrent_growth_staging_cultivars_daymet %>% 
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
cgm_predicted_heading <- concurrent_growth_staging_cultivars_daymet1 %>% 
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

## Save these for plotting later
save("cgm_pred_HD_summary", "cgm_pred_HD_summary", file = file.path(result_dir, "apsim_cultivar_predicted_flowering_results.RData"))

## Find the cultivar that results in the most accuracy prediction
(repr_cultivar <- cgm_pred_HD_summary %>%
  filter(stage == "flowering") %>%
  arrange(rmse, desc(cor)) %>%
  as.data.frame() %>%
  slice(1))

# bass





# Run APSIM crop growth model ------------------------------------------

# File path of the base apsim simulation file
apsim_base <- file.path(data_dir, "barley_base_crop_simulation.apsim")

# Create a temporary directory to store the output files
apsim_dir <- file.path(data_dir, "APSIMfiles")
dir.create(apsim_dir)

# Copy the apsim base file to this directory; change the path of the base apsim file
file.copy(from = apsim_base, to = file.path(apsim_dir, basename(apsim_base)), overwrite = TRUE)
apsim_base_use <- file.path(apsim_dir, basename(apsim_base))

# Edit the apsim base file to switch the cultivar #
# (this is done using the gui)

# Run for each trial
concurrent_growth_staging_daymet <- concurrent_daily_weather_daymet2 %>%
  filter(trial %in% subset(trial_meta1, project2 == "S2MET" | project3 == "S2MET", trial, drop = T)) %>%
  # filter(trial == "S2_MET_CRM16") %>%
  run_apsim(trial.info = ., met.dir = met_dir, apsim.base = apsim_base_use, apsim.exe = apsim_exe) %>%
  ## Add the evapotranspiration to the daily "data"
  mutate(apsim_out = map(apsim_out, "barley1"),
         apsim_out = map(apsim_out, as_tibble), 
         data = map2(data, apsim_out, ~left_join(.x, select(.y, date = Date, pet = eo), by = "date")))


# Remove data from the metfiles directory
unlink(x = list.files(met_dir, full.names = TRUE))




# Historical location growth models ---------------------------------------


# Extract location metadata df
location_meta <- location_historical_daily_weather_daymet %>%
  select(location, latitude, longitude)



# Determine the planting date for each year in each location --------------

# Read in the crop calendar data from http://nelson.wisc.edu/sage/data-and-models/crop-calendar-dataset/index.php
crop_calendar <- read_csv("http://nelson.wisc.edu/sage/data-and-models/crop-calendar-dataset/All_data_with_climate.csv")

## Filter for barley and the U.S.
crop_calendar_barley <- crop_calendar %>%
  filter(Crop == "Barley", Crop.name.in.original.data == "Barley - spring") %>%
  # Select location and relevant climate variables related to planting
  select(Location, contains("plant"))


library(sp)
library(spdep)
library(maptools)
library(maps)
library(rgeos)

## For each location and lat/long, determine the closest state
# Load state data
states <- map("state", fill=TRUE, col="transparent", plot = FALSE)
IDs <- sapply(X = strsplit(x = states$names, split = ":"), "[[", 1)
states_sp <- map2SpatialPolygons(states, IDs = str_to_title(IDs), proj4string = CRS("+proj=longlat +datum=WGS84"))
# subset states in crop calendar
states_sp_use <- states_sp[crop_calendar_barley$Location[-1]]
state_names <- sapply(states_sp_use@polygons, function(x) x@ID)

# Determine closest state
location_closest_state <- location_meta %>%
  distinct(location, latitude, longitude) %>%
  mutate(closest_state = map2_chr(latitude, longitude, ~{
    pointsSP <- SpatialPoints(data.frame(x = .y, y = .x), proj4string = CRS("+proj=longlat +datum=WGS84"))
    state_names[which.min(gDistance(pointsSP, states_sp_use, byid = TRUE)[,1])]
  }))

## Unload packages
detach("package:maptools")
detach("package:spdep")
detach("package:sp")
detach("package:maps")



## Unnest the historical data and group by year
## Add the closest state to the location historical data
location_historical_daily_weather2 <- left_join(location_historical_daily_weather_daymet, location_closest_state) %>%
  unnest(daily_data) %>%
  mutate(date = ymd(paste0(year, "0101")),
         date = `yday<-`(date, yday)) %>%
  rename(mint = tmin, maxt = tmax) %>%
  group_by(location, latitude, longitude, closest_state, year) %>%
  nest() %>%
  mutate(planting_date = as.character(NA))

## Loop over each row
## Determine the planting date in each year for each location
for (i in seq(nrow(location_historical_daily_weather2))) {
  
  ## First pull the state to use with the crop calendar
  crop_calendar_state <- subset(crop_calendar_barley, Location == location_historical_daily_weather2$closest_state[i])
  
  ## Pull out relevant climate thresholds
  tmean <- floor(crop_calendar_state$temp.at.planting[1])
  daylen <- floor(crop_calendar_state$daylength.at.planting[1]) - 1
  
  # Pull out the historical data and determine the planting date
  year_planting_dates <- location_historical_daily_weather2$data[[i]] %>%
    # Calculate daily average temp
    mutate(meant = (mint + maxt) / 2) %>%
    # Filter days for the threshold temperature and daylength
    # Also find a day with little rain
    filter(meant >= tmean - 2,
           daylength >= daylen - 0.5,
           mint > 0, 
           prcp <= 0.5) %>%
    # Find the first day
    subset(., date == min(date), date, drop = TRUE)
  
  # Convert planting data to character
  pd <- as.character(year_planting_dates)
  
  # Return this string; if empty, return NA
  # Add this to the historical data
  location_historical_daily_weather2$planting_date[i] <- ifelse(is_empty(pd), NA, pd)
  
}

location_historical_daily_weather_daymet1 <- location_historical_daily_weather2



# Run apsim for each year of historical location data ---------------------

# File path of the base apsim simulation file
apsim_base <- file.path(data_dir, "barley_base_crop_simulation.apsim")

# Create a temporary directory to store the output files
apsim_dir <- file.path(data_dir, "APSIMfiles")
dir.create(apsim_dir)

# Copy the apsim base file to this directory; change the path of the base apsim file
file.copy(from = apsim_base, to = file.path(apsim_dir, basename(apsim_base)), overwrite = TRUE)
apsim_base_use <- file.path(apsim_dir, basename(apsim_base))


## Add a dummy trial name
## Add date to data
location_historical_daily_weather_daymet2 <- location_historical_daily_weather_daymet1 %>%
  ## Add Dec31 days for leap years (daymet only returns 365 values per year, regardless of year)
  mutate(data = modify_at(data, which(leap_year(year)), ~bind_rows(.x, slice(.x, nrow(.x))) %>% 
                            mutate(yday = seq_len(nrow(.)), date = `yday<-`(date, yday)))) %>%
  mutate(data = map(data, ~rename(., rain = prcp)),
         trial = paste0(location, "_", year),
         apsim_out = list(NULL)) %>%
  arrange(location, year) %>%
  # Filter out trials with missing temperature data
  filter(map_lgl(data, ~mean(is.na(.$mint)) == 0))


## Loop over each year
for (i in seq_len(nrow(location_historical_daily_weather_daymet2))) {
  
  ## Subset the df
  df_to_use <- location_historical_daily_weather_daymet2[i,]
  # Pull out the trial name
  trial_i <- df_to_use$trial[1]
  
  ## Create a met file
  # Use the create_MET function to output MET files
  create_MET(trial.info = df_to_use, data.col = "data", dir = met_dir)
  
  ## Run apsim
  apsim_out <- run_apsim(trial.info = df_to_use, met.dir = met_dir, apsim.base = apsim_base_use, apsim.exe = apsim_exe)
  # Add the results to the df
  location_historical_daily_weather_daymet2$apsim_out[[i]] <- apsim_out$apsim_out[[1]]
  
}

## Rename and save
location_historical_growth_staging_daymet <- location_historical_daily_weather_daymet2 %>%
  select(-trial) %>%
  ## Add evapotranspiration to the "data" df
  mutate(data = map2(data, apsim_out, ~mutate(.x, pet = .y$eo)),
         apsim_out = map(apsim_out, as_tibble) %>% map(1)) %>%
  select(location, latitude, longitude, predicted_planting_date = planting_date,
         daily_data = data, growth_model = apsim_out)



# Remove data from the metfile and apsim file temporary directories; delete the directories
unlink(x = met_dir, recursive = TRUE)
unlink(x = apsim_dir, recursive = TRUE)


# Save the data
save("concurrent_growth_staging_daymet", "location_historical_growth_staging_daymet", "concurrent_growth_staging_cultivars_daymet",
     file = file.path(result_dir, "apsim_growth_model_results.RData"))


