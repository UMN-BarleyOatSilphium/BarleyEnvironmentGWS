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
library(nasapower)

# Significance level
alpha <- 0.05

## Environments to use with corresponding trials
env_trials <- S2_MET_BLUEs %>% 
  distinct(trial, environment) %>%
  filter(environment %in% tp_vp_env,
         str_detect(trial, "S2C1", negate = TRUE))

## Load the location blues
load(file.path(data_dir, "S2MET_Location_BLUEs.RData"))



# ############################
# ### Query historical weather data
# ############################
# 
# # Load functions for querying environmental data or running apsim
# source(file.path(enviro_dir, "enviroPipeline/environmental_data_functions.R"))
# 
# ## Other directories in the EnvironmentalData directories
# growth_stage_dir <- file.path(enviro_dir, "GrowthStaging")
# met_dir <- file.path(growth_stage_dir, "MetFiles")
# apsim_dir <- file.path(growth_stage_dir, "APSIMFiles")
# ec_dir <- file.path(enviro_dir, "EnvironmentalCovariates")
# 
# ## Timeline of historical data to use (in years)
# historical_length <- c(5, 10, 30)
# # Ending year to collect data
# end_year <- min(S2_MET_BLUEs$year) - 1
# 
# ## List of location metadata 
# location_meta <- S2_MET_Loc_BLUEs %>%
#   distinct(location) %>%
#   left_join(., distinct(trial_info, location, latitude, longitude) %>% group_by(location) %>% slice(1) %>% ungroup()) %>%
#   crossing(., time_frame = historical_length) %>%
#   mutate(out = list(NULL))
# 
# 
# 
# 
# 
# 
# ## Parameters to query (others will be included by default)
# pars <- c("tmin" = "TS_MIN", "tmax" = "TS_MAX", "prcp" = "PRECTOT", "relhum" = "RH2M")
# ## Add the needed pars to the parameters list
# # Rename for .met file
# needed_pars <- c("mint" = "TS_MIN", "maxt" = "TS_MAX", "radn" = "ALLSKY_SFC_SW_DWN", "rain" = "PRECTOT")
# # Add names for the parameters (default is to lower)
# pars1 <- setNames(pars, tolower(pars))
# pars2 <- union(needed_pars, pars1)
# names(pars2)[needed_pars %in% pars2] <- names(needed_pars)[needed_pars %in% pars2]
# names(pars2)[pars2 %in% setdiff(pars2, needed_pars)] <- names(pars1)[pars1 %in% setdiff(pars2, needed_pars)]
# 
# # Separate renaming vector (with date)
# pars2_rename <- c("date" = "YYYYMMDD", pars2)
# 
# ## Loop over locations and query data
# for (i in seq(nrow(location_meta))) {
# 
#   meta_i <- location_meta[i,]
#   back_years <- meta_i$time_frame
#   
#   ## 
#   ## Create dates using the end year and historical timeframe
#   dates <- paste0(c((end_year - back_years) + 1, end_year), c("-01-01", "-12-31"))
#   
#   # Get lat/long
#   lonlat <- unlist(meta_i[c("longitude", "latitude")])
#   
#   ## Pull data
#   data_out <- get_power(community = "AG", pars = pars2, temporal_average = "DAILY", lonlat = lonlat, dates = dates)
#   
#   ## Get relevant data and rename
#   data_out1 <- as.data.frame(data_out)[pars2_rename]
#   names(data_out1) <- names(pars2_rename)
#   
#   # Remove missing values
#   data_out1[data_out1 == -99] <- NA
#   
#   
#   # Extract elevation
#   elevation <- as.numeric(str_remove(string = str_extract(string = attr(data_out, "POWER.Elevation"), pattern = "[0-9]*.[0-9]* meters"), 
#                                      pattern = " meters"))
#   
#   ## Output data
#   location_meta$out[[i]] <- list(elevation = elevation, data = data_out1)
#   
#   ## Notify user
#   cat("\nData collected for location", meta_i$location, "going back", meta_i$time_frame, "years.")
#   
# }
# 
# ## Ungroup
# location_historical_daily_weather <- location_meta %>%
#   mutate(data = map(out, "data"),
#          data = map(data, as_tibble)) %>%
#   select(-out)
# 
# 
# ## Calculate daylength
# daily_daylength_data <- location_meta %>%
#   mutate(out = list(NULL))
# 
# 
# ## Loop over trials
# for (i in seq(nrow(daily_daylength_data))) {
#   
#   ## Vector of dates
#   start_end_dates <- ymd(paste0(c((end_year - daily_daylength_data$time_frame[i]) + 1, end_year), c("0101", "1231")))
#   
#   dates <- seq(from = start_end_dates[1], to = start_end_dates[2], by = "day")
#   
#   # query daylength data
#   daylen <- geosphere::daylength(lat = daily_daylength_data$latitude[i], doy = yday(dates))
#   
#   ## Return a df
#   daily_daylength_data$out[[i]] <- data.frame(date = dates, daylength = daylen, stringsAsFactors = FALSE)
# }
# 
# 
# ## Merge with daily weather data
# location_historical_daily_weather1 <- location_historical_daily_weather %>%
#   left_join(., daily_daylength_data) %>%
#   mutate(data = map2(data, out, left_join, by = "date")) %>%
#   select(-out)
# 
# 
# 
# 
# ############################
# ### Determine the planting date for each year
# ############################
# 
# ## First read in crop calendar information
# crop_calendar <- read_csv(file = file.path(growth_stage_dir, "crop_calendar_climate.csv"))
# 
# ## Filter for barley and the U.S.
# crop_calendar_barley <- crop_calendar %>%
#   filter(Crop == "Barley", Crop.name.in.original.data == "Barley - spring") %>%
#   # Select location and relevant climate variables related to planting
#   select(Location, contains("plant"))
# 
# 
# library(sp)
# library(spdep)
# library(maptools)
# library(maps)
# library(rgeos)
# 
# ## For each location and lat/long, determine the closest state
# # Load state data
# states <- map("state", fill=TRUE, col="transparent", plot = FALSE)
# IDs <- sapply(X = strsplit(x = states$names, split = ":"), "[[", 1)
# states_sp <- map2SpatialPolygons(states, IDs = str_to_title(IDs), proj4string = CRS("+proj=longlat +datum=WGS84"))
# # subset states in crop calendar
# states_sp_use <- states_sp[crop_calendar_barley$Location[-1]]
# state_names <- sapply(states_sp_use@polygons, function(x) x@ID)
# 
# # Determine closest state
# location_closest_state <- location_meta %>%
#   distinct(location, latitude, longitude) %>%
#   mutate(closest_state = map2_chr(latitude, longitude, ~{
#     pointsSP <- SpatialPoints(data.frame(x = .y, y = .x), proj4string = CRS("+proj=longlat +datum=WGS84"))
#     state_names[which.min(gDistance(pointsSP, states_sp_use, byid = TRUE)[,1])]
#   }))
# 
# ## Unload packages
# detach("package:maptools")
# detach("package:spdep")
# detach("package:sp")
# detach("package:maps")
# 
# 
# ## Unnest the historical data and group by year
# ## Add the closest state to the location historical data
# location_historical_daily_weather2 <- left_join(location_historical_daily_weather1, location_closest_state) %>%
#   unnest(data) %>%
#   mutate(year = year(date)) %>%
#   filter(time_frame == max(time_frame)) %>%
#   group_by(location, latitude, longitude, time_frame, closest_state, year) %>%
#   nest() %>%
#   mutate(planting_date = as.character(NA))
# 
# ## Loop over each row
# ## Determine the planting date in each year for each location
# for (i in seq(nrow(location_historical_daily_weather2))) {
#   
#   ## First pull the state to use with the crop calendar
#   crop_calendar_state <- subset(crop_calendar_barley, Location == location_historical_daily_weather2$closest_state[i])
#   
#   ## Pull out relevant climate thresholds
#   tmean <- floor(crop_calendar_state$temp.at.planting[1])
#   daylen <- floor(crop_calendar_state$daylength.at.planting[1]) - 1
#   
#   # Pull out the historical data and determine the planting date
#   year_planting_dates <- location_historical_daily_weather2$data[[i]] %>%
#     # Filter days for the threshold temperature and daylength
#     # Also find a day with little rain
#     filter(((mint + maxt) / 2) >= tmean, daylength >= daylen, mint > 0, rain <= 0.5) %>%
#     # Find the first day
#     subset(., date == min(date), date, drop = TRUE)
#   
#   # Add this to the historical data
#   location_historical_daily_weather2$planting_date[i] <- as.character(year_planting_dates)
#   
# }
# 
# 
# 
# 
# 
# ############################
# ### Run APSIM for each year in each location
# ############################
# 
# ## Add a dummy trial name
# location_historical_daily_weather3 <- location_historical_daily_weather2 %>%
#   mutate(trial = paste0(location, "_", year)) %>%
#   mutate(apsim_out = list(NULL)) %>%
#   arrange(location, year)
# 
# 
# ## Loop over each year
# for (i in seq(nrow(location_historical_daily_weather3))) {
#   
#   ## Subset the df
#   df_to_use <- location_historical_daily_weather3[i,]
#   # Pull out the trial name
#   trial_i <- df_to_use$trial[1]
# 
#   ## Create a met file
#   # Use the create_MET function to output MET files
#   create_MET(trial.info = df_to_use, data.col = "data", dir = met_dir)
#   
#   ## Run apsim
#   apsim_out <- run_apsim(trial.info = df_to_use, met.dir = met_dir, apsim.base = apsim_base, apsim.exe = apsim_exe)
#   # Add the results to the df
#   location_historical_daily_weather3$apsim_out[[i]] <- apsim_out$apsim_out[[1]]
#   
#   # Remove the met file
#   invisible(file.remove(list.files(met_dir, pattern = trial_i, full.names = T)))
#   
# }
# 
# ## Rename and save
# location_historical_climate_data <- location_historical_daily_weather3 %>%
#   select(-trial) %>%
#   ## Add evapotranspiration to the "data" df 
#   mutate(data = map2(data, apsim_out, ~mutate(.x, pet = .y$eo)))
# 
# 
# 
# 
# ############################
# ### Summarize covariates
# ############################
# 
# 
# # Use the results of the APSIM models to determine three growth phases:
# # vegetative: emergence-flowering
# # flowering: flowering-start grain fill
# # grain fill: start-end grain fill
# # 
# 
# location_historical_growth_stages <- location_historical_climate_data %>%
#   mutate(growth_stages = imap(apsim_out, ~{
#     
#     df <- .x
#     
#     # print(.y)
#     
#     # Stage ranges
#     stage_ranges <- summarize_at(df, vars(emergence_das, flowering_das, start_grain_fill_das, end_grain_fill_das), ~unique(.[. > 0])) %>% unlist()
#     # List of days in each stage
#     stage_list <- vector("list", length = length(stage_ranges) - 1)
#     names(stage_list) <- c("vegetative", "flowering", "grain_fill")
#     
#     # Iterate
#     for (i in seq(2, length(stage_ranges))) stage_list[[i-1]] <- seq(stage_ranges[i-1], stage_ranges[i]-1)
#     
#     # Create a df
#     stage_df <- tibble(dap = unlist(stage_list), growth_stage = rep(names(stage_list), map_dbl(stage_list, length)))
#     
#     ## Return a df with date, dap, growth stage
#     df %>% 
#       filter(sowing_das == 1) %>%
#       mutate(dap = seq(nrow(.))) %>%
#       select(date, day, dap) %>% 
#       inner_join(., stage_df, by = "dap")
#     
#   }))
# 
# 
# ## Add daily weather observations during each day during a growth stage
# growth_stage_historical_weather <- location_historical_growth_stages %>%
#   unnest(growth_stages) %>%
#   left_join(., unnest(location_historical_growth_stages, data)) %>%
#   ## Calculate water stress as the difference between pet and rain
#   mutate(water_stress = pet - rain)
# 
# 
# ## Summarize covariates in each growth stage
# ## mint, maxt, rh2m get means
# ## rain radn, and water_stress gets sum
# mean_vars <- c("maxt", "mint", "tmean", "trange", "rh2m")
# sum_vars <- c("radn", "rain",  "water_stress")
# 
# growth_stage_historical_covariates <- growth_stage_historical_weather %>%
#   # Calculate diurnal range and mean
#   mutate(tmean = (maxt + mint) / 2,
#          trange = maxt - mint) %>%
#   select(location, year, growth_stage, mean_vars, sum_vars) %>%
#   gather(covariate, value, -location, -year, -growth_stage) %>%
#   group_by(location, year, growth_stage, covariate) %>%
#   summarize_at(vars(value), list(sum = ~sum(., na.rm = TRUE), mean = ~mean(., na.rm = T))) %>%
#   ungroup() %>%
#   mutate(statistic = ifelse(covariate %in% sum_vars, sum, mean)) %>%
#   select(-sum, -mean) %>%
#   spread(covariate, statistic)
# 
# 
# ## Save
# save("growth_stage_historical_covariates", "growth_stage_historical_weather", "location_historical_climate_data", 
#      file = file.path(ec_dir, "s2met_historical_environmental_covariates.RData"))



############################
### Load data
############################


## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_historical_environmental_covariates.RData"))


############################
### Visualization of historical weather data
############################

## Select an example location
historical_weather_select <- growth_stage_historical_weather %>%
  filter(location == "St_Paul") %>%
  select(location, year, dap, growth_stage, mint:water_stress)


#### Example daily temperature for the location

## color vector for growth stage
growth_stage_color <- setNames(neyhart_palette("umn2")[c(4, 5, 2)], names(growth_stage_rename))
growth_stage_color["vegetative"] <- neyhart_palette()[5]

# what covariate to plot
ec_to_plot <- "maxt"

daily_ec_example_plots <- historical_weather_select %>%
  gather(covariate, value, -location, -year, -dap, -growth_stage) %>%
  # Re-order growth stages
  mutate(growth_stage = factor(growth_stage, levels = c("vegetative", "flowering", "grain_fill"))) %>%
  filter(covariate == ec_to_plot) %>%
  mutate(min_dap = min(dap), max_dap = max(dap)) %>%
  ## calculate ranges for each covariate
  group_by(covariate) %>%
  mutate_at(vars(value), list(~min, ~max)) %>%
  group_by(covariate, location) %>%
  do(plot = {
    df <- .
    
    ## Create segments for growth stages
    gs_seg <- df %>% 
      group_by(growth_stage) %>% 
      summarize(start = min(dap) - 1, end = max(dap))
    
    ec_name <- covariate_variable_rename[unique(df$covariate)]
    ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
    
    # x axis limits
    x_limit <- c(df$min_dap[1] - 1, df$max_dap[1])
    
    # y axis limits
    y_limit <- c(df$min[1], df$max[1])
    y_end <- quantile(y_limit, 0.75)
    
    ## Summarize the covariate over dap
    df1 <- df %>%
      group_by(dap) %>%
      summarize(mean_value = mean(value), lower = quantile(value, alpha / 2), upper = quantile(value, 1 - (alpha / 2)),
                min_dap = mean(min_dap), max_dap = mean(max_dap), min = mean(min), max = mean(max))
    
    ## Plot
    ggplot(df1, aes(x = dap, y = mean_value)) +
      geom_segment(data = gs_seg, mapping = aes(x = start, xend = end, y = y_end, yend = y_end, color = growth_stage), lwd = 10) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85", alpha = 0.5) +
      geom_line() +
      scale_y_continuous(name = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), breaks = pretty, limits = y_limit) +
      scale_x_continuous(breaks = pretty, name = "Days after planting", limits = x_limit) +
      scale_color_manual(values = growth_stage_color, name = "Predicted\ngrowth stage",
                         labels = function(x) str_to_title(str_remove_all(x, "_"))) +
      labs(subtitle = unique(df$location)) +
      theme_presentation2() +
      theme(legend.position = "bottom")
      
  })

  





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


