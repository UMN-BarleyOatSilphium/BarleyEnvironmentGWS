## S2MET genotype-environment ecophysiological analysis
## 
## Determine the location covariance matrix by using the long-term
## average of environmental covariables that explain environment
## main effect and GxE
## 
## Author: Jeff Neyhart
## Last modified: 3 March 2020
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
loc_trials <- S2_MET_BLUEs %>% 
  filter(environment %in% tp_vp_env,
         str_detect(trial, "S2C1", negate = TRUE)) %>%
  distinct(trial, location)



############################
### Load data
############################

## Load the location blues
load(file.path(data_dir, "S2MET_Location_BLUEs.RData"))

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/historical_environmental_covariates.RData"))

## Load EC model building
load(file.path(result_dir, "ec_model_building.RData"))
ec_model_building <- unified_ec_models

## Unnest the ec model
ec_model_touse <- ec_model_building %>% 
  unnest(final_model) %>%
  filter(model == "model3_ammi")


############################
### Prepare ECs for modeling
############################


## Create the ECs and select the relevant ones for modeling
ec_select <- growth_stage_covariates %>%
  ## Remove heading as a growth stage
  filter(stage != "heading") %>%
  ## Separate trial into location and year
  mutate(year = as.numeric(str_extract(trial, "[0-9]{4}")),
         location = str_remove(trial, "_[0-9]{4}")) %>%
  select(-trial) %>%
  ## Only use TP environments
  inner_join(., loc_trials, by = "location") %>% 
  # Filter for years before those observed in the S2MET project
  filter(year < min(S2_MET_BLUEs$year)) %>%
  gather(covariate, value, -trial, -stage, -location, -year) %>%
  unite(covariate, stage, covariate, sep = ".") %>%
  ## Remove trange, relhum, and rain (water stress will cover this)
  filter(str_detect(covariate, "radn_mean", negate = TRUE))

# Vector of covariate names
ec_names <- unique(ec_select$covariate)


## Summarize covariates at a location over three time_frames
year_time_frame <- c(5, 15, 30)

ec_select_timeframe <- tibble(
  time_frame = year_time_frame, 
  end_year = max(ec_select$year), 
  start_year = end_year - time_frame + 1) %>%
  mutate(ec_data = map2(start_year, end_year, ~subset(ec_select, between(year, .x, .y))),
         time_frame = paste0("time_frame", time_frame))
         
  
## Summarize
ec_select_timeframe_summary <- ec_select_timeframe %>%
  unnest() %>%
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


## Prepare the covariates for modeling
## Center, but do not scale. Save the mean for later
historical_ec_tomodel_temp <- ec_select_timeframe_summary_wide %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = FALSE, center = TRUE))

historical_ec_tomodel_centers <- historical_ec_tomodel_temp %>%
  map(~summarize_at(.x, vars(-time_frame, -location), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
historical_ec_tomodel_centered <- historical_ec_tomodel_temp %>%
  map(~mutate_at(.x, vars(-time_frame, -location), as.numeric) )

## Center and scale
historical_ec_tomodel_scaled <- ec_select_timeframe_summary_wide %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = TRUE, center = TRUE) %>%
        mutate_at(vars(-time_frame, -location), as.numeric) )

  

# Vector of covariates
environmental_covariates <- ec_names




############################
### Use covariables determined to be significant in the environment-specific model
############################

## Calculate environmental similarity matrices

## Iterate over traits
location_relmat_df <- ec_model_touse %>%
  group_by(trait) %>%
  do({
    row <- .
    
    # Get the model formula
    ec_model_form <- formula(row$object[[1]])
    
    # main_environment_covariates
    main_environment_covariates <- all.vars(ec_model_form) %>% 
      subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
                                           pattern = .))))
    
    # interaction_environment_covariates
    interaction_environment_covariates <- all.vars(ec_model_form) %>% 
      subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
                                           pattern = .)))) %>% setdiff(., "line_name")
    
    ## Subset the ec_tomodel_centered list for these covariates - convert each to a matrix
    main_historical_covariates <- historical_ec_tomodel_scaled %>%
      map(~select(., location, main_environment_covariates) %>%
            as.data.frame() %>%
            column_to_rownames("location") %>%
            as.matrix() )
    
    interaction_historical_covariates <- historical_ec_tomodel_scaled %>%
      map(~select(., location, interaction_environment_covariates) %>%
            as.data.frame() %>%
            column_to_rownames("location") %>%
            as.matrix() )
    
    ## Create relationship matrices
    E_mat_main <- map(main_historical_covariates, ~Env_mat(x = ., method = "Rincent2019"))
    E_mat_int <- map(interaction_historical_covariates, ~Env_mat(x = ., method = "Rincent2019"))
    
    
    ## Output a tibble with the matrices and relationship matrices
    tibble(time_frame = names(main_historical_covariates), main_covariate_mat = main_historical_covariates, 
           interaction_covariate_mat = interaction_historical_covariates, E_mat_main, E_mat_int)
    
    
  }) %>% ungroup()




## Save these results
save("location_relmat_df", "historical_ec_tomodel_centered", "historical_ec_tomodel_centers", 
     "historical_ec_tomodel_scaled", "ec_select_timeframe",
     file = file.path(result_dir, "historical_ec_model_building.RData"))

