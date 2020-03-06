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
library(modelr)
library(broom)
library(pbr)
library(cowplot)
library(pls)
library(car)

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

## Load the fitted ammi models
load(file.path(result_dir, "ammi_model_fit.RData"))

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
  ## Remove "radn_mean" 
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
ec_tomodel_normality <- ec_select_timeframe_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  group_by(time_frame, covariate) %>%
  do(ks_test = ks.test(x = .$value, y = "pnorm", mean = mean(.$value), sd = sd(.$value))) %>%
  ungroup() %>%
  mutate(p_value = map_dbl(ks_test, "p.value"))

(ec_sig_nonnormal <- subset(ec_tomodel_normality, p_value < 0.10))

## Visualize
ec_select_timeframe_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  filter(covariate %in% unique(ec_sig_nonnormal$covariate)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(time_frame ~ covariate)

## Remove this covariate from the df
ec_select_timeframe_summary_wide1 <- ec_select_timeframe_summary_wide %>%
  select(which(! names(.) %in% unique(ec_sig_nonnormal$covariate)) )


## Determine if there is sufficient variation for a covariate
ec_select_timeframe_summary_wide1 %>%
  split(.$time_frame) %>%
  map_df(~summarize_at(., vars(-location, -time_frame), var)) %>%
  as.data.frame()



## Prepare the covariates for modeling
## Center, but do not scale. Save the mean for later
historical_ec_tomodel_temp <- ec_select_timeframe_summary_wide1 %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = FALSE, center = TRUE))

historical_ec_tomodel_centers <- historical_ec_tomodel_temp %>%
  map(~summarize_at(.x, vars(-time_frame, -location), ~attr(., "scaled:center")) %>%
        gather(covariate, center) )

## Convert scaled to numeric
historical_ec_tomodel_centered <- historical_ec_tomodel_temp %>%
  map(~mutate_at(.x, vars(-time_frame, -location), as.numeric) )

## Center and scale
historical_ec_tomodel_scaled <- ec_select_timeframe_summary_wide1 %>%
  split(.$time_frame) %>%
  map(~mutate_at(.x, vars(-time_frame, -location), scale, scale = TRUE, center = TRUE) %>%
        mutate_at(vars(-time_frame, -location), as.numeric) )



## Final model preparation ##


## Prepare data for modelling:
## 1. only use the tp
## 2. add ECs to the df

s2_met_tomodel <- historical_ec_tomodel_centered %>%
  map_df(~left_join(x = filter(S2_MET_BLUEs, line_name %in% tp), y = ., by = "location")) %>%
  # map(~left_join(x = filter(S2_MET_Loc_BLUEs, line_name %in% tp), y = ., by = "location")) %>%
  mutate_at(vars(line_name, environment, location), as.factor)

#  New vector of covariates
environmental_covariates <- s2_met_tomodel %>% 
  select(matches("vegetative|heading|flowering|grain_fill")) %>%
  names()


## Model all covariates without apriori information
historical_ec_by_trait <- tibble(
  trait = sort(traits),
  covariates = list(
    environmental_covariates, # Grain protein
    environmental_covariates, # Grain yield
    str_subset(environmental_covariates, "vegetative"), # Heading date
    str_subset(environmental_covariates, "grain", negate = TRUE), # Plant height
    environmental_covariates # Test weight
  )
)

# New vector of covariates
environmental_covariates_all <- reduce(historical_ec_by_trait$covariates, union)




############################
### Model covariates that are predictive of the location mean
############################

## Get the location effects from the AMMI model (average environmental effects)
base_model_location_effect <- ammiN_fit_location %>%
  unnest(loc_scores) %>% 
  select(trait, location, effect)


## Use stepwise regression to find the best model
location_effect_to_model <- base_model_location_effect %>%
  left_join(., s2_met_tomodel)


# Group by trait
location_effect_models <- location_effect_to_model %>%
  group_by(trait, time_frame) %>%
  nest() %>%
  mutate(model = list(NULL))



for (i in seq(nrow(location_effect_models))) {
  
  df <- location_effect_models$data[[i]]
  tr <- location_effect_models$trait[i]
  
  ## Remove some covariates depending on when the trait is measured
  df1 <- select(df, location, effect, subset(historical_ec_by_trait, trait == tr, covariates, drop = T)[[1]]) %>%
    distinct()
  
  # Minimal model
  min_model <- effect ~ 1
  # Max model
  max_model <- add_predictors(min_model, as.formula(paste0("~", paste0(names(df1)[-1:-2], collapse = " + "))))
  
  # Stepwise regression
  fit_base <- lm(formula = min_model, data = df1, contrasts = "contr.sum")
  fit_step <- step(object = fit_base, scope = max_model, direction = "forward")
  
  
  # Did we run out of dfs? If so, reduce terms
  while (df.residual(fit_step) <= 1) {
    
    ## Drop1
    fit_step_drop1 <- drop1(fit_step) %>%
      as.data.frame() %>%
      rownames_to_column("term")
    
    ## Remove the two covariates that gives the lowest AIC, besides "none"
    to_remove <- fit_step_drop1 %>%
      filter(str_detect(term, "none", negate = T)) %>% 
      top_n(x = ., n = 2, wt = -AIC) %>%
      pull(term) %>% 
      sort() %>%
      first()
    
    ## Reduce the number of parameters by 1
    fit_step_new_formula <- terms(formula(fit_step)) %>% 
      drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% to_remove), keep.response = TRUE) %>%
      formula()
    
    ## Refit the model
    fit_step <- update(object = fit_step, formula = fit_step_new_formula)
    
  }
  
  fit_step1 <- fit_step
  
  
  
  
  ## Remove covariates based on VIF
  # Proceed only if the number of terms is >= 2
  if (length(attr(terms(formula(fit_step1)), "term.labels")) >= 2) {
  
    vif_i <- vif(fit_step1)
    fit_stepi <- fit_step1
    vif_cutoff <- 4
    
    while (any(vif_i > vif_cutoff) & class(vif_i) != "try-error") {
      
      form_i <- terms(formula(fit_stepi)) %>%
        drop.terms(., dropx = which( attr(., "term.labels") %in% names(which.max(vif_i))), keep.response = TRUE ) %>%
        formula()
      
      fit_stepi <- update(object = fit_stepi, formula = form_i)
      vif_i <- try( vif(fit_stepi), silent = TRUE)
    }
    
    # Backward elimination from here
    fit_step_final <- step(object = fit_stepi, direction = "backward")
    
  } else {
    fit_step_final <- fit_step1
    
  }
  
  ## Return the model
  location_effect_models$model[[i]] <- fit_step_final
  
}


## Examine coefficients of covariates
location_effect_models$model %>% 
  map(summary)




## Plot fitted values versus observed environmental mean
location_effect_models %>%
  mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>%
  unnest(model_data) %>%
  ggplot(aes(x = pred, y = effect)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_wrap(~ trait + time_frame, scales = "free", ncol = 3) +
  theme_acs()



############################
### Identify covariates that are correlated with the AMMI location distance matrix
############################


## Using the AMMI output and phi, calculate the distance matrix between environments
## This will be denoted as W_ammi, consistent with Rincent 2019
ammi_dist <- ammiN_fit_location %>%
  # Calculate a location PHI matrix
  mutate(loc_phi = map2(.x = g_scores, .y = loc_scores, ~outer(X = .x$score, Y = .y$score) %>%
                          `dimnames<-`(x = ., value = list(.x$line_name, .y$location)) )) %>%
  mutate(W = map(loc_phi, ~as.matrix(dist(t(.)))),
         W_ammi = map(W, ~1 - (. / max(.)))) %>%
  select(trait, W_ammi)


# Create a matrix of scaled and centered covariates
historical_ec_tomodel_scaled_mat_list <- historical_ec_tomodel_scaled %>% 
  map(~as.data.frame(.) %>%
        column_to_rownames("location") %>% 
        select(-time_frame) %>%
        as.matrix() )

## Add the coviarates assigned to each trait
historical_ec_tomodel_ammi <- left_join(ammi_dist, historical_ec_by_trait, by = "trait") %>%
  crossing(., tibble(time_frame = names(historical_ec_tomodel_scaled_mat_list), 
                     ec_mat = historical_ec_tomodel_scaled_mat_list))


# Tolerance of difference of correlations
cor_tol <- 0.01


## For each trait, use an algorithm to determine the covariates to use
## This is based on the approach of Rincent 2019
historical_ec_ammi_dist <- historical_ec_tomodel_ammi %>%
  group_by(trait, time_frame) %>%
  do({
    
    row <- .
    
    # AMMI dist matrix
    W_ammi <- row$W_ammi[[1]]
    # Vector of covariates
    covariates <- row$covariates[[1]]
    
    ## Subset the EC scaled matrix
    ec_tomodel_scaled_mat_use <- row$ec_mat[[1]][rownames(W_ammi), covariates]
    
    ## Choose the starting covariate as the one whose distance matrix has
    ## the highest correlation with the W_ammi distance matrix
    initial_W_cov <- map(covariates, ~make_dist_mat(ec_tomodel_scaled_mat_use[, ., drop = FALSE])) %>%
      map_dbl(~cor(as.numeric(.), as.numeric(W_ammi)))
    
    starting_covariate <- covariates[which.max(initial_W_cov)]
    # Create two vectors of available and used covariates
    used_covariates <- tested_covariates <- starting_covariate
    available_covariates <- setdiff(covariates, used_covariates)
    
    old_cor <- max(initial_W_cov)
    
    ## Vector of correlations
    i <- 1
    ec_ammi_cor_list <- vector("list", length(covariates))
    ec_ammi_cor_list[[i]] <- tibble(added_covariate = starting_covariate, cor_with_ammi = max(initial_W_cov))
    
    ## For each remaining covariate, find the next one that gives the highest gain (or lowest loss)
    while (length(available_covariates) > 0) {
      
      # Advance i
      i <- i + 1
      
      # Scan the available covariates, add to the matrix, calculate distance, compare to W_ammi
      new_W_cov <- map(available_covariates, ~c(used_covariates, .)) %>%
        map(~make_dist_mat(ec_tomodel_scaled_mat_use[, ., drop = FALSE])) %>%
        map_dbl(~cor(as.numeric(.), as.numeric(W_ammi)))
      
      old_cor <- max(new_W_cov)
      
      # Find the max difference
      which_max_diff <- which.max(new_W_cov - old_cor)
      
      # Get the covariate
      selected_covariate <- available_covariates[which_max_diff]
      
      ## Add to the list
      ec_ammi_cor_list[[i]] <- tibble(added_covariate = selected_covariate, cor_with_ammi = new_W_cov[which_max_diff])
      
      # Determine remaining covariates
      used_covariates <- union(used_covariates, selected_covariate)
      # Determine new available covariate vector
      available_covariates <- setdiff(covariates, used_covariates)
      
    }
    
    ## Collapse the list; calculate difference from previous
    ec_ammi_cor_df <- bind_rows(ec_ammi_cor_list) %>%
      mutate(difference = c(diff(cor_with_ammi), 0),
             stop = which.min(difference >= cor_tol),
             number = seq(nrow(.)))
    
    
    ## List of final covariates
    final_covariates <- subset(ec_ammi_cor_df, number <= stop, added_covariate, drop = TRUE)
    # Final distance matrix
    ec_dist_mat_final <- make_dist_mat(ec_tomodel_scaled_mat_use[, final_covariates, drop = FALSE])
    # final correlation
    final_cor <- subset(ec_ammi_cor_df, number == stop, cor_with_ammi, drop = TRUE)
    
    ## Return these results
    tibble(final_cor = final_cor, final_covariates = list(final_covariates), ec_dist_mat_final = list(ec_dist_mat_final),
           test_results = list(ec_ammi_cor_df))
    
  }) %>% ungroup()





### Create relationship matrices for locations ####
location_relmat_df <- historical_ec_ammi_dist %>%
  # Add location main effect models
  left_join(., select(location_effect_models, -data)) %>%
  # Get a list of main effect models
  mutate(main_environment_covariates = map(model, ~attr(terms(formula(.)), "term.labels"))) %>%
  mutate_at(vars(main_environment_covariates, final_covariates), ~map2(.x = ., .y = time_frame, ~{
    historical_ec_tomodel_scaled_mat_list[[.y]][, .x, drop = FALSE]
  })) %>%
  # Calculate E mat for the main environmental covariates
  mutate(E_mat_main = map(main_environment_covariates, ~{
    if (ncol(.x) == 0) {
      d1 <- diag(x = nrow(.x))
      `dimnames<-`(x = d1, value = replicate(2, row.names(.x), simplify = FALSE))
    } else {
      Env_mat(x = .x, method = "Rincent2")
    } }) ) %>%
  select(trait, time_frame, interaction_covariate_mat = final_covariates,
         main_covariate_mat = main_environment_covariates, E_mat_main, E_mat_int = ec_dist_mat_final)






############################
### Use covariables determined to be significant in the environment-specific model
############################

## Calculate environmental similarity matrices

## Iterate over traits
location_relmat_previous_ecs_df <- ec_model_touse %>%
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
      map(~select(., which(names(.) %in% c("location", main_environment_covariates))) %>%
            as.data.frame() %>%
            column_to_rownames("location") %>%
            as.matrix() )
    
    interaction_historical_covariates <- historical_ec_tomodel_scaled %>%
      map(~select(., which(names(.) %in% c("location", interaction_environment_covariates))) %>%
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
save("location_relmat_df", "location_relmat_previous_ecs_df", "historical_ec_tomodel_centered", 
     "historical_ec_tomodel_centers", "historical_ec_tomodel_scaled", "ec_select_timeframe",
     file = file.path(result_dir, "historical_ec_model_building.RData"))

