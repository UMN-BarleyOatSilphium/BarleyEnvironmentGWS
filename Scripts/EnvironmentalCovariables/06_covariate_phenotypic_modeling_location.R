## S2MET Genomewide Environment Predictions
## 
## Phenotypic modeling with covariates - genotype-location means
## 

# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))

library(modelr)
library(broom)
library(lme4)
library(car)

## significance level
alpha <- 0.05

# data source
source_use <- "daymet"

# Number of cores
n_cores <- 8



# Load data ---------------------------------------------------------------

# Load covariates for environments and historical covariates
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


## Filter BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp,
         environment %in% train_test_env) %>% # only model train/test environments
  mutate_at(vars(line_name, environment), as.factor)


## Ideas
## 1. Calculate locations means and use that as the input (and validation)
## 2. Use all data and fit a model with random year, location:year, and g:year
## 3. Use all data and fit the g + location model


# Modeling ----------------------------------------------------------------

# Historical covariates ===================================================

## Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

# Use just the TP for modeling
S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp)

# Vector of covariate names
covariate_names <- names(historical_ec_tomodel_timeframe_centered[[1]])[-1:-3]




# Identify the most suitable time frame for summarization -----------------


# Implement the stepwise regression algorithm per timeframe to identify the
# timeframe with the lowest RMSEP

# Combine the lists of covariates
historical_ec_tomodel_centered_use <- c(
  historical_ec_tomodel_timeframe_centered %>% subset(., grepl(pattern = source_use, x = names(.))),
  historical_ec_tomodel_window_centered %>% subset(., grepl(pattern = source_use, x = names(.)))
)

historical_ec_tomodel_scaled_use <- c(
  historical_ec_tomodel_timeframe_scaled %>% subset(., grepl(pattern = source_use, x = names(.))),
  historical_ec_tomodel_window_scaled %>% subset(., grepl(pattern = source_use, x = names(.)))
)

# Create a data.frame of covariates per trait
trait_covariate_df <- historical_ec_tomodel_centered_use %>%
  map(~names(.)[-1:-3]) %>%
  reduce(union) %>%
  crossing(trait = traits, covariate = .) %>%
  filter(
    covariate != "awc_range",
    !(str_detect(covariate, "bulk_density")),
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  ) 

# Merge covariates with the phenotypic data
S2_MET_loc_BLUEs_tomodel <- S2_MET_loc_BLUEs_tomodel %>%
  left_join(., bind_rows(historical_ec_tomodel_centered_use), by = "location")


# Create a setup df
historical_timeframe_selection <- S2_MET_loc_BLUEs_tomodel %>%
  group_by(trait, time_frame) %>% 
  group_data() %>% # This creates a pointer to the larger dataframe
  # Output list
  mutate(out = list(NULL)) %>%
  # Adhoc and adhoc no soil
  crossing(., feat_sel_type = c("adhoc", "adhoc_nosoil"))




## Parallelize
historical_timeframe_selection_out <- historical_timeframe_selection %>%
  # Break up by cores
  assign_cores(df = ., n_core = n_cores, split = TRUE) %>%
  # Apply a function over cores
  coreApply(X = ., FUN = function(core_df) {

# # Iterate over rows
# for (i in seq_len(nrow(historical_timeframe_selection))) {

    # Create an output vector
    output <- vector("list", length = nrow(core_df))
    
    # Seq along the output vector
    for (i in seq_along(output)) {
      
      # Get the data row indicator
      df_rows <- core_df$.rows[[i]]
      tr <- core_df$trait[i]
      fst <- core_df$feat_sel_type[i]
    
      # Get the data and factorize line name and location
      
      # Factorize
      df1 <- mutate_at(S2_MET_loc_BLUEs_tomodel[df_rows,,drop = FALSE], vars(line_name, location), ~fct_contr_sum(as.factor(.)))
      
      loo_indices <- df1 %>%
        group_by(location) %>%
        crossv_loo_grouped() %>%
        pull(train) %>%
        map("idx")
      
      # 1. Fit a base model
      base_fit <- lm(value ~ 1 + line_name, data = df1)
      covariates_use <- subset(trait_covariate_df, trait == tr, covariate, drop = TRUE)
      
      # Remove soil variables if called for
      if (fst == "adhoc_nosoil") {
        covariates_use <- str_subset(string = covariates_use, pattern = "soil", negate = TRUE)
      }
      
      
      ## Recursive feature addition ##
      ## Main effect

      # 2. Define the scope
      scope <- list(lower = formula(base_fit), upper = reformulate(c("line_name", covariates_use), response = "value"))
      temp <- rfa_loo(object = base_fit, data = df1, scope = scope, metric = "RMSE",
                      index = loo_indices, env.col = "location")
      
      # Run rfa
      rfa_out <- try(rfa_loo(object = base_fit, data = df1, scope = scope, metric = "RMSE",
                             index = loo_indices, env.col = "location"))

      # If try is error, set rfa_out and rfa_out_int both as null
      if (class(rfa_out) == "try-error") {
        rfa_out <- rfa_out_int <- NULL

        # else proceed
      } else {

        ## Interactions
        # 1. Fit a base model
        base_fit_int <- update(base_fit, formula = reformulate(rfa_out$optVariables, response = "value"))
        # 2. Define the scope
        scope <- list(lower = formula(base_fit_int),
                      upper = reformulate(c(rfa_out$optVariables, paste0("line_name:", covariates_use)), response = "value"))
        # Run rfa
        rfa_out_int <- try(rfa_loo(object = base_fit_int, data = df1, scope = scope, metric = "RMSE",
                                   index = loo_indices, env.col = "location"))

        # Set rfa_out_int to null if error
        if (class(rfa_out_int) == "try-error") rfa_out_int <- NULL

      }

      ## Create a tibble and store in output
      output[[i]] <- tibble(
        model = c("model4", "model5"),
        covariates = list(rfa_out, rfa_out_int),
      )

      # Notify
      cat("Stepwise selection for trait", tr, "with timeframe", core_df$time_frame[i], "complete.\n")
      
    } # Close the loop
    
    # Add output to core_df; return
    core_df %>%
      mutate(out = output) %>%
      select(-.rows, -core) %>%
      unnest(out)
    
  }) %>% bind_rows() # Close the coreApply function


# Save these results
save("historical_timeframe_selection_out", file = file.path(result_dir, "historical_covariate_timeframe_selection.RData"))


## Stop the script if not interactive
if (!interactive()) stop("End of script.")










# Identify the best timeframes using the LASSO ----------------------------

# Calculate environmental means
loc_means <- S2_MET_loc_BLUEs_tomodel %>% 
  distinct(trait, line_name, location, value) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Factorize
    df1 <- df %>%
      droplevels() %>%
      mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.)))
    
    # Fit the model
    fit <- lmer(value ~ (1|line_name) + location, data = df1)
    
    ## Return a df of environmental effects
    fixef(fit) %>% 
      tibble(location = names(.), effect = .) %>% 
      filter(location != "(Intercept)") %>% 
      mutate(location = str_remove(location, "location"),
             mu = fixef(fit)[1]) %>% 
      add_row(location = last(levels(df1$location)), effect = -sum(.$effect), 
              mu = .$mu[1])
    
  }) %>% ungroup()



# Merge with the covariates
S2_MET_loc_means_tomodel <- loc_means %>%
  rename(value = effect) %>%
  left_join(., bind_rows(historical_ec_tomodel_scaled_use), by = "location")
  

# Create a setup df
historical_timeframe_selection_lasso <- S2_MET_loc_means_tomodel %>%
  group_by(trait, time_frame) %>% 
  group_data() %>% # This creates a pointer to the larger dataframe
  # Output list
  mutate(out = list(NULL))

# Iterate over rows
for (i in seq_len(nrow(historical_timeframe_selection_lasso))) {

  row <- historical_timeframe_selection_lasso[i,]

  # Factorize
  df1 <- S2_MET_loc_means_tomodel[row$.rows[[1]], ,drop = FALSE] %>%
    mutate_at(vars(location), ~fct_contr_sum(as.factor(.)))
  
  # Vector of covariates for this trait
  covariates_use <- subset(trait_covariate_df, trait == row$trait, covariate, drop = TRUE)

  ## Include interactions between rainfall and soil variables
  interaction_covariates <- cross(list(str_subset(covariates_use, "water_balance"), str_subset(covariates_use, "soil"))) %>%
    map_chr(~map_chr(., 1) %>% paste0(., collapse = ":"))

  covariates_use <- c(covariates_use, interaction_covariates)

  ## Run variable importance using LASSO
  lasso_out <- select_features_met(data = df1, covariates.use = covariates_use, env.col = "location", search.method = "lasso")

  ## Reformulate for export
  lasso_out1 <- lasso_out %>%
    mutate(adhoc = pmap(list(adhoc, adhoc_RMSE, adhoc_R2), 
                             ~list(optVariables = names(..1[..1 != 0,]), finalResults = c("RMSE" = ..2, "MSE" = ..2^2, R2 = ..3), 
                                   importance = ..1)),
           adhoc_nosoil = pmap(list(adhoc_nosoil, adhoc_nosoil_RMSE, adhoc_nosoil_R2), 
                                    ~list(optVariables = names(..1[..1 != 0,]), finalResults = c("RMSE" = ..2, "MSE" = ..2^2, R2 = ..3), 
                                          importance = ..1))) %>%
    select(model, adhoc, adhoc_nosoil) %>%
    gather(feat_sel_type, covariates, -model)


  ## Create a tibble
  historical_timeframe_selection_lasso$out[[i]] <- lasso_out1

  # Notify
  cat("LASSO analysis for trait", row$trait, "with timeframe", row$time_frame, "complete.\n")

}

## Tidy up
historical_timeframe_selection_lasso_out <- historical_timeframe_selection_lasso %>%
  unnest(out) %>%
  mutate(model = ifelse(model == "model2", "model4", "model5"))



# Analyze the results -----------------------------------------------------


load(file.path(result_dir, "historical_covariate_timeframe_selection.RData"))

# What timeframe returns the lowest RMSE per trait?
historical_timeframe_analysis <- bind_rows(mutate(historical_timeframe_selection_lasso_out, method = "lasso"),
                                           mutate(historical_timeframe_selection_out, method = "stepwise")) %>%
  filter(model == "model5") %>%
  mutate(results = map(covariates, "finalResults") %>% map(as.list) %>% map(as_tibble),
         covariates = map(covariates, "optVariables") %>% map(~setdiff(., "line_name")),
         nCovariates = map_dbl(covariates, length)) %>%
  filter(map_lgl(results, ~ncol(.) > 0)) %>%
  unnest(results) %>%
  # Parse timeframe
  mutate(time_frame_type = str_extract(time_frame, "time_frame|window"),
         time_frame1 = str_remove(time_frame, "time_frame|window")) %>%
  separate(time_frame1, c("length", "start_year", "end_year"), sep = "_") %>%
  mutate_at(vars(length, contains("year")), parse_guess) %>%
  # Round the RMSE
  mutate_at(vars(RMSE, R2), ~round(., 5)) %>%
  arrange(trait, RMSE, length, desc(end_year))


## For each trait and method (stepwise/lasso), choose the timeframe
## with the lowest overall RMSE and the most recent timeframe (end year within
## three years of the trial beginnings) with the lowest RMSE
best_historical_overall_timeframe <- historical_timeframe_analysis %>%
  group_by(trait, method, feat_sel_type) %>%
  top_n(x = ., n = 1, wt = -RMSE) %>%
  ungroup() %>%
  arrange(trait, method, feat_sel_type, desc(time_frame_type)) %>%
  group_by(trait, method, feat_sel_type) %>%
  slice(1) %>% ungroup()

best_historical_recent_timeframe <- historical_timeframe_analysis %>%
  # Filter for qualified timeframes
  filter(end_year %in% seq(min(trial_info$year) - 3, min(trial_info$year) - 1)) %>%
  group_by(trait, method, feat_sel_type) %>%
  top_n(x = ., n = 1, wt = -RMSE) %>%
  ungroup() %>%
  arrange(trait, method, feat_sel_type, desc(time_frame_type)) %>%
  group_by(trait, method, feat_sel_type) %>%
  slice(1) %>% ungroup()

# Combine
best_historical_timeframe <- bind_rows(
  mutate(best_historical_overall_timeframe, selection = "bestOverall"),
  mutate(best_historical_recent_timeframe, selection = "bestRecent")
)

# Create a separate feature selection data.frame using the stepwise covariates
historical_feature_selection <- historical_timeframe_selection_out %>%
  inner_join(., subset(best_historical_timeframe, method == "stepwise", c(trait, time_frame, feat_sel_type, selection))) %>%
  mutate(feat_sel_type = paste0("stepwise_cv_", feat_sel_type))

# Subset the LASSO results
historical_feature_importance <- historical_timeframe_selection_lasso_out %>%
  inner_join(., subset(best_historical_timeframe, method == "lasso", c(trait, time_frame, feat_sel_type, selection))) %>%
  mutate(feat_sel_type = paste0("lasso_cv_", feat_sel_type))



## Use the above timeframes for the apriori covariates


## A priori covariates - quite minimal; each should have a citation
##
## HeadingDate - everything before flowering
##
## PlantHeight - everything before grain fill
##
## Grain yield, grain protein, and test weight
## - Everything for plant height plus...
## - Elevated temperature during grain fill (Passarella et al 2005)
## - Drought during grain fill - (Savin and Nicholas  1996)
##
##
##

apriori_covariate_df <- trait_covariate_df %>%
  filter(str_detect(covariate, "soil", negate = TRUE)) %>%
  {bind_rows(
    filter(., trait %in% c("HeadingDate", "PlantHeight")),
    filter(., trait %in% c("GrainProtein", "GrainYield", "TestWeight"),
           covariate %in% c(subset(., trait == "PlantHeight", covariate, drop = TRUE),
                            "grain_fill.tmean_mean", "grain_fill.water_balance_sum"))
  )}


# Use the timeframes from the feature selection (stepwise) procedure
historical_apriori_feature_selection <- apriori_covariate_df %>% 
  group_by(trait) %>% 
  nest(.key = "covariates") %>% 
  mutate(covariates = map(covariates, "covariate")) %>%
  # Cross with the timeframes
  left_join(., distinct(historical_feature_selection, trait, time_frame, model)) %>%
  # Add feat sel type
  crossing(., feat_sel_type = c("apriori")) %>%
  # Modify covariates for models 4 and 5
  mutate(covariates = modify_if(covariates, model == "model5", ~c(.x, paste0("line_name:", .x))),
         covariates = modify_if(covariates, str_detect(feat_sel_type, "nosoil"), ~str_subset(.x, "soil", negate = TRUE)),
         covariates = map(covariates, ~list(optVariables = .x)))


# Create a data.frame of all covariates using the same timeframes as stepwise
historical_all_features <- trait_covariate_df %>%
  group_by(trait) %>% 
  nest(.key = "covariates") %>% 
  mutate(covariates = map(covariates, "covariate")) %>%
  # Cross with the timeframes
  left_join(., distinct(historical_feature_selection, trait, time_frame, model)) %>%
  # Add feat sel type
  crossing(., feat_sel_type = c("all", "all_nosoil")) %>%
  # Modify covariates for models 4 and 5
  mutate(covariates = modify_if(covariates, model == "model5", ~c(.x, paste0("line_name:", .x))),
         covariates = modify_if(covariates, str_detect(feat_sel_type, "nosoil"), ~str_subset(.x, "soil", negate = TRUE)),
         covariates = map(covariates, ~list(optVariables = .x)))


## Reload the concurrent data and save everything together
concurrent_file_list <- load(file.path(result_dir, "concurrent_feature_selection_results.RData"))

save(list = c("historical_feature_selection", "historical_feature_importance", "historical_apriori_feature_selection",
     "historical_all_features", "historical_timeframe_analysis", concurrent_file_list),
     file = file.path(result_dir, "feature_selection_results.RData"))
