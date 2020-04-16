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
library(car)
library(ggrepel)

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



############################
### Prepare ECs for modeling
############################


## Create the ECs and select the relevant ones for modeling
growth_stage_covariates1 <- growth_stage_covariates %>%
  ## Remove heading as a growth stage
  filter(stage != "heading") %>%
  ## Separate trial into location and year
  mutate(year = as.numeric(str_extract(trial, "[0-9]{4}")),
         location = str_remove(trial, "_[0-9]{4}")) %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca")) %>%
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
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca")) %>%
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


## Summarize covariates at a location over three time_frames
year_time_frame <- seq(1, 30)

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
  mutate(p_value = map_dbl(ks_test, "p.value")) %>%
  split(.$time_frame) %>%
  map_df(~mutate(., p_adj = p.adjust(p_value, method = "bonf")))

(ec_sig_nonnormal <- subset(ec_tomodel_normality, p_adj < 0.10))

# Check covariate by time_frame
# Determine which covariates to exclude
ec_to_remove <- group_by(ec_sig_nonnormal, covariate) %>% 
  summarize(n = n_distinct(time_frame)) %>% 
  filter(n > 20) %>% 
  pull(covariate)




## Visualize
ec_select_timeframe_summary_wide %>%
  gather(covariate, value, environmental_covariates) %>%
  filter(covariate %in% unique(ec_sig_nonnormal$covariate)) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_grid(time_frame ~ covariate)

## Remove this covariate from the df
ec_select_timeframe_summary_wide1 <- ec_select_timeframe_summary_wide %>%
  select(which(! names(.) %in% ec_to_remove) ) %>%
  group_by(time_frame) %>%
  # Inpute using the mean
  mutate_at(vars(-location, -time_frame), impute) %>%
  ungroup()



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
  select(-trial:-time_frame) %>%
  names()


## Model all covariates without apriori information
historical_ec_by_trait <- tibble(
  trait = sort(traits),
  covariates = list(
    environmental_covariates, # Grain protein
    environmental_covariates, # Grain yield
    str_subset(environmental_covariates, "flowering|grain", negate = TRUE), # Heading date
    str_subset(environmental_covariates, "grain", negate = TRUE), # Plant height
    environmental_covariates # Test weight
  )
)

# New vector of covariates
environmental_covariates_all <- reduce(historical_ec_by_trait$covariates, union)







############################
### Identify covariates that are correlated with the L and GL distance matrices
############################


## Calculate distance matrices between environmental main effects and GxE effects
## identified from the AMMI model, consistent with Rincent 2019

dist_matrices <- ammiN_fit_location %>%
  # Calculate a location PHI matrix
  mutate(loc_phi = map2(.x = g_scores, .y = loc_scores, ~outer(X = .x$score, Y = .y$score) %>%
                          `dimnames<-`(x = ., value = list(.x$line_name, .y$location)) )) %>%
  # Convert location effects into a matrix
  mutate(loc_scores = map(loc_scores, ~as.data.frame(.) %>% column_to_rownames(., "location")),
         loc_scores = map(loc_scores, ~as.matrix(select(., effect))),
         loc_phi = map(loc_phi, t)) %>%
  rename(main = loc_scores, int = loc_phi) %>%
  mutate_at(vars(main, int), list(~map(., make_dist_mat))) %>%
  select(trait, main, int)

## Add the coviarates assigned to each trait
historical_ec_tomodel <- left_join(dist_matrices, historical_ec_by_trait, by = "trait") %>%
  gather(term, matrix, -trait, -covariates)


# Create a matrix of scaled and centered covariates
historical_ec_tomodel_scaled_mat_list <- historical_ec_tomodel_scaled %>% 
  map(~as.data.frame(.) %>%
        column_to_rownames("location") %>% 
        select(-time_frame) %>%
        as.matrix() )


# Tolerance of difference of correlations
cor_tol <- 0.005


## For each trait, use an algorithm to determine the covariates to use
## This is based on the approach of Rincent 2019
historical_ec_dist <- historical_ec_tomodel %>%
  # Add timeframe as a character
  crossing(., time_frame = names(historical_ec_tomodel_scaled_mat_list)) %>%
  group_by(trait, term, time_frame) %>%
  do({
    
    row <- .
    
    # Relationship matrix
    relmat <- row$matrix[[1]]
    # Vector of covariates
    covariates <- row$covariates[[1]]
    # timeframe
    tf <- row$time_frame[1]
    
    ## Subset the EC scaled matrix
    ec_tomodel_scaled_mat_use <- historical_ec_tomodel_scaled_mat_list[[tf]][rownames(relmat), covariates]
    
    ## Choose the starting covariate as the one whose distance matrix has
    ## the highest correlation with the relmat distance matrix
    initial_W_cov <- map(covariates, ~make_dist_mat(ec_tomodel_scaled_mat_use[, ., drop = FALSE])) %>%
      map_dbl(~cor(as.numeric(.), as.numeric(relmat)))
    
    starting_covariate <- covariates[which.max(initial_W_cov)]
    # Create two vectors of available and used covariates
    used_covariates <- tested_covariates <- starting_covariate
    available_covariates <- setdiff(covariates, used_covariates)
    
    old_cor <- max(initial_W_cov)
    
    ## Vector of correlations
    i <- 1
    ec_cor_list <- vector("list", length(covariates))
    ec_cor_list[[i]] <- tibble(added_covariate = starting_covariate, cor_with_relmat = max(initial_W_cov))
    
    ## For each remaining covariate, find the next one that gives the highest gain (or lowest loss)
    while (length(available_covariates) > 0) {
      
      # Advance i
      i <- i + 1
      
      # Scan the available covariates, add to the matrix, calculate distance, compare to relmat
      new_W_cov <- map(available_covariates, ~c(used_covariates, .)) %>%
        map(~make_dist_mat(ec_tomodel_scaled_mat_use[, ., drop = FALSE])) %>%
        map_dbl(~cor(as.numeric(.), as.numeric(relmat)))
      
      old_cor <- max(new_W_cov)
      
      # Find the max difference
      which_max_diff <- which.max(new_W_cov - old_cor)
      
      # Get the covariate
      selected_covariate <- available_covariates[which_max_diff]
      
      ## Add to the list
      ec_cor_list[[i]] <- tibble(added_covariate = selected_covariate, cor_with_relmat = new_W_cov[which_max_diff])
      
      # Determine remaining covariates
      used_covariates <- union(used_covariates, selected_covariate)
      # Determine new available covariate vector
      available_covariates <- setdiff(covariates, used_covariates)
      
    }
    
    ## Collapse the list; calculate difference from previous
    ec_cor_df <- bind_rows(ec_cor_list) %>%
      mutate(difference = c(diff(cor_with_relmat), 0),
             stop = which.min(difference >= cor_tol),
             number = seq(nrow(.)))
    
    
    ## List of final covariates
    final_covariates <- subset(ec_cor_df, number <= stop, added_covariate, drop = TRUE)
    # Final distance matrix
    ec_dist_mat_final <- make_dist_mat(ec_tomodel_scaled_mat_use[, final_covariates, drop = FALSE])
    # final correlation
    final_cor <- subset(ec_cor_df, number == stop, cor_with_relmat, drop = TRUE)
    
    ## Return these results
    tibble(relmat = list(relmat), final_cor = final_cor, final_covariates = list(final_covariates), 
           ec_dist_mat_final = list(ec_dist_mat_final), test_results = list(ec_cor_df))
    
  }) %>% ungroup()


# Look at the best time frame per trait
(best_cor <- historical_ec_dist %>%
  group_by(trait, term) %>%
  top_n(x = ., n = 1, wt = final_cor) %>%
  slice(1))

# Compare with a series
best_cor %>%
  select(trait, time_frame, term, final_cor) %>%
  left_join(., historical_ec_dist %>% 
              filter(time_frame %in% paste0("time_frame", c(5, 10, 15, 20, 25, 30))) %>% 
              select(trait, term, time_frame1 = time_frame, final_cor1 = final_cor)) %>%
  as.data.frame()

# It lools like time_frame 5 is ideal


# Get results ready to plot
historical_ec_dist_toplot <- historical_ec_dist %>% 
  unnest(test_results) %>%
  # Find the stopping point and annotate
  mutate(annotation = ifelse(stop == number, stop, ""),
         term = ifelse(term == "int", "GxE", "E"))


## Plot the results - continuous
(g_hist_ec_dist <- historical_ec_dist_toplot %>% 
  mutate(time_frame = parse_number(as.character(time_frame))) %>%
  ggplot(aes(x = number, y = cor_with_relmat, color = time_frame, group = time_frame)) + 
  geom_point(aes(x = stop, y = final_cor)) +
  geom_line(lwd = 0.25) +
  facet_grid(term ~ trait, switch = "y") +
  scale_x_continuous(name = "Number of covariates", breaks = pretty) +
  scale_y_continuous(name = "Correlation with phenotypic relationship matrix", breaks = pretty) +
  scale_color_gradient(name = "Time frame (years)", breaks = pretty) +
  scale_linetype_discrete(name = NULL) +
  theme_light() +
  theme(legend.position = "top"))
ggsave(filename = "ec_locations_correlation__continuous.jpg", plot = g_hist_ec_dist, 
       path = fig_dir, width = 8, height = 5, dpi = 1000)






### Create relationship matrices for locations ####

location_relmat_df <- historical_ec_dist %>%
  select(trait, term, time_frame, K = relmat, EC = ec_dist_mat_final, final_covariates)


## Overlap in main-effect versus interaction ECs
location_relmat_covariates <- location_relmat_df %>%
  select(trait, time_frame, term, final_covariates) %>%
  spread(term, final_covariates) %>%
  mutate(covariate_overlap = map2(int, main, intersect)) %>%
  mutate_at(vars(int, main), list(covariate_unique = ~map2(., covariate_overlap, setdiff))) 

location_relmat_covariates %>%
  mutate_at(vars(-trait, -time_frame), ~map_dbl(., length)) %>%
  filter(time_frame %in% paste0("time_frame", c(5)))
  as.data.frame()

# trait        time_frame    int  main covariate_overlap int_covariate_unique main_covariate_unique
  # 1 GrainProtein time_frame5     4     2                 1                    3                     1
  # 2 GrainYield   time_frame5     1     3                 1                    0                     2
  # 3 HeadingDate  time_frame5     5     4                 1                    4                     3
  # 4 PlantHeight  time_frame5     1     2                 0                    1                     2
  # 5 TestWeight   time_frame5     2     5                 1                    1                     4


## Output a table of covariates for each term
historical_covariate_table <- location_relmat_df %>%
  filter(time_frame == "time_frame5") %>%
  select(trait, term, covariates = final_covariates)  %>% 
  unnest(covariates) %>%
  mutate(term = ifelse(term == "int", "GL", "L")) %>%
  group_by(trait, covariates) %>%
  summarize(term = paste0(term, collapse = "+")) %>%
  spread(trait, term) %>%
  arrange(covariates) %>%
  as.data.frame()
write_csv(x = historical_covariate_table, path = file.path(fig_dir, "historical_covariate_summary_table_TF5.csv"), na = "")


## Save these results
save("location_relmat_df", "historical_ec_tomodel_centered", 
     "historical_ec_tomodel_centers", "historical_ec_tomodel_scaled", "ec_select_timeframe",
     "historical_ec_dist",
     file = file.path(result_dir, "historical_ec_model_building.RData"))












#####################################
## Appendix
#####################################

############################
### Model covariates that are predictive of the location mean
############################

# ## Get the location effects from the AMMI model (average environmental effects)
# base_model_location_effect <- ammiN_fit_location %>%
#   unnest(loc_scores) %>% 
#   select(trait, location, effect)
# 
# 
# ## Use stepwise regression to find the best model
# location_effect_to_model <- base_model_location_effect %>%
#   left_join(., s2_met_tomodel)
# 
# 
# # Group by trait and time_frame and nest
# location_effect_models <- location_effect_to_model %>%
#   group_by(trait, time_frame) %>%
#   nest() %>%
#   mutate(model = list(NULL))
# 
# 
# 
# for (i in seq(nrow(location_effect_models))) {
#   
#   df <- location_effect_models$data[[i]]
#   tr <- location_effect_models$trait[i]
#   
#   ## Remove some covariates depending on when the trait is measured
#   df1 <- select(df, location, effect, subset(historical_ec_by_trait, trait == tr, covariates, drop = T)[[1]]) %>%
#     distinct()
#   
#   # Minimal model
#   min_model <- effect ~ 1
#   # Max model
#   max_model <- add_predictors(min_model, as.formula(paste0("~", paste0(names(df1)[-1:-2], collapse = " + "))))
#   
#   # Stepwise regression
#   fit_base <- lm(formula = min_model, data = df1, contrasts = "contr.sum")
#   fit_step <- step(object = fit_base, scope = max_model, direction = "forward")
#   
#   
#   # Did we run out of dfs? If so, reduce terms
#   while (df.residual(fit_step) <= 2) {
#     
#     ## Drop1
#     fit_step_drop1 <- drop1(fit_step) %>%
#       as.data.frame() %>%
#       rownames_to_column("term")
#     
#     ## Remove the two covariates that gives the lowest AIC, besides "none"
#     to_remove <- fit_step_drop1 %>%
#       filter(str_detect(term, "none", negate = T)) %>% 
#       top_n(x = ., n = 2, wt = -AIC) %>%
#       pull(term) %>% 
#       sort() %>%
#       first()
#     
#     ## Reduce the number of parameters by 1
#     fit_step_new_formula <- terms(formula(fit_step)) %>% 
#       drop.terms(termobj = ., dropx = which(attr(., "term.labels") %in% to_remove), keep.response = TRUE) %>%
#       formula()
#     
#     ## Refit the model
#     fit_step <- update(object = fit_step, formula = fit_step_new_formula)
#     
#   }
#   
#   fit_step1 <- fit_step
#   
#   
#   
#   
#   ## Remove covariates based on VIF
#   # Proceed only if the number of terms is >= 2
#   if (length(attr(terms(formula(fit_step1)), "term.labels")) >= 2) {
#     
#     vif_i <- vif(fit_step1)
#     fit_stepi <- fit_step1
#     vif_cutoff <- 4
#     
#     while (any(vif_i > vif_cutoff) & class(vif_i) != "try-error") {
#       
#       form_i <- terms(formula(fit_stepi)) %>%
#         drop.terms(., dropx = which( attr(., "term.labels") %in% names(which.max(vif_i))), keep.response = TRUE ) %>%
#         formula()
#       
#       fit_stepi <- update(object = fit_stepi, formula = form_i)
#       vif_i <- try( vif(fit_stepi), silent = TRUE)
#     }
#     
#     # Backward elimination from here
#     fit_step_final <- step(object = fit_stepi, direction = "backward")
#     
#   } else {
#     fit_step_final <- fit_step1
#     
#   }
#   
#   ## Return the model
#   location_effect_models$model[[i]] <- fit_step_final
#   
# }
# 
# 
# 
# 
# ## Plot fitted values versus observed environmental mean
# location_effect_models %>%
#   filter(time_frame %in% paste0("time_frame", c(5, 10, 15, 20))) %>%
#   mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>%
#   unnest(model_data) %>%
#   ggplot(aes(x = pred, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   scale_y_continuous(name = "effect", breaks = pretty) +
#   scale_x_continuous(name = "fitted", breaks = pretty) +
#   facet_wrap(~ trait + time_frame, scales = "free", ncol = 4) +
#   theme_acs()
# 
# 
# 
# ## Perform leave-one-out cross-validation
# location_effect_models_looCV <- location_effect_models %>%
#   mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>% 
#   group_by(trait, time_frame) %>%
#   do({
#     row <- .
#     cv_df <- unnest(row, model_data) %>%
#       crossv_loo() %>%
#       mutate_at(vars(test, train), ~map(., as.data.frame))
#     model_fits <- map(cv_df$train, ~update(object = row$model[[1]], data = .))
#     
#     # Predictions
#     predictions <- map2_df(cv_df$test, model_fits, add_predictions)
#     
#     # RMSE
#     rmse_df <- mean(map2_dbl(model_fits, cv_df$test, rmse))
#     
#     # Accuracy
#     accuracy <- cor(predictions$effect, predictions$pred)
#     
#     ## Return summaries
#     tibble(predictions = list(predictions), accuracy = accuracy, rmse = rmse_df)
#     
#   }) %>% ungroup()
# 
# 
# ## Plot predicted values versus observed environmental mean
# location_effect_models_looCV %>%
#   unnest(predictions) %>%
#   filter(time_frame %in% paste0("time_frame", c(5, 10, 15, 20))) %>%
#   ggplot(aes(x = pred, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   scale_y_continuous(name = "effect", breaks = pretty) +
#   scale_x_continuous(name = "prediction", breaks = pretty) +
#   facet_wrap(~ trait + time_frame, scales = "free", ncol = 4) +
#   theme_acs()
# 
# 
# 
# 
# ## Find the best model
# (best_model <- location_effect_models_looCV %>%
#     group_by(trait) %>% 
#     top_n(x = ., n = 1, wt = -rmse)
#   # top_n(x = ., n = 1, wt = accuracy)
# )
# 
# # What about time_frame 5?
# location_effect_models_looCV %>%
#   filter(time_frame == "time_frame5")
# 
# # Plot accuracy from the best model
# best_model %>%
#   unnest(predictions) %>%
#   ggplot(aes(x = pred, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   scale_y_continuous(name = "effect", breaks = pretty) +
#   scale_x_continuous(name = "prediction", breaks = pretty) +
#   facet_wrap(~ trait + time_frame, scales = "free", ncol = 3) +
#   theme_acs()









