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
library(modelr)
library(broom)
library(pbr)
library(cowplot)
library(pls)
library(car)
library(ggrepel)

# Significance level
alpha <- 0.05

## Environments to use with corresponding trials
env_trials <- S2_MET_BLUEs %>% 
  distinct(trial, environment) %>%
  filter(environment %in% tp_vp_env,
         str_detect(trial, "S2C1", negate = TRUE))



############################
### Load data
############################


## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))

## Load the fitted ammi models
load(file.path(result_dir, "ammi_model_fit.RData"))


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
### Prepare ECs for modeling
############################


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






## Prepare data for modelling:
## 1. only use the tp
## 2. add ECs to the df

s2_met_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  left_join(ec_tomodel_centered, by = "environment") %>%
  mutate_at(vars(line_name, environment), as.factor)


## What is the level of coverage in this dataset?
s2_met_tomodel %>% 
  group_by(trait, line_name) %>%
  summarize(n = n_distinct(environment)) %>% 
  mutate(mean_env = n / max(n)) %>%
  top_n(x = ., n = 5, wt = -mean_env)



# New vector of covariates
environmental_covariates <- s2_met_tomodel %>% 
  select(-trial:-std_error) %>%
  names()

## Scatterplot matrix
scatterplotMatrix(distinct(select(s2_met_tomodel, environmental_covariates)), smooth = FALSE)


## Determine if there is sufficient variation for a covariate
distinct(select(s2_met_tomodel, environmental_covariates)) %>%
  map_dbl(var) %>% 
  hist()



## Model all covariates without apriori information
ec_by_trait_all <- tibble(
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
environmental_covariates_all <- reduce(ec_by_trait_all$covariates, union)




############################
### Identify covariates that are correlated with the E and GE distance matrices
############################

## Calculate distance matrices between environmental main effects and GxE effects
## identified from the AMMI model, consistent with Rincent 2019

dist_matrices <- ammiN_fit_location %>%
  # Convert environmental effects into a matrix
  mutate(e_scores = map(e_scores, ~as.matrix(select(., effect))),
         phi = map(phi, t)) %>%
  rename(main = e_scores, int = phi) %>%
  mutate_at(vars(main, int), list(~map(., make_dist_mat))) %>%
  select(trait, main, int)

## Add the coviarates assigned to each trait
ec_tomodel <- left_join(dist_matrices, ec_by_trait_all, by = "trait") %>%
  gather(term, matrix, -trait, -covariates)

# Create a matrix of scaled and centered covariates
ec_tomodel_scaled_mat <- ec_tomodel_scaled %>% 
  as.data.frame() %>%
  column_to_rownames("environment") %>% 
  as.matrix()


# Tolerance of difference of correlations
cor_tol <- 0.005


## For each trait, use an algorithm to determine the covariates to use
## This is based on the approach of Rincent 2019
ec_dist <- ec_tomodel %>%
  group_by(trait, term) %>%
  do({
    
    row <- .
    
    # Relationship matrix
    relmat <- row$matrix[[1]]
    # Vector of covariates
    covariates <- row$covariates[[1]]
    
    ## Subset the EC scaled matrix
    ec_tomodel_scaled_mat_use <- ec_tomodel_scaled_mat[rownames(relmat), covariates]
    
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

# Get results ready to plot
ec_dist_toplot <- ec_dist %>% 
  unnest(test_results) %>%
  # Find the stopping point and annotate
  mutate(annotation = ifelse(stop == number, stop, ""),
         term = ifelse(term == "int", "GxE", "E"))


## Plot the results
(g_ec_dist <- ec_dist_toplot %>% 
  ggplot(aes(x = number, y = cor_with_relmat, lty = term)) + 
  geom_point(aes(x = stop, y = final_cor)) +
  geom_line() +
  geom_text_repel(aes(label = annotation), nudge_y = 0.03, nudge_x = 5) +
  facet_grid(~ trait) +
  scale_x_continuous(name = "Number of covariates", breaks = pretty) +
  scale_y_continuous(name = "Correlation with phenotypic\nrelationship matrix", breaks = pretty) +
  scale_linetype_discrete(name = NULL) +
  theme_light() +
  theme(legend.position = "top") )
# Save
ggsave(filename = "ec_correlation_with_relmat.jpg", plot = g_ec_dist, path = fig_dir, width = 5, height = 3, dpi = 1000)




### Create relationship matrices for environments ####

environmental_relmat_df <- ec_dist %>%
  select(trait, term, K = relmat, EC = ec_dist_mat_final, final_covariates)


## Overlap in main-effect versus interaction ECs
environmental_relmat_covariates <- environmental_relmat_df %>%
  select(trait, term, final_covariates) %>%
  spread(term, final_covariates) %>%
  mutate(covariate_overlap = map2(int, main, intersect)) %>%
  mutate_at(vars(int, main), list(covariate_unique = ~map2(., covariate_overlap, setdiff))) 

environmental_relmat_covariates %>%
  mutate_at(vars(-trait), ~map_dbl(., length))

# trait        int       main      covariate_overlap int_covariate_unique main_covariate_unique
# 1 GrainProtein     2     3                 1                    1                     2
# 2 GrainYield       2     3                 2                    0                     1
# 3 HeadingDate      3     2                 0                    3                     2
# 4 PlantHeight      4     3                 0                    4                     3
# 5 TestWeight       3     4                 1                    2                     3



## Output a table of covariates for each term
ec_covariate_table <- environmental_relmat_df %>%
  select(trait, term, covariates = final_covariates)  %>% 
  unnest(covariates) %>%
  mutate(term = ifelse(term == "int", "GE", "E")) %>%
  group_by(trait, covariates) %>%
  summarize(term = paste0(term, collapse = "+")) %>%
  spread(trait, term) %>%
  arrange(covariates) %>%
  as.data.frame()
write_csv(x = ec_covariate_table, path = file.path(fig_dir, "env_covariate_summary_table.csv"), na = "")



## Save these results
save("environmental_relmat_df", "ec_dist", "ec_tomodel_centered", "ec_tomodel_scaled", "ec_tomodel_centers", 
     file = file.path(result_dir, "ec_model_building.RData"))




































##################################
### Appendix
##################################

############################
### Model covariates that are predictive of the environmental mean
############################

# ## Get the environmental effects from the AMMI model
# base_model_env_effect <- ammiN_fit_location %>%
#   unnest(e_scores) %>% 
#   select(trait, environment, effect)
# 
# 
# 
# ## Use stepwise regression to find the best model
# env_effect_to_model <- base_model_env_effect %>%
#   left_join(s2_met_tomodel)
# 
# # Group by trait
# env_effect_models <- env_effect_to_model %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(model = list(NULL))
# 
# 
# 
# for (i in seq(nrow(env_effect_models))) {
#   
#   df <- env_effect_models$data[[i]]
#   tr <- env_effect_models$trait[i]
#   
#   ## Remove some covariates depending on when the trait is measured
#   df1 <- select(df, environment, effect, subset(ec_by_trait_all, trait == tr, covariates, drop = T)[[1]]) %>%
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
#   while (df.residual(fit_step) <= 1) {
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
#   ## Remove covariates based on VIF
#   vif_i <- vif(fit_step1)
#   fit_stepi <- fit_step1
#   vif_cutoff <- 4
#   
#   while (any(vif_i > vif_cutoff)) {
#     form_i <- terms(formula(fit_stepi)) %>%
#       drop.terms(., dropx = which( attr(., "term.labels") %in% names(which.max(vif_i))), keep.response = TRUE ) %>%
#       formula()
#     
#     fit_stepi <- update(object = fit_stepi, formula = form_i)
#     vif_i <- vif(fit_stepi)
#   }
#   
#   # Backward elimination from here
#   fit_step_final <- step(object = fit_stepi, direction = "backward")
#   
#   ## Return the model
#   env_effect_models$model[[i]] <- fit_step_final
#   
# }
# 
# 
# ## Examine coefficients of covariates
# env_effect_models$model %>% 
#   map(summary)
# 
# 
# ### This seems to have helped.
# 
# 
# 
# 
# ## Plot fitted values versus observed environmental mean
# env_effect_models %>%
#   mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>%
#   unnest(model_data) %>%
#   ggplot(aes(x = pred, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   facet_wrap(~ trait, scales = "free") +
#   theme_acs()
# 
# ## Perform leave-one-out cross-validation
# env_effect_models_looCV <- env_effect_models %>%
#   mutate(model_data = map(model, ~add_predictions(data = model.frame(.x), model = .x))) %>% 
#   group_by(trait) %>%
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
# # Plot
# env_effect_models_looCV %>%
#   unnest(predictions) %>%
#   ggplot(aes(x = pred, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   facet_wrap(~ trait, scales = "free") +
#   theme_acs()
# 
# 
# 
# # trait        predictions       accuracy     rmse
# # 1 GrainProtein <tibble [10 x 8]>    0.853    0.674
# # 2 GrainYield   <tibble [27 x 6]>    0.455 1357.   
# # 3 HeadingDate  <tibble [24 x 4]>    0.805    3.22 
# # 4 PlantHeight  <tibble [27 x 7]>    0.599    7.79 
# # 5 TestWeight   <tibble [12 x 8]>    0.849   19.6  












############################
### Visualization of daily weather data
############################


# ## Look at daily stats
# daily_ec_select <- growth_stage_weather %>%
#   inner_join(., env_trials) %>%
#   select(environment, dap, stage, mint:water_balance)
# 
# 
# ## Summarize max temperatures during grain fill
# grain_fill_maxt_summary <- daily_ec_select %>% 
#   filter(growth_stage == "grain_fill") %>%
#   ## Count number of days with 15-18 max temp,
#   ## < 15, 18-25 (moderate high), 25 - 30 (high), > 30 (very high)
#   mutate(grain_fill_condition = case_when(
#     maxt < 15 ~ "suboptimal",
#     between(maxt, 15, 18) ~ "optimal",
#     between(maxt, 18, 30) ~ "moderate_high",
#     between(maxt, 30, 35) ~ "high",
#     maxt > 35 ~ "very_high"
#   )) %>%
#   group_by(environment, grain_fill_condition) %>%
#   summarize(nDays = n())



# #### Example daily temperature for a drought, non-drought location-year
# 
# ## color vector for growth stage
# growth_stage_color <- setNames(neyhart_palette("umn2")[c(4, 5, 2)], names(growth_stage_rename))
# growth_stage_color["vegetative"] <- neyhart_palette()[5]
# 
# # what covariate to plot
# ec_to_plot <- "maxt"
# 
# daily_ec_example_plots <- daily_ec_select %>%
#   gather(covariate, value, -environment, -dap, -growth_stage) %>%
#   # Re-order growth stages
#   mutate(growth_stage = factor(growth_stage, levels = c("vegetative", "flowering", "grain_fill"))) %>%
#   filter(environment %in% c("EON16", "CRM16"), covariate == ec_to_plot) %>%
#   mutate(min_dap = min(dap), max_dap = max(dap)) %>%
#   ## calculate ranges for each covariate
#   group_by(covariate) %>%
#   mutate_at(vars(value), list(~min, ~max)) %>%
#   group_by(covariate, environment) %>%
#   do(plot = {
#     df <- .
#     
#     ## Create segments for growth stages
#     gs_seg <- df %>% 
#       group_by(growth_stage) %>% 
#       summarize(start = min(dap) - 1, end = max(dap))
#     
#     ec_name <- covariate_variable_rename[unique(df$covariate)]
#     ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
#     
#     # x axis limits
#     x_limit <- c(df$min_dap[1] - 1, df$max_dap[1])
#     
#     # y axis limits
#     y_limit <- c(df$min[1], df$max[1])
#     y_end <- quantile(y_limit, 0.75)
#     
#     ## Plot
#     ggplot(data = df, aes(x = dap, y = value)) +
#       geom_segment(data = gs_seg, mapping = aes(x = start, xend = end, y = y_end, yend = y_end, color = growth_stage), lwd = 10) +
#       geom_line() +
#       scale_y_continuous(name = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), breaks = pretty, limits = y_limit) +
#       scale_x_continuous(breaks = pretty, name = "Days after planting", limits = x_limit) +
#       scale_color_manual(values = growth_stage_color, name = "Predicted\ngrowth stage",
#                          labels = function(x) str_to_title(str_remove_all(x, "_"))) +
#       labs(subtitle = unique(df$environment)) +
#       theme_presentation2() +
#       theme(legend.position = "bottom")
#       
#   })
# 
#   
# # Extract sample plots
# plot_list <- daily_ec_example_plots %>%
#   pull(plot) %>%
#   map(~. + theme(legend.position = "none", axis.title = element_blank()))
# plot_list[[1]] <- plot_list[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# # Merge
# select_ec_plot_merge <- plot_grid(plotlist = plot_list, ncol = 1)
# 
# # y_axis label
# ec_name <- covariate_variable_rename[ec_to_plot]
# ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
# y_label <- grid::textGrob(label = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), rot = 90) 
# 
# x_label <- grid::textGrob(label = "Days after planting")
# 
# select_ec_plot_merge1 <- plot_grid(y_label, select_ec_plot_merge, nrow = 1, rel_widths = c(0.05, 1))  %>%
#   plot_grid(., x_label, get_legend(daily_ec_example_plots$plot[[1]]),ncol = 1, rel_heights = c(1, 0.05, 0.2)) 
# 
# 
# ggsave(filename = "sample_maxt_growth_stage.jpg", plot = select_ec_plot_merge1, path = fig_dir,
#        height = 6, width = 7, dpi = 1000)




# # what covariate to plot
# ec_to_plot <- "water_stress"
# 
# daily_ec_example_plots <- daily_ec_select %>%
#   gather(covariate, value, -environment, -dap, -growth_stage) %>%
#   # Re-order growth stages
#   mutate(growth_stage = factor(growth_stage, levels = c("vegetative", "flowering", "grain_fill"))) %>%
#   filter(environment %in% c("EON16", "CRM16"), covariate == ec_to_plot) %>%
#   mutate(min_dap = min(dap), max_dap = max(dap)) %>%
#   ## calculate ranges for each covariate
#   group_by(covariate) %>%
#   mutate_at(vars(value), list(~min, ~max)) %>%
#   group_by(covariate, environment) %>%
#   do(plot = {
#     df <- .
#     
#     ## Create segments for growth stages
#     gs_seg <- df %>% 
#       group_by(growth_stage) %>% 
#       summarize(start = min(dap) - 1, end = max(dap))
#     
#     ec_name <- covariate_variable_rename[unique(df$covariate)]
#     ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
#     
#     # x axis limits
#     x_limit <- c(df$min_dap[1] - 1, df$max_dap[1])
#     
#     # y axis limits
#     y_limit <- c(df$min[1], df$max[1])
#     y_end <- quantile(y_limit, 0.75)
#     
#     ## Plot
#     ggplot(data = df, aes(x = dap, y = value)) +
#       geom_segment(data = gs_seg, mapping = aes(x = start, xend = end, y = y_end, yend = y_end, color = growth_stage), lwd = 10) +
#       geom_line() +
#       scale_y_continuous(name = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), breaks = pretty, limits = y_limit) +
#       scale_x_continuous(breaks = pretty, name = "Days after planting", limits = x_limit) +
#       scale_color_manual(values = growth_stage_color, name = "Predicted\ngrowth stage",
#                          labels = function(x) str_to_title(str_remove_all(x, "_"))) +
#       labs(subtitle = unique(df$environment)) +
#       theme_presentation2() +
#       theme(legend.position = "bottom")
#     
#   })
# 
# 
# # Extract sample plots
# plot_list <- daily_ec_example_plots %>%
#   pull(plot) %>%
#   map(~. + theme(legend.position = "none", axis.title = element_blank()))
# plot_list[[1]] <- plot_list[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# # Merge
# select_ec_plot_merge <- plot_grid(plotlist = plot_list, ncol = 1)
# 
# # y_axis label
# ec_name <- covariate_variable_rename[ec_to_plot]
# ec_unit <- names(covariate_variable_unit[covariate_variable_unit == ec_name])
# y_label <- grid::textGrob(label = parse(text = paste0(str_to_title(ec_name), "~(", ec_unit, ")")), rot = 90) 
# 
# x_label <- grid::textGrob(label = "Days after planting")
# 
# select_ec_plot_merge1 <- plot_grid(y_label, select_ec_plot_merge, nrow = 1, rel_widths = c(0.05, 1))  %>%
#   plot_grid(., x_label, get_legend(daily_ec_example_plots$plot[[1]]),ncol = 1, rel_heights = c(1, 0.05, 0.2)) 
# 
# 
# ggsave(filename = "sample_maxt_growth_stage.jpg", plot = select_ec_plot_merge1, path = fig_dir,
#        height = 6, width = 7, dpi = 1000)



##########################
## Investiate time-dependent stressors
##########################

# ## Look at "number of days above some threshold temperature x during grainfill" with x = 15, 16, ..., max(T)
# # Determine the threshold that explains the most GxE variance
# temp_thresh <- seq(20, 40)
# 
# grain_fill_temperature_stress <- map(temp_thresh, ~{
#   growth_stage_weather %>%
#     filter(growth_stage == "grain_fill") %>%
#     group_by(environment) %>%
#     summarize(stress_days = sum(maxt > .x), threshold = .x) 
# })
# 
# ## format s2met data for modeling
# s2_met_tomodel <- S2_MET_BLUEs %>%
#   filter(line_name %in% tp) %>%
#   filter(trait == "GrainYield") %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(data = map(data, ~mutate(., line_name = as.factor(line_name),
#                                   line_name = `contrasts<-`(line_name, value = `colnames<-`(contr.sum(levels(line_name)), head(levels(line_name), -1))))
#   ))
# 
# 
# ## fit models
# s2_met_tomodel_stress_days <- s2_met_tomodel %>%
#   crossing(., stress = grain_fill_temperature_stress) %>%
#   mutate(data = map2(data, stress, ~left_join(.x, .y, by = "environment")),
#          threshold = map_dbl(stress, ~unique(.$threshold)))
# 
# 
# ## Extract p-values from anova
# stress_days_fit1 <- s2_met_tomodel_stress_days %>%
#   group_by(trait, threshold) %>%
#   do(fit = lm(value ~ line_name + line_name:stress_days, data = .$data[[1]])) %>%
#   ungroup() %>%
#   mutate(pvalue = map_dbl(fit, ~as.data.frame(anova(.))[2,5]))
# 
# 
# plot(-log10(pvalue) ~ threshold, stress_days_fit1)
# 
# ## Look at regression from lowest pvalue
# best_fit <- stress_days_fit1 %>%
#   filter(pvalue == min(pvalue, na.rm = TRUE)) %>%
#   pull(fit)
# 
# plot(best_fit[[1]])
# 


# ## Plot environment and maxt during grain fill
# daily_ec_select %>% 
#   filter(growth_stage == "grain_fill") %>%
#   mutate(environment = factor(environment, levels = ec_select$environment[order(ec_select$grain_fill_maxt, decreasing = TRUE)])) %>%
#   ggplot(aes(x = environment, y = maxt)) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 35, ymax = 40, fill = "heat stress"), alpha = 0.5) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 35, fill = "high temperature"), alpha = 0.5) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30, fill = "moderate high temperature"), alpha = 0.5) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 15, ymax = 20, fill = "optimal temperature"), alpha = 0.5) +
#   geom_violin(fill = alpha("white", 0.5)) +
#   # geom_boxplot(fill = NA) +
#   scale_y_continuous(breaks = pretty) +
#   scale_fill_manual(values = rev(viridis::inferno(direction = -1, n = 10)[1:4]), name = "Temperature stress level",
#                     guide = guide_legend(title.position = "top")) +
#   labs(subtitle = "Maximum temperature stress during grain fill") +
#   theme_presentation2(base_size = 12) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
# ggsave(filename = "grain_fill_max_temp_stress.jpg", path = fig_dir, width = 8, height = 5, dpi = 300)
# 

 







############################
### Model covariates that explain environment effect
############################

# ## Fit a base, random effect model for all traits
# base_model_fit <- s2_met_tomodel %>%
#   group_by(trait) %>%
#   do(fit = lmer(value ~ 1 + (1|line_name) + (1|environment), data = .)) %>%
#   ungroup()
# 
# 
# 
# 
# ## Get estimates of environmental effect
# base_model_env_effect <- base_model_fit %>%
#   mutate(effect = map(fit, ~ranef(.)$environment %>% rownames_to_column("environment") %>% rename(effect = 2))) %>%
#   unnest(effect)
#   
# ## Histogram
# par(mfrow = c(2,3))
# for (tr in unique(base_model_env_effect$trait)) {
#   hist(subset(base_model_env_effect, trait == tr, effect, drop = T), main = tr)
# }
# par(mfrow = c(1,1))
# 
# 
# # # environments as fixed
# # base_model_env_effect_alt <- s2_met_tomodel %>%
# #   group_by(trait) %>%
# #   do(fit = lmer(value ~ 1 + (1|line_name) + environment, data = ., contrasts = list(environment = "contr.sum"))) %>%
# #   ungroup() %>%
# #   mutate(mf = map(fit, model.frame),
# #          effect = map2(fit, mf, ~tibble(environment = levels(.y$environment), effect = c(fixef(.x)[-1], -sum(fixef(.x)[-1]))))) %>%
# #   unnest(effect)
# # 
# # ## Compare
# # left_join(base_model_env_effect, base_model_env_effect_alt, by = c("trait", "environment")) %>%
# #   qplot(x = effect.x, y = effect.y, data = .) +
# #   facet_wrap(~ trait, scales = "free") +
# #   geom_abline(slope = 1, intercept = 0)
# # 
# # # Good! - proceed with random effect
# 
# 
# 
# 
# ## Use stepwise regression to find the best model
# env_effect_to_model <- base_model_env_effect %>%
#   left_join(s2_met_tomodel)
# 
# # Group by trait
# env_effect_models <- env_effect_to_model %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(model = list(NULL))
# 
# for (i in seq(nrow(env_effect_models))) {
#   
#   df <- env_effect_models$data[[i]]
#   tr <- env_effect_models$trait[i]
#   
#   ## Remove some covariates depending on when the trait is measured
#   df1 <- select(df, environment, effect, subset(ec_by_trait, trait == tr, covariates, drop = T)[[1]])
#   
#   # Minimal model
#   min_model <- effect ~ 1
#   # Max model
#   max_model <- add_predictors(min_model, as.formula(paste0("~", paste0(names(df1)[-1:-2], collapse = " + "))))
#   
#   # Stepwise regression
#   fit_base <- lm(formula = min_model, data = df1)
#   fit_step <- step(object = fit_base, scope = max_model, direction = "both")
#   
#   ## Return the model
#   env_effect_models$model[[i]] <- fit_step
#   
# }
# 
# 
# 
# ## Plot
# env_effect_models %>% 
#   mutate(data = map2(data, model, ~mutate(.x, predicted_effect = predict(.y)))) %>%
#   unnest(data) %>%
#   ggplot(aes(x = predicted_effect, y = effect)) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point() +
#   facet_wrap(~ trait, scales = "free") +
#   theme_acs()






# ############################
# ### Model covariates that explain GxE
# ############################
# 
# 
# 
# #####
# # Model building
# #####
# 
# 
# ## Data.frame to store final models and covariates
# ec_model_building <- s2_met_tomodel %>% 
#   distinct(trait) %>%
#   mutate(all_model = list(NULL), apriori_model = list(NULL), final_model = list(NULL))
# 
# 
# 
# 
# 
# 
# ## 
# ## Grain yield, grain protein, and test weight
# ## 
# ## Start with the following a prior covariates:
# ## 
# ## tmin during flowering
# ## radiation during vegetative
# ## water stress during flowering
# ## grain fill water stress
# ## grain fill tmax
# ## 
# ## 
# ## For heading date and plant height, just use all covariates and backwards elimination
# ## 
# ## 
# 
# 
# ## Control
# control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
# 
# 
# ## Loop over traits
# for (i in seq(nrow(ec_model_building))) {
#   
#   # Get the trait
#   tr <- ec_model_building$trait[i]
#   
#   # Subset the data
#   tr_data <- filter(s2_met_tomodel, trait == tr)
#   
#   # Get the base model
#   tr_base_lmer <- lmer(value ~ 1 + (1|line_name) + environment, data = tr_data)
#   tr_base_lm <- lm(value ~ 1 + line_name + environment, data = tr_data)
#   
#   
#   
#   # ### Apriori model ###
#   # 
#   # ## Assign a priori covariates
#   # apriori_tr_covariates <- subset(ec_by_trait_apriori, trait == tr, covariates, drop = T)[[1]]
#   # 
#   # ## Create a new formula
#   # # tr_covariate_form1 <- as.formula(paste0("~", paste0(apriori_tr_covariates, collapse = " + "), " + ", 
#   # #                                        paste0("line_name:", apriori_tr_covariates, collapse = " + ")))
#   # 
#   # # Diagonal covariance structure
#   # tr_covariate_form2 <- as.formula(paste0("~ ", paste0("(0 +", apriori_tr_covariates, "|line_name)", collapse = " + ")))
#   # # # Unstructured covariance
#   # # tr_covariate_form2 <- as.formula(paste0("~ (0 +", paste0(apriori_tr_covariates, collapse = " + "), "|line_name)"))
#   # 
#   # # gxe_formula1 <- add_predictors(formula(tr_base_model), tr_covariate_form1)
#   # gxe_formula2 <- add_predictors(formula(tr_base_model), tr_covariate_form2)
#   # 
#   # 
#   # # Refit the model
#   # # tr_ec_model_mixed1 <- lmer(formula = gxe_formula1, data = tr_data)
#   # tr_ec_model_mixed2_apriori <- lmer(formula = gxe_formula2, data = tr_data, control = control)
#   # 
#   # 
#   # # Stepwise elimination of extract covariates
#   # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed2_apriori)
#   # # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed, reduce.random = FALSE, alpha.fixed = alpha)
#   # 
#   # # Get the model
#   # tr_ec_model_mixed_apriori_final <- get_model(tr_ec_model_mixed_step)
# 
# 
#   
#   
#   ### All EC model ###
#   # assign covariates
#   all_tr_covariates <- subset(ec_by_trait_all, trait == tr, covariates, drop = T)[[1]]
#   
#   ## Create a new formula
#   # LMER
#   # tr_covariate_form <- as.formula(paste0("~ (0 +", paste0(all_tr_covariates, collapse = " + "), "|line_name)"))
#   tr_covariate_form_lmer <- as.formula(paste0("~ ", paste0("(0 +", all_tr_covariates, "|line_name)", collapse = " + ")))
#   gxe_formula_lmer <- add_predictors(formula(tr_base_lmer), tr_covariate_form_lmer)
#   
#   # Refit the model
#   tr_ec_model_mixed_all <- lmer(formula = gxe_formula_lmer, data = tr_data, control = control)
#   
#   # Stepwise elimination of extra covariates
#   tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed_all)
#   # tr_ec_model_mixed_step <- step(object = tr_ec_model_mixed_all, reduce.random = FALSE, alpha.fixed = alpha)
#   
#   # Get the model
#   tr_ec_model_mixed_all_final <- get_model(tr_ec_model_mixed_step)
#   
#   
#   
#   
#   # LM
#   tr_covariate_form_lm <- as.formula(paste0("~ ", paste0(all_tr_covariates, ":line_name", collapse = " + ")))
#   gxe_formula_lm <- add_predictors(formula(tr_base_lm), tr_covariate_form_lm)
#   
#   ## Stepwise addition of covariates
#   tr_ec_model_fixed_step <- step(object = tr_base_lm, scope = gxe_formula_lm, direction = "both")
#   
#   test <- lm(value ~ 1 + line_name + environment + line_name:grain_fill.radn_sum, tr_data)
#   
#   
#   # Refit the model
#   tr_ec_model_fixed_all <- lm(formula = gxe_formula_lm, data = tr_data)
#   
#   # Stepwise elimination of extra covariates
#   
# 
#   
# 
#   
# 
#   
#   
#   
#   
#   ## Add results to the tibble
#   ec_model_building$apriori_model[i] <- list(tr_ec_model_mixed2_apriori)
#   ec_model_building$final_model[i] <- list(tr_ec_model_mixed_apriori_final)
#   ec_model_building$all_model[i] <- list(tr_ec_model_mixed_all_final)
#   
# }
# 
# 
# ## Tidy the output
# ec_model_building_tidy <- ec_model_building %>%
#   gather(model, object, -trait) %>%
#   mutate(model = str_remove(model, "_model"),
#          # Calculate varprop
#          varprop = map(object, ~as.data.frame(VarCorr(.)) %>% 
#                          mutate(term = case_when(is.na(var1) ~ grp, var1 == "(Intercept)" ~ grp, TRUE ~ var1), 
#                                 variance = vcov, varprop = variance / sum(variance)) ),
#          # Calculate statistics
#          AIC = map_dbl(object, AIC), BIC = map_dbl(object, BIC), logLik = map_dbl(object, logLik))




# 
# 
# ############################
# ### Heritability
# ############################
# 
# 
# ## Fit genomic heritability models for slopes
# # Extract the coefficients for each model
# # Keep final and apriori models
# ec_interaction_coef <- covariate_reg_coefs %>%
#   rename_all(~str_remove(., "_model")) %>%
#   gather(model, out, -trait) %>% 
#   unnest(out)
# 
# 
# # Subset the K matrix
# K_use <- K[tp_geno, tp_geno]
# 
# ## Fit models
# ec_interaction_coef_herit <- ec_interaction_coef %>%
#   group_by(trait, model, covariate) %>%
#   do({
#     df <- .
#     
#     df1 <- subset(df, line_name %in% tp_geno)
#     
#     ## Calculate heritability
#     invisible(capture.output(herit_fit <- marker_h2(data.vector = df1$estimate, geno.vector = df1$line_name, K = K_use, alpha = alpha)))
#     
#     ## Return df
#     tibble(heritability = herit_fit$h2, lower = herit_fit$conf.int1[1], upper = herit_fit$conf.int1[2])
#     
#   })
# 
# ## Write table
# write_csv(x = ec_interaction_coef_herit, path = file.path(fig_dir, "covariate_slope_heritability.csv"))
# 


