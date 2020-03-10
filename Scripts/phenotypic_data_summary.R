## S2MET Phenotypic Adjustment
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## This notebook outlines procedures for calculating adjusted phenotypic means of
## entries in trials that belong to the `S2MET` experiment.
## 

# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))

## Load additional packages
library(Bilinear)
library(modelr)



## significance level
alpha <- 0.05



## Determine the number of environments per trait
S2_MET_BLUEs %>% 
  group_by(trait) %>% 
  summarize_at(vars(location, year, environment), n_distinct)

# trait        location  year environment
# 1 GrainProtein        7     3          10
# 2 GrainYield         13     3          22
# 3 HeadingDate        12     3          25
# 4 PlantHeight        13     3          27
# 5 TestWeight          8     2          12



## Basic Summaries

## Number of lines per environment per trait

## Find the total number of possible line x environment combinations and find
## the proportion that are observed for each trait
## If a trait was not observed in an entire environment, that environment is not
## included in these calculations
(prob_observed <- S2_MET_BLUEs %>% 
    distinct(trait, environment, line_name) %>%
    mutate_at(vars(line_name, environment), as.factor) %>%
    mutate(observed = TRUE) %>%
    split(.$trait) %>%
    map_df(~{
      droplevels(.) %>%
        complete(trait, environment, line_name, fill = list(observed = FALSE)) %>%
        summarize(trait = unique(trait), prop_obs = mean(observed))
    }))

# trait        prop_obs
# 1 GrainProtein    0.985
# 2 GrainYield      0.989
# 3 HeadingDate     0.994
# 4 PlantHeight     0.993
# 5 TestWeight      0.982





## Range of heritabilities
distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  group_by(trait) %>%
  summarize_at(vars(heritability), list(~min, ~max, ~mean))

# trait         min   max  mean
# 1 GrainProtein 0.165 0.838 0.612
# 2 GrainYield   0.180 0.900 0.564
# 3 HeadingDate  0.643 0.970 0.865
# 4 PlantHeight  0.106 0.883 0.525
# 5 TestWeight   0.587 0.945 0.716


## Output a table
intra_env_herit <- distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  mutate(heritability = formatC(heritability, digits = 2, width = 2, format = "g", flag = "#")) %>%
  spread(trait, heritability)

write_csv(x = intra_env_herit, path = file.path(fig_dir, "intra_environment_heritability.csv"))




## Stage-Two analysis


## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs, line_name %in% vp), population = "vp")) %>%
  filter(environment %in% tp_vp_env) %>%
  mutate_at(vars(location, year, line_name), as.factor)


# Boot reps
boot_reps <- 10

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP

#### Note
#### 
#### I am not sure whether to weight the residuals. It seems to result in a more realistic estimate of residual
#### variance, but it may not be correct.
#### 



### Analyze the components of GxE

# Lmer control
lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits <- S2_MET_BLUEs_tomodel %>%
  # filter(location %in% c("St_Paul", "Crookston", "Fargo", "Arlington", "Madison")) %>%
  group_by(trait, population) %>%
  do({
    
    df <- droplevels(.)
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + location + year, data = df)
    
    ## Harmonic means
    # Locations
    harm_loc <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Year
    harm_year <- apply(X = plot_table, MARGIN = c(2,3), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    harm_env <- xtabs(formula = ~ line_name + environment, data = df) %>%
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
      
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2,3), sum) %>% 
      harm_mean()
    
    # Get the weights
    wts <- df$std_error^2
    
    ## fit the full model
    fit_gly <- lmer(formula = value ~ 1 + location + year + location:year + (1|line_name) + (1|line_name:location) +
                  (1|line_name:year) + (1|line_name:location:year), data = df, control = lmer_control, weights = wts)
    
    fit_ge <- lmer(formula = value ~ 1 + environment + (1|line_name) + (1|line_name:environment), 
                   data = df, control = lmer_control, weights = wts)

    
    # ## Likelihood ratio tests
    # lrt_gly <- ranova(fit_gly) %>% 
    #   tidy() %>% 
    #   filter(!str_detect(term, "none")) %>% 
    #   mutate(term = str_remove(term, "\\(1 \\| ") %>% 
    #            str_remove("\\)")) %>% 
    #   select(term, LRT, df, p.value)
    # 
    # lrt_ge <- ranova(fit_ge) %>% 
    #   tidy() %>% 
    #   filter(!str_detect(term, "none")) %>% 
    #   mutate(term = str_remove(term, "\\(1 \\| ") %>% 
    #            str_remove("\\)")) %>% 
    #   select(term, LRT, df, p.value)
    
    
    ## Calculate heritability
    # Expression for heritability
    exp_gly <- "line_name / (line_name + (line_name:location / n_l) + (line_name:year / n_y) +(line_name:location:year / (n_l * n_y)) + (Residual / (n_y * n_l * n_r)))"
    exp_ge <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"
    
    
    ## Use bootstrapping to calculate a confidence interval
    # Generate bootstrapping samples and calculate heritability using the bootMer function
    h2_gly <- herit(object = fit_gly, exp = exp_gly, n_l = harm_loc, n_y = harm_year, n_r = harm_rep)
    h2_ge <- herit(object = fit_ge, exp = exp_ge, n_e = harm_env, n_r = harm_rep)
    
    
    # h2_boot <- bootMer(x = fit, nsim = boot_reps, FUN = function(x) 
    #   herit(object = x, exp = exp, n_l = harm_loc, n_y = harm_year, n_r = harm_rep)$heritability)
    # 
    # # Add the bootstrapping results with a confidence interval
    # h2$heritability <- tidy(h2_boot) %>% 
    #   cbind(., t(quantile(h2_boot$t, probs = c(alpha / 2, 1 - (alpha / 2))))) %>% 
    #   rename_at(vars(4, 5), ~c("lower", "upper"))
    
    
    # Return data_frame
    tibble(fit = list(fit_gly), fit_ge = list(fit_ge), h2 = list(h2_gly), h2_ge = list(h2_ge),
           # lrt = list(lrt), lrt = list(lrt),
           n_l = harm_loc, n_y = harm_year, n_e = harm_env, n_r = harm_rep)
    
  }) %>% ungroup()


# ## Look at LRT results
# stage_two_fits %>% 
#   select(trait, population, lrt) %>% 
#   unnest() %>% 
#   select(-df, -LRT) %>% 
#   spread(term, p.value)



## Create a table for output
stage_two_fits_GYL_varcomp <- stage_two_fits %>%
  mutate(var_comp = map(h2, "var_comp")) %>% 
  unnest(var_comp)


## Plot the proportion of variance from each source
stage_two_fits_GYL_varprop <- stage_two_fits_GYL_varcomp %>% 
  mutate(source = map(source, ~str_split(., pattern = ":", simplify = FALSE) %>% map(str_to_title) %>% .[[1]]) %>% 
           map_chr(~paste(., collapse = " x ")),
         source = str_replace_all(source, "Line_name", "Genotype"),
         source = factor(source, levels = c("Genotype", "Location", "Year", "Location x Year", "Genotype x Location", "Genotype x Year",
                                            "Genotype x Location x Year", "Residual")))

stage_two_fits_GYL_varprop1 <- stage_two_fits_GYL_varprop %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()

## Colors
var_comp_colors <- set_names(umn_palette(3, 8), levels(stage_two_fits_GYL_varprop$source))

## Plot both populations and all traits
g_varprop <- stage_two_fits_GYL_varprop1 %>% 
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of variance") +
  scale_fill_manual(values = var_comp_colors, name = NULL) +
  facet_grid(~ population) + 
  theme_presentation2(base_size = 10) +
  theme(axis.title.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))

# Save
ggsave(filename = "variance_components_expanded.jpg", plot = g_varprop, path = fig_dir, width = 9, height = 6, dpi = 1000)






#### 
#### Heading Date Analysis
#### 
#### Remove potential outlier environments and re-analyze
#### 

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_HD <- S2_MET_BLUEs_tomodel %>%
  filter(trait == "HeadingDate", ! environment %in% c("EON17", "CRM15", "HNY15")) %>%
  group_by(trait, population) %>%
  do({
    
    df <- droplevels(.)
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + location + year, data = df)
    
    ## Harmonic means
    # Locations
    harm_loc <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Year
    harm_year <- apply(X = plot_table, MARGIN = c(2,3), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2,3), sum) %>% 
      harm_mean()
    
    # Get the weights
    wts <- df$std_error^2
    
    ## fit the full model
    fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|location:year) + (1|line_name:location) + 
                  (1|line_name:year) + (1|line_name:location:year),
                data = df, control = lmer_control, weights = wts)
    
    # fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|line_name:location) + (1|line_name:year) + (1|line_name:location:year),
    #             data = df, control = lmer_control)
    
    
    ## Likelihood ratio tests
    lrt <- ranova(fit) %>% 
      tidy() %>% 
      filter(!str_detect(term, "none")) %>% 
      mutate(term = str_remove(term, "\\(1 \\| ") %>% 
               str_remove("\\)")) %>% 
      select(term, LRT, df, p.value)
    
    # Return data_frame
    data_frame(fit = list(fit), lrt = list(lrt), n_l = harm_loc, 
               n_y = harm_year, n_r = harm_rep)
    
  }) %>% ungroup()


## Create a table for output
stage_two_fits_HD_varcomp <- stage_two_fits_HD %>%
  mutate(var_comp = map(fit, ~as.data.frame(VarCorr(.)) %>% select(source = grp, variance = vcov)),
         var_comp = map2(var_comp, lrt, ~full_join(.x, .y, by = c("source" = "term")))) %>% 
  unnest(var_comp) %>%
  mutate(source = map(source, ~str_split(., pattern = ":", simplify = FALSE) %>% map(str_to_title) %>% .[[1]]) %>% 
           map_chr(~paste(., collapse = " x ")),
         source = str_replace_all(source, "Line_name", "Genotype"),
         source = factor(source, levels = c("Genotype", "Location", "Year", "Location x Year", "Genotype x Location", "Genotype x Year",
                                            "Genotype x Location x Year", "Residual"))) %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()








#######
#######
#######



## Look at the heritability and plot
g_herit <- stage_two_fits %>% 
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  # mutate(statistic = statistic + bias) %>%
  mutate(statistic = h2, lower = NA, upper = NA) %>%
  ggplot(aes(x = trait, y = statistic, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  # geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.40, label = round(statistic, 2)), position = position_dodge(0.9), size = 3) +
  # geom_hline(aes(yintercept = unbiased_statistic)) +
  scale_fill_brewer(palette = "Set1", name = "Population") +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = "bottom")

ggsave(filename = "heritability_expanded.jpg", plot = g_herit, path = fig_dir, width = 5, height = 4, dpi = 1000)










### Calculate the proportion of GxE that is due to environmental genetic variance
### heterogeneity versus lack of environmental correlation

# For each environment, calculate the genetic variance via reml
env_varG <- S2_MET_BLUEs_tomodel %>% 
  group_by(population, trait, environment) %>%
  do(varG = {
    df <- .
    wts <- df$std_error^2
    
    # If more than one trial is present, average over trials
    if (n_distinct(df$trial) > 1) {
      formula <- value ~ 1 + trial + (1|line_name)
    } else {
      formula <- value ~ 1 + (1|line_name)
    }
    
    suppressMessages(fit <- lmer(formula = formula, data = df, control = lmer_control, weights = wts) )
    
    as.data.frame(VarCorr(fit))[1,"vcov"]
    
  }) %>% ungroup() %>%
  mutate(varG = unlist(varG))

# Calculate the genetic heterogeneity, V
env_varG_V <- env_varG %>% 
  group_by(population, trait) %>% 
  summarize(V = var(sqrt(varG))) %>%
  ungroup()


## Use the variance components estinated in the previous random model
# prop_varcomp <- stage_two_fits_GE %>%
prop_varcomp <- stage_two_fits %>%
  mutate(varcomp = map(h2, "var_comp")) %>% 
  unnest(varcomp)


# Use the estimate of varGE across all environments to calculate L
env_L <- env_varG_V %>%
  # left_join(., subset(prop_varcomp, source == "line_name:environment", c(population, trait, variance))) %>% 
  left_join(., subset(prop_varcomp, source == "line_name:location:year", c(population, trait, variance))) %>% 
  rename(varGE = variance) %>%
  mutate(L = varGE - V)

# Use the estimate of genetic variance across all environments to calculate the 
# genetic correlation
env_r <- left_join(env_L, subset(prop_varcomp, source == "line_name", c(population, trait, variance)), by = c("population", "trait")) %>% 
  rename(varG = variance) %>%
  mutate(r_G = varG / (varG + L))

env_r %>% 
  select(population, trait, r_G) %>% 
  spread(population, r_G)

# trait          all    tp    vp
# 1 GrainProtein 0.442 0.460 0.300
# 2 GrainYield   0.357 0.362 0.272
# 3 HeadingDate  0.801 0.860 0.600
# 4 PlantHeight  0.401 0.393 0.421
# 5 TestWeight   0.451 0.417 0.459


## What proportion do V and L make up of varGE?
## This is from Li et al 2018 or Cooper and DeLacey 1994
## Add to the variance component table
varGE_components <- env_r %>% 
  mutate_at(vars(V, L), list(prop = ~. / varGE)) 


varGE_components %>%
  mutate(heterogeneity = str_c(round(V, 3), " (", round(V_prop, 2) * 100, "%)"), 
         lackCorrelation = str_c(round(L, 3), " (", round(L_prop, 2) * 100, "%)")) %>% 
  select(population, trait, heterogeneity, lackCorrelation) %>% 
  gather(grp, value, -trait, -population) %>% 
  spread(grp, value) %>%
  arrange(trait, population)

# population trait           heterogeneity  lackCorrelation 
# 1 all        GrainProtein 0.131 (25%)     0.403 (75%)     
# 2 tp         GrainProtein 0.13 (24%)      0.419 (76%)     
# 3 vp         GrainProtein 0.035 (8%)      0.412 (92%)     
# 4 all        GrainYield   84672.705 (33%) 170628.915 (67%)
# 5 tp         GrainYield   88238.353 (33%) 177820.796 (67%)
# 6 vp         GrainYield   51709.965 (25%) 154130.275 (75%)
# 7 all        HeadingDate  3.805 (58%)     2.803 (42%)     
# 8 tp         HeadingDate  3.867 (64%)     2.147 (36%)     
# 9 vp         HeadingDate  2.787 (52%)     2.559 (48%)     
# 10 all        PlantHeight  5.667 (29%)     13.941 (71%)    
# 11 tp         PlantHeight  4.87 (25%)      14.252 (75%)    
# 12 vp         PlantHeight  6.36 (32%)      13.456 (68%)    
# 13 all        TestWeight   129.551 (32%)   273.045 (68%)   
# 14 tp         TestWeight   133.357 (32%)   289.227 (68%)   
# 15 vp         TestWeight   115.139 (34%)   221.69 (66%)



## Output a table of all variance components for all populations
prop_varcomp1 <- prop_varcomp %>% 
  select(trait, population, source, variance) %>% 
  group_by(trait, population) %>% 
  mutate(proportion = variance / sum(variance),
         source = str_replace_all(source, ":", " x "), 
         source = str_replace_all(source, "line_name", "Genotype"), 
         source = str_to_title(source))


varGE_components1 <- varGE_components %>% 
  select(trait, population, V, L) %>% 
  gather(source, variance, V, L) %>%
  group_by(trait, population) %>% 
  mutate(source = ifelse(source == "V", "GeneticHeterogeneity", "LackOfCorrelation"),  
         proportion = variance / sum(variance))

# Component order
comp_order <- unique(stage_two_fits_GYL_varprop1$source) %>% 
  {.[order(str_count(., "x"))]} %>%
  {c(str_subset(., "[^Residual]"), "Genetic Heterogeneity", "Lack Of Correlation", "Residual")}

## Combine and output
var_comp_table <- select(stage_two_fits_GYL_varprop1, trait, population, source, variance, proportion = var_prop) %>%
  bind_rows(., varGE_components1) %>%
  ungroup() %>%
  mutate(significance = "",
         # significance = case_when(p.value < 0.001 ~ "***",
         #                          p.value < 0.01 ~ "**",
         #                          p.value < 0.05 ~ "*",
         #                          TRUE ~ ""),
         variance = signif(variance, 3) %>% formatC(x = ., digits = 3, format = "f") %>% str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         proportion = signif(proportion * 100, 3) %>% formatC(x = ., digits = 3, format = "f") %>% str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         annotation = paste0(variance, significance, " (", proportion, "%)"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, source, annotation) %>%
  mutate(source = str_add_space(source),
         trait = str_add_space(trait),
         source = factor(source, levels = comp_order),
         population = toupper(population)) %>% 
  arrange(population, source) %>%
  rename_all(str_to_title) %>%
  spread(Population, Annotation)

write_csv(x = var_comp_table, path = file.path(fig_dir, "population_variance_components_decomposed.csv"))


####

####






### Calculate the proportion of GxE that is due to location genetic variance
### heterogeneity versus lack of location correlation

## Use sommer
location_year_varG <- S2_MET_BLUEs_tomodel %>%
  group_by(trait, population) %>%
  nest() %>%
  mutate(out = list(NULL))

## Loop over rows
for (i in seq(nrow(location_year_varG))) {
  
  ## Data to model
  df <- location_year_varG$data[[i]]
  
  ## Fit a model per location
  location_varG <- df %>%
    group_by(location) %>%
    do(varG = {
      df1 <- .
      wts <- df1$std_error^2
      
      # If more than one trial is present, average over trials
      if (n_distinct(df1$trial) > 1) {
        formula <- value ~ 1 + trial + (1|line_name)
      } else {
        formula <- value ~ 1 + (1|line_name)
      }
      
      suppressMessages(fit <- lmer(formula = formula, data = df1, control = lmer_control, weights = wts) )
      
      as.data.frame(VarCorr(fit))[1,"vcov"]
      
    }) %>% ungroup() %>% unnest(varG)
  
  
  ## Fit a model per year
  year_varG <- df %>%
    group_by(year) %>%
    do(varG = {
      df1 <- .
      wts <- df1$std_error^2

      # If more than one trial is present, average over trials
      if (n_distinct(df1$trial) > 1) {
        formula <- value ~ 1 + trial + (1|line_name)
      } else {
        formula <- value ~ 1 + (1|line_name)
      }
      
      suppressMessages(fit <- lmer(formula = formula, data = df1, control = lmer_control, weights = wts) )
      
      as.data.frame(VarCorr(fit))[1,"vcov"]
      
    }) %>% ungroup() %>% unnest(varG)
  
  
  ## Add these df to the out list in location_year_varG
  location_year_varG$out[[i]] <- tibble(term = c("location", "year"), varG = list(location_varG, year_varG))
  
}

# Calculate the genetic heterogeneity, V
# do this for both locations and years
varG_V <- location_year_varG %>% 
  unnest(out) %>%
  mutate(V = map_dbl(varG, ~var(sqrt(.$varG))))

## Modify the prop_varcomp for use with varGL and varGY
prop_varcomp1 <- prop_varcomp %>%
  filter(source %in% c("line_name:location", "line_name:year")) %>%
  select(trait, population, term = source, variance) %>%
  mutate(term = str_remove_all(term, "line_name:"))


# Use the estimate of varGL or varGY across all environments to calculate L
varG_L <- varG_V %>%
  left_join(., prop_varcomp1) %>% 
  rename(varG_interaction = variance) %>%
  mutate(L = varG_interaction - V)

# Use the estimate of genetic variance across all environments to calculate the 
# genetic correlation
rG <- left_join(varG_L, subset(prop_varcomp, source == "line_name", c(population, trait, variance)), by = c("population", "trait")) %>% 
  rename(varG_broad = variance) %>%
  mutate(r_G = varG_broad / (varG_broad + L))

rG %>% 
  select(population, trait, term, r_G) %>% 
  spread(population, r_G)

# trait          all    tp    vp
# 1 GrainProtein 0.442 0.460 0.300
# 2 GrainYield   0.357 0.362 0.272
# 3 HeadingDate  0.801 0.860 0.600
# 4 PlantHeight  0.401 0.393 0.421
# 5 TestWeight   0.451 0.417 0.459


## What proportion do V and L make up of varGE?
## This is from Li et al 2018 or Cooper and DeLacey 1994
## Add to the variance component table
varGLY_components <- rG %>% 
  mutate_at(vars(V, L), list(prop = ~. / varG_interaction)) 


varGLY_components %>%
  filter(term == "location") %>%
  mutate(heterogeneity = str_c(round(V, 3), " (", round(V_prop, 2) * 100, "%)"), 
         lackCorrelation = str_c(round(L, 3), " (", round(L_prop, 2) * 100, "%)")) %>% 
  select(population, trait, heterogeneity, lackCorrelation) %>% 
  gather(grp, value, -trait, -population) %>% 
  spread(grp, value) %>%
  arrange(trait, population)

# population trait        heterogeneity   lackCorrelation
# 1 all        GrainProtein 0.113 (96%)     0.004 (4%)     
# 2 tp         GrainProtein 0.084 (77%)     0.025 (23%)    
# 3 vp         GrainProtein 0.048 (38%)     0.078 (62%)    
# 4 all        GrainYield   29415.707 (70%) 12497.899 (30%)
# 5 tp         GrainYield   29964.194 (71%) 12154.738 (29%)
# 6 vp         GrainYield   30732.615 (68%) 14350.573 (32%)
# 7 all        HeadingDate  0.473 (76%)     0.151 (24%)    
# 8 tp         HeadingDate  0.422 (58%)     0.307 (42%)    
# 9 vp         HeadingDate  1.086 (149%)    -0.359 (-49%)  
# 10 all        PlantHeight  1.628 (112%)    -0.174 (-12%)  
# 11 tp         PlantHeight  1.612 (98%)     0.025 (2%)     
# 12 vp         PlantHeight  1.658 (193%)    -0.8 (-93%)    
# 13 all        TestWeight   145.445 (120%)  -23.867 (-20%) 
# 14 tp         TestWeight   145.39 (116%)   -19.957 (-16%) 
# 15 vp         TestWeight   72.335 (88%)    9.514 (12%)



## Output a table of all variance components for all populations
prop_varcomp1 <- prop_varcomp %>% 
  select(trait, population, source, variance) %>% 
  group_by(trait, population) %>% 
  mutate(proportion = variance / sum(variance),
         source = str_replace_all(source, ":", " x "), 
         source = str_replace_all(source, "line_name", "Genotype"), 
         source = str_to_title(source))


varGLY_components1 <- varGLY_components %>% 
  filter(term == "location") %>%
  select(trait, population, V, L) %>% 
  gather(source, variance, V, L) %>%
  group_by(trait, population) %>% 
  mutate(source = ifelse(source == "V", "GeneticHeterogeneity", "LackOfCorrelation"),  
         proportion = variance / sum(variance))

# Component order
comp_order <- unique(stage_two_fits_GYL_varprop1$source) %>% 
  {.[order(str_count(., "x"))]} %>%
  {c(str_subset(., "[^Residual]"), "Genetic Heterogeneity", "Lack Of Correlation", "Residual")}

## Combine and output
var_comp_table <- select(stage_two_fits_GYL_varprop1, trait, population, source, variance, proportion = var_prop) %>%
  bind_rows(., varGLY_components1) %>%
  ungroup() %>%
  mutate(significance = "",
         # significance = case_when(p.value < 0.001 ~ "***",
         #                          p.value < 0.01 ~ "**",
         #                          p.value < 0.05 ~ "*",
         #                          TRUE ~ ""),
         variance = signif(variance, 3) %>% formatC(x = ., digits = 3, format = "f") %>% str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         proportion = signif(proportion * 100, 3) %>% formatC(x = ., digits = 3, format = "f") %>% str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         annotation = paste0(variance, significance, " (", proportion, "%)"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, source, annotation) %>%
  mutate(source = str_add_space(source),
         trait = str_add_space(trait),
         source = factor(source, levels = comp_order),
         population = toupper(population)) %>% 
  arrange(population, source) %>%
  rename_all(str_to_title) %>%
  spread(Population, Annotation)

write_csv(x = var_comp_table, path = file.path(fig_dir, "population_variance_components_gly.csv"))








###########################
## AMMI model fit
###########################

## Amend K to include non-genotyped lines
non_genotyped <- c(setdiff(tp, tp_geno), setdiff(vp, vp_geno))
K_ng <- diag(length(non_genotyped))
dimnames(K_ng) <- list(non_genotyped, non_genotyped)

# Add to K
K1 <- as.matrix(bdiag(K, K_ng))
dimnames(K1) <- replicate(2, c(colnames(K), colnames(K_ng)), simplify = FALSE)
# Sort names
K1 <- K1[order(colnames(K1)), order(colnames(K1))]


# Fit per trait
ammi_fit <- S2_MET_BLUEs_tomodel %>%
  filter(population == "tp") %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Create interaction factor
    df1 <- droplevels(df) %>%
      mutate(environment = as.factor(environment),
             ge = interaction(line_name, environment, sep = ":", drop = FALSE)) %>%
      # Create contrasts
      mutate_at(vars(line_name, environment), fct_contr_sum)
    
    # Create two-way table of genos and enviros
    ge_mat <- df1 %>% 
      select(line_name, environment, value) %>% 
      spread(environment, value) %>% 
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      as.matrix()
    
    ## Calculate correlation matrix between environments
    E_cor <- cor(x = ge_mat, use = "pairwise.complete.obs")
    # Use this to create G-E relationship matrix
    K_ge <- kronecker(K1[levels(df1$line_name), levels(df1$line_name)], E_cor, make.dimnames = TRUE)
    
    ## fit a mixed model to predict all genotype-environment means ##
    # MF
    mf <- model.frame(value ~ line_name + environment, df1)
    y <- model.response(mf)
    X <- model.matrix(~ 1 + line_name + environment, data = mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    
    fit1 <- mixed.solve(y = y, Z = Z, K = K_ge, X = X, method = "REML")
    
    ## Pull out effects
    grand_mean <- fit1$beta[1]
    g_effects <- fit1$beta %>% subset(., str_detect(names(.), "line_name"))
    e_effects <- fit1$beta %>% subset(., str_detect(names(.), "environment"))
    # Add last level
    g_effects <- c(g_effects, -sum(g_effects))
    e_effects <- c(e_effects, -sum(e_effects))
    
    # Create df
    g_effects_df <- tibble(line_name = levels(mf$line_name), effect = g_effects)
    e_effects_df <- tibble(environment = levels(mf$environment), effect = e_effects)
    
    ## Random effects
    rand_eff_df <- fit1$u %>% 
      tibble(term = names(.), effect = .) %>%
      separate(term, c("line_name", "environment"), sep = ":")
    
    ## Combine and calculate y_hat
    y_hat_df <- rand_eff_df %>%
      full_join(., g_effects_df, by = "line_name") %>%
      full_join(., e_effects_df, by = "environment") %>%
      mutate(y_hat = grand_mean + rowSums(select(., contains("effect")))) %>%
      select(-contains("effect"))
    
    ## Calculate prediction accuracy
    acc <- with(left_join(df1, y_hat_df, by = c("environment", "line_name")), cor(value, y_hat))
    
    
    ### Fit a GxE AMMI model ###
    
    # Fit the ammi model using the bilinear package
    # Use Ftest for speed - we don't care about significance
    fit_ammi <- bilinear(x = y_hat_df, G = "line_name", E = "environment", y = "y_hat", 
                         model = "AMMI", alpha = alpha, B = 1)
    
    ## Output a tibble
    tibble(y_hat = list(y_hat_df), model_acc = acc, fit_ammi = list(fit_ammi))
    
  }) %>% ungroup()




## Extract genotype and environment scores
ammi_scores <- ammi_fit %>% 
  mutate(anova = map(fit_ammi, "ANOVA") %>% map(~rownames_to_column(., "term")),
         # Calculate variance explained by PCs using eigenvalues or sums of squares
         varprop = map(fit_ammi, "svdE") %>% map(~tibble(PC = paste0("PC", seq_along(.$d)), eigenvalue = .$d, 
                                                         propvar_eigen = eigenvalue / sum(eigenvalue))),
         varprop = map2(varprop, anova, left_join, by = c("PC" = "term")),
         varprop = map(varprop, ~mutate(.x, propvar_SS = SS / sum(SS, na.rm = TRUE), PC_num = parse_number(PC)) %>% 
                         select(-Df, -MS, -testStatistic)))


## Display proportion of variance explained using normalized eigenvalues
ammi_scores %>%
  unnest(varprop) %>% 
  # Plot
  ggplot(aes(x = PC_num, y = propvar_SS)) +
  geom_point() + 
  facet_wrap(~ trait, scales = "free")

## Determine significant PCs using "elbow" method
# Tolerance for difference in variance explained
tol <- 0.03

ammi_sig_PCs <- ammi_scores %>%
  unnest(varprop) %>%
  #
  mutate(propvar = propvar_SS) %>%
  #
  arrange(trait, PC_num) %>%
  split(.$trait) %>% 
  ## Calculate the difference between steps of adding PCs. Find the first step when the difference is
  ## below the tolerance threshold
  map_df(~mutate(., propvar_diff = c(abs(diff(propvar)), 0), 
                 stop = which.min(propvar_diff >= tol), 
                 nPC = stop - 1))

## Summary df of number of sig PCs
ammi_sig_PCs_summ <- ammi_sig_PCs %>%
  group_by(trait) %>%
  filter(PC_num %in% seq(1, unique(nPC))) %>%
  summarize(total_propvar = sum(propvar), nPC = unique(nPC))


## Fit a model for that number of PCS
ammiN_fit <- ammi_fit %>%
  left_join(., ammi_sig_PCs_summ, by = "trait") %>%
  group_by(trait) %>%
  do({
    
    row <- .
    nPC <- row$nPC
    # Get the fitted ammi model
    fitted_ammi <- row$fit_ammi[[1]]
    
    # Get the environmental and genotypic effects
    g_effects <- fitted_ammi$Geffect %>%
      tibble(line_name = names(.), effect = .)
    e_effects <- fitted_ammi$Eeffect %>%
      tibble(environment = names(.), effect = .)
    
    # Get the environmental and genotypic scores
    g_scores <- fitted_ammi$scores$Gscores
    e_scores <- fitted_ammi$scores$Escores
    
    ## Sum the first nPC scores
    g_scores_sum <- rowSums(g_scores[,seq_len(nPC), drop = FALSE])
    e_scores_sum <- rowSums(e_scores[,seq_len(nPC), drop = FALSE])
    
    ## Combine into DF
    g_effects_scores <- cbind(g_effects, score = g_scores_sum, g_scores)
    e_effects_scores <- cbind(e_effects, score = e_scores_sum, e_scores)
    
    ## Predict y
    # First sum effects
    ge_effect_summ <- outer(
      X = g_effects_scores[,"effect"], 
      Y = e_effects_scores[,"effect"], 
      FUN = "+")
    
    ## Calculate phi, the fitted GxE effect vector using the principal components
    phi <- outer(X = g_scores_sum, Y = e_scores_sum)
    
    # Add the effects to phi, add intercept
    y_hat_mat <- c(fitted_ammi$mu) + ge_effect_summ + phi
    # Convert to df
    y_hat_df <- as.data.frame(y_hat_mat) %>%
      rownames_to_column("line_name") %>%
      gather(environment, y_hat, -line_name)
    
    
    ## Return tibble
    tibble(mu = fitted_ammi$mu, y_hat = list(y_hat_df), g_scores = list(g_effects_scores),
           e_scores = list(e_effects_scores), phi = list(phi))
    
  }) %>% ungroup()


## Calculate average scores by location
ammiN_fit_location <- ammiN_fit %>%
  mutate(loc_scores = map(e_scores, ~{
    # Add location information
    left_join(.x, distinct(trial_info, location, environment), by = "environment") %>%
      # Summarize by location
      group_by(location) %>% 
      summarize_at(vars(-environment), mean)
  }))
    

## Plot environment and location scores
ammi_gplot_list <- ammiN_fit_location %>%
  mutate(plot = map2(e_scores, loc_scores, ~{
    
    # Merge x and y for line segments
    df <- distinct(trial_info, location, environment) %>%
      right_join(., select(.x, environment, effect, score, PC1, PC2) %>% rename_at(vars(-environment), ~paste0("environment_", .)) ) %>%
      left_join(., select(.y, location, effect, score, PC1, PC2) %>% rename_at(vars(-location), ~paste0("location_", .)) ) %>%
      # Rename
      rename(x = location_PC1, xend = environment_PC1, y = location_PC2, yend = environment_PC2,
             x1 = location_effect, xend1 = environment_effect) %>%
      mutate(y1 = x, yend1 = xend)

    # Plot environment scores; then plot location scores
    g_pc1_pc1 <- ggplot(data = .x, aes(x = PC1, y = PC2)) +
      # geom_hline(yintercept = 0) + 
      # geom_vline(xintercept = 0) +
      geom_segment(data = df, aes(x = x, y = y, xend = 0, yend = 0), lwd = 1, color = "red") +
      geom_segment(data = df, aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.5, color = "blue") +
      # geom_point(color = "blue") +
      # geom_point(data = .y, shape = 15, size = 2, color = "red") +
      geom_text(aes(label = environment), color = "blue", size = 2) +
      geom_text(data = .y, aes(label = location), size = 4, color = "red") +
      theme_acs()
    
    # Plot environment scores; then plot location scores
    g_effect_pc1 <- ggplot(data = .x, aes(x = effect, y = PC1)) +
      # geom_hline(yintercept = 0) + 
      # geom_vline(xintercept = 0) +
      geom_segment(data = df, aes(x = x1, y = y1, xend = 0, yend = 0), lwd = 1, color = "red") +
      geom_segment(data = df, aes(x = x1, y = y1, xend = xend1, yend = yend1), lwd = 0.5, color = "blue") +
      # geom_point(color = "blue") +
      # geom_point(data = .y, shape = 15, size = 2, color = "red") +
      geom_text(aes(label = environment), color = "blue", size = 2) +
      geom_text(data = .y, aes(label = location), size = 4, color = "red") +
      theme_acs()
    
    # Return plots
    tibble(plot1 = list(g_pc1_pc1), plot2 = list(g_effect_pc1))
    
  }) ) %>% select(trait, plot)
    





## Save
save("ammi_fit", "ammiN_fit", "ammiN_fit_location", file = file.path(result_dir, "ammi_model_fit.RData"))







































##################################################
### Appendix
##################################################




# ## Visualization of distributions
# env_order <- S2_MET_BLUEs %>% 
#   distinct(environment, location, year) %>%
#   left_join(., data_frame(location = names(colors_use), color = colors_use)) %>%
#   mutate(location = factor(location, levels = loc_order)) %>% 
#   arrange(location, year) %>% 
#   {set_names(x = .$color, .$environment)}
# 
# ## Sort on Lat/Long and year
# S2_MET_BLUEs_toplot <- S2_MET_BLUEs_use %>%
#   mutate(environment = parse_factor(environment, levels = env_order))
# 
# 
# ## Plot
# g_met_dist <- S2_MET_BLUEs_toplot %>%
#   ggplot(aes(x = value, y = environment, fill = environment)) +
#   geom_density_ridges() +
#   facet_grid(. ~ trait, scales = "free_x") +
#   scale_fill_manual(values = env_order, guide = FALSE) +
#   ylab("Environment") +
#   xlab("Phenotypic value") +
#   theme_presentation2(base_size = 10)
# 
# # Save it
# ggsave(filename = "met_trait_dist.jpg", plot = g_met_dist, path = fig_dir, width = 4.5, height = 5, dpi = 1000)
# 
# 
# ## Combine map with distributions
# g_map_dist_combine <- plot_grid(g_map1, g_met_dist, ncol = 1, rel_heights = c(0.65, 1))
# ggsave(filename = "map_and_trait_dist.jpg", plot = g_map_dist_combine, path = fig_dir, width = 5, height = 6, dpi = 1000)
# 


#######
#######
#######


# ## Calculate genotype-location means
# S2_MET_Loc_BLUEs_model <- S2_MET_BLUEs_tomodel %>%
#   filter(population == "all") %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# for (i in seq(nrow(S2_MET_Loc_BLUEs_model))) {
#   
#   df <- S2_MET_Loc_BLUEs_model$data[[i]]
#   wts <- df$std_error^2
# 
#   # print(unique(df$trait))
#     
#   ## Fit a model to estimate genotype-location means
#   fit <- lmer(value ~ line_name + location + (1|year) + (1|line_name:year), data = df, weights = wts)
#   
#   # Get the effects
#   effs <- Effect(focal.predictors = c("line_name", "location"), mod = fit)
#   S2_MET_Loc_BLUEs_model$out[[i]] <- as.data.frame(effs)
#   
# }
# 
# S2_MET_Loc_BLUEs <- unnest(S2_MET_Loc_BLUEs_model, out)
# 
# ## Save
# save("S2_MET_Loc_BLUEs", file = file.path(data_dir, "S2MET_Location_BLUEs.RData"))




 