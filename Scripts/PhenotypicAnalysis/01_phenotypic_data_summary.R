## S2MET Phenotypic Adjustment
## 
## Analyze the phenotypic data across environments in this experiment.
## 
## Specifically, calculate:
## 1. variance components
## 2. correlations between environments
## 3. correlations between traits
## 

# Repository directory
repo_dir <- getwd()
# Source the main project script
source(file.path(repo_dir, "source.R"))

## Load additional packages
library(modelr)
library(broom)


## significance level
alpha <- 0.05



# High-level data summaries -----------------------------------------------




## Determine the number of environments per trait
S2_MET_BLUEs %>% 
  # Assign environment by train/test or external and split
  mutate(env_set = ifelse(environment %in% train_test_env, "trainTest", "external")) %>%
  group_by(trait, env_set) %>%
  summarize_at(vars(location, year, environment), n_distinct)

# trait        env_set   location  year environment
# 1 GrainProtein external         2     1           2
# 2 GrainProtein trainTest        7     3          11
# 3 GrainYield   external         3     2           4
# 4 GrainYield   trainTest       14     3          31
# 5 HeadingDate  trainTest       12     3          28
# 6 PlantHeight  external         1     2           2
# 7 PlantHeight  trainTest       13     3          31
# 8 TestWeight   external         1     2           2
# 9 TestWeight   trainTest        8     2          15



## Basic Summaries



## Range of heritabilities
distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  group_by(trait) %>%
  summarize_at(vars(heritability), list(min = min, max = max, mean = mean))

# trait         min   max  mean
# 1 GrainProtein 0.295 0.840 0.633
# 2 GrainYield   0.169 0.861 0.534
# 3 HeadingDate  0.174 0.980 0.848
# 4 PlantHeight  0.160 0.882 0.548
# 5 TestWeight   0.501 0.945 0.718


## Output a table
intra_env_herit <- distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  mutate(heritability = formatC(heritability, digits = 2, width = 2, format = "g", flag = "#")) %>%
  spread(trait, heritability)

write_csv(x = intra_env_herit, file = file.path(fig_dir, "intra_environment_heritability.csv"))




# Stage-two analysis ------------------------------------------------------


## Combine data
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  # Only look at train/test environments
  filter(environment %in% train_test_env) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter( (population == "all") |
            (population == "tp" & line_name %in% tp) |
            (population == "vp" & line_name %in% vp) ) %>%
  mutate_at(vars(location, year, environment, line_name), as.factor)


# Fit models using environmental means


# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_environment <- S2_MET_BLUEs_tomodel %>%
  group_by(trait, population) %>%
  nest() %>%
  ungroup() %>%
  arrange(trait, population) %>%
  mutate(out = list(NULL))

# First null
first_null <- min(which(sapply(stage_two_environment$out, is.null)))

for (i in seq(first_null, nrow(stage_two_environment))) {

  df <- droplevels(stage_two_environment$data[[i]])
  # Edit factors
  df1 <- df %>%
    droplevels() %>%
    mutate_at(vars(environment, line_name), ~fct_contr_sum(as.factor(.))) %>%
    mutate(wts = std_error^2, 
           gxe = interaction(line_name, environment, drop = TRUE, sep = ":"),
           value1 = scale(value))
  
  ## Model 1
  
  # Fit the model without scaled values
  fit1_out_unscaled <- lmer(formula = value ~ (1|line_name), data = df1, weights = wts,
                           control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  # Fit the model with scaled values
  fit1_out_scaled <- lmer(formula = value1 ~ (1|line_name), data = df1, weights = wts,
                          control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  ## Model 2
  
  # Fit the model without scaled values
  fit2_out_unscaled <- lmer(formula = value ~ (1|line_name) + (1|environment), data = df1, weights = wts,
                            control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  # Fit the model with scaled values
  fit2_out_scaled <- lmer(formula = value1 ~ (1|line_name) + (1|environment), data = df1, weights = wts,
                          control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  ## Model 3
  
  # Fit the model without scaled values
  fit3_out_unscaled <- lmer(formula = value ~ (1|line_name) + (1|environment) + (1|line_name:environment), 
                            data = df1, weights = wts,
                            control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  # Fit the model with scaled values
  fit3_out_scaled <- lmer(formula = value1 ~ (1|line_name) + (1|environment) + (1|line_name:environment), 
                          data = df1, weights = wts,
                          control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  # Combine models, extract variance components, assemble into data.frame
  model_varcomp_df <- ls(pattern = "fit.*_out_") %>% 
    setNames(., .) %>% 
    map(get) %>%
    tibble(model = names(.), fit = .) %>%
    mutate(varcomp = map(fit, ~as.data.frame(VarCorr(.x)))) %>%
    separate(model, c("model", "scaled"), sep = "_out_") %>%
    mutate(model = paste0("model", parse_number(model))) %>%
    select(-fit)
  
  # Get environmental means from model3
  environmental_means <- ranef(fit3_out_unscaled)$environment %>% 
    rownames_to_column("environment") %>%
    rename(effect = 2) %>%
    mutate(mean = effect + fixef(fit3_out_unscaled)[[1]])
  
  
  # ###
  # ### Fit using sommer ###
  # ### 
  # 
  # 
  # ## Construct covariance matrices
  # K1 <- K[levels(df1$line_name), levels(df1$line_name)]
  # ZE1 <- model.matrix(~ -1 + environment, df1)
  # ZG1 <- model.matrix(~ -1 + line_name, df1)
  # 
  # 
  # KE1 <- tcrossprod(ZG1 %*% K1, ZG1) * tcrossprod(ZG1)
  # dimnames(KE1) <- replicate(2, levels(df1$gxe), simplify = FALSE)
  # 
  # # Model 1
  # fit1_mmer <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = K1),
  #                   data = df1, weights = wts, verbose = TRUE)
  # 
  # # Model 2
  # fit1_mmer <- mmer(fixed = value ~ 1, random = ~ vs(line_name, Gu = K1) + environment,
  #                   data = df1, weights = wts, verbose = TRUE)
  # 
  # # Model 3
  # 
  # fit3_mmer <- mmer(fixed = value ~ 1, random = ~ environment + vs(line_name, Gu = K1) + vs(gxe, Gu = KE1),
  #                   data = df1, weights = wts, verbose = TRUE)
  # 
  # ## Get variance components
  # var_comp <- fit_mmer$sigma %>% 
  #   map_df(1) %>% 
  #   gather(term, variance) %>% 
  #   mutate(term = str_remove(term, "u:"),
  #          se = sqrt(diag(fit_mmer$sigmaSE)))
  # 
  # # Get environmental means
  # environmental_means <- fit_mmer$U$environment$value %>% 
  #   tibble(environment = names(.), effect = .) %>%
  #   mutate(environment = str_remove(environment, "environment"),
  #          mean = effect + fit_mmer$Beta$Estimate[1])
  # 
  # 
  
  ## Return the model
  stage_two_environment$out[[i]] <- tibble(var_comp = list(model_varcomp_df), env_mean = list(environmental_means))
  
}

# Unnest
stage_two_varcomp_env <- stage_two_environment %>%
  select(-data) %>%
  unnest(out)



# Fit models using location means


# Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs_tomodel %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  droplevels() %>%
  group_by(trait, population, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

# Group by trait and fit the model
stage_two_fits_loc <- S2_MET_loc_BLUEs %>%
  group_by(trait, population) %>%
  nest() %>%
  ungroup() %>%
  mutate(out = list(NULL))

# First null
first_null <- min(which(sapply(stage_two_fits_loc$out, is.null)))

for (i in seq(first_null, nrow(stage_two_fits_loc))) {
  
  df <- droplevels(stage_two_fits_loc$data[[i]])
  # Edit factors
  df1 <- df %>%
    droplevels() %>%
    mutate_at(vars(location, line_name), ~fct_contr_sum(as.factor(.))) %>%
    mutate(value1 = scale(value), wts = 1)
  
  ## Model 1
  
  # Fit the model without scaled values
  fit1_out_unscaled <- lmer(formula = value ~ (1|line_name), data = df1, weights = wts,
                            control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  # Fit the model with scaled values
  fit1_out_scaled <- lmer(formula = value1 ~ (1|line_name), data = df1, weights = wts,
                          control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  ## Model 2
  
  # Fit the model without scaled values
  fit2_out_unscaled <- lmer(formula = value ~ (1|line_name) + (1|location), data = df1, weights = wts,
                            control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  # Fit the model with scaled values
  fit2_out_scaled <- lmer(formula = value1 ~ (1|line_name) + (1|location), data = df1, weights = wts,
                          control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  ## Model 3
  
  # Fit the model without scaled values
  fit3_out_unscaled <- lmer(formula = value ~ (1|line_name) + (1|location) + (1|line_name:location), 
                            data = df1, weights = wts,
                            control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  # Fit the model with scaled values
  fit3_out_scaled <- lmer(formula = value1 ~ (1|line_name) + (1|location) + (1|line_name:location), 
                          data = df1, weights = wts,
                          control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  # Combine models, extract variance components, assemble into data.frame
  model_varcomp_df <- ls(pattern = "fit.*_out_") %>% 
    setNames(., .) %>% 
    map(get) %>%
    tibble(model = names(.), fit = .) %>%
    mutate(varcomp = map(fit, ~as.data.frame(VarCorr(.x)))) %>%
    separate(model, c("model", "scaled"), sep = "_out_") %>%
    mutate(model = paste0("model", parse_number(model))) %>%
    select(-fit)
  
  # Get location means from model3
  location_means <- ranef(fit3_out_unscaled)$location %>% 
    rownames_to_column("location") %>%
    rename(effect = 2) %>%
    mutate(mean = effect + fixef(fit3_out_unscaled)[[1]])
  
  
  ## Return the model
  stage_two_fits_loc$out[[i]] <- tibble(var_comp = list(model_varcomp_df), loc_mean = list(location_means))
  
}

# Unnest
stage_two_varcomp_loc <- stage_two_fits_loc %>%
  select(-data) %>%
  unnest(out)


# Save
save("stage_two_varcomp_env", "stage_two_varcomp_loc",  
     file = file.path(result_dir, "stage_two_phenotypic_analysis.RData"))















# ## Create a table for output
# stage_two_fits2_varcomp <- stage_two_fits1 %>%
#   mutate(var_comp = map(var_comp, ~select(., term = grp, variance = vcov)),
#          ranova = map(ranova, ~filter(., str_detect(term, "none", negate = TRUE)) %>%
#                         mutate(term = str_trim(str_remove_all(term, "\\(1 \\||\\)"))) ),
#          var_comp2 = map2(var_comp, ranova, full_join, by = "term")) %>%
#   unnest(var_comp2) %>%
#   group_by(trait, population) %>%
#   mutate(variance_prop = variance / sum(variance)) %>%
#   ungroup() %>%
#   mutate(variance = formatC(x = signif(variance, 3), digits = 3, width = 2, format = "fg", flag = "0"),
#          variance_prop = formatC(x = variance_prop * 100, digits = 3, width = 2, format = "fg"),
#          variance1 = paste0(variance, " / ", variance_prop, "%"),
#          p.value = str_trim(formatC(x = p.value, digits = 3, width = 3, format = "E"))) %>%
#   mutate(annotation = ifelse(p.value == "NA", variance1, paste0(variance1, " (P = ", p.value, ")"))) %>%
#   select(trait, population, term, annotation)
# 
# 
# ## Plot the proportion of variance from each source
# stage_two_fits2_varcomp1 <- stage_two_fits2_varcomp %>% 
#   mutate(term = map(term, ~str_split(., pattern = ":", simplify = FALSE) %>% 
#                       map(str_to_title) %>% .[[1]]) %>%  map_chr(~paste(., collapse = " x ")),
#          term = str_replace_all(term, "Line_name", "Genotype"),
#          term = factor(term, levels = c("Genotype", "Environment", "Genotype x Environment", 
#                                         "Residual")))
# 
# ## Save as table
# stage_two_fits2_varcomp1 %>%
#   spread(term, annotation) %>%
#   mutate(trait = str_add_space(trait),
#          population = ifelse(population == "all", "All", toupper(population))) %>%
#   rename_all(str_to_title) %>%
#   write_csv(x = ., path = file.path(fig_dir, "variance_components_table.csv"))
# 
# 
# ## Calculate environmental means
# stage_two_fits_env_mean <- stage_two_fits1 %>%
#   mutate(coef = map(fit_env_fixed, ~as.data.frame(summary(.)$coef) %>% rownames_to_column(., "term") %>%
#                       mutate(environment = str_remove(term, "environment"))),
#          coef = map2(coef, fit_env_fixed, ~add_row(.x, environment = last(levels(model.frame(.y)$environment)),
#                                                    Estimate = -sum(tail(.x$Estimate, -1))) %>%
#                        mutate(mean = Estimate + Estimate[1]) %>% filter(environment != "(Intercept)") %>% 
#                        select(environment, mean)) ) %>%
#   unnest(coef) %>% rename(pheno_mean = mean)
# 
# 
# 
# 
# ## Calculate range for the whole population
# stage_two_fits_env_mean %>%
#   filter(population == "all") %>%
#   group_by(trait) %>%
#   summarize_at(vars(pheno_mean), list(~min, ~mean, ~max))
# 
# 
# # trait            min   mean    max
# # 1 GrainProtein    9.96   12.5   14.7
# # 2 GrainYield   1774.   4407.  9037. 
# # 3 HeadingDate    50.8    61.4   71.0
# # 4 PlantHeight    45.1    72.0  107. 
# # 5 TestWeight    572.    648.   714.
# 
# 
# 
# 
# ### Calculate correlations between environments ###
# 
# env_cor_mat_list <- S2_MET_BLUEs_tomodel %>%
#   select(trait, population, environment, line_name, value) %>%
#   as.data.frame() %>%
#   group_by(trait, population) %>%
#   do(cor_mat = {
#     select(., line_name, environment, value) %>%
#       spread(environment, value) %>%
#       column_to_rownames("line_name") %>%
#       cor(use = "pairwise.complete.obs")
#   }) %>%
#   ungroup() %>%
#   mutate(cor_df = map(cor_mat, ~tidy(as.dist(.))),
#          cor_df = map(cor_df, ~rename_all(., ~c("environment", "environment2", "correlation"))))
# 
# 
# ## Summarize mean and range of the correlation
# cor_summ <- env_cor_mat_list %>%
#   unnest(cor_df) %>%
#   group_by(trait, population) %>%
#   do(summ = summary(.$correlation)) %>%
#   ungroup() %>%
#   mutate(min = map_dbl(summ, min), max = map_dbl(summ, max), mean = map_dbl(summ, "Mean")) %>%
#   select(-summ)
#   
# 
# # trait        population      min   max  mean
# # 1 GrainProtein all         0.113   0.534 0.347
# # 2 GrainProtein tp          0.0798  0.531 0.372
# # 3 GrainProtein vp         -0.0514  0.663 0.262
# # 4 GrainYield   all        -0.255   0.570 0.216
# # 5 GrainYield   tp         -0.307   0.623 0.217
# # 6 GrainYield   vp         -0.272   0.609 0.167
# # 7 HeadingDate  all         0.138   0.887 0.682
# # 8 HeadingDate  tp          0.563   0.898 0.742
# # 9 HeadingDate  vp         -0.122   0.830 0.476
# # 10 PlantHeight  all        -0.0477  0.660 0.331
# # 11 PlantHeight  tp         -0.00618 0.705 0.333
# # 12 PlantHeight  vp         -0.260   0.776 0.344
# # 13 TestWeight   all        -0.103   0.592 0.301
# # 14 TestWeight   tp         -0.0937  0.602 0.284
# # 15 TestWeight   vp         -0.269   0.661 0.303



 