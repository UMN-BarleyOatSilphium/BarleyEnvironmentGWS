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
library(lme4qtl)
library(modelr)
library(broom)


## significance level
alpha <- 0.05


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
# 4 GrainYield   trainTest       14     3          30
# 5 HeadingDate  external         1     2           2
# 6 HeadingDate  trainTest       12     3          27
# 7 PlantHeight  external         2     2           3
# 8 PlantHeight  trainTest       13     3          30
# 9 TestWeight   external         2     2           3
# 10 TestWeight   trainTest        8     2          15



## Basic Summaries



## Range of heritabilities
distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  group_by(trait) %>%
  summarize_at(vars(heritability), list(~min, ~max, ~mean))

# trait         min   max  mean
# 1 GrainProtein 0.295 0.840 0.633
# 2 GrainYield   0.169 0.861 0.539
# 3 HeadingDate  0.124 0.980 0.803
# 4 PlantHeight  0.160 0.882 0.546
# 5 TestWeight   0.501 0.945 0.713


## Output a table
intra_env_herit <- distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  mutate(heritability = formatC(heritability, digits = 2, width = 2, format = "g", flag = "#")) %>%
  spread(trait, heritability)

write_csv(x = intra_env_herit, path = file.path(fig_dir, "intra_environment_heritability.csv"))




## Stage-Two analysis


## Combine data
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  # Only look at train/test environments
  filter(environment %in% train_test_env) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter( (population == "all") |
            (population == "tp" & line_name %in% tp) |
            (population == "vp" & line_name %in% vp) ) %>%
  mutate_at(vars(location, year, environment, line_name), as.factor) %>%
  # Convert yield to Mg ha-1
  mutate_at(vars(value, std_error), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Rescale TW
  mutate_at(vars(value, std_error), ~ifelse(trait == "TestWeight", . / 100, .)) %>%
  # Rescale PH
  mutate_at(vars(value, std_error), ~ifelse(trait == "PlantHeight", . / 10, .))






# Partition variance for environmental means ------------------------------


# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits <- S2_MET_BLUEs_tomodel %>%
  group_by(trait, population) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in which(sapply(stage_two_fits$out, is.null))) {

  df <- droplevels(stage_two_fits$data[[i]])
  # Edit factors
  df1 <- df %>%
    droplevels() %>%
    mutate_at(vars(environment, line_name), ~fct_contr_sum(as.factor(.))) %>%
    mutate(wts = std_error^2, gxe = interaction(line_name, environment, drop = TRUE, sep = ":"))

  ## Construct covariance matrices
  K_use <- K[levels(df1$line_name), levels(df1$line_name)]
  KE_use <- tcrossprod(model.matrix(~ -1 + line_name, df1) %*% K_use, model.matrix(~ -1 + line_name, df1)) *
    tcrossprod(model.matrix(~ -1 + environment, df1))
  dimnames(KE_use) <- replicate(2, levels(df1$gxe), simplify = FALSE)
  
  fit_mmer <- mmer(fixed = value ~ 1, random = ~ environment + vs(line_name, Gu = K_use) + vs(gxe, Gu = KE_use),
                   data = df1, weights = wts, verbose = TRUE)
  
  ## Get variance components
  var_comp <- fit_mmer$sigma %>% 
    map_df(1) %>% 
    gather(term, variance) %>% 
    mutate(term = str_remove(term, "u:"),
           se = sqrt(diag(fit_mmer$sigmaSE)))
  
  # Get environmental means
  environmental_means <- fit_mmer$U$environment$value %>% 
    tibble(environment = names(.), effect = .) %>%
    mutate(environment = str_remove(environment, "environment"),
           mean = effect + fit_mmer$Beta$Estimate[1])
  
  
  ## Return the model
  stage_two_fits$out[[i]] <- tibble(var_comp = list(var_comp), env_mean = list(environmental_means))
  
}

# Unnest
stage_two_varcomp_env <- unnest(stage_two_fits, out) %>%
  select(-data) %>%
  # Rescale test weight
  mutate(var_comp = modify_if(var_comp, trait == "TestWeight", ~mutate(.x, variance = variance * (100^2), se = se * 100)),
         env_mean = modify_if(env_mean, trait == "TestWeight", ~mutate_if(.x, is.numeric, ~. * 100)))



# Partition variance for location means -----------------------------------


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
  mutate(out = list(NULL))

for (i in which(sapply(stage_two_fits_loc$out, is.null))) {
  
  df <- droplevels(stage_two_fits_loc$data[[i]])
  # Edit factors
  df1 <- df %>%
    droplevels() %>%
    mutate_at(vars(location, line_name), ~fct_contr_sum(as.factor(.))) %>%
    mutate(gxe = interaction(line_name, location, drop = TRUE, sep = ":"))
  
  ## Construct covariance matrices
  K_use <- K[levels(df1$line_name), levels(df1$line_name)]
  KE_use <- tcrossprod(model.matrix(~ -1 + line_name, df1) %*% K_use, model.matrix(~ -1 + line_name, df1)) *
    tcrossprod(model.matrix(~ -1 + location, df1))
  dimnames(KE_use) <- replicate(2, levels(df1$gxe), simplify = FALSE)
  
  fit_mmer <- mmer(fixed = value ~ 1, random = ~ location + vs(line_name, Gu = K_use) + vs(gxe, Gu = KE_use),
                   data = df1, verbose = TRUE)
  
  ## Get variance components
  var_comp <- fit_mmer$sigma %>% 
    map_df(1) %>% 
    gather(term, variance) %>% 
    mutate(term = str_remove(term, "u:"),
           se = sqrt(diag(fit_mmer$sigmaSE)))
  
  # Get environmental means
  location_means <- fit_mmer$U$location$value %>% 
    tibble(location = names(.), effect = .) %>%
    mutate(location = str_remove(location, "location"),
           mean = effect + fit_mmer$Beta$Estimate[1])
  
  
  ## Return the model
  stage_two_fits_loc$out[[i]] <- tibble(var_comp = list(var_comp), loc_mean = list(location_means))
  
}

# Unnest
stage_two_varcomp_loc <- unnest(stage_two_fits_loc, out) %>%
  select(-data) %>%
  # Rescale test weight
  mutate(var_comp = modify_if(var_comp, trait == "TestWeight", ~mutate(.x, variance = variance * (100^2), se = se * 100)),
         loc_mean = modify_if(loc_mean, trait == "TestWeight", ~mutate_if(.x, is.numeric, ~. * 100)))

# Save
save("stage_two_varcomp_env", "stage_two_varcomp_loc",file = file.path(result_dir, "stage_two_phenotypic_analysis.RData"))















## Create a table for output
stage_two_fits2_varcomp <- stage_two_fits1 %>%
  mutate(var_comp = map(var_comp, ~select(., term = grp, variance = vcov)),
         ranova = map(ranova, ~filter(., str_detect(term, "none", negate = TRUE)) %>%
                        mutate(term = str_trim(str_remove_all(term, "\\(1 \\||\\)"))) ),
         var_comp2 = map2(var_comp, ranova, full_join, by = "term")) %>%
  unnest(var_comp2) %>%
  group_by(trait, population) %>%
  mutate(variance_prop = variance / sum(variance)) %>%
  ungroup() %>%
  mutate(variance = formatC(x = signif(variance, 3), digits = 3, width = 2, format = "fg", flag = "0"),
         variance_prop = formatC(x = variance_prop * 100, digits = 3, width = 2, format = "fg"),
         variance1 = paste0(variance, " / ", variance_prop, "%"),
         p.value = str_trim(formatC(x = p.value, digits = 3, width = 3, format = "E"))) %>%
  mutate(annotation = ifelse(p.value == "NA", variance1, paste0(variance1, " (P = ", p.value, ")"))) %>%
  select(trait, population, term, annotation)


## Plot the proportion of variance from each source
stage_two_fits2_varcomp1 <- stage_two_fits2_varcomp %>% 
  mutate(term = map(term, ~str_split(., pattern = ":", simplify = FALSE) %>% 
                      map(str_to_title) %>% .[[1]]) %>%  map_chr(~paste(., collapse = " x ")),
         term = str_replace_all(term, "Line_name", "Genotype"),
         term = factor(term, levels = c("Genotype", "Environment", "Genotype x Environment", 
                                        "Residual")))

## Save as table
stage_two_fits2_varcomp1 %>%
  spread(term, annotation) %>%
  mutate(trait = str_add_space(trait),
         population = ifelse(population == "all", "All", toupper(population))) %>%
  rename_all(str_to_title) %>%
  write_csv(x = ., path = file.path(fig_dir, "variance_components_table.csv"))


## Calculate environmental means
stage_two_fits_env_mean <- stage_two_fits1 %>%
  mutate(coef = map(fit_env_fixed, ~as.data.frame(summary(.)$coef) %>% rownames_to_column(., "term") %>%
                      mutate(environment = str_remove(term, "environment"))),
         coef = map2(coef, fit_env_fixed, ~add_row(.x, environment = last(levels(model.frame(.y)$environment)),
                                                   Estimate = -sum(tail(.x$Estimate, -1))) %>%
                       mutate(mean = Estimate + Estimate[1]) %>% filter(environment != "(Intercept)") %>% 
                       select(environment, mean)) ) %>%
  unnest(coef) %>% rename(pheno_mean = mean)




## Calculate range for the whole population
stage_two_fits_env_mean %>%
  filter(population == "all") %>%
  group_by(trait) %>%
  summarize_at(vars(pheno_mean), list(~min, ~mean, ~max))


# trait            min   mean    max
# 1 GrainProtein    9.96   12.5   14.7
# 2 GrainYield   1774.   4407.  9037. 
# 3 HeadingDate    50.8    61.4   71.0
# 4 PlantHeight    45.1    72.0  107. 
# 5 TestWeight    572.    648.   714.




### Calculate correlations between environments ###

env_cor_mat_list <- S2_MET_BLUEs_tomodel %>%
  select(trait, population, environment, line_name, value) %>%
  as.data.frame() %>%
  group_by(trait, population) %>%
  do(cor_mat = {
    select(., line_name, environment, value) %>%
      spread(environment, value) %>%
      column_to_rownames("line_name") %>%
      cor(use = "pairwise.complete.obs")
  }) %>%
  ungroup() %>%
  mutate(cor_df = map(cor_mat, ~tidy(as.dist(.))),
         cor_df = map(cor_df, ~rename_all(., ~c("environment", "environment2", "correlation"))))


## Summarize mean and range of the correlation
cor_summ <- env_cor_mat_list %>%
  unnest(cor_df) %>%
  group_by(trait, population) %>%
  do(summ = summary(.$correlation)) %>%
  ungroup() %>%
  mutate(min = map_dbl(summ, min), max = map_dbl(summ, max), mean = map_dbl(summ, "Mean")) %>%
  select(-summ)
  

# trait        population      min   max  mean
# 1 GrainProtein all         0.113   0.534 0.347
# 2 GrainProtein tp          0.0798  0.531 0.372
# 3 GrainProtein vp         -0.0514  0.663 0.262
# 4 GrainYield   all        -0.255   0.570 0.216
# 5 GrainYield   tp         -0.307   0.623 0.217
# 6 GrainYield   vp         -0.272   0.609 0.167
# 7 HeadingDate  all         0.138   0.887 0.682
# 8 HeadingDate  tp          0.563   0.898 0.742
# 9 HeadingDate  vp         -0.122   0.830 0.476
# 10 PlantHeight  all        -0.0477  0.660 0.331
# 11 PlantHeight  tp         -0.00618 0.705 0.333
# 12 PlantHeight  vp         -0.260   0.776 0.344
# 13 TestWeight   all        -0.103   0.592 0.301
# 14 TestWeight   tp         -0.0937  0.602 0.284
# 15 TestWeight   vp         -0.269   0.661 0.303



 