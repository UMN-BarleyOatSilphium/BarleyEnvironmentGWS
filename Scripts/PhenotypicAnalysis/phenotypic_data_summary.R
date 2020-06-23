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
library(Bilinear)
library(modelr)
library(broom)


## significance level
alpha <- 0.05


## Determine the number of environments per trait
S2_MET_BLUEs %>% 
  # Assign environment by train/test or external and split
  mutate(env_set = ifelse(environment %in% train_test_env, "trainTest", "external")) %>%
  split(.$env_set) %>%
  map(~group_by(., trait) %>% 
      summarize_at(vars(location, year, environment), n_distinct))

# $external
# trait        location  year environment
# 1 GrainProtein        2     1           2
# 2 GrainYield          3     2           4
# 3 HeadingDate         1     2           2
# 4 PlantHeight         2     2           3
# 5 TestWeight          2     2           3
# 
# $trainTest
# trait        location  year environment
# 1 GrainProtein        7     3          11
# 2 GrainYield         14     3          30
# 3 HeadingDate        12     3          27
# 4 PlantHeight        13     3          30
# 5 TestWeight          8     2          15



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
  mutate_at(vars(location, year, environment, line_name), as.factor)


# Get the residuals variance within each environment
S2_MET_varR_tomodel <- s2_metadata %>%
  filter(trial %in% unique(S2_MET_BLUEs_tomodel$trial), trait %in% traits) %>%
  select(trial, trait, varR)


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
  group_by(trait, population) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq_len(nrow(stage_two_fits))) {

  df <- droplevels(stage_two_fits$data[[i]])
  # Edit factors
  df1 <- df %>%
    droplevels() %>%
    mutate_at(vars(location, year, line_name), fct_contr_sum) %>%
    mutate(wts = std_error^2)
  
  # Random model
  fit_lmer <- lmer(formula = value ~ 1 + (1|environment) + (1|line_name) + (1|line_name:environment),
                   data = df1, control = lmer_control, weights = wts)
  
  # Ranova
  fit_ranova <- lmerTest::ranova(model = fit_lmer)
  
  ## Return the model
  stage_two_fits$out[[i]] <- tibble(fit = list(fit_lmer), ranova = list(fit_ranova))
  
}

stage_two_fits1 <- unnest(stage_two_fits, out) %>%
  mutate(var_comp = map(fit, ~as.data.frame(VarCorr(.))), 
         ranova = map(ranova, tidy))


## Create a table for output
stage_two_fits2_varcomp <- stage_two_fits1 %>%
  mutate(var_comp = map(var_comp, ~select(., term = grp, variance = vcov)),
         ranova = map(ranova, ~filter(., str_detect(term, "none", negate = TRUE)) %>%
                        mutate(term = str_trim(str_remove_all(term, "\\(1 \\||\\)"))) ),
         var_comp2 = map2(var_comp, ranova, full_join, by = "term")) %>%
  unnest(var_comp2) %>%
  mutate(variance = formatC(x = signif(variance, 3), digits = 3, width = 3, format = "fg", flag = "0"),
         p.value = formatC(x = p.value, digits = 3, width = 3, format = "E")) %>%
  mutate(annotation = ifelse(p.value == "NA", variance, paste0(variance, " (P = ", p.value, ")"))) %>%
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



### Summarize heritability within environments ###

env_trait_herit %>%
  filter(environment %in% train_test_env) %>%
  group_by(trait) %>%
  summarize(herit_min = min(heritability), herit_max = max(heritability), herit_mean = mean(heritability))

# trait        herit_min herit_max herit_mean
# 1 GrainProtein     0.295     0.840      0.657
# 2 GrainYield       0.179     0.861      0.532
# 3 HeadingDate      0.174     0.980      0.845
# 4 PlantHeight      0.160     0.882      0.547
# 5 TestWeight       0.501     0.945      0.716

## Sort validation environment heritabilities
env_trait_herit %>%
  filter(environment %in% validation_env) %>%
  arrange(trait, desc(heritability))



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



 