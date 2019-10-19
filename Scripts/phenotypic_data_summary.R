## S2MET Phenotypic Adjustment
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## This notebook outlines procedures for calculating adjusted phenotypic means of
## entries in trials that belong to the `S2MET` experiment.
## 

# Load libraries and directories

library(lme4)
library(lmerTest)
library(nlme)
library(modelr)
library(broom)

# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))


## significance level
alpha <- 0.05



## Determine the number of environments per trait
S2_MET_BLUEs %>% 
  group_by(trait) %>% 
  summarize_at(vars(location, year, environment), n_distinct)

# trait       location  year environment
# 1 GrainProtein        7     3          11
# 2 GrainYield         13     3          23
# 3 HeadingDate        12     3          26
# 4 PlantHeight        13     3          28
# 5 TestWeight          8     2          14



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
# 1 GrainProtein    0.982
# 2 GrainYield      0.986
# 3 HeadingDate     0.992
# 4 PlantHeight     0.991
# 5 TestWeight      0.980





## Range of heritabilities
distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  group_by(trait) %>%
  summarize_at(vars(heritability), list(~min, ~max, ~mean))

# trait         min   max  mean
# 1 GrainProtein 0.165 0.838 0.626
# 2 GrainYield   0.180 0.900 0.551
# 3 HeadingDate  0.643 0.970 0.869
# 4 PlantHeight  0.106 0.883 0.531
# 5 TestWeight   0.587 0.945 0.727


## Output a table
intra_env_herit <- distinct(S2_MET_BLUEs, trait, environment) %>%
  left_join(., env_trait_herit) %>%
  mutate(heritability = formatC(heritability, digits = 2, width = 2, format = "g", flag = "#")) %>%
  spread(trait, heritability)

write_csv(x = intra_env_herit, path = file.path(fig_dir, "intra_environment_heritability.csv"))














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










## Stage-Two analysis


## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs, line_name %in% vp), population = "vp")) %>%
  filter(environment %in% tp_vp_env)


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
  mutate_at(vars(location, year, line_name), as.factor) %>%
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
    
    fit_ge <- lmer(formula = value ~ 1 + environment + (1|line_name) + (1|line_name:environment), data = df, control = lmer_control, weights = wts)

    
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
  mutate_at(vars(location, year, line_name), as.factor) %>%
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



## Look at genetic components
stage_two_fits_GYL_varprop2 <- stage_two_fits_GYL_varprop %>%
  filter(str_detect(source, "Genotype")) %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()



## Plot both populations and all traits
g_varprop_gen <- stage_two_fits_GYL_varprop2 %>% 
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of genetic-related variance") +
  scale_fill_manual(values = var_comp_colors, name = NULL) +
  facet_grid(~ population) + 
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_genetic_expanded.jpg", plot = g_varprop_gen, path = fig_dir, width = 8, height = 4, dpi = 1000)





## Look at the heritability and plot
g_herit <- stage_two_fits_GYL %>% 
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
prop_varcomp <- stage_two_fits_GYL %>%
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
# 1 GrainProtein 1.11  1.03  0.745
# 2 GrainYield   0.712 0.762 0.530
# 3 HeadingDate  0.889 0.923 0.670
# 4 PlantHeight  1.03  0.981 1.31 
# 5 TestWeight   0.594 0.597 0.469


## What proportion do V and L make up of varGE?
## This is from Li et al 2018 or Cooper and DeLacey 1994
## Add to the variance component table
varGE_components <- env_r %>% 
  mutate_at(vars(V, L), funs(prop = . / varGE)) 


varGE_components %>%
  mutate(heterogeneity = str_c(round(V, 3), " (", round(V_prop, 2) * 100, "%)"), 
         lackCorrelation = str_c(round(L, 3), " (", round(L_prop, 2) * 100, "%)")) %>% 
  select(population, trait, heterogeneity, lackCorrelation) %>% 
  gather(grp, value, -trait, -population) %>% 
  spread(grp, value)

# population trait           heterogeneity  lackCorrelation 
# 1 all        GrainYield      22853.898 (9%) 233134.396 (91%)
# 2 all        HeadingDate     1.057 (15%)    6.127 (85%)     
# 3 all        HeadingDateAGDD 1301.759 (15%) 7646.485 (85%)  
# 4 all        PlantHeight     2.034 (11%)    16.958 (89%)    
# 5 tp         GrainYield      19983.738 (8%) 239570.652 (92%)
# 6 tp         HeadingDate     1.271 (21%)    4.78 (79%)      
# 7 tp         HeadingDateAGDD 1444.957 (19%) 6300.711 (81%)  
# 8 tp         PlantHeight     2 (11%)        16.474 (89%)    
# 9 vp         GrainYield      19347.096 (9%) 190971.903 (91%)
# 10 vp         HeadingDate     1.181 (22%)    4.135 (78%)     
# 11 vp         HeadingDateAGDD 1205.73 (18%)  5435.906 (82%)  
# 12 vp         PlantHeight     2.41 (13%)     15.976 (87%)


# Plot
g_varGE_comp <- varGE_components %>% 
  select(population, trait, GeneticHeterogen = V_prop, LackCorrelation = L_prop) %>% 
  gather(group, proportion, -population, -trait) %>% 
  mutate(group = factor(group, levels = rev(unique(group))),
         population = str_to_upper(population)) %>%
  ggplot(data = ., aes(x = trait, y = proportion, fill = group)) + 
  geom_col() +
  geom_text(aes(y = proportion / 2, label = round(proportion, 2))) +
  scale_fill_brewer(name = NULL, palette = "Set2") +
  facet_grid(~ population) +
  ylab(expression("Proportion of"~sigma[GE]^2)) +
  xlab("Trait") +
  theme_acs() +
  theme(legend.position = c(0.87, 0.75))

# ggsave(filename = "varGE_components.jpg", plot = g_varGE_comp, width = 8, height = 3, path = fig_dir, dpi = 1000)
ggsave(filename = "varGE_components_withAGDD.jpg", plot = g_varGE_comp, width = 8, height = 3, path = fig_dir, dpi = 1000)



## Combine the plots for all variance components with the GxE variance components
# Plot
g_varGE_comp_tp <- varGE_components %>% 
  filter(population == "tp") %>%
  select(population, trait, GeneticHeterogen = V_prop, LackCorrelation = L_prop) %>% 
  gather(group, proportion, -population, -trait) %>% 
  mutate(group = factor(group, levels = rev(unique(group))),
         population = str_to_upper(population)) %>%
  ggplot(data = ., aes(x = trait, y = proportion, fill = group)) + 
  geom_col() +
  # geom_text(aes(y = proportion / 2, label = round(proportion, 2))) +
  scale_fill_manual(name = NULL, values = c(neyhart_palette("umn1")[3], neyhart_palette("umn2", 8)[8])) +
  ylab(expression("Proportion of"~sigma[GE]^2)) +
  xlab("Trait") +
  theme_acs() +
  theme(legend.position = "bottom")

g_varcomp_combine <- plot_grid(g_varprop_tp + theme(legend.position = "top"), 
                               g_varGE_comp_tp, ncol = 1, rel_heights = c(1, 0.4), axis = "lr", align = "v")

ggsave(filename = "variance_components_combined.jpg", plot = g_varcomp_combine, path = fig_dir,
       height = 8, width = 5, dpi = 1000)


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
var_comp_table <- select(stage_two_fits_GYL_varprop1, trait, population, source, variance, proportion = var_prop, p.value) %>%
  bind_rows(., varGE_components1) %>%
  ungroup() %>%
  mutate(significance = case_when(p.value < 0.001 ~ "***",
                                  p.value < 0.01 ~ "**",
                                  p.value < 0.05 ~ "*",
                                  TRUE ~ ""),
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


write_csv(x = var_comp_table, path = file.path(fig_dir, "population_variance_components_decomposed_withAGDD.csv"))

# Remove AGDD
write_csv(x = filter(var_comp_table, Trait %in% str_add_space(traits)), path = file.path(fig_dir, "population_variance_components_decomposed.csv"))


## Just the TP
var_comp_table %>% 
  select(Trait, Source, TP) %>% 
  spread(Trait, TP) %>% 
  write_csv(x = ., path = file.path(fig_dir, "tp_population_variance_components_decomposed_withAGDD.csv"))

var_comp_table %>% 
  filter(Trait %in% str_add_space(traits)) %>%
  select(Trait, Source, TP) %>% 
  spread(Trait, TP) %>% 
  write_csv(x = ., path = file.path(fig_dir, "tp_population_variance_components_decomposed.csv"))




##


