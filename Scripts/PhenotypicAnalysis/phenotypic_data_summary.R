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
library(msm)



## significance level
alpha <- 0.05



## Determine the number of environments per trait
S2_MET_BLUEs %>% 
  group_by(trait) %>% 
  summarize_at(vars(location, year, environment), n_distinct)

# trait        location  year environment
# 1 GrainProtein        7     3          10
# 2 GrainYield         13     3          27
# 3 HeadingDate        11     3          24
# 4 PlantHeight        12     3          27
# 5 TestWeight          7     2          12



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
# 1 GrainProtein 0.295 0.840 0.647
# 2 GrainYield   0.179 0.861 0.541
# 3 HeadingDate  0.656 0.980 0.863
# 4 PlantHeight  0.160 0.882 0.530
# 5 TestWeight   0.501 0.945 0.689


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
  group_by(trait, population) %>%
  do({

    df <- droplevels(.)
    # Edit factors
    df1 <- df %>%
      droplevels() %>%
      mutate_at(vars(location, year, line_name), fct_contr_sum) %>%
      mutate(weights = std_error^2)


    ## fit the full model
    # fit_gly <- lmer(formula = value ~ 1 + location + year + location:year + (1|line_name) + (1|line_name:location) +
    #               (1|line_name:year) + (1|line_name:location:year), data = df, control = lmer_control)
    #
    # fit_ge <- lmer(formula = value ~ 1 + environment + (1|line_name) + (1|line_name:environment),
    #                data = df, control = lmer_control)


    # Random models
    fit_gly <- lmer(formula = value ~ 1 + (1|location) + (1|year) + (1|location:year) + (1|line_name) + (1|line_name:location) +
                      (1|line_name:year) + (1|line_name:location:year), data = df1, control = lmer_control)

    fit_ge <- lmer(formula = value ~ 1 + (1|environment) + (1|line_name) + (1|line_name:environment),
                   data = df1, control = lmer_control)


    # fit_mmer <- mmer(value ~ 1, random = ~ line_name + location + year, data = df1)


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
     
    
    # VCOV matrix of variance components
    vcov_rand_gly <- fit_gly@vcov_varpar
    vcov_rand_ge <- fit_ge@vcov_varpar


    ## Get variance components
    var_comp_gly <- as.data.frame(VarCorr(fit_gly)) %>%
      select(term = grp, variance = vcov) %>%
      mutate(std_error = sqrt(diag(vcov_rand_gly)))
    var_comp_ge <- as.data.frame(VarCorr(fit_ge)) %>%
      select(term = grp, variance = vcov) %>%
      mutate(std_error = sqrt(diag(vcov_rand_ge)))
      

    # Return data_frame
    tibble(model = c("gly", "ge"), fit = list(fit_gly, fit_ge),
           h2 = c(herit2(fit_gly), herit2(fit_ge)), var_comp = list(var_comp_gly, var_comp_ge))

  }) %>% ungroup()



################################
# Use ASREML to fit the models
# 
# 
# Import the results
load(file.path(result_dir, "stage_two_fits_asreml.RData"))


################################




## Create a table for output
stage_two_fits_asreml_GYL_varcomp <- stage_two_fits_asreml %>%
  mutate(var_comp = map(var_comp, ~rownames_to_column(., "term"))) %>%
  unnest(var_comp) %>%
  rename(variance = component) %>%
  filter(model == "gly") %>%
  mutate(term = ifelse(term == "units!R", "Residual", term))

## Create a table for output
stage_two_fits_lmer_GYL_varcomp <- stage_two_fits %>%
  unnest(var_comp) %>%
  filter(model == "gly")

## Choose the data.frame to use
stage_two_fits_GYL_varcomp <- stage_two_fits_lmer_GYL_varcomp
# stage_two_fits_GYL_varcomp <- stage_two_fits_asreml_GYL_varcomp




## Plot the proportion of variance from each source
stage_two_fits_GYL_varprop <- stage_two_fits_GYL_varcomp %>% 
  mutate(term = map(term, ~str_split(., pattern = ":", simplify = FALSE) %>% map(str_to_title) %>% .[[1]]) %>% 
           map_chr(~paste(., collapse = " x ")),
         term = str_replace_all(term, "Line_name", "Genotype"),
         term = factor(term, levels = c("Genotype", "Location", "Year", "Location x Year", 
                                            "Genotype x Location", "Genotype x Year",
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
  ggplot(aes(x = trait, y = var_prop, fill = term)) + 
  geom_col() +
  ylab("Proportion of variance") +
  scale_fill_manual(values = var_comp_colors, name = NULL) +
  facet_grid(~ population) + 
  theme_presentation2(base_size = 10) +
  theme(axis.title.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))

# Save
ggsave(filename = "variance_components_expanded.jpg", plot = g_varprop, path = fig_dir, width = 9, height = 6, dpi = 1000)




#######
#######
#######



## Look at the heritability and plot
g_herit <- stage_two_fits_GYL_varcomp %>% 
  distinct(trait, population, h2) %>%
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



## Output a table of all variance components for all populations
prop_varcomp1 <- stage_two_fits_GYL_varcomp %>% 
  select(trait, population, term, variance) %>% 
  group_by(trait, population) %>% 
  mutate(proportion = variance / sum(variance),
         term = str_replace_all(term, ":", " x "), 
         term = str_replace_all(term, "line_name", "Genotype"), 
         term = str_to_title(term))


# Component order
comp_order <- unique(stage_two_fits_GYL_varprop1$term) %>% 
  {.[order(str_count(., "x"))]} %>%
  {c(str_subset(., "[^Residual]"), "Residual")}

## Combine and output
var_comp_table <- select(stage_two_fits_GYL_varprop1, trait, population, term, variance, proportion = var_prop) %>%
  # bind_rows(., varGE_components1) %>%
  ungroup() %>%
  mutate(significance = "",
         # significance = case_when(p.value < 0.001 ~ "***",
         #                          p.value < 0.01 ~ "**",
         #                          p.value < 0.05 ~ "*",
         #                          TRUE ~ ""),
         variance = signif(variance, 3) %>% formatC(x = ., digits = 3, format = "f") %>%
           str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         proportion = signif(proportion * 100, 3) %>% formatC(x = ., digits = 3, format = "f") %>%
           str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         annotation = paste0(variance, significance, " (", proportion, "%)"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, term, annotation) %>%
  mutate(term = str_add_space(term),
         trait = str_add_space(trait),
         term = factor(term, levels = comp_order),
         population = toupper(population)) %>% 
  arrange(population, term) %>%
  rename_all(str_to_title) %>%
  spread(Population, Annotation)

write_csv(x = var_comp_table, path = file.path(fig_dir, "population_variance_components_decomposed.csv"))


## Combine and output
var_comp_table1 <- select(stage_two_fits_GYL_varprop1, trait, population, term, variance, std_error) %>%
  ungroup() %>%
  mutate(significance = "",
         # significance = case_when(p.value < 0.001 ~ "***",
         #                          p.value < 0.01 ~ "**",
         #                          p.value < 0.05 ~ "*",
         #                          TRUE ~ ""),
         variance = signif(variance, 3) %>% formatC(x = ., digits = 3, format = "f") %>%
           str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         std_error = signif(std_error, 3) %>% formatC(x = ., digits = 3, format = "f"),
         annotation = paste0(variance, significance, " (", std_error, ")"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, term, annotation) %>%
  mutate(term = str_add_space(term),
         trait = str_add_space(trait),
         term = factor(term, levels = comp_order),
         population = toupper(population)) %>% 
  arrange(population, term) %>%
  rename_all(str_to_title) %>%
  spread(Population, Annotation)


write_csv(x = var_comp_table1, path = file.path(fig_dir, "population_variance_components_decomposed1.csv"))




####

####






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
mixed_model_fit <- S2_MET_BLUEs_tomodel %>%
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
    
    # ## fit a mixed model to predict all genotype-environment means ##
    # # MF
    # mf <- model.frame(value ~ line_name + environment, df1)
    # y <- model.response(mf)
    # X <- model.matrix(~ 1 + line_name + environment, data = mf)
    # Z <- model.matrix(~ -1 + line_name:environment, mf)
    # 
    # fit1 <- mixed.solve(y = y, Z = Z, K = K_ge, X = X, method = "REML")
    
    ## Fit the model using sommer
    fit1 <- mmer(fixed = value ~ line_name + environment, random = ~ vs(line_name:environment, Gu = K_ge), data = df1,
                 date.warning = FALSE, weights = std_error)
    
    ## Calculate degrees of freedom
    df_fit1 <- nrow(df1) - n_distinct(df1$line_name) - 1 - n_distinct(df1$environment) - 1 + 1
    # Extract the residual variance
    varR <- c(fit1$sigma$units)
    
    
    ## Pull out effects
    # grand_mean <- fit1$beta[1]
    grand_mean <- subset(fit1$Beta, Effect == "(Intercept)", Estimate, drop = TRUE)
    
    # g_effects <- fit1$beta %>% subset(., str_detect(names(.), "line_name"))
    # e_effects <- fit1$beta %>% subset(., str_detect(names(.), "environment"))
    # # Add last level
    # g_effects <- c(g_effects, -sum(g_effects))
    # e_effects <- c(e_effects, -sum(e_effects))
    
    g_effects_df <- subset(fit1$Beta, str_detect(Effect, "line_name")) %>%
      select(line_name = Effect, effect = Estimate) %>%
      mutate(line_name = str_remove_all(line_name, "line_name")) %>%
      add_row(line_name = last(levels(df1$line_name)), effect = -sum(.$effect))
    e_effects_df <- subset(fit1$Beta, str_detect(Effect, "environment")) %>%
      select(environment = Effect, effect = Estimate) %>%
      mutate(environment = str_remove_all(environment, "environment")) %>%
      add_row(environment = last(levels(df1$environment)), effect = -sum(.$effect))
    
 
    ## Random effects
    # rand_eff_df <- fit1$u %>% 
    rand_eff_df <- fit1$U$`u:line_name:environment`$value %>% 
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
    
    ## Output a tibble
    tibble(y_hat = list(y_hat_df), model_acc = acc, df = df_fit1, varR = varR)
    
  }) %>% ungroup()

## Calculate accuracy and RMSE of model fit
mixed_model_fit %>%
  unnest(y_hat) %>%
  left_join(., filter(S2_MET_BLUEs_tomodel, population == "tp")) %>%
  group_by(trait) %>%
  summarize(acc = cor(y_hat, value, use = "complete.obs"),
            rmse = sqrt(mean((y_hat - value)^2, na.rm = TRUE)))


## Create a new df for fitting the AMMI model
pheno_tomodel_ammi <- mixed_model_fit %>%
  unnest(y_hat) %>%
  left_join(., filter(S2_MET_BLUEs_tomodel, population == "tp")) %>%
  select(trial, environment, location, year, trait, line_name, y_hat, value, std_error)



## Fit ammi models using the predictions from the first model:
ammi_fit <- pheno_tomodel_ammi %>%
  group_by(trait) %>%
  do({
    df <- .
    # Set factors
    df1 <- droplevels(df) %>%
      mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.)))
    
    ## Create an incomplete two-way table
    two_way_missing <- df1 %>% 
      select(line_name, environment, value) %>% 
      as.data.frame() %>%
      spread(environment, value) %>% 
      column_to_rownames("line_name") %>%
      as.matrix()
    # Imput
    two_way_complete <- em(Y = two_way_missing, model = "AMMI", maxiter = 200, k = 1)
    
    # Convert back into df
    df2 <- two_way_complete %>% 
      as.data.frame() %>% 
      rownames_to_column("line_name") %>% 
      gather(environment, value_imputed, -line_name) %>%
      left_join(df1, .) %>%
      select(trial, environment, location, year, trait, line_name, mixed_model = y_hat, imputed = value_imputed) %>%
      gather(type, value, mixed_model, imputed)
    
    ## Split and fit the AMMI models
    ammi_fit_list <- df2 %>%
      split(., .$type) %>%
      map(~{
        
        dat <- .x
        # Fit the AMMI model
        fit_ammi <- bilinear(x = dat, G = "line_name", E = "environment", y = "value", 
                             model = "AMMI", alpha = alpha, test = "Ftest")
        
        ## Calculate scores
        svd <- fit_ammi$svdE
        # Get the singular values
        sing_val <- svd$d
        
        # Calculate scores
        g_scores <- svd$u %*% diag(sqrt(sing_val))
        dimnames(g_scores) <- list(names(fit_ammi$Geffect), paste0("PC", seq_len(ncol(g_scores))))
        e_scores <- svd$v %*% diag(sqrt(sing_val))
        dimnames(e_scores) <- list(names(fit_ammi$Eeffect), paste0("PC", seq_len(ncol(e_scores))))
        
        ## Add these scores to the fit_ammi object
        fit_ammi$scores <- list(Gscores = g_scores, Escores = e_scores)
        
        # Output this object
        fit_ammi
        
      })
    
    ## Output tibble
    tibble(type = names(ammi_fit_list), fit_ammi = ammi_fit_list)
    
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

## Stick with the imputed type
ammi_scores <- filter(ammi_scores, type == "imputed")


## Display proportion of variance explained using normalized eigenvalues
(g_pc_ss <- ammi_scores %>%
  unnest(varprop) %>%
  # Plot
  ggplot(aes(x = PC_num, y = propvar_SS)) +
  geom_point() + 
  scale_x_continuous(name = "PC", breaks = pretty) +
  facet_grid(type ~ trait, scales = "free_x") +
  theme_light() )
# Save
ggsave(filename = "ammi_propvar_SS.jpg", plot = g_pc_ss, path = fig_dir, height = 4, width = 6, 
       dpi = 1000)

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
                 nPC = which.min(propvar_diff >= tol)))

## Summary df of number of sig PCs
ammi_sig_PCs_summ <- ammi_sig_PCs %>%
  group_by(trait) %>%
  filter(PC_num %in% seq(1, unique(nPC))) %>%
  summarize(total_propvar = sum(propvar), nPC = unique(nPC))

## Output a table
write_csv(x = ammi_sig_PCs_summ, path = file.path(fig_dir, "sig_ammi_PC_propvar.csv"))


## Fit a model for that number of PCS
ammiN_fit <- ammi_fit %>%
  filter(type == "imputed") %>%
  left_join(., ammi_sig_PCs_summ, by = "trait") %>%
  group_by(trait) %>%
  do({
    
    row <- .
    nPC <- row$nPC
    # Get the fitted ammi model
    fitted_ammi <- row$fit_ammi[[1]]
    
    # print(row$trait)
    
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
save("mixed_model_fit", "ammi_fit", "ammiN_fit", "ammiN_fit_location", "ammi_gplot_list", 
     file = file.path(result_dir, "ammi_model_fit.RData"))

#





































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




 