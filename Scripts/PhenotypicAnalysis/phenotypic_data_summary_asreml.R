## S2MET Phenotypic Adjustment
## 
## Author: Jeff Neyhart
## Last updated: 10 June 2019
## 
## This notebook outlines procedures for calculating adjusted phenotypic means of
## entries in trials that belong to the `S2MET` experiment.
## 

# Run on a local machine
load(file.path(repo_dir, "Data/asreml_model_environment.RData"))

setwd("G:/AGRO/BARLEY_LAB/Jeff/Projects/S2MET_Predictions_Models/")
repo_dir <- getwd()
source(file.path(repo_dir, "source_functions.R"))

library(tibble)
library(modelr)
library(purrr)
library(dplyr)
library(tidyr)
library(asreml)

## Stage-Two analysis


## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs, line_name %in% vp), population = "vp")) %>%
  filter(environment %in% tp_vp_env) %>%
  mutate_at(vars(location, year, line_name), as.factor)


# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP


# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_asreml <- S2_MET_BLUEs_tomodel %>%
  group_by(trait, population) %>%
  do({
    
    df <- droplevels(.)
    # Edit factors
    df1 <- df %>%
      droplevels() %>%
      mutate_at(vars(location, year, line_name, environment), as.factor) %>%
      mutate_at(vars(location, year, line_name, environment), fct_contr_sum) %>%
      mutate(weights = std_error^2)
    
    
    ## fit the full model
    # fit_gly <- lmer(formula = value ~ 1 + location + year + location:year + (1|line_name) + (1|line_name:location) +
    #               (1|line_name:year) + (1|line_name:location:year), data = df, control = lmer_control)
    # 
    # fit_ge <- lmer(formula = value ~ 1 + environment + (1|line_name) + (1|line_name:environment), 
    #                data = df, control = lmer_control)
    
    
    # Fit a model
    fit_gly_asreml <- asreml(fixed = value ~ 1, 
                             random = ~ line_name + location + year + location:year + line_name:location +
                               line_name:year + line_name:location:year,
                             residual = ~ units,
                             data = df1, 
                             weights = weights,
                             fail = "soft",
                             workspace = "1gb")
    
    fit_ge_asreml <- asreml(fixed = value ~ 1, 
                             random = ~ line_name + environment + line_name:environment,
                             residual = ~ units,
                             data = df1, 
                             weights = weights,
                             fail = "soft",
                             workspace = "1gb")
    
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
    
    ## Calculate heritability ##
    # Genetic variance component
    fit_h2_list <- lapply(X = list(gly = fit_gly_asreml, ge = fit_ge_asreml), FUN = function(fit) {
      
      fit_summary <- summary(fit)
      vc.g <- fit_summary$varcomp["line_name","component"]
      
      # Mean variance of a difference of two genotypic BLUPs
      vdBLUP.mat <- predict(fit, classify="line_name",
                            sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
      vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
      h2 <- 1 - (vdBLUP.avg / 2 / vc.g)
      
      list(h2 = h2, var_comp = fit_summary$varcomp)
      
    })
    
    
    # Return data_frame
    tibble(model = c("gly", "ge"),
           fit = list(fit_gly_asreml, fit_ge_asreml), 
           h2 = sapply(fit_h2_list, "[[", "h2"),
           var_comp = lapply(fit_h2_list, "[[", "var_comp"))
    
  }) %>% ungroup()


## Save
save("stage_two_fits", file = file.path(repo_dir, "Results/stage_two_fits_asreml.RData"))



