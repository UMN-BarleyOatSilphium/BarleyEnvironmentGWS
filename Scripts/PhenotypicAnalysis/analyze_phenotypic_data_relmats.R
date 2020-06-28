## S2MET Prediction Models
## 
## Fit full models to estimate variance components
## 
## Author: Jeff Neyhart
## Last modified: 8 January 2020
## 

# This script will analyze the phenotypic data using covariates and markers
# 

# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions_Models/"
source(file.path(repo_dir, "source_MSI.R"))



# Other packages
library(modelr)
library(broom)
# library(lme4qtl)

load(file = file.path(result_dir, "feature_selection_results.RData"))
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

# Set the number of cores
n_cores <- 8

# time_frame to use for the location relationship matrix
time_frame_use <- "time_frame5_2010_2014"

# Simple function to return input.
f_return <- function(x) x


# Environment-specific ----------------------------------------------------



## Data.frame to use to fit models
data_to_model <- S2_MET_BLUEs %>% 
  filter(environment %in% train_test_env, line_name %in% c(tp_geno, vp_geno)) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter(!(population == "tp" & line_name %in% vp_geno), !(population == "vp" & line_name %in% tp_geno)) %>%
  # Copy columns for modeling covariates
  droplevels() %>%
  mutate_at(vars(line_name, environment), ~fct_contr_sum(as.factor(.))) %>%
  mutate(gxe = interaction(line_name, environment, drop = TRUE, sep = ":")) %>%
  mutate_at(vars(line_name, environment, gxe), list(cov = f_return))

## Assemble covariates to model
covariates_to_model <- concurrent_fact_reg_feature_selection %>% 
  select(trait, model, source, apriori) %>% 
  gather(feat_sel_type, features, apriori) %>%
  bind_rows(., rename(concurrent_feature_selection, features = adhoc)) %>%
  mutate(features = map(features, "optVariables") %>% map(~setdiff(x = ., y = "line_name"))) %>%
  filter(source == "daymet", model == "model3") %>%
  select(trait, feat_sel_type, features) %>%
  bind_rows(., tibble(trait = traits, feat_sel_type = "all", features = list(names(ec_tomodel_scaled$daymet)[-1:-2])))


# Matrix of covariates
covariate_mat <- ec_tomodel_scaled$daymet %>%
  as.data.frame() %>%
  select(-source) %>%
  column_to_rownames("environment") %>%
  as.matrix()
  

## Create a results df
data_to_model_split <- data_to_model %>% 
  group_by(trait, population) %>% 
  nest() %>%
  ungroup() %>%
  full_join(., covariates_to_model, by = "trait") %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE)

## Iterate over cores
pheno_variance_analysis_out <- coreApply(X = data_to_model_split, FUN = function(core_df) {
  
  # Output list
  out <- vector("list", length = nrow(core_df))
  
  # Iterate over rows
  for (i in seq_along(out)) {
  
    row <- core_df[i,]
    
    # Trait
    tr <- row$trait
    data <- droplevels(row$data[[1]])
    genotypes <- levels(data$line_name)
    environments <- levels(data$environment)
    covariates <- row$features[[1]]
    
    # Separate covariates into main effect or interaction
    covariate_list <- tibble(term = unique(covariates)) %>% 
      mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
             covariate = str_remove(term, "line_name:")) %>%
      split(.$class) %>% 
      map("covariate")
    
    
    # Create or subset relationship matrices
    K_use <- K[genotypes, genotypes]
    Zg <- model.matrix(~ -1 + line_name, data)
    Emain <- Eint2 <- Env_mat(x = covariate_mat[environments, covariate_list$main, drop = FALSE], method = "Jarq")
    Eint <- Env_mat(x = covariate_mat[environments, covariate_list$int, drop = FALSE], method = "Jarq")
    # Replace with diagonal matrix if NA
    if (all(is.na(Eint))) { 
      Eint[] <- 0
      diag(Eint) <- 1 
      no_gxe_cov <- TRUE
      
    } else {
      no_gxe_cov <- FALSE
      
    }
    
    Ze <- model.matrix(~ -1 + environment, data)
    KE <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% Eint, Ze)
    dimnames(KE) <- replicate(2, data$gxe, simplify = FALSE)
    
    # # Fit the model
    # fit <- relmatLmer(value ~ (1|line_name_cov) + (1|line_name) + (1|environment_cov) + (1|environment) + 
    #                     (1|gxe), 
    #                   data = data, relmat = list(line_name_cov = K_use, environment_cov = Emain))
    
    ## Sommer
    model_stdout <- capture.output({ 
      model_try <- try( {
        fit_mmer <- mmer(fixed = value ~ 1, 
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + environment + vs(environment_cov, Gu = Emain) +
                           gxe + vs(gxe_cov, Gu = KE),
                         data = data, verbose = TRUE) 
      }, silent = TRUE) })
    
    # Retry with fewer iterations if singular
    if (class(model_try) == "try-error") {
      
      # Parse the monitor for iterations
      model_monitor <- read_table(head(model_stdout, -1))
      
      # Select 1 iteration fewer
      iter_max <- max(model_monitor$iteration) - 1
      
      # If only one iteration, then error out
      if (iter_max == 0) {
        varcomp_df <- tibble(term = "Error in singularlity.")
        
      } else {
      
        # refit the model
        fit_mmer <- mmer(fixed = value ~ 1, 
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + environment + vs(environment_cov, Gu = Emain) +
                           gxe + vs(gxe_cov, Gu = KE),
                         data = data, verbose = TRUE, iter = iter_max)
        
        # Get the variance components; asemble into a tibble
        varcomp_df <- fit_mmer$sigma %>% 
          map_dbl(~.) %>% 
          tibble(term = names(.), variance = .) %>% 
          mutate(term = str_remove_all(term, "u:"), 
                 source = str_remove_all(term, "_cov")) %>%
          group_by(source) %>%
          mutate(total_source_variance = sum(variance)) %>%
          ungroup() %>%
          mutate(variance_prop = variance / total_source_variance,
                 note = list(NULL))
        
      }
      
      # Give an error message if empty
    } else if (is_empty(fit_mmer)) {
      varcomp_df <- tibble(term = "Error in singularlity.")
      
      
      # else proceed normally
    } else {
      
      # Get the variance components; asemble into a tibble
      varcomp_df <- fit_mmer$sigma %>% 
        map_dbl(~.) %>% 
        tibble(term = names(.), variance = .) %>% 
        mutate(term = str_remove_all(term, "u:"), 
               source = str_remove_all(term, "_cov")) %>%
        group_by(source) %>%
        mutate(total_source_variance = sum(variance)) %>%
        ungroup() %>%
        mutate(variance_prop = variance / total_source_variance,
               note = list(NULL))
      
    }
    
    # Return
    out[[i]] <- varcomp_df
    
  } # Close loop
  
  # Add out to core_df and return
  core_df %>%
    mutate(results = out) %>%
    select(-data, -core)
  
})


# Rename
environment_pheno_variance_analysis <- pheno_variance_analysis_out



## If working locally, load the results to complete unfitted models
if (startsWith(repo_dir, "C:/")) {
  pheno_variance_analysis1 <- bind_rows(environment_pheno_variance_analysis)
  
  # Find the results with 1 row- these are errors
  which_errors <- map_lgl(pheno_variance_analysis1$results, ~nrow(.) == 1)
  
  # If any, refit them
  if (any(which_errors)) {
    
    pheno_variance_analysis_refit <- pheno_variance_analysis1[which_errors,,drop = FALSE] %>%
      left_join(., select(bind_rows(data_to_model_split), -features)) %>%
      mutate(results = list(NULL))
    
    # Iterate over rows
    for (i in seq_len(nrow(pheno_variance_analysis_refit))) {
      
      row <- pheno_variance_analysis_refit[i,]
      
      # Trait
      tr <- row$trait
      data <- droplevels(row$data[[1]])
      genotypes <- levels(data$line_name)
      environments <- levels(data$environment)
      covariates <- row$features[[1]]
      
      # Separate covariates into main effect or interaction
      covariate_list <- tibble(term = unique(covariates)) %>%
        mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
               covariate = str_remove(term, "line_name:")) %>%
        split(.$class) %>%
        map("covariate")
      
      
      # Create or subset relationship matrices
      K_use <- K[genotypes, genotypes]
      Zg <- model.matrix(~ -1 + line_name, data)
      Emain <- Eint2 <- Env_mat(x = covariate_mat[environments, covariate_list$main, drop = FALSE], method = "Jarq")
      Eint <- Env_mat(x = covariate_mat[environments, covariate_list$int, drop = FALSE], method = "Jarq")
      # Replace with diagonal matrix if NA
      if (all(is.na(Eint))) {
        Eint[] <- 0
        diag(Eint) <- 1
        no_gxe_cov <- TRUE
        
      } else {
        no_gxe_cov <- FALSE
        
      }
      
      Ze <- model.matrix(~ -1 + environment, data)
      KE <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% Eint, Ze)
      dimnames(KE) <- replicate(2, data$gxe, simplify = FALSE)
      
      ## Sommer
      model_stdout <- capture.output({
        model_try <- try( {
          fit_mmer <- mmer(fixed = value ~ 1,
                           random = ~ line_name + vs(line_name_cov, Gu = K_use) + environment + vs(environment_cov, Gu = Emain) +
                             gxe + vs(gxe_cov, Gu = KE),
                           data = data, verbose = TRUE)
        }, silent = TRUE) })
      
      
      
      # Retry with fewer iterations if singular
      if (class(model_try) == "try-error" | any(grepl(pattern = "singular", x = model_stdout))) {
        
        # Parse the monitor for iterations
        model_monitor <- read_table(head(model_stdout, -1))
        
        # Select 1 iteration fewer
        iter_max <- max(model_monitor$iteration) - 1
        
          
        # refit the model
        fit_mmer <- mmer(fixed = value ~ 1, 
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + environment + vs(environment_cov, Gu = Emain) +
                           gxe + vs(gxe_cov, Gu = KE),
                         data = data, verbose = TRUE, iter = iter_max)
        
      }
          
      
      # Get the variance components; asemble into a tibble
      varcomp_df <- fit_mmer$sigma %>%
        map_dbl(~.) %>%
        tibble(term = names(.), variance = .) %>%
        mutate(term = str_remove_all(term, "u:"),
               source = str_remove_all(term, "_cov")) %>%
        group_by(source) %>%
        mutate(total_source_variance = sum(variance)) %>%
        ungroup() %>%
        mutate(variance_prop = variance / total_source_variance,
               note = list(NULL))
      
      # Add to the df
      pheno_variance_analysis_refit$results[[i]] <- varcomp_df
      
    } # end the loop
    
    # Add the results to the original df
    pheno_variance_analysis2 <- bind_rows(
      pheno_variance_analysis1[!which_errors,,drop = FALSE],
      pheno_variance_analysis_refit
    ) %>% select(-data, -core)
    
    
  }
  
  environment_pheno_variance_analysis <- pheno_variance_analysis2
  
}
















# Location-specific -------------------------------------------------------


# Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

## Pull out location covariates to use
loc_ec_tomodel_scaled <-  historical_ec_tomodel_timeframe_scaled %>%
  subset(., grepl(pattern = time_frame_use, x = names(.)))


data_to_model <- S2_MET_loc_BLUEs %>% 
  filter(location %in% train_test_loc, line_name %in% c(tp_geno, vp_geno)) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter(!(population == "tp" & line_name %in% vp_geno), !(population == "vp" & line_name %in% tp_geno)) %>%
  # Copy columns for modeling covariates
  droplevels() %>%
  mutate_at(vars(line_name, location), ~fct_contr_sum(as.factor(.))) %>%
  mutate(gxl = interaction(line_name, location, drop = TRUE, sep = ":")) %>%
  mutate_at(vars(line_name, location, gxl), list(cov = f_return))

# Convert model 2/3 to 4/5
historical_feature_selection <- historical_feature_selection %>% 
  mutate(model = str_replace_all(model, c("model2" = "model4", "model3" = "model5")))

## Assemble covariates to model
covariates_to_model <- historical_fact_reg_feature_selection %>% 
  select(trait, model, source, apriori) %>% 
  gather(feat_sel_type, features, apriori) %>%
  bind_rows(., rename(historical_feature_selection, features = adhoc)) %>%
  mutate(features = map(features, "optVariables") %>% map(~setdiff(x = ., y = "line_name"))) %>%
  filter(source == "daymet", model == "model5") %>%
  select(trait, feat_sel_type, features) %>%
  bind_rows(., tibble(trait = traits, 
                      feat_sel_type = "all", features = list(names(loc_ec_tomodel_scaled$daymet.time_frame5_2010_2014)[-1:-3])))


# Matrix of covariates
covariate_mat <- loc_ec_tomodel_scaled$daymet.time_frame5_2010_2014 %>%
  as.data.frame() %>%
  select(-source, -time_frame) %>%
  column_to_rownames("location") %>%
  as.matrix()


## Create a results df
data_to_model_split <- data_to_model %>% 
  group_by(trait, population) %>% 
  nest() %>%
  ungroup() %>%
  full_join(., covariates_to_model, by = "trait") %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE)

## Iterate over cores
pheno_variance_analysis_out <- coreApply(X = data_to_model_split, FUN = function(core_df) {
  
  # Output list
  out <- vector("list", length = nrow(core_df))
  
  # Iterate over rows
  for (i in seq_along(out)) {
    
    row <- core_df[i,]
    
    # Trait
    tr <- row$trait
    data <- droplevels(row$data[[1]])
    genotypes <- levels(data$line_name)
    locations <- levels(data$location)
    covariates <- row$features[[1]]
    
    # Separate covariates into main effect or interaction
    covariate_list <- tibble(term = unique(covariates)) %>% 
      mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
             covariate = str_remove(term, "line_name:")) %>%
      split(.$class) %>% 
      map("covariate")
    
    
    # Create or subset relationship matrices
    K_use <- K[genotypes, genotypes]
    Zg <- model.matrix(~ -1 + line_name, data)
    Emain <- Eint2 <- Env_mat(x = covariate_mat[locations, covariate_list$main, drop = FALSE], method = "Jarq")
    Eint <- Env_mat(x = covariate_mat[locations, covariate_list$int, drop = FALSE], method = "Jarq")
    # Replace with diagonal matrix if NA
    if (all(is.na(Eint))) { 
      Eint[] <- 0
      diag(Eint) <- 1 
      no_gxe_cov <- TRUE
      
    } else {
      no_gxe_cov <- FALSE
      
    }
    
    Ze <- model.matrix(~ -1 + location, data)
    KE <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% Eint, Ze)
    dimnames(KE) <- replicate(2, data$gxl, simplify = FALSE)
    
    ## Sommer
    model_stdout <- capture.output({ 
      model_try <- try( {
        fit_mmer <- mmer(fixed = value ~ 1, 
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + location + vs(location_cov, Gu = Emain) +
                           gxl + vs(gxl_cov, Gu = KE),
                         data = data, verbose = TRUE) 
      }, silent = TRUE) })
    
    # Retry with fewer iterations if singular
    if (class(model_try) == "try-error") {
      
      # Parse the monitor for iterations
      model_monitor <- read_table(head(model_stdout, -1))
      
      # Select 1 iteration fewer
      iter_max <- max(model_monitor$iteration) - 1
      
      # If only one iteration, then error out
      if (iter_max == 0) {
        varcomp_df <- tibble(term = "Error in singularlity.")
        
      } else {
        
        # refit the model
        fit_mmer <- mmer(fixed = value ~ 1, 
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + location + vs(location_cov, Gu = Emain) +
                           gxl + vs(gxl_cov, Gu = KE),
                         data = data, verbose = TRUE, iter = iter_max)
        
        # Get the variance components; asemble into a tibble
        varcomp_df <- fit_mmer$sigma %>% 
          map_dbl(~.) %>% 
          tibble(term = names(.), variance = .) %>% 
          mutate(term = str_remove_all(term, "u:"), 
                 source = str_remove_all(term, "_cov")) %>%
          group_by(source) %>%
          mutate(total_source_variance = sum(variance)) %>%
          ungroup() %>%
          mutate(variance_prop = variance / total_source_variance,
                 note = list(NULL))
        
      }
      
      # Give an error message if empty
    } else if (is_empty(fit_mmer)) {
      varcomp_df <- tibble(term = "Error in singularlity.")
      
      
      # else proceed normally
    } else {
      
      # Get the variance components; asemble into a tibble
      varcomp_df <- fit_mmer$sigma %>% 
        map_dbl(~.) %>% 
        tibble(term = names(.), variance = .) %>% 
        mutate(term = str_remove_all(term, "u:"), 
               source = str_remove_all(term, "_cov")) %>%
        group_by(source) %>%
        mutate(total_source_variance = sum(variance)) %>%
        ungroup() %>%
        mutate(variance_prop = variance / total_source_variance,
               note = list(NULL))
      
    }
    
    # Return
    out[[i]] <- varcomp_df
    
  } # Close loop
  
  # Add out to core_df and return
  core_df %>%
    mutate(results = out) %>%
    select(-data, -core)
  
})


# Rename
location_pheno_variance_analysis <- pheno_variance_analysis_out


## If working locally, load the results to complete unfitted models
if (startsWith(repo_dir, "C:/")) {
  pheno_variance_analysis1 <- bind_rows(location_pheno_variance_analysis)
  
  # Find the results with 1 row- these are errors
  which_errors <- map_lgl(pheno_variance_analysis1$results, ~nrow(.) == 1)
  
  # If any, refit them
  if (any(which_errors)) {
    
    pheno_variance_analysis_refit <- pheno_variance_analysis1[which_errors,,drop = FALSE] %>%
      left_join(., select(bind_rows(data_to_model_split), -features)) %>%
      mutate(results = list(NULL))
    
    # Iterate over rows
    for (i in seq_len(nrow(pheno_variance_analysis_refit))) {
      
      row <- pheno_variance_analysis_refit[i,]
      
      # Trait
      tr <- row$trait
      data <- droplevels(row$data[[1]])
      genotypes <- levels(data$line_name)
      locations <- levels(data$location)
      covariates <- row$features[[1]]
      
      # Separate covariates into main effect or interaction
      covariate_list <- tibble(term = unique(covariates)) %>%
        mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
               covariate = str_remove(term, "line_name:")) %>%
        split(.$class) %>%
        map("covariate")
      
      
      # Create or subset relationship matrices
      K_use <- K[genotypes, genotypes]
      Zg <- model.matrix(~ -1 + line_name, data)
      Emain <- Eint2 <- Env_mat(x = covariate_mat[locations, covariate_list$main, drop = FALSE], method = "Jarq")
      Eint <- Env_mat(x = covariate_mat[locations, covariate_list$int, drop = FALSE], method = "Jarq")
      # Replace with diagonal matrix if NA
      if (all(is.na(Eint))) {
        Eint[] <- 0
        diag(Eint) <- 1
        no_gxe_cov <- TRUE
        
      } else {
        no_gxe_cov <- FALSE
        
      }
      
      Ze <- model.matrix(~ -1 + location, data)
      KE <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% Eint, Ze)
      dimnames(KE) <- replicate(2, data$gxl, simplify = FALSE)
      
      ## Sommer
      model_stdout <- capture.output({
        model_try <- try( {
          fit_mmer <- mmer(fixed = value ~ 1,
                           random = ~ line_name + vs(line_name_cov, Gu = K_use) + location + vs(location_cov, Gu = Emain) +
                             gxl + vs(gxl_cov, Gu = KE),
                           data = data, verbose = TRUE)
        }, silent = TRUE) })
      
      
      
      # Retry with fewer iterations if singular
      if (class(model_try) == "try-error" | any(grepl(pattern = "singular", x = model_stdout))) {
        
        # Parse the monitor for iterations
        model_monitor <- read_table(head(model_stdout, -1))
        
        # Select 1 iteration fewer
        iter_max <- max(model_monitor$iteration) - 1
        
        ## Stop if iter_max == 0
        if (iter_max == 0) {
          # Add to the df
          pheno_variance_analysis_refit$results[[i]] <- tibble(term = "Convergency error")
          next
          
        }
        
        
        # refit the model
        fit_mmer <- mmer(fixed = value ~ 1,
                         random = ~ line_name + vs(line_name_cov, Gu = K_use) + location + vs(location_cov, Gu = Emain) +
                           gxl + vs(gxl_cov, Gu = KE),
                         data = data, verbose = TRUE, iter = iter_max)
        
      }
      
      
      # Get the variance components; asemble into a tibble
      varcomp_df <- fit_mmer$sigma %>%
        map_dbl(~.) %>%
        tibble(term = names(.), variance = .) %>%
        mutate(term = str_remove_all(term, "u:"),
               source = str_remove_all(term, "_cov")) %>%
        group_by(source) %>%
        mutate(total_source_variance = sum(variance)) %>%
        ungroup() %>%
        mutate(variance_prop = variance / total_source_variance,
               note = list(NULL))
      
      # Add to the df
      pheno_variance_analysis_refit$results[[i]] <- varcomp_df
      
    } # end the loop
    
    # Add the results to the original df
    pheno_variance_analysis2 <- bind_rows(
      pheno_variance_analysis1[!which_errors,,drop = FALSE],
      pheno_variance_analysis_refit
    ) %>% select(-data, -core)
    
    
  }
  
  location_pheno_variance_analysis <- pheno_variance_analysis2
  
}


# Save the df
save("environment_pheno_variance_analysis", "location_pheno_variance_analysis", 
     file = file.path(result_dir, "full_model_variance_analysis.RData"))
