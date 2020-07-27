## S2MET Prediction Models
## 
## Fit full models to estimate variance components
## 
## Author: Jeff Neyhart
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
# library(BGLR)

load(file = file.path(result_dir, "feature_selection_results.RData"))
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

# Set the number of cores
n_cores <- 8


# Simple function to return input.
f_return <- function(x) x


# Fit models ----------------------------------------------------



## Data.frame to use to fit models
# Environment means
env_data_to_model <- S2_MET_BLUEs %>% 
  filter(environment %in% train_test_env, line_name %in% c(tp_geno, vp_geno)) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter(!(population == "tp" & line_name %in% vp_geno), !(population == "vp" & line_name %in% tp_geno)) %>%
  # Copy columns for modeling covariates
  droplevels() %>%
  # Rename site as environment
  rename(site = environment) %>%
  mutate_at(vars(line_name, site), ~fct_contr_sum(as.factor(.))) %>%
  mutate(gxs = interaction(line_name, site, drop = TRUE, sep = ":"),
         analysis = "environment") %>%
  mutate_at(vars(line_name, site, gxs), list(cov = f_return)) %>%
  ## Nest the data
  select(which(! names(.) %in% c("trial", "location", "year", "std_error"))) %>%
  group_by(analysis, trait, population) %>%
  nest() %>% ungroup()


# Remove soil variable from all variables
concurrent_all_features <- concurrent_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(concurrent_all_features, .)

# Assemble covariates to model
env_covariates_to_model <- bind_rows(concurrent_fact_reg_feature_selection, concurrent_feature_selection, 
                                     concurrent_all_features)  %>%
  mutate(features = map(covariates, "optVariables") %>% map(~setdiff(x = ., y = "line_name"))) %>%
  filter(source == "daymet", model == "model3", str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  select(trait, feat_sel_type, features) %>%
  mutate(analysis = "environment")




# Location means
# Calculate location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()

## Pull out location covariates to use
loc_ec_tomodel_scaled <-  historical_ec_tomodel_timeframe_scaled %>%
  subset(., grepl(pattern = "daymet", x = names(.)))


loc_data_to_model <- S2_MET_loc_BLUEs %>% 
  filter(location %in% train_test_loc, line_name %in% c(tp_geno, vp_geno)) %>%
  crossing(., population = c("all", "tp", "vp")) %>%
  filter(!(population == "tp" & line_name %in% vp_geno), !(population == "vp" & line_name %in% tp_geno)) %>%
  # Copy columns for modeling covariates
  droplevels() %>%
  # Rename site as location
  rename(site = location) %>%
  mutate_at(vars(line_name, site), ~fct_contr_sum(as.factor(.))) %>%
  mutate(gxs = interaction(line_name, site, drop = TRUE, sep = ":"),
         analysis = "location") %>%
  mutate_at(vars(line_name, site, gxs), list(cov = f_return)) %>%
  ## Nest the data
  select(which(! names(.) %in% c("trial", "location", "year", "std_error"))) %>%
  group_by(analysis, trait, population) %>%
  nest() %>% ungroup()

# Remove soil variable from all variables
historical_all_features <- historical_all_features %>%
  mutate(covariates = map(covariates, ~modify_at(.x, "optVariables", ~str_subset(.x, "soil", negate = TRUE))), 
         feat_sel_type = "all_nosoil") %>%
  bind_rows(historical_all_features, .)


# Combine feature selection df
loc_covariates_to_model <- bind_rows(historical_fact_reg_feature_selection, mutate(historical_feature_selection, source = "daymet"),
                                  historical_all_features) %>%
  select(-contains("time_frame"))  %>%
  mutate(features = map(covariates, "optVariables") %>% map(~setdiff(x = ., y = "line_name"))) %>%
  filter(source == "daymet", model == "model5", str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  select(trait, feat_sel_type, features) %>%
  mutate(analysis = "location")


## First estimate genetic variance compoments; these will be held constant
genetic_variance_components <- bind_rows(env_data_to_model, loc_data_to_model) %>%
  group_by(trait, population, analysis) %>%
  do(varG_hat = {
    
    row <- .
    
    # Trait
    data <- droplevels(row$data[[1]])
    genotypes <- levels(data$line_name)
    analysis <- row$analysis
    
    # Create or subset relationship matrices
    K_use <- K[genotypes, genotypes]
    
    ## Mixed solve
    ms_out <- mixed.solve(y = data$value, K = K_use, Z = model.matrix(~ -1 + line_name, data))
    
    # Return the estimate variance components
    ms_out$Vu
    
  }) %>% ungroup() %>% unnest()



# Combine the 'data_to_model' dataframes
data_to_model_split <- bind_rows(env_data_to_model, loc_data_to_model) %>% 
  full_join(., bind_rows(env_covariates_to_model, loc_covariates_to_model), by = c("analysis", "trait")) %>%
  left_join(., genetic_variance_components, by = c("analysis", "trait", "population")) %>%
  assign_cores(df = ., n_core = n_cores, split = TRUE)


## Iterate over cores
pheno_variance_analysis_out <- coreApply(X = data_to_model_split, FUN = function(core_df) {
  
  # Output list
  out <- vector("list", length = nrow(core_df))
  
  # Create a prefix depending on the core number
  bglr_prefix <- paste0("bglr_core", unique(core_df$core))
  
  # Iterate over rows
  for (i in seq_along(out)) {
    
    row <- core_df[i,]
    
    # Trait
    tr <- row$trait
    data <- droplevels(row$data[[1]]) %>%
      # Center and scale
      mutate(value_scaled = as.numeric(scale(value)))
    genotypes <- levels(data$line_name)
    sites <- levels(data$site)
    covariates <- row$features[[1]]
    analysis <- row$analysis
    
    # Get the genetic variance component
    varG_hat <- row$varG_hat
    
    # Separate covariates into main effect or interaction
    covariate_list <- tibble(term = unique(covariates)) %>% 
      mutate(class = ifelse(str_detect(term, ":"), "interaction", "main"),
             covariate = str_remove(term, "line_name:")) %>%
      split(.$class) %>% 
      map("covariate")
    
    ## Create covariate matrix depending on analysis type
    if (analysis == "environment") {
      
      # Create a covariate matrix
      covariate_mat <- ec_tomodel_scaled$daymet %>%
        as.data.frame() %>%
        select(-source) %>%
        column_to_rownames("environment") %>%
        as.matrix()
      
    } else {
    
      # Get the time_frame corresponding to this trait
      time_frame_use <- unique(subset(historical_feature_selection, trait == tr, time_frame, drop = TRUE))
      
      covariate_mat <- loc_ec_tomodel_scaled[[which(grepl(pattern = time_frame_use, x = names(loc_ec_tomodel_scaled)))]] %>%
        as.data.frame() %>%
        select(-source, -time_frame) %>%
        column_to_rownames("location") %>%
        as.matrix()
      
    }
    
    
    # Create or subset relationship matrices
    K_use <- K[genotypes, genotypes]
    Zg <- model.matrix(~ -1 + line_name, data)
    Emain <- Env_mat(x = covariate_mat[sites, covariate_list$main, drop = FALSE], method = "Jarq")
    # Emain; but with zero off-diagonal
    Emain1 <- Emain
    Emain1[] <- 0
    diag(Emain1) <- diag(Emain)
    
    Eint <- Env_mat(x = covariate_mat[sites, covariate_list$int, drop = FALSE], method = "Jarq")
    # Replace with diagonal matrix if NA
    if (all(is.na(Eint))) { 
      Eint[] <- 0
      diag(Eint) <- 1 
      no_gxe_cov <- TRUE
      
    } else {
      no_gxe_cov <- FALSE
      
    }
    
    Ze <- model.matrix(~ -1 + site, data)
    Zge <- model.matrix(~ -1 + gxs, data)
    KE <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% Eint, Ze)
    dimnames(KE) <- replicate(2, data$gxs, simplify = FALSE)
    
    # Another matrix with unrelated environments
    KE1 <- tcrossprod(Zg %*% K_use, Zg) * tcrossprod(Ze %*% diag(ncol(Eint)), Ze)
    dimnames(KE1) <- dimnames(KE)
    
    # Reset the flag 
    parse_model_out <- TRUE
    
    # Matrix specifying constraints
    nem <- as.matrix(3) # Not to be estimated but fixed
    varGt <- as.matrix(varG_hat)
    cm <- as.matrix(1)
    
    ## Create formula objects
    fixed_form <- value ~ 1
    
    # If no_gxe_cov is TRUE, drop gxs_cov
    if (no_gxe_cov) {
      rand_form <- ~ vs(line_name_cov, Gu = K_use, Gt = varGt, Gtc = nem) + 
        vs(site, Gu = Emain1) + vs(site_cov, Gu = Emain) +
        vs(gxs, Gu = KE1)
      
    } else {
      rand_form <- ~ vs(line_name_cov, Gu = K_use, Gt = varGt, Gtc = nem) + 
        vs(site, Gu = Emain1) + vs(site_cov, Gu = Emain) +
        vs(gxs, Gu = KE1) + vs(gxs_cov, Gu = KE)
      
    }
    
    
    ## Sommer
    model_stdout <- capture.output({
      model_try <- try( fit_mmer <- mmer(fixed = fixed_form, random = rand_form, data = data, verbose = TRUE), silent = TRUE) 
    })
    
    # Parse the output
    stdout_parsed <- suppressWarnings(read_table2(model_stdout))
    # Max LL iteration
    maxLL_iter <- which.max(as.numeric(stdout_parsed$LogLik))
    
    # If the output contains singular, refit the model using the iteration that maximized the LL
    if (any(str_detect(model_stdout, "singular")) | maxLL_iter != max(as.numeric(stdout_parsed$iteration))) {

      # If maxLL_iter is empty, just skip; else refit
      if (is_empty(maxLL_iter)) {
        varcomp_df <- NULL
        parse_model_out <- FALSE # Use this as a flag
        
      } else {

        ## Sommer
        model_stdout <- capture.output({
          model_try <- try( fit_mmer <- mmer(fixed = fixed_form, random = rand_form, data = data, verbose = TRUE, 
                                             iters = maxLL_iter), 
                            silent = TRUE) 
        })
        
        
      }
      
    }
    
    if (parse_model_out) {
    
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
    
    
    # ## BGLR
    # ## Try BGLR
    # fit_bglr <- BGLR(y = data$value_scaled, verbose = FALSE, burnIn = 2000, nIter = 12000, saveAt = bglr_prefix,
    #                  ETA = list( genoK = list(X = Zg, K = K_use, model = "BRR"), genoI = list(X = Zg, K = diag(ncol(Zg)), model = "BRR"),
    #                              siteK = list(X = Ze, K = Emain, model = "BRR"), siteI = list(X = Ze, K = diag(ncol(Ze)), model = "BRR"),
    #                              gxsK = list(X = Zge, K = KE, model = "BRR"), gxsI = list(X = Zge, K = diag(ncol(Zge)), model = "BRR")
    #                  ) )
    # 
    # ## Extract variance components and tidy
    # varcomp_df <- c(map_dbl(fit_bglr$ETA, "varB"), Residuals = fit_bglr$varE) %>% 
    #   tibble(term = names(.), variance = .) %>%
    #   mutate(source = str_remove_all(term, "[KI]$")) %>%
    #   group_by(source) %>%
    #   mutate(total_source_variance = sum(variance)) %>%
    #   ungroup() %>%
    #   mutate(variance_prop = variance / total_source_variance,
    #          note = list(NULL))
    # 
    
      

    
    
    # Return
    out[[i]] <- varcomp_df
    
    # Print a message
    cat("\nVariance components estimated for", unlist(select(row, analysis, trait, population, feat_sel_type)) %>% 
          paste0(names(.), ": ", .) %>% paste0(., collapse = ", "))
    
  } # Close loop
  
  # Add out to core_df and return
  core_df %>%
    mutate(results = out) %>%
    select(-data, -core)
  
})





# Save the df
save("pheno_variance_analysis_out", file = file.path(result_dir, "full_model_variance_analysis.RData"))
