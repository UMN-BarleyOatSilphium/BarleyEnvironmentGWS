## S2MET Predictions Models
## 
## Scripts to generate figures
## 
## 

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(broom)
library(patchwork)


### Plot Environmental Covariables ####

## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))

# Rename
ec_model_building <- unified_ec_models %>%
  unnest(final_model) %>%
  filter(model == "model3_ammi")


## For each model, list the random and fixed effect covariates in order
## of descending var comp or regression coefficient
ec_model_table <- ec_model_building %>%
  mutate(fixefs = map(object, ~fixef(.x)[-1] %>%  tibble(term = names(.), coef = .) %>% 
                        arrange(desc(abs(coef))) ),
         ranefs = map(object, ~as.data.frame(VarCorr(.x)) %>% filter(var1 != "(Intercept)", is.na(var2)) %>%
                        arrange(desc(vcov)) %>% select(term = var1, variance = vcov) ),
         # ranefs = map(object, ~as.data.frame(VarCorr(.x)) %>% filter(! grp %in% c("line_name", "Residual")) %>%
         #                arrange(desc(vcov)) %>% select(term = var1, variance = vcov) ),
         effects = map2(ranefs, fixefs, bind_rows) ) %>%
  unnest(effects) %>%
  select(trait, term, variance, coef)

## Output table
write_csv(x = ec_model_table, path = file.path(result_dir, "ec_variance_coef.csv"))


## Plot distribution of random effects and coefficients of fixed effects
# First get the predicted random effects
ec_pred_ranef <- ec_model_building %>% 
  mutate(pred_ranef = map(object, ~ranef(.x)[[1]] %>% rownames_to_column("level") %>%
                            gather(term, effect, -level) %>% filter(term != "(Intercept)") )) %>%
  left_join(., group_by(ec_model_table, trait) %>% nest(.key = "term_var_coef"))


## Plot
ec_model_plots <- ec_pred_ranef %>%
  group_by(trait) %>%
  do(plot = {
    row <- .
    # Get the fitted model object
    fit <- row$object[[1]]
    
    ## Plot 1 shows the main regression coefs and CI
    fit_coefs <- summary(fit)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      left_join(., confint.merMod(object = fit, parm = "beta_", method = "Wald") %>% 
                  as.data.frame() %>% rownames_to_column("term") ) %>%
      filter(term != "(Intercept)") %>%
      rename_all(make.names) %>%
      rename_all(tolower) %>%
      # Reorder
      arrange(desc(abs(estimate))) %>%
      mutate(term = factor(term, levels = term))
    
    # Create plot 1
    plot1 <- fit_coefs %>%
      ggplot(aes(x = term, y = estimate, ymin = x2.5.., ymax = x97.5..)) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_pointrange(shape = 1, color = "blue") +
      coord_flip() +
      theme_minimal(base_size = 8) +
      labs(subtitle = "Main effect coefficient")
    
    
    ## Plot 2 shows histograms of random effects
    plot2 <- unnest(row, pred_ranef) %>%
      ggplot(aes(x = effect)) +
      geom_histogram() +
      facet_wrap(~ term, scales = "free_x") +
      theme_minimal(base_size = 8) +
      labs(subtitle = "Random interaction effects")
    
    
    ## Combine plots and return
    plot1 + plot2 + 
      plot_layout(widths = c(0.75, 1))
    
    #
    
  }) %>% ungroup()


## Loop and save figures
for (i in seq(nrow(ec_model_plots))) {
  tr <- ec_model_plots$trait[i]
  filename <- paste0("ec_model_effects_", tr, ".jpg")
  
  ggsave(filename = filename, plot = ec_model_plots$plot[[i]], path = fig_dir,
         height = 4, width = 8, dpi = 300)
  
}







#####################
## Covariate identification
#####################

## Plot the results
g_hist_ec_ammi <- historical_ec_ammi_dist %>% 
  unnest(test_results) %>% 
  ggplot(aes(x = number, y = cor_with_ammi, color = trait, lty = time_frame)) +
  # geom_point() +
  geom_line() +
  facet_grid(~ trait, scales = "free_x") +
  theme_acs()
ggsave(filename = "ec_locations_correlation_with_AMMI.jpg", plot = g_hist_ec_ammi, 
       path = fig_dir, width = 8, height = 4, dpi = 1000)











