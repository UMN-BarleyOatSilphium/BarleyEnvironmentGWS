## S2MET Predictions
## 
## Script for analyzing cross-validation results
## 
## Author: Jeff Neyhart
## Last modified: May 14, 2019
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# # Load some packages
library(lubridate)
library(effects) # For ls means / marginal means
library(ggforce)
library(ggridges)
library(gridExtra)
library(broom)
library(car)
library(lmerTest)


## Load the validation results
file_list <- list.files(result_dir, "results.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))


## Significant level
alpha <- 0.05


## A vector to rename models
model_replace <- c("M1" = "M1 (Main effect)", "M2" = "M2 (GxE)")


## combine all data and assign the scheme number
pov00_predictions <- rename(pov00_predictions, environment = val_environment)

results_df <- map(object_list, get) %>%
  map_df(ungroup) %>%
  select(-core) %>%
  mutate(scheme_number = str_extract(scheme, "[0-9]{1,2}")) %>%
  mutate_at(vars(model, scheme), as.factor)


## 
## Fit a model
## 

model_results <- results_df %>% 
  group_by(trait, scheme_number) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(model_results))) {
  
  df <- model_results$data[[i]]
  
  if (model_results$scheme_number[i] == "00") {
    next
    
  } else {
    
    fit <- lmer(accuracy ~ model + scheme + model:scheme + (1|environment) + (1|environment:model) + 
                  (1|environment:scheme)+ (1|environment:model:scheme), 
                data = df)
  }
  
  out <- data_frame(model = list(fit),
                    effects = list(as.data.frame(Effect(focal.predictors = c("model", "scheme"), fit))),
                    ranova = list(tidy(ranova(fit))),
                    anova = list(tidy(anova(fit))))
  
  model_results$out[[i]] <- out
  
}
                    

## Unnest
model_results1 <- model_results %>%
  filter(!map_lgl(out, is.null)) %>%
  unnest(out)


## Plot
model_results1 %>% 
  unnest(effects) %>%
  # filter(scheme_number == "2") %>%
  ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
  geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
  geom_point(position = position_dodge(0.7)) +
  facet_grid(trait ~ scheme_number, scales = "free_x") +
  theme_presentation2()




## Model within scheme
## 
## POV00
## 
pov00_analysis <- pov00_predictions %>% 
  mutate(zscore = ztrans(accuracy)) %>%
  group_by(trait, scheme) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(pov00_analysis))) {
  
  df <- pov00_analysis$data[[i]]
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment), data = df, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  ## Rescale accuracy and return
  pov00_analysis$out[[i]] <- data_frame(
    model = list(fit),
    effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
    ranova = list(tidy(ranova(fit))),
    anova = list(tidy(anova(fit)))
  )
  
}

pov00_analysis1 <- unnest(pov00_analysis, out)

## Plot
g_pov00 <- pov00_analysis1 %>% 
  unnest(effects) %>%
  mutate(model = str_replace_all(model, model_replace),
         model = factor(model, levels = model_replace)) %>%
  ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
  geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
  geom_point(position = position_dodge(0.7), size = 3) +
  # facet_grid(trait ~ scheme_number, scales = "free_x") +
  facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
  scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
  scale_x_discrete(name = "Validation scheme", labels = toupper) +
  scale_color_discrete(name = "Model") +
  theme_presentation2() +
  theme(legend.position = "bottom")

ggsave(filename = "pov00_accuracy_analysis.jpg", plot = g_pov00, path = fig_dir, width = 7, height = 4, dpi = 1000)





## Model within scheme
## 
## POV1
## 

# Histogram
qplot(x = accuracy, data = pov1_predictions, geom = "density", facets = ~trait)

## Needs transformation

pov1_analysis <- pov1_predictions %>% 
  mutate(zscore = ztrans(accuracy),
         model = as.factor(model)) %>%
  group_by(trait, scheme) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(pov1_analysis))) {
  
  df <- pov1_analysis$data[[i]]
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment), data = df, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
  
  ## Rescale accuracy and return
  pov1_analysis$out[[i]] <- data_frame(
    model = list(fit),
    effects = list(Effect(fit, focal.predictors = "model") %>% as.data.frame() %>% mutate_at(vars(-model), zexp)),
    ranova = list(tidy(ranova(fit))),
    anova = list(tidy(anova(fit)))
  )
  
}

pov1_analysis1 <- unnest(pov1_analysis, out)

## Plot
g_pov1 <- pov1_analysis1 %>% 
  unnest(effects) %>%
  mutate(model = str_replace_all(model, model_replace),
         model = factor(model, levels = model_replace)) %>%
  ggplot(aes(x = scheme, y = fit, ymax = upper, ymin = lower, color = model, group = model)) +
  geom_errorbar(color = "black", position = position_dodge(0.7), width = 0.5) +
  geom_point(position = position_dodge(0.7), size = 3) +
  # facet_grid(trait ~ scheme_number, scales = "free_x") +
  facet_grid(~ trait, scales = "free_x", labeller = labeller(trait = str_add_space)) +
  scale_y_continuous(name = "Predicton accuracy", breaks = pretty) +
  scale_x_discrete(name = "Validation scheme", labels = toupper) +
  scale_color_discrete(name = "Model") +
  theme_presentation2() +
  theme(legend.position = "bottom")

ggsave(filename = "pov1_accuracy_analysis.jpg", plot = g_pov1, path = fig_dir, width = 7, height = 4, dpi = 1000)


