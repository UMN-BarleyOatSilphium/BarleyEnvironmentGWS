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



