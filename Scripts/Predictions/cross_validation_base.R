## S2MET cross-validation
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: January 17, 2019
## 
## This is the base script from which other scripts will drawn upon
## Cross-validation is initiated by drawing 50 replicates of 5-folds
## 
## This script should be run once, and then the results will be used as starting material
## for other scripts.
## 


# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)



## Some parameters
# Number of folds
nFolds <- 5

# Number of CV iterations
nCV <- 50







## Sample environments to keep
set.seed(1122)
sample_envs <- trial_info %>% 
  filter(environment %in% reduce(tp_vp_env_trait, intersect)) %>%
  distinct(environment, year) %>% 
  sample_n(size = 10) %>%
  arrange(year)



## Data to use for CV
cv_data <- S2_MET_BLUEs %>%
  ## Sample
  filter(environment %in% sample_envs$environment) %>%
  filter(line_name %in% tp_geno,
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.)))


## Generate CV randomizations
# First generate training/test lines
set.seed(1171039)
train_test <- replicate(n = nCV, crossv_kfold(data = distinct(cv_data, line_name), k = nFolds), simplify = FALSE) %>% 
  map2_df(.x = ., .y = seq_along(.), ~mutate(.x, .id = paste0(.y, "_", .id))) %>%
  mutate_at(vars(train, test), funs(map(., as.data.frame) %>% map("line_name"))) %>%
  crossing(trait = traits, .)

# CV1 - predict totally untested genotypes using all environments
cv1_rand <- train_test %>%
  mutate_at(vars(train, test), funs(map2(.x = ., .y = trait, ~filter(cv_data, line_name %in% .x, trait == .y))))

# CV2 - predict partially tested genotypes
# This randomly assigns observations in environments to a fivefold partition
cv2_rand <- cv_data %>%
  group_by(trait) %>%
  # do(crossv_mc(data = ., n = nCV, test = pTest)) %>%
  do({
    replicate(n = nCV, crossv_kfold(data = ., k = nFolds), simplify = FALSE) %>%
      map2_df(.x = ., .y = seq_along(.), ~mutate(.x, .id = paste0(.y, "_", .id))) }) %>%
  ungroup() %>%
  mutate_at(vars(train, test), funs(map(., as.data.frame)))
  

## Predict unobserved environments
## Leave-one-out predictions

# CV0 - predict previously tested genotypes in unobserved environments
cv0_rand_loeo <- cv_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    distinct(df, environment) %>% 
      mutate(train = map(environment, ~filter(df, environment != .)), 
             test = map(environment, ~filter(df, environment == .))) 
    
  }) %>% ungroup() %>%
  mutate(.id = "01")

# CV00 - predict totally untested genotypes in totally untested environments
cv00_rand_loeo <- full_join(train_test, select(cv0_rand_loeo, -.id), by = "trait") %>%
  mutate(train = map2(train.x, train.y, ~filter(.y, line_name %in% .x)),
         test = map2(test.x, test.y, ~filter(.y, line_name %in% .x))) %>%
  select(trait, environment, .id, train, test) %>%
  arrange(trait, environment)


## Predict 2017 using 2015 and 2016
# CV0 - predict previously tested genotypes in unobserved environments
cv0_rand_future <- cv_data %>%
  group_by(trait) %>%
  do(data_frame(train = list(filter(., year != 2017)), test = list(filter(., year == 2017)))) %>% 
  ungroup() %>%
  mutate(.id = "01")

# CV00 - predict totally untested genotypes in totally untested environments
cv00_rand_future <- full_join(train_test, select(cv0_rand_future, -.id), by = "trait") %>%
  mutate(train = map2(train.x, train.y, ~filter(.y, line_name %in% .x)),
         test = map2(test.x, test.y, ~filter(.y, line_name %in% .x))) %>%
  select(trait, .id, train, test) %>%
  arrange(trait)



## Combine CV1 and CV2
CV12 <- bind_rows(
  mutate(cv1_rand, cv = "cv1"),
  mutate(cv2_rand, cv = "cv2")
)

CV_zero_loeo <- bind_rows(
  mutate(cv0_rand_loeo, cv = "cv0"),
  mutate(cv00_rand_loeo, cv = "cv00")
)

CV_zero_future <- bind_rows(
  mutate(cv0_rand_future, cv = "cv0"),
  mutate(cv00_rand_future, cv = "cv00")
)





# ### Parent-offspring cross-validation
# 
# ## Data to use for CV
# pocv_data <- S2_MET_BLUEs %>% 
#   filter(environment %in% unique(cv_data$environment)) %>%
#   filter(line_name %in% c(tp_geno, vp_geno),
#          trait %in% traits,
#          environment %in% tp_vp_env) %>%
#   mutate(id = seq(nrow(.)))
# 
# 
# ## Generate CV randomizations
# # First generate training/test lines
# train_test1 <- train_test %>%
#   mutate(test = list(vp_geno))
# 
# 
# # POCV1 - predict the VP using all environments
# pocv1_rand <- train_test1 %>%
#   mutate_at(vars(train, test), funs(map2(.x = ., .y = trait, ~filter(pocv_data, line_name %in% .x, trait == .y))))
# 
# # POCV2 - predict partially tested VP
# pocv2_rand <- pocv_data %>%
#   group_by(trait) %>%
#   # do(crossv_mc(data = ., n = nCV, test = pTest)) %>%
#   do({
#     replicate(n = nCV, crossv_kfold(data = ., k = nFolds), simplify = FALSE) %>%
#       map2_df(.x = ., .y = seq_along(.), ~mutate(.x, .id = paste0(.y, "_", .id))) }) %>%
#   ungroup()
# 
# # POCV0 - predict totally tested VP in unobserved environments
# pocv0_rand_loeo <- pocv_data %>%
#   group_by(trait) %>%
#   do({
#     df <- .
#     
#     distinct(., environment) %>% 
#       mutate(train = map(environment, ~filter(df, environment != .)), 
#              test = map(environment, ~filter(df, environment == .))) }) %>% 
#   ungroup() %>%
#   mutate(.id = "01")
# 
# # POCV00 - predict totally untested VP in totally untested environments
# pocv00_rand_loeo <- full_join(train_test1, select(pocv0_rand_loeo, -.id), by = "trait") %>%
#   mutate(train = map2(train.x, train.y, ~filter(.y, line_name %in% .x)), # Take the available training data and subset the TP
#          test = map2(test.x, test.y, ~filter(.y, line_name %in% .x))) %>% # Take the available testing data and subset the VP
#   select(trait, environment, .id, train, test) %>%
#   arrange(trait, environment)
# 
# 
# ## Predict 2017 using 2015 and 2016
# # POCV0 - predict previously tested genotypes in unobserved environments
# pocv0_rand_future <- pocv_data %>%
#   group_by(trait) %>%
#   do(data_frame(train = list(filter(., year != 2017)), test = list(filter(., year == 2017)))) %>% 
#   ungroup() %>%
#   mutate(.id = "01")
# 
# # POCV00 - predict totally untested genotypes in totally untested environments
# pocv00_rand_future <- full_join(train_test, select(pocv0_rand_future, -.id), by = "trait") %>%
#   mutate(train = map2(train.x, train.y, ~filter(.y, line_name %in% .x)),
#          test = map2(test.x, test.y, ~filter(.y, line_name %in% .x))) %>%
#   select(trait, .id, train, test) %>%
#   arrange(trait)
# 
# 
# 
# 
# ## Combine POCV1 and POCV2
# POCV12 <- bind_rows(
#   mutate(pocv1_rand, cv = "pocv1"),
#   mutate(pocv2_rand, cv = "pocv2")
# )
# 
# POCV_zero_loeo <- bind_rows(
#   mutate(pocv0_rand_loeo, cv = "pocv0"),
#   mutate(pocv00_rand_loeo, cv = "pocv00")
# )
# 
# POCV_zero_future <- bind_rows(
#   mutate(pocv0_rand_future, cv = "pocv0"),
#   mutate(pocv00_rand_future, cv = "pocv00")
# )


## Save
if (all(tp_vp_env %in% unique(cv_data$environment))) {
  save("CV12", "CV_zero_future", "CV_zero_loeo", file = file.path(proj_dir, "Scripts/Predictions/CrossValidation/cv_starting_material.RData"))
  
} else {
  save("CV12", "CV_zero_future", "CV_zero_loeo", file = file.path(proj_dir, "Scripts/Predictions/CrossValidation/cv_starting_material_sample.RData"))
  
}


