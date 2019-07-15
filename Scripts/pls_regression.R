## S2MET Prediction Models
## 
## Selection of ECs via PLS
## 
## This script will look at the weighing of ECs by PLS regression
## 
## Author: Jeff Neyhart
## Last modified: June 4, 2019
## 

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(pls)


## Load the distance methods
load(file.path(alt_proj_dir, "Results/distance_method_results.RData"))
# Load the environmental variables
load(file.path(alt_proj_dir, "Data/environmental_data_compiled.RData"))
# Load the AMMI results
load(file.path(alt_proj_dir, "Results/genotype_environment_phenotypic_analysis.RData"))


## Calculate GxE residuals for each trait
training_set_gxe <- training_sets_twoway %>%
  filter(set %in% c("complete")) %>%
  filter(trait %in% traits) %>%
  group_by(trait) %>%
  do({
    row <- .
    
  })





### Example using vargas dataset

library(agridat)


data(vargas.wheat1.covs)
data(vargas.wheat1.traits)

# Tidy
df <- vargas.wheat1.traits %>% 
  group_by(year, gen) %>%
  summarize(yield = mean(yield)) %>%
  ungroup() %>%
  mutate_at(vars(gen, year), as.factor)

# Fit a g + y main effect model
fit <- lm(yield ~ gen + year, df, contrasts = list(gen = "contr.sum", year = "contr.sum"))

## Add residuals to df
df2 <- add_residuals(data = df, model = fit)

# Create GxE matrix
Y_gxe <- scale(xtabs(resid ~ year + gen, df2))

## Prepare covariates
# Yield as a function of environment covariates
cov_mat <- vargas.wheat1.covs %>%
  column_to_rownames("year") %>%
  as.matrix() %>%
  scale()

## Run PLS regressions
pls_fit <- plsr(Y_gxe ~ cov_mat)
# Get the coefficients
ec_beta <- loadings(pls_fit)[,1]
# Number of coefficients
m <- ncol(cov_mat)

# Create a weight matrix
D <- diag(1 / (m * abs(ec_beta)))
E <- cov_mat %*% D %*% t(cov_mat)

## The alternative is simply to divide by m
E1 <- tcrossprod(cov_mat) / m

## Lastly, the correlation among environments
E2 <- cor(xtabs(yield ~ gen + year, df2))

## Relationship matrix of individuals
K <- diag(ncol(Y_gxe)); dimnames(K) <- list(colnames(Y_gxe), colnames(Y_gxe))
Kge <- kronecker(E, K, make.dimnames = T)
Kge1 <- kronecker(E1, K, make.dimnames = T)
Kge2 <- kronecker(E2, K, make.dimnames = T)


## Factorial regression using ECs in order of descending PLS coefficient
ec_order <- names(ec_beta)[order(abs(ec_beta), decreasing = TRUE)]
# Only use the first n ECs, where n is the number of environments
ec_order_use <- ec_order[seq_along(levels(df_tomodel$year))]

## Combine ECs with df
df_tomodel <- left_join(df, mutate(vargas.wheat1.covs, year = as.factor(year)))

# Define the maximum model for forward stepwise regression
max_model <- add_predictors(f = formula(fit), as.formula(paste0("~", paste0("gen:", ec_order_use, collapse = " + "))))
# update the model with new data
fit1 <- update(object = fit, data = df_tomodel)

# Fit a g + y + beta h with forward stepwise regression
fit_step <- step(object = fit1, direction = "both", scope = max_model, data = df_tomodel)



## Fit a model
library(lmerTest)
library(lme4qtl)

# Prepare for modeling
df3 <- df %>%
  mutate(year2 = year, gen2 = gen)

mm_fit1 <- relmatLmer(yield ~ (1|gen) + (1|year), df3)
VarProp(mm_fit1)

# With year matrix
mm_fit2 <- relmatLmer(yield ~ (1|gen) + (1|year2), df3, relmat = list(year2 = E))
VarProp(mm_fit2)


## What matrix leads to the highest predictions?
# Leave-one-year-out CV
uniq_year <- unique(df2$year)
year_acc <- numeric(length(uniq_year))

for (i in seq_along(year_acc)) {
  # Subset data.frame
  mf <- model.frame(yield ~ gen + year, df2, subset = year != uniq_year[i])
  y <- model.response(mf)
  X <- model.matrix(~ 1, mf)
  Zg <- model.matrix(~ -1 + gen, mf); colnames(Zg) <-  levels(mf$gen)
  Ze <- model.matrix(~ -1 + year, mf); colnames(Ze) <-  levels(mf$year)
  Zge <- model.matrix(~ -1 + gen:year, mf); colnames(Zge) <- colnames(Kge)

  
  # ## Fit models
  # fit_somm1 <- sommer::mmer(Y = y, X = X, Z = list(gen = list(Z = Zg, K = K), year = list(Z = Ze, K = E)), silent = T)
  # fit_somm2 <- sommer::mmer(Y = y, X = X, Z = list(gen = list(Z = Zg, K = K), year = list(Z = Ze, K = E1)), silent = T)
  # fit_somm3 <- sommer::mmer(Y = y, X = X, Z = list(gen = list(Z = Zg, K = K), year = list(Z = Ze, K = E2)), silent = T)
  
  # Fit models
  fit_somm1 <- sommer::mmer(Y = y, X = X, silent = T,
                            Z = list(gen = list(Z = Zg, K = K), year = list(Z = Ze, K = E), int = list(Z = Zge, K = Kge)))
  fit_somm2 <- sommer::mmer(Y = y, X = X, silent = T,
                            Z = list(gen = list(Z = Zg, K = K), year = list(Z = Ze, K = E1), int = list(Z = Zge, K = Kge1)))
  fit_somm3 <- sommer::mmer(Y = y, X = X, silent = T,
                            Z = list(gen = list(Z = Zg, K = K), year = list(Z = Ze, K = E2), int = list(Z = Zge, K = Kge2)))
  
  ## List and calculate phenotypes
  predictions_list <- list(fit_somm1, fit_somm2, fit_somm3) %>%
    map(~map(.$u.hat, ~as.data.frame(.) %>% rownames_to_column("term")) %>% 
          map2(.x = ., .y = names(.), ~{names(.) <- c(.y, "value"); return(.)}) %>%
          do.call("crossing", .) %>% 
          separate(col = int, into = c("year1", "gen1"), sep = ":") %>%
          filter(gen == gen1, year == year1) %>%
          mutate(yield_hat = value + value1 + value2) )
  
  ## Calculate accuracy
  year_acc[i] <- predictions_list %>%
    setNames(object = ., nm = c("scaled", "unscaled", "cor")) %>%
    map_df(~left_join(subset(df2, year == uniq_year[i]), .) %>% summarize(acc = cor(yield, yield_hat)))
  
}





