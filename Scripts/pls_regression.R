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
library(broom)


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
  mutate_at(vars(gen, year), as.factor)

# Fit a g + y main effect model
fit <- lm(yield ~ gen + year + gen:year, df, contrasts = list(gen = "contr.sum", year = "contr.sum"))

# Create GxE matrix
Y_gxe <- model.tables(aov(fit), type = "effects", cterms = c("gen:year"))
Y_gxe <- Y_gxe$tables$`gen:year` %>% 
  as.matrix() %>%
  t()

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


## Factorial regression using ECs in order of descending PLS coefficient
ec_order <- names(ec_beta)[order(abs(ec_beta), decreasing = TRUE)]

## Combine ECs with df
df_tomodel <- left_join(df, mutate(vargas.wheat1.covs, year = as.factor(year)))

# Only use the first n ECs, where n is the number of environments
ec_order_use <- ec_order[seq_along(levels(df_tomodel$year))]

# Define the maximum model for forward stepwise regression
# First drop gxe term
formula_noGxE <- formula(drop.terms(termobj = terms(formula(fit)), dropx = 3, keep.response = T))
  
max_model <- formula_noGxE %>%
  # Add ECs
  add_predictors(as.formula(paste0("~", paste0("gen:", ec_order_use, collapse = " + ")))) %>%
  # Add GxE
  add_predictors(~ gen:year)


# update the model with new data
fit1 <- update(object = fit, data = df_tomodel)
fit_noGxE <- lm(formula_noGxE, df_tomodel)

## Regular model anova
anova_base <- anova(fit1)

## Fit a g + y + beta h model with the first EC with the largest PLS coefficient
formula_ec1 <- yield ~ gen + year + gen:SHF + gen:year
fit_ec1 <- lm(formula = formula_ec1, data = df_tomodel)

## Proper F test is EC versus GxE
anova_ec1 <- tidy(anova(fit_ec1)) %>%
  mutate(gxe_meansq = list(subset(anova_ec1, term == "gen:year", c(df, meansq)))) %>%
  unnest() %>%
  mutate(fstat = ifelse(str_detect(term, ":"), meansq / meansq1, NA),
         pvalue = pf(q = fstat, df1 = df, df2 = df1, lower.tail = FALSE))

fit_ecAll <- lm(formula = max_model, data = df_tomodel)




# Fit a g + y + beta h with forward stepwise regression
fit_step <- step(object = fit_noGxE, direction = "forward", scope = max_model, data = df_tomodel)






