## S2MET GWAS
## 
## Map response to environmental conditions
## 
## Author: Jeff Neyhart
## Last modified: 14 October 2019
## 

# Run on a local machine
repo_dir <- getwd()
# Other packages
library(modelr)
library(broom)
library(sommer)

source(file.path(repo_dir, "source.R"))



## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))

K_use <- K[tp_geno, tp_geno]


# ## Convert markers for mapping
# geno <- snp_info %>% 
#   rename(marker = 1) %>% 
#   select(marker, chrom, pos) %>%
#   as.data.frame() %>%
#   cbind(., t(s2_imputed_mat_use)) %>%
#   select(marker:pos, tp_geno)
   
##

## Load BOPA/GBS SNP data
bopa_geno <- read_tsv(file = "C:/GoogleDrive/BarleyLab/Breeding/GenotypicData/Combined_GBS_BOPA_Genotype_Data/S2TP_final_discrete_multi_genos_hmp.txt")

bopa_geno2 <- bopa_geno %>%
  rename(marker = 1) %>%
  select(marker, chrom, pos, tp) %>%
  as.data.frame()

## Impute
bopa_impute <- bopa_geno2 %>%
  select(tp) %>%
  t() %>%
  `colnames<-`(., bopa_geno2$marker) %>%
  A.mat(X = ., min.MAF = 15 / nrow(.), max.missing = 0.8, impute.method = "EM", return.imputed = TRUE)

# Get the K matrix
K_use <- bopa_impute$A
# Get imputed markers
bopa_impute1 <- bopa_impute$imputed

# Convert to data.frame
geno <- t(bopa_impute1) %>% 
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  right_join(select(bopa_geno2, marker, chrom, pos), .)


## Organize phenotypes for mapping
pheno <- covariate_reg_coefs %>% 
  unnest(apriori_model) %>%
  unite(map_trait, c("trait", "covariate")) %>%
  spread(map_trait, estimate) %>%
  as.data.frame() %>%
  filter(line_name %in% names(geno))

# # Map
# map_out <- GWAS(pheno = pheno, geno = geno, K = K_use, min.MAF = 0)

# Organize phenotypes for GWAS 
pheno1 <- unnest(covariate_reg_coefs, apriori_model) %>%
  filter(line_name %in% names(geno))

# Marker data to use
# M <- s2_imputed_mat_use[tp_geno,, drop = T]
M <- t(select(geno, -marker, -chrom, -pos)) %>%
  `colnames<-`(., geno$marker)

# Group by trait and covariate and map
sommer_gwas <- pheno1 %>%
  group_by(trait, covariate) %>%
  do({
    df <- .
    
    ##
    out <- sommer::GWAS(fixed = estimate ~ 1, random = ~ vs(line_name, Gu = K_use), 
                        data = df, M = M, gTerm = "u:line_name", verbose = FALSE, 
                        min.MAF = 0)
    
    ## 

    marker_scores <- t(out$scores) %>% 
      as.data.frame() %>% 
      rownames_to_column("marker") %>%
      rename_all(~str_remove_all(., "estimate "))
    
    ## Get FDR
    fdr_level <- sommer:::fdr(p = 10^-marker_scores$score)$fdr.10
    
    # Get variance components
    var_comp <- setNames(out$sigma[,,], c("line_name", "Residuals"))
    
    tibble(marker_scores = list(marker_scores), var_comp = list(var_comp), fdr = fdr_level)
    
  })



# # Group by trait and covariate and map
# # Single-step GWAS
# sommer_ssgwas <- pheno1 %>%
#   group_by(trait, covariate) %>%
#   do({
#     df <- .
#     
#     # Model frame
#     mf <- model.frame(estimate ~ line_name, df)
#     y <- model.response(mf)
#     X <- model.matrix(~ 1, mf)
#     
#     ##
#     fit <- EMMREML::emmreml(y = y, X = X, Z = M, K = diag(ncol(M)), varbetahat = TRUE, PEVuhat = TRUE,
#                             varuhat = TRUE, test = FALSE)
#         
#     # Get marker scores
#     
#     
#     
#     ## 
#     
#     marker_scores <- t(out$scores) %>% 
#       as.data.frame() %>% 
#       rownames_to_column("marker") %>%
#       rename_all(~str_remove_all(., "estimate "))
#     
#     ## Get FDR
#     fdr_level <- sommer:::fdr(p = 10^-marker_scores$score)$fdr.10
#     
#     # Get variance components
#     var_comp <- setNames(out$sigma[,,], c("line_name", "Residuals"))
#     
#     tibble(marker_scores = list(marker_scores), var_comp = list(var_comp), fdr = fdr_level)
#     
#   })
# 
# 


## Plot - only those with significant QTL
sommer_gwas %>%
  unnest(marker_scores) %>%
  inner_join(select(geno, marker, chrom, pos), .) %>%
  ggplot(aes(x = pos, y = score)) +
  geom_point() +
  facet_grid(trait + covariate ~ chrom, scales = "free", space = "free_x") +
  theme_manhattan()


## Select a few examples
sommer_gwas %>%
  unnest(marker_scores) %>%
  inner_join(select(geno, marker, chrom, pos), .) %>%
  filter(trait == "GrainYield" & covariate == "flowering_mint" |
         trait == "GrainProtein" & covariate == "grain_fill_maxt") %>%
  ggplot(aes(x = pos, y = score)) +
  geom_point() +
  facet_grid(trait + covariate ~ chrom, scales = "free_x", space = "free_x", switch = "both") +
  theme_presentation2(base_size = 10) +
  theme(strip.background.x = element_rect())




## Plot























