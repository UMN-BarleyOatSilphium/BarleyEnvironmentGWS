## S2MET GWAS
## 
## Map response to environmental conditions
## 
## Author: Jeff Neyhart
## Last modified: 14 October 2019
## 

# Directories and other packages
repo_dir <- getwd()

# Other packages
library(modelr)
library(broom)
library(sommer)
library(heritability)

source(file.path(repo_dir, "source.R"))

## significance level
alpha <- 0.05
# fdr level
fdr_level <- 0.10

## Load the GBS/BOPA markers
bopa_snps <- load("C:/GoogleDrive/BarleyLab/Breeding/GenotypicData/Combined_GBS_BOPA_Genotype_Data/S2TP_multi_genos.RData")

## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))


# ## Convert markers for mapping
# ## GBS
# geno <- snp_info %>% 
#   rename(marker = 1) %>% 
#   select(marker, chrom, pos) %>%
#   as.data.frame() %>%
#   cbind(., t(s2_imputed_mat_use)) %>%
#   select(marker:pos, tp_geno)

## BOPA/GBS
geno <- select(s2tp_genos_imputed_hmp, -alleles, cM_pos)
# K matrix
K1 <- A.mat(X = s2tp_genos_imputed, min.MAF = 0, max.missing = 1)

## Reorganize the covariate coefficients
ec_interaction_coef <- covariate_reg_coefs %>%
  rename_all(~str_remove(., "_model")) %>%
  gather(model, out, -trait) %>% 
  unnest(out)

## Determine heritability
ec_interaction_coef_herit <- ec_interaction_coef %>%
  group_by(trait, model, covariate) %>%
  do({
    df <- .
    
    df1 <- subset(df, line_name %in% tp_geno)
    
    ## Calculate heritability
    invisible(capture.output(herit_fit <- marker_h2(data.vector = df1$estimate, geno.vector = df1$line_name, K = K1, alpha = alpha)))
    
    ## Return df
    tibble(heritability = herit_fit$h2, lower = herit_fit$conf.int1[1], upper = herit_fit$conf.int1[2])
    
  })


## Organize phenotypes for mapping
pheno <- ec_interaction_coef %>%
  unite(map_trait, c("trait", "covariate", "model"), sep = ":") %>%
  spread(map_trait, estimate) %>%
  as.data.frame() %>%
  filter(line_name %in% tp_geno)

K1_use <- K1[unique(pheno$line_name), unique(pheno$line_name)]

geno_use <- select(geno, marker, chrom, pos, unique(pheno$line_name))

# Map
map_out <- GWAS(pheno = pheno, geno = geno_use, K = K1_use, min.MAF = 0)



## Reorganize scores, determine fdr threshold
gwas_scores <- map_out %>% 
  gather(map_trait, marker_score, -marker:-pos) %>% 
  separate(map_trait, c("trait", "covariate", "model"), sep = "\\.") %>%
  group_by(trait, covariate, model) %>%
  mutate(fdr = sommer:::fdr(p = 10^-marker_score, fdr.level = fdr_level)$fdr.10,
         fdr = -log10(fdr)) %>%
  ungroup()


## Save
save("gwas_scores", file = file.path(result_dir, "ec_covariate_coef_gwas.RData"))


## Select a few examples
chrom_colors <- setNames(rep(c("black", "grey"), length.out = 7), seq(7))
# Chromosome position breaks
chrom_breaks <- seq(0, 800, by = 200)

gwas_scores %>%
  filter(model == "final") %>%
  mutate(chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos, y = marker_score, color = chrom)) +
  geom_hline(data = distinct(gwas_scores, trait, covariate, model, fdr), aes(yintercept = fdr), lty = 2) +
  geom_point() +
  facet_grid(trait + covariate ~ chrom, scales = "free_x", space = "free_x", switch = "both") +
  scale_color_manual(values = chrom_colors, guide = FALSE) +
  theme_presentation2(base_size = 10) +
  theme(panel.border = element_blank(), panel.spacing.x = unit(0, "line"))



g_gwas_grain <- gwas_scores %>%
  filter(model == "final") %>%
  filter(trait %in% c("GrainYield", "GrainProtein")) %>%
  mutate(chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos / 1e6, y = marker_score, color = chrom)) +
  geom_point() +
  facet_grid(trait + covariate ~ chrom, scales = "free_x", space = "free_x", switch = "both") +
  scale_color_manual(values = chrom_colors, guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression(-log[10]~'(p)')) +
  scale_x_continuous(breaks = chrom_breaks, name = "Position (Mbp)") +
  theme_presentation2(base_size = 10) +
  theme(panel.border = element_blank(), panel.spacing.x = unit(0.1, "line"),
        axis.text.x = element_text(angle = 45, hjust = 1))


## 



#




















# # Organize phenotypes for GWAS 
# pheno1 <- unnest(covariate_reg_coefs, apriori_model) %>%
#   filter(line_name %in% tp_geno)
# 
# # Marker data to use
# M <- s2_imputed_mat_use[tp_geno,, drop = T]
# 
# 
# # Group by trait and covariate and map
# sommer_gwas <- pheno1 %>%
#   group_by(trait, covariate) %>%
#   do({
#     df <- .
#     
#     ##
#     out <- sommer::GWAS(fixed = estimate ~ 1, random = ~ vs(line_name, Gu = K_use), 
#                         data = df, M = M, gTerm = "u:line_name", verbose = FALSE, 
#                         date.warning = FALSE)
#     
#     ## 
# 
#     marker_scores <- t(out$scores) %>% 
#       as.data.frame() %>% 
#       rownames_to_column("marker") %>%
#       rename_all(~str_remove_all(., "estimate "))
#   
#     
#     # Get variance components
#     var_comp <- setNames(out$sigma[,,], c("line_name", "Residuals"))
#     
#     tibble(marker_scores = list(marker_scores), var_comp = list(var_comp), fdr = fdr_level)
#     
#   })
# 
# 
# ## Get FDR level
# fdr_level <- 0.10
# sommer_gwas1 <- sommer_gwas %>%
#   ungroup() %>%
#   mutate(fdr = map_dbl(marker_scores, ~ sommer:::fdr(p = .$score, fdr.level = fdr_level)$fdr.10))
# 
# ## Plot - only those with significant QTL
# sommer_gwas %>%
#   unnest(marker_scores) %>%
#   inner_join(select(geno, marker, chrom, pos), .) %>%
#   ggplot(aes(x = pos, y = score)) +
#   geom_point() +
#   facet_grid(trait + covariate ~ chrom, scales = "free", space = "free_x") +
#   theme_manhattan()





































