## S2MET Prediction Models
## 
## GWAS of environment-specific response
## Genomewide prediction of responses
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
library(cowplot)

source(file.path(repo_dir, "source.R"))

## significance level
alpha <- 0.05
# fdr level
fdr_level <- 0.10

## Load environmental covariables
load(file = file.path(result_dir, "ec_model_building.RData"))








#####

# Organize data
s2_covariate_response <- S2_MET_BLUEs %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  group_by(trait, pop) %>%
  nest() %>% 
  crossing(., model = c("apriori_model", "final_model")) %>%
  mutate(response = list(NULL), mean = list(NULL))

## Fit models to determine genotypic response to covariates.
## Separate by population
for (i in seq(nrow(s2_covariate_response))) {
  
  # Get the df
  df <- s2_covariate_response$data[[i]]
  tr <- s2_covariate_response$trait[i]
  model_type <- s2_covariate_response$model[i]
  
  ## Get the relevant covariates for the trait
  ecs <- subset(covariate_reg_coefs, trait == tr, model_type, drop = TRUE)[[1]] %>%
    distinct(covariate) %>%
    pull(covariate)
  
  # Edit factors
  df1 <- df %>%
    mutate_at(vars(line_name, environment), as.factor) %>%
    mutate_at(vars(line_name, environment), droplevels) %>%
    mutate_at(vars(line_name, environment), ~`contrasts<-`(., value = `colnames<-`(contr.sum(levels(.)), head(levels(.), -1)))) %>%
    ## Add covariate information
    left_join(., select(ec_tomodel, environment, ecs), by = "environment")
  
  ## Refit the larger model
  ec_model_i <- lmer(formula = formula(subset(ec_model_building, trait == tr, model_type, drop = TRUE)[[1]]),
                     data = df1)

  ## Get the coefficients
  covariate_reg_coefs_i <- fixef(ec_model_i) %>%
    tibble(term = names(.), estimate = .) %>% 
    filter(str_detect(term, "line_name")) %>% 
    mutate(term = str_remove_all(term, "line_name")) %>% 
    separate(term, c("covariate", "line_name"), sep = ":") %>%
    spread(covariate, estimate) %>%
    add_row(line_name = tail(levels(df1$line_name), 1))
  
  covariate_reg_coefs_i[nrow(covariate_reg_coefs_i),-1] <- map_dbl(covariate_reg_coefs_i[-1], ~-sum(., na.rm = TRUE))
  
  ## Add intercepts
  covariate_intercepts <- fixef(ec_model_i)[ecs] %>% 
    t() %>% 
    as.data.frame()
  
  
  ## Compare to previous coefficients
  # subset(covariate_reg_coefs, trait == tr, final_model, drop = TRUE)[[1]] %>% 
  #   spread(covariate, estimate) %>% 
  #   left_join(., covariate_reg_coefs_i, by = "line_name") %>% 
  #   summarize(cor = cor(grain_fill_water_stress.x, grain_fill_water_stress.y))
  #   
  #   # Correlation of 1
  
  s2_covariate_response$response[[i]] <- covariate_reg_coefs_i
  s2_covariate_response$mean[[i]] <- covariate_intercepts
  
}


## Heritability

## Reorganize the covariate coefficients
s2_covariate_response1 <- s2_covariate_response %>%
  mutate(response = map(response, ~gather(., covariate, estimate, -line_name)),
         mean = map(mean, ~gather(., covariate, mean)),
         response = map2(response, mean, left_join, by = "covariate")) %>%
  unnest(response)



s2_covariate_heritability <- s2_covariate_response1 %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  group_by(trait, covariate, model, pop) %>%
  # group_by(trait, covariate) %>%
  do({
    df <- .
    
    ## Calculate heritability
    invisible(capture.output(herit_fit <- marker_h2(data.vector = df$estimate, geno.vector = df$line_name, K = K, alpha = alpha)))
    
    ## Return df
    tibble(heritability = herit_fit$h2, lower = herit_fit$conf.int1[1], upper = herit_fit$conf.int1[2])
     
  })  





## Plot distributions
s2_covariate_response1 %>%
  ggplot(aes(x = estimate)) +
  geom_histogram() +
  facet_wrap(~ pop + model + trait + covariate, scales = "free_x") +
  theme_presentation2(base_size = 10)
ggsave(filename = "all_covariate_response_histogram.jpg", path = fig_dir, height = 20, width = 20,
       dpi = 1000)

## Select relevant parameters and plot
s2_covariate_response_hist <- s2_covariate_response1 %>%
  filter(pop == "tp") %>%
  left_join(., s2_covariate_heritability) %>% # add heritability estimates
  mutate(response = estimate + mean) %>% # add mean to estimate
  group_by(trait, model, covariate) %>%
  do(plot = {
    df <- .
    # Get the new name and the units
    ec_name <- covariate_rename[unique(df$covariate)]
    ec_unit <- covariate_units[unique(df$covariate)]
    
    ## Combine a trait unit with ec unit
    unit_use <- parse(text = paste0(trait_units[unique(df$trait)], "~(", ec_unit, ")^-1"))
    
    ## Create a subtitle
    subtitle <- parse(text = paste0("Reaction~of~", str_replace(str_add_space(unique(df$trait)), " ", "~"), "~to~", ec_name))
    # Get the heritability
    h2 <- round(unique(df$heritability), 2)
    
    ## Plot, include heritability
    ggplot(data = df, aes(x = response)) +
      geom_histogram(bins = 20, fill = "grey85", color = "black") +
      annotate(geom = "text", x = -Inf, y = Inf, label = paste0("h^2==", h2), parse = TRUE,
                hjust = -0.3, vjust = 1) + 
      scale_x_continuous(breaks = pretty, name = unit_use) +
      scale_y_continuous(breaks = pretty, name = NULL) +
      labs(subtitle = subtitle) +
      theme_acs(base_size = 10) +
      theme(panel.border = element_blank(), axis.line = element_line())
    
  }) %>% ungroup()


## Select some plots of interest
# Designate an expression for filtering
# interest_exp <- '(trait == "GrainYield" & covariate == "grain_fill_water_stress") | (trait == "GrainProtein" & covariate == "grain_fill_maxt") | (trait == "HeadingDate" & covariate == "vegetative_tmean") '
interest_exp <- '(trait == "GrainYield" & covariate %in% c("flowering_mint", "grain_fill_water_stress") | (trait == "GrainProtein" & covariate == "grain_fill_maxt")) '

plot_list <- s2_covariate_response_hist %>%
  filter(model == "apriori_model") %>%
  filter(eval(parse(text  = interest_exp))) %>%
  pull(plot)

## Merge plots
covariate_plot_merge <- plot_grid(plotlist = plot_list, nrow = 1)

ggsave(filename = "select_covariate_response_histogram.jpg", plot = covariate_plot_merge,
       path = fig_dir, height = 3.5, width = 10, dpi = 1000)






## Plot reaction norms with heritability and histograms
s2_covariate_reaction_norms <- s2_covariate_response1 %>% 
  filter(pop == "tp") %>%
  filter(eval(parse(text = interest_exp))) %>%
  group_by(trait, covariate) %>%
  {
    bind_rows(top_n(x = ., n = 5, wt = estimate),
              top_n(x = ., n = 5, wt = -estimate))
  } %>% 
  ungroup() %>%
  left_join(., S2_MET_BLUEs) %>% 
  left_join(., gather(ec_tomodel, covariate, covariate_value, -environment))

s2_covariate_reaction_norm_plots <- s2_covariate_reaction_norms %>%
  group_by(trait, covariate) %>%
  do(plot = {
    df <- .
    
    # Get the new name and the units
    ec_name <- covariate_rename[unique(df$covariate)]
    ec_unit <- covariate_units[unique(df$covariate)]
    
    ## Combine a trait unit with ec unit
    trait_name <- str_add_space(unique(df$trait))
    trait_unit <- trait_units[unique(df$trait)]
    
    ## Create a subtitle
    subtitle <- parse(text = paste0("Reaction~of~", str_replace(str_add_space(unique(df$trait)), " ", "~"), "~to~", ec_name))
    y_title <- parse(text = paste0(str_replace(str_add_space(unique(df$trait)), " ", "~"), "~(", trait_unit, ")"))
    x_title <- parse(text = paste0(ec_name, "~(", ec_unit, ")"))
    
    ## Calculate line_means
    df1 <- df %>%
      group_by(line_name) %>%
      do(fit = lm(value ~ covariate_value, data = .)) %>%
      ungroup() %>% 
      mutate(geno_mean = map_dbl(fit, ~coef(.)[1])) %>%
      left_join(., df)
    
    ggplot(data = df1, aes(x = covariate_value, y = value, group = line_name)) +
      # geom_point() +
      # geom_smooth(method = "lm", se = FALSE) +
      geom_blank() +
      geom_abline(aes(intercept = geno_mean, slope = estimate + mean, group = line_name)) + 
      scale_x_continuous(breaks = pretty, name = x_title) +
      scale_y_continuous(breaks = pretty, name = y_title) +
      labs(subtitle = subtitle) +
      theme_acs(base_size = 10) +
      theme(panel.border = element_blank(), axis.line = element_line())
  })


plot_list <- s2_covariate_reaction_norm_plots$plot

## Merge plots
covariate_reaction_plot_merge <- plot_grid(plotlist = plot_list, nrow = 1)

ggsave(filename = "select_covariate_response_reaction_norm.jpg", plot = covariate_reaction_plot_merge,
       path = fig_dir, height = 3.5, width = 10, dpi = 1000)

## Merge reaction norms and histograms
s2_covariate_reaction_norm_hist_plots <- left_join(s2_covariate_reaction_norm_plots, s2_covariate_response_hist,
                                                   by = c("trait", "covariate")) %>%
  group_by(trait, covariate) %>%
  do(plot_new = {
    row <- .
    plot_grid(row$plot.x[[1]], row$plot.y[[1]] + labs(subtitle = NULL), ncol = 1, align = "hv")
  })


covariate_reaction_plot_merge1 <- plot_grid(plotlist = s2_covariate_reaction_norm_hist_plots$plot_new, nrow = 1)

ggsave(filename = "select_covariate_response_histogram_reaction_norm.jpg", plot = covariate_reaction_plot_merge1,
       path = fig_dir, height = 6, width = 10, dpi = 1000)




## Save the new responses
save("s2_covariate_response1", file = file.path(result_dir, "s2_ec_response.RData"))









###############################
### Genomewide prediction
###############################


## Cross-validation ##
# Subset K
K_cv <- K[tp_geno, tp_geno]

s2_covariate_cv <- s2_covariate_response1 %>%
  filter(line_name %in% c(tp_geno)) %>%
  group_by(trait, covariate, model) %>%
  do({
    df <- .
    
    # Refactor
    df1 <- mutate(df, line_name = factor(line_name, levels = tp_geno))
    
    # List lines
    cv_lines <- as.character(unique(df1$line_name))
    
    # Map over the lines
    cv_predictions <- map_dbl(cv_lines, ~{
      mf <- model.frame(estimate ~ line_name, data = df1, subset = line_name != .)
      y <- model.response(mf)
      Z <- model.matrix(~ -1 + line_name, mf)
      # Return prediction
      mixed.solve(y = y, Z = Z, K = K_cv)$u[.]
    })
    
    # Calculate accuracy
    tibble(line_name = cv_lines, prediction = cv_predictions) %>%
      inner_join(., df1, by = "line_name")
    
  }) %>% ungroup()


## Calculate acc


## Plot predicted versus observed responses
s2_covariate_cv_plot <- s2_covariate_cv %>%
  filter(eval(parse(text  = interest_exp))) %>%
  group_by(trait, covariate) %>%
  do(plot = {
    df <- .
    
    ## Combine a trait unit with ec unit
    trait_name <- str_add_space(unique(df$trait))
    trait_unit <- trait_units[unique(df$trait)]
    
    # Get the new name and the units
    ec_name <- covariate_rename[unique(df$covariate)]
    ec_unit <- covariate_units[unique(df$covariate)]
    
    ## Combine a trait unit with ec unit
    unit_use <- parse(text = paste0(trait_units[unique(df$trait)], "~(", ec_unit, ")^-1"))
    
    ## Calculate accuracy
    accuracy <- paste0("r==", round(cor(df$estimate, df$prediction), 2))
    
    df %>%
      mutate(group = paste0(str_replace_all(str_add_space(trait_name), " ", "~"), ":~", ec_name)) %>%
      ggplot(aes(x = prediction, y = estimate)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      annotate(geom = "text", x = -Inf, y = Inf, label = accuracy, parse = T, 
               hjust = -0.3, vjust = 1.5) + 
      scale_x_continuous(name = "Prediction", breaks = pretty) +
      scale_y_continuous(name = "Observation", breaks = pretty) + 
      facet_wrap(~ group, scales = "free", labeller = labeller(group = label_parsed)) +
      theme_presentation2(base_size = 12)
    
  }) %>% ungroup()


# Extract sample plots
plot_list <- s2_covariate_cv_plot %>%
  pull(plot) %>%
  map(~ . + theme(axis.title = element_blank()))

# Merge
select_ec_plot_merge <- plot_grid(plotlist = plot_list, nrow = 1)

# y_axis label
y_label <- grid::textGrob(label = "Observation", rot = 90) 
x_label <- grid::textGrob(label = "Prediction")

select_ec_plot_merge1 <- plot_grid(y_label, select_ec_plot_merge, nrow = 1, rel_widths = c(0.05, 1))  %>%
  plot_grid(., x_label, ncol = 1, rel_heights = c(1, 0.05)) 


ggsave(filename = "sample_covariate_cv.jpg", plot = select_ec_plot_merge1, path = fig_dir,
       height = 4, width = 10, dpi = 1000)







## Plot
s2_covariate_cv %>%
  ggplot(aes(x = "acc", y = accuracy, fill = model)) +
  geom_col(position = "dodge") +
  facet_grid(~trait + covariate)


## Parnet-offspring validation ##
# Subset K
K_pov <- K

s2_covariate_pov <- s2_covariate_response1 %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  group_by(trait, covariate, model) %>%
  do({
    df <- .
    
    # Refactor
    df1 <- mutate(df, line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
    
    # Fit the model
    mf <- model.frame(estimate ~ line_name, data = df1, subset = ! line_name %in% vp_geno)
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name, mf)
    # Return prediction
    fit <- mixed.solve(y = y, Z = Z, K = K_pov)
    
    # Calculate accuracy
    fit$u %>%
      tibble(line_name = names(.), prediction = .) %>%
      inner_join(., df1, by = "line_name") %>%
      filter(line_name %in% vp_geno)
    
  })


## Plot predicted versus observed responses
s2_covariate_pov_plot <- s2_covariate_pov %>%
  filter(eval(parse(text  = interest_exp))) %>%
  group_by(trait, covariate) %>%
  do(plot = {
    df <- .
    
    ## Combine a trait unit with ec unit
    trait_name <- str_add_space(unique(df$trait))
    trait_unit <- trait_units[unique(df$trait)]
    
    # Get the new name and the units
    ec_name <- covariate_rename[unique(df$covariate)]
    ec_unit <- covariate_units[unique(df$covariate)]
    
    ## Combine a trait unit with ec unit
    unit_use <- parse(text = paste0(trait_units[unique(df$trait)], "~(", ec_unit, ")^-1"))
    
    ## Calculate accuracy
    accuracy <- paste0("r==", round(cor(df$estimate, df$prediction), 2))
    
    df %>%
      mutate(group = paste0(str_replace_all(str_add_space(trait_name), " ", "~"), ":~", ec_name)) %>%
      ggplot(aes(x = prediction, y = estimate)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      annotate(geom = "text", x = -Inf, y = Inf, label = accuracy, parse = T, 
               hjust = -0.3, vjust = 1.5) + 
      scale_x_continuous(name = "Prediction", breaks = pretty) +
      scale_y_continuous(name = "Observation", breaks = pretty) + 
      facet_wrap(~ group, scales = "free", labeller = labeller(group = label_parsed)) +
      theme_presentation2(base_size = 12)
    
  }) %>% ungroup()


# Extract sample plots
plot_list <- s2_covariate_pov_plot %>%
  pull(plot) %>%
  map(~ . + theme(axis.title = element_blank()))

# Merge
select_ec_plot_merge <- plot_grid(plotlist = plot_list, nrow = 1)

# y_axis label
y_label <- grid::textGrob(label = "Observation", rot = 90) 
x_label <- grid::textGrob(label = "Prediction")

select_ec_plot_merge1 <- plot_grid(y_label, select_ec_plot_merge, nrow = 1, rel_widths = c(0.05, 1))  %>%
  plot_grid(., x_label, ncol = 1, rel_heights = c(1, 0.05)) 


ggsave(filename = "sample_covariate_pov.jpg", plot = select_ec_plot_merge1, path = fig_dir,
       height = 4, width = 10, dpi = 1000)




## Convert markers for mapping
## GBS
geno <- snp_info %>%
  rename(marker = 1) %>%
  select(marker, chrom, pos) %>%
  as.data.frame() %>%
  cbind(., t(s2_imputed_mat_use)) %>%
  select(marker:pos, tp_geno)

K1 <- K


# ## BOPA/GBS
# geno <- select(s2tp_genos_imputed_hmp, -alleles, cM_pos)
# # K matrix
# K1 <- A.mat(X = s2tp_genos_imputed, min.MAF = 0, max.missing = 1)
# 



## Organize phenotypes for mapping
pheno <- s2_covariate_response1 %>%
  unite(map_trait, c("trait", "covariate", "model"), sep = ":") %>%
  spread(map_trait, estimate) %>%
  as.data.frame() %>%
  filter(line_name %in% tp_geno)

K1_use <- K1[unique(pheno$line_name), unique(pheno$line_name)]

geno_use <- select(geno, marker, chrom, pos, unique(pheno$line_name))


## SS-GWAS
ss_geno <- geno_use %>%
  select(-marker:-pos) %>%
  t()
  

ss_map_out <- s2_covariate_response1 %>%
  filter(line_name %in% tp_geno) %>%
  mutate(line_name = factor(line_name, levels = tp_geno)) %>%
  group_by(trait, covariate, model) %>%
  do({
    df <- .
    
    # print(unique(df$covariate))
    
    ## Model matrices
    mf <- model.frame(estimate ~ line_name, data = df)
    y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    Z <- ss_geno
    
    ## EMMREML to estimate marker effects
    fit <- EMMREML::emmreml(y = y, X = X, Z = Z, K = diag(ncol(Z)), varbetahat = TRUE, varuhat = TRUE, PEVuhat = TRUE, test = TRUE)
    
    ## Calculate scores for each marker
    stat <- c(fit$uhat) / sqrt(fit$varuhat)
    df <- length(y) - 1
    pvalue <- 2 * pt(q = abs(stat), df = df, lower.tail = FALSE)
    score <- -log10(pvalue)
    
    ## Output the scores
    snp_info %>%
      mutate(score = score) %>%
      rename(marker = 1)
    
  }) %>% ungroup()


## reorganize, add fdr
ss_map_out1 <- ss_map_out %>%
  rename(marker_score = score) %>%
  group_by(trait, covariate, model) %>%
  mutate(fdr = sommer:::fdr(p = 10^-marker_score, fdr.level = fdr_level)$fdr.10,
         fdr = -log10(fdr)) %>%
  ungroup() 



# ## Reorganize scores, determine fdr threshold
# gwas_scores <- map_out %>% 
#   gather(map_trait, marker_score, -marker:-pos) %>% 
#   # separate(map_trait, c("trait", "covariate", "model"), sep = "\\.") %>%
#   separate(map_trait, c("trait", "covariate"), sep = "\\.") %>%
#   # group_by(trait, covariate, model) %>%
#   group_by(trait, covariate) %>%
#   mutate(fdr = sommer:::fdr(p = 10^-marker_score, fdr.level = fdr_level)$fdr.10,
#          fdr = -log10(fdr)) %>%
#   ungroup() %>%
#   left_join(snp_info, .) %>%
#   arrange(trait, covariate, chrom, pos)
# 
# 
# ## Compare pvalues
# inner_join(ss_map_out1, gwas_scores, by = c("trait", "covariate", "marker", "allele", "chrom", "pos", "cM_pos")) %>%
#   group_by(trait, covariate) %>% 
#   summarize(compare = cor(marker_score.x, marker_score.y))
# 



## Save
save("ss_map_out1", file = file.path(result_dir, "ec_covariate_coef_gwas.RData"))


## Select a few examples
chrom_colors <- setNames(rep(c("black", "grey"), length.out = 7), seq(7))
# Chromosome position breaks
chrom_breaks <- seq(0, 800, by = 200)

gwas_scores %>%
  # filter(model == "final") %>%
  mutate(chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos, y = marker_score, color = chrom)) +
  geom_hline(data = distinct(gwas_scores, trait, covariate, model, fdr), aes(yintercept = fdr), lty = 2) +
  geom_point() +
  facet_grid(trait + covariate ~ chrom, scales = "free_x", space = "free_x", switch = "both") +
  scale_color_manual(values = chrom_colors, guide = FALSE) +
  theme_presentation2(base_size = 10) +
  theme(panel.border = element_blank(), panel.spacing.x = unit(0, "line"))


## Which trait/covariate/model have an FDR threshold
ss_map_out1 %>%
  filter(!is.na(fdr)) %>%
  distinct(trait, covariate, model)



# g_gwas_grain <- gwas_scores %>%
g_gwas_select <- ss_map_out1 %>%
  filter(eval(parse(text  = interest_exp))) %>%
  mutate(chrom = as.factor(chrom),
         covariate = str_replace_all(covariate, covariate_rename),
         group = paste0(str_replace_all(str_add_space(trait), " ", "~"), "~-~", covariate)) %>%
  {
    ggplot(data = ., aes(x = pos / 1e6, y = marker_score, color = chrom)) +
      # ggplot(aes(x = cM_pos, y = marker_score, color = chrom)) +
      geom_point() +
      geom_hline(data = distinct(., trait, covariate, model, fdr), aes(yintercept = fdr), lty = 2) +
      facet_grid(group ~ chrom, scales = "free_x", space = "free_x", switch = "both",
                 labeller = labeller(group = label_parsed) ) +
      scale_color_manual(values = chrom_colors, guide = FALSE) +
      scale_y_continuous(breaks = pretty, name = expression(-log[10]~'(p)'), limits = c(0, 5)) +
      scale_x_continuous(breaks = chrom_breaks, name = "Position (Mbp)") +
      theme_presentation2(base_size = 16) +
      theme(panel.border = element_blank(), panel.spacing.x = unit(0.1, "line"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
  }

## Save
ggsave(filename = "select_covariate_response_gwas.jpg", plot = g_gwas_select, path = fig_dir,
       height = 8, width = 10, dpi = 1000)


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





































