## S2MET Predictions Models
## 
## Scripts to generate manuscript figures
## 
## 

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(broom)
library(patchwork)
library(paletteer)
library(cowplot)
library(ggrepel)
library(ggsn)
library(lubridate)
library(ggdendro)

# Set subfigure labels
subfigure_labels <- LETTERS

# Set resolution
dpi_use <- 2000


# Define heat colors
heat_colors <- wesanderson::wes_palette("Zissou1", n = 5)[c(1,3,5)]



# Figure 1: map, crop model output, generation of covariates --------------


## Create the map with number of environments per trait ##

# Compute number of environments and locations per trait
trait_experiment_summary <- S2_MET_BLUEs %>% 
  group_by(trait) %>% 
  summarize_at(vars(environment, location), n_distinct) %>%
  mutate(trait = str_add_space(trait)) %>%
  rename_all(str_to_title) %>%
  rename_at(vars(-Trait), ~paste0(., "s"))



# Get the map data for canada
canada <- map_data("world", "Canada")

# Download map data for US by county
usa_county <- map_data(map = "county")
# Download state data
usa_state <- map_data(map = "state")

# Adjust the groups in the states
usa_state <- usa_state %>%
  mutate(group = group + max(canada$group))

# Adjust the groups in the counties
usa_county <- usa_county %>%
  mutate(group = group + max(usa_state$group))

# Tidy and combine
north_america <- bind_rows(usa_state, usa_county, canada)


## Now create labels for each location (as a combination of all environments)
use_loc_info_toplot <-  trial_info %>% 
  # Correct GRV spelling mistake
  mutate(location = ifelse(location == "Grande_rhonde_valley", "Grande_ronde_valley", location)) %>%
  # Add trial designator
  mutate(set = ifelse(environment %in% train_test_env, "Train/test locations", "Holdout locations"),
         set = fct_relevel(set, "Holdout locations", after = Inf)) %>%
  filter(environment %in% unique(S2_MET_BLUEs$environment)) %>%
  group_by(location, set, latitude, longitude) %>%
  summarize(n_years = n_distinct(year)) %>%
  ungroup() %>%
  mutate(n_years = factor(n_years, levels = sort(unique(n_years))),
         location = str_to_title(str_replace_all(location, "_", " ")))


## Collapse Ithaca
use_loc_info_toplot1 <- use_loc_info_toplot %>% 
  mutate(location = str_remove(location, "[0-9]{1}"), 
         n_years = parse_number(as.character(n_years))) %>% 
  group_by(location) %>% 
  mutate(nExperiments = sum(n_years)) %>% 
  slice(1) %>%
  ungroup()


# Coordinate limits
long_limit <- c(-117, -63)
lat_limit <- c(39, 49)

## Different version of the map - grey
g_map_alt <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "grey85") +
  geom_polygon(data = canada, fill = NA, color = "white", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "white", lwd = 0.3) +
  geom_point(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, shape = set), size = 3.5) +
  geom_text_repel(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = location),
                  size = 2, hjust = 0.5, nudge_x = ifelse(use_loc_info_toplot1$location == "Bozeman", 3, -1), segment.size = 0.2, 
                  point.padding = unit(2, "pt"), min.segment.length = 1) +
  geom_text(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = nExperiments), size = 2, 
            color = ifelse(use_loc_info_toplot1$location == "Arlington", "black", "white")) +
  coord_fixed(ratio = 1.5, xlim = long_limit, ylim = lat_limit) +
  scale_shape_manual(name = NULL, values = c(16, 15), guide = guide_legend(label.position = "left", override.aes = list(size = 2))) +
  scale_x_continuous(breaks = NULL, name = NULL, labels = NULL) + 
  scale_y_continuous(breaks = NULL, name = NULL, labels = NULL) +
  theme_classic(base_size = 8) +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = alpha("white", 0)), axis.line = element_blank(),
        legend.position = c(0.90, 0.20), legend.text.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.5, "lines"), legend.key.width = unit(0.25, "lines"))


## Add a simple bar chart of the number of environments/locations per trait

## Determine the number of environments per trait
trait_env_breakdown <- S2_MET_BLUEs %>% 
  # Assign environment by train/test or external and split
  mutate(env_set = ifelse(environment %in% train_test_env, "Train/test", "Holdout"),
         env_set = as.factor(env_set)) %>%
  arrange(trait, env_set) %>%
  group_by(trait, env_set) %>%
  summarize_at(vars(location, year, environment), n_distinct) %>%
  mutate(y = ifelse(env_set == "Holdout", last(environment) + (environment / 2), environment / 2))

# Plot
g_trait_env_count <- trait_env_breakdown %>%
  ggplot(aes(x = trait, y = environment, fill = env_set)) +
  geom_col() +
  geom_text(aes(y = y, label = environment), size = 2) +
  scale_fill_manual(name = NULL, values = neyhart_palette("umn2")[c(8, 9)]) +
  scale_x_discrete(labels = function(x) abbreviate(str_add_space(x), 2), name = NULL) +
  scale_y_continuous(breaks = pretty, name = "Environments") +
  theme_genetics(base_size = 6) +
  theme(legend.position = c(0.10, 0.90))




## Plot examples of crop model outputs

# Load concurrent and historical crop model output
load(file.path(enviro_dir, "GrowthStaging/apsim_s2met_model_results_daymet.RData"))
load(file.path(enviro_dir, "GrowthStaging/apsim_historical_growth_model_results.RData"))


# Assign units to variables
covariate_rename_abbr <- c("mint" = "T[min]", "maxt" = "T[max]", "tmean" = "T[mean]", "gdd" = "GDD",
                           "water_balance" = "W", "radn" = "R")
                           
covariate_rename <- c("mint" = "Min. temperature", "maxt" = "Max. temperature", "tmean" = "Mean temperature",
                      "gdd" = "Growing deg. days", "water_balance" = "Water balance", "radn" = "Solar radiation")

# Units
covariate_variable_unit <- setNames(c(rep("degree*C", 4), "mm", "MJ~m^-2"), names(covariate_rename))

## Function for renaming a covariate
f_covariate_replace <- function(x) paste0("'", covariate_rename[x], "'~(", covariate_variable_unit[x], ")")
# Function for renaming timeframe
f_timeframe_replace <- function(x) c("concurrent" = "Environment-concurrent", "historical" = "Historical location average")[x]


## Plot Bozeman and Columbus (Wooster)
cgm_locations_plot <- c("Bozeman", "Wooster")
cgm_year_plot <- 2017
# Assign colors to location
cgm_locations_color <- setNames(paletteer_d("ggsci", palette = "nrc_npg")[c(1,3)], cgm_locations_plot)

## Replacements and colors for growth stages

growth_stage_color <- setNames(object = paletteer_d("ggsci", palette = "default_jco", n = 5),
                               nm = c("early_vegetative", "late_vegetative", "heading", "flowering", "grain_fill"))


concurrent_cgm_toplot <- apsim_s2met_out %>%
  mutate(location = ifelse(location == "Columbus", "Wooster", location)) %>%
  filter(location %in% cgm_locations_plot, year == cgm_year_plot)
historical_cgm_toplot <- location_historical_growth_staging_daymet %>%
  mutate(location = ifelse(location == "Columbus", "Wooster", location)) %>%
  filter(location %in% cgm_locations_plot) %>%
  mutate(predicted_planting_date = parse_date(predicted_planting_date, format = "%Y-%m-%d"),
         data = map2(daily_data, growth_model, ~left_join(.x, select(.y, yday = day, pet = eo), by = "yday")))

## Combine
cgm_toplot <- bind_rows(
  select(concurrent_cgm_toplot, location, planting_date, latitude, longitude, data, growth_model = apsim_out) %>% 
    mutate(timeframe = "concurrent"),
  select(historical_cgm_toplot, location, planting_date = predicted_planting_date, latitude, longitude, data, growth_model) %>% 
    mutate(timeframe = "historical")
) %>% mutate(year = year(planting_date))


## Assign growth stage
# Use the results of the APSIM models to determine growth stages
# Z10 - Z30: early vegetative
# Z30 - Z49: late vegetative
# Z50 - Z59: heading
# Z60 - Z69: flowering
# Z70 - Z91: grain fill
# 

cgm_toplot_growth_stages <- cgm_toplot %>%
  mutate(growth_stages = map(growth_model, ~{
    
    df <- .
    
    ## Determine stages
    stages <- mutate(df, stage = case_when(
      between(zadok_stage, 10, 30) ~ "early_vegetative",
      between(zadok_stage, 30, 50) ~ "late_vegetative",
      # between(zadok_stage, 50, 60) ~ "heading",
      between(zadok_stage, 50, 70) ~ "flowering",
      between(zadok_stage, 70, 91) ~ "grain_fill"),
      stage = fct_inorder(stage))
    
    
    ## Return a df with date, dap, growth stage
    stages %>% 
      filter(sowing_das == 1) %>% 
      mutate(dap = seq(nrow(.))) %>% 
      select(date = Date, day, dap, stage)
    
  }))




## Add daily weather observations during for each day during a growth stage
cgm_toplot_growth_stage_weather <- cgm_toplot_growth_stages %>%
  unnest(growth_stages) %>%
  left_join(., unnest(cgm_toplot_growth_stages, data)) %>%
  # Add GDD info
  left_join(., unnest(cgm_toplot_growth_stages, growth_model) %>% select(location, year, timeframe, date = Date, tt = TT)) %>%
  ## Calculate water stress as the difference between pet and rain
  mutate(water_balance = rain - pet) %>%
  rename(gdd = tt) %>%
  ## Gather
  gather(variable, value, mint:water_balance) %>%
  # Skip NA stages
  filter(!is.na(stage))

## Create bars for growth stages
cgm_toplot_growth_stage_times <- cgm_toplot_growth_stage_weather %>% 
  group_by(location, timeframe, year, stage) %>% 
  filter(!is.na(stage)) %>% 
  summarize_at(vars(dap), list(~min, ~max)) %>%
  ungroup() %>%
  mutate(y = as.numeric(as.factor(location)), dap = min,
         type = "Growth\nstage") # Dummy variable for growth stage facet


## Plot 2-3 covariates
covariates_plot <- c("mint", "radn", "water_balance")


## Plot both concurrent and historical data ##
# First make some edits and then summarize over years
cgm_toplot_growth_stage_weather_historical <- cgm_toplot_growth_stage_weather %>%
  filter(variable %in% covariates_plot, timeframe == "historical", 
         # 5 years of data
         between(year, cgm_year_plot - 5, cgm_year_plot - 1)) %>%
  mutate(group = paste0(location, "_", year))

cgm_toplot_growth_stage_weather_historical_mean <- cgm_toplot_growth_stage_weather_historical %>% 
  group_by(timeframe, location, dap, variable) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>% 
  mutate(group = location) %>%
  ungroup()


## Combine to plot
cgm_toplot_growth_stage_weather_summary <- bind_rows(
  filter(cgm_toplot_growth_stage_weather, timeframe == "concurrent") %>% mutate(group = location),
  cgm_toplot_growth_stage_weather_historical_mean
) %>% filter(variable %in% covariates_plot) %>%
  mutate(variable = f_covariate_replace(variable))

# Simple subplot label df
g_weather_subplot_label <- tibble(timeframe = c("concurrent", "historical"), variable = f_covariate_replace(covariates_plot[1]), 
                                  label = subfigure_labels[2:3], dap = min(cgm_toplot_growth_stage_weather_historical$dap))


## Plot of both concurrent and historical mean (plus individual year) data
g_weather_concurrent_historical <- cgm_toplot_growth_stage_weather_historical %>%
  filter(variable %in% covariates_plot) %>%
  mutate(variable = f_covariate_replace(variable)) %>%
  ggplot(aes(x = dap, y = value, color = location, group = group)) +
  geom_line(alpha = 0.2, lwd = 0.15) +
  # Add average line
  geom_line(data = cgm_toplot_growth_stage_weather_summary, lwd = 0.3) +
  # Add text for plot labels
  geom_text(data = tibble(location = names(cgm_locations_color), variable = f_covariate_replace("mint"),
                          timeframe = "concurrent", x = 10, y = c(23, 20)),
            aes(x = x, y = y, color = location, label = location), inherit.aes = FALSE, hjust = 0, size = 2) +
  facet_grid(variable ~ timeframe, scales = "free", space = "free_x", switch = "y", 
             labeller = labeller(variable = label_parsed, timeframe = f_timeframe_replace)) +
  scale_y_continuous(name = NULL, breaks = pretty) +
  scale_x_continuous(name = "Days after planting date", breaks = pretty) +
  scale_color_manual(values = cgm_locations_color, name = NULL) +
  theme_genetics(base_size = 6) +
  theme(strip.placement = "outside", strip.background = element_blank(), panel.spacing.y = unit(0.5, "line"),
        legend.position = "none")





## Create summary data.frames for historical and concurrent growth stages
# First summarize the entire length of the growing season per location
cgm_toplot_growing_season <- cgm_toplot_growth_stage_times %>%
  split(.$timeframe) %>%
  modify_at(., "historical", ~filter(., between(year, cgm_year_plot - 5, cgm_year_plot - 1))) %>%
  bind_rows() %>%
  group_by(type, timeframe, location) %>% 
  summarize(min = min(min), max = max(max), y = mean(y)) %>%
  ungroup()


## Calculate average growth stages
cgm_toplot_growth_stage_times_summary <- cgm_toplot_growth_stage_times %>%
  filter(timeframe == "historical", between(year, cgm_year_plot - 5, cgm_year_plot - 1)) %>%
  group_by(type, timeframe, location, stage) %>%
  summarize_at(vars(min, max, y, dap), mean) %>%
  ungroup() %>%
  bind_rows(., filter(cgm_toplot_growth_stage_times, timeframe == "concurrent"))


# Plot segments for growth stages
# both concurrent and historical
g_growth_stages <- cgm_toplot_growth_stage_times %>%
  filter(timeframe == "historical", 
         # 3 years of data
         between(year, cgm_year_plot - 3, cgm_year_plot - 1)) %>%
  ggplot(aes(x = dap, color = stage)) +
  # Individual-year stages
  geom_segment(aes(x = min, xend = max, y = y, yend = y), lwd = 4, alpha = 0.05) +
  # Average stages
  geom_segment(data = cgm_toplot_growth_stage_times_summary, 
               aes(x = min, xend = max, y = y, yend = y), lwd = 2) +
  # Lines for the entire growing season per location
  geom_segment(data = cgm_toplot_growing_season, aes(x = min, xend = max, y = y, yend = y),
               color = cgm_locations_color[cgm_toplot_growing_season$location], lwd = 0.5, inherit.aes = FALSE) +
  # Annotation for growth stages
  geom_text(data = subset(cgm_toplot_growth_stage_times_summary, y == 1),
            aes(x = (min + max) / 2, y = y - 1, label = f_growth_stage_replace(as.character(stage))), size = 1.5) +
  facet_grid(type ~ timeframe, labeller = labeller(type = str_add_space, timeframe = f_timeframe_replace), 
             switch = "y", scales = "free_x", space = "free_x") +
  scale_color_manual(values = growth_stage_color, guide = FALSE) +
  scale_y_continuous(limits = c(-0.5, 2.5)) +
  scale_x_continuous(breaks = pretty, name = "Days after planting date") +
  theme_genetics(base_size = 6) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank(),
        axis.title.y = element_blank(), axis.line.y = element_blank())


## Combine weather and growth stage plot
# g_weather_stages <- plot_grid(
#   g_weather_concurrent_historical + theme(legend.position = "none", axis.title.x = element_blank(),
#                                           axis.text.x = element_blank()), 
#   g_growth_stages + theme(strip.text.x = element_blank()),
#   ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.25), labels = subfigure_labels[2:3], label_size = 8)

## Combine weather and growth stage plot
g_weather_stages <- plot_grid(
  g_growth_stages + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
  g_weather_concurrent_historical + theme(legend.position = "none", strip.text.x = element_blank()), 
  ncol = 1, align = "v", axis = "lr", rel_heights = c(0.20, 1), labels = subfigure_labels[2:3], label_size = 8,
  label_y = c(1, 1.03), label_x = 0.01)



## Highlight locations with color in the map
map_location_colors <- setdiff(g_map_alt$layers[[4]]$data$location, names(cgm_locations_color)) %>%
  setNames(rep("black", length(.)), .) %>%
  c(., cgm_locations_color)
# New map plot
g_map_alt1 <- g_map_alt
g_map_alt1$layers <- g_map_alt1$layers[-5]
g_map_alt1 <- g_map_alt1 +
  geom_text_repel(data = use_loc_info_toplot1, 
                  aes(x = longitude, y = latitude, group = location, label = location, color = location),
                  size = 2, hjust = 0.5, nudge_x = ifelse(use_loc_info_toplot1$location == "Bozeman", 3, -1), segment.size = 0.2, 
                  point.padding = unit(2, "pt"), min.segment.length = 1) +
  scale_color_manual(values = map_location_colors, guide = FALSE) +
  theme(legend.position = c(0.88, 0.20))
  
## Save map and environment count separately
# Save the map
ggsave(filename = "figure1_partA_draft1.jpg", plot = g_map_alt1, path = fig_dir,
       width = 8.7, height = 0.35 * 11, unit = "cm", dpi = dpi_use)

ggsave(filename = "figure1_partB_draft1.jpg", plot = g_trait_env_count, path = fig_dir,
       width = 2.7, height = 0.35 * 11, unit = "cm", dpi = dpi_use)



g_fig1 <- plot_grid(g_map_alt1, g_weather_stages, ncol = 1, labels = subfigure_labels[1], label_size = 8,
                    label_y = 0.93, label_x = 0.01, rel_heights = c(0.35, 1))

# Save
ggsave(filename = "figure1_map_growth_stage_example.jpg", plot = g_fig1, path = fig_dir, 
       width = 8.7, height = 11, units = "cm", dpi = dpi_use)













 # Figure 2. LOEO prediction accuracy ---------------------------------------


# Load the compiled prediction results
load(file.path(result_dir, "prediction_accuracy_compiled.RData"))


# Color scheme for models
model_colors <- c(neyhart_palette("umn1")[1], rep(c(rep(neyhart_palette("umn1")[3], 2), rep(neyhart_palette("umn1")[4], 3)), 2))
names(model_colors) <- names(model_replace)

# Vector of models to report
models_use <- str_subset(string = names(model_colors), pattern = "cov1", negate = TRUE)


## Part A - plot predicted/observed phenotypes for each trait and by target population
## for model 3 with stepwise covariates - LOEO

loo_pred_obs_df <- predictions_df %>%
  filter(type == "loeo") %>%
  filter(model %in% str_subset(models_use, "3|5"), selection == "rfa_cv_adhoc") %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  # mutate(trait = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"))
  mutate(trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))

# Split by trait and Plot
g_loo_pred_obs_list <- loo_pred_obs_df %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.2) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


# Combine the plots
g_loo_pred_obs <- g_loo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = 6)), ., rel_widths = c(0.03, 1))


# Save
ggsave(filename = "figure2_partA_draft.jpg", plot = g_loo_pred_obs, path = fig_dir,
       height = 120, width = 80, units = "mm", dpi = dpi_use)


## Part B - accuracy within environments - LOEO

## Filter for relevant cases
accuracy_within_env_toplot <- within_environment_prediction_accuracy %>%
  filter(type == "loeo") %>%
  filter(model %in% str_subset(models_use, "id", negate = TRUE), selection %in% c("none", "rfa_cv_adhoc")) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), "\n(n = ", nGroup+1, ")")) %>%
  unite(group, trait, model, pop, remove = FALSE)
  
accuracy_within_env_toplot2 <- accuracy_within_env_toplot %>%
  group_by(type, group, trait, trait1, model, pop, selection) %>%
  do({
    df <- .
    # Get ranges from boxplot
    bp <- boxplot(x = df$accuracy, plot = FALSE)
    # Return a tibble
    bind_cols(summarize(df, accuracy_mean = mean(accuracy)), 
              as_tibble(setNames(object = as.list(bp$stats), nm = c("lower", "q25", "median", "q75", "upper"))) )
  }) %>% ungroup()



# Plot point and line range
g_accuracy_within_env <- accuracy_within_env_toplot %>%
  ggplot(aes(x = trait1, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = accuracy), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = accuracy_within_env_toplot2, aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = accuracy_within_env_toplot2, aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = accuracy_within_env_toplot2, aes(y = accuracy_mean), position = position_dodge(0.9), size = 0.8) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = , breaks = pretty) +
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)) +
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)) +
  facet_grid(type ~ pop, switch = "y", labeller = labeller(pop = f_validation_replace)) +
  theme_genetics(base_size = 6) +
  theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
        legend.text.align = 1, legend.position = c(0.89, 0.95), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line")) +
  theme(axis.title = element_blank(), strip.text.y = element_text(color = "white"))

# Add text grob
g_accuracy_within_env1 <- plot_grid(textGrob(label = expression('Within-environment'~italic(r[MG])), 
                                             rot = 90, gp = gpar(fontsize = 6)), 
                                    g_accuracy_within_env, rel_widths = c(0.03, 1))  

## Save
ggsave(filename = "figure2_partB_draft.jpg", plot = g_accuracy_within_env1, path = fig_dir,
       width = 80, height = 40, dpi = dpi_use, units = "mm")

## Combine plots
g_figure2 <- plot_grid(g_loo_pred_obs, g_accuracy_within_env1, 
                       ncol = 1, rel_heights = c(1, 0.35), labels = subfigure_labels[1:2], label_size = 8)

## Save
ggsave(filename = "figure2_loeo_predictions.jpg", plot = g_figure2, path = fig_dir,
       width = 8.0, height = 14, dpi = dpi_use, units = "cm")


## Toy example of eliminating undesirable lines/environments based on GPC
loo_pred_obs_df %>%
  filter(trait == "GrainProtein", pop == "tp") %>%
  # Calculate number of environments/lines
  mutate(nEnvAll = n_distinct(leave_one_group), nLineAll = n_distinct(line_name),
         nObsAll = n()) %>%
  filter(pred_complete <= 12) %>%
  mutate(nEnvSel = n_distinct(leave_one_group), nLineSel = n_distinct(line_name), nObsSel = n()) %>%
  distinct_at(vars(contains("All"), contains("Sel")))






# Figure 3. LOLO predictions ---------------------------------

## Parts A and B: timeframe search results

# Load the search results
load(file.path(result_dir, "historical_covariate_timeframe_selection.RData"))
# Load the feature selection results
load(file.path(result_dir, "feature_selection_results.RData"))

# Extract results
historical_timeframe_selection_out1 <- historical_timeframe_selection_out %>%
  mutate(stepwise_results = map(adhoc, "finalResults") %>% map(as.list) %>% map(as_tibble)) %>%
  unnest(stepwise_results) %>%
  # Parse timeframe
  mutate(time_frame_type = str_extract(time_frame, "time_frame|window"),
         time_frame1 = str_remove(time_frame, "time_frame|window")) %>%
  separate(time_frame1, c("length", "start_year", "end_year"), sep = "_") %>%
  mutate_at(vars(length, contains("year")), parse_guess) %>%
  # Scale the RMSE for better plotting
  group_by(trait, model) %>%
  mutate(RMSE_scale = RMSE / max(RMSE)) %>%
  ungroup()


# Colors for traits
trait_colors <- setNames(neyhart_palette("umn2", 5), sort(traits))

# Get a range for scaling by trait
y_trait_range <- historical_timeframe_selection_out1 %>% 
  filter(model == "model3") %>%
  group_by(trait) %>%
  do(breaks = pretty(.$RMSE_scale, 3)) %>%
  ungroup()

# Plot modifier
g_mod <- list(
  scale_color_manual(values = trait_colors, guide = FALSE),
  facet_grid(trait ~ ., switch = "y", scales = "free_y", 
             labeller = labeller(trait = function(x) abbreviate(str_add_space(x), 2))),
  theme_genetics(base_size = 6),
  theme(strip.placement = "outside", strip.background = element_blank())
)


## Plot timeframe results
historical_timeframe_selection_out1_toplot <- historical_timeframe_selection_out1 %>%
  filter(time_frame_type == "time_frame", model == "model3")

# Separate df for the selected timeframe (to add as points and an annotation)
timeframe_used_df <- inner_join(distinct(historical_feature_selection, trait, time_frame), historical_timeframe_selection_out1_toplot)

# Plot together
g_timeframe_analysis <- historical_timeframe_selection_out1_toplot %>%
  ggplot(aes(x = length, y = RMSE_scale, color = trait)) +
  geom_line(lwd = 0.25) +
  geom_point(data = timeframe_used_df, size = 0.3) +
  # Add an annotation
  annotate(geom = "curve", x = 10, y = 0.35, xend = min(timeframe_used_df$length) + 1, yend = min(timeframe_used_df$RMSE_scale),
           curvature = 0.25, arrow = arrow(angle = 30, length = unit(0.25, "line"), ends = "last"), lwd = 0.25) +
  annotate(geom = "text", label = "Chosen timeframe", x = 10.5, y = 0.39, hjust = 1, size = 2) +
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Scaled RMSEP") +
  scale_x_continuous(name = "Historical data years (before 2015)", trans = "reverse") +
  scale_color_manual(values = trait_colors, name = NULL, labels = function(x) abbreviate(str_add_space(x), 2),
                     guide = guide_legend(nrow = 2)) +
  theme_genetics(base_size = 6) +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.3, 0.15),
        legend.direction = "horizontal", legend.background = element_rect(fill = alpha("white", 0)),
        legend.key.width = unit(0.5, "line"), legend.key.height = unit(0.5, "line"))

  
## Plot window results
## 
## Make this supplemental?
## 
historical_window_selection_out1_toplot <- historical_timeframe_selection_out1 %>%
  filter(time_frame_type == "window", model == "model3") %>%
  mutate(length = fct_inseq(as.factor(length))) %>%
  mutate_at(vars(contains("year")), ~ymd(paste0(., "0101")))


# Plot together
g_window_analysis <- historical_window_selection_out1_toplot %>%
  ggplot(aes(x = end_year, y = RMSE_scale, color = trait, lty = length)) +
  geom_line(lwd = 0.25) +
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Scaled RMSEP") +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y", name = "End year of window") +
  scale_color_manual(values = trait_colors, guide = FALSE) +
  scale_linetype_discrete(name = "Window length (yrs)", guide = guide_legend(title.position = "top")) +
  theme_genetics(base_size = 6) +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.3, 0.15),
        legend.direction = "horizontal", legend.key.width = unit(0.5, "line"), legend.key.height = unit(0.5, "line"))



# Combine plot
g_time_frame_combine <- plot_grid(
  g_timeframe_analysis, 
  g_window_analysis + theme(axis.title.y = element_blank(), axis.text.y = element_blank()), 
  labels = subfigure_labels[1:2], label_size = 8, label_x = c(0, -0.05), rel_widths = c(1, 0.95))

# Save
ggsave(filename = "figure3_partA_draft1.jpg", plot = g_time_frame_combine, path = fig_dir,
       height = 4, width = 10, units = "cm", dpi = dpi_use)




## Part C - accuracy within environments - LOLO

## Filter for relevant cases
accuracy_within_loc_toplot <- predictive_ability %>%
  ## Adjust model names
  mutate(model = str_replace_all(model, c("model2" = "model4", "model3" = "model5"))) %>%
  filter(type == "lolo") %>%
  filter(model %in% str_subset(models_use, "id", negate = TRUE), selection %in% c("none", "rfa_cv_adhoc")) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), "\n(n = ", nGroup+1, ")"),
         pop = paste0(f_validation_replace(pop), " (", f_pop_replace(pop), ")")) %>%
  unite(group, trait, model, pop, remove = FALSE)

accuracy_within_loc_toplot2 <- accuracy_within_loc_toplot %>%
  group_by(type, group, trait, trait1, model, pop, selection) %>%
  do({
    df <- .
    # Get ranges from boxplot
    bp <- boxplot(x = df$ability, plot = FALSE)
    # Return a tibble
    bind_cols(summarize(df, ability_mean = mean(ability)), 
              as_tibble(setNames(object = as.list(bp$stats), nm = c("lower", "q25", "median", "q75", "upper"))) )
  }) %>% ungroup()



# Plot point and line range
g_accuracy_within_loc <- accuracy_within_loc_toplot %>%
  ggplot(aes(x = trait1, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = ability), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = accuracy_within_loc_toplot2, aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = accuracy_within_loc_toplot2, aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = accuracy_within_loc_toplot2, aes(y = ability_mean), position = position_dodge(0.9), size = 0.8) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = expression('Within-location'~italic(r[MP])), breaks = pretty) +
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)) +
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title = NULL)) +
  facet_grid(. ~ pop, switch = "y") +
  theme_genetics(base_size = 6) +
  theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
        legend.text.align = 1, legend.position = c(0.89, 0.95), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line")) +
  theme(axis.title.x = element_blank(), strip.text.y = element_text(color = "white"))

## Save
ggsave(filename = "figure3_partc_draft.jpg", plot = g_accuracy_within_loc, path = fig_dir,
       width = 80, height = 40, dpi = dpi_use, units = "mm")



## Part D - plot predicted/observed phenotypes for each trait and by target population
## for model 3 with stepwise covariates - LOLO

lolo_pred_obs_df <- predictions_df %>%
  filter(type == "lolo") %>%
  filter(model == "model3_cov", selection == "rfa_cv_adhoc") %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  # mutate(trait = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"))
  mutate(trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))

# Split by trait and Plot
g_lolo_pred_obs_list <- lolo_pred_obs_df %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.2) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(pop ~ trait1, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_pop_replace)) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


# Combine the plots
g_lolo_pred_obs <- g_lolo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.y = element_blank(), axis.title.y = element_blank())) %>%
  map(~. + theme(axis.title.x = element_blank())) %>%
  plot_grid(plotlist = ., nrow = 1, align = "h", rel_widths = c(1, rep(0.90, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1)


# Save
ggsave(filename = "figure3_partD_draft.jpg", plot = g_lolo_pred_obs, path = fig_dir,
       height = 5, width = 16, units = "cm", dpi = dpi_use)



## Combine plots
top_plot <- plot_grid(g_time_frame_combine, g_accuracy_within_loc, nrow = 1, 
                      labels = c("", subfigure_labels[3]), label_size = 8, rel_widths = c(1, 0.8))

g_figure3 <- plot_grid(top_plot, g_lolo_pred_obs, ncol = 1, rel_heights = c(0.75, 1), 
                       labels = c("", subfigure_labels[4]), label_size = 8)

## Save
ggsave(filename = "figure3_lolo_predictions.jpg", plot = g_figure3, path = fig_dir,
       width = 16, height = 8.5, dpi = dpi_use, units = "cm")











# Figure 4. External environment prediction accuracy ---------------------------------------


## Part A - plot predicted/observed phenotypes for each trait and by target population
## for model 3 with stepwise covariates

external_pred_obs_df <- predictions_df %>%
  filter(str_detect(type, "ext")) %>%
  filter(model == "model3_cov", selection == "rfa_cv_adhoc") %>%
  # Add annotations - prediction accuracy across and within environments/locations
  left_join(., predictive_ability) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  # mutate(trait = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"))
  mutate(trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))

# Split by trait and Plot
g_external_pred_obs_list <- external_pred_obs_df %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(select(.x, trait, pop, type,  leave_one_group, ability)) %>% 
      mutate(leave_one_group1 =  ifelse(str_detect(leave_one_group, "[0-9]{2}"), leave_one_group, 
                                                   str_to_title(str_replace_all(leave_one_group, "_", " "))),
             ann1 = paste0("r[MP]==", formatC(signif(ability, 2), digits = 2, format = "f", flag = "0")),
             annotation = paste0(abbreviate(leave_one_group1, 5), "~", ann1),
             x = ifelse(str_detect(type, "env") & trait %in% c("PlantHeight", "TestWeight"), -Inf, Inf), y = -Inf)
    
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.2) +
      geom_text_repel(data = ann_df, aes(x = x, y = y, label = annotation, color = leave_one_group), inherit.aes = FALSE,
                      parse = TRUE, hjust = 1, size = 1, direction = "y", segment.alpha = 0, box.padding = 0.01) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(type ~ trait1, switch = "y", labeller = labeller(trait1 = label_parsed, type = f_type_replace)) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


## Combine and save
g_holdout_pred_obs <- g_external_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.y = element_blank(), axis.title.y = element_blank())) %>%
  map(~. + theme(axis.title.x = element_blank())) %>%
  plot_grid(plotlist = ., nrow = 1, align = "h", rel_widths = c(1, rep(0.80, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1)

# Save
ggsave(filename = "figure4_partA_draft1.jpg", plot = g_holdout_pred_obs, path = fig_dir,
       height = 5, width = 12, units = "cm", dpi = 2000)




## Supplementary figures and tables ##

## Number and overlap of covariates in the analyses ##

# Load the data
load(file = file.path(result_dir, "feature_selection_results.RData"))
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


## Concurrent

## Combine concurrent feature selection df
concurrent_features <- bind_rows(
  concurrent_fact_reg_feature_selection %>% select(source, trait, model, apriori, stepAIC_adhoc = adhoc) %>% 
    gather(feat_sel_type, features, apriori, stepAIC_adhoc),
  gather(concurrent_feature_selection, feat_subset, features, adhoc, adhoc_nosoil) %>%
    unite(feat_sel_type, feat_sel_type, feat_subset, sep = "_")
)

## combine model2 and model3 covariates
concurrent_features1 <- concurrent_features %>%
  filter(str_detect(feat_sel_type, "adhoc_nosoil", negate = TRUE)) %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model2, model3, union)) %>%
  ## Add all covariates
  bind_rows(.,
            tibble(source = names(ec_tomodel_centered), features = map(ec_tomodel_centered, ~names(.)[-1:-2])) %>% 
              mutate(features = map(features, ~c(., paste0("line_name:", .)))) %>% 
              crossing(., trait = traits, feat_sel_type = "all")
  ) %>%
  mutate(features = map(features, ~setdiff(., "line_name"))) %>%
  select(-contains("model")) %>%
  mutate(interaction_features = map(features, ~str_subset(., ":")),
         main_features = map2(features, interaction_features, setdiff))


## Calculate the number of covariates in each contingency
concurrent_features_count <- concurrent_features1 %>% 
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  filter(feature_class != "features")



## Compare trait, selections, etc. for covariates

# These will be polished plots for the manuscript
concurrent_features_count1 <- concurrent_features_count %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc")

concurrent_features2 <- concurrent_features1 %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc")


## Modify the feature count plot to remove nasapower
g_concurrent_features_count2 <- concurrent_features_count1 %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  filter(source == "daymet") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_count1$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_count1$feature_class == "main_features", 25, 10), size = 3) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "top", strip.placement = "outside")

ggsave(filename = "concurrent_features_count_daymet.jpg", plot = g_concurrent_features_count2,
       path = fig_dir, height = 8, width = 3.5, dpi = 1000)



## Overlap between non-all feature selection types
concurrent_features_overlap_featsel <- concurrent_features2 %>%
  select(-features) %>%
  filter(! feat_sel_type %in% c("all")) %>%
  full_join(., ., by = c("trait", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(concurrent_features1$feat_sel_type), m = 2))), 
                        ~paste0("feat_sel_type", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap"))



## Plot
g_concurrent_features_overlap_featsel <- concurrent_features_overlap_featsel %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("feat_sel")), f_ec_selection_replace) %>%
  unite(feat_sel_type_pair, feat_sel_type.x, feat_sel_type.y, sep = "-") %>%
  ggplot(aes(x = trait, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_overlap_featsel$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_overlap_featsel$feature_class == "main_features", 0, -2), size = 3) +
  facet_grid(~ feat_sel_type_pair, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "concurrent_features_overlap_featsel_daymet.jpg", plot = g_concurrent_features_overlap_featsel,
       path = fig_dir, height = 3, width = 4, dpi = 1000)


## Overlap between traits
concurrent_features_overlap_trait <- concurrent_features2 %>%
  select(-features) %>%
  full_join(., ., by = c("feat_sel_type", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(concurrent_features2$trait), m = 2))), 
                        ~paste0("trait", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap")) %>%
  # We don't need to see apriori or all
  filter(! feat_sel_type %in% c("all"))



## Plot
g_concurrent_features_overlap_trait <- concurrent_features_overlap_trait %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("trait")), ~abbreviate(str_add_space(.), 2)) %>%
  unite(trait_pair, trait.x, trait.y, sep = ":") %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_overlap_trait$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_overlap_trait$feature_class == "main_features", 7, 2), size = 2) +
  facet_grid(feat_sel_type ~ ., switch = "y", labeller = labeller(feat_sel_type = f_ec_selection_replace)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "concurrent_features_overlap_trait_daymet.jpg", plot = g_concurrent_features_overlap_trait,
       path = fig_dir, height = 3, width = 3.5, dpi = 1000)



## Count the number of times a particular covariate is select

concurrent_features3 <- concurrent_features2 %>%
  select(-features) %>%
  gather(feature_class, covariates, contains("features")) %>%
  unnest() %>%
  mutate(covariates = str_remove(covariates, "line_name:"))


# Counts by source and feat_sel_type
concurrent_indiv_feature_counts <- concurrent_features3 %>%
  # Remove all and apriori feat selections
  filter(! feat_sel_type %in% c("all", "apriori")) %>%
  group_by(covariates, source, feature_class, ) %>%
  summarize(n = n()) %>%
  mutate(nTotal = sum(n)) %>%
  ungroup() %>%
  mutate(covariates = fct_reorder(covariates, nTotal, .fun = max, .desc = TRUE))


## Plot
g_concurrent_indiv_feature_counts <- concurrent_indiv_feature_counts %>%
  ggplot(aes(x = covariates, y = n, fill = feature_class)) +
  geom_col() +
  scale_fill_discrete(labels = c("main_features" = "Environment covariate", "interaction_features" = "G x E covariate"), name = NULL) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(legend.position = c(0.75, 0.75), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Save
ggsave(filename = "concurrent_indiv_feature_counts_daymet.jpg", plot = g_concurrent_indiv_feature_counts, 
       path = fig_dir, width = 4, height = 5, dpi = 1000)


## Stitch these plots together - patchwork
left_plot <- plot_grid(g_concurrent_features_count2 + theme(legend.position = "none"), g_concurrent_features_overlap_trait, 
                       ncol = 1, rel_heights = c(1, 0.45), align = "hv", axis = "lr", labels = letters[1:2],
                       label_size = 10)  
right_plot <- plot_grid(g_concurrent_features_overlap_featsel, g_concurrent_indiv_feature_counts, ncol = 1,
                        rel_heights = c(0.6, 1), align = "hv", axis = "lr", labels = letters[3:4], label_size = 10)  

merged_plot <- (left_plot | right_plot) + plot_layout(widths = c(1, 1))

# Save
ggsave(filename = "concurrent_features_comparison_merged_daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = 1000)





# Supplementary fig XX - environmental relationships ----------------------

# Matrix of covariates
ec_tomodel_scaled_mat <- ec_tomodel_scaled$daymet %>%
  select(-source) %>%
  filter(environment %in% c(train_test_env, validation_env)) %>%
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()

## Calculate environmental relationship matrices based on these features
concurrent_features_env_relmat <- concurrent_features1 %>%
  filter(feat_sel_type == "rfa_cv_adhoc", source == "daymet") %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate(interaction_features_ec_mat = map(interaction_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE]),
         main_features_ec_mat = map(main_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE])) %>%
  mutate_at(vars(ends_with("ec_mat")), ~map(., ~Env_mat(x = .x, method = "Jarq")))

concurrent_features_env_relmat1 <- concurrent_features_env_relmat %>%
  select(trait, feat_sel_type, interaction = interaction_features_ec_mat, main = main_features_ec_mat) %>%
  gather(covariate_type, Emat, main, interaction) %>%
  ## Remove lines where matrices are all NA
  filter(map_lgl(Emat, ~!any(is.na(.))))


## For each trait and covariate type (interaction/main), plot a heatmap and
## a dendrogram
concurrent_features_heatmap_plots <- concurrent_features_env_relmat1 %>%
  group_by(trait, feat_sel_type, covariate_type) %>%
  do(plot = {
    row <- .
    
    # First perform clustering on the relationship matrix
    env_clust <- hclust(dist(row$Emat[[1]], method = "euclidean"), method = "average")
    # get the data
    clust_data <- dendro_data(model = env_clust)
    # Label data
    clust_lab_data <- mutate(clust_data$labels, group = ifelse(label %in% train_test_env, "Train/test", "Holdout"),
                             group = fct_relevel(group, "Holdout", after = Inf))

    # Factor order of environments
    env_order <- levels(clust_lab_data$label)
    
    # Create the plotting data.frame
    dat <- row$Emat[[1]] %>%
      as.data.frame(.) %>% 
      rownames_to_column(., "environment") %>%
      gather(environment2, relationship, -environment) %>%
      # Refactor the environments
      mutate_at(vars(contains("environment")), ~factor(., levels = env_order)) %>%
      mutate_at(vars(contains("environment")), list(group = ~ifelse(. %in% train_test_env, "training", "external"))) %>%
      mutate(environment = fct_rev(environment))
    
    # Plot
    g_heat <- dat %>%
      ggplot(aes(x = environment, y = environment2, fill = relationship)) +
      geom_tile() +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3]) +
      scale_y_discrete(position = "right") +
      labs(subtitle = paste(str_add_space(row$trait), str_to_title(row$covariate_type), sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = "none",
            axis.text.y = element_blank())
    
    ## Plot dendrogram
    g_dendro <- ggplot(segment(clust_data)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_text(data = clust_lab_data, aes(x = x, y = y, label = label, color = group), hjust = 1, 
                size = 2.5, nudge_y = -0.1) +
      coord_flip() + 
      scale_y_continuous(expand=c(0.2, 0)) + 
      scale_color_manual(name = NULL, values = c("black", "red")) +
      theme_genetics(base_size = 8) +
      theme_dendro() +
      theme(legend.position = c(0.8, 0.90))

    # Combine plots
    plot_grid(g_heat, g_dendro, align = "v", axis = "tb", rel_widths = c(1, 0.7))
    
  }) %>% ungroup()

## Combine all plots
g_heat_dendro_all <- plot_grid(plotlist = concurrent_features_heatmap_plots$plot, 
                               nrow = ceiling(nrow(concurrent_features_heatmap_plots)/2))
# Save
ggsave(filename = "figure_sXX_environmental_covariate_heatmap.jpg", plot = g_heat_dendro_all, 
       path = fig_dir, width = 12, height = 15, dpi = 1000)





## Historical

## Combine historical feature selection df
historical_features <- bind_rows(
  historical_fact_reg_feature_selection %>% select(source, trait, model, apriori, stepAIC_adhoc = adhoc) %>% 
    gather(feat_sel_type, features, apriori, stepAIC_adhoc),
  historical_feature_selection %>% 
    mutate(model = ifelse(model == "model2", "model4", "model5"), source = "daymet") %>%
    gather(feat_subset, features, adhoc) %>%
    unite(feat_sel_type, feat_sel_type, feat_subset, sep = "_")
)

## combine model2 and model3 covariates
historical_features1 <- historical_features %>%
  filter(str_detect(feat_sel_type, "adhoc_nosoil", negate = TRUE)) %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model4, model5, union)) %>%
  ## Add all covariates
  bind_rows(.,
            tibble(source = names(ec_tomodel_centered), 
                   features = list(names(historical_ec_tomodel_timeframe_centered[[1]])[-1:-3])) %>% 
              mutate(features = map(features, ~c(., paste0("line_name:", .)))) %>% 
              crossing(., trait = traits, feat_sel_type = "all")
  ) %>%
  mutate(features = map(features, ~setdiff(., "line_name"))) %>%
  select(-contains("model")) %>%
  mutate(interaction_features = map(features, ~str_subset(., ":")),
         main_features = map2(features, interaction_features, setdiff)) %>%
  ## Use the time_frame selected by the stepwise search
  select(-time_frame) %>%
  left_join(., distinct(historical_feature_selection, trait, time_frame))


## Calculate the number of covariates in each contingency
historical_features_count <- historical_features1 %>% 
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  filter(feature_class != "features")



## Compare trait, selections, etc. for covariates

# These will be polished plots for the manuscript
historical_features_count1 <- historical_features_count %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc")

historical_features2 <- historical_features1 %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc")


## Modify the feature count plot to remove nasapower
g_historical_features_count1 <- historical_features_count1 %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  filter(source == "daymet") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_count1$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_count1$feature_class == "main_features", 25, 10), size = 3) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "top", strip.placement = "outside")

ggsave(filename = "historical_features_count_daymet.jpg", plot = g_historical_features_count1,
       path = fig_dir, height = 8, width = 3.5, dpi = 1000)



## Overlap between non-all feature selection types
historical_features_overlap_featsel <- historical_features2 %>%
  select(-features) %>%
  filter(! feat_sel_type %in% c("all")) %>%
  full_join(., ., by = c("trait", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(historical_features1$feat_sel_type), m = 2))), 
                        ~paste0("feat_sel_type", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap"))



## Plot
g_historical_features_overlap_featsel <- historical_features_overlap_featsel %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("feat_sel")), f_ec_selection_replace) %>%
  unite(feat_sel_type_pair, feat_sel_type.x, feat_sel_type.y, sep = "-") %>%
  ggplot(aes(x = trait, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_overlap_featsel$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_overlap_featsel$feature_class == "main_features", 0, -2), size = 3) +
  facet_grid(~ feat_sel_type_pair, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "historical_features_overlap_featsel_daymet.jpg", plot = g_historical_features_overlap_featsel,
       path = fig_dir, height = 3, width = 4, dpi = 1000)


## Overlap between traits
historical_features_overlap_trait <- historical_features2 %>%
  select(-features) %>%
  full_join(., ., by = c("feat_sel_type", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(historical_features2$trait), m = 2))), 
                        ~paste0("trait", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap")) %>%
  # We don't need to see apriori or all
  filter(! feat_sel_type %in% c("all"))



## Plot
g_historical_features_overlap_trait <- historical_features_overlap_trait %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("trait")), ~abbreviate(str_add_space(.), 2)) %>%
  unite(trait_pair, trait.x, trait.y, sep = ":") %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_overlap_trait$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_overlap_trait$feature_class == "main_features", 7, 2), size = 2) +
  facet_grid(feat_sel_type ~ ., switch = "y", labeller = labeller(feat_sel_type = f_ec_selection_replace)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "historical_features_overlap_trait_daymet.jpg", plot = g_historical_features_overlap_trait,
       path = fig_dir, height = 3, width = 3.5, dpi = 1000)



## Count the number of times a particular covariate is select

historical_features3 <- historical_features2 %>%
  select(-features) %>%
  gather(feature_class, covariates, contains("features")) %>%
  unnest() %>%
  mutate(covariates = str_remove(covariates, "line_name:"))


# Counts by source and feat_sel_type
historical_indiv_feature_counts <- historical_features3 %>%
  # Remove all and apriori feat selections
  filter(! feat_sel_type %in% c("all", "apriori")) %>%
  group_by(covariates, source, feature_class, ) %>%
  summarize(n = n()) %>%
  mutate(nTotal = sum(n)) %>%
  ungroup() %>%
  mutate(covariates = fct_reorder(covariates, nTotal, .fun = max, .desc = TRUE))


## Plot
g_historical_indiv_feature_counts <- historical_indiv_feature_counts %>%
  ggplot(aes(x = covariates, y = n, fill = feature_class)) +
  geom_col() +
  scale_fill_discrete(labels = c("main_features" = "Environment covariate", "interaction_features" = "G x E covariate"), name = NULL) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(legend.position = c(0.75, 0.75), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Save
ggsave(filename = "historical_indiv_feature_counts_daymet.jpg", plot = g_historical_indiv_feature_counts, 
       path = fig_dir, width = 4, height = 5, dpi = 1000)


## Stitch these plots together - patchwork
left_plot <- plot_grid(g_historical_features_count1 + theme(legend.position = "none"), g_historical_features_overlap_trait, 
                       ncol = 1, rel_heights = c(1, 0.45), align = "hv", axis = "lr", labels = letters[1:2],
                       label_size = 10)  
right_plot <- plot_grid(g_historical_features_overlap_featsel, g_historical_indiv_feature_counts, ncol = 1,
                        rel_heights = c(0.6, 1), align = "hv", axis = "lr", labels = letters[3:4], label_size = 10)  

merged_plot <- (left_plot | right_plot) + plot_layout(widths = c(1, 1))

# Save
ggsave(filename = "historical_features_comparison_merged_daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = 1000)




# Supplementary fig XX - location relationships ----------------------

# list of covariate matrices
historical_ec_tomodel_list <- historical_ec_tomodel_timeframe_centered %>%
  subset(., str_detect(names(.), "daymet")) %>%
  bind_rows() %>%
  select(-source) %>%
  filter(time_frame %in% unique(historical_feature_selection$time_frame)) %>%
  split(.$time_frame) %>%
  map(~select(., -time_frame) %>% as.data.frame() %>% column_to_rownames("location") %>% as.matrix())

## Calculate environmental relationship matrices based on these features
historical_features_env_relmat <- historical_features1 %>%
  filter(feat_sel_type %in% c("all", "rfa_cv_adhoc"), source == "daymet") %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate(interaction_features_ec_mat = map2(interaction_features, time_frame, ~historical_ec_tomodel_list[[.y]][, .x, drop = FALSE]),
         main_features_ec_mat = map2(main_features, time_frame, ~historical_ec_tomodel_list[[.y]][, .x, drop = FALSE])) %>%
  mutate_at(vars(ends_with("ec_mat")), ~map(., ~Env_mat(x = .x, method = "Jarq")))

historical_features_env_relmat1 <- historical_features_env_relmat %>%
  select(trait, feat_sel_type, interaction = interaction_features_ec_mat, main = main_features_ec_mat) %>%
  gather(covariate_type, Emat, main, interaction) %>%
  ## Remove lines where matrices are all NA
  filter(map_lgl(Emat, ~!any(is.na(.))))


## For each trait and covariate type (interaction/main), plot a heatmap and
## a dendrogram
historical_features_heatmap_plots <- historical_features_env_relmat1 %>%
  # No need for all and interaction covariates
  filter(!(feat_sel_type == "all" & covariate_type == "interaction")) %>%
  group_by(trait, feat_sel_type, covariate_type) %>%
  do(plot = {
    row <- .
    
    # First perform clustering on the relationship matrix
    env_clust <- hclust(dist(row$Emat[[1]], method = "euclidean"), method = "average")
    # get the data
    clust_data <- dendro_data(model = env_clust)
    # Label data
    clust_lab_data <- mutate(clust_data$labels, group = ifelse(label %in% train_test_loc, "Train/test", "Holdout"),
                             group = fct_relevel(group, "Holdout", after = Inf),
                             label1 = str_to_title(str_replace_all(label, "_", " ")),
                             label1 = ifelse(str_detect(label1, " ") & nchar(label1) > 10, str_replace(label1, " ", "\n"), label1))
    # Segment data
    clust_seg_data <- segment(clust_data)
    # Factor order of environments
    env_order <- levels(clust_lab_data$label)
    
    # Create the plotting data.frame
    dat <- row$Emat[[1]] %>%
      as.data.frame(.) %>% 
      rownames_to_column(., "location") %>%
      gather(location2, relationship, -location) %>%
      # Refactor the environments
      mutate_at(vars(contains("location")), ~factor(., levels = env_order)) %>%
      mutate(location = fct_rev(location))
    
    # Plot
    g_heat <- dat %>%
      ggplot(aes(x = location, y = location2, fill = relationship)) +
      geom_tile() +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3]) +
      scale_y_discrete(position = "right") +
      scale_x_discrete(labels = function(x) str_to_title(str_replace_all(x, "_", " "))) +
      labs(subtitle = paste(str_add_space(row$trait), str_to_title(row$covariate_type), sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = "none",
            axis.text.y = element_blank())
    
    x_nudge <- 0.01 * max(clust_seg_data$y)
    
    ## Plot dendrogram
    g_dendro <- ggplot(data = clust_seg_data) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_text(data = clust_lab_data, aes(x = x, y = y, label = label1, color = group), hjust = 1, 
                size = 2.5, nudge_y = -x_nudge) +
      coord_flip() + 
      scale_y_continuous(expand = c(0.5, 0)) + 
      scale_color_manual(name = NULL, values = c("black", "red")) +
      theme_genetics(base_size = 8) +
      theme_dendro() +
      theme(legend.position = c(0.8, 0.90))
    
    # Combine plots
    plot_grid(g_heat, g_dendro, align = "h", axis = "tblr", rel_widths = c(1, 0.7), nrow = 1)
    
  }) %>% ungroup()

## Combine all plots
g_heat_dendro_all <- historical_features_heatmap_plots %>%
  filter(feat_sel_type == "all") %>%
  pull(plot) %>%
  plot_grid(plotlist = ., nrow = ceiling(length(.)/2))

# Save
ggsave(filename = "figure_sXX_location_covariate_heatmap_all.jpg", plot = g_heat_dendro_all, 
       path = fig_dir, width = 14, height = 14, dpi = 1000)

## Combine stepwise covariate plots
g_heat_dendro_stepwise <- historical_features_heatmap_plots %>%
  filter(feat_sel_type == "rfa_cv_adhoc") %>%
  pull(plot) %>%
  plot_grid(plotlist = ., nrow = ceiling(length(.)/2))

# Save
ggsave(filename = "figure_sXX_location_covariate_heatmap_stepwise.jpg", plot = g_heat_dendro_stepwise, 
       path = fig_dir, width = 14, height = 14, dpi = 1000)











# Supplemental Figure XX - LOEO predictive ability across environments ----------

## Plot predicted vs observed phenotypic values for LOEO for all models, traits,
## EC selections, etc.

# First create a data.frame from which to plot all scenarios

loo_pred_obs_df <- predictions_df %>%
  # Filter for relevant models
  filter(model %in% models_use) %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # Add environment/location specific predictive abilities
  left_join(., select(predictive_ability, trait:ability)) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  mutate(trait1 = paste0("'", str_add_space(trait), "'~(", trait_units1[trait], ")")) %>%
  # Group by trait and add a dummy faceting variable
  unite(dummy, model, selection, sep = "_", remove = FALSE) %>%
  mutate(dummy = as.factor(dummy)) %>%
  group_by(trait, type) %>%
  mutate(model = f_model_replace(model),
         model = ifelse(dummy == first(levels(dummy)), paste0("Model: ", model), model),
         selection = f_ec_selection_replace(selection),
         selection = ifelse(dummy == first(levels(dummy)), paste0("Covariate Set: ", selection), selection)) %>%
  ungroup() %>%
  mutate_at(vars(model, selection), fct_inorder)



## Plot all
g_loeo_all_list <- loo_pred_obs_df %>%
  filter(type == "loeo") %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, model, selection, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(data = .x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.2) +
      geom_point(size = 0.2) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(pop ~ model + selection, switch = "y", scales = "free_x", labeller = labeller(pop = f_pop_replace)) +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_loeo_all <- plot_grid(plotlist = g_loeo_all_list, ncol = 1, align = "v", axis = "lr",
                        labels = subfigure_labels[seq_along(g_loeo_all_list)], label_size = 8)


# Save
ggsave(filename = "figure_sXX_loeo_across_site_prediction_accuracy.jpg", plot = g_loeo_all, path = fig_dir,
       height = 20, width = 22, units = "cm", dpi = 500)



# Supplemental Figure XX - LOEO prediction accuracy and bias within environments ----------

## Filter for relevant cases
accuracy_bias_within_env_toplot <- within_environment_prediction_accuracy %>%
  filter(type == "loeo") %>%
  filter(model %in% models_use) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         selection = fct_inorder(selection),
         model = str_replace(model, "_id", "_cov"),
         bias = bias * 100) %>%
  unite(group, trait, model, pop, remove = FALSE) %>%
  select(group:selection, trait1, accuracy, bias) %>%
  gather(statistic, value, accuracy, bias)

accuracy_bias_within_env_toplot2 <- accuracy_bias_within_env_toplot %>%
  group_by(type, group, trait, trait1, model, pop, selection, statistic) %>%
  do({
    df <- .
    # Get ranges from boxplot
    bp <- boxplot(x = df$value, plot = FALSE)
    # Return a tibble
    bind_cols(summarize(df, value_mean = mean(value)), 
              as_tibble(setNames(object = as.list(bp$stats), nm = c("lower", "q25", "median", "q75", "upper"))) )
  }) %>% ungroup()


## Common plot modifier
g_mod <- list(
  scale_x_discrete(name = "Covariate set", labels = f_ec_selection_replace),
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)),
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title = NULL)),
  facet_grid(type ~ pop + trait1, switch = "y", labeller = labeller(pop = f_pop_replace)),
  theme_genetics(base_size = 6),
  theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
        legend.text.align = 1, legend.position = c(0.89, 0.95), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line"),
        strip.text.y = element_text(color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)


# Plot point and line range for accuracy
g_accuracy_within_env <- accuracy_bias_within_env_toplot %>%
  filter(statistic == "accuracy") %>%
  ggplot(aes(x = selection, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = filter(accuracy_bias_within_env_toplot2, statistic == "accuracy"), 
                 aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = filter(accuracy_bias_within_env_toplot2, statistic == "accuracy"), 
                 aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = filter(accuracy_bias_within_env_toplot2, statistic == "accuracy"), 
             aes(y = value_mean), position = position_dodge(0.9), size = 0.8) +
  scale_y_continuous(name = expression('Within-environment'~italic(r[MG])), breaks = pretty) +
  g_mod
  

# Plot point and line range for bias
g_bias_within_env <- accuracy_bias_within_env_toplot %>%
  filter(statistic == "bias") %>%
  ggplot(aes(x = selection, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = filter(accuracy_bias_within_env_toplot2, statistic == "bias"), 
                 aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = filter(accuracy_bias_within_env_toplot2, statistic == "bias"), 
                 aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = filter(accuracy_bias_within_env_toplot2, statistic == "bias"), 
             aes(y = value_mean), position = position_dodge(0.9), size = 0.8) +
  scale_y_continuous(name = "Within-environment bias (%)", breaks = pretty) +
  g_mod

# Combine plots
g_accuracy_bias_within_env <- plot_grid(g_accuracy_within_env, g_bias_within_env, align = "hv",
                                        labels = subfigure_labels[1:2], label_size = 8)


## Save
ggsave(filename = "figure_sXX_loeo_accuracy_bias_within_environments.jpg", plot = g_accuracy_bias_within_env, 
       path = fig_dir, width = 12, height = 4, dpi = 500)




# Supplemental Figure XX - LOLO predictive ability across locations ----------

## Plot predicted vs observed phenotypic values for LOLO for all models, traits,
## EC selections, etc.

## Plot all
g_lolo_all_list <- loo_pred_obs_df %>%
  filter(type == "lolo") %>%
  filter(str_detect(selection, "Concurrent|AIC", negate = TRUE)) %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, model, selection, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(data = .x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.2) +
      geom_point(size = 0.2) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(pop ~ model + selection, switch = "y", scales = "free_x", labeller = labeller(pop = f_pop_replace)) +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_lolo_all <- plot_grid(plotlist = g_lolo_all_list, ncol = 1, align = "v", axis = "lr",
                        labels = subfigure_labels[seq_along(g_lolo_all_list)], label_size = 8)

# Save
ggsave(filename = "figure_sXX_lolo_across_site_prediction_accuracy.jpg", plot = g_lolo_all, path = fig_dir,
       height = 20, width = 22, units = "cm", dpi = 500)



# Supplemental Figure XX - LOLO prediction accuracy and bias within locations ----------

## Filter for relevant cases
accuracy_bias_within_loc_toplot <- predictive_ability %>%
  filter(type == "lolo", model %in% models_use, str_detect(selection, "concurrent|AIC", negate = TRUE)) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         selection = fct_inorder(selection),
         model = str_replace(model, "_id", "_cov"),
         bias = bias * 100) %>%
  unite(group, trait, model, pop, remove = FALSE) %>%
  select(group:selection, trait1, accuracy = ability, bias) %>%
  gather(statistic, value, accuracy, bias)

accuracy_bias_within_loc_toplot2 <- accuracy_bias_within_loc_toplot %>%
  group_by(type, group, trait, trait1, model, pop, selection, statistic) %>%
  do({
    df <- .
    # Get ranges from boxplot
    bp <- boxplot(x = df$value, plot = FALSE)
    # Return a tibble
    bind_cols(summarize(df, value_mean = mean(value)), 
              as_tibble(setNames(object = as.list(bp$stats), nm = c("lower", "q25", "median", "q75", "upper"))) )
  }) %>% ungroup()


# Plot point and line range for accuracy
g_accuracy_within_loc <- accuracy_bias_within_loc_toplot %>%
  filter(statistic == "accuracy") %>%
  ggplot(aes(x = selection, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = filter(accuracy_bias_within_loc_toplot2, statistic == "accuracy"), 
                 aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = filter(accuracy_bias_within_loc_toplot2, statistic == "accuracy"), 
                 aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = filter(accuracy_bias_within_loc_toplot2, statistic == "accuracy"), 
             aes(y = value_mean), position = position_dodge(0.9), size = 0.8) +
  scale_y_continuous(name = expression('Within-location'~italic(r[MP])), breaks = pretty) +
  g_mod


# Plot point and line range for bias
g_bias_within_loc <- accuracy_bias_within_loc_toplot %>%
  filter(statistic == "bias") %>%
  ggplot(aes(x = selection, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = filter(accuracy_bias_within_loc_toplot2, statistic == "bias"), 
                 aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = filter(accuracy_bias_within_loc_toplot2, statistic == "bias"), 
                 aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = filter(accuracy_bias_within_loc_toplot2, statistic == "bias"), 
             aes(y = value_mean), position = position_dodge(0.9), size = 0.8) +
  scale_y_continuous(name = "Within-location bias (%)", breaks = pretty) +
  g_mod

# Combine plots
g_accuracy_bias_within_loc <- plot_grid(g_accuracy_within_loc, g_bias_within_loc, align = "hv",
                                        labels = subfigure_labels[1:2], label_size = 8)


## Save
ggsave(filename = "figure_sXX_lolo_accuracy_bias_within_locations.jpg", plot = g_accuracy_bias_within_loc, 
       path = fig_dir, width = 12, height = 4, dpi = 500)



# Supplemental Figure XX - external predictive ability across environments ----------

## Plot all
g_ext_env_all_list <- loo_pred_obs_df %>%
  filter(type == "env_external") %>%
  filter(str_detect(selection, "Concurrent|AIC", negate = TRUE)) %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(select(.x, trait, pop, type,  leave_one_group, ability, model, selection)) %>% 
      mutate(leave_one_group1 =  ifelse(str_detect(leave_one_group, "[0-9]{2}"), leave_one_group, 
                                        str_to_title(str_replace_all(leave_one_group, "_", " "))),
             ann1 = paste0("r[MP]==", formatC(signif(ability, 2), digits = 2, format = "f", flag = "0")),
             annotation = paste0(abbreviate(leave_one_group1, 5), "~", ann1)) %>%
      mutate(x = case_when(trait == "TestWeight" & selection == "StepwiseCV" ~ -Inf,
                           TRUE ~ Inf),
             y = case_when(trait == "GrainYield" & str_detect(selection, "None", negate = TRUE) ~ Inf,
                           TRUE ~ -Inf))
    
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.2) +
      geom_label_repel(data = ann_df, aes(x = x, y = y, label = annotation, color = leave_one_group), inherit.aes = FALSE,
                       parse = TRUE, hjust = 1, size = 1, direction = "y", segment.alpha = 0, box.padding = 0.01, 
                       label.size = NA, label.padding = 0.01, fill = alpha("white", 0.5)) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(. ~ model + selection, switch = "y", scales = "free_x") +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_ext_env_all <- plot_grid(plotlist = g_ext_env_all_list, ncol = 1, align = "v", axis = "lr",
                           labels = subfigure_labels[seq_along(g_lolo_all_list)], label_size = 8)

# Save
ggsave(filename = "figure_sXX_external_environment_across_site_prediction_accuracy.jpg", plot = g_ext_env_all, path = fig_dir,
       height = 16, width = 22, units = "cm", dpi = 500)





# Supplemental Figure XX - external predictive ability across locations ----------

## Plot all
g_ext_loc_all_list <- loo_pred_obs_df %>%
  filter(type == "loc_external") %>%
  filter(str_detect(selection, "Concurrent|AIC", negate = TRUE)) %>%
  split(.$trait) %>%
  map(~{
    
    # Extract annotations
    ann_df <- distinct(select(.x, trait, pop, type,  leave_one_group, ability, model, selection)) %>% 
      mutate(leave_one_group1 =  ifelse(str_detect(leave_one_group, "[0-9]{2}"), leave_one_group, 
                                        str_to_title(str_replace_all(leave_one_group, "_", " "))),
             ann1 = paste0("r[MP]==", formatC(signif(ability, 2), digits = 2, format = "f", flag = "0")),
             annotation = paste0(abbreviate(leave_one_group1, 5), "~", ann1)) %>%
      mutate(x = case_when(trait == "GrainYield" & selection %in% c("All", "Literature") ~ -Inf,
                           TRUE ~ Inf),
             y = case_when(trait == "GrainYield" & selection == "All" ~ Inf,
                           TRUE ~ -Inf))
    
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.2) +
      geom_label_repel(data = ann_df, aes(x = x, y = y, label = annotation, color = leave_one_group), inherit.aes = FALSE,
                       parse = TRUE, hjust = 1, size = 1, direction = "y", segment.alpha = 0, box.padding = 0.01, 
                       label.size = NA, label.padding = 0.01, fill = alpha("white", 0.5)) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_paletteer_d(package = "ggsci", palette = "default_igv", guide = FALSE) +
      facet_grid(. ~ model + selection, switch = "y", scales = "free_x") +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_ext_loc_all <- plot_grid(plotlist = g_ext_loc_all_list, ncol = 1, align = "v", axis = "lr",
                           labels = subfigure_labels[seq_along(g_lolo_all_list)], label_size = 8)

# Save
ggsave(filename = "figure_sXX_external_location_across_site_prediction_accuracy.jpg", plot = g_ext_loc_all, path = fig_dir,
       height = 16, width = 22, units = "cm", dpi = 500)







# Supplemental Table XX - full model phenotypic variance analysis ---------

# Load the full model analysis results
load(file.path(result_dir, "full_model_variance_analysis.RData"))

# Collapse if a list
if (!is.data.frame(environment_pheno_variance_analysis)) {
  environment_pheno_variance_analysis <- environment_pheno_variance_analysis %>% 
    subset(., sapply(., is.data.frame)) %>%
    bind_rows()
}

# Unnest results
environment_pheno_variance_analysis1 <- environment_pheno_variance_analysis %>%
  # Make a note as to whether any GxE covariates were included
  mutate(any_gxe_cov = map_lgl(features, ~any(str_detect(., "line_name:"))) | feat_sel_type == "all") %>%
  unnest(results)

## Calculate the total variance and calculate the proportion explained by each source
environment_pheno_variance_analysis2 <- environment_pheno_variance_analysis1 %>%
  group_by(trait, population, feat_sel_type) %>%
  mutate(total_variance = sum(variance), prop_total_variance = variance / total_variance) %>%
  rename(prop_source_variance = variance_prop) %>%
  ungroup()

## Clean up for a table
environment_pheno_variance_analysis_table <- environment_pheno_variance_analysis2 %>%
  mutate(feat_sel_type = ifelse(feat_sel_type == "rfa_cv", "rfa_cv_adhoc", feat_sel_type)) %>%
  select(trait, population, covariate_set = feat_sel_type, any_gxe_cov, source, term, prop_source_variance, prop_total_variance) %>%
  mutate(trait = str_add_space(trait), population = f_pop_replace(population), covariate_set = f_ec_selection_replace(covariate_set),
         source = str_replace_all(source, c("line_name" = "G", "environment" = "E", "gxe" = "G x E", "units" = "Residuals")),
         term = str_replace_all(term, c("line_name" = "G", "environment" = "E", "gxe" = "G x E", "units" = "Residuals")),
         term = str_replace_all(term, c("G_cov" = "Markers", "G x E_cov" = "Markers x Covariates",  "E_cov" = "Covariates"))) %>%
  mutate_at(vars(contains("prop")), ~formatC(x = signif(., 2), digits = 2, width = 2, format = "f")) %>%
  arrange(trait, population, covariate_set) %>%
  rename_all(~str_to_title(str_replace_all(., "_", " ")))


# Calculate difference in variance explained
environment_pheno_variance_analysis_diff <- environment_pheno_variance_analysis_table %>%
  rename_all(make.names) %>%
  select(-Any.Gxe.Cov, -Prop.Total.Variance) %>%
  filter(Population == "FP", str_detect(Term, "Covariates")) %>%
  {full_join(x = filter(., Covariate.Set == "StepwiseCV"), y = filter(., Covariate.Set != "StepwiseCV"),
             by = c("Trait", "Population", "Source", "Term"))} %>%
  mutate_at(vars(contains("Prop")), parse_number) %>%
  mutate(prop_diff = Prop.Source.Variance.x - Prop.Source.Variance.y)


# Output table
write_csv(x = environment_pheno_variance_analysis_table, path = file.path(fig_dir, "table_sXX_full_model_phenotypic_variance_analysis.csv"))


## Subset for inclusion in main text
environment_pheno_variance_analysis_diff %>% 
  filter(Population == "FP", Term == "Markers x Covariates") %>% 
  as.data.frame() %>%
  arrange(Covariate.Set.y, desc(prop_diff))








# Supplemental Table XX - full model phenotypic variance analysis for locations ---------

# Collapse if a list
if (!is.data.frame(location_pheno_variance_analysis)) {
  location_pheno_variance_analysis <- location_pheno_variance_analysis %>% 
    subset(., sapply(., is.data.frame)) %>%
    bind_rows()
}

# Unnest results
location_pheno_variance_analysis1 <- location_pheno_variance_analysis %>%
  # Make a note as to whether any GxE covariates were included
  mutate(any_gxe_cov = map_lgl(features, ~any(str_detect(., "line_name:"))) | feat_sel_type == "all") %>%
  unnest(results)

## Calculate the total variance and calculate the proportion explained by each source
location_pheno_variance_analysis2 <- location_pheno_variance_analysis1 %>%
  group_by(trait, population, feat_sel_type) %>%
  mutate(total_variance = sum(variance), prop_total_variance = variance / total_variance) %>%
  rename(prop_source_variance = variance_prop) %>%
  ungroup()

## Clean up for a table
location_pheno_variance_analysis_table <- location_pheno_variance_analysis2 %>%
  mutate(feat_sel_type = ifelse(feat_sel_type == "rfa_cv", "rfa_cv_adhoc", feat_sel_type)) %>%
  select(trait, population, covariate_set = feat_sel_type, any_gxe_cov, source, term, prop_source_variance, prop_total_variance) %>%
  mutate(trait = str_add_space(trait), population = f_pop_replace(population), covariate_set = f_ec_selection_replace(covariate_set),
         source = str_replace_all(source, c("line_name" = "G", "location" = "L", "gxl" = "G x L", "units" = "Residuals")),
         term = str_replace_all(term, c("line_name" = "G", "location" = "L", "gxl" = "G x L", "units" = "Residuals")),
         term = str_replace_all(term, c("G_cov" = "Markers", "G x L_cov" = "Markers x Covariates",  "L_cov" = "Covariates"))) %>%
  mutate_at(vars(contains("prop")), ~formatC(x = signif(., 2), digits = 2, width = 2, format = "f")) %>%
  arrange(trait, population, covariate_set) %>%
  rename_all(~str_to_title(str_replace_all(., "_", " ")))


# Calculate difference in variance explained
location_pheno_variance_analysis_table %>%
  rename_all(make.names) %>%
  select(-Any.Gxe.Cov, -Prop.Total.Variance) %>%
  filter(Population == "FP", str_detect(Term, "Covariates")) %>%
  {full_join(x = filter(., Covariate.Set == "Stepwise"), y = filter(., Covariate.Set != "Stepwise"),
             by = c("Trait", "Population", "Source", "Term"))} %>%
  mutate_at(vars(contains("Prop")), parse_number) %>%
  mutate(prop_diff = Prop.Source.Variance.x - Prop.Source.Variance.y)


# Output table
write_csv(x = location_pheno_variance_analysis_table, 
          path = file.path(fig_dir, "full_model_phenotypic_variance_analysis_location.csv"))


























