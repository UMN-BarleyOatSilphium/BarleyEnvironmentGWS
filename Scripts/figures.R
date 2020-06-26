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
  mutate(set = ifelse(environment %in% train_test_env, "Train/test locations", "Validation locations")) %>%
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

# Save the map
ggsave(filename = "map_of_experiments.jpg", plot = g_map_alt, path = fig_dir,
       height = 3, width = 4.5, dpi = 500)

# ## Add an adjacent table
# g_map_alt1 <- g_map_alt +
#   gridExtra::tableGrob(d = trait_experiment_summary, rows = NULL, 
#                        theme = gridExtra::ttheme_minimal()) +
#   plot_layout(widths = c(0.5, 0.25))
# 
# plot_grid(g_map_alt, gridExtra::tableGrob(d = trait_experiment_summary, rows = NULL, 
#                                           theme = gridExtra::ttheme_minimal()),
#           rel_widths = c(1, 0.25)





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
cgm_locations_color <- setNames(paletteer_d("ggsci", palette = "nrc_npg", n = 2), cgm_locations_plot)

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
                                  label = letters[2:3], dap = min(cgm_toplot_growth_stage_weather_historical$dap))


## Plot of both concurrent and historical mean (plus individual year) data
g_weather_concurrent_historical <- cgm_toplot_growth_stage_weather_historical %>%
  filter(variable %in% covariates_plot) %>%
  mutate(variable = f_covariate_replace(variable)) %>%
  ggplot(aes(x = dap, y = value, color = location, group = group)) +
  geom_line(alpha = 0.2, lwd = 0.15) +
  # Add average line
  geom_line(data = cgm_toplot_growth_stage_weather_summary, lwd = 0.3) +
  # Add text for plot labels
  # geom_text(data = g_weather_subplot_label, aes(x = dap, y = Inf, vjust = 1, label = label), 
  #           fontface = 2, size = 3, inherit.aes = FALSE) +
  facet_grid(variable ~ timeframe, scales = "free", space = "free_x", switch = "y", 
             labeller = labeller(variable = label_parsed, timeframe = f_timeframe_replace)) +
  scale_y_continuous(name = NULL, breaks = pretty) +
  scale_x_continuous(name = "Days after planting date", breaks = pretty) +
  scale_color_manual(values = cgm_locations_color) +
  theme_genetics(base_size = 6) +
  theme(strip.placement = "outside", strip.background = element_blank(), panel.spacing.y = unit(0.5, "line"))





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
         # 5 years of data
         between(year, cgm_year_plot - 5, cgm_year_plot - 1)) %>%
  ggplot(aes(x = dap, color = stage)) +
  # Individual-year stages
  geom_segment(aes(x = min, xend = max, y = y, yend = y), lwd = 4, alpha = 0.01) +
  # Average stages
  geom_segment(data = cgm_toplot_growth_stage_times_summary, 
               aes(x = min, xend = max, y = y, yend = y), lwd = 2) +
  # Lines for the entire growing season per location
  geom_segment(data = cgm_toplot_growing_season, aes(x = min, xend = max, y = y, yend = y),
               color = cgm_locations_color[cgm_toplot_growing_season$location], lwd = 0.5, inherit.aes = FALSE) +
  # Annotation for growth stages
  geom_text(data = subset(cgm_toplot_growth_stage_times_summary, y == 1),
            aes(x = (min + max) / 2, y = y - 1, label = f_growth_stage_replace(as.character(stage))), size = 1.5) +
  facet_grid(type ~ timeframe, labeller = labeller(type = str_add_space), switch = "y", scales = "free_x", space = "free_x") +
  scale_color_manual(values = growth_stage_color, guide = FALSE) +
  scale_y_continuous(limits = c(-0.5, 2.5)) +
  scale_x_continuous(breaks = pretty, name = "Days after planting date") +
  theme_genetics(base_size = 6) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank(),
        axis.title.y = element_blank(), axis.line.y = element_blank())


## Combine weather and growth stage plot
g_weather_stages <- plot_grid(
  g_weather_concurrent_historical + theme(legend.position = "none", axis.title.x = element_blank(),
                                          axis.text.x = element_blank()), 
  g_growth_stages + theme(strip.text.x = element_blank()),
  ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.25), labels = letters[2:3], label_size = 8)



## Highlight locations with color in the map
map_location_colors <- setdiff(g_map_alt$layers[[4]]$data$location, names(cgm_locations_color)) %>%
  setNames(rep("black", length(.)), .) %>%
  c(., cgm_locations_color)
# New map plot
g_map_alt1 <- g_map_alt
g_map_alt1$layers <- g_map_alt1$layers[-5]
g_map_alt1 <- g_map_alt1 +
  geom_text_repel(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = location, color = location),
                  size = 2, hjust = 0.5, nudge_x = ifelse(use_loc_info_toplot1$location == "Bozeman", 3, -1), segment.size = 0.2, 
                  point.padding = unit(2, "pt"), min.segment.length = 1) +
  scale_color_manual(values = map_location_colors, guide = FALSE) +
  theme(legend.position = c(0.88, 0.20))
  

g_fig1 <- plot_grid(g_map_alt1, g_weather_stages, ncol = 1, labels = letters[1], label_size = 8,
                    rel_heights = c(0.4, 1))

# Save
ggsave(filename = "figure1_draft.jpg", plot = g_fig1, path = fig_dir, 
       width = 88, height = 115, units = "mm", dpi = 1000)










# Figure 3: predicted versus observed phenotypic values -------------------


## Load prediction data
file_list <- list.files(result_dir, pattern = "fact_reg.RData$", full.names = TRUE)
object_list <- unlist(lapply(file_list, load, envir = .GlobalEnv))

## A vector to rename models
model_replace <- c("model1" = "y = g", "model2_id" = "y = g + E", "model2_cov" = "y = g + e", 
                   "model3_id" = "y = g + e + gE", "model3_cov" = "y = g + e + ge")

f_model_replace <- function(x) model_replace[x]
f_pop_replace <- function(x) str_replace_all(x, c("tp" = "CV", "vp" = "POV"))
# Replace type
f_type_replace <- function(x) c("loeo" = "New environment", "lolo" = "New location", "loyo" = "New year")[x]
# Replace ec selection
f_ec_selection_replace <- function(x) 
  c("adhoc" = "italic(ad~hoc)", "adhoc_nosoil" = "italic(ad~hoc)~(no~soil)", "apriori" = "italic(a~priori)")[x]

# Color scheme for models
model_colors <- c(neyhart_palette("umn1")[1], neyhart_palette("umn3")[3], neyhart_palette("umn1")[3],
                  neyhart_palette("umn3")[4], neyhart_palette("umn1")[4])
# names(model_colors) <- grep(pattern = "_", x = names(model_replace), value = TRUE, invert = TRUE)
names(model_colors) <- names(model_replace)



## Grab the prediction outputs
loo_prediction_list <- map(set_names(object_list, object_list), get) %>%
  # Bind rows if necessary
  modify_if(is.list, bind_rows) %>%
  subset(., map_lgl(., ~nrow(.) > 1)) %>%
  map(~unnest(., out)) %>%
  set_names(x = ., nm = str_extract(names(.), "^[a-z]{4}"))


## Combine data.frames and mutate columns
loo_predictions_df <- loo_prediction_list %>%
  imap(~unnest(.x, prediction) %>% mutate(type = .y) ) %>%
  map(~rename_at(.x, vars(which(names(.x) %in% c("env", "loc"))),
                 ~str_replace_all(., c("loc" = "location", "env" = "environment")))) %>%
  map_df(~mutate_if(., is.character, parse_guess) %>%
           mutate_if(is.factor, ~parse_guess(as.character(.)))) %>%
  mutate(pop = ifelse(line_name %in% tp, "tp", "vp")) %>%
  select(-which(names(.) %in% c(".id", "core", "trait1"))) %>%
  # Coalesce columns
  mutate(leave_one_group = coalesce(environment, location),
         nGroup = coalesce(nEnv, nLoc)) %>%
  select(-which(names(.) %in% c("environment", "location", "nLoc", "nEnv", "loc1", "env1")))


## Calculate accuracy and bias per train group, model, and population
loo_predictive_ability <- loo_predictions_df %>%
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  # First calculate accuracy per environment
  mutate(ability = cor(pred_complete, value), 
         bias = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  group_by(trait, model, pop, type, selection) %>%
  # Next calculate accuracy across all environments
  mutate(ability_all = cor(pred_complete, value), 
         bias_all = mean((pred_complete - value) / value)) %>% # Bias as percent deviation from observed
  # Now summarize across all
  group_by(trait, model, pop, type, leave_one_group, selection) %>%
  summarize_at(vars(ability, bias, ability_all, bias_all, nObs, nGroup), mean) %>%
  ungroup()



# First create annotation df
loo_prediction_accuracy_annotation <- loo_predictive_ability %>% 
  distinct(trait, model, pop, type, selection, ability_all, bias_all) %>%
  mutate_at(vars(contains("ability")), ~formatC(., width = 3, digits = 2, format = "f")) %>%
  mutate(ability_all_annotation = paste0("r[MP]==", ability_all))


# Plot predicted versus observed value using the best model and selection for each trait
best_results <- loo_prediction_accuracy_annotation %>% 
  filter(model == "model3_cov", pop == "tp") %>%
  group_by(trait, type) %>%
  top_n(x = ., n = 1, wt = ability_all) %>%
  ungroup() %>%
  distinct(trait, model, type, selection) %>%
  left_join(., subset(loo_prediction_accuracy_annotation, model == "model3_cov"))


## Plot - combine by trait
loo_pred_obs_tp_list <- loo_predictions_df %>%
  filter(pop == "tp") %>%
  inner_join(., best_results) %>%
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  mutate(trait = paste0("'", str_add_space(trait), "'~(", trait_units1[trait], ")")) %>%
  split(.$trait) %>%
  map(~{
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = distinct(.x, trait, type, ability_all_annotation),
                aes(x = Inf, y = -Inf, label = ability_all_annotation), inherit.aes = FALSE, 
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(trait ~ type, switch = "y", labeller = labeller(trait = label_parsed, type = f_type_replace)) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank(), panel.spacing.x = unit(1, "line"))
  })

# Combine the plots
loo_pred_obs_tp <- loo_pred_obs_tp_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "hv") %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6)

# Save
ggsave(filename = "figure2_draft.jpg", plot = loo_pred_obs_tp, path = fig_dir,
       height = 120, width = 88, units = "mm", dpi = 1000)

# Alternative
loo_pred_obs_tp1 <- loo_pred_obs_tp_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  .[-4] %>%
  plot_grid(plotlist = ., ncol = 1, align = "hv") %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6)

ggsave(filename = "figure2_draft1.jpg", plot = loo_pred_obs_tp1, path = fig_dir,
       height = 120, width = 88, units = "mm", dpi = 1000)



## VP ##

## Plot - combine by trait
loo_pred_obs_vp_list <- loo_predictions_df %>%
  filter(pop == "vp") %>%
  inner_join(., best_results) %>%
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  mutate(trait = paste0("'", str_add_space(trait), "'~(", trait_units1[trait], ")")) %>%
  split(.$trait) %>%
  map(~{
    ggplot(.x, aes(x = pred_complete, y = value, color = leave_one_group)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text(data = distinct(.x, trait, type, ability_all_annotation),
                aes(x = Inf, y = -Inf, label = ability_all_annotation), inherit.aes = FALSE, 
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_color_discrete(guide = FALSE) +
      facet_grid(trait ~ type, switch = "y", labeller = labeller(trait = label_parsed, type = f_type_replace)) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank(), panel.spacing.x = unit(1, "line"))
  })

# Combine the plots
loo_pred_obs_vp <- loo_pred_obs_vp_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "hv") %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6)

# Save
ggsave(filename = "figure3_draft.jpg", plot = loo_pred_obs_vp, path = fig_dir,
       height = 120, width = 88, units = "mm", dpi = 1000)

# Alternative
loo_pred_obs_vp1 <- loo_pred_obs_vp_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  .[-4] %>%
  plot_grid(plotlist = ., ncol = 1, align = "hv") %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6)

ggsave(filename = "figure3_draft1.jpg", plot = loo_pred_obs_vp1, path = fig_dir,
       height = 120, width = 88, units = "mm", dpi = 1000)



## Plot VP within TP
loo_pred_obs_combine_list <- loo_predictions_df %>%
  inner_join(., best_results) %>%
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  mutate(trait = paste0("'", str_add_space(trait), "'~(", trait_units1[trait], ")")) %>%
  split(.$trait) %>%
  map(~{
    
    ## Determine x and y axis breaks and limits
    xlims <- range(.x$pred_complete)
    ylims <- range(.x$value)
    xbreaks <- pretty(xlims)
    ybreaks <- pretty(ylims)
    
    ## Define the limits for the plot-in-plot (topleft)
    pp_scale <- 0.33
    pp_xs <- c(xlims[1] * 1.01, (xlims[1] * 1.01) + diff(xlims) * pp_scale)
    pp_ys <- rev(c(ylims[2] * 0.99, (ylims[2] * 0.99) - diff(ylims) * pp_scale))
    
    # Plot modifier
    g_mod <- list(
      scale_x_continuous(name = "Predicted phenotypic value", breaks = xbreaks, limits = xlims),
      scale_y_continuous(name = "Observed phenotypic value", breaks = ybreaks, limits = ylims),
      scale_color_discrete(guide = FALSE),
      theme_genetics(base_size = 6),
      theme(strip.placement = "outside", strip.background = element_blank(), axis.title = element_blank())
    )
    
    ## Split by type
    g_list_type <- split(.x, .x$type) %>%
      map(.x = ., .f = function(x1) {
        
        # Plot modifier
        
        # Plot the VP plot-in-plot window
        g_vp <- filter(x1, pop == "vp") %>%
          ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
          geom_point(size = 0.5, alpha = 0.5) +
          geom_text(data = filter(x1, pop == "vp") %>% distinct(trait, type, ability_all_annotation),
                    aes(x = Inf, y = -Inf, label = ability_all_annotation), inherit.aes = FALSE, 
                    parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
          g_mod
        
        ## Plot tp
        g_tp <- filter(x1, pop == "tp") %>%
          ggplot(aes(x = pred_complete, y = value, color = leave_one_group)) +
          geom_point(size = 0.5, alpha = 0.5) +
          geom_text(data = distinct(x1, trait, type, ability_all_annotation),
                    aes(x = Inf, y = -Inf, label = ability_all_annotation), inherit.aes = FALSE, 
                    parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
          facet_grid(trait ~ type, switch = "y", labeller = labeller(trait = label_parsed, type = f_type_replace)) +
          g_mod
        
        ## Add custom annotation
        g_tp + annotation_custom(grob = ggplotGrob(g_vp), xmin = pp_xs[1], xmax = pp_xs[2], 
                                 ymin = pp_ys[1], ymax = pp_ys[2])
        
      })
    
    ## Combine these plots
    # First modify so the second does not have a y strip or axis line
    g_combine <- g_list_type %>%
      modify_at(2, ~. + theme(strip.text.y = element_blank(), axis.line.y = element_blank(), 
                              axis.ticks.y = element_blank(), axis.text.y = element_blank()))  %>%
      plot_grid(plotlist = ., nrow = 1, align = "hv")
    
    # Return
    g_combine
    
  })
        
 
# Combine the plots
loo_pred_obs_combine <- loo_pred_obs_combine_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "hv") %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6)

# Save
ggsave(filename = "figure2_draft2.jpg", plot = loo_pred_obs_combine, path = fig_dir,
       height = 120, width = 88, units = "mm", dpi = 1000)

# Alternative
loo_pred_obs_combine1 <- loo_pred_obs_combine_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  .[-4] %>%
  plot_grid(plotlist = ., ncol = 1, align = "hv") %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6)

ggsave(filename = "figure2_draft2b.jpg", plot = loo_pred_obs_combine1, path = fig_dir,
       height = 130, width = 88, units = "mm", dpi = 1000)





# Figure 4: prediction accuracy within environments and locations ---------


## Summarize the best selection predictions
loo_predictive_ability_summ <- loo_predictive_ability %>%
  inner_join(., distinct(best_results, trait, type, selection)) %>%
  group_by(trait, model, pop, type, selection) %>%
  mutate_at(vars(ability, bias), list(mean = ~mean, sd = ~sd, n = ~n())) %>%
  ungroup() %>%
  mutate_at(vars(contains("sd")), list(se = ~. / sqrt(ability_n))) %>%
  select(trait, model, pop, type, leave_one_group, selection, ability, ability_mean, 
         ability_se = ability_sd_se) %>%
  unite(group, model, pop, remove = FALSE)


## Determine the LSD between model-selection groups
loo_prediction_accuracy_summ1 <- loo_predictive_ability %>%
  filter(!is.na(ability)) %>%
  mutate(type = str_extract(type, "l[a-z]{3}")) %>%
  group_by(type, trait, pop) %>%
  do({
    dat <- .
    fit <- lm(ability ~ model + selection + model:selection, data = dat)
    n <- max(as.numeric(xtabs(~ model + selection, dat)))
    MS_within <- sigma(fit)^2 # Within-group error (residuals)
    df <- df.residual(fit)
    LSD <- qt(p = 1 - (0.05 / 2), df = df) * sqrt(MS_within * (2 / n))
    
    # Return mean and LSD of prediction accuracy
    aggregate(ability ~ model + selection, dat, mean) %>% 
      mutate(LSD = LSD)
  }) %>% ungroup()





## Plot
g_accuracy_within_env <- loo_predictive_ability_summ %>%
  distinct(trait, model, pop, type, ability_mean, ability_se) %>%
  unite(group, model, pop, remove = FALSE) %>%
  ggplot(aes(x = trait, y = ability_mean, color = model, shape = pop, group = group)) +
  geom_point(position = position_dodge(0.9)) +
  scale_x_discrete(name = NULL, labels = str_add_space) +
  scale_y_continuous(name = expression('Predictive ability'~r[MP]), breaks = pretty) +
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace) +
  scale_shape_discrete(name = "Validation\nscheme", labels = f_pop_replace) +
  facet_grid(~ type, labeller = labeller(type = f_type_replace)) +
  theme_genetics(base_size = 6) +
  theme(strip.background = element_blank(), panel.border = element_rect(fill = alpha("white", 0)))


## Only plot models with covariates
g_accuracy_within_env2 <- g_accuracy_within_env %>% 
  modify_at("data", ~filter(., str_detect(model, "id", negate = TRUE)))

## Only plot model3_cov
g_accuracy_within_env3 <- g_accuracy_within_env %>% 
  modify_at("data", ~filter(., model == "model3_cov"))




g_accuracy_within_env <- loo_predictive_ability_summ %>%
  distinct(trait, model, pop, type, group, ability_mean, ability_se) %>%
  ggplot(aes(x = trait, y = ability_mean, group = group)) +
  # add points with jitter
  geom_jitter(data = loo_predictive_ability_summ, aes(y = ability), position = position_dodge(0.9),
              size = 0.15, alpha = 0.75, shape = 1) +
  geom_col(aes(fill = model), position = position_dodge(0.9), color = "black", lwd = 0.25, alpha = 0.75) +
  geom_errorbar(aes(ymin = ability_mean - ability_se, ymax = ability_mean + ability_se), 
                position = position_dodge(0.9), width = 0.5) +
  scale_x_discrete(name = NULL, labels = function(x) str_replace_all(str_add_space(x), " ", "\n")) +
  scale_y_continuous(name = expression('Predictive ability'~r[MP]), breaks = pretty) +
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)) +
  scale_shape_discrete(name = "Validation\nscheme", labels = f_pop_replace) +
  facet_grid(type ~ pop, switch = "y", labeller = labeller(type = f_type_replace, pop = f_pop_replace)) +
  theme_genetics(base_size = 6) +
  theme(strip.background = element_blank(), panel.border = element_rect(fill = alpha("white", 0)),
        strip.placement = "outside", axis.line = element_blank(),
        legend.text.align = 1, legend.position = c(0.89, 0.94), legend.key.size = unit(0.6, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5))


## Only plot models with covariates
g_accuracy_within_env2 <- g_accuracy_within_env %>% 
  modify_at("data", ~filter(., str_detect(model, "id", negate = TRUE))) 
g_accuracy_within_env2$layers[[1]]$data <- g_accuracy_within_env2$layers[[1]]$data %>%
  filter(., str_detect(model, "id", negate = TRUE))

## Save
ggsave(filename = "figure4_draft2.jpg", plot = g_accuracy_within_env2, path = fig_dir,
       width = 88, height = 90, dpi = 1000, units = "mm")










# Supplemental figures -------------------


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


## Plot
g_concurrent_features_count <- concurrent_features_count %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_count$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_count$feature_class == "main_features", 25, 10), size = 3) +
  facet_grid(trait ~ source, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Feature selection method") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "top", strip.placement = "outside")

ggsave(filename = "concurrent_features_count.jpg", plot = g_concurrent_features_count,
       path = fig_dir, height = 8, width = 4, dpi = 1000)


## Save as a table
concurrent_features_count %>% 
  spread(feature_class, n) %>% 
  arrange(source, feat_sel_type, trait) %>%
  write_csv(x = ., path = file.path(fig_dir, "concurrent_environmental_covariate.csv"))

# Polish as a supplemental table





#### Compare among data sources ####



## Overlap between sources for each trait and feature selection type
concurrent_features_overlap_source <- concurrent_features1 %>%
  select(-features) %>%
  full_join(., ., by = c("trait", "feat_sel_type")) %>%
  left_join(rename_all(as_tibble(t(combn(x = unique(concurrent_features1$source), m = 2))), ~paste0("source", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), trait, feat_sel_type, contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap"))


## Plot
g_concurrent_features_overlap_source <- concurrent_features_overlap_source %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  unite(source_pair, source.x, source.y, sep = ":") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_overlap_source$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_overlap_source$feature_class == "main_features", 30, 10), size = 3) +
  facet_grid(trait ~ source_pair, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Feature selection method") +
  scale_y_continuous(name = "Number of overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside")

ggsave(filename = "concurrent_features_overlap_source.jpg", plot = g_concurrent_features_overlap_source,
       path = fig_dir, height = 8, width = 3, dpi = 1000)


## Overlap between non-apriori / all feature selection types
concurrent_features_overlap_featsel <- concurrent_features1 %>%
  select(-features) %>%
  filter(! feat_sel_type %in% c("apriori", "all")) %>%
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
  unite(feat_sel_type_pair, feat_sel_type.x, feat_sel_type.y, sep = ":") %>%
  ggplot(aes(x = feat_sel_type_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_overlap_featsel$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_overlap_featsel$feature_class == "main_features", -8, -10), size = 3) +
  facet_grid(trait ~ source, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Feature selection pair") +
  scale_y_continuous(name = "Number of overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside")

ggsave(filename = "concurrent_features_overlap_featsel.jpg", plot = g_concurrent_features_overlap_featsel,
       path = fig_dir, height = 8, width = 2.5, dpi = 1000)


## Overlap between traits
concurrent_features_overlap_trait <- concurrent_features1 %>%
  select(-features) %>%
  full_join(., ., by = c("feat_sel_type", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(concurrent_features1$trait), m = 2))), 
                        ~paste0("trait", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap")) %>%
  # We don't need to see apriori or all
  filter(! feat_sel_type %in% c("apriori", "all"))



## Plot
g_concurrent_features_overlap_trait <- concurrent_features_overlap_trait %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("trait")), ~abbreviate(str_add_space(.), 2)) %>%
  unite(trait_pair, trait.x, trait.y, sep = ":") %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(concurrent_features_overlap_trait$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(concurrent_features_overlap_trait$feature_class == "main_features", -7, -10), size = 2) +
  facet_grid(feat_sel_type ~ source, switch = "y", labeller = labeller(feat_sel_type = f_ec_selection_replace)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Number of overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "concurrent_features_overlap_trait.jpg", plot = g_concurrent_features_overlap_trait,
       path = fig_dir, height = 3, width = 8, dpi = 1000)


## Stitch these plots together - patchwork
left_plot <- (g_concurrent_features_count + theme(legend.position = "none")) / g_concurrent_features_overlap_trait +
  plot_layout(heights = c(1, 0.3), guides = "collect")
right_plot <- g_concurrent_features_overlap_source | g_concurrent_features_overlap_featsel

merged_plot <- (left_plot | right_plot) + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1, 0.80))

# Save
ggsave(filename = "concurrent_features_comparison_merged.jpg", plot = merged_plot,
       path = fig_dir, height = 10, width = 10, dpi = 1000)


## Count the number of times a particular covariate is select

concurrent_features2 <- concurrent_features1 %>%
  select(-features) %>%
  gather(feature_class, covariates, contains("features")) %>%
  unnest() %>%
  mutate(covariates = str_remove(covariates, "line_name:"))
  

# Counts by source and feat_sel_type
concurrent_indiv_feature_counts <- concurrent_features2 %>%
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
  facet_grid(. ~ source, switch = "y") +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(legend.position = "top", strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Save
ggsave(filename = "concurrent_indiv_feature_counts.jpg", plot = g_concurrent_indiv_feature_counts, 
       path = fig_dir, width = 8, height = 4, dpi = 1000)






#### Do not compare among data sources ####

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
ggsave(filename = "concurrent_features_comparison_merged+daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = 1000)








#### Historical ####


## Combine concurrent feature selection df
historical_features <- bind_rows(
  historical_fact_reg_feature_selection %>% select(source, trait, model, apriori, stepAIC_adhoc = adhoc) %>% 
    gather(feat_sel_type, features, apriori, stepAIC_adhoc),
  gather(historical_feature_selection, feat_subset, features, adhoc, adhoc_nosoil) %>%
    unite(feat_sel_type, feat_sel_type, feat_subset, sep = "_")
) %>% mutate(model = case_when(model == "model2" ~ "model4", model == "model3" ~ "model5", TRUE ~ model))

## combine model4 and model5 covariates
historical_features1 <- historical_features %>%
  filter(str_detect(feat_sel_type, "adhoc_nosoil", negate = TRUE)) %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model4, model5, union)) %>%
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
historical_features_count <- historical_features1 %>% 
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  filter(feature_class != "features")


## Plot
g_historical_features_count <- historical_features_count %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_count$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_count$feature_class == "main_features", 25, 10), size = 3) +
  facet_grid(trait ~ source, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Feature selection method") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "top", strip.placement = "outside")

ggsave(filename = "historical_features_count.jpg", plot = g_historical_features_count,
       path = fig_dir, height = 8, width = 4, dpi = 1000)





## Overlap between sources for each trait and feature selection type
historical_features_overlap_source <- historical_features1 %>%
  select(-features) %>%
  full_join(., ., by = c("trait", "feat_sel_type")) %>%
  left_join(rename_all(as_tibble(t(combn(x = unique(historical_features1$source), m = 2))), ~paste0("source", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), trait, feat_sel_type, contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap"))


## Plot
g_historical_features_overlap_source <- historical_features_overlap_source %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  unite(source_pair, source.x, source.y, sep = ":") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_overlap_source$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_overlap_source$feature_class == "main_features", 30, 10), size = 3) +
  facet_grid(trait ~ source_pair, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Feature selection method") +
  scale_y_continuous(name = "Number of overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside")

ggsave(filename = "historical_features_overlap_source.jpg", plot = g_historical_features_overlap_source,
       path = fig_dir, height = 8, width = 3, dpi = 1000)


## Overlap between non-apriori / all feature selection types
historical_features_overlap_featsel <- historical_features1 %>%
  select(-features) %>%
  filter(! feat_sel_type %in% c("apriori", "all")) %>%
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
  unite(feat_sel_type_pair, feat_sel_type.x, feat_sel_type.y, sep = ":") %>%
  ggplot(aes(x = feat_sel_type_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_overlap_featsel$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_overlap_featsel$feature_class == "main_features", 3, 1), size = 3) +
  facet_grid(trait ~ source, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Feature selection pair") +
  scale_y_continuous(name = "Number of overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside")

ggsave(filename = "historical_features_overlap_featsel.jpg", plot = g_historical_features_overlap_featsel,
       path = fig_dir, height = 8, width = 2.5, dpi = 1000)


## Overlap between traits
historical_features_overlap_trait <- historical_features1 %>%
  select(-features) %>%
  full_join(., ., by = c("feat_sel_type", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(historical_features1$trait), m = 2))), 
                        ~paste0("trait", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap")) %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap")) %>%
  # We don't need to see apriori or all
  filter(! feat_sel_type %in% c("apriori", "all"))



## Plot
g_historical_features_overlap_trait <- historical_features_overlap_trait %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("trait")), ~abbreviate(str_add_space(.), 2)) %>%
  unite(trait_pair, trait.x, trait.y, sep = ":") %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(y = 2*max(historical_features_overlap_trait$n), label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)), 
            nudge_y = ifelse(historical_features_overlap_trait$feature_class == "main_features", -3, -4), size = 2) +
  facet_grid(feat_sel_type ~ source, switch = "y", labeller = labeller(feat_sel_type = f_ec_selection_replace)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Number of overlapping covariates", breaks = pretty) +
  theme_genetics() +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "historical_features_overlap_trait.jpg", plot = g_historical_features_overlap_trait,
       path = fig_dir, height = 3, width = 8, dpi = 1000)


## Stitch these plots together - patchwork
left_plot <- (g_historical_features_count + theme(legend.position = "none")) / g_historical_features_overlap_trait +
  plot_layout(heights = c(1, 0.3), guides = "collect")
right_plot <- g_historical_features_overlap_source | g_historical_features_overlap_featsel

merged_plot <- (left_plot | right_plot) + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1, 0.80))

# Save
ggsave(filename = "historical_features_comparison_merged.jpg", plot = merged_plot,
       path = fig_dir, height = 10, width = 10, dpi = 1000)



## Count the number of times a particular covariate is select

historical_features2 <- historical_features1 %>%
  select(-features) %>%
  gather(feature_class, covariates, contains("features")) %>%
  unnest() %>%
  mutate(covariates = str_remove(covariates, "line_name:"))


# Counts by source and feat_sel_type
historical_indiv_feature_counts <- historical_features2 %>%
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
  facet_grid(. ~ source, switch = "y") +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(legend.position = "top", strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Save
ggsave(filename = "historical_indiv_feature_counts.jpg", plot = g_historical_indiv_feature_counts, 
       path = fig_dir, width = 8, height = 4, dpi = 1000)

















## Calculate environmental relationship matrices based on these features
concurrent_features_env_relmat <- concurrent_features1 %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate(interaction_features_ec_mat = map(interaction_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE]),
         main_features_ec_mat = map(main_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE])) %>%
  mutate_at(vars(ends_with("ec_mat")), ~map(., ~Env_mat(x = .x, method = "Jarq")))

concurrent_features_env_relmat1 <- concurrent_features_env_relmat %>%
  select(trait, feat_sel_type, interaction = interaction_features_ec_mat, main = main_features_ec_mat) %>%
  gather(covariate_type, Emat, main, interaction) %>%
  ## If a matrix is all NA, convert to diagonal
  mutate(Emat = modify_if(Emat, ~all(is.na(.)), ~`diag<-`(ifelse(is.na(.), 0, 1), 1)))


## Compare the relationship between training/test and external environments for 
## each trait and feature selection type
concurrent_features_env_relmat1 %>%
  mutate(train_val_relat = map_dbl(Emat, ~mean(.[train_test_env, validation_env])),
         train_test_relat = ) %>% 
  arrange(covariate_type, trait, feat_sel_type) %>% 
  select(-Emat) %>% 
  View



# Define heat colors
heat_colors <- wesanderson::wes_palette("Zissou1", n = 5)[c(1,3,5)]


## Plot heatmaps of relationship matrices
# Separate plots by main/int covariate types
concurrent_features_heatmap_plots <- concurrent_features_env_relmat1 %>%
  group_by(trait, feat_sel_type, covariate_type) %>%
  do(plot = {
    row <- .
    
    # Create the heatmap to order the environments
    row_heat <- heatmap(x = row$Emat[[1]])
    
    # Factor order of environments
    env_order <- fct_inorder(row.names(row$Emat[[1]])[row_heat$rowInd])
    
    # Create the plotting data.frame
    dat <- row$Emat[[1]] %>%
      as.data.frame(.) %>% 
      rownames_to_column(., "environment") %>%
      gather(environment2, relationship, -environment) %>%
      # Refactor the environments
      mutate_at(vars(contains("environment")), ~factor(., levels = levels(env_order))) %>%
      mutate_at(vars(contains("environment")), list(group = ~ifelse(. %in% train_test_env, "training", "external"))) %>%
      mutate(environment_group = factor(environment_group, levels = c("training", "external")))
    
    # Plot
    g_heat <- dat %>%
      ggplot(aes(x = environment, y = environment2, fill = relationship)) +
      geom_tile() +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3]) +
      facet_grid(environment2_group ~ environment_group, scales = "free", space = "free") +
      labs(subtitle = paste(str_add_space(row$trait), row$feat_sel_type, sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = "none")
    
    # Return the plot
    g_heat
    
  }) %>% ungroup()

## Create plots by trait
for (plotList in split(concurrent_features_heatmap_plots, concurrent_features_heatmap_plots$trait)) {
  # Create plot
  plot_to_save <- plot_grid(plotlist = plotList$plot, nrow = n_distinct(plotList$feat_sel_type),
                            labels = c("Int.", "Main"), label_x = 0)
  
  # File name
  filename <- paste0("environment_covariate_relationship_", unique(plotList$trait), ".jpg")
  ggsave(filename = filename, plot = plot_to_save, path = fig_dir, 
         width = 8, height = 14, dpi = 1000)
  
}








# Supplemental Table XX - full model phenotypic variance analysis ---------

# Load the full model analysis results
load(file.path(result_dir, "full_model_variance_analysis.RData"))

# Collapse if a list
if (is.list(pheno_variance_analysis)) {
  pheno_variance_analysis <- pheno_variance_analysis %>% 
    subset(., sapply(., is.data.frame)) %>%
    bind_rows()
}

# Unnest results
pheno_variance_analysis1 <- pheno_variance_analysis %>%
  # Make a note as to whether any GxE covariates were included
  mutate(any_gxe_cov = map_lgl(features, ~any(str_detect(., "line_name:"))) | feat_sel_type == "all") %>%
  unnest(results)

## Calculate the total variance and calculate the proportion explained by each source
pheno_variance_analysis2 <- pheno_variance_analysis1 %>%
  group_by(trait, population, feat_sel_type) %>%
  mutate(total_variance = sum(variance), prop_total_variance = variance / total_variance) %>%
  rename(prop_source_variance = variance_prop) %>%
  ungroup()

## Clean up for a table
pheno_variance_analysis_table <- pheno_variance_analysis2 %>%
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
pheno_variance_analysis_table %>%
  rename_all(make.names) %>%
  select(-Any.Gxe.Cov, -Prop.Total.Variance) %>%
  filter(Population == "FP", str_detect(Term, "Covariates")) %>%
  {full_join(x = filter(., Covariate.Set == "Stepwise"), y = filter(., Covariate.Set != "Stepwise"),
             by = c("Trait", "Population", "Source", "Term"))} %>%
  mutate_at(vars(contains("Prop")), parse_number) %>%
  mutate(prop_diff = Prop.Source.Variance.x - Prop.Source.Variance.y)


# Output table
write_csv(x = pheno_variance_analysis_table, path = file.path(fig_dir, "full_model_phenotypic_variance_analysis.csv"))



























