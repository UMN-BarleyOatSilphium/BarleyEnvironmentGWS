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
# Subfigure label size
subfigure_label_size <- 9

# Set resolution
dpi_use <- 2000


# Define heat colors
heat_colors <- wesanderson::wes_palette("Zissou1", n = 5)[c(1,3,5)]

# Base font size
base_font_size <- 6
# Base geom text size
base_geom_text_size <- base_font_size * (5/14)

# 1 column figure width
figwidth_onecol <- 5.5
figwidth_twocol <- 12.0


## Create location-based color palettes
location_colors <- paletteer_d(package = "ggsci", palette = "default_igv", n = n_distinct(trial_info$location) + 2)
location_colors <- setNames(object = (location_colors[-c(4, 18)]), unique(trial_info$location))



# Use shading to designate environment-specific colors
environment_colors <- distinct(trial_info, location, environment) %>% 
  left_join(., tibble(location = names(location_colors), color = location_colors)) %>%
  arrange(location, environment) %>%
  group_by(location) %>%
  mutate(alpha = seq(n())) %>%
  ungroup() %>%
  mutate(alpha = alpha - 1, alpha = 1 - (alpha * 0.1), 
         color = alpha(color, alpha = alpha)) %>% 
  {setNames(pull(., color), pull(., environment))}



###
### Figures and Supplemental information is created in order of reference 
### in the manuscript 
### 



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

canada <- rnaturalearth::ne_states(country = "canada") %>%
  tidy(x = ., region = "name_en") %>%
  mutate(group = as.numeric(as.factor(group)))

# Download state data
usa_state <- as_tibble(map_data(map = "state"))

# Adjust the groups in the states
usa_state <- usa_state %>%
  mutate(group = group + max(canada$group))


# Tidy and combine
north_america <- bind_rows(usa_state, canada)


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
g_map_v1 <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "grey85") +
  geom_polygon(data = canada, fill = NA, color = "white", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "white", lwd = 0.3) +
  geom_point(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, shape = set), size = 3.5) +
  geom_text_repel(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = location),
                  size = 1.75, hjust = 0.5, nudge_x = ifelse(use_loc_info_toplot1$location == "Bozeman", 3, -1), segment.size = 0.2, 
                  point.padding = unit(2, "pt"), min.segment.length = 1) +
  geom_text(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = nExperiments), size = 2, 
            color = ifelse(use_loc_info_toplot1$location == "Arlington", "black", "white")) +
  coord_map(xlim = long_limit, ylim = lat_limit) +
  scale_shape_manual(name = NULL, values = c(16, 15), guide = guide_legend(label.position = "left", override.aes = list(size = 2))) +
  scale_x_continuous(breaks = NULL, name = NULL, labels = NULL) + 
  scale_y_continuous(breaks = NULL, name = NULL, labels = NULL) +
  theme_classic(base_size = 8) +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = alpha("white", 0)), axis.line = element_blank())





## Add visualization of genetic relationships using principal coordinate analysis

# First perform multi-dimensional scaling
K_mds <- cmdscale(d = dist(K), k = 2, eig = TRUE)

# Get the variance explained by each PCo
K_pco_varprop <- K_mds$eig %>% 
  {. / sum(.)} %>% 
  {. * 100} %>% 
  round(2) %>% 
  formatC(., digits = 2, format = "f") %>%
  str_c(str_c("PCo", seq_along(.)), " (", ., "%)") %>%
  set_names(str_c("PCo", seq_along(.)))

# Convert to DF and annotate
K_pco_df <- K_mds$points %>% 
  as.data.frame() %>% 
  `row.names<-`(., row.names(K)) %>%
  `colnames<-`(., str_c("PCo", seq(ncol(.)))) %>%
  rownames_to_column("line_name") %>% 
  left_join(., select(entry_list, line_name = Line, class = Class, program = Program, parent = `Parent?`)) %>%
  arrange(desc(class), program, line_name) %>%
  mutate(class = ifelse(class == "S2TP", "tp", "vp"), program = fct_inorder(program),
         parent = case_when(parent == "TRUE" ~ "parent", parent == "FALSE" ~ "nonparent", TRUE ~ "offspring"))

# Color scheme for populations
pop_colors <- setNames(object = c(neyhart_palette("umn1", n = 4), neyhart_palette("umn2")[c(5, 8)]), nm = levels(K_pco_df$program))

# Create the visualization
g_popstr <- K_pco_df %>% 
  ggplot(aes(x = PCo1, y = PCo2, fill = program, color = program)) + 
  # Plot non-parents and offspring
  geom_point(data = filter(K_pco_df, parent != "parent"), 
             shape = ifelse(subset(K_pco_df, parent != "parent")$parent == "nonparent", 16, 17), size = 0.7) +
  # Plot parents
  geom_point(data = filter(K_pco_df, parent == "parent"), shape = 21, color = "black", size = 0.7) +
  # scale_x_continuous(name = K_pco_varprop[1], breaks = pretty, limits = c(min(K_pco_df$PCo1)*1.65, max(K_pco_df$PCo1))) +
  scale_x_continuous(name = K_pco_varprop[1], breaks = pretty) +
  scale_y_continuous(name = K_pco_varprop[2], breaks = pretty, position = "right") +
  scale_color_manual(name = NULL, values = pop_colors, drop = FALSE) +
  scale_fill_manual(guide = FALSE, values = pop_colors, drop = FALSE) +
  theme_genetics(base_size = base_font_size) +
  theme(legend.position = "none", plot.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)),
        panel.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)))

# Save this
ggsave(filename = "figure1_partB_population_structure.jpg", plot = g_popstr, path = fig_dir,
       height = 2.5, width = 3, dpi = dpi_use)


# New limits
long_limit <- c(-123, -63)
lat_limit <- c(30, 50)

## Create a df of states with population colors
usa_state1 <- usa_state %>%
  mutate(fill = case_when(
    region == "idaho" ~ pop_colors["AB"],
    region == "montana" ~ pop_colors["MT"],
    region == "north dakota" ~ pop_colors["N2"],
    region == "washington" ~ pop_colors["WA"],
  ))


## Different projection
g_map_v2 <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, fill = NA, color = "grey85", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state1, aes(x = long, y = lat, group = group), fill = alpha(usa_state1$fill, 0.75), color = "grey85", lwd = 0.3) +
  geom_point(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, shape = set), size = 2.75) +
  geom_text_repel(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = location),
                  size = base_geom_text_size, hjust = 0.5, nudge_x = ifelse(use_loc_info_toplot1$location == "Buffalo County", -3, 0), 
                  segment.size = 0.2) +
  geom_text(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = nExperiments), 
            size = base_geom_text_size, color = ifelse(use_loc_info_toplot1$location == "Arlington", "black", "white")) +
  coord_map(projection = "bonne", lat0 = 50, xlim = long_limit, ylim = lat_limit) +
  scale_shape_manual(name = NULL, values = c(16, 15), guide = guide_legend(label.position = "left", override.aes = list(size = 2))) +
  scale_x_continuous(breaks = NULL, name = NULL, labels = NULL) + 
  scale_y_continuous(breaks = NULL, name = NULL, labels = NULL) +
  theme_classic(base_size = base_font_size) +
  theme(panel.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)), panel.grid = element_blank(),
        panel.border = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)), axis.line = element_blank(),
        plot.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)), legend.position = c(0.20, 0.20),
        legend.text.align = 1,
        legend.text = element_text(size = 5), legend.key.height = unit(0.5, "lines"), legend.key.width = unit(0.25, "lines"))

# Save this
ggsave(filename = "figure1_partA_draft2.jpg", plot = g_map_v2, path = fig_dir,
       height = 5, width = 8.7, units = "cm", dpi = dpi_use)
    
          
## Combine plots
layout <- c(
  area(t = 1, l = 1, b = 10, r = 8),
  area(t = 3, l = 7, b = 10, r = 9)
)

g_pop_map <- plot_grid(g_map_v2, labels = subfigure_labels[1], label_size = subfigure_label_size) + 
  plot_grid(g_popstr, labels = subfigure_labels[2], label_size = subfigure_label_size, label_x = 0.9)+ 
  plot_layout(design = layout)

# Save this
ggsave(filename = "figure1_partAB_draft1.jpg", plot = g_pop_map, path = fig_dir,
       height = 6, width = 12, units = "cm", dpi = dpi_use)





## Plot examples of crop model outputs

# Load concurrent and historical crop model output
load(file.path(enviro_dir, "GrowthStaging/apsim_s2met_model_results_daymet.RData"))
load(file.path(enviro_dir, "GrowthStaging/apsim_historical_growth_model_results.RData"))


# Assign units to variables
covariate_rename_abbr <- c("mint" = "T[min]", "maxt" = "T[max]", "tmean" = "T[mean]", "gdd" = "GDD",
                           "water_balance" = "H[2]O Bal.", "radn" = "Radn.")
                           
covariate_rename <- c("mint" = "Min. temp.", "maxt" = "Max. temp.", "tmean" = "Mean temp.",
                      "gdd" = "Growing deg. days", "water_balance" = "Water balance", "radn" = "Solar radiation")

# Units
covariate_variable_unit <- setNames(c(rep("degree*C", 4), "mm", "MJ~m^-2"), names(covariate_rename_abbr))

## Function for renaming a covariate
f_covariate_replace2 <- function(x) paste0("'", covariate_rename[x], "'~(", covariate_variable_unit[x], ")")
f_covariate_replace <- function(x) paste0("'", covariate_rename[x], "'~(", covariate_variable_unit[x], ")")
# Function for renaming timeframe
f_timeframe_replace <- function(x) c("concurrent" = "In-season", "historical" = "Historical location average")[x]


## Plot Bozeman and Columbus (Wooster)
cgm_locations_plot <- c("Bozeman", "Wooster")
cgm_year_plot <- 2017
# Assign colors to location
cgm_locations_color <- setNames(location_colors[1:2], cgm_locations_plot)

## Replacements and colors for growth stages

growth_stage_color <- setNames(object = paletteer_d(package = "ggsci", palette = "default_igv")[c(28, 41, 1, 47, 34)],
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
covariates_plot <- c("mint", "water_balance")


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
  theme_genetics(base_size = base_font_size) +
  theme(strip.placement = "outside", strip.background = element_blank(), panel.spacing.y = unit(0.5, "line"),
        legend.position = "none", strip.text.x = element_text(size = 8))




# Yearly growth stage information; subsetted properly for historical data
cgm_toplot_growth_stage_yearly <- cgm_toplot_growth_stage_times %>%
  split(.$timeframe) %>%
  modify_at(., "historical", ~filter(., between(year, min(trial_info$year) - 3, min(trial_info$year) - 1))) %>%
  bind_rows() 


## Create summary data.frames for historical and concurrent growth stages
# First summarize the entire length of the growing season per location
cgm_toplot_growing_season <- cgm_toplot_growth_stage_yearly %>%
  group_by(type, timeframe, location) %>% 
  summarize(min = min(min), max = max(max), y = mean(y)) %>%
  ungroup()

## Calculate average growth stages
cgm_toplot_growth_stage_times_summary <- cgm_toplot_growth_stage_yearly %>%
  filter(timeframe == "historical") %>%
  group_by(type, timeframe, location, stage) %>%
  summarize_at(vars(min, max, y, dap), mean) %>%
  ungroup() %>%
  bind_rows(., filter(cgm_toplot_growth_stage_times, timeframe == "concurrent"))





# Plot segments for growth stages
# both concurrent and historical
g_growth_stages <- cgm_toplot_growth_stage_yearly %>%
  filter(timeframe == "historical") %>%
  ggplot(aes(x = dap, color = stage)) +
  # Individual-year stages
  geom_segment(aes(x = min, xend = max, y = y, yend = y), lwd = 4, alpha = 0.15) +
  # Average stages
  geom_segment(data = cgm_toplot_growth_stage_times_summary, 
               aes(x = min, xend = max, y = y, yend = y), lwd = 2) +
  # Lines for the entire growing season per location
  geom_segment(data = cgm_toplot_growing_season, aes(x = min, xend = max, y = y, yend = y),
               color = cgm_locations_color[cgm_toplot_growing_season$location], lwd = 0.5, inherit.aes = FALSE) +
  # Annotation for growth stages
  geom_text(data = subset(cgm_toplot_growth_stage_times_summary, y == 1),
            aes(x = (min + max) / 2, y = y - 1, label = f_growth_stage_replace(as.character(stage))), size = base_geom_text_size) +
  facet_grid(type ~ timeframe, labeller = labeller(type = str_add_space, timeframe = f_timeframe_replace), 
             switch = "y", scales = "free_x", space = "free_x") +
  scale_color_manual(values = growth_stage_color, guide = FALSE) +
  scale_y_continuous(limits = c(-0.5, 2.5)) +
  scale_x_continuous(breaks = pretty, name = "Days after planting date", position = "top") +
  theme_genetics(base_size = base_font_size) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank(),
        axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.length.x = unit(-base_font_size / 4, "pt"),
        strip.text.x = element_text(size = 6))







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
  ncol = 1, align = "v", axis = "lr", rel_heights = c(0.35, 1), labels = subfigure_labels[3:4], label_size = subfigure_label_size,
  label_y = c(1, 1.03), label_x = 0.01)


## Combine plots
g_fig1 <- ( g_pop_map / g_weather_stages ) + plot_layout(heights = c(0.75, 1))

# Save
ggsave(filename = "figure1_map_growth_stage_example.jpg", plot = g_fig1, path = fig_dir, 
       width = 12, height = 12, units = "cm", dpi = dpi_use)

# Save as vector-based
ggsave(filename = "figure1_map_growth_stage_example.svg", plot = g_fig1, path = fig_dir, 
       width = 12, height = 12, units = "cm", dpi = 300)












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
  mutate(trait2 = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"),
         trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))


## Filter for relevant cases
accuracy_within_env_toplot <- within_environment_prediction_accuracy %>%
  filter(type == "loeo") %>%
  filter(model %in% str_subset(models_use, "id", negate = TRUE), selection %in% c("none", "rfa_cv_adhoc")) %>%
  # Add individual points
  # mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), "\n(n = ", nGroup+1, ")")) %>%
  mutate(trait1 = abbreviate(str_add_space(trait), 2)) %>%
  unite(group, trait, model, pop, remove = FALSE) %>%
  # Add a test name for faceting
  mutate(type_facet_x = "'Within-environment'~italic(r[MG])")

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





# Split by trait and Plot
g_loo_pred_obs_list <- loo_pred_obs_df %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = environment_colors, guide = FALSE) +
      facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


# Combine the plots
g_loo_pred_obs <- g_loo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
                     rel_widths = c(0.05, 1))


# Save
ggsave(filename = "figure2_partA_draft.jpg", plot = g_loo_pred_obs, path = fig_dir,
       height = 120, width = 80, units = "mm", dpi = dpi_use)



# ## Special empty plot for presentations
# g_loo_pred_obs1 <- g_loo_pred_obs_list %>% 
#   map(~. + geom_blank()) %>% # Add blank geom
#   map(~modify_at(.x, "layers", ~.[c(1,4)])) %>%
#   modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
#   map(~. + theme(axis.title = element_blank())) %>%
#   plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
#   add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
#   plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
#             rel_widths = c(0.05, 1))
# 
# # Save
# ggsave(filename = "figure2_partA_blank.jpg", plot = g_loo_pred_obs1, path = fig_dir,
#        height = 120, width = 80, units = "mm", dpi = dpi_use)
#        
#  # Blank out the progeny validation
# g_loo_pred_obs1 <- g_loo_pred_obs_list %>% 
#   map(~. + geom_blank()) %>%
#   map(~{
#     x1 <- modify_at(.x, "data", ~filter(., pop == "tp"))
#     x1$layers[[3]]$data <- mutate(x1$layers[[3]]$data, annotation = ifelse(pop == "vp", "", annotation))
#     return(x1)
#   }) %>%
#   modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
#   map(~. + theme(axis.title = element_blank())) %>%
#   plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
#   add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
#   plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
#             rel_widths = c(0.05, 1))
# 
# # Save
# ggsave(filename = "figure2_partA_blank2.jpg", plot = g_loo_pred_obs1, path = fig_dir,
#        height = 120, width = 80, units = "mm", dpi = dpi_use)






## Part B - accuracy within environments - LOEO



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
g_accuracy_within_env1 <- plot_grid(
  textGrob(label = expression('Within-environment'~italic(r[MG])), rot = 90, gp = gpar(fontsize = base_font_size)), 
  g_accuracy_within_env, rel_widths = c(0.05, 1))  

## Save
ggsave(filename = "figure2_partB_draft.jpg", plot = g_accuracy_within_env1, path = fig_dir,
       width = 80, height = 40, dpi = dpi_use, units = "mm")



## Plot both accuracy across environments and within environments per trait

# calculate the character width of the y axis
y_axis_width <- max(nchar(pretty(loo_pred_obs_df$value))) 

# Split by trait and Plot
g_plotlist2 <- loo_pred_obs_df %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    # Predicted/observed phenotypes
    p1 <- ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty, labels = function(x) formatC(x, width = y_axis_width)) +
      scale_fill_manual(values = environment_colors, guide = FALSE) +
      facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.background = element_blank(), axis.title = element_blank())
    
    within_env_df <- subset(accuracy_within_env_toplot, trait == unique(.x$trait))
    within_env_df2 <- subset(accuracy_within_env_toplot2, trait == unique(.x$trait))
    
    # Add annotation for number of environments
    nE <- paste0("n = ", unique(within_env_df$nGroup) + 1)
    
      
    # Accuracy within environments
    p2 <- within_env_df %>%
      ggplot(aes(x = pop, y = accuracy, color = model, fill = model, group = group)) +
      # jitter points
      geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
                  size = 0.2, color = "grey85") +
      # boxplot
      geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
      annotate(geom = "text", label = nE, x = -Inf, y = Inf, vjust = 1, hjust = -0.25, size = base_geom_text_size*0.8) +
      scale_x_discrete(name = NULL, labels = f_validation_replace) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                        guide = guide_legend(label.position = "left", title = NULL)) +
      scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                         guide = guide_legend(label.position = "left", title = NULL)) +
      facet_grid(. ~ type_facet_x, labeller = label_parsed) +
      theme_genetics(base_size = 6) +
      theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
            legend.text.align = 1, legend.position = c(0.80, 0.95), legend.key.size = unit(0.4, "line"),
            legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
            legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line")) +
      theme(axis.title = element_blank(), strip.text.y = element_text(color = "white"))
    
    # Return a list of plots
    list(plot1 = p1, plot2 = p2)
    
  })

## Edit plots
# Remove strip text from all plots except the first; do the same for axis x text for plot2
g_figure2 <- g_plotlist2 %>%
  modify_at(-1, ~{
    modify_at(.x, "plot1", ~. + theme(strip.text.x = element_blank())) %>%
      modify_at("plot2", ~. + theme(legend.position = "none", strip.text.x = element_blank())) }) %>%
  modify_at(-length(.), ~{modify_at(.x, "plot2", ~.+theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_blank()))}) %>%
  unlist(recursive = FALSE) %>%
  # Combine plots
  plot_grid(plotlist = ., ncol = 2, align = "v", axis = "lr", rel_heights = c(1, rep(0.85, length(.) - 1)), rel_widths = c(1, 0.65),
            labels = subfigure_labels[1:2], label_size = subfigure_label_size) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, hjust = 1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.03, 1))


# ## Combine plots
# g_figure2 <- plot_grid(plotlist = g_plotlist2, 
#                        ncol = 1, rel_heights = c(1, 0.35), labels = subfigure_labels[1:2], label_size = 8)

## Save
ggsave(filename = "figure2_loeo_predictions.jpg", plot = g_figure2, path = fig_dir,
       width = figwidth_twocol, height = 12, dpi = dpi_use, units = "cm")

## Save as vector image
ggsave(filename = "figure2_loeo_predictions.svg", plot = g_figure2, path = fig_dir,
       width = figwidth_onecol, height = 12, dpi = dpi_use, units = "cm")







# Figure 3. LOLO predictions ---------------------------------

## Parts A and B: timeframe search results

# Load the search results
load(file.path(result_dir, "historical_covariate_timeframe_selection.RData"))
# Load the feature selection results
load(file.path(result_dir, "feature_selection_results.RData"))
## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


# Extract results
historical_timeframe_selection_out1 <- historical_timeframe_selection_out %>%
  filter(feat_sel_type == "adhoc") %>%
  mutate(stepwise_results = map(covariates, "finalResults") %>% map(as.list) %>% map(as_tibble)) %>%
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
  theme_genetics(base_size = base_font_size),
  theme(strip.placement = "outside", strip.background = element_blank())
)


## Plot timeframe results
historical_timeframe_selection_out1_toplot <- historical_timeframe_selection_out1 %>%
  filter(time_frame_type == "time_frame", model == "model5")

# Separate df for the selected timeframe (to add as points and an annotation)
timeframe_used_df <- inner_join(distinct(historical_feature_selection, trait, time_frame), historical_timeframe_selection_out1_toplot)

# Plot together
g_timeframe_analysis <- historical_timeframe_selection_out1_toplot %>%
  ggplot(aes(x = length, y = RMSE_scale, color = trait)) +
  geom_line(lwd = 0.25) +
  geom_point(data = timeframe_used_df, size = 0.3) +
  # # Add an annotation
  # annotate(geom = "curve", x = 10, y = 0.35, xend = min(timeframe_used_df$length) + 1, yend = min(timeframe_used_df$RMSE_scale),
  #          curvature = 0.25, arrow = arrow(angle = 30, length = unit(0.25, "line"), ends = "last"), lwd = 0.25) +
  # annotate(geom = "text", label = "Chosen timeframe", x = 10.5, y = 0.39, hjust = 1, size = 2) +
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Scaled RMSEP") +
  scale_x_continuous(name = "Historical weather data years (before 2015)", trans = "reverse") +
  scale_color_manual(values = trait_colors, name = NULL, labels = function(x) abbreviate(str_add_space(x), 2),
                     guide = guide_legend(nrow = 2)) +
  theme_genetics(base_size = base_font_size) +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.3, 0.15),
        legend.direction = "horizontal", legend.background = element_rect(fill = alpha("white", 0)),
        legend.key.width = unit(0.5, "line"), legend.key.height = unit(0.5, "line"))

  
# ## Plot window results
# ## 
# ## Make this supplemental?
# ## 
# historical_window_selection_out1_toplot <- historical_timeframe_selection_out1 %>%
#   filter(time_frame_type == "window", model == "model3") %>%
#   mutate(length = fct_inseq(as.factor(length))) %>%
#   mutate_at(vars(contains("year")), ~ymd(paste0(., "0101")))
# 
# 
# # Plot together
# g_window_analysis <- historical_window_selection_out1_toplot %>%
#   ggplot(aes(x = end_year, y = RMSE_scale, color = trait, lty = length)) +
#   geom_line(lwd = 0.25) +
#   scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Scaled RMSEP") +
#   scale_x_date(date_breaks = "5 year", date_labels = "%Y", name = "End year of window") +
#   scale_color_manual(values = trait_colors, guide = FALSE) +
#   scale_linetype_discrete(name = "Window length (yrs)", guide = guide_legend(title.position = "top")) +
#   theme_genetics(base_size = 6) +
#   theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.3, 0.15),
#         legend.direction = "horizontal", legend.key.width = unit(0.5, "line"), legend.key.height = unit(0.5, "line"))
# 





## Part B - accuracy within locations - LOLO

## Filter for relevant cases
accuracy_within_loc_toplot <- predictive_ability %>%
  ## Adjust model names
  mutate(model = str_replace_all(model, c("model2" = "model4", "model3" = "model5"))) %>%
  filter(type == "lolo") %>%
  filter(model %in% str_subset(models_use, "id", negate = TRUE), selection %in% c("none", "rfa_adhoc")) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), "\n(n = ", nGroup+1, ")")) %>%
  # mutate(trait1 = abbreviate(str_add_space(trait), 2)) %>%
  mutate(pop = paste0(f_validation_replace(pop), " (", f_pop_replace(pop), ")")) %>%
  unite(group, trait, model, pop, remove = FALSE)


# Plot point and line range
g_accuracy_within_loc <- accuracy_within_loc_toplot %>%
  ggplot(aes(x = trait1, y = ability, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.2, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = expression('Within-location'~italic(r[MP])), breaks = pretty) +
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)) +
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title = NULL)) +
  facet_grid(. ~ pop, switch = "y") +
  theme_genetics(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
        legend.text.align = 1, legend.position = c(0.89, 0.95), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line")) +
  theme(axis.title.x = element_blank(), strip.text.y = element_text(color = "white"))

## Save
ggsave(filename = "figure3_partB_draft.jpg", plot = g_accuracy_within_loc, path = fig_dir,
       width = 80, height = 40, dpi = dpi_use, units = "mm")



## Part C - plot predicted/observed phenotypes for each trait and by target population
## for model 3 with stepwise covariates - LOLO

lolo_pred_obs_df <- predictions_df %>%
  filter(type == "lolo") %>%
  filter(model == "model3_cov", selection == "rfa_adhoc") %>%
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
    
    ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = location_colors, guide = FALSE) +
      facet_grid(pop ~ trait1, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_pop_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


# Combine the plots
g_lolo_pred_obs <- g_lolo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.y = element_blank(), axis.title.y = element_blank())) %>%
  map(~. + theme(axis.title.x = element_blank())) %>%
  plot_grid(plotlist = ., nrow = 1, align = "h", rel_widths = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = base_font_size, vjust = -1)


# Save
ggsave(filename = "figure3_partC_draft.jpg", plot = g_lolo_pred_obs, path = fig_dir,
       height = 5, width = 16, units = "cm", dpi = dpi_use)



## Combine plots
# Align the left-most top plot with the bottom plot
plots <- align_plots(g_timeframe_analysis, g_lolo_pred_obs, align = "h", axis = "tblr")
# Build the top row
top_row <- plot_grid(plots[[1]], g_accuracy_within_loc, nrow = 1, align = "h", axis = "tb",
                     labels = subfigure_labels[1:2], label_size = subfigure_label_size, rel_widths = c(0.65, 1))
# Then combine with the bottom row
g_figure3 <- plot_grid(top_row, plots[[2]], ncol = 1, rel_heights = c(0.75, 1),
                       labels = c("", subfigure_labels[3]), label_size = subfigure_label_size)



## Save
ggsave(filename = "figure3_lolo_predictions.jpg", plot = g_figure3, path = fig_dir,
       width = figwidth_twocol, height = 10, dpi = dpi_use, units = "cm")

## Save as vector image
ggsave(filename = "figure3_lolo_predictions.svg", plot = g_figure3, path = fig_dir,
       width = figwidth_twocol, height = 10, dpi = dpi_use, units = "cm")











# Figure 4. External environment prediction accuracy ---------------------------------------


## Part A - plot predicted/observed phenotypes for each trait and by target population
## for model 3 with stepwise covariates

external_pred_obs_df <- predictions_df %>%
  filter(str_detect(type, "ext")) %>%
  filter(model == "model3_cov", str_detect(selection, "rfa"), str_detect(selection, "nosoil", negate = TRUE)) %>%
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
    ann_df <- distinct(select(.x, trait, pop, type,  test_group, ability)) %>% 
      mutate(test_group1 =  ifelse(str_detect(test_group, "[0-9]{2}"), test_group, 
                                                   str_to_title(str_replace_all(test_group, "_", " "))),
             ann1 = paste0("r[MP]==", formatC(signif(ability, 2), digits = 2, format = "f", flag = "0")),
             annotation = paste0(abbreviate(test_group1, 5), "~", ann1),
             x = ifelse(str_detect(type, "env") & trait %in% c("PlantHeight", "TestWeight"), -Inf, Inf), y = -Inf)
    
    ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
      geom_text_repel(data = ann_df, aes(x = x, y = y, label = annotation, color = test_group), inherit.aes = FALSE,
                      parse = TRUE, hjust = 1, size = 1, direction = "y", segment.alpha = 0, box.padding = 0.01) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = c(location_colors, environment_colors), guide = FALSE) +
      scale_color_manual(values = c(location_colors, environment_colors), guide = FALSE) +
      facet_grid(type ~ trait1, switch = "y", labeller = labeller(trait1 = label_parsed, type = f_type_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


## Combine and save
g_holdout_pred_obs <- g_external_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.y = element_blank(), axis.title.y = element_blank())) %>%
  map(~. + theme(axis.title.x = element_blank())) %>%
  plot_grid(plotlist = ., nrow = 1, align = "h", rel_widths = c(1, rep(0.80, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = base_font_size, vjust = -1)

## Vertical
g_holdout_pred_obs_vert <- g_external_pred_obs_list %>%
  map(~. + facet_grid(trait1 ~ type, switch = "y", labeller = labeller(trait1 = label_parsed, type = f_type_replace))) %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.05, 1))

# Save
ggsave(filename = "figure4_exernal_validation.jpg", plot = g_holdout_pred_obs_vert, path = fig_dir,
       height = 10, width = 5.5, units = "cm", dpi = dpi_use)

# Save as vector image
ggsave(filename = "figure4_exernal_validation.svg", plot = g_holdout_pred_obs_vert, path = fig_dir,
       height = 10, width = 5.5, units = "cm", dpi = dpi_use)












# Supplementary figures ---------------------------------------------------


# Figure S1: validation of predicted flowering time from crop model -------

# The code for this picture is in the 'define_environmental_covariables.R' script




# Figure S2: Comparison of environment-specific EC sets -------------------

## Number and overlap of covariates in the analyses ##


## Concurrent

## Combine concurrent feature selection df
concurrent_features <- bind_rows( concurrent_fact_reg_feature_selection, concurrent_feature_selection, concurrent_all_features ) %>%
  rename(features = covariates) %>% select(-direction)

## combine model2 and model3 covariates
concurrent_features1 <- concurrent_features %>%
  filter(str_detect(feat_sel_type, "adhoc_nosoil", negate = TRUE)) %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model2, model3, union)) %>%
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
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc", n > 0) %>%
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n))

concurrent_features2 <- concurrent_features1 %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc")


## Modify the feature count plot to remove nasapower
g_concurrent_features_count2 <- concurrent_features_count1 %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  filter(source == "daymet") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = "top", strip.placement = "outside")




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
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(~ feat_sel_type_pair, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))


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
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(feat_sel_type ~ ., switch = "y", labeller = labeller(feat_sel_type = f_ec_selection_replace)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))



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
  mutate(feature_class = fct_inorder(feature_class)) %>%
  ggplot(aes(x = covariates, y = n, fill = feature_class)) +
  geom_col() +
  scale_fill_discrete(labels = c("main_features" = "Environment covariate", "interaction_features" = "G x E covariate"), name = NULL) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count") +
  theme_genetics(base_size = base_font_size) +
  theme(legend.position = c(0.75, 0.75), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## Stitch these plots together - patchwork
left_plot <- plot_grid(g_concurrent_features_count2 + theme(legend.position = "none"), g_concurrent_features_overlap_trait, 
                       ncol = 1, rel_heights = c(1, 0.45), align = "hv", axis = "lr", labels = subfigure_labels[1:2],
                       label_size = 10)  
right_plot <- plot_grid(g_concurrent_features_overlap_featsel, g_concurrent_indiv_feature_counts, ncol = 1,
                        rel_heights = c(0.6, 1), align = "hv", axis = "lr", labels = subfigure_labels[3:4], label_size = 10)  

merged_plot <- (left_plot | right_plot) + plot_layout(widths = c(1, 1))

# Save
ggsave(filename = "figure_S2_concurrent_features_comparison_merged_daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = dpi_use)






# Figure S3: Comparison of location-specific EC sets ----------------------

# Edit historical_feature_selection
historical_feature_selection1 <- historical_feature_selection %>%
  mutate(source = "daymet", feat_sel_type = str_replace(feat_sel_type, "rfa_adhoc", "rfa_cv_adhoc"))


## Combine historical feature selection df
historical_features <- bind_rows( historical_fact_reg_feature_selection, historical_all_features,  historical_feature_selection1 ) %>%
  rename(features = covariates) %>% select(-direction, -time_frame)
  
  

## combine model2 and model3 covariates
historical_features1 <- historical_features %>%
  filter(str_detect(feat_sel_type, "adhoc_nosoil", negate = TRUE)) %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model4, model5, union)) %>%
  mutate(features = map(features, ~setdiff(., "line_name"))) %>%
  select(-contains("model")) %>%
  mutate(interaction_features = map(features, ~str_subset(., ":")),
         main_features = map2(features, interaction_features, setdiff))


## Calculate the number of covariates in each contingency
historical_features_count <- historical_features1 %>% 
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  filter(feature_class != "features")



## Compare trait, selections, etc. for covariates

# These will be polished plots for the manuscript
historical_features_count1 <- historical_features_count %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc", n > 0) %>%
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n))

historical_features2 <- historical_features1 %>%
  filter(source == "daymet", feat_sel_type != "stepAIC_adhoc")


## Modify the feature count plot to remove nasapower
g_historical_features_count1 <- historical_features_count1 %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  filter(source == "daymet") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = "top", strip.placement = "outside")


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
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(~ feat_sel_type_pair, switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))



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
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(feat_sel_type ~ ., switch = "y", labeller = labeller(feat_sel_type = f_ec_selection_replace)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Overlapping covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))


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
  mutate(feature_class = fct_inorder(feature_class)) %>%
  ggplot(aes(x = covariates, y = n, fill = feature_class)) +
  geom_col() +
  scale_fill_discrete(labels = c("main_features" = "Environment covariate", "interaction_features" = "G x E covariate"), name = NULL) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count") +
  theme_genetics(base_size = base_font_size) +
  theme(legend.position = c(0.75, 0.75), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## Stitch these plots together - patchwork
left_plot <- plot_grid(g_historical_features_count1 + theme(legend.position = "none"), g_historical_features_overlap_trait, 
                       ncol = 1, rel_heights = c(1, 0.45), align = "hv", axis = "lr", labels = subfigure_labels[1:2],
                       label_size = 10)  
right_plot <- plot_grid(g_historical_features_overlap_featsel, g_historical_indiv_feature_counts, ncol = 1,
                        rel_heights = c(0.6, 1), align = "hv", axis = "lr", labels = subfigure_labels[3:4], label_size = 10)  

merged_plot <- (left_plot | right_plot) + plot_layout(widths = c(1, 1))

# Save
ggsave(filename = "figure_S3_historical_features_comparison_merged_daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = dpi_use)







# Supplementary fig XX - environmental relationships ----------------------


## Calculate environmental correlation

env_phenoCor <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  group_by(trait) %>%
  do(phenoCor = {
    select(., line_name, environment, value) %>%
      as.data.frame() %>%
      spread(environment, value) %>%
      column_to_rownames("line_name") %>%
      cor(use = "pairwise.complete.obs")
  }) %>% ungroup() %>%
  mutate(phenoCor_df = map(phenoCor, ~{
    as.data.frame(.x) %>%
      rownames_to_column("environment") %>%
      gather(environment2, phenoCorrelation, -environment)
  }))



## Calculate environmental means

## Fit models to calculate environmental means
env_means <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Factorize
    df1 <- mutate_at(df, vars(line_name, environment), ~fct_contr_sum(as.factor(.))) %>%
      mutate(weight = std_error^2)
    
    # Fit the model
    fit <- lmer(value ~ (1|line_name) + environment, data = df1, weights = weight)
    
    ## Return a df of environmental effects
    fixef(fit) %>% 
      tibble(environment = names(.), effect = .) %>% 
      filter(environment != "(Intercept)") %>% 
      mutate(environment = str_remove(environment, "environment"),
             mu = fixef(fit)[1]) %>% 
      add_row(environment = last(levels(df1$environment)), effect = -sum(.$effect), 
              mu = .$mu[1])
    
  }) %>% ungroup()

# Calculate the distance between environmental means per trait
env_mean_dist <- env_means %>%
  mutate(mean = mu + effect) %>%
  group_by(trait) %>%
  do(dist = as.matrix(dist(setNames(.$mean, .$environment)))) %>%
  ungroup()


# Matrix of covariates
ec_tomodel_scaled_mat <- ec_tomodel_scaled$daymet %>%
  select(-source) %>%
  filter(environment %in% c(train_test_env, validation_env)) %>%
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()

## Calculate environmental relationship matrices based on these features
concurrent_features_env_relmat <- concurrent_features1 %>%
  filter(source == "daymet") %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate(interaction_features_ec_mat = map(interaction_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE]),
         main_features_ec_mat = map(main_features, ~ec_tomodel_scaled_mat[, .x, drop = FALSE])) %>%
  mutate_at(vars(ends_with("ec_mat")), ~map(., ~Env_mat(x = .x, method = "Jarq")))

concurrent_features_env_relmat1 <- concurrent_features_env_relmat %>%
  select(trait, feat_sel_type, interaction = interaction_features_ec_mat, main = main_features_ec_mat) %>%
  # Set to 0 any that are completely NA
  mutate_at(vars(interaction, main), ~map(., ~if (all(is.na(.x))) {.x[] <- 0; .x} else .x)) %>%
  # Sum the covariance matrices
  mutate(covMat = map2(main, interaction, `+`))


## Calculate correlation matrices
concurrent_features_env_cormat <- concurrent_features_env_relmat1 %>%
  rename(corMat = covMat) %>%
  mutate_at(vars(corMat, main, interaction), ~map(., ~{
    D <- diag((diag(.x)^-0.5))
    .x1 <- D %*% .x %*% D
    `dimnames<-`(.x1, dimnames(.x))
  }))



## Plot the distance between environmental means versus the correlation due to covariates
env_mean_concurrent_features <- concurrent_features_env_cormat %>%
  left_join(., env_mean_dist) %>%
  # mutate(corMat = map2(corMat, dist, ~.x[row.names(.y), row.names(.y)]))
  mutate(corMat = map2(main, dist, ~.x[row.names(.y), row.names(.y)]))

env_mean_cor_concurrent_features <- env_mean_concurrent_features %>%
  mutate(corMat_envMean_df = map2(corMat, dist, ~{
    .x[lower.tri(.x, diag = TRUE)] <- NA
    .y[lower.tri(.y, diag = TRUE)] <- NA
    
    .x1 <- as.data.frame(.x) %>% 
      rownames_to_column("environment") %>% 
      gather(environment2, correlation, -environment) %>% 
      filter(!is.na(correlation))
    
    .y1 <- as.data.frame(.y) %>% 
      rownames_to_column("environment") %>% 
      gather(environment2, env_mean_dist, -environment) %>% 
      filter(!is.na(env_mean_dist))
    
    full_join(.x1, .y1, by = c("environment", "environment2"))
    
  })) %>% ungroup()

## Plot
g_env_mean_cor_concurrent_features <- env_mean_cor_concurrent_features %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(corMat_envMean_df) %>%
  mutate(env_mean_dist = ifelse(trait == "GrainYield", env_mean_dist / 1000, env_mean_dist)) %>%
  ggplot(aes(x = correlation, y = env_mean_dist)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Difference in environmental mean", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "concurrent_ec_correlation_vs_env_mean_dist.jpg", plot = g_env_mean_cor_concurrent_features,
       path = fig_dir, height = 5, width = 5, dpi = 1000)


# For each environment, find the environment with which it has the highest correlation
# magnitude
# Compare this to their difference in environmental means
env_mean_cor_concurrent_features1 <- env_mean_cor_concurrent_features %>%
  mutate(corMat_envMean_df = map(corMat_envMean_df, ~{
    envs <- unique(unlist(select(.x, 1:2)))
    map(envs, function(e) filter_at(.x, vars(contains("environment")), any_vars(. == e))) %>%
      map_df(~top_n(x = ., n = 1, wt = correlation))
  }))

# Plot
g_env_mean_cor_concurrent_features2 <- env_mean_cor_concurrent_features1 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(corMat_envMean_df) %>%
  mutate(env_mean_dist = ifelse(trait == "GrainYield", env_mean_dist / 1000, env_mean_dist)) %>%
  ggplot(aes(x = correlation, y = env_mean_dist)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Difference in environmental mean", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "concurrent_ec_correlation_vs_env_mean_dist_highestCor.jpg", plot = g_env_mean_cor_concurrent_features2,
       path = fig_dir, height = 5, width = 5, dpi = 1000)



## Plot the phenotypic correlation versus the correlation due to covariates
phenCor_concurrent_features <- concurrent_features_env_cormat %>%
  left_join(., env_phenoCor) %>%
  mutate(corMat = map2(corMat, phenoCor, ~.x[row.names(.y), row.names(.y)]))

phenCor_concurrent_features1 <- phenCor_concurrent_features %>%
  mutate(phenoCor_corMat = map2(corMat, phenoCor, ~{
    .x[lower.tri(.x, diag = TRUE)] <- NA
    .y[lower.tri(.y, diag = TRUE)] <- NA
    
    .x1 <- as.data.frame(.x) %>% 
      rownames_to_column("environment") %>% 
      gather(environment2, correlation, -environment) %>% 
      filter(!is.na(correlation))
    
    .y1 <- as.data.frame(.y) %>% 
      rownames_to_column("environment") %>% 
      gather(environment2, phenoCor, -environment) %>% 
      filter(!is.na(phenoCor))
    
    full_join(.x1, .y1, by = c("environment", "environment2"))
    
  })) %>% ungroup()

## Plot
g_phenCor_concurrent_features1 <- phenCor_concurrent_features1 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(phenoCor_corMat) %>%
  ggplot(aes(x = correlation, y = phenoCor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "concurrent_ec_correlation_vs_phenoCor.jpg", plot = g_phenCor_concurrent_features1,
       path = fig_dir, height = 5, width = 5, dpi = 1000)


# For each environment, find the environment with which it has the highest correlation
# magnitude
# Compare this to their difference in environmental means
phenCor_concurrent_features2 <- phenCor_concurrent_features1 %>%
  mutate(phenoCor_corMat = map(phenoCor_corMat, ~{
    envs <- unique(unlist(select(.x, 1:2)))
    map(envs, function(e) filter_at(.x, vars(contains("environment")), any_vars(. == e))) %>%
      map_df(~top_n(x = ., n = 1, wt = correlation))
  }))

# Plot
g_phenCor_concurrent_features2 <- phenCor_concurrent_features2 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(phenoCor_corMat) %>%
  ggplot(aes(x = correlation, y = phenoCor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "concurrent_ec_correlation_vs_phenoCor_highestCor.jpg", plot = g_phenCor_concurrent_features2,
       path = fig_dir, height = 5, width = 5, dpi = 1000)




## For each trait and covariate type (interaction/main), plot a heatmap and
## a dendrogram
concurrent_features_heatmap_plots <- concurrent_features_env_cormat %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  group_by(trait, feat_sel_type) %>%
  do(plot = {
    row <- .
    
    corMat <- row$corMat[[1]]

    # First perform clustering on the relationship matrix
    env_clust <- hclust(dist(corMat, method = "euclidean"), method = "average")
    # get the data
    clust_data <- dendro_data(model = env_clust)
    # Label data
    clust_lab_data <- mutate(clust_data$labels, group = ifelse(label %in% train_test_env, "Train/test", "Holdout"),
                             group = fct_relevel(group, "Holdout", after = Inf))

    # Factor order of environments
    env_order <- levels(clust_lab_data$label)
    
    # Re-orient the correlation matrix
    corMat <- corMat1 <- corMat[env_order, env_order]
    # Set the lower triangle to NA
    corMat1[upper.tri(corMat)] <- NA

    # Create the plotting data.frame
    dat1 <- corMat1 %>%
      as.data.frame(.) %>%
      rownames_to_column(., "environment") %>%
      gather(environment2, correlation, -environment)
    
    ## Add correlation annotation for the lower triangle
    dat <- corMat %>%
      as.data.frame(.) %>%
      rownames_to_column(., "environment") %>%
      gather(environment2, correlation, -environment) %>%
      mutate(annotation = format_numbers(correlation, 2))

    # Merge
    dat2 <- left_join(dat1, select(dat, -correlation), by = c("environment", "environment2")) %>%
      # Refactor the environments
      mutate_at(vars(contains("environment")), ~factor(., levels = env_order)) %>%
      mutate_at(vars(contains("environment")), list(group = ~ifelse(. %in% train_test_env, "training", "external"))) %>%
      mutate(annotation = ifelse(is.na(correlation), annotation, NA))
    

    # Plot
    g_heat <- dat2 %>%
      ggplot(aes(x = environment, y = environment2)) +
      # geom_tile(aes(fill = correlation), color = ifelse(is.na(dat2$correlation), "black", NA), lwd = 0.1) +
      geom_tile(aes(fill = correlation), lwd = 0.1) +
      # geom_text(aes(label = annotation), size = 1) +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3], na.value = "white",
                           midpoint = 0, limits = c(-1.01, 1.01), breaks = pretty,
                           name = "Estimated environmental\ncorrelation", guide = guide_colorbar(title.position = "top")) +
      scale_y_discrete(position = "right") +
      labs(subtitle = paste(str_add_space(row$trait), f_ec_selection_replace(row$feat_sel_type), sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = c(0.25, 0.75),
            axis.text.y = element_blank(), legend.direction = "horizontal")

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
      theme(legend.position = c(0.8, 0.90), legend.key.height = unit(5, "pt"))

    # Combine plots
    plot_grid(g_heat, g_dendro, align = "v", axis = "tb", rel_widths = c(1, 0.7))

  }) %>% ungroup()

## Combine all plots
g_heat_dendro_all <- plot_grid(plotlist = concurrent_features_heatmap_plots$plot,
                               ncol = n_distinct(concurrent_features_heatmap_plots$feat_sel_type))
# Save
ggsave(filename = "figure_sXX_environmental_covariate_heatmap.jpg", plot = g_heat_dendro_all,
       path = fig_dir, width = 20, height = 20, dpi = 1000)






# Supplementary fig XX - location relationships ----------------------


# Calculate genotype-location means

# Data to use - location means
S2_MET_loc_BLUEs <- S2_MET_BLUEs %>%
  # Remove irrigated trials - these will eventually be included
  filter(!str_detect(environment, "HTM|BZI|AID")) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup()


# Calculate location phenotypic correlation

loc_phenoCor <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp) %>%
  group_by(trait) %>%
  do(phenoCor = {
    select(., line_name, location, value) %>%
      as.data.frame() %>%
      spread(location, value) %>%
      column_to_rownames("line_name") %>%
      cor(use = "pairwise.complete.obs")
  }) %>% ungroup() %>%
  mutate(phenoCor_df = map(phenoCor, ~{
    as.data.frame(.x) %>%
      rownames_to_column("location") %>%
      gather(location2, phenoCorrelation, -location)
  }))



## Calculate environmental means

## Fit models to calculate environmental means
loc_means <- S2_MET_loc_BLUEs %>%
  filter(line_name %in% tp) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Factorize
    df1 <- mutate_at(df, vars(line_name, location), ~fct_contr_sum(as.factor(.))) 
    
    # Fit the model
    fit <- lmer(value ~ (1|line_name) + location, data = df1)
    
    ## Return a df of environmental effects
    fixef(fit) %>% 
      tibble(location = names(.), effect = .) %>% 
      filter(location != "(Intercept)") %>% 
      mutate(location = str_remove(location, "location"),
             mu = fixef(fit)[1]) %>% 
      add_row(location = last(levels(df1$location)), effect = -sum(.$effect), 
              mu = .$mu[1])
    
  }) %>% ungroup()

# Calculate the distance between environmental means per trait
loc_mean_dist <- loc_means %>%
  mutate(mean = mu + effect) %>%
  group_by(trait) %>%
  do(dist = as.matrix(dist(setNames(.$mean, .$location)))) %>%
  ungroup()


# Matrix of covariates
loc_tomodel_scaled_mat_list <- historical_ec_tomodel_timeframe_scaled %>%
  subset(., grepl(pattern = "daymet", x = names(.))) %>%
  setNames(object = ., nm = gsub(pattern = "daymet\\.", replacement = "", x = names(.))) %>%
  map(~select(.x, -source, -time_frame) %>%
        filter(location %in% c(train_test_loc, validation_loc)) %>%
        as.data.frame() %>%
        column_to_rownames("location") %>%
        as.matrix() )

## Calculate environmental relationship matrices based on these features
historical_features_loc_relmat <- historical_features1 %>%
  left_join(., distinct(historical_feature_selection, trait, time_frame)) %>%
  filter(source == "daymet") %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate(interaction_features_ec_mat = map2(interaction_features, time_frame, ~{
    loc_tomodel_scaled_mat_list[[.y]][, .x, drop = FALSE] }),
    main_features_ec_mat = map2(main_features, time_frame, ~{
      loc_tomodel_scaled_mat_list[[.y]][, .x, drop = FALSE] })) %>%
  mutate_at(vars(ends_with("ec_mat")), ~map(., ~Env_mat(x = .x, method = "Jarq")))

historical_features_loc_relmat1 <- historical_features_loc_relmat %>%
  select(trait, feat_sel_type, interaction = interaction_features_ec_mat, main = main_features_ec_mat) %>%
  # Set to 0 any that are completely NA
  mutate_at(vars(interaction, main), ~map(., ~if (all(is.na(.x))) {.x[] <- 0; .x} else .x)) %>%
  # Sum the covariance matrices
  mutate(covMat = map2(main, interaction, `+`))


## Calculate correlation matrices
historical_features_loc_cormat <- historical_features_loc_relmat1 %>%
  rename(corMat = covMat) %>%
  mutate_at(vars(corMat, main, interaction), ~map(., ~{
    D <- diag((diag(.x)^-0.5))
    .x1 <- D %*% .x %*% D
    `dimnames<-`(.x1, dimnames(.x))
  }))



## Plot the distance between environmental means versus the correlation due to covariates
loc_mean_historical_features <- historical_features_loc_cormat %>%
  left_join(., loc_mean_dist) %>%
  # mutate(corMat = map2(corMat, dist, ~.x[row.names(.y), row.names(.y)]))
  mutate(corMat = map2(main, dist, ~.x[row.names(.y), row.names(.y)]))

loc_mean_cor_historical_features <- loc_mean_historical_features %>%
  mutate(corMat_locMean_df = map2(corMat, dist, ~{
    .x[lower.tri(.x, diag = TRUE)] <- NA
    .y[lower.tri(.y, diag = TRUE)] <- NA
    
    .x1 <- as.data.frame(.x) %>% 
      rownames_to_column("location") %>% 
      gather(location2, correlation, -location) %>% 
      filter(!is.na(correlation))
    
    .y1 <- as.data.frame(.y) %>% 
      rownames_to_column("location") %>% 
      gather(location2, loc_mean_dist, -location) %>% 
      filter(!is.na(loc_mean_dist))
    
    full_join(.x1, .y1, by = c("location", "location2"))
    
  })) %>% ungroup()

## Plot
g_loc_mean_cor_historical_features <- loc_mean_cor_historical_features %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(corMat_locMean_df) %>%
  mutate(loc_mean_dist = ifelse(trait == "GrainYield", loc_mean_dist / 1000, loc_mean_dist)) %>%
  ggplot(aes(x = correlation, y = loc_mean_dist)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Difference in location mean", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "historical_ec_correlation_vs_loc_mean_dist.jpg", plot = g_loc_mean_cor_historical_features,
       path = fig_dir, height = 5, width = 5, dpi = 1000)


# For each location, find the location with which it has the highest correlation
# magnitude
# Compare this to their difference in environmental means
loc_mean_cor_historical_features1 <- loc_mean_cor_historical_features %>%
  mutate(corMat_locMean_df = map(corMat_locMean_df, ~{
    locs <- unique(unlist(select(.x, 1:2)))
    map(locs, function(e) filter_at(.x, vars(contains("location")), any_vars(. == e))) %>%
      map_df(~top_n(x = ., n = 1, wt = correlation))
  }))

# Plot
g_loc_mean_cor_historical_features2 <- loc_mean_cor_historical_features1 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(corMat_locMean_df) %>%
  mutate(loc_mean_dist = ifelse(trait == "GrainYield", loc_mean_dist / 1000, loc_mean_dist)) %>%
  ggplot(aes(x = correlation, y = loc_mean_dist)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Difference in location mean", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))


# Save
ggsave(filename = "historical_ec_correlation_vs_loc_mean_dist_highestCor.jpg", plot = g_loc_mean_cor_historical_features2,
       path = fig_dir, height = 5, width = 5, dpi = 1000)



## Plot the phenotypic correlation versus the correlation due to covariates
phenCor_historical_features <- historical_features_loc_cormat %>%
  left_join(loc_phenoCor) %>%
  mutate(corMat = map2(corMat, phenoCor, ~.x[row.names(.y), row.names(.y)]))

phenCor_historical_features1 <- phenCor_historical_features %>%
  mutate(phenoCor_corMat = map2(corMat, phenoCor, ~{
    .x[lower.tri(.x, diag = TRUE)] <- NA
    .y[lower.tri(.y, diag = TRUE)] <- NA
    
    .x1 <- as.data.frame(.x) %>% 
      rownames_to_column("location") %>% 
      gather(location2, correlation, -location) %>% 
      filter(!is.na(correlation))
    
    .y1 <- as.data.frame(.y) %>% 
      rownames_to_column("location") %>% 
      gather(location2, phenoCor, -location) %>% 
      filter(!is.na(phenoCor))
    
    full_join(.x1, .y1, by = c("location", "location2"))
    
  })) %>% ungroup()

## Plot
g_phenCor_historical_features1<- phenCor_historical_features1 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(phenoCor_corMat) %>%
  ggplot(aes(x = correlation, y = phenoCor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "historical_ec_correlation_vs_phenoCor.jpg", plot = g_phenCor_historical_features1,
       path = fig_dir, height = 5, width = 5, dpi = 1000)


# For each environment, find the environment with which it has the highest correlation
# magnitude
# Compare this to their difference in environmental means
phenCor_historical_features2 <- phenCor_historical_features1 %>%
  mutate(phenoCor_corMat = map(phenoCor_corMat, ~{
    envs <- unique(unlist(select(.x, 1:2)))
    map(envs, function(e) filter_at(.x, vars(contains("location")), any_vars(. == e))) %>%
      map_df(~top_n(x = ., n = 1, wt = correlation))
  }))

# Plot
g_phenCor_historical_features2 <- phenCor_historical_features2 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(phenoCor_corMat) %>%
  ggplot(aes(x = correlation, y = phenoCor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Correlation estimated from environmental covariates", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)))

# Save
ggsave(filename = "historical_ec_correlation_vs_phenoCor_highestCor.jpg", plot = g_phenCor_historical_features2,
       path = fig_dir, height = 5, width = 5, dpi = 1000)




## For each trait and covariate type (interaction/main), plot a heatmap and
## a dendrogram
historical_features_heatmap_plots <- historical_features_loc_cormat %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  group_by(trait, feat_sel_type) %>%
  do(plot = {
    row <- .
    
    corMat <- row$corMat[[1]]
    
    # First perform clustering on the relationship matrix
    env_clust <- hclust(dist(corMat, method = "euclidean"), method = "average")
    # get the data
    clust_data <- dendro_data(model = env_clust)
    # Segment data
    clust_seg_data <- segment(clust_data)
    # Label data
    clust_lab_data <- mutate(clust_data$labels, group = ifelse(label %in% train_test_loc, "Train/test", "Holdout"),
                             group = fct_relevel(group, "Holdout", after = Inf),
                             label1 = str_to_title(str_replace_all(label, "_", " ")),
                             label1 = ifelse(str_detect(label1, " ") & nchar(label1) > 10, str_replace(label1, " ", "\n"), label1))
    
    # Factor order of environments
    env_order <- levels(clust_lab_data$label)
    
    # Re-orient the correlation matrix
    corMat <- corMat1 <- corMat[env_order, env_order]
    # Set the lower triangle to NA
    corMat1[upper.tri(corMat)] <- NA
    
    # Create the plotting data.frame
    dat1 <- corMat1 %>%
      as.data.frame(.) %>%
      rownames_to_column(., "location") %>%
      gather(location2, correlation, -location)
    
    ## Add correlation annotation for the lower triangle
    dat <- corMat %>%
      as.data.frame(.) %>%
      rownames_to_column(., "location") %>%
      gather(location2, correlation, -location) %>%
      mutate(annotation = format_numbers(correlation, 2))
    
    # Factor order of environments
    env_order <- clust_lab_data$label1
    
    # Merge
    dat2 <- left_join(dat1, select(dat, -correlation), by = c("location", "location2")) %>%
      left_join(., select(clust_lab_data, label, label1), c("location" = "label")) %>%
      left_join(., select(clust_lab_data, label, label1), c("location2" = "label")) %>%
      select(location = label1.x, location2 = label1.y, correlation, annotation) %>%
      # Refactor the environments
      mutate_at(vars(contains("location")), ~factor(., levels = env_order)) %>%
      mutate_at(vars(contains("location")), list(group = ~ifelse(. %in% train_test_env, "training", "external"))) %>%
      mutate(annotation = ifelse(is.na(correlation), annotation, NA))


    
    
    # Plot
    g_heat <- dat2 %>%
      ggplot(aes(x = location, y = location2)) +
      # geom_tile(aes(fill = correlation), color = ifelse(is.na(dat2$correlation), "black", NA), lwd = 0.1) +
      geom_tile(aes(fill = correlation), lwd = 0.1) +
      # geom_text(aes(label = annotation), size = 1) +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3], na.value = "white",
                           midpoint = 0, limits = c(-1.01, 1.01), breaks = pretty,
                           name = "Estimated location\ncorrelation", guide = guide_colorbar(title.position = "top")) +
      scale_y_discrete(position = "right") +
      labs(subtitle = paste(str_add_space(row$trait), f_ec_selection_replace(row$feat_sel_type), sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = c(0.25, 0.75),
            axis.text.y = element_blank(), legend.direction = "horizontal")
    
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
      theme(legend.position = c(0.8, 0.90), legend.key.height = unit(5, "pt"))
    
    # Combine plots
    plot_grid(g_heat, g_dendro, align = "h", axis = "tb", rel_widths = c(1, 0.7))
    
  }) %>% ungroup()

## Combine all plots
g_heat_dendro_all <- plot_grid(plotlist = historical_features_heatmap_plots$plot,
                               ncol = n_distinct(historical_features_heatmap_plots$feat_sel_type))
# Save
ggsave(filename = "figure_sXX_location_covariate_heatmap.jpg", plot = g_heat_dendro_all,
       path = fig_dir, width = 20, height = 20, dpi = 1000)
























# Figure S4 - LOEO predictive ability across environments ----------

## Plot predicted vs observed phenotypic values for LOEO for all models, traits,
## EC selections, etc.

# First create a data.frame from which to plot all scenarios

loo_pred_obs_df <- predictions_df %>%
  # Filter for relevant models and remove nonsoil selections
  filter(model %in% models_use, str_detect(selection, "nosoil", negate = TRUE)) %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # Add environment/location specific predictive abilities
  left_join(., select(predictive_ability, trait:ability)) %>%
  ## Adjust rfa_cv selections
  mutate(selection = str_replace(selection, "rfa_adhoc", "rfa_cv_adhoc")) %>%
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
    
    ggplot(data = .x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.2) +
      geom_point(size = 0.1, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = environment_colors, guide = FALSE) +
      facet_grid(pop ~ model + selection, switch = "y", scales = "free_x", labeller = labeller(pop = f_pop_replace)) +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_loeo_all <- plot_grid(plotlist = g_loeo_all_list, ncol = 1, align = "v", axis = "lr",
                        labels = subfigure_labels[seq_along(g_loeo_all_list)], label_size = 8)


# Save
ggsave(filename = "figure_S4_loeo_across_site_prediction_accuracy.jpg", plot = g_loeo_all, path = fig_dir,
       height = 20, width = 22, units = "cm", dpi = 1000)



















# Figure S5 - LOEO prediction accuracy and bias within environments ----------

## Filter for relevant cases
accuracy_bias_within_env_toplot <- within_environment_prediction_accuracy %>%
  filter(type == "loeo") %>%
  filter(model %in% models_use, str_detect(selection, "nosoil", negate = TRUE)) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         selection = fct_inorder(selection),
         model = str_replace(model, "_id", "_cov"),
         bias = bias * 100) %>%
  unite(group, trait, model, pop, remove = FALSE) %>%
  select(-rmse, -varR, -heritability, -ability) %>%
  gather(statistic, value, accuracy, bias)

accuracy_bias_within_env_toplot2 <- accuracy_bias_within_env_toplot %>%
  group_by(group, type, trait, trait1, model, pop, selection, statistic) %>%
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
ggsave(filename = "figure_S5_loeo_accuracy_bias_within_environments.jpg", plot = g_accuracy_bias_within_env, 
       path = fig_dir, width = 12, height = 4, dpi = 1000)




# Figure S6 - LOLO predictive ability across locations ----------

## Plot predicted vs observed phenotypic values for LOLO for all models, traits,
## EC selections, etc.

## Plot all
g_lolo_all_list <- loo_pred_obs_df %>%
  filter(type == "lolo") %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, model, selection, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(data = .x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.2) +
      geom_point(size = 0.1, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = location_colors, guide = FALSE) +
      facet_grid(pop ~ model + selection, switch = "y", scales = "free_x", labeller = labeller(pop = f_pop_replace)) +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_lolo_all <- plot_grid(plotlist = g_lolo_all_list, ncol = 1, align = "v", axis = "lr",
                        labels = subfigure_labels[seq_along(g_lolo_all_list)], label_size = 8)

# Save
ggsave(filename = "figure_S6_lolo_across_site_prediction_accuracy.jpg", plot = g_lolo_all, path = fig_dir,
       height = 20, width = 22, units = "cm", dpi = 1000)



# Figure S7 - LOLO prediction accuracy and bias within locations ----------

## Filter for relevant cases
accuracy_bias_within_loc_toplot <- predictive_ability %>%
  filter(type == "lolo") %>%
  filter(model %in% models_use, str_detect(selection, "nosoil", negate = TRUE)) %>%
  # Add individual points
  mutate(selection = str_replace(selection, "rfa_adhoc", "rfa_cv_adhoc"),
         trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         selection = fct_inorder(selection),
         model = str_replace(model, "_id", "_cov"),
         bias = bias * 100) %>%
  unite(group, trait, model, pop, remove = FALSE) %>%
  select(-contains("all"), -rmse) %>%
  gather(statistic, value, ability, bias)

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
  filter(statistic == "ability") %>%
  ggplot(aes(x = selection, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  # Pointrange
  geom_linerange(data = filter(accuracy_bias_within_loc_toplot2, statistic == "ability"), 
                 aes(ymin = lower, ymax = upper), position = position_dodge(0.9), lwd = 0.25) +
  geom_linerange(data = filter(accuracy_bias_within_loc_toplot2, statistic == "ability"), 
                 aes(ymin = q25, ymax = q75), position = position_dodge(0.9), lwd = 0.6) +
  geom_point(data = filter(accuracy_bias_within_loc_toplot2, statistic == "ability"), 
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
ggsave(filename = "figure_S7_lolo_accuracy_bias_within_locations.jpg", plot = g_accuracy_bias_within_loc, 
       path = fig_dir, width = 12, height = 4, dpi = 1000)



# Figure S8 - external predictive ability across environments ----------

## Plot all
g_ext_env_all_list <- loo_pred_obs_df %>%
  filter(type == "env_external") %>%
  filter(str_detect(selection, "Concurrent|AIC", negate = TRUE)) %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(select(.x, trait, pop, type,  test_group, ability, model, selection)) %>% 
      mutate(test_group1 =  ifelse(str_detect(test_group, "[0-9]{2}"), test_group, 
                                        str_to_title(str_replace_all(test_group, "_", " "))),
             ann1 = paste0("r[MP]==", formatC(signif(ability, 2), digits = 2, format = "f", flag = "0")),
             annotation = paste0(abbreviate(test_group1, 5), "~", ann1)) %>%
      mutate(x = case_when(trait == "TestWeight" & selection == "StepwiseCV" ~ -Inf,
                           TRUE ~ Inf),
             y = case_when(trait == "GrainYield" & str_detect(selection, "None", negate = TRUE) ~ Inf,
                           TRUE ~ -Inf))
    
    ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.1, shape = 21, color = alpha("white", 0)) +
      geom_label_repel(data = ann_df, aes(x = x, y = y, label = annotation, color = test_group), inherit.aes = FALSE,
                       parse = TRUE, hjust = 1, size = 1, direction = "y", segment.alpha = 0, box.padding = 0.01, 
                       label.size = NA, label.padding = 0.01, fill = alpha("white", 0.5)) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = environment_colors, guide = FALSE) +
      scale_color_manual(values = environment_colors, guide = FALSE) +
      facet_grid(. ~ model + selection, switch = "y", scales = "free_x") +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_ext_env_all <- plot_grid(plotlist = g_ext_env_all_list, ncol = 1, align = "v", axis = "lr",
                           labels = subfigure_labels[seq_along(g_lolo_all_list)], label_size = 8)

# Save
ggsave(filename = "figure_S8_external_environment_across_site_prediction_accuracy.jpg", plot = g_ext_env_all, path = fig_dir,
       height = 16, width = 22, units = "cm", dpi = 1000)





# Figure S9 - external predictive ability across locations ----------

## Plot all
g_ext_loc_all_list <- loo_pred_obs_df %>%
  filter(type == "loc_external") %>%
  filter(str_detect(selection, "Concurrent|AIC", negate = TRUE)) %>%
  split(.$trait) %>%
  map(~{
    
    # Extract annotations
    ann_df <- distinct(select(.x, trait, pop, type,  test_group, ability, model, selection)) %>% 
      mutate(test_group1 =  ifelse(str_detect(test_group, "[0-9]{2}"), test_group, 
                                        str_to_title(str_replace_all(test_group, "_", " "))),
             ann1 = paste0("r[MP]==", formatC(signif(ability, 2), digits = 2, format = "f", flag = "0")),
             annotation = paste0(abbreviate(test_group1, 5), "~", ann1)) %>%
      mutate(x = case_when(trait == "GrainYield" & selection %in% c("All", "Literature") ~ -Inf,
                           TRUE ~ Inf),
             y = case_when(trait == "GrainYield" & selection == "All" ~ Inf,
                           TRUE ~ -Inf))
    
    ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.1, shape = 21, color = alpha("white", 0)) +
      geom_label_repel(data = ann_df, aes(x = x, y = y, label = annotation, color = test_group), inherit.aes = FALSE,
                       parse = TRUE, hjust = 1, size = 1, direction = "y", segment.alpha = 0, box.padding = 0.01, 
                       label.size = NA, label.padding = 0.01, fill = alpha("white", 0.5)) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = location_colors, guide = FALSE) +
      scale_color_manual(values = location_colors, guide = FALSE) +
      facet_grid(. ~ model + selection, switch = "y", scales = "free_x") +
      labs(subtitle = parse(text = unique(.x$trait1))) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())
    
  })

# Combine the plots
g_ext_loc_all <- plot_grid(plotlist = g_ext_loc_all_list, ncol = 1, align = "v", axis = "lr",
                           labels = subfigure_labels[seq_along(g_lolo_all_list)], label_size = 8)

# Save
ggsave(filename = "figure_S9_external_location_across_site_prediction_accuracy.jpg", plot = g_ext_loc_all, path = fig_dir,
       height = 16, width = 22, units = "cm", dpi = 1000)




# Figure S10 - comparison of genomic relationships between populations ----------

# Tidy the K matrix
genomic_relationship_df <- tidy(as.dist(K)) %>%
  rename_all(~c("line_name1", "line_name2", "genomic_rel"))

## Subset the K matrix for FP-FP comparisons and FP-OP comparisons
fp_fp_genomic_relationship <- filter_at(genomic_relationship_df, vars(contains("line_name")), all_vars(. %in% tp)) %>%
  mutate(comparison = "fp_fp")

fp_op_genomic_relationship <- filter(genomic_relationship_df, line_name1 %in% vp, line_name2 %in% tp) %>%
  mutate(comparison = "fp_op")

## Plot
bind_rows(fp_fp_genomic_relationship, fp_op_genomic_relationship) %>%
  ggplot(aes(x = genomic_rel, fill = comparison)) +
  geom_density(alpha = 0.5)




# Supplementary fig XX - 75-25 environment cross-validation ---------------

# Filter across-environment accuracy for LOEO and CV
environment_pred_acc_compare <- predictive_ability %>%
  filter(type %in% c("cv_env", "loeo"), str_detect(selection, "soil", negate = TRUE)) %>%
  distinct_at(vars(trait, type, cv_rep, pop, model, selection, contains("all"))) %>%
  # Modify the names
  mutate(type = str_replace_all(type, c("loeo" = "Leave-one-environment-out", "cv_env" = "75-25 environment CV"))) %>%
  # Summarize over reps
  group_by(trait, type, pop, model, selection) %>%
  summarize(ability_all_mean = mean(ability_all), ability_all_sd = sd(ability_all), n = n()) %>%
  ungroup() %>%
  mutate_at(vars(ends_with("sd")), list(se = ~. / sqrt(n)))


## Compare 75-25 cross-validation results to LOEO results
g_env_cv_loeo_compare <- environment_pred_acc_compare %>%
  # Filter for relevant models and variable selections
  filter(selection == "rfa_cv_adhoc", model %in% c("model2_cov", "model3_cov")) %>%
  ggplot(aes(x = pop, y = ability_all_mean, fill = type)) +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = ability_all_mean - se, ymax = ability_all_mean + se), position = position_dodge(0.9), width = 0.5) + 
  scale_fill_discrete(name = NULL) +
  scale_y_continuous(name = "Predictive ability", breaks = pretty) +
  scale_x_discrete(name = NULL, labels = f_validation_replace) +
  facet_grid(trait ~ model + selection, switch = "y", scales = "free_y",
             labeller = labeller(model = f_model_replace, selection = f_ec_selection_replace)) +
  facet_grid(trait ~ model, switch = "y", labeller = labeller(model = f_model_replace, trait = str_add_space)) +
  theme_genetics(base_size = 8) +
  theme(legend.position = "top")

ggsave(filename = "concurrent_compare_crossValidation_across_environment.jpg", plot = g_env_cv_loeo_compare,
       path = fig_dir, height = 6, width = 4, dpi = 1000)







# Supplementary tables ----------------------------------------




# Table S1: Variance component breakdown ----------------------------------

# Load the initial variance component partition
load(file.path(result_dir, "stage_two_phenotypic_analysis.RData"))

## Create a df of the range in environmental means
ge_env_mean_range <- stage_two_varcomp_env %>% 
  mutate(env_mean_range = map_chr(env_mean, ~paste(format_numbers(range(.$mean)), collapse = "-"))) %>%
  select(trait, population, env_mean_range)

## Create a table for the partition of genotype-environmental means
ge_var_comp <- stage_two_varcomp_env %>% 
  unnest(var_comp) %>%
  # Calculate proportions of variance
  group_by(trait, population) %>%
  mutate(prop_var = (variance / sum(variance)) * 100) %>%
  ungroup() %>%
  mutate_if(is.numeric, format_numbers) %>% 
  mutate(annotation = paste0(variance, " (", se, ") / ", prop_var, "%")) %>% 
  select(trait, population, term, annotation) %>% 
  spread(term, annotation) %>%
  ## Add environmental mean range
  left_join(., ge_env_mean_range)

## Adjust names and save the table
ge_var_comp %>%
  mutate(trait = str_add_space(trait), population = f_pop_replace(population)) %>%
  select(Trait = trait, Population = population, Genotype = line_name,
         Environment = environment, GxE = gxe, Residual = units, `Environmental Mean Range` = env_mean_range) %>%
  write_csv(x = ., path = file.path(fig_dir, "table_s1_environmental_variance_partition.csv"))


## Do the same thing for location

## Create a df of the range in environmental means
ge_loc_mean_range <- stage_two_varcomp_loc %>% 
  mutate(loc_mean_range = map_chr(loc_mean, ~paste(format_numbers(range(.$mean)), collapse = "-"))) %>%
  select(trait, population, loc_mean_range)

## Create a table for the partition of genotype-environmental means
ge_var_comp <- stage_two_varcomp_loc %>% 
  unnest(var_comp) %>%
  # Calculate proportions of variance
  group_by(trait, population) %>%
  mutate(prop_var = (variance / sum(variance)) * 100) %>%
  ungroup() %>%
  mutate_if(is.numeric, format_numbers) %>% 
  mutate(annotation = paste0(variance, " (", se, ") / ", prop_var, "%")) %>% 
  select(trait, population, term, annotation) %>% 
  spread(term, annotation) %>%
  ## Add environmental mean range
  left_join(., ge_loc_mean_range)

## Adjust names and save the table
ge_var_comp %>%
  mutate(trait = str_add_space(trait), population = f_pop_replace(population)) %>%
  select(Trait = trait, Population = population, Genotype = line_name,
         Location = location, GxL = gxe, Residual = units, `Location Mean Range` = loc_mean_range) %>%
  write_csv(x = ., path = file.path(fig_dir, "table_s1_location_variance_partition.csv"))



# Table S2: Name, definition, and source of variables -------------------

## Load the environmental covariates
load(file.path(enviro_dir, "EnvironmentalCovariates/s2met_environmental_covariates.RData"))

covariate_names <- bind_rows(
  crossing(distinct(filter(growth_stage_covariates, stage != "heading"), stage), 
           variable = str_subset(string = names(growth_stage_covariates)[-1:-3], pattern = "radn_mean", negate = TRUE)) %>% 
    unite(environmental_covariate, sep = ".") %>% mutate(type = "weather", source = "Daymet"),
  tibble(type = "soil", environmental_covariate = names(soil_covariates)[-1], source = "HWSD"))
  
## Save table
write_csv(x = variable_names, path = file.path(fig_dir, "table_s2_variable_names_definitions_draft1.csv"))




# Table S3: ECs selected per trait - concurrent and historical -----------------

selected_features <- bind_rows(
  mutate(concurrent_features, analysis = "concurrent"),
  mutate(historical_features, analysis = "historical")
) %>%
  mutate(features = map(features, "optVariables")) %>%
  filter(model %in% c("model3", "model5"), source == "daymet",
         str_detect(feat_sel_type, "apriori|rfa"), str_detect(feat_sel_type, "nosoil", negate = TRUE)) %>%
  select(-source, -model) %>%
  unnest(features) %>%
  filter(features != "line_name") %>%
  mutate(feature_type = ifelse(str_detect(features, "line_name"), "main", "gxe"),
         feature_type = fct_rev(feature_type),
         features = str_remove(features, "line_name:"),
         obs = "X")

# Output a table
selected_features %>%
  select(analysis, trait, CovariateSet = feat_sel_type, names(.)) %>%
  spread(feature_type, obs) %>%
  mutate(analysis = f_timeframe_replace(analysis), trait = str_add_space(trait),
         CovariateSet = f_ec_selection_replace(CovariateSet)) %>%
  rename_at(vars(matches("^[a-z]")), str_to_title) %>%
  rename(GxE = Gxe) %>%
  write_csv(x = ., path = file.path(fig_dir, "table_s3_ecs_per_trait.csv"), na = "")


# Table S4: within-trait reliability ------------------------------------------






 

























# Supplemental Table XX - full model phenotypic variance analysis ---------

# Load the full model analysis results
load(file.path(result_dir, "full_model_variance_analysis.RData"))

pheno_variance_analysis1 <- pheno_variance_analysis_out %>%
  # Make a note as to whether any interaction covariates were included
  mutate(any_int_cov = map_lgl(features, ~any(str_detect(., "line_name:"))) | feat_sel_type == "all") %>%
  unnest(results)

# Calculate proportions of total phenotypic variance
pheno_variance_analysis2 <- pheno_variance_analysis1 %>%
  group_by(trait, population, analysis, feat_sel_type) %>%
  mutate(total_variance = sum(variance)) %>%
  ungroup() %>%
  mutate(source_prop_total_variance = total_source_variance / total_variance,
         prop_total_variance  = variance / total_variance)
  


## Environment ##
# Clean up for a table
env_pheno_variance_analysis <- pheno_variance_analysis2 %>%
  filter(analysis == "environment") %>%
  mutate(feat_sel_type = ifelse(feat_sel_type == "rfa_cv", "rfa_cv_adhoc", feat_sel_type)) %>%
  select(trait, population, covariate_set = feat_sel_type, any_int_cov, source, term, 
         prop_source_variance = variance_prop, prop_total_variance, source_prop_total_variance) %>%
  mutate(trait = str_add_space(trait), population = f_pop_replace(population), covariate_set = f_ec_selection_replace(covariate_set),
         source = str_replace_all(source, c("geno" = "G", "site" = "E", "gxs" = "G x E", "units" = "Residuals")),
         term = str_replace_all(term, c("geno" = "G", "site" = "E", "gxs" = "G x E", "units" = "Residuals")),
         term = str_replace_all(term, c("GK" = "Markers", "G x EK" = "Markers x Covariates",  "EK" = "Covariates", "I" = "")))  %>%
  # convert to factors
  mutate_at(vars(term, source), fct_inorder)

## Plot overall contributions to variance
env_pheno_variance_analysis %>%
  filter(covariate_set == "All") %>%
  distinct(trait, population, source, source_prop_total_variance) %>%
  # Create annotation
  mutate_at(vars(contains("prop")), list(annotation = ~paste0(format_numbers(. * 100), "%"))) %>%
  # Plot
  ggplot(aes(x = population, y = source_prop_total_variance, fill = source)) +
  geom_col() +
  geom_text(aes(y = source_prop_total_variance, label = annotation), size = base_geom_text_size,
            position = position_stack(vjust = 0.5)) +
  scale_y_continuous(breaks = pretty, labels = scales::percent, name = "Percent variance explained") +
  facet_grid(~ trait) +
  theme_genetics(base_size = base_font_size)
  

## Plot variance explained by markers/covariates/etc.
env_pheno_variance_analysis %>%
  filter(term != "Residuals") %>%
  ggplot(aes(x = source, y = prop_source_variance, fill = term)) +
  geom_col() +
  facet_grid(trait ~ population + covariate_set) +
  theme_genetics(8)
    

# Create a table
env_pheno_variance_analysis_table <- env_pheno_variance_analysis %>%
  mutate_at(vars(contains("prop")), ~formatC(x = signif(., 2), digits = 2, width = 2, format = "fg")) %>%
  arrange(trait, population, covariate_set) %>%
  rename_all(~str_to_title(str_replace_all(., "_", " ")))


# Calculate difference in variance explained
environment_pheno_variance_analysis_diff <- env_pheno_variance_analysis_table %>%
  rename_all(make.names) %>%
  select(-Any.Int.Cov, -Prop.Total.Variance, -Source.Prop.Total.Variance) %>%
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


# Clean up for a table
loc_pheno_variance_analysis <- pheno_variance_analysis2 %>%
  filter(analysis == "location") %>%
  mutate(feat_sel_type = ifelse(feat_sel_type == "rfa_cv", "rfa_cv_adhoc", feat_sel_type)) %>%
  select(trait, population, covariate_set = feat_sel_type, any_int_cov, source, term, 
         prop_source_variance = variance_prop, prop_total_variance, source_prop_total_variance) %>%
  mutate(trait = str_add_space(trait), population = f_pop_replace(population), covariate_set = f_ec_selection_replace(covariate_set),
         source = str_replace_all(source, c("geno" = "G", "site" = "L", "gxs" = "G x L", "units" = "Residuals")),
         term = str_replace_all(term, c("geno" = "G", "site" = "L", "gxs" = "G x L", "units" = "Residuals")),
         term = str_replace_all(term, c("GK" = "Markers", "G x LK" = "Markers x Covariates",  "LK" = "Covariates", "I" = "")))  %>%
  # convert to factors
  mutate_at(vars(term, source), fct_inorder)

## Plot overall contributions to variance
loc_pheno_variance_analysis %>%
  filter(covariate_set == "All") %>%
  distinct(trait, population, source, source_prop_total_variance) %>%
  # Create annotation
  mutate_at(vars(contains("prop")), list(annotation = ~paste0(format_numbers(. * 100), "%"))) %>%
  # Plot
  ggplot(aes(x = population, y = source_prop_total_variance, fill = source)) +
  geom_col() +
  geom_text(aes(y = source_prop_total_variance, label = annotation), size = base_geom_text_size,
            position = position_stack(vjust = 0.5)) +
  scale_y_continuous(breaks = pretty, labels = scales::percent, name = "Percent variance explained") +
  facet_grid(~ trait) +
  theme_genetics(base_size = base_font_size)




## Plot variance explained by markers/covariates/etc.
loc_pheno_variance_analysis %>%
  filter(term != "Residuals") %>%
  ggplot(aes(x = source, y = prop_source_variance, fill = term)) +
  geom_col() +
  facet_grid(trait ~ population + covariate_set) +
  theme_genetics(8)


# Create a table
loc_pheno_variance_analysis_table <- loc_pheno_variance_analysis %>%
  mutate_at(vars(contains("prop")), ~formatC(x = signif(., 2), digits = 2, width = 2, format = "fg")) %>%
  arrange(trait, population, covariate_set) %>%
  rename_all(~str_to_title(str_replace_all(., "_", " ")))


# Calculate difference in variance explained
loc_pheno_variance_analysis_table %>%
  rename_all(make.names) %>%
  select(-Any.Int.Cov, -Prop.Total.Variance, -Source.Prop.Total.Variance) %>%
  filter(Population == "FP", str_detect(Term, "Covariates")) %>%
  {full_join(x = filter(., Covariate.Set == "StepwiseCV"), y = filter(., Covariate.Set != "StepwiseCV"),
             by = c("Trait", "Population", "Source", "Term"))} %>%
  mutate_at(vars(contains("Prop")), parse_number) %>%
  mutate(prop_diff = Prop.Source.Variance.x - Prop.Source.Variance.y)


# Output table
write_csv(x = location_pheno_variance_analysis_table, 
          path = file.path(fig_dir, "full_model_phenotypic_variance_analysis_location.csv"))


