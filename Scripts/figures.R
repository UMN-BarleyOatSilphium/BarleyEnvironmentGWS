## BarleyEnvironmentGWS
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
### Figures and Supplemental information are created in order of reference 
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
load(file.path(result_dir, "apsim_growth_model_results.RData"))
load(file.path(result_dir, "feature_selection_results.RData"))


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


concurrent_cgm_toplot <- concurrent_growth_staging_daymet %>%
  filter(location %in% cgm_locations_plot, year == cgm_year_plot)

historical_cgm_toplot <- location_historical_growth_staging_daymet %>%
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
# Z50 - Z69: heading
# Z70 - Z91: grain fill
# 

cgm_toplot_growth_stages <- cgm_toplot %>%
  mutate(growth_stages = map(growth_model, ~{
    
    df <- .
    
    ## Determine stages
    stages <- mutate(df, stage = case_when(
      between(zadok_stage, 10, 30) ~ "early_vegetative",
      between(zadok_stage, 30, 50) ~ "late_vegetative",
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
# Historical timeframe to highlight
time_frame_ex <- subset(historical_feature_selection, trait == "GrainYield" & feat_sel_type == "stepwise_cv_adhoc" & model == "model5",
                        time_frame, drop = TRUE)
# Extract the years
time_frame_ex <- seq_range(as.numeric(str_extract_all(string = time_frame_ex, pattern = "[0-9]{4}")[[1]]), by = 1)

# First make some edits and then summarize over years
cgm_toplot_growth_stage_weather_historical <- cgm_toplot_growth_stage_weather %>%
  filter(variable %in% covariates_plot, timeframe == "historical", year %in% time_frame_ex) %>%
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
  modify_at(., "historical", ~filter(., year %in% time_frame_ex)) %>%
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




# Figure 2: Outline of prediction approach --------------------------------









# Figure 3. LOEO prediction accuracy ---------------------------------------


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
  filter(model %in% str_subset(models_use, "3|5"), selection == "stepwise_cv_adhoc") %>%
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
  filter(model %in% str_subset(models_use, "id", negate = TRUE), selection %in% c("none", "stepwise_cv_adhoc")) %>%
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
ggsave(filename = "figure3_partA_draft.jpg", plot = g_loo_pred_obs, path = fig_dir,
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
ggsave(filename = "figure3_partB_draft.jpg", plot = g_accuracy_within_env1, path = fig_dir,
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
      mutate(annotation = paste0("r[MP]=='", ability_all, "'"))
    
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
ggsave(filename = "figure3_loeo_predictions.jpg", plot = g_figure2, path = fig_dir,
       width = figwidth_twocol, height = 12, dpi = dpi_use, units = "cm")

## Save as vector image
ggsave(filename = "figure3_loeo_predictions.svg", plot = g_figure2, path = fig_dir,
       width = figwidth_twocol, height = 12, dpi = dpi_use, units = "cm")





# Figure 4. LOLO predictions ---------------------------------


## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))

## Part A - accuracy across locations - LOLO

lolo_pred_obs_df <- predictions_df %>%
  filter(type == "lolo") %>%
  filter(model == "model3_cov", selection == "stepwise_cv_adhoc", time_frame_selection == "bestOverall") %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  # mutate(trait = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"))
  mutate(trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))


## Filter for relevant within-location cases 
accuracy_within_loc_toplot <- predictive_ability %>%
  ## Adjust model names
  mutate(model = str_replace_all(model, c("model2" = "model4", "model3" = "model5"))) %>%
  filter(type == "lolo",
         model %in% str_subset(models_use, "id", negate = TRUE), 
         selection %in% c("none", "stepwise_cv_adhoc"),
         (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), "\n(n=", nGroup+1, ")")) %>%
  unite(group, trait, model, pop, remove = FALSE) %>%
  # Add a test name for faceting
  mutate(type_facet_x = "'Within-location'~italic(r[MP])")



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
      facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.background = element_blank())
  })


# Combine the plots
g_lolo_pred_obs <- g_lolo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.05, 1))


# Save
ggsave(filename = "figure4_partA_draft.jpg", plot = g_lolo_pred_obs, path = fig_dir,
       height = 120, width = 80, units = "mm", dpi = dpi_use)


## Part B - accuracy within locations - LOLO

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
                    guide = guide_legend(label.position = "left", title.position = "left")) +
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title.position = "left")) +
  facet_grid(. ~ pop, switch = "y") +
  theme_genetics(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
        legend.text.align = 1, legend.position = c(0.95, 0.95), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line"),
        axis.title.x = element_blank(), strip.text.y = element_text(color = "white"),
        legend.justification = c(1, 0.5))

## Save
ggsave(filename = "figure4_partB_draft.jpg", plot = g_accuracy_within_loc, path = fig_dir,
       width = 80, height = 40, dpi = dpi_use, units = "mm")



## Plot both accuracy across environments and within environments per trait

# calculate the character width of the y axis
y_axis_width <- max(nchar(pretty(lolo_pred_obs_df$value))) 

# Split by trait and Plot
g_plotlist2 <- lolo_pred_obs_df %>%
  split(.$trait) %>%
  map(~{
    # Extract annotations
    ann_df <- distinct(.x, trait, pop, ability_all) %>%
      mutate(annotation = paste0("r[MP]=='", ability_all, "'"))
    
    # Predicted/observed phenotypes
    p1 <- ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty, labels = function(x) formatC(x, width = y_axis_width)) +
      scale_fill_manual(values = location_colors, guide = FALSE) +
      facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.background = element_blank(), axis.title = element_blank())
    
    within_env_df <- subset(accuracy_within_loc_toplot, trait == unique(.x$trait))

    # Add annotation for number of environments
    nE <- paste0("n = ", unique(within_env_df$nGroup) + 1)
    
    
    # Accuracy within environments
    p2 <- within_env_df %>%
      ggplot(aes(x = pop, y = ability, color = model, fill = model, group = group)) +
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

## Save
ggsave(filename = "figure4_lolo_predictions.jpg", plot = g_figure2, path = fig_dir,
       width = figwidth_twocol, height = 12, dpi = dpi_use, units = "cm")

## Save as vector image
ggsave(filename = "figure4_lolo_predictions.svg", plot = g_figure2, path = fig_dir,
       width = figwidth_twocol, height = 12, dpi = dpi_use, units = "cm")
















# Supplementary figures ---------------------------------------------------


# Figure S1: validation of predicted flowering time from crop model -------

# Load the results of the cultivar analysis
load(file.path(result_dir, "apsim_cultivar_predicted_flowering_results.RData"))

# Create an annotation df
ann_df <- cgm_pred_HD_summary %>%
  mutate(annotation = paste0("R^2==", format_numbers(cor^2), "*';'~RMSE==", format_numbers(rmse)),
         cultivar = fct_reorder(cultivar, rmse))


## Plot the results of each cultivar
g_cultivar_growth_modelA <- pred_obs_heading %>%
  mutate(cultivar = factor(cultivar, levels = levels(ann_df$cultivar))) %>%
  ggplot(aes(x = pred_HD, y = obs_HD)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 0.3) +
  facet_wrap(~ cultivar, ncol = 7) +
  scale_x_continuous(name = "Predicted heading date (days)", breaks = pretty) +
  scale_y_continuous(name = "Mean environment heading date (days)", breaks = pretty) +
  theme_genetics() +
  theme(panel.background = element_rect(fill = alpha("white", 0), color = "grey85"))
  
# Subset the cultivars used in the analysis
g_cultivar_growth_modelB <- g_cultivar_growth_modelA %>%
  modify_at("data", ~filter(., cultivar == "bass")) +
  geom_point(size = 1) +
  geom_text(data = subset(ann_df, cultivar == "bass"), aes(x = Inf, y = -Inf, label = annotation), size = base_geom_text_size, parse = TRUE,
            vjust = -0.5, hjust = 1)

## Merge plots
g_combined <- plot_grid(g_cultivar_growth_modelA, g_cultivar_growth_modelB, nrow = 1, rel_widths = c(1, 0.5),
                        labels = subfigure_labels[1:2], label_size = subfigure_label_size)

# Save
ggsave(filename = "figure_s1_growth_model_cultivar_comparison.jpg", plot = g_combined, path = fig_dir,
       width = 8, height = 4, dpi = 1000)






## Find the cultivar that results in the most accuracy prediction
(repr_cultivar <- cgm_pred_HD_summary %>%
    filter(stage == "flowering") %>%
    arrange(rmse, desc(cor)) %>%
    as.data.frame() %>%
    slice(1))





# Predict heading date from the crop growth model
predicted_heading_date_cultivars <- concurrent_growth_staging_cultivars_daymet %>% 
  unnest(apsim_out) %>% 
  filter(between(zadok_stage, 50, 70)) %>% 
  mutate(dap = as.numeric(Date - planting_date)) %>% 
  group_by(trial, location, year, environment, simulation_name) %>% 
  summarize(first_dap = first(dap), mean_dap = mean(dap))
  






# Figure S2: Comparison of environment-specific EC sets -------------------

## Number and overlap of covariates in the analyses ##


## Concurrent

## Combine concurrent feature selection df
concurrent_features <- bind_rows( concurrent_apriori_feature_selection, concurrent_stepwise_feature_selection, concurrent_all_features ) %>%
  rename(features = covariates) %>%
  mutate(feat_sel_type = str_replace(feat_sel_type, "rfa_cv", "stepwise_cv"))

## combine model2 and model3 covariates
concurrent_features1 <- concurrent_features %>%
  filter(str_detect(feat_sel_type, "nosoil", negate = TRUE)) %>%
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
  filter(source == "daymet", n > 0) %>%
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n))

concurrent_features2 <- concurrent_features1 %>%
  filter(source == "daymet")


## Modify the feature count plot to remove nasapower
g_concurrent_features_count2 <- concurrent_features_count1 %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  filter(source == "daymet") %>%
  ggplot(aes(x = feat_sel_type, y = n, fill = feature_class)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main environment covariate", "interaction_features" = "GxE covariates"), 
                      name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = c(0.50, 0.97), strip.placement = "outside",
        legend.justification = "left")




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
  scale_y_continuous(name = "Number of\noverlapping covariates", breaks = pretty) +
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
  scale_y_continuous(name = "Number of\noverlapping covariates", breaks = pretty) +
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
  group_by(covariates, source, feature_class) %>%
  summarize(n = n()) %>%
  mutate(nTotal = sum(n)) %>%
  ungroup() %>%
  mutate(covariates = fct_reorder(covariates, nTotal, .fun = max, .desc = TRUE))


## Plot
g_concurrent_indiv_feature_counts <- concurrent_indiv_feature_counts %>%
  mutate(feature_class = fct_inorder(feature_class)) %>%
  ggplot(aes(x = covariates, y = n, fill = feature_class)) +
  geom_col() +
  scale_fill_discrete(guide = FALSE) +
  scale_x_discrete(name = "Covariate") +
  scale_y_continuous(name = "Count") +
  theme_genetics(base_size = base_font_size) +
  theme(legend.position = c(0.75, 0.75), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Combine the plots
right_plot <- plot_grid(g_concurrent_features_overlap_trait, g_concurrent_features_overlap_featsel, g_concurrent_indiv_feature_counts, 
                        ncol = 1, rel_heights = c(0.5, 0.5, 1), align = "hv", axis = "lr", labels = subfigure_labels[-1], 
                        label_size = subfigure_label_size)  

merged_plot <- plot_grid(g_concurrent_features_count2, right_plot, nrow = 1, align = "v", axis = "tb",
                         labels = c(subfigure_labels[1], ""), label_size = subfigure_label_size)  

# Save
ggsave(filename = "figure_S2_concurrent_features_comparison_merged_daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = dpi_use)






# Figure S3: Comparison of location-specific EC sets ----------------------

# Edit historical_feature_selection
historical_feature_selection1 <- historical_feature_selection %>%
  filter(selection == "bestOverall", feat_sel_type == "stepwise_cv_adhoc")

historical_apriori_feature_selection1 <- historical_apriori_feature_selection %>%
  inner_join(., distinct(historical_feature_selection1, trait, time_frame))

historical_all_features1 <- historical_all_features %>%
  inner_join(., distinct(historical_feature_selection1, trait, time_frame))



## Combine historical feature selection df
historical_features <- bind_rows( historical_apriori_feature_selection1, historical_all_features1,  historical_feature_selection1 ) %>%
  rename(features = covariates) %>% select(-time_frame)
  
  

## combine model2 and model3 covariates
historical_features1 <- historical_features %>%
  filter(str_detect(feat_sel_type, "nosoil", negate = TRUE)) %>%
  spread(model, features) %>%
  mutate_at(vars(contains("model")), ~map(., "optVariables")) %>%
  mutate(features = map2(model4, model5, union)) %>%
  mutate(features = map(features, ~setdiff(., "line_name"))) %>%
  select(-contains("model")) %>%
  mutate(interaction_features = map(features, ~str_subset(., ":")),
         main_features = map2(features, interaction_features, setdiff),
         source = "daymet")


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
  scale_fill_discrete(labels = c("main_features" = "Main location covariate", "interaction_features" = "GxL covariates"), 
                      name = NULL) +  
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = c(0.50, 0.97), strip.placement = "outside",
        legend.justification = "left")


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
  scale_y_continuous(name = "Number of\noverlapping covariates", breaks = pretty) +
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
  scale_y_continuous(name = "Number of\noverlapping covariates", breaks = pretty) +
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
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Combine the plots
right_plot <- plot_grid(g_historical_features_overlap_trait, g_historical_features_overlap_featsel, g_historical_indiv_feature_counts, 
                        ncol = 1, rel_heights = c(0.5, 0.5, 1), align = "hv", axis = "lr", labels = subfigure_labels[-1], 
                        label_size = subfigure_label_size)  

merged_plot <- plot_grid(g_historical_features_count1, right_plot, nrow = 1, align = "v", axis = "tb",
                         labels = c(subfigure_labels[1], ""), label_size = subfigure_label_size)  

# Save
ggsave(filename = "figure_S3_historical_features_comparison_merged_daymet.jpg", plot = merged_plot,
       path = fig_dir, height = 8, width = 6, dpi = dpi_use)






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
  filter(type == "loeo", selection != "LASSO") %>%
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
  filter(type == "loeo", model %in% models_use, str_detect(selection, "nosoil", negate = TRUE),
         str_detect(selection, "lasso", negate = TRUE)) %>%
  # Add individual points
  mutate(trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         selection = fct_inorder(selection),
         model = str_replace(model, "_id", "_cov"),
         bias = bias * 100) %>%
  unite(group, trait, model, pop, selection, remove = FALSE) %>%
  select(-rmse, -varR, -heritability, -ability) %>%
  gather(statistic, value, accuracy, bias)


## Common plot modifier
g_mod <- list(
  scale_x_discrete(name = "Covariate set", labels = f_ec_selection_replace, guide = guide_axis(check.overlap = TRUE, n.dodge = 2)),
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)),
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title = NULL)),
  facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(pop = f_pop_replace), scales = "free_y"),
  theme_genetics(base_size = 6),
  theme(strip.background = element_blank(), strip.placement = "outside", axis.line.x = element_blank(),
        legend.text.align = 1, legend.position = c(0.89, 0.95), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line"))
)


# Plot point and line range for accuracy
g_accuracy_within_env <- accuracy_bias_within_env_toplot %>%
  filter(statistic == "accuracy") %>%
  ggplot(aes(x = selection, y = value, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = expression('Within-environment'~italic(r[MG])), breaks = pretty) +
  g_mod


# Plot point and line range for bias
g_bias_within_env <- accuracy_bias_within_env_toplot %>%
  filter(statistic == "bias") %>%
  ggplot(aes(x = selection, y = value, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = "Within-environment bias (%)", breaks = pretty) +
  g_mod +
  theme(legend.position = "none")

# Combine plots
g_accuracy_bias_within_env <- plot_grid(g_accuracy_within_env, g_bias_within_env, align = "hv",
                                        labels = subfigure_labels[1:2], label_size = subfigure_label_size)


## Save
ggsave(filename = "figure_S5_loeo_accuracy_bias_within_environments.jpg", plot = g_accuracy_bias_within_env, 
       path = fig_dir, width = 6, height = 5, dpi = 1000)




# Figure S6 - Historical weather data timeframe analysis -------------------


## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))


# Extract results
historical_timeframe_selection_out1 <- historical_timeframe_analysis %>%
  filter(feat_sel_type == "adhoc", method == "stepwise") %>%
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
  filter(model == "model5") %>%
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

# Separate df for the selected timeframe (to add as points and an annotation)
timeframe_used_df <- historical_feature_selection %>% 
  filter(feat_sel_type == "stepwise_cv_adhoc", model == "model5", selection == "bestRecent") %>%
  distinct(trait, time_frame) %>% 
  inner_join(., historical_timeframe_selection_out1)

# Plot together
g_timeframe_analysis <- historical_timeframe_selection_out1 %>%
  ggplot(aes(x = length, y = RMSE_scale, color = trait)) +
  geom_line(lwd = 0.25) +
  geom_point(data = timeframe_used_df, size = 0.3, shape = 15) +
  # # Add an annotation
  # annotate(geom = "curve", x = 10, y = 0.35, xend = min(timeframe_used_df$length) + 1, yend = min(timeframe_used_df$RMSE_scale),
  #          curvature = 0.25, arrow = arrow(angle = 30, length = unit(0.25, "line"), ends = "last"), lwd = 0.25) +
  # annotate(geom = "text", label = "Chosen timeframe", x = 10.5, y = 0.39, hjust = 1, size = 2) +
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Scaled RMSEP") +
  scale_x_continuous(name = "Years of weather data", trans = "reverse") +
  scale_color_manual(values = trait_colors, name = NULL, labels = function(x) abbreviate(str_add_space(x), 2),
                     guide = guide_legend(nrow = 2)) +
  theme_genetics(base_size = base_font_size) +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.5, 0.10),
        legend.direction = "horizontal", legend.background = element_rect(fill = alpha("white", 0)),
        legend.key.width = unit(0.25, "line"), legend.key.height = unit(0.5, "line"))


## Plot window results
##
historical_window_selection_out1_toplot <- historical_timeframe_selection_out1 %>%
  filter(time_frame_type == "window", model == "model5") %>%
  mutate(length = fct_inseq(as.factor(length))) %>%
  mutate_at(vars(contains("year")), ~ymd(paste0(., "0101")))

# Separate df for the selected timeframe (to add as points and an annotation)
timeframe_used_df <- historical_feature_selection %>% 
  filter(feat_sel_type == "stepwise_cv_adhoc", model == "model5", selection == "bestOverall", str_detect(time_frame, "window")) %>%
  distinct(trait, time_frame) %>% 
  inner_join(., historical_window_selection_out1_toplot)



# Plot together
g_window_analysis <- historical_window_selection_out1_toplot %>%
  ggplot(aes(x = end_year, y = RMSE_scale, color = trait, lty = length)) +
  geom_line(data = subset(historical_window_selection_out1_toplot, length == 5), lwd = 0.25, alpha = 1) +
  geom_line(data = subset(historical_window_selection_out1_toplot, length == 10), lwd = 0.25, alpha = 0.5) +
  geom_line(data = subset(historical_window_selection_out1_toplot, length == 15), lwd = 0.25, alpha = 0.5) +
  geom_point(data = timeframe_used_df, size = 0.3) +
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Scaled RMSEP") +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y", name = "End year of window") +
  scale_color_manual(values = trait_colors, guide = FALSE) +
  scale_linetype_manual(name = "Window length (yrs)", values = c("5" = 1, "10" = 2, "15" = 3), breaks = c(5, 10, 15),
                        guide = guide_legend(title.position = "top")) +
  theme_genetics(base_size = base_font_size) +
  theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.5, 0.10),
        legend.direction = "horizontal", legend.key.width = unit(0.5, "line"), legend.key.height = unit(0.5, "line"),
        legend.background = element_rect(fill = alpha("white", 0)))


# Combine plots
g_timeframe_window <- plot_grid(g_timeframe_analysis, g_window_analysis, nrow = 1,
                                labels = subfigure_labels[1:2], label_size = subfigure_label_size)

# Save
ggsave(filename = "figure_S6_timeframe_window_analysis.jpg", plot = g_timeframe_window, path = fig_dir,
       width = 3, height = 1.5,  dpi = 1000)







# Figure S7 - LOLO predictive ability across locations ----------

## Plot predicted vs observed phenotypic values for LOLO for all models, traits,
## EC selections, etc.

## Plot all
g_lolo_all_list <- loo_pred_obs_df %>%
  filter(type == "lolo", selection != "LASSO", (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
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
ggsave(filename = "figure_S7_lolo_across_site_prediction_accuracy.jpg", plot = g_lolo_all, path = fig_dir,
       height = 20, width = 22, units = "cm", dpi = 1000)



# Figure S8 - LOLO prediction accuracy and bias within locations ----------

## Filter for relevant cases
accuracy_bias_within_loc_toplot <- predictive_ability %>%
  filter(type == "lolo", model %in% models_use, str_detect(selection, "nosoil", negate = TRUE),
         str_detect(selection, "lasso", negate = TRUE), (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
  # Add individual points
  mutate(selection = str_replace(selection, "rfa_adhoc", "rfa_cv_adhoc"),
         trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         selection = fct_inorder(selection),
         model = str_replace(model, "_id", "_cov"),
         bias = bias * 100) %>%
  unite(group, trait, model, pop, selection, remove = FALSE) %>%
  select(-contains("all"), -rmse) %>%
  gather(statistic, value, ability, bias)


# Plot point and line range for accuracy
g_accuracy_within_loc <- accuracy_bias_within_loc_toplot %>%
  filter(statistic == "ability") %>%
  ggplot(aes(x = selection, y = value, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = expression('Within-location'~italic(r[MP])), breaks = pretty) +
  g_mod


# Plot point and line range for bias
g_bias_within_loc <- accuracy_bias_within_loc_toplot %>%
  filter(statistic == "bias") %>%
  ggplot(aes(x = selection, y = value, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = "Within-location bias (%)", breaks = pretty) +
  g_mod +
  theme(legend.position = "none")

# Combine plots
g_accuracy_bias_within_loc <- plot_grid(g_accuracy_within_loc + theme(legend.position = c(0.10, 0.85)), 
                                        g_bias_within_loc, align = "hv",
                                        labels = subfigure_labels[1:2], label_size = subfigure_label_size)
  


## Save
ggsave(filename = "figure_S8_lolo_accuracy_bias_within_locations.jpg", plot = g_accuracy_bias_within_loc, 
       path = fig_dir, width = 6, height = 5, dpi = 1000)


# Figure S9 - external predictive ability across environments ----------

## Plot all
g_ext_env_all_list <- loo_pred_obs_df %>%
  filter(type == "env_external", selection != "LASSO") %>%
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
             y = case_when(trait == "GrainYield" & str_detect(selection, "None|StepwiseCV", negate = TRUE) ~ Inf,
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
ggsave(filename = "figure_S9_external_environment_across_site_prediction_accuracy.jpg", plot = g_ext_env_all, path = fig_dir,
       height = 16, width = 22, units = "cm", dpi = 1000)





# Figure S10 - external predictive ability across locations ----------

## Plot all
g_ext_loc_all_list <- loo_pred_obs_df %>%
  filter(type == "loc_external", selection != "LASSO", (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
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
ggsave(filename = "figure_S10_external_location_across_site_prediction_accuracy.jpg", plot = g_ext_loc_all, path = fig_dir,
       height = 16, width = 22, units = "cm", dpi = 1000)




# Figure S11 - Comparison of accuracy when using different timeframes ----

# Plot LOLO results when using longer-term historical data
loc_longterm_pred_obs_df <- loo_pred_obs_df %>%
  # Extract data from the GxE model and both the recent and best timeframes
  filter(type %in% c("lolo", "loc_external"), selection == "StepwiseCV", model == "g + e + (ge)") %>%
  mutate(time_frame = str_extract(time_frame, "[0-9]{4}_[0-9]{4}") %>% str_replace(., "_", "-"),
         time_frame_facet = paste0(time_frame_selection, "Timeframe: ", time_frame),
         group = paste0(f_validation_replace(pop), ", ", f_type_replace(type)),
         group = fct_reorder2(group, type, pop, .desc = FALSE) %>% fct_relevel(., "Untested offspring, Holdout location", after = Inf))

g_loc_longterm_list <- loc_longterm_pred_obs_df %>%
  split(.$trait) %>%
  map(~{

    # Extract annotations
    ann_df <- distinct(.x, trait1, group, model, selection, ability_all, time_frame_facet) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))

    left_join(complete(distinct(.x, trait1, group, time_frame_facet), trait1, group, time_frame_facet), .x) %>%
      ggplot(aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.2) +
      geom_point(size = 0.1, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = location_colors, guide = FALSE) +
      facet_grid(trait1 ~ group + time_frame_facet, switch = "y", scales = "free_x", labeller = labeller(trait1 = label_parsed)) +
      theme_genetics(base_size = 6) +
      theme(strip.placement = "outside", strip.background = element_blank())

  })

# Combine the plots
g_loc_longterm_all <- plot_grid(plotlist = g_loc_longterm_list, ncol = 1, align = "hv", axis = "lr",
                                labels = subfigure_labels[seq_along(g_loc_longterm_list)], label_size = subfigure_label_size)


# Save
ggsave(filename = "figure_S11_location_longterm_ec_predictions.jpg", plot = g_loc_longterm_all, path = fig_dir,
       height = 18, width = 18, units = "cm", dpi = 1000)





# Figure S12 - environmental relationships ----------------------


# Calculate environmental means
# 
# Fit models to calculate environmental means
# 
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
  mutate_at(vars(ends_with("_features")),
            list(ec_mat = ~modify_if(., sapply(., length) > 0, ~Env_mat(x = ec_tomodel_scaled_mat, terms = .x, method = "Jarq")))) %>%
  mutate_at(vars(ends_with("_ec_mat")), 
            ~modify_if(., sapply(., inherits, "character"), ~tcrossprod(ec_tomodel_scaled_mat) - tcrossprod(ec_tomodel_scaled_mat))) %>%
  # Subset environments
  left_join(., nest(group_by(distinct(S2_MET_BLUEs, trait, environment), trait), .key = "environment")) %>%
  mutate_at(vars(ends_with("_ec_mat")), ~map2(., environment, ~.x[.y$environment, .y$environment]))


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
    
  })) %>% ungroup() %>%
  # Calculate correlations and annotate
  mutate(cor_test = map(corMat_envMean_df, ~cor.test(x = .x$correlation, y = .x$env_mean_dist)),
         cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
         annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))

## Plot
g_env_mean_cor_concurrent_features <- env_mean_cor_concurrent_features %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(corMat_envMean_df) %>%
  mutate(env_mean_dist = ifelse(trait == "GrainYield", env_mean_dist / 1000, env_mean_dist)) %>%
  ggplot(aes(x = correlation, y = env_mean_dist)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5, color = "blue") +
  geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y", 
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Covariate-estimated\nenvironmental correlation", breaks = pretty) +
  scale_y_continuous(name = "Difference in environmental mean", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")

# Save
ggsave(filename = "concurrent_ec_correlation_vs_env_mean_dist.jpg", plot = g_env_mean_cor_concurrent_features,
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
    
  })) %>% ungroup() %>%
  # Calculate correlations and annotate
  mutate(cor_test = map(phenoCor_corMat, ~cor.test(x = .x$correlation, y = .x$phenoCor)),
         cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
         annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))

## Plot
g_phenCor_concurrent_features1 <- phenCor_concurrent_features1 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(phenoCor_corMat) %>%
  ggplot(aes(x = correlation, y = phenoCor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Covariate-estimated\nenvironmental correlation", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation between environments", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")

# Save
ggsave(filename = "concurrent_ec_correlation_vs_phenoCor.jpg", plot = g_phenCor_concurrent_features1,
       path = fig_dir, height = 5, width = 5, dpi = 1000)


## For each trait and covariate type (interaction/main), plot a heatmap and
## a dendrogram
concurrent_features_heatmap_plots <- concurrent_features_env_cormat %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  group_by(feat_sel_type) %>%
  mutate(n = seq(n())) %>%
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
    segment_data <- segment(clust_data)
    
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
      mutate(annotation = ifelse(is.na(correlation), annotation, NA)) %>%
      # Add continuous x/y coordinates
      left_join(., clust_lab_data, by = c("environment" = "label")) %>%
      left_join(., clust_lab_data, by = c("environment2" = "label")) %>%
      select(x_center = x.x, y_center = x.y, correlation)
    
    # Generate common axis limit
    env_axis_limits <- with(clust_lab_data, c(min(x - 0.5), max(x + 0.5))) + 0.1 * c(-1, 1)
    
    
    # Plot
    g_heat <- dat2 %>%
      ggplot(aes(x = x_center, y = y_center)) +
      # geom_tile(aes(fill = correlation), color = ifelse(is.na(dat2$correlation), "black", NA), lwd = 0.1) +
      geom_tile(aes(fill = correlation), lwd = 0.1) +
      # geom_text(aes(label = annotation), size = 1) +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3], na.value = "white",
                           midpoint = 0, limits = c(-1.01, 1.01), breaks = pretty,
                           name = "Estimated environmental\ncorrelation", guide = guide_colorbar(title.position = "top")) +
      scale_x_continuous(breaks = clust_lab_data$x, labels = clust_lab_data$label, limits = env_axis_limits, expand = c(0, 0)) +
      scale_y_continuous(breaks = clust_lab_data$x, labels = clust_lab_data$label, expand = c(0, 0), position = "right") +
      # facet_grid(row$trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
      # labs(subtitle = paste(str_add_space(row$trait), f_ec_selection_replace(row$feat_sel_type), sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = c(0.25, 0.75),
            axis.text.y = element_blank(), legend.direction = "horizontal")
    
    # Remove legend if n != 1
    if (row$n != 1) {
      g_heat <- g_heat + theme(legend.position = "none")
    }
    
    ## Plot dendrogram
    g_dendro <- ggplot(data = segment_data) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_x_continuous(breaks = clust_lab_data$x, labels = clust_lab_data$label, limits = env_axis_limits, expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_color_manual(name = NULL, values = c("black", "red")) +
      theme_genetics(base_size = 8) +
      theme(legend.key.height = unit(5, "pt"), axis.title = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(),
            axis.ticks.x = element_blank())
    
    # Combine plots
    plot_grid(g_heat, g_dendro, align = "h", rel_widths = c(1, 0.8))
    
  }) %>% ungroup()



## Combine plots for use in the supplemental figure

# First edit each plot
g_env_mean_cor_concurrent_features1 <- g_env_mean_cor_concurrent_features %>%
  modify_at("data", ~filter(.x, feat_sel_type == "stepwise_cv_adhoc")) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y", labeller = labeller(trait = str_add_space))

g_phenCor_concurrent_features2 <- g_phenCor_concurrent_features1 %>%
  modify_at("data", ~filter(.x, feat_sel_type == "stepwise_cv_adhoc")) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y", labeller = labeller(trait = str_add_space))

# Merge heatmaps
g_heat_dendro_select <- subset(concurrent_features_heatmap_plots, feat_sel_type == "stepwise_cv_adhoc", plot, drop = TRUE) %>%
  modify_at(-1, ~. + theme(legend.position = "none")) %>%
  plot_grid(plotlist = ., ncol = 1)

# Merge all plots
g_heatmap_cor_compare <- plot_grid(g_env_mean_cor_concurrent_features1, 
                                   g_phenCor_concurrent_features2 + theme(strip.text = element_blank()),
                                   plot_spacer() + theme(panel.background = element_rect(fill = alpha("white", 0))), 
                                   g_heat_dendro_select,
                                   nrow = 1, align = "h", axis = "tblr", rel_widths = c(0.5, 0.5, 0.05, 1),
                                   labels = subfigure_labels[c(1,2,NA,3)], label_size = subfigure_label_size)


# Save
ggsave(filename = "figure_s12_environmental_covariate_heatmap.jpg", plot = g_heatmap_cor_compare,
       path = fig_dir, width = 10, height = 12, dpi = 1000)






# Figure S13 - location relationships ----------------------


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
loc_tomodel_scaled_mat_list <- c(historical_ec_tomodel_timeframe_scaled, historical_ec_tomodel_window_scaled) %>%
  subset(., grepl(pattern = "daymet", x = names(.))) %>%
  setNames(object = ., nm = gsub(pattern = "daymet\\.", replacement = "", x = names(.))) %>%
  map(~select(.x, -source, -time_frame) %>%
        filter(location %in% c(train_test_loc, validation_loc)) %>%
        as.data.frame() %>%
        column_to_rownames("location") %>%
        as.matrix() )

## Calculate environmental relationship matrices based on these features
historical_features_loc_relmat <- historical_features1 %>%
  left_join(., subset(historical_feature_selection, selection == "bestOverall" & feat_sel_type == "stepwise_cv_adhoc") %>%
              distinct(trait, time_frame)) %>%
  mutate(interaction_features = map(interaction_features, ~str_remove(., "line_name:"))) %>%
  mutate_at(vars(ends_with("_features")), ~modify_if(., sapply(., length) == 0, ~c("1"))) %>%
  mutate_at(vars(ends_with("_features")),
            list(ec_mat = ~map2(., time_frame, ~Env_mat(x = loc_tomodel_scaled_mat_list[[.y]], terms = .x, method = "Jarq")))) %>%
  # Subset environments
  left_join(., nest(group_by(distinct(S2_MET_BLUEs, trait, location), trait), .key = "location")) %>%
  mutate_at(vars(ends_with("_ec_mat")), ~map2(., location, ~.x[row.names(.x) %in% .y$location, colnames(.x) %in% .y$location]))
  

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
env_mean_concurrent_features <- historical_features_loc_cormat %>%
  left_join(., loc_mean_dist) %>%
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
    
  })) %>% ungroup() %>%
  # Calculate correlations and annotate
  mutate(cor_test = map(corMat_envMean_df, ~cor.test(x = .x$correlation, y = .x$env_mean_dist)),
         cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
         annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))

## Plot
g_env_mean_cor_concurrent_features <- env_mean_cor_concurrent_features %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(corMat_envMean_df) %>%
  mutate(env_mean_dist = ifelse(trait == "GrainYield", env_mean_dist / 1000, env_mean_dist)) %>%
  ggplot(aes(x = correlation, y = env_mean_dist)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5, color = "blue") +
  geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y", 
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Covariate-estimated\nlocation correlation", breaks = pretty) +
  scale_y_continuous(name = "Difference in location mean", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")


## Plot the phenotypic correlation versus the correlation due to covariates
phenCor_concurrent_features <- historical_features_loc_cormat %>%
  left_join(., loc_phenoCor) %>%
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
    
  })) %>% ungroup() %>%
  # Calculate correlations and annotate
  mutate(cor_test = map(phenoCor_corMat, ~cor.test(x = .x$correlation, y = .x$phenoCor)),
         cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
         annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))

## Plot
g_phenCor_concurrent_features1 <- phenCor_concurrent_features1 %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  unnest(phenoCor_corMat) %>%
  ggplot(aes(x = correlation, y = phenoCor)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
  facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
             labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
  scale_x_continuous(name = "Covariate-estimated\nlocation correlation", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation between locations", breaks = pretty) +
  theme_genetics(base_size = 8) +
  theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")

# Location abbreviations
lop_abbr <- trial_info %>%
  distinct(location) %>%
  mutate(loc = str_replace_all(location, "_", " ") %>% str_to_title() %>% 
           str_remove_all(" ") %>% abbreviate() )

## For each trait and covariate type (interaction/main), plot a heatmap and
## a dendrogram
features_heatmap_plots <- historical_features_loc_cormat %>%
  filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
  group_by(feat_sel_type) %>%
  mutate(n = seq(n())) %>%
  group_by(trait, feat_sel_type) %>%
  do(plot = {
    row <- .
    
    corMat <- row$corMat[[1]]
    
    # First perform clustering on the relationship matrix
    env_clust <- hclust(dist(corMat, method = "euclidean"), method = "average")
    # get the data
    clust_data <- dendro_data(model = env_clust)
    # Label data
    clust_lab_data <- mutate(clust_data$labels, group = ifelse(label %in% train_test_loc, "Train/test", "Holdout"),
                             group = fct_relevel(group, "Holdout", after = Inf),
                             location = label) %>%
      left_join(., lop_abbr, by = c("location")) %>%
      mutate(location = label, label = loc, loc = location)
    segment_data <- segment(clust_data)
    
    # Factor order of environments
    env_order <- levels(clust_lab_data$loc)
    
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
      mutate(annotation = ifelse(is.na(correlation), annotation, NA)) %>%
      # Add continuous x/y coordinates
      left_join(., clust_lab_data, by = c("environment" = "loc")) %>%
      left_join(., clust_lab_data, by = c("environment2" = "loc")) %>%
      select(x_center = x.x, y_center = x.y, correlation)
    
    # Generate common axis limit
    env_axis_limits <- with(clust_lab_data, c(min(x - 0.5), max(x + 0.5))) + 0.1 * c(-1, 1)
    
    
    # Plot
    g_heat <- dat2 %>%
      ggplot(aes(x = x_center, y = y_center)) +
      # geom_tile(aes(fill = correlation), color = ifelse(is.na(dat2$correlation), "black", NA), lwd = 0.1) +
      geom_tile(aes(fill = correlation), lwd = 0.1) +
      # geom_text(aes(label = annotation), size = 1) +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3], na.value = "white",
                           midpoint = 0, limits = c(-1.01, 1.01), breaks = pretty,
                           name = "Estimated location\ncorrelation", guide = guide_colorbar(title.position = "top")) +
      scale_x_continuous(breaks = clust_lab_data$x, labels = clust_lab_data$label, limits = env_axis_limits, expand = c(0, 0)) +
      scale_y_continuous(breaks = clust_lab_data$x, labels = clust_lab_data$label, expand = c(0, 0), position = "right") +
      # facet_grid(row$trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
      # labs(subtitle = paste(str_add_space(row$trait), f_ec_selection_replace(row$feat_sel_type), sep = ", ")) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.25, "line"),
            strip.placement = "outside", axis.title = element_blank(), legend.position = c(0.25, 0.75),
            axis.text.y = element_blank(), legend.direction = "horizontal")
    
    # Remove legend if n != 1
    if (row$n != 1) {
      g_heat <- g_heat + theme(legend.position = "none")
    }
    
    ## Plot dendrogram
    g_dendro <- ggplot(data = segment_data) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_x_continuous(breaks = clust_lab_data$x, labels = clust_lab_data$label, limits = env_axis_limits, expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_color_manual(name = NULL, values = c("black", "red")) +
      theme_genetics(base_size = 8) +
      theme(legend.key.height = unit(5, "pt"), axis.title = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(),
            axis.ticks.x = element_blank())
    
    # Combine plots
    plot_grid(g_heat, g_dendro, align = "h", rel_widths = c(1, 0.8))
    
  }) %>% ungroup()



## Combine plots for use in the supplemental figure

# First edit each plot
g_mean_cor_features1 <- g_env_mean_cor_concurrent_features %>%
  modify_at("data", ~filter(.x, feat_sel_type == "stepwise_cv_adhoc")) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y", labeller = labeller(trait = str_add_space))

g_phenCor_features2 <- g_phenCor_concurrent_features1 %>%
  modify_at("data", ~filter(.x, feat_sel_type == "stepwise_cv_adhoc")) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y", labeller = labeller(trait = str_add_space))

# Merge heatmaps
g_heat_dendro_select <- subset(features_heatmap_plots, feat_sel_type == "stepwise_cv_adhoc", plot, drop = TRUE) %>%
  modify_at(-1, ~. + theme(legend.position = "none")) %>%
  plot_grid(plotlist = ., ncol = 1)

# Merge all plots
g_heatmap_cor_compare <- plot_grid(g_mean_cor_features1, 
                                   g_phenCor_features2 + theme(strip.text = element_blank()),
                                   plot_spacer() + theme(panel.background = element_rect(fill = alpha("white", 0))), 
                                   g_heat_dendro_select,
                                   nrow = 1, align = "h", axis = "tblr", rel_widths = c(0.5, 0.5, 0.05, 1),
                                   labels = subfigure_labels[c(1,2,NA,3)], label_size = subfigure_label_size)


# Save
ggsave(filename = "figure_s13_location_covariate_heatmap.jpg", plot = g_heatmap_cor_compare,
       path = fig_dir, width = 10, height = 12, dpi = 1000)






# # Figure S14 - 75-25 environment cross-validation ---------------
# 
# # Filter across-environment accuracy for LOEO and CV
# loo_cv_compare <- predictive_ability %>%
#   filter(type %in% c("cv_env", "cv_loc", "loeo", "lolo"), str_detect(selection, "soil", negate = TRUE)) %>%
#   distinct_at(vars(trait, type, cv_rep, pop, model, selection, contains("all"))) %>%
#   # Modify the names
#   mutate(target = ifelse(type %in% c("loeo", "cv_env"), "Environment", "Location"),
#          type = ifelse(type %in% c("loeo", "lolo"), "Leave-one-out", "75%-25% train-test"), 
#          selection = ifelse(selection == "rfa_adhoc", "rfa_cv_adhoc", selection)) %>%
#   # Summarize over reps
#   group_by(trait, type, target, pop, model, selection) %>%
#   mutate(n = n()) %>%
#   summarize_at(vars(ability_all, rmse_all, n), list(~mean, ~sd)) %>%
#   ungroup() %>%
#   mutate_at(vars(ends_with("sd")), list(se = ~. / sqrt(n_mean)))
# # 
# # 
# # ## Compare 75-25 cross-validation results to LOO results
# # g_cv_loeo_compare_acc_list <- loo_cv_compare %>%
# #   # Filter for relevant models and variable selections
# #   filter(selection %in% c("rfa_cv_adhoc", "none"), model %in% c("model1", "model2_cov", "model3_cov")) %>%
# #   split(.$target) %>%
# #   map(~{
# #     ggplot(data = .x, aes(x = pop, y = ability_all_mean, color = type)) +
# #       geom_point(position = position_dodge(0.9), size = 0.75) +
# #       geom_errorbar(aes(ymin = ability_all_mean - ability_all_sd_se, ymax = ability_all_mean + ability_all_sd_se), 
# #                     position = position_dodge(0.9), width = 0.5, lwd = 0.25) + 
# #       scale_color_discrete(name = NULL) +
# #       scale_y_continuous(name = "Predictive ability", breaks = pretty) +
# #       scale_x_discrete(name = NULL, labels = function(x) str_replace(f_validation_replace(x), " ", "\n")) +
# #       facet_grid(trait ~ model + selection, switch = "y", scales = "free_y",
# #                  labeller = labeller(model = f_model_replace, selection = f_ec_selection_replace)) +
# #       facet_grid(trait ~ model, switch = "y", labeller = labeller(model = f_model_replace, trait = str_add_space)) +
# #       theme_genetics(base_size = 8) +
# #       theme(legend.position = "top", panel.border = element_rect(colour = "grey85", fill = alpha("white", 0)))
# #   })
# # 
# # 
# # # Combine
# # g_cv_loeo_compare_acc <- plot_grid(plotlist = map(g_cv_loeo_compare_acc_list, ~.+theme(legend.position = "none")), ncol = 1,
# #                                labels = subfigure_labels, label_size = subfigure_label_size)
# # g_cv_loeo_compare_acc1 <- plot_grid(get_legend(g_cv_loeo_compare_acc_list[[1]]), g_cv_loeo_compare_acc, ncol = 1, rel_heights = c(0.05, 1))
# # 
# # ggsave(filename = "figure_S13_compare_loo_and_cv.jpg", plot = g_cv_loeo_compare_acc1,
# #        path = fig_dir, height = 7, width = 4, dpi = 1000)
# # 
# 
# g_cv_loeo_compare_rmse_list <- loo_cv_compare %>%
#   # Filter for relevant models and variable selections
#   filter(selection %in% c("rfa_cv_adhoc", "none"), model %in% c("model1", "model2_cov", "model3_cov")) %>%
#   split(.$target) %>%
#   map(~{
#     ggplot(data = .x, aes(x = pop, y = rmse_all_mean, fill = type)) +
#       geom_col(position = position_dodge(0.9)) +
#       geom_errorbar(aes(ymin = rmse_all_mean - rmse_all_sd_se, ymax = rmse_all_mean + rmse_all_sd_se), 
#                     position = position_dodge(0.9), width = 0.5, lwd = 0.25) + 
#       scale_fill_paletteer_d(package = "wesanderson", palette = "Royal1", name = NULL) +
#       scale_y_continuous(name = "RMSE", breaks = pretty) +
#       scale_x_discrete(name = NULL, labels = function(x) str_replace(f_validation_replace(x), " ", "\n")) +
#       facet_grid(trait ~ model, switch = "y", scales = "free_y", labeller = labeller(model = f_model_replace, trait = str_add_space)) +
#       labs(subtitle = paste0("Prediction target: ", unique(.x$target), "s")) +
#       theme_genetics(base_size = 8) +
#       theme(legend.position = "top", panel.border = element_rect(colour = "grey85", fill = alpha("white", 0)))
#   })
# 
# 
# # Combine
# g_cv_loeo_compare_rmse <- plot_grid(plotlist = map(g_cv_loeo_compare_rmse_list, ~.+theme(legend.position = "none")), nrow = 1,
#                                    labels = subfigure_labels, label_size = subfigure_label_size)
# g_cv_loeo_compare_rmse1 <- plot_grid(get_legend(g_cv_loeo_compare_rmse_list[[1]]), g_cv_loeo_compare_rmse, ncol = 1, rel_heights = c(0.05, 1))
# 
# ggsave(filename = "figure_S14_compare_loo_and_cv.jpg", plot = g_cv_loeo_compare_rmse1,
#        path = fig_dir, height = 4, width = 8, dpi = 1000)














# Supplementary tables ----------------------------------------




# Table S1: Variance component breakdown ----------------------------------

# Load the initial variance component partition
load(file.path(result_dir, "stage_two_phenotypic_analysis.RData"))

## Create a table for the partition of environment-specific variance components
ge_var_comp <- stage_two_varcomp_env %>% 
  unnest(var_comp) %>%
  # Calculate proportions of variance
  group_by(trait, population) %>%
  mutate(prop_var = (variance / sum(variance)) * 100) %>%
  ungroup() %>%
  mutate_if(is.numeric, format_numbers) %>% 
  # mutate(annotation = paste0(variance, " (", se, ") / ", prop_var, "%")) %>% 
  mutate(annotation = paste0(variance, " (", prop_var, "%)")) %>% 
  select(trait, population, term, annotation) %>% 
  spread(trait, annotation) %>%
  mutate(analysis = "environment")


## Create a table for the partition of location-specific variance components
gl_var_comp <- stage_two_varcomp_loc %>% 
  unnest(var_comp) %>%
  # Calculate proportions of variance
  group_by(trait, population) %>%
  mutate(prop_var = (variance / sum(variance)) * 100) %>%
  ungroup() %>%
  mutate_if(is.numeric, format_numbers) %>% 
  # mutate(annotation = paste0(variance, " (", se, ") / ", prop_var, "%")) %>% 
  mutate(annotation = paste0(variance, " (", prop_var, "%)")) %>% 
  select(trait, population, term, annotation) %>% 
  spread(trait, annotation) %>%
  mutate(analysis = "location")

## Combine amd save
var_comp_print <- bind_rows(ge_var_comp, gl_var_comp) %>%
  mutate(population = f_pop_replace(population),
         term = str_replace_all(term, c("environment" = "Environment", "location" = "Location", "gxe" = "GxE", 
                                        "line_name" = "Genotype", "units" = "Residual")),
         term = ifelse(analysis == "location" & term == "GxE", "GxL", term),
         term = fct_relevel(term, "Genotype", "Location"),
         analysis = str_to_title(analysis)) %>%
  arrange(analysis, population, term) %>%
  rename_at(vars(matches("^[A-Z]", ignore.case = FALSE)), str_add_space) %>%
  select(analysis, Population = population, Term = term, names(.))

# Save
write_csv(x = var_comp_print, path = file.path(fig_dir, "table_s1_variance_partition.csv"))




# Table S2: Name, definition, and source of variables -------------------

## Load the environmental covariates
load(file.path(result_dir, "environmental_covariates.RData"))

variable_names <- bind_rows(
  crossing(distinct(filter(concurrent_growth_stage_covariates), stage), 
           variable = str_subset(string = names(concurrent_growth_stage_covariates)[-1:-3], pattern = "radn_mean", negate = TRUE)) %>% 
    unite(environmental_covariate, sep = ".") %>% mutate(type = "weather", source = "Daymet"),
  tibble(type = "soil", environmental_covariate = names(soil_covariates)[-1:-2], source = "HWSD")) %>%
  separate(environmental_covariate, c("stage", "variable"), sep = "\\.", fill = "left") %>%
  distinct(type, source, variable) %>%
  arrange(desc(type), variable)
  
## Save table
write_csv(x = variable_names, path = file.path(fig_dir, "table_s2_variable_names_definitions_draft1.csv"))




# Table S3: ECs selected per trait - concurrent and historical -----------------

selected_features <- bind_rows(
  mutate(concurrent_features, analysis = "concurrent"),
  mutate(historical_features, analysis = "historical")
) %>%
  mutate(covariate = map(features, "optVariables")) %>%
  filter(model %in% c("model3", "model5"),
         str_detect(feat_sel_type, "apriori|stepwise"), str_detect(feat_sel_type, "nosoil", negate = TRUE)) %>%
  select(-source, -model, -selection, -features) %>%
  unnest(covariate) %>%
  filter(covariate != "line_name") %>%
  mutate(feature_type = ifelse(str_detect(covariate, "line_name"), "gxe", "main"),
         feature_type = fct_rev(feature_type),
         covariate = str_remove(covariate, "line_name:"),
         obs = "X")

# Output a table
selected_features %>%
  select(analysis, trait, names(.)) %>%
  spread(feature_type, obs) %>%
  mutate(analysis = f_timeframe_replace(analysis), trait = str_add_space(trait),
         feat_sel_type = f_ec_selection_replace(feat_sel_type)) %>%
  rename_at(vars(matches("^[a-z]")), str_to_title) %>%
  rename(Interaction = Gxe, CovariateSet = Feat_sel_type) %>%
  write_csv(x = ., path = file.path(fig_dir, "table_s3_ecs_per_trait.csv"), na = "")


# Table S4: within-trait reliability ------------------------------------------

# Add location and year information to the heritability estimates;
# clean-up and save as a table
env_herit_table_plot <- env_trait_herit %>% 
  select(-varR) %>% 
  left_join(distinct(trial_info, environment, location, year)) %>% 
  select(environment, location, year, trait, heritability) %>%
  mutate(heritability = format_numbers(x = heritability)) %>%
  rename_all(str_to_title) %>%
  spread(Trait, Heritability) %>%
  arrange(Environment)


write_csv(x = env_herit_table_plot, path = file.path(fig_dir, "table_s4_env_trait_herit.csv"), na = "")





## Other notes to reference in the manuscript

# Number of covariates per selection type and trait
concurrent_features %>% 
  unnest(features) %>% 
  filter(model == "model2", str_detect(feat_sel_type, "soil", negate = TRUE)) %>% 
  mutate(nFeatures = sapply(features, length)) %>%
  group_by(trait, feat_sel_type) %>% 
  summarize_at(vars(nFeatures), list(~min, ~max)) %>%
  arrange(feat_sel_type)


# Prediction accuracy within environments
accuracy_within_env_toplot %>%
  group_by(trait, model, selection, type, pop) %>%
  summarize_at(vars(ability, accuracy, bias, rmse), mean) %>%
  ungroup() %>%
  arrange(type, model, selection, pop, trait) %>%
  as.data.frame()

accuracy_within_loc_toplot %>%
  group_by(trait, model, selection, type, pop) %>%
  summarize_at(vars(ability, bias, rmse), mean) %>%
  ungroup() %>%
  arrange(type, model, selection, pop, trait) %>%
  as.data.frame()


# External environment predictions
predictive_ability %>%
  filter(type == "env_external", selection == "stepwise_cv_adhoc", (is.na(time_frame_selection) | time_frame_selection == "bestOverall"),
         model == "model3_cov") %>%
  distinct(trait, site, model, pop, ability, bias, rmse) %>% 
  left_join(., env_trait_herit, by = c("trait", "site" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability)) %>%
  group_by(trait, pop, model) %>%
  summarize_at(vars(ability, accuracy, bias, rmse), mean) %>%
  ungroup() %>%
  arrange(model, pop, trait)

predictive_ability %>%
  filter(type == "loc_external", selection == "stepwise_cv_adhoc", (is.na(time_frame_selection) | time_frame_selection == "bestOverall"),
         model == "model3_cov") %>%
  distinct(trait, site, model, pop, ability, bias, rmse) %>% 
  group_by(trait, pop, model) %>%
  summarize_at(vars(ability, bias, rmse), mean) %>%
  ungroup() %>%
  arrange(model, pop, trait)


# Number of trials per trait
S2_MET_BLUEs %>%
  group_by(trait) %>% 
  summarize(nEnv = n_distinct(environment))
