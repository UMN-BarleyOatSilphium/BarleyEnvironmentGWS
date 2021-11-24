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
base_font_size <- 8
# Base geom text size
base_geom_text_size <- base_font_size * (5/14)

# 1 column figure width (inch)
figwidth_onecol <- 3.5
figwidth_twocol <- 7.0


## Create location-based color palettes
location_colors <- paletteer_d("ggsci::default_igv", n = n_distinct(trial_info$location) + 2)
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
         parent = case_when(parent == "TRUE" ~ "parent", parent == "FALSE" ~ "nonparent", TRUE ~ "offspring"),
         parent = fct_relevel(parent, "parent"))

# Color scheme for populations
pop_colors <- setNames(object = c(neyhart_palette("umn1", n = 4), neyhart_palette("umn2")[c(5, 8)]), nm = levels(K_pco_df$program))

# Function to rename the "parent" column
f_parent_rename <- function(x) c("parent" = "FP (parent)", "nonparent" = "FP", "offspring" = "OP")[x]

# Create the visualization
g_popstr <- K_pco_df %>% 
  ggplot(aes(x = PCo1, y = PCo2)) + 
  # Parents
  geom_point(data = subset(K_pco_df, parent == "parent"), aes(fill = program, shape = parent), size = 1) +
  # Non-parents and offspring
  geom_point(data = subset(K_pco_df, parent != "parent"), aes(color = program, shape = parent), size = 1) +
  scale_x_continuous(name = K_pco_varprop[1], breaks = pretty) +
  scale_y_continuous(name = K_pco_varprop[2], breaks = pretty, position = "right") +
  scale_color_manual(values = pop_colors, drop = FALSE, guide = FALSE) +
  scale_fill_manual(values = pop_colors, drop = FALSE, guide = FALSE) +
  scale_shape_manual(name = NULL, values = c(parent = 21, nonparent = 16, offspring = 17), 
                     labels = f_parent_rename, drop = FALSE) +
  theme_genetics(base_size = base_font_size) +
  # theme(legend.position = "none", plot.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)),
  #       panel.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)))
  theme(legend.key.height = unit(0.5, "line"), legend.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)),
        # legend.position = c(1,1), legend.direction = "vertical", legend.justification = c(1,1),
        legend.position = "top", legend.direction = "vertical", legend.justification = "right",
        plot.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0)),
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
            size = floor(base_geom_text_size), color = ifelse(use_loc_info_toplot1$location == "Arlington", "black", "white")) +
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

## Combine plots
layout <- c(
  area(t = 1, l = 1, b = 10, r = 8),
  area(t = 2.5, l = 7, b = 10, r = 10)
)

g_pop_map <- plot_grid(g_map_v2, labels = subfigure_labels[1], label_size = subfigure_label_size) + 
  plot_grid(g_popstr, labels = subfigure_labels[2], label_size = subfigure_label_size, label_x = 0.9)+ 
  plot_layout(design = layout)

# Save this
ggsave(filename = "figure1_map_and_popstr.jpg", plot = g_pop_map, path = fig_dir,
       height = 3, width = 6, dpi = dpi_use)








# Figure 2: Using crop growth model to generate environmental covariates --------------------------------


## Plot examples of crop model outputs

# Load concurrent and historical crop model output
load(file.path(result_dir, "apsim_growth_model_results.RData"))
## Load covariate data
load(file.path(result_dir, "concurrent_historical_covariables.RData"))
load(file.path(result_dir, "feature_selection_results.RData"))


# Assign units to variables
covariate_rename_abbr <- c("mint" = "T[min]", "maxt" = "T[max]", "tmean" = "T[mean]", "gdd" = "GDD",
                           "water_balance" = "H[2]O Bal.", "radn" = "Radn.")

covariate_rename <- c("mint" = "Min. temp.", "maxt" = "Max. temp.", "tmean" = "Mean temp.",
                      "gdd" = "Growing deg. days", "water_balance" = "Water bal.", "radn" = "Solar radiation")

# Units
covariate_variable_unit <- setNames(c(rep("degree*C", 4), "mm", "MJ~m^-2"), names(covariate_rename_abbr))

## Function for renaming a covariate
f_covariate_replace2 <- function(x) paste0("atop('", covariate_rename[x], "',(", covariate_variable_unit[x], "))")
f_covariate_replace <- function(x) paste0("'", covariate_rename[x], "'~(", covariate_variable_unit[x], ")")
# Function for renaming timeframe
f_timeframe_replace <- function(x) c("concurrent" = "In-season", "historical" = "Historical location average")[x]


## Plot Bozeman and Columbus (Wooster)
cgm_locations_plot <- c("Bozeman", "Wooster")
cgm_year_plot <- 2017
# Assign colors to location
cgm_locations_color <- setNames(location_colors[1:2], cgm_locations_plot)

## Replacements and colors for growth stages

growth_stage_color <- setNames(object = paletteer_d("ggsci::default_igv")[c(28, 41, 1, 47, 34)],
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
  ## Combine data, growth model output, and growth stages
  mutate(data_use = pmap(select(., data, growth_model, growth_stages), ~{
    left_join(x = inner_join(..1, select(..2, date = Date, yday = day, tt = TT), by = c("yday", "date")), 
              y = ..3, by = c("date", "yday" = "day")) })) %>%
  unnest(data_use, names_repair = tidyr_legacy) %>%
  ## Calculate water stress as the difference between pet and rain
  mutate(water_balance = rain - pet) %>%
  rename(gdd = tt) %>%
  ## Gather
  gather(variable, value, rain, swe, maxt, mint, vp, radn, daylength, pet, elevation, gdd, water_balance) %>%
  # Skip NA stages
  filter(!is.na(stage))

## Create bars for growth stages
cgm_toplot_growth_stage_times <- cgm_toplot_growth_stage_weather %>% 
  group_by(location, timeframe, year, stage) %>% 
  filter(!is.na(stage)) %>% 
  summarize_at(vars(dap), list(min = min, max = max)) %>%
  ungroup() %>%
  mutate(y = as.numeric(as.factor(location)), dap = min,
         type = "Growth\nstage") # Dummy variable for growth stage facet




## Plot 2-3 covariates
covariates_plot <- c("mint", "water_balance")
covariates_plot1 <- paste0(covariates_plot, c("_mean", "_sum"))


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
  facet_grid(variable ~ timeframe, scales = "free_y", switch = "y", 
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
             switch = "y") +
  scale_color_manual(values = growth_stage_color, guide = FALSE) +
  scale_y_continuous(limits = c(-0.5, 2.5)) +
  scale_x_continuous(breaks = pretty, name = "Days after planting date", position = "bottom") +
  theme_genetics(base_size = base_font_size) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.background = element_blank(),
        axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.length.x = unit(base_font_size / 4, "pt"),
        strip.text.x = element_text(size = 6))


# Calculate averages of each variable within growth stages, plot below
concurrent_ec_plot <- ec_tomodel_centered[[1]] %>%
  gather(covariate, deviation, -source, -environment) %>%
  left_join(ec_tomodel_centers$daymet) %>%
  mutate(value = center + deviation) %>%
  separate(covariate, c("stage", "variable"), sep = "\\.", fill = "left") %>%
  mutate(variable = str_remove(variable, "_sum|_mean")) %>%
  filter(variable %in% covariates_plot) %>%
  inner_join(., select(concurrent_cgm_toplot, location, environment)) %>%
  # Add timeframes
  left_join(subset(cgm_toplot_growth_stage_yearly, timeframe == "concurrent"), .)

# Do the same for historical ECs
historical_ec_plot <- historical_ec_tomodel_timeframe_centered$daymet.time_frame17_1998_2014 %>%
  gather(covariate, deviation, -source, -time_frame, -location) %>%
  left_join(historical_ec_tomodel_timeframe_centers$daymet.time_frame17_1998_2014) %>%
  mutate(value = center + deviation) %>%
  separate(covariate, c("stage", "variable"), sep = "\\.", fill = "left") %>%
  mutate(variable = str_remove(variable, "_sum|_mean")) %>%
  filter(variable %in% covariates_plot) %>%
  inner_join(., select(concurrent_cgm_toplot, location)) %>%
  # Add timeframes
  left_join(subset(cgm_toplot_growth_stage_times_summary, timeframe == "historical"), .)

# Combine
ec_stages_plot <- bind_rows(concurrent_ec_plot, historical_ec_plot) %>%
  mutate(variable = f_covariate_replace(variable))



g_growth_stage_ecs <- ec_stages_plot %>%
  ggplot(aes(x = dap, color = stage)) +
  # Segments for variables in growth stages
  geom_segment(aes(x = min, xend = max, y = value, yend = value), lwd = 4, alpha = 0.15) + 
  # Segments for location colors
  # Lines for the entire growing season per location
  geom_segment(data = ec_stages_plot, aes(x = min, xend = max, y = value, yend = value), 
                                              color = cgm_locations_color[ec_stages_plot$location], 
                                              lwd = 0.5, inherit.aes = FALSE) +
  facet_grid(variable ~ timeframe, labeller = labeller(variable = label_parsed, timeframe = f_timeframe_replace), 
             switch = "y", scales = "free_y") +
  scale_color_manual(values = growth_stage_color, guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = NULL) +
  scale_x_continuous(breaks = pretty, name = "Days after planting date", position = "bottom") +
  theme_genetics(base_size = base_font_size) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        axis.title.y = element_blank())



## Combine weather and growth stage plot
# g_weather_stages <- plot_grid(
#   g_weather_concurrent_historical + theme(legend.position = "none", axis.title.x = element_blank(),
#                                           axis.text.x = element_blank()), 
#   g_growth_stages + theme(strip.text.x = element_blank()),
#   ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.25), labels = subfigure_labels[2:3], label_size = 8)

## Combine weather and growth stage plot
g_weather_stages <- plot_grid(
  g_weather_concurrent_historical + theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank()), 
  g_growth_stages + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                          strip.text.x = element_blank()),
  g_growth_stage_ecs + theme(strip.text.x = element_blank()),
  ncol = 1, align = "hv", axis = "lr", rel_heights = c(1, 0.35, 0.75), labels = subfigure_labels[1:3], label_size = subfigure_label_size,
  label_y = c(1, 1.03, 1.03), label_x = 0.01)

# Save
ggsave(filename = "figure2_gcm_to_ecs.jpg", plot = g_weather_stages, path = fig_dir, 
       width = figwidth_onecol, height = 5, dpi = dpi_use)







# Figure 3. LOEO accuracy within environments ---------------------------------------


# Load the compiled prediction results
load(file.path(result_dir, "prediction_accuracy_compiled.RData"))


# Color scheme for models
model_colors <- c(neyhart_palette("umn1")[1], rep(c(rep(neyhart_palette("umn1")[3], 2), rep(neyhart_palette("umn1")[4], 2)), 2))
names(model_colors) <- names(model_replace)

# Vector of models to report
models_use <- str_subset(string = names(model_colors), pattern = "cov1", negate = TRUE)


## Filter for relevant cases
accuracy_bias_within_env_toplot <- within_environment_prediction_accuracy %>%
  filter(type == "loeo", model %in% models_use, str_detect(selection, "nosoil", negate = TRUE),
         str_detect(selection, "lasso", negate = TRUE)) %>%
  # Add individual points
  mutate(selection = fct_inorder(selection),
         # trait1 = paste0(str_add_space(trait), " (n = ", nGroup+1, ")"),
         trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         trait2 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"),
         model = str_replace(model, "_id", "_cov"),
         pop = paste0("Prediction target: ", pop),
         bias = bias * 100,
         # Convert RMSE for GY to Mg ha
         rmse = ifelse(trait == "GrainYield", rmse / 1000, rmse)) %>%
  unite(group, trait, model, pop, selection, remove = FALSE) %>%
  select(-varR, -heritability, -ability) %>%
  gather(statistic, value, accuracy, bias, rmse)


## Common plot modifier
g_mod <- list(
  scale_x_discrete(name = "Covariate set", labels = f_ec_selection_replace,
                   guide = guide_axis(check.overlap = TRUE, n.dodge = 2)),
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)),
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title = NULL)),
  facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(pop = f_pop_replace), scales = "free_y"),
  theme_genetics(base_size = base_font_size),
  theme(strip.placement = "outside", axis.line.x = element_blank(),
        # panel.border = element_rect(color = "grey85", fill = alpha("white", 0)),
        panel.spacing.y = unit(0.5, "line"),
        panel.grid.major.y = element_line(color = "grey85"),
        legend.text.align = 1, legend.position = c(0.89, 0.97), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(color = "grey85", fill = alpha("white", 1)),
        legend.margin = margin(0, 2, 2, 2),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line"))
)


# Plot boxplot for accuracy
g_accuracy_within_env <- accuracy_bias_within_env_toplot %>%
  filter(statistic == "accuracy") %>%
  ggplot(aes(x = selection, y = value, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = expression('Within-environment'~italic(r[MG])), breaks = pretty) +
  g_mod



## Save
ggsave(filename = "figure3_loeo_accuracy_within_environments.jpg", plot = g_accuracy_within_env, 
       path = fig_dir, width = figwidth_onecol, height = 5, dpi = dpi_use)



# Figure 4. LOEO accuracy across environments ---------------------------------------


## plot predicted/observed phenotypes for each trait and by target population
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



# # Split by trait and Plot
# g_loo_pred_obs_list <- loo_pred_obs_df %>%
#   split(.$trait) %>%
#   map(~{
#     # Extract annotations
#     ann_df <- distinct(.x, trait, pop, ability_all) %>%
#       mutate(annotation = paste0("r[MP]==", ability_all))
#     
#     ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
#       geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
#       geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
#       geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
#                 parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
#       scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
#       scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
#       scale_fill_manual(values = environment_colors, guide = FALSE) +
#       facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
#       theme_genetics(base_size = base_font_size) +
#       theme(strip.placement = "outside", strip.background = element_blank())
#   })
# 
# 
# # Combine the plots
# g_loo_pred_obs <- g_loo_pred_obs_list %>%
#   modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
#   map(~. + theme(axis.title = element_blank())) %>%
#   plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
#   add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
#   plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
#                      rel_widths = c(0.05, 1))


# Split by trait and Plot
g_loo_pred_obs_list <- loo_pred_obs_df %>%
  mutate(pop = paste0("Prediction target: ", pop)) %>%
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
      facet_grid(pop~ trait1, switch = "y", 
                 labeller = labeller(trait1 = label_parsed, pop = f_pop_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside")
  })


# Combine the plots
g_loo_pred_obs <- g_loo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.y = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., nrow = 1, align = "h", rel_widths = c(1, rep(0.90, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = base_font_size, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.03, 1))



## Save
ggsave(filename = "figure4_loeo_predictions_across_environment.jpg", plot = g_loo_pred_obs, path = fig_dir,
       width = figwidth_twocol, height = 3, dpi = dpi_use)



# Figure 5. LOLO accuracy within locations ---------------------------------------

## Filter for relevant within-location cases 
accuracy_within_loc_toplot <- predictive_ability %>%
  ## Adjust model names
  mutate(model = str_replace_all(model, c("model2" = "model4", "model3" = "model5"))) %>%
  filter(type == "lolo",
         model %in% models_use, 
         str_detect(selection, "nosoil", negate = TRUE),
         str_detect(selection, "lasso", negate = TRUE),
         (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
  # Add individual points
  mutate(selection = fct_inorder(selection),
         # trait1 = paste0(str_add_space(trait), " (n = ", nGroup+1, ")"),
         trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         model = str_replace(model, "_id", "_cov"),
         pop = paste0("Prediction target: ", pop)) %>%
  unite(group, trait, model, pop, selection, remove = FALSE) %>%
  # Add a test name for faceting
  mutate(type_facet_x = "'Within-location'~italic(r[MP])")

  

## Common plot modifier
g_mod <- list(
  scale_x_discrete(name = "Covariate set", labels = f_ec_selection_replace,
                   guide = guide_axis(check.overlap = TRUE, n.dodge = 2)),
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace,
                    guide = guide_legend(label.position = "left", title = NULL)),
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace,
                     guide = guide_legend(label.position = "left", title = NULL)),
  facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(pop = f_pop_replace), scales = "free_y"),
  theme_genetics(base_size = base_font_size),
  theme(strip.placement = "outside", axis.line.x = element_blank(),
        # panel.border = element_rect(color = "grey85", fill = alpha("white", 0)),
        panel.spacing.y = unit(0.5, "line"),
        panel.grid.major.y = element_line(color = "grey85"),
        legend.text.align = 1, legend.position = c(0.10, 0.86), legend.key.size = unit(0.4, "line"),
        legend.title.align = 1, legend.background = element_rect(color = "grey85", fill = alpha("white", 1)),
        legend.margin = margin(0, 2, 2, 2),
        legend.text = element_text(size = 5), legend.key.height = unit(0.3, "line"))
)


# Plot boxplot for accuracy
g_accuracy_within_loc <- accuracy_within_loc_toplot %>%
  ggplot(aes(x = selection, y = ability, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = expression('Within-location'~italic(r[MP])), breaks = pretty) +
  g_mod



## Save
ggsave(filename = "figure5_lolo_accuracy_within_locations.jpg", plot = g_accuracy_within_loc, 
       path = fig_dir, width = figwidth_onecol, height = 5, dpi = dpi_use)



# Figure 6. LOLO accuracy across locations ---------------------------------------

loo_pred_obs_df <- predictions_df %>%
  filter(type == "lolo") %>%
  filter(model %in% str_subset(models_use, "3|5"), selection == "stepwise_cv_adhoc",
         (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  mutate(trait2 = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"),
         trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))



# # Split by trait and Plot
# g_loo_pred_obs_list <- loo_pred_obs_df %>%
#   split(.$trait) %>%
#   map(~{
#     # Extract annotations
#     ann_df <- distinct(.x, trait, pop, ability_all) %>%
#       mutate(annotation = paste0("r[MP]==", ability_all))
#     
#     ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
#       geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
#       geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
#       geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
#                 parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
#       scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
#       scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
#       scale_fill_manual(values = environment_colors, guide = FALSE) +
#       facet_grid(trait1 ~ pop, switch = "y", labeller = labeller(trait1 = label_parsed, pop = f_validation_replace)) +
#       theme_genetics(base_size = base_font_size) +
#       theme(strip.placement = "outside", strip.background = element_blank())
#   })
# 
# 
# # Combine the plots
# g_loo_pred_obs <- g_loo_pred_obs_list %>%
#   modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
#   map(~. + theme(axis.title = element_blank())) %>%
#   plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
#   add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
#   plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
#                      rel_widths = c(0.05, 1))


# Split by trait and Plot
g_loo_pred_obs_list <- loo_pred_obs_df %>%
  mutate(pop = paste0("Prediction target: ", pop)) %>%
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
      facet_grid(pop~ trait1, switch = "y", 
                 labeller = labeller(trait1 = label_parsed, pop = f_pop_replace)) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside")
  })


# Combine the plots
g_loo_pred_obs <- g_loo_pred_obs_list %>%
  modify_at(-1, ~. + theme(strip.text.y = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., nrow = 1, align = "h", rel_widths = c(1, rep(0.90, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = base_font_size, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.03, 1))



## Save
ggsave(filename = "figure6_lolo_predictions_across_locations.jpg", plot = g_loo_pred_obs, path = fig_dir,
       width = figwidth_twocol, height = 3, dpi = dpi_use)













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
  unnest(apsim_out, names_repair = tidyr_legacy) %>% 
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
  # geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Covariates for main\nenvironmental effect", 
                                 "interaction_features" = "Covariates for GxE effect"), 
                      name = NULL) +
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = c(0.40, 0.97), strip.placement = "outside",
        legend.justification = "left")




## Overlap between non-all feature selection types
concurrent_features_overlap <- concurrent_features2 %>%
  select(-features) %>%
  filter(! feat_sel_type %in% c("all")) %>%
  full_join(., ., by = c("trait", "source")) %>%
  inner_join(rename_all(as_tibble(t(combn(x = unique(concurrent_features1$feat_sel_type), m = 2))), 
                        ~paste0("feat_sel_type", c(".x", ".y"))), .) %>%
  mutate(interaction_features_overlap = map2(interaction_features.x, interaction_features.y, intersect),
         main_features_overlap = map2(main_features.x, main_features.y, intersect)) %>%
  select(contains("source"), contains("trait"), contains("feat_sel_type"), contains("overlap"))

# Examine the overlapping features
set_names(concurrent_features_overlap$main_features_overlap, concurrent_features_overlap$trait)

concurrent_features_overlap_featsel <- concurrent_features_overlap %>%
  mutate_at(vars(contains("features")), ~map_dbl(., n_distinct)) %>%
  gather(feature_class, n, contains("features")) %>%
  mutate(feature_class = str_remove(feature_class, "_overlap"))



## Plot
g_concurrent_features_overlap_featsel <- concurrent_features_overlap_featsel %>%
  mutate(feature_class = fct_rev(feature_class)) %>%
  mutate_at(vars(contains("feat_sel")), ~f_ec_selection_replace(., parse = FALSE)) %>%
  unite(feat_sel_type_pair, feat_sel_type.x, feat_sel_type.y, sep = " : ") %>%
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait, y = n, fill = feature_class)) +
  geom_col() +
  # geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(~ feat_sel_type_pair, switch = "y", labeller = labeller(feat_sel_type_pair = label_parsed)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = "Number of\noverlapping covariates") +
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
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n),
         feat_sel_type = f_ec_selection_replace(feat_sel_type, parse = FALSE)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  # geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(feat_sel_type ~ ., switch = "y", labeller = labeller(feat_sel_type = label_parsed)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Number of\noverlapping covariates") +
  theme_genetics(base_font_size) +
  theme(legend.position = "none", strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1))



## Count the number of times a particular covariate is selected

concurrent_features3 <- concurrent_features2 %>%
  select(-features) %>%
  gather(feature_class, covariates, contains("features")) %>%
  unnest(cols = c(covariates)) %>%
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
  scale_y_continuous(name = "Number of traits for which covariate was selected") +
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
  # geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(trait ~ ., switch = "y", labeller = labeller(trait = str_add_space)) +
  scale_fill_discrete(labels = c("main_features" = "Main location covariate", "interaction_features" = "GxL covariates"), 
                      name = NULL) +  
  scale_x_discrete(labels = f_ec_selection_replace, name = "Covariate set") +
  scale_y_continuous(name = "Number of covariates", breaks = pretty) +
  theme_genetics(base_font_size) +
  theme(legend.position = c(0.40, 0.97), strip.placement = "outside",
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
  mutate_at(vars(contains("feat_sel")), ~f_ec_selection_replace(., parse = FALSE)) %>%
  unite(feat_sel_type_pair, feat_sel_type.x, feat_sel_type.y, sep = " : ") %>%
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait, y = n, fill = feature_class)) +
  geom_col() +
  # geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(~ feat_sel_type_pair, switch = "y", labeller = labeller(feat_sel_type_pair = label_parsed)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait", labels = str_add_space) +
  scale_y_continuous(name = "Number of\noverlapping covariates") +
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
  mutate(label = paste0(toupper(abbreviate(feature_class, 1)), ": ", n),
         feat_sel_type = f_ec_selection_replace(feat_sel_type, parse = FALSE)) %>%
  filter(n > 0) %>%
  ggplot(aes(x = trait_pair, y = n, fill = feature_class)) +
  geom_col() +
  # geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = base_geom_text_size) +
  facet_grid(feat_sel_type ~ ., switch = "y", labeller = labeller(feat_sel_type = label_parsed)) +
  scale_fill_discrete(labels = c("main_features" = "Main", "interaction_features" = "Interaction"), name = NULL) +
  scale_x_discrete(name = "Trait pair") +
  scale_y_continuous(name = "Number of\noverlapping covariates") +
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
  scale_y_continuous(name = "Number of traits for which covariate was selected") +
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
  # left_join(., select(predictive_ability, trait:ability)) %>%
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
         selection = f_ec_selection_replace(selection, parse = FALSE),
         selection = ifelse(dummy == first(levels(dummy)), paste0("Covariate~Set: ", selection), selection)) %>%
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
      facet_grid(pop ~ model + selection, switch = "y", scales = "free_x", 
                 labeller = labeller(pop = f_pop_replace, selection = label_parsed)) +
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



















# Figure S5 - LOEO RMSEP within environments ----------



## Common plot modifier
g_mod <- list(
  scale_x_discrete(name = "Covariate set", labels = f_ec_selection_replace),
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace),
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace),
  facet_grid(trait2 ~ pop, switch = "y", labeller = labeller(pop = f_pop_replace, trait2 = label_parsed), 
             scales = "free_y"),
  theme_genetics(base_size = base_font_size),
  theme(strip.placement = "outside", axis.line.x = element_blank(),
        # panel.border = element_rect(color = "grey85", fill = alpha("white", 0)),
        panel.spacing.y = unit(0.5, "line"), panel.grid.major.y = element_line(color = "grey85"),
        legend.position = c(0.89, 0.97),  legend.title.align = 1, 
        legend.margin = margin(0, 2, 2, 2))
)



# Plot point and line range for bias
g_rmsep_within_env <- accuracy_bias_within_env_toplot %>%
  filter(statistic == "rmse") %>%
  mutate(pop = paste0("Prediction target: ", pop)) %>%
  ggplot(aes(x = selection, y = value, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = "Within-environment RMSEP", breaks = pretty) +
  g_mod +
  theme(legend.position = "top")



## Save
ggsave(filename = "figure_S5_loeo_bias_within_environments.jpg", plot = g_rmsep_within_env, 
       path = fig_dir, width = figwidth_onecol + 1, height = 5, dpi = dpi_use)







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
  geom_point(data = timeframe_used_df, size = 0.75, shape = 15) +
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
  geom_point(data = timeframe_used_df, size = 0.75) +
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
       width = 5, height = 2,  dpi = 1000)




# Figure S7 - LOLO predictive ability across locations ----------

## Plot predicted vs observed phenotypic values for LOLO for all models, traits,
## EC selections, etc.

# For predictions that used the variable selection procedure, only choose the correct
# time_frame
loo_pred_obs_df_stepwise <- historical_feature_selection %>%
  filter(str_detect(feat_sel_type, "nosoil", negate = TRUE), selection == "bestOverall") %>%
  select(trait, time_frame, selection = feat_sel_type) %>%
  mutate(selection = f_ec_selection_replace(selection, parse = FALSE)) %>%
  distinct() %>%
  inner_join(loo_pred_obs_df, .)

# Just the model without ECs
loo_pred_obs_df_noECs <- loo_pred_obs_df %>%
  filter(str_detect(selection, "None"))

# All other predictions should have time_frame == NA
loo_pred_obs_df_other <- loo_pred_obs_df %>%
  filter(selection != "LASSO", str_detect(selection, "stepwise", negate = TRUE)) %>%
  inner_join(., distinct(loo_pred_obs_df_stepwise, trait, time_frame))

loo_pred_obs_df1 <- bind_rows(loo_pred_obs_df_stepwise, loo_pred_obs_df_other, loo_pred_obs_df_noECs) %>%
  filter(type == "lolo") %>%
  arrange(trait, model, selection, pop)


## Plot all
g_lolo_all_list <- loo_pred_obs_df1 %>%
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
      facet_grid(pop ~ model + selection, switch = "y", scales = "free_x", 
                 labeller = labeller(pop = f_pop_replace, selection = label_parsed)) +
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



## Example figure for presentations
## Show base model vs full model for cross-validation and parent-offspring

## Part A - accuracy across locations - LOLO

lolo_pred_obs_df1 <- predictions_df %>%
  filter(type == "lolo", trait %in% c("GrainProtein", "GrainYield", "HeadingDate")) %>%
  filter( (model == "model3_cov" & selection == "stepwise_cv_adhoc" & time_frame_selection == "bestOverall") |
            (model == "model1" & selection == "none") ) %>%
  # Add annotations
  left_join(., across_site_prediction_accuracy_annotation) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Add trait units
  # mutate(trait = paste0("atop('", str_add_space(trait), "',(", trait_units1[trait], "))"))
  mutate(trait1 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"))

# Split by trait and pop
g_lolo_pred_obs_list <- lolo_pred_obs_df1 %>%
  group_by(trait, pop) %>%
  do(plot = {
    .x <- .
    # Extract annotations
    ann_df <- distinct(.x, model, trait, pop, ability_all) %>%
      mutate(annotation = paste0("r[MP]==", ability_all))
    
    ggplot(.x, aes(x = pred_complete, y = value, fill = test_group)) +
      geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
      geom_point(size = 0.15, shape = 21, color = alpha("white", 0)) +
      geom_text(data = ann_df, aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE,
                parse = TRUE, vjust = -1, hjust = 1.2, size = 1.5) +
      scale_x_continuous(name = "Predicted phenotypic value", breaks = pretty) +
      scale_y_continuous(name = "Observed phenotypic value", breaks = pretty) +
      scale_fill_manual(values = location_colors, guide = FALSE) +
      facet_grid(trait1 ~ model, switch = "y", 
                 labeller = labeller(trait1 = label_parsed, model = function(x) 
                   c("model1" = "Genomic markers", "model3_cov" = "Genomic markers +\nEnvironmental data")[x])) +
      theme_genetics(base_size = base_font_size) +
      theme(strip.placement = "outside", strip.text.y = element_blank(), strip.background = element_blank())
    
  }) %>% ungroup()


# Combine cross-validation plots
g_lolo_pred_obs_tp <- g_lolo_pred_obs_list %>%
  filter(pop == "tp") %>%
  pull(plot) %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.05, 1))

# Save
ggsave(filename = "lolo_prediction_example_tp.jpg", plot = g_lolo_pred_obs_tp, path = fig_dir,
       height = 3, width = 2.5, dpi = dpi_use)

# Combine parent-offspring validation plots
g_lolo_pred_obs_vp <- g_lolo_pred_obs_list %>%
  filter(pop == "vp") %>%
  pull(plot) %>%
  modify_at(-1, ~. + theme(strip.text.x = element_blank())) %>%
  map(~. + theme(axis.title = element_blank())) %>%
  plot_grid(plotlist = ., ncol = 1, align = "v", rel_heights = c(1, rep(0.85, length(.) - 1))) %>%
  add_sub(plot = ., label = "Predicted phenotypic value", size = 6, vjust = -1) %>%
  plot_grid(textGrob(label = "Observed phenotypic value", rot = 90, gp = gpar(fontsize = base_font_size)), .,
            rel_widths = c(0.05, 1))

# Save
ggsave(filename = "lolo_prediction_example_vp.jpg", plot = g_lolo_pred_obs_vp, path = fig_dir,
       height = 3, width = 2.5, dpi = dpi_use)











# Figure S8 - LOLO RMSEP within locations ----------

## Filter for relevant within-location cases 
accuracy_within_loc_toplot <- predictive_ability %>%
  ## Adjust model names
  mutate(model = str_replace_all(model, c("model2" = "model4", "model3" = "model5"))) %>%
  filter(type == "lolo",
         model %in% models_use, 
         str_detect(selection, "nosoil", negate = TRUE),
         str_detect(selection, "lasso", negate = TRUE),
         (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
  # Add individual points
  mutate(selection = fct_inorder(selection),
         # trait1 = paste0(str_add_space(trait), " (n = ", nGroup+1, ")"),
         trait1 = paste0(abbreviate(str_add_space(trait), 2), " (n = ", nGroup+1, ")"),
         model = str_replace(model, "_id", "_cov"),
         trait2 = paste0("'", abbreviate(str_add_space(trait), 2), "'~(", trait_units1[trait], ")"),
         pop = paste0("Prediction target: ", pop),
         bias = bias * 100,
         # Convert RMSE for GY to Mg ha
         rmse = ifelse(trait == "GrainYield", rmse / 1000, rmse)) %>%
  unite(group, trait, model, pop, selection, remove = FALSE) %>%
  # Add a test name for faceting
  mutate(type_facet_x = "'Within-location'~italic(r[MP])")



## Common plot modifier
g_mod <- list(
  scale_x_discrete(name = "Covariate set", labels = f_ec_selection_replace),
  scale_fill_manual(name = "Model", values = model_colors, labels = f_model_replace),
  scale_color_manual(name = "Model", values = model_colors, labels = f_model_replace),
  facet_grid(trait2 ~ pop, switch = "y", labeller = labeller(pop = f_pop_replace, trait2 = label_parsed), 
             scales = "free_y"),
  theme_genetics(base_size = base_font_size),
  theme(strip.placement = "outside", axis.line.x = element_blank(),
        # panel.border = element_rect(color = "grey85", fill = alpha("white", 0)),
        panel.spacing.y = unit(0.5, "line"), panel.grid.major.y = element_line(color = "grey85"),
        legend.position = c(0.89, 0.97),  legend.title.align = 1, 
        legend.margin = margin(0, 2, 2, 2))
)


# Plot boxplot for accuracy
g_rmsep_within_loc <- accuracy_within_loc_toplot %>%
  ggplot(aes(x = selection, y = rmse, color = model, fill = model, group = group)) +
  # jitter points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.1, color = "grey85") +
  geom_boxplot(fill = alpha("white", 0), position = position_dodge(0.9), outlier.shape = NA, lwd = 0.35) +
  scale_y_continuous(name = "Within-location RMSEP", breaks = pretty) +
  g_mod +
  theme(legend.position = "top")



## Save
ggsave(filename = "figure_S8_lolo_bias_within_locations.jpg", plot = g_rmsep_within_loc, 
       path = fig_dir, width = figwidth_onecol + 1, height = 5, dpi = dpi_use)


# Figure S9 - external predictive ability across environments ----------

## Plot all
g_ext_env_all_list <- loo_pred_obs_df %>%
  filter(type == "env_external", selection != "LASSO") %>%
  split(.$trait) %>%
  map(~{
    
    # Get the annotation for across environment prediction
    across_env_ann <- distinct(select(.x, trait, pop, type, ability_ann, model, selection)) %>% 
      mutate(annotation = paste0("'All:'~", ability_ann),
             x = case_when(trait == "TestWeight" & selection == "StepwiseCV" ~ -Inf,
                           TRUE ~ Inf),
             y = case_when(trait == "GrainYield" & str_detect(selection, "None|StepwiseCV", negate = TRUE) ~ Inf,
                           TRUE ~ -Inf),
             test_group = "All")
             
                           
    # Annotation for accuracy within environments
    within_env_ann <- distinct(within_environment_prediction_accuracy, trait, pop, type, selection, model, test_group, ability) %>% 
      mutate(model = f_model_replace(model), 
             selection = f_ec_selection_replace(selection, parse = FALSE),
             selection = ifelse(selection == "None" & model == "g", "Covariate~Set: None", selection),
             model = ifelse(model == "g", "Model: g", model)) %>%
      left_join(distinct(.x, trait, pop, type, selection, model, test_group), .) %>%
      mutate(annotation = paste0("'", test_group, ":'~r[MP]==", format_numbers(ability, 2)))
    
    # Combine annotations
    ann_df <- select(across_env_ann, -annotation) %>% 
      full_join(., within_env_ann) %>%
      bind_rows(., across_env_ann)

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
      facet_grid(. ~ model + selection, switch = "y", scales = "free_x",
                 labeller = labeller(selection = label_parsed)) +
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
    ann_df <- distinct(select(.x, trait, pop, type,  test_group, ability_ann, model, selection)) %>% 
      mutate(test_group1 =  ifelse(str_detect(test_group, "[0-9]{2}"), test_group, 
                                        str_to_title(str_replace_all(test_group, "_", " "))),
             annotation = paste0(abbreviate(test_group1, 5), "~", ability_ann)) %>%
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
      facet_grid(. ~ model + selection, switch = "y", scales = "free_x",
                 labeller = labeller(selection = label_parsed)) +
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
  filter(type %in% c("lolo", "loc_external"), str_detect(selection, "stepwise"), model == "g + e + (ge)") %>%
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

# Display the means
env_means %>%
  mutate(mean = mu + effect) %>%
  group_by(trait) %>%
  summarize_at(vars(mean), list(min = min, max = max))

# trait            min    max
# 1 GrainProtein    9.98   14.7
# 2 GrainYield   1612.   9056. 
# 3 HeadingDate    50.9    70.5
# 4 PlantHeight    45.2   107. 
# 5 TestWeight    568.    713


## Calculate environmental correlation
## Use all available lines

env_phenoCor <- S2_MET_BLUEs %>%
  # filter(line_name %in% tp) %>%
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




# Average train/holdout environment phenotypic correlations
env_phenoCor %>%
  mutate(phenoCor = map(phenoCor, ~as.dist(.) %>% broom::tidy())) %>%
  unnest(phenoCor) %>%
  select(-phenoCor_df) %>%
  rename_at(vars(contains("item")), ~str_replace(., "item", "environment")) %>%
  filter(environment1 %in% train_test_env) %>%
  mutate(environment2_group = ifelse(environment2 %in% train_test_env, "train", "holdout")) %>%
  group_by(trait, environment2_group) %>%
  summarize(env1_env2_dist = mean(distance))




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



# ## Plot the distance between environmental means versus the correlation due to covariates
# env_mean_concurrent_features <- concurrent_features_env_cormat %>%
#   left_join(., env_mean_dist) %>%
#   # mutate(corMat = map2(corMat, dist, ~.x[row.names(.y), row.names(.y)]))
#   mutate(corMat = map2(main, dist, ~.x[row.names(.y), row.names(.y)]))
# 
# env_mean_cor_concurrent_features <- env_mean_concurrent_features %>%
#   mutate(corMat_envMean_df = map2(corMat, dist, ~{
#     .x[lower.tri(.x, diag = TRUE)] <- NA
#     .y[lower.tri(.y, diag = TRUE)] <- NA
#     
#     .x1 <- as.data.frame(.x) %>% 
#       rownames_to_column("environment") %>% 
#       gather(environment2, correlation, -environment) %>% 
#       filter(!is.na(correlation))
#     
#     .y1 <- as.data.frame(.y) %>% 
#       rownames_to_column("environment") %>% 
#       gather(environment2, env_mean_dist, -environment) %>% 
#       filter(!is.na(env_mean_dist))
#     
#     full_join(.x1, .y1, by = c("environment", "environment2"))
#     
#   })) %>% ungroup() %>%
#   # Calculate correlations and annotate
#   mutate(cor_test = map(corMat_envMean_df, ~cor.test(x = .x$correlation, y = .x$env_mean_dist)),
#          cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
#          annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))
# 
# ## Plot
# g_env_mean_cor_concurrent_features <- env_mean_cor_concurrent_features %>%
#   filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
#   unnest(corMat_envMean_df) %>%
#   mutate(env_mean_dist = ifelse(trait == "GrainYield", env_mean_dist / 1000, env_mean_dist)) %>%
#   ggplot(aes(x = correlation, y = env_mean_dist)) +
#   geom_point(size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 0.5, color = "blue") +
#   geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
#   facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y", 
#              labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
#   scale_x_continuous(name = "Covariate-estimated\nenvironmental correlation", breaks = pretty) +
#   scale_y_continuous(name = "Difference in environmental mean", breaks = pretty) +
#   theme_genetics(base_size = 8) +
#   theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")
# 
# # Save
# ggsave(filename = "concurrent_ec_correlation_vs_env_mean_dist.jpg", plot = g_env_mean_cor_concurrent_features,
#        path = fig_dir, height = 5, width = 5, dpi = 1000)



# ## Plot the phenotypic correlation versus the correlation due to covariates
# phenCor_concurrent_features <- concurrent_features_env_cormat %>%
#   left_join(., env_phenoCor) %>%
#   mutate(corMat = map2(corMat, phenoCor, ~.x[row.names(.y), row.names(.y)]))
# 
# phenCor_concurrent_features1 <- phenCor_concurrent_features %>%
#   mutate(phenoCor_corMat = map2(corMat, phenoCor, ~{
#     .x[lower.tri(.x, diag = TRUE)] <- NA
#     .y[lower.tri(.y, diag = TRUE)] <- NA
#     
#     .x1 <- as.data.frame(.x) %>% 
#       rownames_to_column("environment") %>% 
#       gather(environment2, correlation, -environment) %>% 
#       filter(!is.na(correlation))
#     
#     .y1 <- as.data.frame(.y) %>% 
#       rownames_to_column("environment") %>% 
#       gather(environment2, phenoCor, -environment) %>% 
#       filter(!is.na(phenoCor))
#     
#     full_join(.x1, .y1, by = c("environment", "environment2"))
#     
#   })) %>% ungroup() %>%
#   # Calculate correlations and annotate
#   mutate(cor_test = map(phenoCor_corMat, ~cor.test(x = .x$correlation, y = .x$phenoCor)),
#          cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
#          annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))
# 
# ## Plot
# g_phenCor_concurrent_features1 <- phenCor_concurrent_features1 %>%
#   filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
#   unnest(phenoCor_corMat, names_repair = tidyr_legacy) %>%
#   ggplot(aes(x = correlation, y = phenoCor1)) +
#   geom_point(size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
#   geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
#   facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
#              labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
#   scale_x_continuous(name = "Covariate-estimated\nenvironmental correlation", breaks = pretty) +
#   scale_y_continuous(name = "Phenotypic correlation between environments", breaks = pretty) +
#   theme_genetics(base_size = 8) +
#   theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")
# 
# # Save
# ggsave(filename = "concurrent_ec_correlation_vs_phenoCor.jpg", plot = g_phenCor_concurrent_features1,
#        path = fig_dir, height = 5, width = 5, dpi = 1000)




# Plot a joint heatmap of the estimate phenotypic correlation versus the correlation
# between environments estimated by the EC-covariance matrix
phenoCor_ec_heatmaps <- concurrent_features_env_cormat %>%
  left_join(., env_phenoCor) %>%
  group_by(trait, feat_sel_type) %>%
  do(plot = {
    row <- .
    
    # Find the common set of environments
    common_envs <- intersect(row.names(row$phenoCor[[1]]), row.names(row$main[[1]]))
    
    # Subset each correlation matrix, convert to a df
    pheno_cor_df <- row$phenoCor[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("environment1") %>%
      gather(environment2, phenoCor, -environment1) %>%
      filter_at(vars(contains("environment")), all_vars(. %in% common_envs))
    
    # Calculate the average phenotypic correlation of each environment with all others; each with 
    # only training environments
    pheno_cor_df1 <- pheno_cor_df %>% 
      mutate(environment2_training = environment2 %in% train_test_env)
    
    avg_phenocor <- bind_rows(
      mutate(aggregate(phenoCor ~ environment1, data = pheno_cor_df1, FUN = mean), scope = "all_env"),
      mutate(aggregate(phenoCor ~ environment1, data = pheno_cor_df1, FUN = mean, subset = environment2_training), scope = "train_env")
    ) %>%
      rename(avg_phenoCor = phenoCor)
    
    ec_cor_df <- row$main[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("environment1") %>%
      gather(environment2, ecCor, -environment1) %>%
      filter_at(vars(contains("environment")), all_vars(. %in% common_envs))
    
    
    # Cluster based on phenotypic correlations
    clust <- hclust(as.dist(-row$phenoCor[[1]]))
    # Order of factors
    env_fact_order <- factor(common_envs, levels = clust$labels[clust$order])
    
    # Upper triangle df
    upper_tri_df <- row$corMat[[1]][levels(env_fact_order), levels(env_fact_order)] %>%
      upper.tri() %>% 
      as.data.frame() %>% 
      rownames_to_column("environment1") %>% 
      gather(environment2, upper_tri, -environment1) %>%
      mutate_at(vars(contains("environment")), ~parse_number(.) %>% {levels(env_fact_order)[.]})
    
    # Combine df
    cor_df <- full_join(pheno_cor_df, ec_cor_df, by = c("environment1", "environment2")) %>%
      # Decide what elements will go on the upper versus lower diagonal
      full_join(., upper_tri_df, by = c("environment1", "environment2")) %>%
      # Add the average phenotypic correlations
      left_join(., subset(avg_phenocor, scope == "train_env"), by = c("environment1")) %>%
      mutate(heatvalue = ifelse(upper_tri, phenoCor, ecCor),
             heatvalue = ifelse(environment1 == environment2, as.numeric(NA), heatvalue),
             label = ifelse(environment1 == environment2, format_numbers(avg_phenoCor, 2), as.character(NA))) %>%
      # Set environments as factors with levels in the following order: first holdout environments, then 
      # based on phenocor clustering
      mutate_at(vars(contains("environment")), ~factor(., levels = levels(env_fact_order))) %>%
      mutate(x = "Phenotypic correlation", y = "EC-correlation")

    # Bold font for validation environments
    font_face <- ifelse(levels(env_fact_order) %in% validation_env, "bold", "plain")
    
    # plot the heatmap
    g_heat <- cor_df %>%
      ggplot(aes(x = environment1, y = environment2)) +
      # geom_tile(aes(fill = correlation), color = ifelse(is.na(dat2$correlation), "black", NA), lwd = 0.1) +
      geom_tile(aes(fill = heatvalue)) +
      # Include text
      geom_text(aes(label = label), size = 2) +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3], na.value = "white",
                           midpoint = 0, limits = c(-1.01, 1.01), breaks = pretty,
                           name = "Correlation", guide = guide_colorbar(title.position = "top")) +
      facet_grid(x ~ y, switch = "both") +
      labs(subtitle = str_add_space(row$trait)) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = font_face), 
            axis.text.y = element_text(face = font_face), axis.title = element_blank(),
            strip.background = element_blank())
    
    # Output the heatmap
    g_heat
    
  }) %>% ungroup()




## Combine plots for use in the supplemental figure

g_heatmap_combined <- plot_grid(plotlist = subset(phenoCor_ec_heatmaps, feat_sel_type == "stepwise_cv_adhoc", plot, drop = TRUE),
                                ncol = 2)


# Save
ggsave(filename = "figure_s12_environmental_covariate_heatmap.jpg", plot = g_heatmap_combined,
       path = fig_dir, width = 12, height = 12, dpi = 1000)






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
  # filter(line_name %in% tp) %>%
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

# Average train/holdout environment phenotypic correlations
loc_phenoCor %>%
  mutate(phenoCor = map(phenoCor, ~as.dist(.) %>% broom::tidy())) %>%
  unnest(phenoCor) %>%
  select(-phenoCor_df) %>%
  rename_at(vars(contains("item")), ~str_replace(., "item", "environment")) %>%
  filter(environment1 %in% train_test_loc) %>%
  mutate(environment2_group = ifelse(environment2 %in% train_test_loc, "train", "holdout")) %>%
  group_by(trait, environment2_group) %>%
  summarize(env1_env2_dist = mean(distance))




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



# ## Plot the distance between environmental means versus the correlation due to covariates
# env_mean_concurrent_features <- historical_features_loc_cormat %>%
#   left_join(., loc_mean_dist) %>%
#   # mutate(corMat = map2(corMat, dist, ~.x[row.names(.y), row.names(.y)]))
#   mutate(corMat = map2(main, dist, ~.x[row.names(.y), row.names(.y)]))
# 
# env_mean_cor_concurrent_features <- env_mean_concurrent_features %>%
#   mutate(corMat_envMean_df = map2(corMat, dist, ~{
#     .x[lower.tri(.x, diag = TRUE)] <- NA
#     .y[lower.tri(.y, diag = TRUE)] <- NA
# 
#     .x1 <- as.data.frame(.x) %>%
#       rownames_to_column("environment") %>%
#       gather(environment2, correlation, -environment) %>%
#       filter(!is.na(correlation))
# 
#     .y1 <- as.data.frame(.y) %>%
#       rownames_to_column("environment") %>%
#       gather(environment2, env_mean_dist, -environment) %>%
#       filter(!is.na(env_mean_dist))
# 
#     full_join(.x1, .y1, by = c("environment", "environment2"))
# 
#   })) %>% ungroup() %>%
#   # Calculate correlations and annotate
#   mutate(cor_test = map(corMat_envMean_df, ~cor.test(x = .x$correlation, y = .x$env_mean_dist)),
#          cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
#          annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))
# 
# ## Plot
# g_env_mean_cor_concurrent_features <- env_mean_cor_concurrent_features %>%
#   filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
#   unnest(corMat_envMean_df) %>%
#   mutate(env_mean_dist = ifelse(trait == "GrainYield", env_mean_dist / 1000, env_mean_dist)) %>%
#   ggplot(aes(x = correlation, y = env_mean_dist)) +
#   geom_point(size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 0.5, color = "blue") +
#   geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
#   facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
#              labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
#   scale_x_continuous(name = "Covariate-estimated\nlocation correlation", breaks = pretty) +
#   scale_y_continuous(name = "Difference in location mean", breaks = pretty) +
#   theme_genetics(base_size = 8) +
#   theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")
# 
# 
# ## Plot the phenotypic correlation versus the correlation due to covariates
# phenCor_concurrent_features <- historical_features_loc_cormat %>%
#   left_join(., loc_phenoCor) %>%
#   mutate(corMat = map2(corMat, phenoCor, ~.x[row.names(.y), row.names(.y)]))
# 
# phenCor_concurrent_features1 <- phenCor_concurrent_features %>%
#   mutate(phenoCor_corMat = map2(corMat, phenoCor, ~{
#     .x[lower.tri(.x, diag = TRUE)] <- NA
#     .y[lower.tri(.y, diag = TRUE)] <- NA
# 
#     .x1 <- as.data.frame(.x) %>%
#       rownames_to_column("environment") %>%
#       gather(environment2, correlation, -environment) %>%
#       filter(!is.na(correlation))
# 
#     .y1 <- as.data.frame(.y) %>%
#       rownames_to_column("environment") %>%
#       gather(environment2, phenoCor, -environment) %>%
#       filter(!is.na(phenoCor))
# 
#     full_join(.x1, .y1, by = c("environment", "environment2"))
# 
#   })) %>% ungroup() %>%
#   # Calculate correlations and annotate
#   mutate(cor_test = map(phenoCor_corMat, ~cor.test(x = .x$correlation, y = .x$phenoCor)),
#          cor_estimate = map_dbl(cor_test, "estimate"), p_value = map_dbl(cor_test, "p.value"),
#          annotation = paste0("r = ", format_numbers(x = cor_estimate, 2), "; P = ", formatC(x = p_value, digits = 2, width = 2, format = "g")))
# 
# ## Plot
# g_phenCor_concurrent_features1 <- phenCor_concurrent_features1 %>%
#   filter(str_detect(feat_sel_type, "AIC", negate = TRUE)) %>%
#   unnest(phenoCor_corMat) %>%
#   ggplot(aes(x = correlation, y = phenoCor)) +
#   geom_point(size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
#   geom_text(aes(x = Inf, y = Inf, label = annotation), color = "blue", size = base_geom_text_size, hjust = 1.05, vjust = 1.2) +
#   facet_grid(trait ~ feat_sel_type, scales = "free_y", switch = "y",
#              labeller = labeller(feat_sel_type = f_ec_selection_replace, trait = str_add_space)) +
#   scale_x_continuous(name = "Covariate-estimated\nlocation correlation", breaks = pretty) +
#   scale_y_continuous(name = "Phenotypic correlation between locations", breaks = pretty) +
#   theme_genetics(base_size = 8) +
#   theme(panel.border = element_rect(color = "grey85", fill = alpha("white", 0)), strip.placement = "outside")


# Function to abbreviate location names
f_abbreviate_locations <- function(x) {
  locations_short <- distinct(trial_info, location, environment) %>%
    mutate(environment = str_sub(environment, 1, 3)) %>% 
    distinct() %>%
    {setNames(.$environment, .$location)}
  locations_short[x]
}



# Plot a joint heatmap of the estimate phenotypic correlation versus the correlation
# between environments estimated by the EC-covariance matrix
loc_phenoCor_ec_heatmaps <- historical_features_loc_cormat %>%
  left_join(., loc_phenoCor) %>%
  group_by(trait, feat_sel_type) %>%
  do(plot = {
    row <- .
    
    # Find the common set of environments
    common_envs <- intersect(row.names(row$phenoCor[[1]]), row.names(row$main[[1]]))
    
    # Subset each correlation matrix, convert to a df
    pheno_cor_df <- row$phenoCor[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("environment1") %>%
      gather(environment2, phenoCor, -environment1) %>%
      filter_at(vars(contains("environment")), all_vars(. %in% common_envs))
    
    # Calculate the average phenotypic correlation of each environment with all others; each with 
    # only training environments
    pheno_cor_df1 <- pheno_cor_df %>% 
      mutate(environment2_training = environment2 %in% train_test_loc)
    
    avg_phenocor <- bind_rows(
      mutate(aggregate(phenoCor ~ environment1, data = pheno_cor_df1, FUN = mean), scope = "all_env"),
      mutate(aggregate(phenoCor ~ environment1, data = pheno_cor_df1, FUN = mean, subset = environment2_training), scope = "train_env")
    ) %>%
      rename(avg_phenoCor = phenoCor)
    
    ec_cor_df <- row$main[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("environment1") %>%
      gather(environment2, ecCor, -environment1) %>%
      filter_at(vars(contains("environment")), all_vars(. %in% common_envs))
    
    
    # Cluster based on phenotypic correlations
    clust <- hclust(as.dist(-row$phenoCor[[1]]))
    # Order of factors
    env_fact_order <- factor(common_envs, levels = clust$labels[clust$order])
    
    # Upper triangle df
    upper_tri_df <- row$corMat[[1]][levels(env_fact_order), levels(env_fact_order)] %>%
      upper.tri() %>% 
      as.data.frame() %>% 
      rownames_to_column("environment1") %>% 
      gather(environment2, upper_tri, -environment1) %>%
      mutate_at(vars(contains("environment")), ~parse_number(.) %>% {levels(env_fact_order)[.]})
    
    # Combine df
    cor_df <- full_join(pheno_cor_df, ec_cor_df, by = c("environment1", "environment2")) %>%
      # Decide what elements will go on the upper versus lower diagonal
      full_join(., upper_tri_df, by = c("environment1", "environment2")) %>%
      # Add the average phenotypic correlations
      left_join(., subset(avg_phenocor, scope == "train_env"), by = c("environment1")) %>%
      mutate(heatvalue = ifelse(upper_tri, phenoCor, ecCor),
             heatvalue = ifelse(environment1 == environment2, as.numeric(NA), heatvalue),
             label = ifelse(environment1 == environment2, format_numbers(avg_phenoCor, 2), as.character(NA))) %>%
      # Set environments as factors with levels in the following order: first holdout environments, then 
      # based on phenocor clustering
      mutate_at(vars(contains("environment")), ~factor(., levels = levels(env_fact_order))) %>%
      mutate(x = "Phenotypic correlation", y = "EC-correlation")
    
    # Bold font for validation environments
    font_face <- ifelse(levels(env_fact_order) %in% validation_loc, "bold", "plain")
    
    # plot the heatmap
    g_heat <- cor_df %>%
      ggplot(aes(x = environment1, y = environment2)) +
      # geom_tile(aes(fill = correlation), color = ifelse(is.na(dat2$correlation), "black", NA), lwd = 0.1) +
      geom_tile(aes(fill = heatvalue)) +
      # Include text
      geom_text(aes(label = label), size = 2) +
      scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[2], high = heat_colors[3], na.value = "white",
                           midpoint = 0, limits = c(-1.01, 1.01), breaks = pretty,
                           name = "Correlation", guide = guide_colorbar(title.position = "top")) +
      scale_x_discrete(labels = f_abbreviate_locations) +
      scale_y_discrete(labels = f_abbreviate_locations) +
      facet_grid(x ~ y, switch = "both") +
      labs(subtitle = str_add_space(row$trait)) +
      theme_genetics(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = font_face), 
            axis.text.y = element_text(face = font_face), axis.title = element_blank(),
            strip.background = element_blank())
    
    # Output the heatmap
    g_heat
    
  }) %>% ungroup()



## Combine plots for use in the supplemental figure

g_heatmap_combined <- plot_grid(plotlist = subset(loc_phenoCor_ec_heatmaps, feat_sel_type == "stepwise_cv_adhoc", plot, drop = TRUE),
                                ncol = 2)


# Save
ggsave(filename = "figure_s13_location_covariate_heatmap.jpg", plot = g_heatmap_combined,
       path = fig_dir, width = 12, height = 12, dpi = 1000)






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


# Table S4: Holdout environment/location accuracy --------------------------


# External environment predictions
env_external_pred_accuracy <- predictive_ability %>%
  filter(type == "env_external", (is.na(time_frame_selection) | time_frame_selection == "bestOverall"),
         str_detect(selection, "lasso|nosoil", negate = TRUE)) %>%
  distinct(trait, site, selection, model, pop, ability, bias, rmse) %>% 
  left_join(., env_trait_herit, by = c("trait", "site" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability))

as.data.frame(env_external_pred_accuracy)

(env_external_pred_accuracy_summary <- env_external_pred_accuracy %>%
    group_by(trait, pop, selection, model) %>%
    summarize_at(vars(ability, accuracy, bias, rmse), mean) %>%
    ungroup() %>%
    arrange(model, pop, trait) %>%
    as.data.frame())

(env_external_pred_accuracy1 <- env_external_pred_accuracy %>%
    group_by(trait, pop, selection, model, site) %>%
    summarize_at(vars(ability, accuracy, bias, rmse), mean) %>%
    ungroup() %>%
    arrange(model, pop, trait, site) %>%
    as.data.frame())



# Create a table to save
table_towrite <- env_external_pred_accuracy1 %>%
  select(trait, selection, model, site, accuracy, rmse) %>%
  # Add summary
  bind_rows(., select(mutate(env_external_pred_accuracy_summary, site = "Mean"), trait, selection, model, site, accuracy, rmse)) %>%
  mutate(selection = f_ec_selection_replace(selection, parse = FALSE) %>% fct_inorder(),
         model = f_model_replace(model)) %>%
  mutate_at(vars(accuracy, rmse), ~format_numbers(.) %>% str_remove(., "\\.$")) %>%
  mutate(annotation = paste0(accuracy, " (", rmse, ")")) %>%
  select(-accuracy, -rmse) %>%
  # Rename traits
  mutate(trait = str_add_space(trait)) %>%
  spread(trait, annotation) %>%
  # mutate(scenario = "Holdout environments") %>%
  select(model, covariate_set = selection, names(.)) %>%
  arrange(model, covariate_set) %>%
  group_by(model, covariate_set) %>%
  # This code keeps only the first element in a group
  do({
    df <- .
    df1 <- mutate(df, model = "", covariate_set = "")
    df1$model[1] <- as.character(unique(df$model))
    df1$covariate_set[1] <- as.character(unique(df$covariate_set))
    df1
  }) %>%
  # rename(prediction_scenario = scenario, EC_set = covariate_set) %>%
  rename_all(~str_replace_all(., "_", " ")) %>%
  rename_at(vars(-matches(paste0(str_add_space(traits), collapse = "|"))), str_to_sentence)

ft1 <- flextable(data = table_towrite) %>% 
  autofit() %>%
  bold(i = which(table_towrite$Site == "Mean"), j = seq_len(ncol(table_towrite)), bold = TRUE) %>%
  border(i = which(table_towrite$Site == "Mean"), j = seq_len(ncol(table_towrite)),
         border.bottom = fp_border())


# Same thing for locations
# External location predictions
loc_external_pred_accuracy <- predictive_ability %>%
  filter(type == "loc_external", (is.na(time_frame_selection) | time_frame_selection == "bestOverall"),
         str_detect(selection, "lasso|nosoil", negate = TRUE)) %>%
  distinct(trait, site, selection, model, pop, ability, bias, rmse)

as.data.frame(loc_external_pred_accuracy)

(loc_external_pred_accuracy_summary <- loc_external_pred_accuracy %>%
    group_by(trait, pop, selection, model) %>%
    summarize_at(vars(ability, ability, bias, rmse), mean) %>%
    ungroup() %>%
    arrange(model, pop, trait) %>%
    as.data.frame())

(loc_external_pred_accuracy1 <- loc_external_pred_accuracy %>%
    group_by(trait, pop, selection, model, site) %>%
    summarize_at(vars(ability, ability, bias, rmse), mean) %>%
    ungroup() %>%
    arrange(model, pop, trait, site) %>%
    as.data.frame())

# Create a table to save
table_towrite2 <- loc_external_pred_accuracy1 %>%
  select(trait, selection, model, site, ability, rmse) %>%
  # Rename sites
  left_join(., distinct(mutate(trial_info, location_code = str_sub(environment, 1, 3)), site = location, location_code)) %>%
  mutate(site = location_code) %>% select(-location_code) %>%
  # Add summary
  bind_rows(., select(mutate(loc_external_pred_accuracy_summary, site = "Mean"), trait, selection, model, site, ability, rmse)) %>%
  mutate(selection = f_ec_selection_replace(selection, parse = FALSE) %>% fct_inorder(),
         model = f_model_replace(model)) %>%
  mutate_at(vars(ability, rmse), ~format_numbers(.) %>% str_remove(., "\\.$")) %>%
  mutate(annotation = paste0(ability, " (", rmse, ")")) %>%
  select(-ability, -rmse) %>%
  # Rename traits
  mutate(trait = str_add_space(trait)) %>%
  # Edit model
  mutate(model = str_replace_all(model, "e", "l")) %>%
  spread(trait, annotation) %>%
  select(model, covariate_set = selection, names(.)) %>%
  arrange(model, covariate_set) %>%
  group_by(model, covariate_set) %>%
  # This code keeps only the first element in a group
  do({
    df <- .
    df1 <- mutate(df, model = "", covariate_set = "")
    df1$model[1] <- as.character(unique(df$model))
    df1$covariate_set[1] <- as.character(unique(df$covariate_set))
    df1
  }) %>%
  # rename(prediction_scenario = scenario, EC_set = covariate_set) %>%
  rename_all(~str_replace_all(., "_", " ")) %>%
  rename_at(vars(-matches(paste0(str_add_space(traits), collapse = "|"))), str_to_sentence)


ft2 <- flextable(data = table_towrite2) %>% 
  autofit() %>%
  bold(i = which(table_towrite2$Site == "Mean"), j = seq_len(ncol(table_towrite2)), bold = TRUE) %>%
  border(i = which(table_towrite2$Site == "Mean"), j = seq_len(ncol(table_towrite2)),
         border.bottom = fp_border())

# Render the tables
save_as_docx("Table S4" = ft1, "Table S5" = ft2, path = file.path(fig_dir, "supplemental_tables_s4_s5.docx"))

# Table S5: Holdout environment/location accuracy --------------------------


## Compare average accuracy in the LOLO or LEOE set versus the holdout environments/locations
loeo_prediction_accuracy_summary <- predictive_ability %>%
  filter(type == "loeo", (is.na(time_frame_selection) | time_frame_selection == "bestOverall"),
         str_detect(selection, "lasso|nosoil", negate = TRUE)) %>%
  distinct(trait, site, selection, model, pop, ability, bias, rmse) %>% 
  left_join(., env_trait_herit, by = c("trait", "site" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability)) %>%
  group_by(trait, pop, selection, model) %>%
  summarize_at(vars(ability, accuracy, bias, rmse), mean) %>%
  ungroup() %>%
  arrange(model, pop, trait) %>%
  as.data.frame()

loeo_envExternal_accuracy_summary <- full_join(loeo_prediction_accuracy_summary, env_external_pred_accuracy_summary,
                                               by = c("trait", "pop", "selection", "model")) %>%
  filter(pop == "vp") %>%
  mutate(ability_cv_ts_diff = ability.x - ability.y, ability_cv_ts_diff_per = ability_cv_ts_diff / ability.x,
         accuracy_cv_ts_diff = accuracy.x - accuracy.y, accuracy_cv_ts_diff_per = accuracy_cv_ts_diff / accuracy.x,
         rmse_cv_ts_diff = rmse.x - rmse.y, rmse_cv_ts_diff_per = rmse_cv_ts_diff / rmse.x) %>%
  # Summarize by trait
  group_by(trait) %>%
  summarize_at(vars(contains("cv_ts_diff")), list(mean = mean, min = min, max = max)) %>%
  as.data.frame()

loeo_envExternal_accuracy_summary


lolo_prediction_accuracy_summary <- predictive_ability %>%
  filter(type == "lolo", (is.na(time_frame_selection) | time_frame_selection == "bestOverall"),
         str_detect(selection, "lasso|nosoil", negate = TRUE)) %>%
  distinct(trait, site, selection, model, pop, ability, bias, rmse) %>% 
  left_join(., env_trait_herit, by = c("trait", "site" = "environment")) %>%
  mutate(accuracy = ability / sqrt(heritability)) %>%
  group_by(trait, pop, selection, model) %>%
  summarize_at(vars(ability, accuracy, bias, rmse), mean) %>%
  ungroup() %>%
  arrange(model, pop, trait) %>%
  as.data.frame()

lolo_envExternal_accuracy_summary <- full_join(lolo_prediction_accuracy_summary, loc_external_pred_accuracy_summary,
                                               by = c("trait", "pop", "selection", "model")) %>%
  filter(pop == "vp") %>%
  mutate(ability_cv_ts_diff = ability.x - ability.y, ability_cv_ts_diff_per = ability_cv_ts_diff / ability.x,
         # accuracy_cv_ts_diff = accuracy.x - accuracy.y, accuracy_cv_ts_diff_per = accuracy_cv_ts_diff / accuracy.x,
         rmse_cv_ts_diff = rmse.x - rmse.y, rmse_cv_ts_diff_per = rmse_cv_ts_diff / rmse.x) %>%
  # Summarize by trait
  group_by(trait) %>%
  summarize_at(vars(contains("cv_ts_diff")), list(mean = mean, min = min, max = max)) %>%
  as.data.frame()





# Other notes to reference in the manuscript -------

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
  arrange(trait, model, pop, desc(ability)) %>%
  as.data.frame() %>% filter(model == "model5_cov", selection == "stepwise_cv_adhoc") %>% arrange(pop)


predictive_ability %>%
  filter(type == "lolo", selection %in% c("none", "stepwise_cv_adhoc"))





# Summarize best model by trait


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


# How do LOLO predictions compared to LOEO averages over years?
# 
# First calculate LOEO predictions averaged over years
loeo_predictions_per_location <- predictions_df %>%
  filter(type == "loeo", model %in% models_use, str_detect(selection, "nosoil|lasso", negate = TRUE)) %>%
  # Convert GY to Mg ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  # Remove all NA columns
  select_if(map_lgl(., ~any(!is.na(.)))) %>%
  # add location information
  left_join(., distinct(select(trial_info, site = environment, location))) %>%
  group_by(trait, location, model, selection, pop, line_name) %>%
  summarize_at(vars(contains("pred")), mean, na.rm = TRUE) %>%
  ungroup()

# New df to merge
loo_pred_obs_df1 <- predictions_df %>%
  filter(type == "lolo", model %in% models_use, str_detect(selection, "nosoil|lasso", negate = TRUE),
         (is.na(time_frame_selection) | time_frame_selection == "bestOverall")) %>%
  # convert grain yield to Mg/ha
  mutate_at(vars(value, pred_complete), ~ifelse(trait == "GrainYield", . / 1000, .)) %>%
  select_if(map_lgl(., ~any(!is.na(.))))

# Add and trim the correct time frames for each trait and model
loo_pred_obs_df2 <- loo_pred_obs_df1 %>%
  filter(selection %in% c("none", "stepwise_cv_adhoc")) %>%
  distinct(trait, model, time_frame) %>%
  # Inner join
  inner_join(loo_pred_obs_df1, .)
  
loeo_average_location_predictions <- loeo_predictions_per_location %>%
  left_join(., select(loo_pred_obs_df2, trait, location = site, model, selection, time_frame, 
                      time_frame_selection, pop, line_name, loc_mean_value = value)) %>%
  # This removes Aberdeen as a location - dryland
  filter(!is.na(loc_mean_value))
  
loeo_average_location_predictions_accuracy_byLoc <- loeo_average_location_predictions %>%
  group_by(trait, model, selection, time_frame, time_frame_selection, pop, location) %>%
  summarize(ability = cor(pred_complete, loc_mean_value)) %>%
  # Take the average accuracy of all locations
  summarize_at(vars(ability), list(mean = mean, min = min, max = max)) %>%
  ungroup()

## Compare average predictive abilities

loeo_average_location_predictions_accuracy_byLoc1 <- predictive_ability %>% 
  filter(type == "lolo") %>%
  group_by(trait, model, selection, time_frame, pop) %>%
  summarize_at(vars(ability), list(lolo_mean = mean, lolo_min = min, lolo_max = max)) %>%
  ungroup() %>%
  left_join(loeo_average_location_predictions_accuracy_byLoc, .)

# Plot
loeo_average_location_predictions_accuracy_byLoc1 %>%
  ggplot(aes(x = ))
  
  
loeo_average_location_predictions_accuracy_acrossLoc <- loeo_average_location_predictions %>%
  group_by(trait, model, selection, time_frame, time_frame_selection, pop) %>%
  summarize(ability = cor(pred_complete, loc_mean_value), .groups = "drop")



# Examine outlier environments for rmsep for yield
outlier_env_yield <- predictive_ability %>%
  filter(trait == "GrainYield", type == "loeo", model == "model3_cov", selection == "stepwise_cv_adhoc", pop == "tp") %>%
  filter(rmse == max(rmse))

# Plot pred vs obs values
predictions_df %>%
  left_join(outlier_env_yield, .) %>%
  plot(value ~ pred_complete, .)

predictions_df %>%
  left_join(select(outlier_env_yield, trait, type, model, selection, pop), .) %>%
  filter(str_detect(site, "WLI")) %>%
  ggplot(aes(x = pred_complete, y = value, color = site))+
  geom_point()

left_join(outlier_env_yield, within_environment_prediction_accuracy) %>%
  as.data.frame()


env_means %>%
  filter(str_detect(environment, "WLI")) %>%
  mutate(mean = mu + effect)
  
ec_tomodel_centered$daymet %>% 
  filter(str_detect(environment, "WLI")) %>% 
  gather(covariate, value, -source, -environment) %>%
  spread(environment, value) %>%
  as.data.frame()
