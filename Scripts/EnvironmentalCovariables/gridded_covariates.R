## Example nasapower data for a region

# Repository directory
repo_dir <- getwd()
# Source the main project script
source(file.path(repo_dir, "source.R"))

# Other packages
library(nasapower)


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
  filter(environment %in% unique(S2_MET_BLUEs$environment)) %>%
  group_by(location, latitude, longitude) %>%
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
long_limit <- c(-111, -72)
lat_limit <- c(39, 49)



## Different version of the map - grey
g_map_alt <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "grey85") +
  geom_polygon(data = canada, fill = NA, color = "white", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "white", lwd = 0.3) +
  coord_fixed(ratio = 1.5, xlim = long_limit, ylim = lat_limit) +
  scale_x_continuous(breaks = NULL, name = NULL, labels = NULL) + 
  scale_y_continuous(breaks = NULL, name = NULL, labels = NULL) +
  theme_classic() +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = alpha("white", 0)), axis.line = element_blank())


## Break up the longitude and latitude limits into 4.5-degree bins
long_bin <- seq.default(min(long_limit), max(long_limit), by = 4.5)
long_bin <- if (max(long_bin) == max(long_limit)) long_bin else c(long_bin, max(long_limit))
lat_bin <- seq.default(min(lat_limit), max(lat_limit), by = 4.5)
lat_bin <- if (max(lat_bin) == max(lat_limit)) lat_bin else c(lat_bin, max(lat_limit))

## Create a list of grids
grid_boxes <- list()
# p = number of list elements
p <- 1

# Iterate over lat (i or y)
for (i in seq(1, length(lat_bin) - 1)) {
  # Iterate over long (j or x)
  for (j in seq(1, length(long_bin) - 1)) {
    
    # Create a vector of xmin, ymin, xmax, ymax
    grid_boxes[[p]] <- c(xmin = long_bin[j], ymin = lat_bin[i], 
                         xmax = long_bin[j+1], ymax = lat_bin[i+1])
    # Increase p
    p <- p + 1
    
  }
}

## Combine into a df
grid_boxes_df <- as_tibble(do.call("rbind", grid_boxes))

# Draw the rectangles
g_map_alt + 
  geom_rect(data = grid_boxes_df, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), 
            inherit.aes = FALSE, fill = alpha("white", 0), color = "black")



## Get the monthly average of parameters within each grid space
# List of parameter
pars <- c("T2M_MIN", "T2M_MAX", "PRECTOT")

# A empty list for storing the data
climate_data_grid_list <- grid_boxes

# Iterate over the list of grids
for (i in seq_along(grid_boxes)) {
  # Get the data and save
  grid_climate_data <- get_power(community = "AG", pars = pars, temporal_average = "CLIMATOLOGY", 
                                 lonlat = grid_boxes[[i]])
  climate_data_grid_list[[i]] <- as_tibble(grid_climate_data)
}


## Combine list into a data.frame
climate_data_grid_df <- bind_rows(climate_data_grid_list) %>%
  rename_all(tolower) %>%
  distinct()
  
  
  

## Plot annual precipitation over the map
g_map_alt +
  geom_tile(data = subset(climate_data_grid_df, parameter == "T2M_MAX"), 
            aes(x = lon, y = lat, fill = ann), alpha = 0.6, inherit.aes = FALSE) +
  scale_fill_viridis_c()
  



