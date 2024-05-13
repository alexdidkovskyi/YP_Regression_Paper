#Load Libraries-----
# Load required libraries
library(openxlsx)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(geosphere)
library(sp)
library(rgdal)
library(proj4)
library(raster)
library(igraph)

sort_by_row <- function(df) {
  df %>%
    rowwise() %>%
    mutate(across(everything(), sort, decreasing = TRUE)) %>%
    ungroup()
}


#Process pebbles' movement data-----
### Preprocess pebbles' locations. 2016-2017
path_to_16_17_data <- 'data/initial/Rilievo2016'
movement_16_17_files <- list.files(path_to_16_17_data)

locations_16_17 <-
  lapply(movement_16_17_files, function(movement_16_17_file)
    read_csv(paste0(
      path_to_16_17_data, '/', movement_16_17_file
    ))) %>%
  do.call('rbind', .)

###  Preprocess  pebbles' characteristics. 2016-2017
characteristics_columns <-
  c('IDREF', 'Weight', 'a_axis', 'b_axis', 'c_axis')

characteristics_16_17 <-
  read.xlsx('data/initial/Pebbles_16_17.xlsx') %>%
  dplyr::select(all_of(characteristics_columns))
characteristics_16_17[, -1] %<>% mutate_if(is.character, as.numeric)


### Read pebbles' locations. 2018
locations_18 <- read.xlsx('data/initial/MERGED2018_XLS.xlsx')

CRS.n <-
  CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

locations_18[, c("X_Start", "Y_Start")] <-
  data.frame(project(locations_18[, c("X_Start", "Y_Start")], CRS.n, inverse =
                       TRUE))
locations_18[, c("X_End", "Y_End")] <-
  data.frame(project(locations_18[, c("X_End", "Y_End")], CRS.n, inverse =
                       TRUE))
### Preprocess  pebbles' characteristics. 2018


characteristics_18 <-
  read.xlsx('data/initial/DB 2018.xlsx', sheet = 'Foglio1') %>%
  set_colnames(characteristics_columns) %>%
  mutate(
    a_axis = ifelse(a_axis > 200, 200, a_axis),
    # fix outliers
    b_axis = ifelse(b_axis > 135, 135, b_axis)
  )


###  Aggregate and filter pebbles' locations and characteristics. 2018
df_loc <- locations_16_17 %>% bind_rows(locations_18)
df_ch  <- characteristics_16_17 %>% bind_rows(characteristics_18)

good_IDREF <- intersect(df_loc$IDREF, df_ch$IDREF)


df_loc %<>% dplyr::filter(IDREF %in% good_IDREF) %>% filter(Y_Start >= 45, X_Start >= 9,
                                                            Y_End >= 45, X_End >= 9)
df_ch %<>% dplyr::filter(IDREF %in% good_IDREF)



df_ch[, c("a_axis", "b_axis", "c_axis")] <-
  sort_by_row(df_ch[, c("a_axis", "b_axis", "c_axis")])

df_ch %<>% mutate(
  Nominal_diameter = (a_axis * b_axis * c_axis) ^ (1 / 3),
  Elongation = b_axis / a_axis,
  Platyness  = c_axis / a_axis,
  Sphericity = (b_axis * c_axis / a_axis ^ 2) ^ (1 / 3)
)




### add missing events for the cases when
### - a pebble didn't move (covered distance is 0) and
### - EventoStart, EventoEnd are not consecutive

# Generate a sequence of event pairs given the start and end values
generate_event_sequence <- function(start, end) {
  data.frame(EventoStart = start:(end - 1),
             EventoEnd = (start + 1):end)
}

# Expand each row of a data frame into multiple rows based on the event sequence
expand_events <- function(df) {
  df %>%
    rowwise() %>%
    # Generate the event sequence for each row
    mutate(event_sequence = list(generate_event_sequence(EventoStart, EventoEnd))) %>%
    # Unnest the event_sequence column
    unnest(event_sequence, names_sep = '_') %>%
    # Update the original EventoStart and EventoEnd columns with the values from the unnested columns
    mutate(EventoStart = event_sequence_EventoStart, EventoEnd = event_sequence_EventoEnd) %>%
    # Select only the required columns
    dplyr::select(names(df),-starts_with("event_sequence"))
}

# Process the location data frame to handle cases where EventoEnd - EventoStart >= 2 and Distance_m == 0
process_loc_data <- function(df_loc) {
  # Filter the data frame to get rows where EventoEnd - EventoStart >= 2 and Distance_m == 0
  df_loc_temp <- df_loc %>%
    filter(EventoEnd - EventoStart >= 2, Distance_m == 0)
  
  
  # Expand the events for each row in df_loc_temp
  df_loc_expanded <- df_loc_temp %>%
    expand_events()
  
  # Filter the original data frame to exclude rows meeting the condition
  # and bind the expanded rows to the filtered data frame
  df_loc %>%
    filter(!(EventoEnd - EventoStart >= 2 & Distance_m == 0)) %>%
    bind_rows(df_loc_expanded)
}

df_loc_processed <- process_loc_data(df_loc)



### Add pebble PC1
df_ch %<>% bind_cols(data.frame((
  df_ch %>% dplyr::select(a_axis, b_axis, c_axis) %>%
    princomp()
)$scores[, 1:2]) %>%
  set_colnames(c('pebble_PC1', 'pebble_PC2')))

#Add polygon id and recalculate movement distance-----
# Read pool shape data
read_pool_shape_data <- function(dsn, layer) {
  pools_shape <- rgdal::readOGR(dsn = dsn, layer = layer)
  pools_shape_polygon_union <- rgeos::gUnaryUnion(polygons(pools_shape))
  return(list(pools_shape = pools_shape, pools_shape_polygon_union = pools_shape_polygon_union))
}

# Extract pool border coordinates
extract_pool_border <- function(pools_shape_polygon_union) {
  pool_border <- pools_shape_polygon_union@polygons[[1]]@Polygons[[1]]@coords %>%
    as.data.frame() %>%
    distinct(round(., 8))
  return(pool_border)
}

# Extract pool intra coordinates
extract_pool_intra <- function(pools_shape, pool_border) {
  pool_intra <- lapply(pools_shape@polygons, function(x) x@Polygons[[1]]@coords) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    distinct(round(., 8)) %>%
    set_names(c('x', 'y')) %>%
    anti_join(pool_border, by = c('x', 'y'))
  return(pool_intra)
}

# Extract pebble start and end locations
extract_pebble_locations <- function(pebble_coords, coords_for_triang) {
  pebble_locs <- pebble_coords %>%
    pivot_longer(cols = everything(), names_to = c(".value", "coord"), names_sep = "_") %>%
    set_names(c('x', 'y', 'coord')) %>%
    split(.$coord) %>%
    purrr::map(~ anti_join(., as.data.frame(coords_for_triang)))
  return(pebble_locs)
}

# Triangulate domain and calculate distance on the graph
triangulate_and_calculate_distance <- function(coords_for_triang, pool_border) {
  m_border <- cbind(1:nrow(pool_border), c(2:nrow(pool_border), 1))
  data_plsg <- RTriangle::pslg(P = coords_for_triang, S = m_border)
  data_traingl <- RTriangle::triangulate(data_plsg, a = 1e-11)
  
  edge_list <- data_traingl$E
  graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
  
  edge_weights <- geosphere::distHaversine(data_traingl$P[edge_list[, 1], 1:2],
                                           data_traingl$P[edge_list[, 2], 1:2])
  
  triangl_dist <- igraph::distances(graph, weights = edge_weights, algorithm = "dijkstra")
  
  return(triangl_dist)
}

# Calculate graph distances for pebble movements
calculate_graph_dist <- function(pebble_coords, distances, dsn, layer) {
  pool_data <- read_pool_shape_data(dsn, layer)
  pools_shape <- pool_data$pools_shape
  pools_shape_polygon_union <- pool_data$pools_shape_polygon_union
  
  pool_border <- extract_pool_border(pools_shape_polygon_union)
  pool_intra <- extract_pool_intra(pools_shape, pool_border)
  
  coords_for_triang <- unique(rbind(as.matrix(pool_border), as.matrix(pool_intra)))
  pebble_locs <- extract_pebble_locations(pebble_coords, coords_for_triang)
  coords_for_triang <- unique(rbind(coords_for_triang, as.matrix(pebble_locs$Start), as.matrix(pebble_locs$End)))
  
  triangl_dist <- triangulate_and_calculate_distance(coords_for_triang, pool_border)
  
  pebble_ids <- apply(pebble_coords, 1, function(x) {
    start_id <- which(coords_for_triang[, 1] == x[1] & coords_for_triang[, 2] == x[2])[1]
    end_id <- which(coords_for_triang[, 1] == x[3] & coords_for_triang[, 2] == x[4])[1]
    c(start_id, end_id)
  })
  
  graph_dist <- triangl_dist[pebble_ids]
  graph_dist[is.infinite(graph_dist)] <- distances[is.infinite(graph_dist)]
  
  return(graph_dist)
}

# Determine polygon ID for each pebble location
get_polygon_id <- function(pebble_points, dsn, layer) {
  pools_shape <- read_pool_shape_data(dsn, layer)$pools_shape
  
  polygon_id <- apply(pebble_points, 1, function(p) {
    which.max(sapply(pools_shape@polygons, function(polygon) {
      sp::point.in.polygon(p[1], p[2], polygon@Polygons[[1]]@coords[, 1], polygon@Polygons[[1]]@coords[, 2])
    }))
  })
  
  return(polygon_id)
}

distances <- df_loc$Distance_m
df_loc$graph_dist <- calculate_graph_dist(df_loc[, c("X_Start", "Y_Start", "X_End", "Y_End")], distances, "data/initial/shp", "MORFOLOGIA_OKSHP_WGS84")
df_loc$polygon_id <- get_polygon_id(df_loc[, c("X_Start", "Y_Start")], "data/initial/shp", "MORFOLOGIA_OKSHP_WGS84")

rm(list = setdiff(ls(), c('df_loc', 'df_ch')))


### Preprocess weather data
preprocess_weather_data <- function(file_path, start_sheet = 3, end_sheet = 30) {
  # Generate the sheet indices
  sheet_indices <- start_sheet:end_sheet
  
  # Read and preprocess the data from multiple sheets using map_dfr()
  df_weather <- purrr::map_dfr(sheet_indices, function(i) {
    readxl::read_excel(
      path = file_path,
      sheet = i,
      skip = 2,
      col_types = "guess"
    ) %>%
      dplyr::select(2:5) %>%
      set_names(gsub("\\.+", "_", names(.))) %>%
      mutate(EventoStart = i - start_sheet + 1)
  }) %>%
    set_colnames(c("Data","h_CP_cm_","h_CP_m_","Q_CP_m3_s_","EventoStart")) %>%
    dplyr::select(-h_CP_m_)
  
  # Exclude the unwanted column
  #df_weather <- dplyr::select(df_weather, -h_CP_m_)
  
  return(df_weather)
}

df_weather <- preprocess_weather_data("data/initial/Dati idrometeo.xlsx")

saveRDS(
  list(
    df_ch = df_ch,
    df_loc = df_loc,
    df_weather = df_weather
  ),
  file = 'data/preprocessed/preprocessed_observations.RDS',
  compress = F
)

rm(list = ls())
