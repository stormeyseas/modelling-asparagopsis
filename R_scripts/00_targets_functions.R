suppressMessages(suppressWarnings({
  library(tidyr)
  library(dplyr)
  library(units)
  library(qs)
  library(qs2)
  library(magrittr)
  library(stringr)
  library(raster)
  library(sp)
  library(units)
  library(streamMetabolizer)
  library(lubridate)
}))

fixnum <- function(n, digits = 4) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}
maxna <- function(x) {max(x, na.rm = T)}
minna <- function(x) {min(x, na.rm = T)}
meanna <- function(x) {mean(x, na.rm = T)}
medna <- function(x) {median(x, na.rm = T)}
sdna <- function(x) {sd(x, na.rm = T)}

adj_params <- function(params, focus_param, factor) {
  params[focus_param] <- params[focus_param] * factor
  return(params)
}

buffer_cells <- function(cell_coords, coastline, dist) {
  cell_coords$in_buffer <- NA
  
  for (i in 1:nrow(cell_coords)) {
    tmp <- geosphere::distHaversine(
      p1 = as.matrix(coastline), 
      p2 = cbind(cell_coords$lon[i], cell_coords$lat[i])
    ) * 10^-3
    if(any(tmp <= dist)) {
      cell_coords$in_buffer[i] <- 1
    }
  }
  return(cell_coords)
}

define_coastline <- function(coast_data) {
  statelist <- list()
  for (i in 1:length(coast_data)) {
    statelist[[i]] <- as.data.frame(coast_data[[i]][[1]])
  }
  coast_data <- statelist %>% 
    bind_rows()
  return(coast_data)
}

positive_end <- function(df) {
  df <- df %>% 
    filter(output %in% c("Nf", "Ns")) %>% 
    mutate(plant_date = min(date)) %>% 
    rename(harv_date = date) %>% 
    group_by(plant_date, harv_date, state, cell_ID, species, batch) %>% 
    reframe(value = sum(value))
  df <- df[df$value == max(df$value), ]
  return(df)
}

get_total_N_long <- function(df) {
  df <- df[df$output %in% c("Nf", "Ns"), ] 
  df <- dplyr::group_by(df, across(all_of(colnames(df)[!colnames(df) %in% c("output", "value")]))) 
  df <- dplyr::reframe(df, total_N = sum(value, na.rm = F))
  return(df)
}

map_cells <- function(cell_df, data_df, try_dist) {
  # Distance of the data cell to every map cell centre
  dist <- geosphere::distHaversine(
    p1 = c(cell_df$lon, cell_df$lat), 
    p2 = as.matrix(cbind(data_df$lon, data_df$lat))
  ) * 10^-3
  dist <- dist[which(dist <= max(try_dist))]
  
  if (any(dist <= try_dist[1])) {
    dat <- data.frame(data_cell_ID = data_df$cell_ID[which(dist <= try_dist[1])],
                      dist = dist[which(dist <= try_dist[1])])
  } else if (any(dist <= try_dist[2])) {
    dat <- data.frame(data_cell_ID = data_df$cell_ID[which(dist <= try_dist[2])],
                      dist = dist[which(dist <= try_dist[2])])
  } else if (any(dist <= try_dist[3])) {
    dat <- data.frame(data_cell_ID = data_df$cell_ID[which(dist <= try_dist[3])],
                      dist = dist[which(dist <= try_dist[3])])
  } else if (any(dist <= try_dist[4])) {
    dat <- data.frame(data_cell_ID = data_df$cell_ID[which(dist <= try_dist[4])],
                      dist = dist[which(dist <= try_dist[4])])
  } else {
    dat <- data.frame(data_cell_ID = NA,
                      dist = NA)
  }
  dat$map_cell_ID <- cell_df$cell_ID
  
  return(dat)
}

match_location <- function(data_coords, cell_coords, choose = 2) {
  data_coords$dist <- geosphere::distHaversine(
    p = c(cell_coords$lon, cell_coords$lat), 
    p2 = as.matrix(cbind(data_coords$lon, data_coords$lat))
  )
  if (!is.na(choose)) {
    data_coords <- arrange(data_coords, dist) %>% slice_head(n = choose)
  }
  return(data_coords)
}

combine_outfall_data <- function(conc_df, vol_df, dist_df) {
  if (nrow(conc_df) != 0 & nrow(vol_df) != 0) {
    df <- merge(conc_df, vol_df, by = c("name", "indicator", "month", "cell_ID"), all = T) %>%
      mutate(value = outfall_conc_mgm3 * prop_vol) %>% 
      dplyr::select(-c(outfall_conc_mgm3, prop_vol)) %>% 
      merge(dist_df, by = c("cell_ID", "name"))
  } else {
    df <- dist_df %>% 
      mutate(indicator = NA, month = NA, value = NA) %>% 
      relocate(dist, .after = value) %>%
      relocate(cell_ID, .before = name)
    colnames(df) <- c("cell_ID", "name", "indicator", "month", "value", "dist")
  }
  return(df)
}

consecutive <- function(df, measure, threshold) {
  colnames(df) <- colnames(df) %>% stringr::str_replace(measure, "result")
  df <- df %>% 
    mutate(thresh_met = case_when(result >= threshold ~ 1, TRUE ~ 0),
           consec = case_when(t == 1 & thresh_met == 1 ~ 1, TRUE ~ 0))
  
  for (i in 2:nrow(df)) {
    if (df$thresh_met[i] == 1 & df$thresh_met[i-1] == 1) {
      df$consec[i] <- df$consec[i-1] + 1
    } 
  }
  return(df)
}

# Function to fill NAs with weighted averages
fill_nas_weighted <- function(x, window_size = 8, decay_factor = 0.65) {
  n <- length(x)
  result <- x
  
  # Find positions of NAs
  na_positions <- which(is.na(x))
  
  for (i in na_positions) {
    # Define window range around missing value
    start_idx <- max(1, i - window_size)
    end_idx <- min(n, i + window_size)
    
    # Get non-NA values in window
    window_values <- x[start_idx:end_idx]
    window_positions <- start_idx:end_idx
    valid_indices <- which(!is.na(window_values))
    
    if (length(valid_indices) > 0) {
      valid_values <- window_values[valid_indices]
      valid_positions <- window_positions[valid_indices]
      
      # Calculate distances from the NA value
      distances <- abs(valid_positions - i)
      
      # Calculate weights based on distances (closer values have higher weights)
      weights <- decay_factor^distances
      
      # Calculate weighted average
      result[i] <- sum(valid_values * weights) / sum(weights)
    }
  }
  
  # Any remaining NAs (e.g., if entire window was NA)
  still_na <- which(is.na(result))
  if (length(still_na) > 0) {
    # Use global mean for any remaining NAs
    result[still_na] <- mean(x, na.rm = TRUE)
  }
  
  return(result)
}

