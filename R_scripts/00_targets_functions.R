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

do_grow_macroalgae <- function(start, grow_days, temperature, salinity, light, velocity, nitrate, ammonium, ni_uptake, am_uptake, site_params, spec_params, initials) {
  name_cols <- c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")

  # if (any(is.na(temperature))) {
  #   df <- cbind(data.frame(date = start), 
  #               as.data.frame(matrix(NA, nrow = 1, ncol = 22)),
  #               data.frame(note = "Error in temperature vector"))
  # } else if (any(is.na(light))) {
  #   df <- cbind(data.frame(date = start), 
  #               as.data.frame(matrix(NA, nrow = 1, ncol = 22)),
  #               data.frame(note = "Error in light vector"))
  # } else 
  if (site_params['hz'] <= (site_params['d_top'] + site_params['hc'])) {
    df <- cbind(date = start, matrix(NA, nrow = 1, ncol = 22, dimnames = list(NULL, name_cols)))
  } else {
    df <- #tryCatch(
    #   expr = {
        grow_macroalgae(
          start = start, 
          grow_days = grow_days,
          temperature = temperature, 
          salinity = salinity, 
          light = light, 
          velocity = velocity,
          nitrate = nitrate, 
          ammonium = ammonium, 
          ni_uptake = ni_uptake, 
          am_uptake = am_uptake,
          site_params = site_params, 
          spec_params = spec_params, 
          initials = initials
        )
      #   }, 
      # warning = function(w){
      #   df <- cbind(date = start, matrix(NA, nrow = 1, ncol = 22, dimnames = list(NULL, name_cols)))
      # }, 
      # error = function(e){
      #   df <- cbind(date = start, matrix(NA, nrow = 1, ncol = 22, dimnames = list(NULL, name_cols)))
      # }
    # )
  }
  # colnames(df) <- name_cols
  return(df)
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

# make_long_base <- function(df, keep_cols, state_levs, species_levs, depth_levs, param_levs, factor_levs) {
#   df$cell_ID <- NULL
#   if (sum(df$growth_rate) == 0) {
#     df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#     colnames(df) <- c("date", "state", "species", "cult_dep", "param", "factor", "measure", "value")
#   } else {
#     cols <- setdiff(colnames(df), keep_cols)
#     df <- tidyr::pivot_longer(df, names_to = "measure", values_to = "value", cols = all_of(cols))
#   }
#   df$date <- as.Date(df$date)
#   df$state <- factor(df$state, levels = state_levs)
#   df$species <- factor(df$species, levels = species_levs)
#   df$cult_dep <- factor(df$cult_dep, levels = depth_levs)
#   df$param <- factor(df$param, levels = param_levs)
#   df$factor <- factor(df$factor, levels = factor_levs)
#   df$measure <- as.factor(df$measure)
#   df$value <- as.numeric(df$value)
#   return(df)
# }

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

try_Slim <- function(date_range_partial_lims, S_input_cell_av, yday_range_partial_lims, species_data_3){
  if (any(is.na(S_input_cell_av$value[yday_range_partial_lims]))) {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      salinity = S_input_cell_av$value[yday_range_partial_lims], 
      S_lim = rep(NA, length(date_range_partial_lims))
      )
  } else {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      salinity = S_input_cell_av$value[yday_range_partial_lims], 
      S_lim = sapply(
        X = S_input_cell_av$value[yday_range_partial_lims], 
        FUN = S_lim,
        spec_params = unlist(species_data_3)
      ))
  }
}

try_loss <- function(date_range_partial_lims, UV_input_cell_av, yday_range_partial_lims, species_data_3){
  if (any(is.na(UV_input_cell_av$value[yday_range_partial_lims]))) {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      UV_velocity = UV_input_cell_av$value[yday_range_partial_lims], 
      loss = rep(NA, length(date_range_partial_lims))
      )
  } else {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      UV_velocity = UV_input_cell_av$value[yday_range_partial_lims], 
      loss = sapply(
        X = UV_input_cell_av$value[yday_range_partial_lims], 
        FUN = loss,
        turbulence = 0,
        spec_params = unlist(species_data_3)
      ))
  }
}

try_Ni_uptake <- function(date_range_partial_lims, Ni_input, species_data_3, spec_ni_uptake_3){
  if (any(is.na(Ni_input))) {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      conc = Ni_input, 
      uptake = rep(NA, length(date_range_partial_lims))
    )
  } else {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      conc = Ni_input, 
      uptake = sapply(
        X = Ni_input, 
        FUN = get_uptake,
        uptake_shape = spec_ni_uptake_3,
        Nform_abbr = "ni",
        spec_params = unlist(species_data_3)
      ))
  }
}

try_Am_uptake <- function(date_range_partial_lims, Am_input_cell_av, yday_range_partial_lims, species_data_3, spec_am_uptake_3) {
  if (any(is.na(Am_input_cell_av$value[yday_range_partial_lims]))) {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      conc = Am_input_cell_av$value[yday_range_partial_lims], 
      uptake = rep(NA, length(date_range_partial_lims))
    )
  } else {
    data.frame(
      t = 1:length(date_range_partial_lims),
      date = date_range_partial_lims,
      conc = Am_input_cell_av$value[yday_range_partial_lims], 
      uptake = sapply(
        X = Am_input_cell_av$value[yday_range_partial_lims], 
        FUN = get_uptake,
        uptake_shape = spec_am_uptake_3,
        Nform_abbr = "am",
        spec_params = unlist(species_data_3)
      ))
  }
}

modify_theo_inputs <- function(input, start, time, mod){
  if (any(is.na(input))) {
    input
  } else {
    input * mod
  }
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

