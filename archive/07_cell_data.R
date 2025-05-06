suppressMessages(suppressWarnings({
  library(targets)
  library(tarchetypes)
  library(crew)
  library(tidyr)
  library(conflicted)
  library(qs)
  library(qs2)
  library(magrittr)
  library(stringr)
  library(geotargets)
  library(raster)
  library(sp)
  library(lubridate)
}))
conflicts_prefer(dplyr::filter(), dplyr::select())

tar_option_set(
  packages = c("dplyr", "macrogrow", "geosphere", "stringr"),
  format = "qs", 
  controller = crew_controller_local(workers = 10, seconds_idle = 120),
  error = "continue",
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

base_path <- file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "modelling-asparagopsis")
data_path <- file.path(base_path, "data", "processed-spatial-cells")

this_stat <- "TAS"

# For outfall volume processing
units::remove_unit("ML")
units::install_unit(symbol = "ML", def = "1000000 L")
cell_vol <- units::set_units(12, "km") * units::set_units(12, "km")
cell_vol <- cell_vol * units::set_units(10, "m") 
cell_vol <- units::set_units(cell_vol, "ML")

load(file = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata"))

list(
  # Species data --------------------------------------------------------------------------------------------------
  tar_target(
    species_data, 
    list(
      targets::tar_read(a_armata_gametophyte, store = "targets_outputs/_species"),
      targets::tar_read(a_taxiformis_gametophyte, store = "targets_outputs/_species")
    )
  ),
    tar_target(
    species_names, 
    c("A. armata", "A. taxiformis")
  ),
  
  # Big data rasters ----------------------------------------------------------------------------------------------
  tar_terra_rast(BARRA_C2_land_rast, terra::rast(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_land_rast.tif"))),
  
  # These have already been projected (but not masked)
  tar_terra_rast(BARRA_C2_rsdsdir_raster, terra::rast(file.path("data_raw", "BARRA-C2", "BARRA_C2_rsdsdir_data.tif")), 
                 memory = "persistent", deployment = "main"),
  tar_terra_rast(BARRA_C2_ts_raster, terra::rast(file.path("data_raw", "BARRA-C2", "BARRA_C2_ts_data.tif")), 
                 memory = "persistent", deployment = "main"),
  tar_terra_rast(BRAN_salt_raster, terra::rast(file.path("data_raw", "BRAN2023", "BRAN_salt_data.tif")), 
                 memory = "persistent", deployment = "main"),
  tar_terra_rast(BRAN_u_raster, terra::rast(file.path("data_raw", "BRAN2023", "BRAN_u_data.tif")), 
                 memory = "persistent", deployment = "main"),
  tar_terra_rast(BRAN_v_raster, terra::rast(file.path("data_raw", "BRAN2023", "BRAN_v_data.tif")), 
                 memory = "persistent", deployment = "main"),
  
  # Cell info -----------------------------------------------------------------------------------------------------
  tar_target(
    states_bbox,
    list(
      AUS = c(lonmin = 110, lonmax = 157, latmin = -44.5, latmax = -7.5), # All Australia
      SAU = c(lonmin = 128.94, lonmax = 140.97, latmin = -39, latmax = -31), # All South Australia
      QLD = c(lonmin = 138.05, lonmax = 155, latmin = -28.20, latmax = -8), # All Queensland
      # WAU = c(lonmin = 111, lonmax = 128.94, latmin = -36, latmax = -10), # All Western Australia (not in use!)
      WAS = c(lonmin = 111, lonmax = 128.94, latmin = -36, latmax = -24.35), # South Western Australia
      WAN = c(lonmin = 111, lonmax = 128.94, latmin = -24.35, latmax = -10), # North Western Australia
      TAS = c(lonmin = 143, lonmax = 149.4, latmin = -44.5, latmax = -39.4), # All Tasmania
      VIC = c(lonmin = 140.97, lonmax = 151, latmin = -39.4, latmax = -37.50), # All Victoria
      NSW = c(lonmin = 148, lonmax = 155, latmin = -37.50, latmax = -28.20), # All NSW
      NTE = c(lonmin = 128.94, lonmax = 138.05, latmin = -17, latmax = -10) # All Northern Territory
    )
  ),
  tar_target(states, names(states_bbox)[names(states_bbox) != "AUS"]),
  
  tar_target(
    BARRA_C2_cell_coords, 
    terra::extract(x = BARRA_C2_land_rast, y = 1:(terra::ncell(BARRA_C2_land_rast)), xy = TRUE) %>% 
      rename(longitude = x, latitude = y) %>% 
      mutate(cell_no = 1:n()) %>% 
      filter(layer < 95) %>% 
      mutate(state = case_when(
        longitude > states_bbox[["SAU"]]["lonmin"] & longitude <= states_bbox[["SAU"]]["lonmax"] & 
          latitude > states_bbox[["SAU"]]["latmin"] & latitude <= states_bbox[["SAU"]]["latmax"]  ~ "SAU",
        longitude > states_bbox[["QLD"]]["lonmin"] & longitude <= states_bbox[["QLD"]]["lonmax"] & 
          latitude > states_bbox[["QLD"]]["latmin"] & latitude <= states_bbox[["QLD"]]["latmax"]  ~ "QLD",
        longitude > states_bbox[["WAS"]]["lonmin"] & longitude <= states_bbox[["WAS"]]["lonmax"] & 
          latitude > states_bbox[["WAS"]]["latmin"] & latitude <= states_bbox[["WAS"]]["latmax"]  ~ "WAS",
        longitude > states_bbox[["WAN"]]["lonmin"] & longitude <= states_bbox[["WAN"]]["lonmax"] & 
          latitude > states_bbox[["WAN"]]["latmin"] & latitude <= states_bbox[["WAN"]]["latmax"]  ~ "WAN",
        longitude > states_bbox[["VIC"]]["lonmin"] & longitude <= states_bbox[["VIC"]]["lonmax"] & 
          latitude > states_bbox[["VIC"]]["latmin"] & latitude <= states_bbox[["VIC"]]["latmax"]  ~ "VIC",
        longitude > states_bbox[["NSW"]]["lonmin"] & longitude <= states_bbox[["NSW"]]["lonmax"] & 
          latitude > states_bbox[["NSW"]]["latmin"] & latitude <= states_bbox[["NSW"]]["latmax"]  ~ "NSW",
        longitude > states_bbox[["NTE"]]["lonmin"] & longitude <= states_bbox[["NTE"]]["lonmax"] & 
          latitude > states_bbox[["NTE"]]["latmin"] & latitude <= states_bbox[["NTE"]]["latmax"]  ~ "NTE",
        longitude > states_bbox[["TAS"]]["lonmin"] & longitude <= states_bbox[["TAS"]]["lonmax"] & 
          latitude > states_bbox[["TAS"]]["latmin"] & latitude <= states_bbox[["TAS"]]["latmax"]  ~ "TAS",
        TRUE ~ NA
      )) %>% 
      mutate(state = case_when(
        # Adjusting the border between tas and vic
        cell_no %in% c(885685:885690, 886855:886866, 888025:888042, 889194:889218, 890363:890394, 891534:891570, 892705:892746) ~ "TAS",
        cell_no %in% c(893937:893980, 895114:895154, 896290:896328, 897466:897502, 898642:898676, 899818:899850, 900995:901024, 902171:902198, 903347:903372, 904524:904546, 905701:905719, 906879:906891) ~ "VIC",
        cell_no %in% c(892807:892844, 891633:891668, 890459:890492, 889285:889316, 888111:888140, 886937:886964, 885763:885787, 884589:884610, 883419:883432) ~ "TAS", 
        # Adjusting border between vic and nsw
        cell_no %in% c(837645) ~ "VIC",
        cell_no %in% c(838826:838844, 840003:840019, 841180:841194, 842357:842369, 843534:843544, 844711:844718, 845888:845893, 847065:847067, 848242) ~ "NSW",
        cell_no %in% c(950, 2125, 3300, 4475, 5650, 6825) ~ NA, 
        TRUE ~ state
      )) %>%
      filter(!is.na(state))
  ),
  
  tar_target(
    BARRA_C2_cell_nos, 
    BARRA_C2_cell_coords %>% 
      filter(state == this_state) %>%
      select(cell_no) %>% unlist() %>% unname()
  ),
  
  # Nitrogen prioritisation ---------------------------------------------------------------------------------------
  # Nutrients -----------------------------------------------------------------------------------------------------
  ## Refstation matching ------------------------------------------------------------------------------------------
  # The refstation interpolations serve as a "baseline" for nitrogen concentrations
  tar_target(refstation_coords_file, file.path("data/nitrogen/refstation_locations.parquet"), format = "file"),
  tar_target(
    refstation_coords, 
    refstation_coords_file %>% 
      arrow::read_parquet() %>% 
      filter(!code %in% c("WAU_ningaloo", "WAU_esper", "NSW_bonney"))
  ),
  tar_target(refstation_data_files, list.files(file.path("data/nitrogen"), full.names = T) %>% str_subset(".csv"), format = "file"),
  tar_target(num_refstation_data_files, 1:length(refstation_data_files)),
  tar_target(
    refstation_data, 
    refstation_data_files[num_refstation_data_files] %>% 
      read.csv() %>% 
      mutate(
        SampleDate = lubridate::yday(as.Date(SampleDate)),
        depth_weight = case_when(depth_band == "5-10m" ~ 0.75,
                                 depth_band == "10-20m" ~ 0.5, 
                                 TRUE ~ 1)
      ) %>% 
      merge(refstation_coords, by = "StationName") %>% 
      group_by(code, StationName, SampleDate, measure) %>% 
      reframe(value = weighted.mean(value, depth_weight, na.rm = T)),
    pattern = num_refstation_data_files
  ),
  tar_target(
    refstation_cell_match,
    description = "Match each cell to its nearest 2 refstations",
    match_location(
      data_coords = refstation_coords,
      cell_coords = BARRA_C2_cell_coords[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos, ]
    ),
    pattern = BARRA_C2_cell_nos
  ),
  tar_target(
    Ni_refstation_data,
    description = "Load distance-weighted refstation nitrate data applicable to each cell",
    command = {
      df1 <- refstation_data %>% 
        filter(code == refstation_cell_match$code[1] & measure == "Nitrate_mgm3") %>% 
        mutate(weight = 1/refstation_cell_match$dist[1])
      df2 <- refstation_data %>% 
        filter(code == refstation_cell_match$code[2] & measure == "Nitrate_mgm3") %>% 
        mutate(weight = 1/refstation_cell_match$dist[2])
      df <- rbind(df1, df2) %>% mutate(cell_no = BARRA_C2_cell_nos)
      df$weight <- df$weight/sum(unique(df$weight))
      df %>% 
        mutate(data_source = "refstation")
    },
    pattern = map(refstation_cell_match, BARRA_C2_cell_nos)
  ),
  tar_target(
    Am_refstation_data,
    description = "Load distance-weighted refstation ammonium data applicable to each cell",
    command = {
      df1 <- refstation_data %>% 
        filter(code == refstation_cell_match$code[1] & measure == "Ammonium_mgm3") %>% 
        mutate(weight = 1/refstation_cell_match$dist[1])
      df2 <- refstation_data %>% 
        filter(code == refstation_cell_match$code[2] & measure == "Ammonium_mgm3") %>% 
        mutate(weight = 1/refstation_cell_match$dist[2])
      df <- rbind(df1, df2) %>% mutate(cell_no = BARRA_C2_cell_nos)
      df$weight <- df$weight/sum(unique(df$weight))
      df %>% 
        mutate(data_source = "refstation")
    },
    pattern = map(refstation_cell_match, BARRA_C2_cell_nos)
  ),
  
  ## Outfall nutrients --------------------------------------------------------------------------------------------
  tar_target(
    outfall_data_file,
    file.path("data_raw/national-outfall-database/data-output/outflow_results.parquet"),
    format = "file"
  ),
  tar_target(
    outfall_locations_file,
    file.path("data_raw/national-outfall-database/data-output/outflow_site_locations.parquet"),
    format = "file"
  ),
  tar_target(
    outfall_data,
    arrow::read_parquet(outfall_data_file) %>%
      dplyr::select(c("name", "indicator", "month", "outfall_conc_mgL", "outflow_vol_ML", "quality")) %>%
      mutate(
        outfall_conc_mgL = units::set_units(outfall_conc_mgL, "mg L-1"),
        outfall_conc_mgm3 = units::set_units(outfall_conc_mgL, "mg m-3"),
        outfall_conc_mgm3 = units::drop_units(outfall_conc_mgm3),
        prop_vol = units::set_units(outflow_vol_ML, "ML") / cell_vol,
        prop_vol = units::drop_units(prop_vol)
      ) %>%
      filter(outfall_conc_mgm3 < 355000) %>% # ~5000 [mg/L]
      dplyr::select(-c(outfall_conc_mgL, outflow_vol_ML))
  ),   tar_target(
    outfall_locations,
    arrow::read_parquet(outfall_locations_file) %>%
      filter(name %in% unique(outfall_data$name)) %>%
      mutate(droplevels(name)) %>% 
      rename(lon = longitude, lat = latitude) %>% 
      mutate(cell_no = as.character(name))
  ),
  tar_target(
    outfall_sites_cellpaired,
    description = "Determine which outfall sites should influence each cell",
    match_location(
      data_coords = outfall_locations,
      cell_coords = BARRA_C2_cell_coords[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos, ],
      choose = NA
    ) %>% 
      filter(dist <= 48*10^3),
    pattern = BARRA_C2_cell_nos
  ),
  tar_target(
    outfall_mapped_data,
    command = {
      df <- outfall_data %>% 
        filter(name %in% outfall_sites_cellpaired$name) %>% 
        merge(outfall_sites_cellpaired, by = "name") %>% 
        mutate(weight = 1/dist,
               value = outfall_conc_mgm3 * prop_vol)
      df$weight <- df$weight/sum(unique(df$weight))
      df %>% 
        mutate(data_source = "outfall")
    },
    pattern = outfall_sites_cellpaired
  ),
  
  ## Prioritisation -----------------------------------------------------------------------------------------------
  tar_target(
    Ni_data_prioritised,
    command = {
      df_outfall <- outfall_mapped_data %>% 
        filter(indicator == "nitrate_nitrite") %>% 
        mutate(yday = lubridate::yday(lubridate::make_date("2023", month, 15)),
               weight = weight * 0.5) %>% 
        dplyr::select(yday, weight, value, data_source)
      df_refstat <- Ni_refstation_data %>% 
        rename(yday = SampleDate) %>% 
        dplyr::select(yday, weight, value, data_source)
      df <- rbind(df_outfall, df_refstat)
      # ggplot(df, aes(x = yday, y = value, colour = data_source)) + geom_point()
      
      mean_v <- weighted.mean(df$value, df$weight, na.rm = T)
      coefs <- tryCatch(
        expr = {
          fit <- nls(
            formula = as.formula(value ~ mean_v + b * sin((yday * pi + c) / 182.5)),
            data = df,
            weights = df$weight,
            start = c(b = mean_v * 0.1, c = df$yday[which(df$value == max(df$value))]),
            lower = c(b = 0, c = -365),
            upper = c(b = mean_v, c = 365),
            algorithm = "port"
          )
          coefficients(fit)
        },
        error = function(e) {c(b = mean_v, c = df$yday[which(df$value == max(df$value))])}
      )
      
      list(
        data = df %>% mutate(cell_no = BARRA_C2_cell_nos),
        a = mean_v,
        b = coefs["b"],
        c = as.integer(coefs["c"])
      )
    },
    pattern = map(BARRA_C2_cell_nos, outfall_mapped_data, Ni_refstation_data),
    iteration = "list"
  ),
  
  tar_target(
    Am_data_prioritised,
    command = {
      df_outfall <- outfall_mapped_data %>% 
        filter(indicator == "ammonia") %>% 
        mutate(yday = lubridate::yday(lubridate::make_date("2023", month, 15)),
               weight = weight * 0.5) %>% 
        dplyr::select(yday, weight, value, data_source)
      df_refstat <- Am_refstation_data %>% 
        rename(yday = SampleDate) %>% 
        dplyr::select(yday, weight, value, data_source)
      df <- rbind(df_outfall, df_refstat)
      # ggplot(df, aes(x = yday, y = value, colour = data_source)) + geom_point()
      
      mean_v <- weighted.mean(df$value, df$weight, na.rm = T)
      coefs <- tryCatch(
        expr = {
          fit <- nls(
            formula = as.formula(value ~ mean_v + b * sin((yday * pi + c) / 182.5)),
            data = df,
            weights = df$weight,
            start = c(b = mean_v * 0.1, c = df$yday[which(df$value == max(df$value))]),
            lower = c(b = 0, c = -365),
            upper = c(b = mean_v*0.9, c = 365),
            algorithm = "port"
          )
          coefficients(fit)
        },
        error = function(e) {c(b = mean_v, c = df$yday[which(df$value == max(df$value))])}
      )
      
      list(
        data = df %>% mutate(cell_no = BARRA_C2_cell_nos),
        a = mean_v,
        b = unname(coefs["b"]),
        c = as.integer(coefs["c"])
      )
    },
    pattern = map(BARRA_C2_cell_nos, outfall_mapped_data, Am_refstation_data),
    iteration = "list"
  ),
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    cell_input_timeseries,
    data.frame(
      cell_no = BARRA_C2_cell_nos,
      state = this_state,
      yday = 1:730,
      I_input  = terra::extract(BARRA_C2_rsdsdir_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
      T_input  = terra::extract(BARRA_C2_ts_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
      S_input  = terra::extract(BRAN_salt_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
      U_input  = terra::extract(BRAN_u_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
      V_input  = terra::extract(BRAN_v_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname()
      ) %>% 
      mutate(
        Ni_input = Ni_data_prioritised[["a"]] + Ni_data_prioritised[["b"]] * sin((yday * pi + Ni_data_prioritised[["c"]]) / 182.5),
        Am_input = Am_data_prioritised[["a"]] + Am_data_prioritised[["b"]] * sin((yday * pi + Am_data_prioritised[["c"]]) / 182.5),
      ) %>%
      mutate(
        I_input = streamMetabolizer::convert_SW_to_PAR(I_input),
        T_input = T_input %>% units::set_units("kelvin") %>% units::set_units("degree_Celsius") %>% units::drop_units(),
        UV_input = case_when(
          is.na(U_input) & is.na(V_input) ~ NA, 
          is.na(U_input) ~ sqrt(V_input^2 + V_input^2),
          is.na(V_input) ~ sqrt(U_input^2 + U_input^2),
          T ~ sqrt(U_input^2 + V_input^2)
          )
      ) %>% 
      select(-U_input, -V_input),
    pattern = map(BARRA_C2_cell_nos, Ni_data_prioritised, Am_data_prioritised),
    deployment = "main"
  )
)



