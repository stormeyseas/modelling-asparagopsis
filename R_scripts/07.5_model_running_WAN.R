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
  # library(lubridate)
}))
conflicts_prefer(dplyr::filter(), dplyr::select())

tar_option_set(
  packages = c("dplyr", "macrogrow", "geosphere", "stringr"),
  format = "qs", 
  controller = crew_controller_local(workers = 20, seconds_idle = 120),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

this_state <- "WAN"

# For outfall volume processing
units::remove_unit("ML")
units::install_unit(symbol = "ML", def = "1000000 L")
cell_vol <- units::set_units(4, "km") * units::set_units(4, "km")
cell_vol <- cell_vol * units::set_units(12.5, "m") 
cell_vol <- units::set_units(cell_vol, "ML")

Kd_scale_factor <- 0.000199999994947575
Kd_offset <- 0

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
  tar_terra_rast(BathyTopo_raster, terra::rast(file.path("data_raw", "AusBathyTopo 2024", "bathy_projected.tif")),
                 memory = "persistent", deployment = "main"),
  tar_terra_rast(AquaMODIS_Kd_raster, terra::rast(file.path("D:", "FRDC-Seaweed-Raw-Data", "Aqua_MODIS_KD", "AquaMODIS_Kd_490_data.tif")),
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
      # slice_sample(n = 500) %>% 
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
  ),   
  tar_target(
    outfall_locations,
    arrow::read_parquet(outfall_locations_file) %>%
      filter(name %in% unique(outfall_data$name)) %>%
      mutate(name = droplevels(name)) %>% 
      rename(lon = longitude, lat = latitude)
  ),
  tar_target(
    outfall_sites_cellpaired,
    description = "Determine which outfall sites should influence each cell",
    command = {
      data_coords <- outfall_locations
      cell_coords <- BARRA_C2_cell_coords[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos, ]
      data_coords$dist <- geosphere::distHaversine(
        p = c(cell_coords$longitude, cell_coords$latitude), 
        p2 = as.matrix(cbind(data_coords$lon, data_coords$lat))
      ) * 10^-3
      data_coords <- data_coords[data_coords$dist <= 48,]
      nms <- colnames(data_coords)
      if (nrow(data_coords) == 0) {
        data_coords <- rbind(data_coords, as.data.frame(matrix(NA, nrow = 1, ncol = ncol(data_coords))))
        colnames(data_coords) <- nms
      }
      # data_coords$cell_no <- BARRA_C2_cell_nos
      data_coords
    },
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
        dplyr::filter(indicator == "nitrate_nitrite") %>% 
        mutate(yday = lubridate::yday(lubridate::make_date("2023", month, sample(10:20, 1)))) %>% 
        dplyr::select(yday, weight, value, data_source)
      df_refstat <- Ni_refstation_data %>% 
        rename(yday = SampleDate) %>% 
        dplyr::select(yday, weight, value, data_source)
      df <- rbind(df_outfall, df_refstat)
      
      mean_v <- weighted.mean(df$value, df$weight, na.rm = T)
      sd_v <- sqrt(Hmisc::wtd.var(df$value, df$weight, na.rm = T))
      
      coefs <- tryCatch(
        expr = {
          fit <- nls(
            formula = as.formula(value ~ a + b * sin((yday * pi + c) / 182.5)),
            data = df,
            weights = df$weight,
            start = c(a = mean_v,     b = mean_v*0.1,  c = -450),
            lower = c(a = mean_v*0.5, b = 0,           c = -Inf),
            upper = c(a = mean_v*1.5, b = mean_v*0.95, c = Inf),
            algorithm = "port"
          )
          coefficients(fit)
        },
        error = function(e) {c(a = mean_v, b = mean_v*0.1, c = -450)}
      )
      
      # df$curve <- mean_v + min(sd_v, mean_v) * sin((df$yday * pi - 450) / 182.5)
      # qus_v <- cNORM::weighted.quantile(df$value, probs = c(0.25, 0.5, 0.75), weights = df$weight, type = "Harrell-Davis")
      # df$fit <- coefs['a'] + coefs['b'] * sin((df$yday * pi + coefs['c']) / 182.5)
      # ggplot(df, aes(x = yday, y = value, colour = data_source, size = weight)) + 
      #   geom_point() +
      #   geom_hline(yintercept = mean_v, linetype = "dashed") +
      #   geom_hline(yintercept = qus_v, linetype = "dotted") +
      #   geom_line(aes(y = curve), size = 0.5) +
      #   geom_line(aes(y = fit), size = 1)
      
      list(
        data = df %>% mutate(cell_no = BARRA_C2_cell_nos),
        mean = mean_v, 
        sd = sd_v,
        a = unname(coefs["a"]),
        b = unname(coefs["b"]),
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
        mutate(yday = lubridate::yday(lubridate::make_date("2023", month, sample(10:20, 1)))) %>% 
        dplyr::select(yday, weight, value, data_source)
      df_refstat <- Am_refstation_data %>% 
        rename(yday = SampleDate) %>% 
        dplyr::select(yday, weight, value, data_source)
      df <- rbind(df_outfall, df_refstat)
      
      mean_v <- weighted.mean(df$value, df$weight, na.rm = T)
      sd_v <- sqrt(Hmisc::wtd.var(df$value, df$weight, na.rm = T))
      
      coefs <- tryCatch(
        expr = {
          fit <- nls(
            formula = as.formula(value ~ a + b * sin((yday * pi + c) / 182.5)),
            data = df,
            weights = df$weight,
            start = c(a = mean_v,     b = mean_v*0.1,  c = -450),
            lower = c(a = mean_v*0.5, b = 0,           c = -Inf),
            upper = c(a = mean_v*1.5, b = mean_v*0.95, c = Inf),
            algorithm = "port"
          )
          coefficients(fit)
        },
        error = function(e) {c(a = mean_v, b = mean_v*0.1, c = -450)}
      )
      
      list(
        data = df %>% mutate(cell_no = BARRA_C2_cell_nos),
        mean = mean_v, 
        sd = sd_v,
        a = unname(coefs["a"]),
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
    command = {
      df <- data.frame(
        cell_no = BARRA_C2_cell_nos,
        state = this_state,
        yday = 1:730,
        I_input  = terra::extract(BARRA_C2_rsdsdir_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
        T_input  = terra::extract(BARRA_C2_ts_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
        S_input  = terra::extract(BRAN_salt_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
        U_input  = terra::extract(BRAN_u_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
        V_input  = terra::extract(BRAN_v_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
        Kd_490   = Kd_offset + Kd_scale_factor * fill_nas_weighted(
          terra::extract(AquaMODIS_Kd_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname()
        )
      ) %>%
        # Construct nutrient curves
        mutate(
          Ni_input = Ni_data_prioritised[["a"]] + Ni_data_prioritised[["b"]] * sin((yday * pi + Ni_data_prioritised[["c"]]) / 182.5),
          Ni_input = Ni_input + rnorm(730, 0, 0.1)* Ni_data_prioritised[["sd"]],
          Ni_input = case_when(Ni_input < 0 ~ 0, T ~ Ni_input),
          Am_input = Am_data_prioritised[["a"]] + Am_data_prioritised[["b"]] * sin((yday * pi + Am_data_prioritised[["c"]]) / 182.5),
          Am_input = Am_input + rnorm(730, 0, 0.1)* Am_data_prioritised[["sd"]],
          Am_input = case_when(Am_input < 0 ~ 0, T ~ Am_input),
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
        dplyr::select(-U_input, -V_input)
      
      # Add some variation to nutrient curves
      df$Ni_input <- df$Ni_input + rnorm(nrow(df), 0, 0.05) * Ni_data_prioritised[["sd"]]
      # df$Am_input <- df$Am_input + rnorm(nrow(df), 0, 0.025) * Am_data_prioritised[["sd"]]
    },
    pattern = map(BARRA_C2_cell_nos, Ni_data_prioritised, Am_data_prioritised),
    deployment = "main"
  ),
  tar_target(
    state_input_timeseries,
    cell_input_timeseries %>% 
      group_by(state, yday) %>% 
      reframe(
        I_input = mean(I_input, na.rm = T),
        T_input = mean(T_input, na.rm = T),
        S_input = mean(S_input, na.rm = T),
        UV_input = mean(UV_input, na.rm = T),
        Kd_490 = mean(Kd_490, na.rm = T),
        Ni_input = mean(Ni_input, na.rm = T),
        Am_input = mean(Am_input, na.rm = T)
      )
  ),
  
  # Run model -----------------------------------------------------------------------------------------------------
  ## Partial lims (per cell) --------------------------------------------------------------------------------------
  # tar_target(start_date_PL,  lubridate::make_date(year = 2023, month = 1, day = 1)),
  # tar_target(end_date_PL, lubridate::make_date(year = 2023, month = 12, day = 31)),
  # tar_target(date_range_PL, seq(start_date_PL, end_date_PL, by = 'days')),
  # tar_target(yday_range_PL, lubridate::yday(date_range_PL)),
  # 
  tar_target(init_biomass_PL, 0.005 %>% units::set_units("g L-1") %>% units::set_units("mg m-3") %>% units::drop_units()),
  tar_target(
    init_state_PL,
    c(
      Nf = biomass_to_Nf(biomass = init_biomass_PL, Q_rel = 0.5, spec_params = unlist(species_data), dry = T),
      Q_rel = 0.5
    ),
    pattern = species_data,
    iteration = "list"
  ),
  # 
  ### Temperature -------------------------------------------------------------------------------------------------
  # tar_target(
  #   Tlim_cell_PL,
  #   data.frame(
  #     t = 1:length(date_range_PL),
  #     yday= lubridate::yday(date_range_PL),
  #     temperature = cell_input_timeseries$T_input[yday_range_PL],
  #     T_lim = sapply(
  #       X = cell_input_timeseries$T_input[yday_range_PL],
  #       FUN = T_lim,
  #       spec_params = unlist(species_data)
  #     )) %>%
  #     mutate(
  #       state = this_state,
  #       cell_no = BARRA_C2_cell_nos,
  #       species = species_names
  #     ), 
  #   pattern = cross(map(species_data, species_names), map(BARRA_C2_cell_nos, cell_input_timeseries))
  # ),
  # 
  # tar_target(
  #   Tlim_cell_consec,
  #   command = {
  #     df <- Tlim_cell_PL %>% 
  #       mutate(consec_0.75 = NA, 
  #              consec_0.50 = NA)
  #     df$consec_0.75[1] <- ifelse(df$T_lim[1] >= 0.75, 1, 0)
  #     df$consec_0.50[1]  <- ifelse(df$T_lim[1] >= 0.50, 1, 0)
  #     for (r in 2:nrow(df)) {
  #       df$consec_0.75[r] <- ifelse(df$T_lim[r] >= 0.75, df$consec_0.75[r-1]+1, 0)
  #       df$consec_0.50[r] <- ifelse(df$T_lim[r] >= 0.50, df$consec_0.50[r-1]+1, 0)
  #     }
  #     df %>% 
  #       group_by(state, cell_no, species) %>% 
  #       reframe(consec_0.75 = max(consec_0.75),
  #               consec_0.50 = max(consec_0.50))
  #   },
  #   pattern = Tlim_cell_PL
  # ),
  # 
  ### Irradiance --------------------------------------------------------------------------------------------------
  # tar_target(d_top_PL, (c(1, 2.5, 5))),
  # tar_target(
  #   Ilim_cell_PL, 
  #   data.frame(
  #     t = 1:length(date_range_PL),
  #     yday = lubridate::yday(date_range_PL),
  #     irradiance = cell_input_timeseries$I_input[yday_range_PL],
  #     I_lim = sapply(
  #       X = cell_input_timeseries$I_input[yday_range_PL],
  #       FUN = I_lim,
  #       Nf = unlist(init_state_PL[[1]]["Nf"]),
  #       spec_params = unlist(species_data[[1]]),
  #       site_params = c(d_top = d_top_PL, kW = Secchi_to_Kd(10))
  #     )
  #   ) %>%
  #     mutate(
  #       state = this_state,
  #       cell_no = BARRA_C2_cell_nos,
  #       depth = d_top_PL,
  #     ),     
  #   pattern = cross(d_top_PL, map(BARRA_C2_cell_nos, cell_input_timeseries))
  # ),
  # 
  # tar_target(
  #   Ilim_cell_consec,
  #   command = {
  #     df <- Ilim_cell_PL %>% 
  #       mutate(consec_0.75 = NA, 
  #              consec_0.50 = NA)
  #     df$consec_0.75[1] <- ifelse(df$I_lim[1] >= 0.75, 1, 0)
  #     df$consec_0.50[1]  <- ifelse(df$I_lim[1] >= 0.50, 1, 0)
  #     for (r in 2:nrow(df)) {
  #       df$consec_0.75[r] <- ifelse(df$I_lim[r] >= 0.75, df$consec_0.75[r-1]+1, 0)
  #       df$consec_0.50[r] <- ifelse(df$I_lim[r] >= 0.50, df$consec_0.50[r-1]+1, 0)
  #     }
  #     df %>% 
  #       group_by(state, cell_no, depth) %>% 
  #       reframe(consec_0.75 = max(consec_0.75),
  #               consec_0.50 = max(consec_0.50))
  #   },
  #   pattern = Ilim_cell_PL
  # ),
  # 
  ### Salinity ----------------------------------------------------------------------------------------------------
  # tar_target(
  #   Slim_cell_PL,
  #   command = {
  #     if(any(is.na(cell_input_timeseries$S_input[yday_range_PL]))) {
  #       salinity <- state_input_timeseries$S_input[yday_range_PL]
  #     } else {
  #       salinity <- cell_input_timeseries$S_input[yday_range_PL]
  #     }
  #     data.frame(
  #       t = 1:length(date_range_PL),
  #       yday= lubridate::yday(date_range_PL),
  #       salinity = salinity,
  #       S_lim = sapply(
  #         X = salinity,
  #         FUN = S_lim,
  #         spec_params = unlist(species_data[[1]])
  #       )) %>%
  #       mutate(
  #         state = this_state,
  #         cell_no = BARRA_C2_cell_nos
  #       )
  #   }, 
  #   pattern = map(BARRA_C2_cell_nos, cell_input_timeseries)
  # ),
  # 
  # tar_target(
  #   Slim_cell_consec,
  #   command = {
  #     df <- Slim_cell_PL %>% 
  #       mutate(consec_0.75 = NA, 
  #              consec_0.50 = NA)
  #     df$consec_0.75[1] <- ifelse(df$S_lim[1] >= 0.75, 1, 0)
  #     df$consec_0.50[1]  <- ifelse(df$S_lim[1] >= 0.50, 1, 0)
  #     for (r in 2:nrow(df)) {
  #       df$consec_0.75[r] <- ifelse(df$S_lim[r] >= 0.75, df$consec_0.75[r-1]+1, 0)
  #       df$consec_0.50[r] <- ifelse(df$S_lim[r] >= 0.50, df$consec_0.50[r-1]+1, 0)
  #     }
  #     df %>% 
  #       group_by(state, cell_no) %>% 
  #       reframe(consec_0.75 = max(consec_0.75),
  #               consec_0.50 = max(consec_0.50))
  #   },
  #   pattern = Slim_cell_PL
  # ),
  # 
  ### Nf/Ns loss --------------------------------------------------------------------------------------------------
  # tar_target(
  #   N_loss_PL,
  #   command = {
  #     if(any(is.na(cell_input_timeseries$UV_input[yday_range_PL]))) {
  #       velocity <- state_input_timeseries$UV_input[yday_range_PL]
  #     } else {
  #       velocity <- cell_input_timeseries$UV_input[yday_range_PL]
  #     }
  #     
  #     data.frame(
  #       t = 1:length(date_range_PL),
  #       yday= lubridate::yday(date_range_PL),
  #       velocity = velocity,
  #       N_loss = sapply(
  #         X = velocity,
  #         FUN = loss,
  #         spec_params = unlist(species_data[[1]])
  #       )) %>%
  #       mutate(
  #         state = this_state,
  #         cell_no = BARRA_C2_cell_nos
  #       )
  #   }, 
  #   pattern = map(BARRA_C2_cell_nos, cell_input_timeseries)
  # ),
  
  # Total growth --------------------------------------------------------------------------------------------------
  tar_target(months_growth, 1:12),
  tar_target(start_dates_growth, lubridate::make_date(year = 2023, month = months_growth, day = 24), pattern = months_growth),
  tar_target(start_ydays_growth, lubridate::yday(start_dates_growth), pattern = start_dates_growth),
  tar_target(
    site_params,
    c(
      hc = 5,
      farmA = 50 * 50,
      # kW = Secchi_to_Kd(10), # possibility to add raster of Secchi depth?
      turbulence = NA,
      d_top = 1.5,
      hz = terra::extract(BathyTopo_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname() %>% abs()
    ),
    pattern = BARRA_C2_cell_nos
  ),
  tar_target(
    total_cell_growth,
    command = {
      if (site_params['hz'] > (site_params['d_top'] + site_params['hc'])) {
        if(any(is.na(cell_input_timeseries$S_input[start_ydays_growth:(start_ydays_growth + 42)]))) {
          salinity <- state_input_timeseries$S_input[start_ydays_growth:(start_ydays_growth + 42)]
        } else {
          salinity <- cell_input_timeseries$S_input[start_ydays_growth:(start_ydays_growth + 42)]
        }
        if(any(is.na(cell_input_timeseries$UV_input[start_ydays_growth:(start_ydays_growth + 42)]))) {
          velocity <- state_input_timeseries$UV_input[start_ydays_growth:(start_ydays_growth + 42)]
        } else {
          velocity <- cell_input_timeseries$UV_input[start_ydays_growth:(start_ydays_growth + 42)]
        }
        
        mat <- grow_macroalgae(
          start = start_ydays_growth,
          grow_days = 42,
          temperature =  cell_input_timeseries$T_input[start_ydays_growth:(start_ydays_growth + 42)],
          salinity =     salinity,
          light =        cell_input_timeseries$I_input[start_ydays_growth:(start_ydays_growth + 42)],
          velocity =     velocity,
          nitrate =      cell_input_timeseries$Ni_input[start_ydays_growth:(start_ydays_growth + 42)],
          ammonium =     cell_input_timeseries$Am_input[start_ydays_growth:(start_ydays_growth + 42)],
          ni_uptake = "linear",
          am_uptake = "MM",
          site_params = c(site_params, 'kW' = Kd_490), 
          spec_params = unlist(species_data),
          initials = init_state_PL
        ) 
        # %>% as.data.frame() %>%
        #   mutate(N_loss = Ns_loss + Nf_loss,
        #          TN = (Nf + Ns)*hm*site_params['hc']) %>%
        #   select(-c(Ns_to_Nf, Ns_loss, Nf_loss, Q_int, Nf, Ns, up_Ot, Q_rel))
        
      } else { # blank df
        mat <- matrix(NA, nrow = 1, ncol = 24)
      }
      mat <- cbind(
        mat,
        matrix(data = mat[,'Ns_loss'] + mat[,'Nf_loss'], nrow = nrow(mat), ncol = 1),
        matrix(data = (mat[,'Nf'] + mat[,'Ns']) * mat[,'hm'] * site_params['hc'], nrow = nrow(mat), ncol = 1),
        matrix(data = BARRA_C2_cell_nos, nrow = nrow(mat), ncol = 1)
      )
      colnames(mat) <- c("t", "Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "conc_other", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c", "N_loss", "TN", "cell_no")
      
      mat
      
      # df %>% 
      #   tidyr::pivot_longer(
      #     names_to = "output", names_transform = list(output = as.factor),
      #     values_to = "value", values_transform = list(value = as.numeric),
      #     cols = c("growth_rate", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "T_lim", "S_lim", "I_top", "I_lim", "u_c", "N_loss", "TN")
      #   ) %>% 
      #   mutate(state = factor(state, levels = states))
    },
    pattern = cross(
      map(species_data, init_state_PL), 
      map(start_dates_growth, start_ydays_growth), 
      map(BARRA_C2_cell_nos, cell_input_timeseries, site_params)
    ),
  ),
  tar_target(
    growth_lims,
    get_limiting(total_cell_growth, spec_params = unlist(species_data)),
    pattern = map(total_cell_growth, cross(
      map(species_data, init_state_PL),
      map(start_dates_growth, start_ydays_growth),
      map(BARRA_C2_cell_nos, cell_input_timeseries, site_params)
    ))
  ),
  tar_target(
    total_cell_growth_end_TN,
    total_cell_growth %>% 
      filter(output == "TN" & date == max(date))
  )
)



