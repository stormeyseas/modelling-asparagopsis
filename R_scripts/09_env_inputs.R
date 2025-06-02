suppressMessages(suppressWarnings({
  library(targets)
  library(tarchetypes)
  library(crew)
  library(mirai)
  library(tidyr)
  library(conflicted)
  library(qs)
  library(qs2)
  library(magrittr)
  library(stringr)
  library(geotargets)
  library(raster)
  library(sp)
}))
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("dplyr", "macrogrow", "geosphere", "stringr"),
  format = "qs",
  controller = crew::crew_controller_local(workers = 16, seconds_idle = 120),
  workspace_on_error = T,
  garbage_collection = 500
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

# For outfall volume processing
units::remove_unit("ML")
units::install_unit(symbol = "ML", def = "1000000 L")
cell_vol <- units::set_units(4, "km") * units::set_units(4, "km")
cell_vol <- cell_vol * units::set_units(12.5, "m") 
cell_vol <- units::set_units(cell_vol, "ML")

Kd_scale_factor <- 0.000199999994947575
Kd_offset <- 0

units::remove_unit("molN")
units::install_unit(symbol = "molN", def = "14.007 g")

chunk_size <- 150

load(file = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata"))

list(
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
  tar_terra_rast(AquaMODIS_Kd_raster, terra::rast(file.path("D:", "FRDC-Seaweed-Raw-Data", "Aqua_MODIS_KD", "AquaMODIS_Kd_490_data.tif")),
                 memory = "persistent", deployment = "main"),
  tar_terra_rast(BathyTopo_raster, terra::rast(file.path("data_raw", "AusBathyTopo 2024", "bathy_projected.tif")),
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
        cell_no %in% c(885685:885690, 886855:886866, 888025:888042, 889194:889218, 890363:890394, 891534:891570, 
                       892705:892746) ~ "TAS",
        cell_no %in% c(893937:893980, 895114:895154, 896290:896328, 897466:897502, 898642:898676, 899818:899850, 
                       900995:901024, 902171:902198, 903347:903372, 904524:904546, 905701:905719, 
                       906879:906891) ~ "VIC",
        cell_no %in% c(892807:892844, 891633:891668, 890459:890492, 889285:889316, 888111:888140, 886937:886964, 
                       885763:885787, 884589:884610, 883419:883432) ~ "TAS", 
        # Adjusting border between vic and nsw
        cell_no %in% c(837645) ~ "VIC",
        cell_no %in% c(838826:838844, 840003:840019, 841180:841194, 842357:842369, 843534:843544, 844711:844718, 
                       845888:845893, 847065:847067, 848242) ~ "NSW",
        cell_no %in% c(950, 2125, 3300, 4475, 5650, 6825) ~ NA, 
        TRUE ~ state
      )) %>%
      filter(!is.na(state))
  ),
  
  tar_target(
    BARRA_C2_cell_nos, 
    BARRA_C2_cell_coords %>% 
      # slice_sample(n = 2000) %>% 
      select(cell_no) %>% unlist() %>% unname()
  ),
  
  tar_target(
    BARRA_C2_cell_nos_chunked, 
    split(BARRA_C2_cell_nos, ceiling(seq_along(BARRA_C2_cell_nos)/chunk_size))
  ),
  
  
  # Nitrogen prioritisation ---------------------------------------------------------------------------------------
  # Nutrients -----------------------------------------------------------------------------------------------------
  ## Refstation matching ------------------------------------------------------------------------------------------
  # The refstation interpolations serve as a "baseline" for nitrogen concentrations
  tar_target(refstation_coords_file, file.path("data/nitrogen/refstation_locations.parquet"), format = "file"),
  tar_target(refstation_data_files, list.files(file.path("data/nitrogen"), full.names = T) %>% str_subset(".csv"), format = "file"),
  tar_target(
    refstation_coords, 
    refstation_coords_file %>% 
      arrow::read_parquet() %>% 
      filter(!code %in% c("WAU_ningaloo", "WAU_esper", "NSW_bonney"))
  ),
  tar_target(
    refstation_data, 
    refstation_data_files %>% 
      purrr::map_dfr(read.csv) %>% 
      mutate(
        SampleDate = lubridate::yday(as.Date(SampleDate)),
        depth_weight = case_when(depth_band == "5-10m" ~ 0.75,
                                 depth_band == "10-20m" ~ 0.5, 
                                 TRUE ~ 1)
      ) %>% 
      merge(refstation_coords, by = "StationName") %>% 
      group_by(code, StationName, SampleDate, measure) %>% 
      reframe(value = weighted.mean(value, depth_weight, na.rm = T))
  ),
  tar_target(
    N_refstation_data,
    command = {
      purrr::map_dfr(BARRA_C2_cell_nos_chunked[[1]], function(cell_no) {
        refstation_cellpaired <- match_location(
          data_coords = refstation_coords,
          cell_coords = BARRA_C2_cell_coords[BARRA_C2_cell_coords$cell_no == cell_no, ]
        )
        df <- refstation_data %>% 
          filter(code %in% unique(refstation_cellpaired$code)) %>% 
          merge(refstation_cellpaired, by = c("code", "StationName")) %>% 
          mutate(weight = 1/dist,
                 cell_no = cell_no,
                 code = as.factor(code),
                 StationName = as.factor(StationName),
                 data_source = as.factor("refstation"),
                 measure = as.factor(measure))
        df$weight <- df$weight/sum(unique(df$weight))
        df
      })
    },
    pattern = BARRA_C2_cell_nos_chunked
  ),
  
  ## Outfall nutrients --------------------------------------------------------------------------------------------
  tar_target(
    outfall_locations_file,
    file.path("data_raw/national-outfall-database/data-output/outflow_site_locations.parquet"),
    format = "file"
  ),
  tar_target(
    outfall_locations,
    arrow::read_parquet(outfall_locations_file) %>%
      rename(lon = longitude, lat = latitude)
  ),
  tar_target(
    outfall_data_file,
    file.path("data_raw/national-outfall-database/data-output/outflow_results.parquet"),
    format = "file"
  ),
  tar_target(
    N_outfall_data,
    command = {
      purrr::map_dfr(BARRA_C2_cell_nos_chunked[[1]], function(cell_no) {
        cell_coords <- BARRA_C2_cell_coords[BARRA_C2_cell_coords$cell_no == cell_no, ]
        outfall_locs <- outfall_locations
        outfall_locs$dist <- geosphere::distHaversine(
          p = c(cell_coords$longitude, cell_coords$latitude), 
          p2 = as.matrix(cbind(outfall_locs$lon, outfall_locs$lat))
        ) * 10^-3
        outfall_locs <- outfall_locs[outfall_locs$dist <= 48,]
        if (nrow(outfall_locs) == 0) {
          outfall_locs <- rbind(
            outfall_locs, data.frame(state = factor(NA), name = factor(NA), lat = as.numeric(NA), lon = as.numeric(NA), dist = as.numeric(NA)))
        }
        outfall_locs$cell_no <- cell_no
        
        df <- arrow::read_parquet(outfall_data_file) %>%
          dplyr::filter(name %in% outfall_locs$name) %>% 
          dplyr::select(c("name", "indicator", "month", "outfall_conc_mgL", "outflow_vol_ML", "quality")) %>%
          mutate(
            outfall_conc_mgL = units::set_units(outfall_conc_mgL, "mg L-1"),
            outfall_conc_mgm3 = units::set_units(outfall_conc_mgL, "mg m-3"),
            outfall_conc_mgm3 = units::drop_units(outfall_conc_mgm3),
            prop_vol = units::set_units(outflow_vol_ML, "ML") / cell_vol,
            prop_vol = units::drop_units(prop_vol),
            name = as.factor(name)
          ) %>%
          dplyr::filter(outfall_conc_mgm3 < 355000) %>% # ~5000 [mg/L]
          dplyr::select(-c(outfall_conc_mgL, outflow_vol_ML)) %>% 
          merge(outfall_locs, by = "name") %>% 
          mutate(weight = 1/dist,
                 value = outfall_conc_mgm3 * prop_vol)
        df$weight <- df$weight/sum(unique(df$weight))
        df %>% 
          mutate(data_source = as.factor("outfall"))
      })
    },
    pattern = BARRA_C2_cell_nos_chunked
  ),
  
  ## Prioritisation -----------------------------------------------------------------------------------------------
  
  tar_target(
    N_data_combined,
    description = "Combine outfall and refstation data into single dataframe",
    command = {
      df <- rbind(
        N_outfall_data %>% 
          mutate(yday = lubridate::yday(lubridate::make_date("2023", month, sample(10:20, 1)))) %>% 
          dplyr::select(yday, weight, value, data_source, indicator, cell_no) %>% 
          mutate(weight = weight/2), # Outfalls are weighted half
        N_refstation_data %>% 
          rename(yday = SampleDate, indicator = measure) %>% 
          dplyr::select(yday, weight, value, data_source, indicator, cell_no) %>% 
          mutate(indicator = factor(indicator, levels = c("Ammonium_mgm3", "Nitrate_mgm3"), labels = c("ammonia", "nitrate_nitrite")))
      )
      df$weight <- df$weight/sum(unique(df$weight))
      df
    },
    pattern = map(N_outfall_data, N_refstation_data)
  ),
      
  tar_target(
    N_data_prioritised,
    command = {
      purrr::map(BARRA_C2_cell_nos_chunked[[1]], function(cell_no) {
        # Split combined data by measure, get means and sds
        df_1 <- N_data_combined %>% 
          dplyr::filter(cell_no == cell_no & indicator == "nitrate_nitrite")
        df_2 <- N_data_combined %>% 
          dplyr::filter(cell_no == cell_no & indicator == "ammonia")
        means <- c(Ni = weighted.mean(df_1$value, df_1$weight, na.rm = T),
                   Am = weighted.mean(df_2$value, df_2$weight, na.rm = T))
        sds <- c(Ni = sqrt(Hmisc::wtd.var(df_1$value, df_1$weight, na.rm = T)),
                 Am = sqrt(Hmisc::wtd.var(df_2$value, df_2$weight, na.rm = T)))
        
        # Fit curves for each measure
        coefs_1 <- tryCatch(
          expr = {
            fit <- nls(
              formula = as.formula(value ~ a + b * sin((yday * pi + c) / 182.5)),
              data = df_1,
              weights = df_1$weight,
              start = c(a = unname(means["Ni"]), b = unname(means["Ni"]), c = -450),
              lower = c(a = unname(means["Ni"])*0.5, b = unname(means["Ni"])*0.5, c = -Inf),
              upper = c(a = unname(means["Ni"])*1.5, b = unname(sds["Ni"]), c = Inf),
              algorithm = "port"
            )
            coefficients(fit)
          },
          error = function(e) {c(a = unname(means["Ni"]), b = unname(means["Ni"])*0.1, c = -450)}
        ) %>% magrittr::set_names(c("a", "b", "c"))
        
        # df_1 %>% 
        #   filter(cell_no == sample(df_1$cell_no, size = 1)) %>% 
        #   mutate(curve = (coefs_1["a"] + coefs_1["b"] * sin((yday * pi + coefs_1["c"]) / 182.5))) %>% 
        #   ggplot(aes(x = yday, y = value, colour = data_source, size = weight)) +
        #   geom_point() +
        #   geom_line(aes(x = yday, y = curve), inherit.aes = F) +
        #   geom_hline(yintercept = unname(means["Ni"])) +
        #   geom_hline(yintercept = unname(sds["Ni"]), linetype = "dashed")
        
        coefs_2 <- tryCatch(
          expr = {
            fit <- nls(
              formula = as.formula(value ~ a + b * sin((yday * pi + c) / 182.5)),
              data = df_2,
              weights = df_2$weight,
              start = c(a = unname(means["Am"]), b = unname(means["Am"]), c = -450),
              lower = c(a = unname(means["Am"])*0.5, b = unname(means["Am"])*0.5, c = -Inf),
              upper = c(a = unname(means["Am"])*1.5, b = unname(sds["Am"]), c = Inf),
              algorithm = "port"
            )
            coefficients(fit)
          },
          error = function(e) {c(a = unname(means["Am"]), b = unname(means["Am"])*0.1, c = -450)}
        ) %>% magrittr::set_names(c("a", "b", "c"))
        
        # df_2 %>% 
        #   mutate(curve = (coefs_2["a"] + coefs_2["b"] * sin((yday * pi + coefs_2["c"]) / 182.5))) %>% 
        #   ggplot(aes(x = yday, y = value, colour = data_source, size = weight)) +
        #   geom_point() +
        #   geom_line(aes(x = yday, y = curve), inherit.aes = F) +
        #   geom_hline(yintercept = unname(means["Am"])) +
        #   geom_hline(yintercept = unname(sds["Am"]), linetype = "dashed")
        
        dat <- data.frame(
          mean = means, 
          sd = sds,
          a = c(Ni = unname(coefs_1["a"]), Am = unname(coefs_2["a"])),
          b = c(Ni = unname(coefs_1["b"]), Am = unname(coefs_2["b"])),
          c = c(Ni = as.integer(coefs_1["c"]), Am = as.integer(coefs_2["c"]))
        ) %>% 
          mutate(cell_no = cell_no)
        dat$measure <- as.factor(rownames(dat))
        dat %>% tibble::remove_rownames()
      }) %>% purrr::list_rbind()
    },
    pattern = map(BARRA_C2_cell_nos_chunked, N_data_combined)
  ),
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    cell_inputs_chunked,
    command = {
      purrr::map(BARRA_C2_cell_nos_chunked[[1]], function(cell_no) {
        I_input <- terra::extract(BARRA_C2_rsdsdir_raster, cell_no) %>% unlist() %>% unname()
        T_input <- terra::extract(BARRA_C2_ts_raster, cell_no) %>% unlist() %>% unname()
        S_input <- terra::extract(BRAN_salt_raster, cell_no) %>% unlist() %>% unname()
        U_input <- terra::extract(BRAN_u_raster, cell_no) %>% unlist() %>% unname()
        V_input <- terra::extract(BRAN_v_raster, cell_no) %>% unlist() %>% unname()
        Kd_490 <- Kd_offset + Kd_scale_factor * fill_nas_weighted(
          terra::extract(AquaMODIS_Kd_raster, cell_no) %>% unlist() %>% unname()
        )
        hz <- terra::extract(BathyTopo_raster, cell_no) %>% unlist() %>% unname() %>% abs()
        
        N_data <- N_data_prioritised[N_data_prioritised$cell_no == cell_no, ]
        Ni_input <- N_data$a[N_data$measure == "Ni"] + N_data$b[N_data$measure == "Ni"] * 
          sin((1:400 * pi + N_data$c[N_data$measure == "Ni"]) / 182.5)
        Am_input <- N_data$a[N_data$measure == "Am"] + N_data$b[N_data$measure == "Am"] * 
          sin((1:400 * pi + N_data$c[N_data$measure == "Am"]) / 182.5)
        
        df <- data.frame(yday = 1:400,
                   cell_no = cell_no,
                   state = as.factor(BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_no == cell_no])) %>%
          mutate(
            I_input = c(I_input, I_input[1:35]),
            T_input = c(T_input, T_input[1:35]),
            S_input = c(S_input, S_input[1:35]),
            U_input = c(U_input, U_input[1:35]),
            V_input = c(V_input, V_input[1:35]),
            Kd_490 = c(Kd_490, Kd_490[1:35]),
            Ni_input = Ni_input,
            Am_input = Am_input,
            hz = hz
          )
        
        # Add a small amount of random variation
        df$Ni_input <- df$Ni_input + rnorm(nrow(df), 0, 0.1) * N_data$sd[N_data$measure == "Ni"]
        df$Am_input <- df$Am_input + rnorm(nrow(df), 0, 0.1) * N_data$sd[N_data$measure == "Am"]
        df$Ni_input[df$Ni_input < 0] <- 0
        df$Am_input[df$Am_input < 0] <- 0
        df
      }) %>% purrr::list_rbind() %>% 
        # Cleanup and convert to correct units
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
    },
    pattern = map(BARRA_C2_cell_nos_chunked, N_data_prioritised)
  ),
  tar_target(
    state_input_timeseries,
    cell_inputs_chunked %>% 
      group_by(state, yday) %>% 
      reframe(
        I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
        T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
        S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
        UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input),
        Kd_490_mean = meanna(Kd_490), Kd_490_min = minna(Kd_490), Kd_490_max = maxna(Kd_490), Kd_490_sd = sdna(Kd_490),
        Ni_input_mean = meanna(Ni_input), Ni_input_min = minna(Ni_input), Ni_input_max = maxna(Ni_input), Ni_input_sd = sdna(Ni_input),
        Am_input_mean = meanna(Am_input), Am_input_min = minna(Am_input), Am_input_max = maxna(Am_input), Am_input_sd = sdna(Am_input)
      )
  ),
  tar_target(
    cell_input_chunked_gapfilled,
    command = {
      purrr::map(BARRA_C2_cell_nos_chunked[[1]], function(cell_no) {
        cia <- cell_inputs_chunked[cell_inputs_chunked$cell_no == cell_no, ]
        sia <- state_input_timeseries[state_input_timeseries$state == unique(cia$state), ]
        if (any(is.na(cia$I_input))) {cia$I_input <- sia$I_input_mean}
        if (any(is.na(cia$T_input))) {cia$T_input <- sia$T_input_mean}
        if (any(is.na(cia$S_input))) {cia$S_input <- sia$S_input_mean}
        if (any(is.na(cia$Kd_490))) {cia$Kd_490 <- sia$Kd_490_mean}
        if (any(is.na(cia$Ni_input))) {cia$Ni_input <- sia$Ni_input_mean}
        if (any(is.na(cia$Am_input))) {cia$Am_input <- sia$Am_input_mean}
        if (any(is.na(cia$UV_input))) {cia$UV_input <- sia$UV_input_mean}
        cia
      })
    },
    pattern = map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked)
  )
)



