library(targets)
library(tarchetypes)
library(crew)
library(mirai)
library(tidyr)
library(conflicted)
library(qs, quietly = T)
library(qs2, quietly = T)
library(magrittr)
library(stringr)
library(geotargets)
library(raster)
library(sp, quietly = T)
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("dplyr", "macrogrow", "geosphere", "stringr"),
  format = "qs",
  controller = crew::crew_controller_local(workers = 4, seconds_idle = 120),
  workspace_on_error = T,
  garbage_collection = 500
)

tar_source(
  files = c(file.path("R_scripts", "00_functions.R"))
)

# For outfall volume processing
units::remove_unit("ML")
units::install_unit(symbol = "ML", def = "1000000 L")
cell_vol <- units::set_units(12, "km") * units::set_units(12, "km")
cell_vol <- cell_vol * units::set_units(10, "m") 
cell_vol <- units::set_units(cell_vol, "ML")

Kd_scale_factor <- 0.000199999994947575
Kd_offset <- 0

units::remove_unit("molN")
units::install_unit(symbol = "molN", def = "14.007 g")

chunk_size <- 100

load(file = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2", "BARRA_R2_cell_mask.Rdata"))

list(
  # Big data rasters ----------------------------------------------------------------------------------------------
  tar_terra_rast(BARRA_R2_land_rast, terra::rast(
    file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2", "BARRA_R2_land_rast.tif")
  )), 
  
  tar_terra_rast(
    BARRA_R2_rsdsdir_raster,
    terra::rast(
      file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2", "BARRA_R2_rsdsdir_data.tif")
    ),
    memory = "persistent",
    deployment = "main"
  ), 
  
  tar_terra_rast(
    BARRA_R2_ts_raster,
    terra::rast(
      file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2", "BARRA_R2_ts_data.tif")
    ),
    memory = "persistent",
    deployment = "main"
  ), 
  
  tar_terra_rast(
    BRAN_salt_raster,
    terra::rast(
      file.path("D:", "FRDC-Seaweed-Raw-Data", "BRAN2023", "BRAN_ocean_salt.tif")
    ),
    memory = "persistent",
    deployment = "main"
  ), 
  
  tar_terra_rast(
    BRAN_uv_raster, 
    terra::rast(
      file.path("D:", "FRDC-Seaweed-Raw-Data", "BRAN2023", "BRAN_ocean_uv.tif")
    ),
    memory = "persistent",
    deployment = "main"
  ), 
  
  tar_terra_rast(
    AquaMODIS_Kd_raster,
    terra::rast(
      file.path("D:", "FRDC-Seaweed-Raw-Data", "Aqua_MODIS_KD", "AquaMODIS_Kd_490_data.tif")
    ),
    memory = "persistent",
    deployment = "main"
  ), 
  
  tar_terra_rast(
    BathyTopo_raster,
    terra::rast(
      file.path("data_raw", "AusBathyTopo 2024", "bathy_projected.tif")
    ),
    memory = "persistent",
    deployment = "main"
  ),
  
  # Cell info -----------------------------------------------------------------------------------------------------
  tar_target(
    states_bbox,
    list(
      AUS = c(lonmin = 110, lonmax = 157, latmin = -44.5, latmax = -7.5), # All Australia
      SAU = c(lonmin = 128.94, lonmax = 140.97, latmin = -39, latmax = -31), # All South Australia
      QLD = c(lonmin = 138.05, lonmax = 155, latmin = -28.20, latmax = -8), # All Queensland
      WAS = c(lonmin = 111, lonmax = 128.94, latmin = -36, latmax = -24.35), # South Western Australia
      WAN = c(lonmin = 111, lonmax = 128.94, latmin = -24.35, latmax = -10), # North Western Australia
      TAS = c(lonmin = 143, lonmax = 149.4, latmin = -44.5, latmax = -39.4), # All Tasmania
      VIC = c(lonmin = 140.97, lonmax = 151, latmin = -39.4, latmax = -37.50), # All Victoria
      NSW = c(lonmin = 148, lonmax = 155, latmin = -37.50, latmax = -28.20), # All NSW
      NTE = c(lonmin = 128.94, lonmax = 138.05, latmin = -17, latmax = -10) # All Northern Territory
    )
  ),

  tar_target(
    BARRA_R2_cell_coords, 
    terra::extract(x = BARRA_R2_land_rast, y = 1:(terra::ncell(BARRA_R2_land_rast)), xy = TRUE) %>% 
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
        cell_no %in% c(118170:118173, 119025, 118598, 118599, 118185:118187, 117758:117759, 118620:118632, 119049:119059, 119477:119486, 119905:119913, 118192:118193, 116932:116939, 117358:117367, 117785:117795, 118218:118223) ~ NA, 
        # Adjusting the border between tas and vic
        cell_no %in% c(118174:118184, 118207:118217, 117782:117784) ~ "TAS",
        cell_no %in% c(118626:118632) ~ "VIC",
        # Adjusting border between vic and nsw
        cell_no %in% c(110975) ~ "VIC",
        cell_no %in% c(110976:110982, 111405:111409, 111835:111836) ~ "NSW",
        TRUE ~ state
      )) %>%
      filter(!is.na(state))
  ),
  
  tar_target(
    BARRA_R2_cell_nos, 
    BARRA_R2_cell_coords %>% 
      # slice_sample(n = 250) %>%
      select(cell_no) %>% unlist() %>% unname()
  ),
  
  tar_target(
    BARRA_R2_cell_nos_chunked, 
    split(BARRA_R2_cell_nos, ceiling(seq_along(BARRA_R2_cell_nos)/chunk_size))
  ),
  
  
  # Nitrogen prioritisation ---------------------------------------------------------------------------------------
  # Nutrients -----------------------------------------------------------------------------------------------------
  ## Location matching --------------------------------------------------------------------------------------------
  tar_target(
    refstation_coords_file, 
    file.path("data", "nitrogen", "refstation_locations.parquet"), 
    format = "file"
  ),
  tar_target(
    refstation_coords, 
    refstation_coords_file %>% 
      arrow::read_parquet() %>% 
      filter(!code %in% c("WAU_ningaloo", "WAU_esper", "NSW_bonney"))
  ),
  tar_target(
    outfall_locations_file,
    file.path("data_raw", "national-outfall-database", "data-output", "outflow_site_locations.parquet"),
    format = "file"
  ),
  tar_target(
    outfall_locations,
    arrow::read_parquet(outfall_locations_file) %>%
      rename(lon = longitude, lat = latitude)
  ),
  
  tar_target(
    locations_matched,
    description = "Match each cell two its two nearest refstations and any outfalls within 48 km, weight by distance",
    command = {
      cell_coords <- BARRA_R2_cell_coords %>% 
        filter(cell_no %in% BARRA_R2_cell_nos_chunked[[1]])
      cell_coords <- split(cell_coords, cell_coords$cell_no)
      
      purrr::map2(BARRA_R2_cell_nos_chunked[[1]], cell_coords, function(cell_no, coords) {
        refstation_cellpaired <- match_location(
          data_coords = refstation_coords,
          cell_coords = coords
        ) %>% 
          mutate(
            cell_no = cell_no,
            state = as.factor(coords$state),
            data_source = "refstation"
          ) %>% 
          dplyr::select(-c(lat, lon, StationName)) %>% 
          rename(name = code)
        
        outfall_cellpaired <- outfall_locations %>% 
          mutate(
            dist = geosphere::distHaversine(
              p = c(coords$longitude, coords$latitude), 
              p2 = as.matrix(cbind(outfall_locations$lon, outfall_locations$lat))
              ) * 10^-3
            ) %>% 
          filter(dist <= 48)%>% 
          dplyr::select(-c(state, lat, lon))
        
        outfall_cellpaired %>% 
          mutate(
            cell_no = cell_no,
            state = as.factor(coords$state),
            data_source = "outfall"
          ) %>% 
          rbind(refstation_cellpaired) %>% 
          mutate(
            raw_weight = case_when(data_source == "refstation" ~ 0.5 + 1/dist, 
                                   T ~ 1/dist)
          ) %>%
          group_by(cell_no) %>%
          mutate(
            weight = raw_weight/sum(raw_weight),
            data_source = as.factor(data_source)
            ) %>%
          ungroup() %>%
          dplyr::select(-raw_weight)
      }) %>% 
        bind_rows()
    },
    pattern = BARRA_R2_cell_nos_chunked
  ),
  
  ## Load in data -------------------------------------------------------------------------------------------------
  tar_target(
    refstation_data_files, 
    list.files(file.path("data", "nitrogen"), full.names = T) %>% 
      str_subset(".csv"), 
    format = "file"
  ),
  tar_target(
    outfall_data_file,
    file.path("data_raw", "national-outfall-database", "data-output", "outflow_results.parquet"),
    format = "file"
  ),
  
  tar_target(
    refstation_data, 
    description = "Read in refstation data",
    command = {
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
        reframe(
          value = weighted.mean(value, depth_weight, na.rm = T),
          measure = as.factor(measure),
          StationName = as.factor(StationName),
          code = as.factor(code)) %>% 
        filter(value > 0)
    }
  ),
  tar_target(
    outfall_data,
    description = "Read in outfall data",
    command = {
      arrow::read_parquet(outfall_data_file) %>%
        dplyr::select(c("name", "indicator", "month", "outfall_conc_mgL", "outflow_vol_ML", "quality")) %>%
        mutate(
          outfall_conc_mgL = units::set_units(outfall_conc_mgL, "mg L-1"),
          outfall_conc_mgm3 = units::set_units(outfall_conc_mgL, "mg m-3"),
          outfall_conc_mgm3 = units::drop_units(outfall_conc_mgm3),
          prop_vol = units::set_units(outflow_vol_ML, "ML") / cell_vol,
          prop_vol = units::drop_units(prop_vol),
          name = as.factor(name)
        ) %>%
        dplyr::filter(outfall_conc_mgm3 < 355000 & outfall_conc_mgm3 > 0) %>% # ~5000 [mg/L]
        dplyr::select(-c(outfall_conc_mgL, outflow_vol_ML))
    }
  ),
  
  ## Prioritise and create curves ---------------------------------------------------------------------------------
  tar_target(
    N_data_combined,
    description = "Combine outfall and refstation data into single dataframe, with weights",
    command = {
      refdata <- refstation_data %>% 
        merge(filter(locations_matched, data_source == "refstation"),
              by.x = "code", by.y = "name") %>% 
        rename(yday = SampleDate,
               name = code) %>% 
        mutate(indicator = factor(
          measure,
          levels = levels(refstation_data$measure),
          labels = levels(outfall_data$indicator))
        ) %>%
        relocate(indicator, .after = yday) %>% 
        relocate(cell_no, .before = name) %>% 
        relocate(state, .before = cell_no) %>% 
        relocate(data_source, .after = yday) %>% 
      dplyr::select(-c(measure, StationName, dist))
        
      outdata <- outfall_data %>% 
        merge(filter(locations_matched, data_source == "outfall"),
              by = "name") %>% 
        mutate(
          # Give the outfall reading a day somewhere within the month
          yday = lubridate::yday(lubridate::make_date("2023", month, sample(5:25, 1))),
          value = outfall_conc_mgm3 * prop_vol
        ) %>% 
        relocate(yday, .after = name) %>% 
        relocate(cell_no, .before = name) %>% 
        relocate(state, .before = cell_no) %>% 
        relocate(data_source, .after = yday) %>% 
        relocate(value, .before = weight) %>% 
        dplyr::select(-c(month, quality, outfall_conc_mgm3, prop_vol, dist))
      
      rbind(refdata, outdata) %>% 
        mutate(data_source = as.factor(data_source))
    }
  ),
      
  tar_target(
    N_data_prioritised,
    description = "Create model input curves from the combined N data",
    command = {
      N_data <- filter(N_data_combined, cell_no %in% BARRA_R2_cell_nos_chunked[[1]])
      N_data <- split(N_data, N_data$cell_no)
      
      purrr::map2(BARRA_R2_cell_nos_chunked[[1]], N_data, function(cell_no, N_dat) {
        # Split combined data by measure, get means and sds
        df_Ni <- N_dat %>% dplyr::filter(indicator == "nitrate_nitrite")
        df_Am <- N_dat %>% dplyr::filter(indicator == "ammonia")
        means <- c(Ni = weighted.mean(df_Ni$value, df_Ni$weight, na.rm = T),
                   Am = weighted.mean(df_Am$value, df_Am$weight, na.rm = T))
        sds <- c(Ni = sqrt(Hmisc::wtd.var(df_Ni$value, df_Ni$weight, na.rm = T)),
                 Am = sqrt(Hmisc::wtd.var(df_Am$value, df_Am$weight, na.rm = T)))
        
        # Fit curves for each measure
        coefs_Ni <- tryCatch(
          expr = {
            fit <- nls(
              formula = as.formula(value ~ a + b * sin((yday * pi + c) / 182.5)),
              data = df_Ni,
              weights = df_Ni$weight,
              start = c(a = unname(means["Ni"])*2, b = unname(sds["Ni"]), c = -90),
              lower = c(a = unname(means["Ni"])*0.75, b = -Inf, c = -365),
              upper = c(a = unname(means["Ni"])*2.5, b = Inf, c = 365),
              algorithm = "port"
            )
            coefficients(fit)
          },
          error = function(e) {c(a = unname(means["Ni"]), b = unname(means["Ni"])*0.1, c = -90)}
          ) %>% 
          magrittr::set_names(c("a", "b", "c"))
        
        if (unname(coefs_Ni["b"] > coefs_Ni["a"])) {
          coefs_Ni["b"] <- coefs_Ni["a"] * 0.95
        }
        
        # df_Ni %>%
        #   mutate(curve = (coefs_Ni["a"] + coefs_Ni["b"] * sin((yday * pi + coefs_Ni["c"]) / 182.5))) %>%
        #   ggplot(aes(x = yday, y = (value), colour = data_source, size = weight)) +
        #   geom_point() +
        #   geom_line(aes(x = yday, y = (curve)), inherit.aes = F) +
        #   geom_hline(yintercept = (unname(means["Ni"]))) +
        #   geom_hline(yintercept = (unname(sds["Ni"]+means["Ni"])), linetype = "dashed")
        
        coefs_Am <- tryCatch(
          expr = {
            fit <- nls(
              formula = as.formula(value ~ a + b * sin((yday * pi + c) / 182.5)),
              data = df_Am,
              weights = df_Am$weight,
              start = c(a = unname(means["Am"])*2, b = unname(sds["Am"]), c = -90),
              lower = c(a = unname(means["Am"])*0.75, b = -Inf, c = -365),
              upper = c(a = unname(means["Am"])*2.5, b = Inf, c = 365),
              algorithm = "port"
            )
            coefficients(fit)
          },
          error = function(e) {c(a = unname(means["Am"]), b = unname(means["Am"])*0.1, c = -90)}
        ) %>% magrittr::set_names(c("a", "b", "c"))
        
        if (unname(coefs_Am["b"] > coefs_Am["a"])) {
          coefs_Am["b"] <- coefs_Am["a"] * 0.95
        }
        
        # df_Am %>%
        #   mutate(curve = (coefs_Am["a"] + coefs_Am["b"] * sin((yday * pi + coefs_Am["c"]) / 182.5))) %>%
        #   ggplot(aes(x = yday, y = value, colour = data_source, size = weight)) +
        #   geom_point() +
        #   geom_line(aes(x = yday, y = curve), inherit.aes = F) +
        #   geom_hline(yintercept = unname(means["Am"])) +
        #   geom_hline(yintercept = unname(sds["Am"]), linetype = "dashed")
        
        dat <- data.frame(
          mean = means, 
          sd = sds,
          a = c(Ni = unname(coefs_Ni["a"]), Am = unname(coefs_Am["a"])),
          b = c(Ni = unname(coefs_Ni["b"]), Am = unname(coefs_Am["b"])),
          c = c(Ni = as.integer(coefs_Ni["c"]), Am = as.integer(coefs_Am["c"]))
        ) %>% 
          mutate(cell_no = cell_no)
        dat$measure <- as.factor(rownames(dat))
        dat %>% tibble::remove_rownames()
      }) %>% purrr::list_rbind()
    },
    pattern = BARRA_R2_cell_nos_chunked
  ),
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    gapfilling_inputs_chunked,
    command = {
      coords <- BARRA_R2_cell_coords %>% filter(cell_no %in% BARRA_R2_cell_nos_chunked[[1]])
      coords <- split(coords, coords$cell_no)

      # cell_no <- BARRA_R2_cell_nos_chunked[[1]][1]
      # cell_coords <- coords[[1]]
      purrr::map2(BARRA_R2_cell_nos_chunked[[1]], coords, function(cell_no, cell_coords) {
        # Get initial outputs from vectors (1 cell)
        I_input <- terra::extract(BARRA_R2_rsdsdir_raster, cell_no) %>% unlist() %>% unname()
        T_input <- terra::extract(BARRA_R2_ts_raster, cell_no) %>% unlist() %>% unname()

        S_input <- extract_stepwise_buffer(BRAN_salt_raster, c(cell_coords$lon, cell_coords$lat))
        UV_input <- extract_stepwise_buffer(BRAN_uv_raster, c(cell_coords$lon, cell_coords$lat))
        
        Kd_490 <- extract_stepwise_buffer(AquaMODIS_Kd_raster, c(cell_coords$lon, cell_coords$lat))
        Kd_490 <- Kd_offset + Kd_scale_factor * fill_nas_weighted(Kd_490)
        hz <- terra::extract(BathyTopo_raster, y = as.matrix(cell_coords[1,1:2]), method = "bilinear") %>% unlist() %>% unname() %>% abs()

        df <- data.frame(
          yday = 1:400,
          cell_no = cell_no,
          state = as.factor(BARRA_R2_cell_coords$state[BARRA_R2_cell_coords$cell_no == cell_no])
        ) %>%
          mutate(
            I_input = c(I_input, I_input[1:35]),
            T_input = c(T_input, T_input[1:35]),
            S_input = c(S_input, S_input[1:35]),
            UV_input = c(UV_input, UV_input[1:35]),
            Kd_490 = c(Kd_490, Kd_490[1:35]),
            hz = hz
          )
      }) %>% purrr::list_rbind() %>% 
        # Cleanup and convert to correct units
        mutate(
          I_input = streamMetabolizer::convert_SW_to_PAR(I_input),
          T_input = T_input %>% units::set_units("kelvin") %>% units::set_units("degree_Celsius") %>% units::drop_units()
        )
    },
    pattern = BARRA_R2_cell_nos_chunked,
    deployment = "main"
  ),

  tar_target(
    missing_data,
    command = {
      gapfilling_inputs_chunked %>% 
        group_by(state, cell_no) %>% 
        reframe(
          hz = sum(is.na(hz))/400,
          I_input = sum(is.na(I_input))/400,
          Kd_490 = sum(is.na(Kd_490))/400,
          S_input = sum(is.na(S_input))/400,
          T_input = sum(is.na(T_input))/400,
          UV_input = sum(is.na(UV_input))/400
        )
    },
    pattern = gapfilling_inputs_chunked
  ),

  tar_target(
    cell_inputs_chunked,
    command = {
      N_data <- split(N_data_prioritised, N_data_prioritised$cell_no)
      inputs_data <- split(gapfilling_inputs_chunked, gapfilling_inputs_chunked$cell_no)
      purrr::map2(inputs_data, N_data, function(inputs, N_dat) {
        Ni_input <- N_dat$a[N_dat$measure == "Ni"] + N_dat$b[N_dat$measure == "Ni"] * 
          sin((1:400 * pi + N_dat$c[N_dat$measure == "Ni"]) / 182.5)
        Am_input <- N_dat$a[N_dat$measure == "Am"] + N_dat$b[N_dat$measure == "Am"] * 
          sin((1:400 * pi + N_dat$c[N_dat$measure == "Am"]) / 182.5)
        
        # Add a small amount of random variation
        Ni_input <- Ni_input + rnorm(length(Ni_input), 0, 0.1) * N_dat$sd[N_dat$measure == "Ni"]
        Am_input <- Am_input + rnorm(length(Am_input), 0, 0.1) * N_dat$sd[N_dat$measure == "Am"]
        Ni_input[Ni_input < 0] <- 0
        Am_input[Am_input < 0] <- 0

        inputs %>% 
          mutate(
            Ni_input = Ni_input,
            Am_input = Am_input
          )
      }) %>% purrr::list_rbind()
    },
    pattern = map(gapfilling_inputs_chunked, N_data_prioritised)
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
     ins_list <- split(cell_inputs_chunked, cell_inputs_chunked$cell_no)
     purrr::map(ins_list, function(ins) {
       state <- unique(ins$state) %>% as.character()
       sia <- state_input_timeseries[state_input_timeseries$state == state, ]
       if (any(is.na(ins$I_input))) {ins$I_input <- sia$I_input_mean}
       if (any(is.na(ins$T_input))) {ins$T_input <- sia$T_input_mean}
       if (any(is.na(ins$S_input))) {ins$S_input <- sia$S_input_mean}
       if (any(is.na(ins$Kd_490))) {ins$Kd_490 <- sia$Kd_490_mean}
       if (any(is.na(ins$Ni_input))) {ins$Ni_input <- sia$Ni_input_mean}
       if (any(is.na(ins$Am_input))) {ins$Am_input <- sia$Am_input_mean}
       if (any(is.na(ins$UV_input))) {ins$UV_input <- sia$UV_input_mean}
       ins
     }) %>% purrr::list_rbind()
   },
   pattern = map(BARRA_R2_cell_nos_chunked, cell_inputs_chunked)
  )
)



