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
  error = "abridge",
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

base_path <- file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "modelling-asparagopsis")
data_path <- file.path(base_path, "data", "processed-spatial-cells")

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
      # slice_sample(n = 5726) %>% 
      select(cell_no) %>% unlist() %>% unname()
  ),
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    cell_input_timeseries_files, 
    list.files(data_path, full.names = T) %>% 
      str_subset("cell_input_timeseries") %>% 
      str_subset(states),
    format = "file",
    pattern = states
  ),
  tar_target(
    cell_input_timeseries,
    arrow::read_parquet(cell_input_timeseries_files),
    pattern = cell_input_timeseries_files,
    memory = "persistent"
  ), 
  
  tar_terra_rast(BathyTopo_raster, terra::rast(file.path("data_raw", "AusBathyTopo 2024", "bathy_projected.tif"))),
  tar_target(
    cell_input_static, 
    terra::extract(BathyTopo_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname(),
    pattern = BARRA_C2_cell_nos,
    deployment = "main"
  ),
  
  ## State averages -----------------------------------------------------------------------------------------------
  tar_target(
    state_input_timeseries_files, 
    list.files(data_path, full.names = T) %>% 
      str_subset("state_input_timeseries") %>% 
      str_subset(states),
    format = "file",
    pattern = states
  ),
  tar_target(
    state_input_timeseries,
    arrow::read_parquet(state_input_timeseries_files),
    pattern = state_input_timeseries_files,
    memory = "persistent"
  ), 
  
  # Run model -----------------------------------------------------------------------------------------------------
  ## Partial lims (per cell) --------------------------------------------------------------------------------------
  tar_target(start_date,  lubridate::make_date(year = 2023, month = 1, day = 1)),
  tar_target(end_date, lubridate::make_date(year = 2023, month = 12, day = 31)),
  tar_target(date_range, seq(start_date, end_date, by = 'days')),
  tar_target(yday_range, lubridate::yday(date_range)),
  
  tar_target(init_biomass, 0.005 %>% units::set_units("g L-1") %>% units::set_units("mg m-3") %>% units::drop_units()),

  tar_target(
    name = total_growth_arma,
    command = {
      site_params <- list(
        hz = terra::extract(BathyTopo_raster, BARRA_C2_cell_nos) %>% unlist() %>% unname()
        hc = 2.5,
        d_top = culture_depths,
        kW = Kd_from_Secchi(20)
      )
      if (site_params['hz'] <= (site_params['d_top'] + site_params['hc'])) {
        cbind(date = as.Date(parse_date_time(x = paste0("2022-", starts), orders = "yj")), 
              matrix(NA, nrow = 1, ncol = 22, dimnames = list(NULL, c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c"))))
      } else {
        grow_macroalgae(
          start = as.Date(parse_date_time(x = paste0("2022-", starts), orders = "yj")), 
          grow_days = 42,
          temperature = cell_input_timeseries$T_input[cell_input_timeseries$yday %in% starts:(starts + 42)], 
          salinity = cell_input_timeseries$S_input[cell_input_timeseries$yday %in% starts:(starts + 42)], 
          light = cell_input_timeseries$I_input[cell_input_timeseries$yday %in% starts:(starts + 42)], 
          velocity = cell_input_timeseries$U_input[cell_input_timeseries$yday %in% starts:(starts + 42)],
          nitrate = cell_input_timeseries$Ni_input[cell_input_timeseries$yday %in% starts:(starts + 42)], 
          ammonium = cell_input_timeseries$Am_input[cell_input_timeseries$yday %in% starts:(starts + 42)], 
          ni_uptake = "linear", 
          am_uptake = "MM",
          site_params = site_params,
          spec_params = unlist(species_data),
          initials = list(
            Nf = biomass_to_Nf(biomass = init_biomass, Q_rel = 0.5, spec_params = unlist(species_data), dry = T),
            Q_rel = 0.5
          )
        )
      }
      },
    pattern = cross(BARRA_C2_cell_nos, starts, species_data)
  ),
  
)



