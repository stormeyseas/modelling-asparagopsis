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
  controller = crew_controller_local(workers = 12, seconds_idle = 120),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

base_path <- file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "modelling-asparagopsis")
data_path <- file.path(base_path, "data", "processed_cell_data")

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
  tar_target(states, targets::tar_read(states, store = "targets_outputs/_cell_data")),
  tar_target(BARRA_C2_cell_coords, targets::tar_read(BARRA_C2_cell_coords, store = "targets_outputs/_cell_data")),
  tar_target(
    BARRA_C2_cell_nos,
    targets::tar_read(BARRA_C2_cell_nos, store = "targets_outputs/_cell_data") %>% 
      sample(size = 1000, replace = F)
  ), 
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    cell_input_timeseries_files, 
    file.path(data_path, "cell_input_timeseries") %>% 
      list.files(full.names = T) %>% 
      str_subset(fixnum(BARRA_C2_cell_nos, 8)),
    format = "file",
    pattern = BARRA_C2_cell_nos
  ),
  tar_target(
    cell_input_timeseries,
    arrow::read_parquet(cell_input_timeseries_files),
    pattern = cell_input_timeseries_files
  ), 
  
  # Run model -----------------------------------------------------------------------------------------------------
  ## Partial lims (per cell) --------------------------------------------------------------------------------------
  tar_target(start_date_PL,  lubridate::make_date(year = 2023, month = 1, day = 1)),
  tar_target(end_date_PL, lubridate::make_date(year = 2023, month = 12, day = 31)),
  tar_target(date_range_PL, seq(start_date_PL, end_date_PL, by = 'days')),
  tar_target(yday_range_PL, lubridate::yday(date_range_PL)),
  
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
    
  ### Temperature -------------------------------------------------------------------------------------------------
  tar_target(
    Tlim_cell_PL,
    data.frame(
      t = 1:length(date_range_PL),
      date = date_range_PL,
      temperature = cell_input_timeseries$T_input[yday_range_PL],
      T_lim = sapply(
        X = cell_input_timeseries$T_input[yday_range_PL],
        FUN = T_lim,
        spec_params = unlist(species_data)
    )) %>%
      mutate(
        state = BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos],
        cell_no = BARRA_C2_cell_nos,
        species = species_names
      ), 
    pattern = cross(map(species_data, species_names), map(BARRA_C2_cell_nos, cell_input_timeseries))
  ),

  ### Irradiance --------------------------------------------------------------------------------------------------
  tar_target(d_top_PL, (c(1, 2.5, 5, 10))),
  tar_target(
    Ilim_cell_PL, 
    data.frame(
      t = 1:length(date_range_PL),
      date = date_range_PL,
      irradiance = cell_input_timeseries$I_input[yday_range_PL],
      I_lim = sapply(
        X = cell_input_timeseries$I_input[yday_range_PL],
        FUN = I_lim,
        Nf = unlist(init_state_PL["Nf"]),
        spec_params = unlist(species_data),
        site_params = c(d_top = d_top_PL, kW = Kd_from_Secchi(20))
      )
    ) %>%
      mutate(
        state = BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos],
        cell_no = BARRA_C2_cell_nos,
        depth = d_top_PL,
        species = species_names
      ),     
    pattern = cross(map(species_data, species_names, init_state_PL), d_top_PL, map(BARRA_C2_cell_nos, cell_input_timeseries))
  ),
  
  ### Salinity ----------------------------------------------------------------------------------------------------
  tar_target(
    Slim_cell_PL,
    data.frame(
      t = 1:length(date_range_PL),
      date = date_range_PL,
      salinity = cell_input_timeseries$S_input[yday_range_PL],
      S_lim = sapply(
        X = cell_input_timeseries$S_input[yday_range_PL],
        FUN = S_lim,
        spec_params = unlist(species_data)
      )) %>%
      mutate(
        state = BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos],
        cell_no = BARRA_C2_cell_nos,
        species = species_names
      ), 
    pattern = cross(map(species_data, species_names), map(BARRA_C2_cell_nos, cell_input_timeseries))
  ),

  ### Nf/Ns loss --------------------------------------------------------------------------------------------------
  tar_target(
    N_loss_PL,
    data.frame(
      t = 1:length(date_range_PL),
      date = date_range_PL,
      velocity = cell_input_timeseries$UV_input[yday_range_PL],
      N_loss = sapply(
        X = cell_input_timeseries$UV_input[yday_range_PL],
        FUN = loss,
        spec_params = unlist(species_data)
      )) %>%
      mutate(
        state = BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_no == BARRA_C2_cell_nos],
        cell_no = BARRA_C2_cell_nos,
        species = species_names
      ), 
    pattern = cross(map(species_data, species_names), map(BARRA_C2_cell_nos, cell_input_timeseries))
  )
)



