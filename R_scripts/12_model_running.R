suppressMessages(suppressWarnings({
  library(targets)
  library(tarchetypes)
  library(crew)
  library(conflicted)
  library(units)
  library(qs)
  library(qs2)
}))
# devtools::install_github("https://github.com/stormeyseas/macrogrow.git", quiet = T)

tar_option_set(
  packages = c("dplyr", "tidyr", "magrittr", "lubridate", "ggplot2", "arrow", "macrogrow", "units", "geosphere", "qs", "qs2", "stringr"),
  format = "qs", 
  controller = crew_controller_local(workers = 10, seconds_idle = 60),
  workspace_on_error = T,
  error = "abridge"
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R")) 
)

base_path <- file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "modelling-asparagopsis")
spat_path <- file.path(base_path, "data", "processed-spatial-cells")
part_path <- file.path(base_path, "data", "processed-partial-lims")

list(
  # Species data --------------------------------------------------------------------------------------------------
  tar_target(
    name = species_data, 
    command = list(
      targets::tar_read(a_armata_gametophyte, store = "targets_outputs/_species"),
      targets::tar_read(a_taxiformis_gametophyte, store = "targets_outputs/_species"),
    )
  ),
    tar_target(name = species_names, command = c("A. armata", "A. taxiformis")),

  # Cell info -----------------------------------------------------------------------------------------------------
  tar_target(BARRA_C2_cell_nos, targets::tar_read(BARRA_C2_cell_nos, store = "targets_outputs/_spatial_cells")#[seq(1,11677, 3)]
  ),
  tar_target(BARRA_C2_cell_coords, targets::tar_read(BARRA_C2_cell_coords, store = "targets_outputs/_spatial_cells")),
  tar_target(states, targets::tar_read(states, store = "targets_outputs/_spatial_cells")),
  
  # Spatial data files --------------------------------------------------------------------------------------------
  tar_target(
    name = cell_bathy_files, 
    command = file.path(base_path, "data_raw", "AusBathyTopo 2024", "cell_vals") %>% list.files(full.names = T),
    format = "file"
  ),
  tar_target(
    name = I_input_cell_files, 
    command = file.path(base_path, "data_raw", "BARRA-C2", "cell_vals") %>% list.files(full.names = T) %>% str_subset("rsdsdir"),
    format = "file"
  ),
  tar_target(
    name = T_input_cell_files, 
    command = file.path(base_path, "data_raw", "BARRA-C2", "cell_vals") %>% list.files(full.names = T) %>% str_subset("ts"),
    format = "file"
  ),
  tar_target(
    name = S_input_cell_files, 
    command = file.path(base_path, "data_raw", "BARRA2023", "cell_vals") %>% list.files(full.names = T) %>% str_subset("ocean_salt"),
    format = "file"
  ),
  tar_target(
    name = U_input_cell_files, 
    command = file.path(base_path, "data_raw", "BARRA2023", "cell_vals") %>% list.files(full.names = T) %>% str_subset("ocean_u"),
    format = "file"
  ),
  tar_target(
    name = V_input_cell_files, 
    command = file.path(base_path, "data_raw", "BARRA2023", "cell_vals") %>% list.files(full.names = T) %>% str_subset("ocean_v"),
    format = "file"
  ),
  
  ## Cell inputs --------------------------------------------------------------------------------------------------
  tar_target(
    name = cell_inputs,
    command = data.frame(
      yday = 1:730,
      I_input  = c(qread(str_subset(I_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(I_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8)))), 
      T_input  = c(qread(str_subset(T_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(T_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8)))), 
      S_input  = c(qread(str_subset(S_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(S_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8)))), 
      U_input  = c(qread(str_subset(U_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(U_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8)))), 
      V_input  = c(qread(str_subset(V_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(V_input_cell_files,  fixnum(BARRA_C2_cell_nos, digits = 8)))), 
      Ni_input = c(qread(str_subset(Ni_input_cell_files, fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(Ni_input_cell_files, fixnum(BARRA_C2_cell_nos, digits = 8)))), 
      Am_input = c(qread(str_subset(Am_input_cell_files, fixnum(BARRA_C2_cell_nos, digits = 8))), 
                   qread(str_subset(Am_input_cell_files, fixnum(BARRA_C2_cell_nos, digits = 8))))
    ),
    pattern = BARRA_C2_cell_nos
  ),
  
  # Run model -----------------------------------------------------------------------------------------------------
  ## Total growth and remediation (per cell) ----------------------------------------------------------------------
  tar_target(name = months_growth, command = 1:12),
  tar_target(
    name = start_dates_growth,
    command = make_date(year = 2023, month = months_growth, day = 1),
    pattern = months_growth
  ),
  tar_target(
    name = start_ydays_growth,
    command = yday(start_dates_growth),
    pattern = start_dates_growth
  ),

  tar_target(name = growtime_growth, command = c(42, 42, 90, 30)),
  # tar_target(name = init_biomass, command = c(0.005, 0.005, 0.005, 0.005)),
  tar_target(
    name = init_state_growth,
    command = c(
      Nf = biomass_to_Nf(
        biomass = set_units(0.0075, "g L-1") %>%
          set_units("mg m-3") %>%
          drop_units(),
        Q_rel = 0.5,
        spec_params = unlist(species_data),
        dry = T
      ),
      Q_rel = 0.5
    ),
    pattern = species_data,
    iteration = "list"
  ),

  ## Growth -------------------------------------------------------------------------------------------------------
  tar_target(
    name = total_growth_arma,
    command = cbind(
      do_grow_macroalgae(
        start = start_dates_growth,
        grow_days = 42,
        temperature =  cell_inputs$T_input[start_ydays_growth:(start_ydays_growth + 42)],
        salinity =     cell_inputs$S_input[start_ydays_growth:(start_ydays_growth + 42)],
        light =        cell_inputs$I_input[start_ydays_growth:(start_ydays_growth + 42)],
        velocity =     cell_inputs$UV_input[start_ydays_growth:(start_ydays_growth + 42)],
        nitrate =      cell_inputs$Ni_input[start_ydays_growth:(start_ydays_growth + 42)],
        ammonium =     cell_inputs$Am_input[start_ydays_growth:(start_ydays_growth + 42)],
        ni_uptake = "linear",
        am_uptake = "MM",
        site_params = c(hc = 5, farmA = 50 * 50, kW = 0.6, turbulence = NA, d_top = 0.5, hz = abs(cell_bathy$bathy[cell_bathy$cell_ID == BARRA_C2_cell_nos])),
        spec_params = unlist(species_data[[1]]),
        initials = unlist(init_state_growth[[1]])
      ),
      data.frame(
        start_date = rep(format(start_dates_growth, "%Y-%m-%d"), 42+1),
        state = rep(BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_ID == BARRA_C2_cell_nos], 42+1),
        cell_ID = rep(BARRA_C2_cell_nos, 42+1)
      )
    ) %>% 
      mutate(date = format(date, "%Y-%m-%d")) %>% 
      pivot_longer(
        names_to = "output", names_transform = list(output = as.factor),
        values_to = "value", values_transform = list(value = as.numeric),
        cols = c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")
      ) %>% 
      filter(output != "up_Ot"),
    pattern = cross(map(start_dates_growth, start_ydays_growth), map(BARRA_C2_cell_nos, cell_inputs))
  ),

  tar_target(
    name = total_growth_taxi,
    command = cbind(
      do_grow_macroalgae(
        start = start_dates_growth,
        grow_days = 42,
        temperature =  cell_inputs$T_input[start_ydays_growth:(start_ydays_growth + 42)],
        salinity =     cell_inputs$S_input[start_ydays_growth:(start_ydays_growth + 42)],
        light =        cell_inputs$I_input[start_ydays_growth:(start_ydays_growth + 42)],
        velocity =     cell_inputs$UV_input[start_ydays_growth:(start_ydays_growth + 42)],
        nitrate =      cell_inputs$Ni_input[start_ydays_growth:(start_ydays_growth + 42)],
        ammonium =     cell_inputs$Am_input[start_ydays_growth:(start_ydays_growth + 42)],
        ni_uptake = "linear",
        am_uptake = "MM",
        site_params = c(hc = 5, farmA = 50 * 50, kW = 0.6, turbulence = NA, d_top = 0.5, hz = abs(cell_bathy$bathy[cell_bathy$cell_ID == BARRA_C2_cell_nos])),
        spec_params = unlist(species_data[[2]]),
        initials = unlist(init_state_growth[[2]])
      ),
      data.frame(
        start_date = rep(format(start_dates_growth, "%Y-%m-%d"), 42+1),
        state = rep(BARRA_C2_cell_coords$state[BARRA_C2_cell_coords$cell_ID == BARRA_C2_cell_nos], 42+1),
        cell_ID = rep(BARRA_C2_cell_nos, 42+1)
      )
    ) %>% 
      mutate(date = format(date, "%Y-%m-%d")) %>% 
      pivot_longer(
        names_to = "output", names_transform = list(output = as.factor),
        values_to = "value", values_transform = list(value = as.numeric),
        cols = c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")
      ) %>% 
      filter(output != "up_Ot"),
    pattern = cross(map(start_dates_growth, start_ydays_growth), map(BARRA_C2_cell_nos, cell_inputs))
  )
)



