library(targets)
library(crew)
# devtools::install_github("https://github.com/stormeyseas/macrogrow.git", quiet = T)

tar_option_set(
  packages = c("dplyr", "tidyr", "magrittr", "lubridate", "ggplot2", "arrow", "macrogrow", "units", "geosphere", "qs", "qs2", "stringr"),
  format = "qs", 
  controller = crew_controller_local(workers = 14, seconds_idle = 60),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R")) 
)

base_path <- file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "FRDC-seaweed")
data_path <- file.path(base_path, "data", "processed-spatial-cells")
part_path <- file.path(base_path, "data", "processed-partial-lims")

list(
  # Species data --------------------------------------------------------------------------------------------------
  tar_target(
    name = species_data, 
    command = list(
      targets::tar_read(a_armata_gametophyte, store = "targets_outputs/_species"),
      targets::tar_read(a_taxiformis_gametophyte, store = "targets_outputs/_species"),
      targets::tar_read(ecklonia, store = "targets_outputs/_species"),
      targets::tar_read(ulva, store = "targets_outputs/_species")
    )
  ),
  tar_target(
    name = species_names, 
    command = c("A. armata", "A. taxiformis", "Ecklonia", "Ulva")
  ),
  tar_target(
    name = spec_ni_uptake, 
    command = c("linear", "linear", "MM", "MM")
  ),
  tar_target(
    name = spec_am_uptake, 
    command = c("MM", "MM", "MM", "MM")
  ),

  # State average inputs (doubled) --------------------------------------------------------------------------------
  tar_target(
    name = states,
    command = targets::tar_read(states, store = "targets_outputs/_spatial_cells")
  ),

  tar_target(name = I_input_state_av_file, command = file.path(data_path, "I_input_state_av.parquet"), format = "file"),
  tar_target(
    name = I_input_state_av, 
    command = rbind(I_input_state_av_file %>% read_parquet() %>% filter(state == states), 
                    I_input_state_av_file %>% read_parquet() %>% filter(state == states)), 
    pattern = states
  ),
  
  tar_target(name = T_input_state_av_file, command = file.path(data_path, "T_input_state_av.parquet"), format = "file"),
  tar_target(
    name = T_input_state_av, 
    command = rbind(T_input_state_av_file %>% read_parquet() %>% filter(state == states), 
                    T_input_state_av_file %>% read_parquet() %>% filter(state == states)), 
    pattern = states
  ),
  
  tar_target(name = S_input_state_av_file, command = file.path(data_path, "S_input_state_av.parquet"), format = "file"),
  tar_target(
    name = S_input_state_av, 
    command = rbind(S_input_state_av_file %>% read_parquet() %>% filter(state == states), 
                    S_input_state_av_file %>% read_parquet() %>% filter(state == states)), 
    pattern = states
  ),
  
  tar_target(name = UV_input_state_av_file, command = file.path(data_path, "UV_input_state_av.parquet"), format = "file"),
  tar_target(
    name = UV_input_state_av, 
    command = rbind(UV_input_state_av_file %>% read_parquet() %>% filter(state == states), 
                    UV_input_state_av_file %>% read_parquet() %>% filter(state == states)), 
    pattern = states
  ),
  
  tar_target(name = Ni_input_state_av_file, command = file.path(data_path, "Ni_input_state_av.parquet"), format = "file"),
  tar_target(
    name = Ni_input_state_av, 
    command = rbind(Ni_input_state_av_file %>% read_parquet() %>% filter(state == states), 
                    Ni_input_state_av_file %>% read_parquet() %>% filter(state == states)), 
    pattern = states
  ),
  
  tar_target(name = Am_input_state_av_file, command = file.path(data_path, "Am_input_state_av.parquet"), format = "file"),
  tar_target(
    name = Am_input_state_av, 
    command = rbind(Am_input_state_av_file %>% read_parquet() %>% filter(state == states), 
                    Am_input_state_av_file %>% read_parquet() %>% filter(state == states)), 
    pattern = states
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
  tar_target(name = init_biomass, command = c(0.005, 0.005, 0.005, 0.005)),
  tar_target(
    name = species_ideal_depth,
    command = c(0.5, 0.5, 0.5, 0.5)
  ),

  ## Theoretical locations/scenarios (per state) ------------------------------------------------------------------
  tar_target(
    name = theo_scens,
    command = tibble::tribble(
      ~name,            ~T_mod, ~bathy, ~kW_mod, ~I_mod, ~S_mod, ~Ni_add,  ~Am_add,  ~UV_mod,
      "base",           1,      50,     1,       1,        1,    0,        0,        1,
      "fish_farm",      1,      30,     1,       1,        1,    49.02345, 91.04355, 1,
      "deep_water",     0.9,    100,    1.2,     1,        1,    0,        0,        4,
      "shallow_bay",    1.1,    15,     0.8,     1,        1.1,  49.02345, 0,        0.25,
      "estuary_mouth",  1,      30,     0.8,     1,        0.9,  49.02345, 0,        1
    )
  ),
  tar_target(
    name = theo_scen_names,
    command = unique(theo_scens$name)
  ),
  tar_target(
    name = theo_scen_Nf,
    command = list(
      Nf = biomass_to_Nf(biomass = set_units(init_biomass, "g L-1") %>% set_units("mg m-3") %>% drop_units(),
                         Q_rel = 0.5,
                         spec_params = unlist(species_data),
                         dry = T),
      Q_rel = 0.5
    ),
    pattern = map(init_biomass, species_data),
    iteration = "list"
  ),
  tar_target(
    name = theo_scen_site_params,
    command = c(
      kW = 0.6 * theo_scens$kW_mod[theo_scens$name == theo_scen_names],
      hc = 5,
      hz = theo_scens$bathy[theo_scens$name == theo_scen_names],
      farmA = 50*50
    ),
    pattern = theo_scen_names
  ),
  
  tar_target(
    name = theo_T_input,
    command = T_input_state_av %>% 
      mutate(value = value * theo_scens$T_mod[theo_scens$name == theo_scen_names],
             theo_scen = theo_scen_names,
             state = states),
    pattern = cross(map(T_input_state_av, states), theo_scen_names)
  ),
  tar_target(
    name = theo_S_input,
    command = S_input_state_av %>% 
      mutate(value = value * theo_scens$S_mod[theo_scens$name == theo_scen_names],
             theo_scen = theo_scen_names,
             state = states),
    pattern = cross(map(S_input_state_av, states), theo_scen_names)
  ),
  tar_target(
    name = theo_I_input,
    command = I_input_state_av %>% 
      mutate(value = value * theo_scens$I_mod[theo_scens$name == theo_scen_names],
             theo_scen = theo_scen_names,
             state = states),
    pattern = cross(map(I_input_state_av, states), theo_scen_names)
  ),
  tar_target(
    name = theo_UV_input,
    command = UV_input_state_av %>% 
      mutate(value = value * theo_scens$UV_mod[theo_scens$name == theo_scen_names],
             theo_scen = theo_scen_names,
             state = states),
    pattern = cross(map(UV_input_state_av, states), theo_scen_names)
  ),
  tar_target(
    name = theo_Ni_input,
    command = Ni_input_state_av %>% 
      mutate(value = value + theo_scens$Ni_add[theo_scens$name == theo_scen_names],
             theo_scen = theo_scen_names,
             state = states),
    pattern = cross(map(Ni_input_state_av, states), theo_scen_names)
  ),
  tar_target(
    name = theo_Am_input,
    command = Am_input_state_av %>% 
      mutate(value = value + theo_scens$Am_add[theo_scens$name == theo_scen_names],
             theo_scen = theo_scen_names,
             state = states),
    pattern = cross(map(Am_input_state_av, states), theo_scen_names)
  ),
  
  tar_target(
    name = theo_scen_grow_arma,
    command = do_grow_macroalgae(
      start = start_dates_growth,
      grow_days = growtime_growth[1],
      temperature = theo_T_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[1])],
      salinity = theo_S_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[1])],
      light = theo_I_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[1])],
      velocity = theo_UV_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[1])],
      nitrate = theo_Ni_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[1])],
      ammonium = theo_Am_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[1])],
      ni_uptake = spec_ni_uptake[1],
      am_uptake = spec_am_uptake[1],
      site_params = c(unlist(theo_scen_site_params), d_top = species_ideal_depth[1]),
      spec_params = unlist(species_data[[1]]),
      initials = unlist(theo_scen_Nf[[1]])
    ) %>%
      mutate(
        value = Nf + Ns,
        start_date = format(start_dates_growth, "%Y-%m-%d"),
        start_value = unname(unlist(theo_scen_Nf[[4]])["Nf"]),
        species = "A. armata",
        state = states,
        scen = theo_scen_names,
        date = format(date, "%Y-%m-%d")
      ) %>% 
      select(-c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")),
    pattern = cross(
      map(start_dates_growth, start_ydays_growth),
      map(map(theo_T_input, theo_S_input, theo_UV_input, theo_Ni_input, theo_Am_input), cross(states, map(theo_scen_names, theo_scen_site_params)))
    )
  ),
  
  tar_target(
    name = theo_scen_grow_taxi,
    command = do_grow_macroalgae(
      start = start_dates_growth,
      grow_days = growtime_growth[2],
      temperature = theo_T_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[2])],
      salinity = theo_S_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[2])],
      light = theo_I_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[2])],
      velocity = theo_UV_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[2])],
      nitrate = theo_Ni_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[2])],
      ammonium = theo_Am_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[2])],
      ni_uptake = spec_ni_uptake[2],
      am_uptake = spec_am_uptake[2],
      site_params = c(unlist(theo_scen_site_params), d_top = species_ideal_depth[2]),
      spec_params = unlist(species_data[[2]]),
      initials = unlist(theo_scen_Nf[[2]])
    ) %>%
      mutate(
        value = Nf + Ns,
        start_date = format(start_dates_growth, "%Y-%m-%d"),
        start_value = unname(unlist(theo_scen_Nf[[2]])["Nf"]),
        species = "A. taxiformis",
        state = states,
        scen = theo_scen_names,
        date = format(date, "%Y-%m-%d")
      ) %>% 
      select(-c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")),
    pattern = cross(
      map(start_dates_growth, start_ydays_growth),
      map(map(theo_T_input, theo_S_input, theo_UV_input, theo_Ni_input, theo_Am_input), cross(states, map(theo_scen_names, theo_scen_site_params)))
    )
  ),
  
  tar_target(
    name = theo_scen_grow_eckl,
    command = do_grow_macroalgae(
      start = start_dates_growth,
      grow_days = growtime_growth[3],
      temperature = theo_T_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[3])],
      salinity = theo_S_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[3])],
      light = theo_I_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[3])],
      velocity = theo_UV_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[3])],
      nitrate = theo_Ni_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[3])],
      ammonium = theo_Am_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[3])],
      ni_uptake = spec_ni_uptake[3],
      am_uptake = spec_am_uptake[3],
      site_params = c(unlist(theo_scen_site_params), d_top = species_ideal_depth[3]),
      spec_params = unlist(species_data[[3]]),
      initials = unlist(theo_scen_Nf[[3]])
    ) %>%
      mutate(
        value = Nf + Ns,
        start_date = format(start_dates_growth, "%Y-%m-%d"),
        start_value = unname(unlist(theo_scen_Nf[[3]])["Nf"]),
        species = "Ecklonia",
        state = states,
        scen = theo_scen_names,
        date = format(date, "%Y-%m-%d")
      ) %>% 
      select(-c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")),
    pattern = cross(
      map(start_dates_growth, start_ydays_growth),
      map(map(theo_T_input, theo_S_input, theo_UV_input, theo_Ni_input, theo_Am_input), cross(states, map(theo_scen_names, theo_scen_site_params)))
    )
  ),
  
  tar_target(
    name = theo_scen_grow_ulva,
    command = do_grow_macroalgae(
      start = start_dates_growth,
      grow_days = growtime_growth[4],
      temperature = theo_T_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[4])],
      salinity = theo_S_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[4])],
      light = theo_I_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[4])],
      velocity = theo_UV_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[4])],
      nitrate = theo_Ni_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[4])],
      ammonium = theo_Am_input$value[start_ydays_growth:(start_ydays_growth + growtime_growth[4])],
      ni_uptake = spec_ni_uptake[4],
      am_uptake = spec_am_uptake[4],
      site_params = c(unlist(theo_scen_site_params), d_top = species_ideal_depth[4]),
      spec_params = unlist(species_data[[4]]),
      initials = unlist(theo_scen_Nf[[4]])
    ) %>%
      mutate(
        value = Nf + Ns,
        start_date = format(start_dates_growth, "%Y-%m-%d"),
        start_value = unname(unlist(theo_scen_Nf[[4]])["Nf"]),
        species = "Ulva",
        state = states,
        scen = theo_scen_names,
        date = format(date, "%Y-%m-%d")
      ) %>% 
      select(-c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "up_Ot", "T_lim", "S_lim", "I_top", "I_lim", "u_c")),
    pattern = cross(
      map(start_dates_growth, start_ydays_growth),
      map(map(theo_T_input, theo_S_input, theo_UV_input, theo_Ni_input, theo_Am_input), cross(states, map(theo_scen_names, theo_scen_site_params)))
    )
  )
)



