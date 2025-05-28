# nolint start

# Setup -----------------------------------------------------------------------------------------------------------
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
}))
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("dplyr", "macrogrow", "geosphere", "stringr"),
  format = "qs", 
  controller = crew_controller_local(workers = 16, seconds_idle = 120),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

units::remove_unit("molN")
units::install_unit(symbol = "molN", def = "14.007 g")
Ni_add <- 3.5 %>% units::set_units("umolN L-1") %>% units::set_units("mg m-3") %>% units::drop_units()
Am_add <- 6.5 %>% units::set_units("umolN L-1") %>% units::set_units("mg m-3") %>% units::drop_units()

chunk_size <- 100

list(
  # Species data --------------------------------------------------------------------------------------------------
  tar_target(
    species_data, 
    list(targets::tar_read(a_armata_gametophyte, store = "targets_outputs/_species"),
         targets::tar_read(a_taxiformis_gametophyte, store = "targets_outputs/_species"))
  ),
  tar_target(species_names, c("armata", "taxiformis")),

  tar_target(BARRA_C2_cell_coords_file, "data/processed_env_inputs/BARRA_C2_cell_coords.parquet", format = "file"),
  tar_target(BARRA_C2_cell_coords, arrow::read_parquet(BARRA_C2_cell_coords_file)),
  tar_target(states, unique(BARRA_C2_cell_coords$state)),
  
  tar_target(BARRA_C2_cell_nos_file, "data/processed_env_inputs/BARRA_C2_cell_nos.qs", format = "file"),
  tar_target(bathy_file, "data/processed_env_inputs/cell_bathy.parquet", format = "file"),
  tar_target(cell_bathy, arrow::read_parquet(bathy_file)),
  tar_target(
    cells_to_omit, 
    cell_bathy %>% 
      dplyr::filter(hz < 8) %>% 
      pull(cell_no)
  ),
  tar_target(
    BARRA_C2_cell_nos, 
    qs::qread(BARRA_C2_cell_nos_file) %>% 
      setdiff(cells_to_omit) #%>% head(1000)
  ),
  tar_target(
    BARRA_C2_cell_nos_chunked, 
    split(BARRA_C2_cell_nos, ceiling(seq_along(BARRA_C2_cell_nos)/chunk_size))
  ),
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    cell_inputs_files, list.files("data/processed_env_inputs", full.names = T) %>% 
      str_subset("cell_input_all"), 
    format = "file"
  ),
  tar_target(
    cell_inputs_chunked, 
    command = {
      df <- purrr::map_dfr(cell_inputs_files, function(file) {
        arrow::read_parquet(file) %>% 
          dplyr::filter(cell_no %in% BARRA_C2_cell_nos_chunked[[1]])
      }) %>% 
        mutate(state = as.factor(state))
      split(df, df$cell_no)
    },
    pattern = BARRA_C2_cell_nos_chunked,
    iteration = "list"
  ), 
  
  tar_target(state_inputs_file, "data/processed_env_inputs/state_input_timeseries.parquet", format = "file"),
  tar_target(
    state_inputs, 
    state_inputs_file %>% 
      purrr::map_dfr(arrow::read_parquet) %>% 
      dplyr::filter(yday <= 400 & state == states),
    pattern = states
  ),
  
  # Run model -----------------------------------------------------------------------------------------------------
  tar_target(init_biomass, 0.005 %>% units::set_units("g L-1") %>% 
               units::set_units("mg m-3") %>% units::drop_units()),
  tar_target(
    init_state,
    c(Nf = biomass_to_Nf(biomass = init_biomass, Q_rel = 0.5, spec_params = unlist(species_data), dry = T),
      Q_rel = 0.5),
    pattern = species_data,
    iteration = "list"
  ),
  
  # Total growth --------------------------------------------------------------------------------------------------
  tar_target(months_growth, 1:12),
  tar_target(
    start_dates_growth, 
    lubridate::ceiling_date(lubridate::make_date(2023, months_growth, 15), 'month') - lubridate::days(42+1), 
    pattern = months_growth
  ),
  tar_target(start_ydays_growth, lubridate::yday(start_dates_growth), pattern = start_dates_growth),
  tar_target(
    cell_growth,
    command = {
      purrr::map2_dfr(BARRA_C2_cell_nos_chunked[[1]], cell_inputs_chunked, function(cell, input) {
        mat <- grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    42,
          temperature =  input$T_input[start_ydays_growth:(start_ydays_growth + 42)],
          salinity =     input$S_input[start_ydays_growth:(start_ydays_growth + 42)],
          light =        input$I_input[start_ydays_growth:(start_ydays_growth + 42)],
          kW =           input$Kd_490[start_ydays_growth:(start_ydays_growth + 42)],
          velocity =     input$UV_input[start_ydays_growth:(start_ydays_growth + 42)],
          nitrate =      input$Ni_input[start_ydays_growth:(start_ydays_growth + 42)],
          ammonium =     input$Am_input[start_ydays_growth:(start_ydays_growth + 42)],
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(
            hc = 5,
            farmA = 50 * 50,
            turbulence = NA,
            d_top = 2.5,
            hz = cell_bathy$hz[cell_bathy$cell_no == cell]
          ), 
          spec_params =  species_data[[1]], #unlist(species_data),
          initials =     init_state
        )
        mat %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(cell_no = cell)
      })
    },
    pattern = cross(map(species_data, init_state), start_ydays_growth, map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked))
  ),
  tar_target(
    cell_growth_lims,
    command = {
      spec_params <- unlist(species_data)
      cell_growth %>% 
        mutate(
          lim = case_when(
            growth_rate == spec_params['mu'] ~ "No_limit",
            T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
            I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
            S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
            Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
            Nf_loss >= Ns_to_Nf ~ "Nf_loss",
            Ns_loss >= (up_Am+up_Ni) ~ "Ns_loss",
            T ~ NA
          ),
          lim = as.factor(lim),
          species = species_names,
          start = min(start_ydays_growth),
        )
    },
    pattern = map(cell_growth, cross(map(species_data, species_names), start_ydays_growth, BARRA_C2_cell_nos_chunked))
  ),
  tar_target(
    cell_growth_end,
    command = {
      spec_params <- unlist(species_data)
      all <- cell_growth %>% 
        mutate(t_start = start_ydays_growth)
      start <- all %>% 
        filter(t == t_start) %>% 
        mutate(TN_hm_start = (Ns + Nf) * hm * 5,
               TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw_start.mg = B_dw.mg) %>% 
        dplyr::select(cell_no, TN_hm_start, TN_start, B_dw_start.mg)
      end <- all %>% 
        filter(t == max(t)) %>% 
        rename(t_end = t) %>% 
        mutate(gr_mean = mean(growth_rate)/unname(spec_params['mu']),
               TN_hm = (Ns + Nf) * hm * 5,
               TN = (Ns + Nf) * 5) %>% 
        dplyr::select(t_start, t_end, cell_no, TN, TN_hm, gr_mean, B_dw.mg)
      merge(start, end, by = "cell_no") %>% 
        mutate(success = case_when(TN > TN_start ~ T, T ~ F),
               species = species_names)
      },
    pattern = map(cell_growth, cross(map(species_data, species_names), start_ydays_growth, BARRA_C2_cell_nos_chunked))
  ),
  
  # Supplemented growth -------------------------------------------------------------------------------------------
  tar_target(
    cell_growth_extra,
    command = {
      purrr::map2_dfr(BARRA_C2_cell_nos_chunked[[1]], cell_inputs_chunked, function(cell, input) {
        mat <- grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    42,
          temperature =  input$T_input[start_ydays_growth:(start_ydays_growth + 42)],
          salinity =     input$S_input[start_ydays_growth:(start_ydays_growth + 42)],
          light =        input$I_input[start_ydays_growth:(start_ydays_growth + 42)],
          kW =           input$Kd_490[start_ydays_growth:(start_ydays_growth + 42)],
          velocity =     input$UV_input[start_ydays_growth:(start_ydays_growth + 42)],
          nitrate =      input$Ni_input[start_ydays_growth:(start_ydays_growth + 42)] + Ni_add,
          ammonium =     input$Am_input[start_ydays_growth:(start_ydays_growth + 42)] + Am_add,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(
            hc = 5,
            farmA = 50 * 50,
            turbulence = NA,
            d_top = 2.5,
            hz = cell_bathy$hz[cell_bathy$cell_no == cell]
          ), 
          spec_params =  species_data[[1]], #unlist(species_data),
          initials =     init_state
        )
        mat %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(cell_no = cell)
      })
    },
    pattern = cross(map(species_data, init_state), start_ydays_growth, map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked))
  ),
  tar_target(
    cell_growth_extra_lims,
    command = {
      spec_params <- unlist(species_data)
      cell_growth_extra %>% 
        mutate(
          lim = case_when(
            growth_rate == spec_params['mu'] ~ "No_limit",
            T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
            I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
            S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
            Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
            Nf_loss >= Ns_to_Nf ~ "Nf_loss",
            Ns_loss >= (up_Am+up_Ni) ~ "Ns_loss",
            T ~ NA
          ),
          lim = as.factor(lim),
          species = species_names,
          start = min(start_ydays_growth),
        )
    },
    pattern = map(cell_growth_extra, cross(map(species_data, species_names), start_ydays_growth, BARRA_C2_cell_nos_chunked))
  ),
  tar_target(
    cell_growth_extra_end,
    command = {
      spec_params <- unlist(species_data)
      all <- cell_growth_extra %>% 
        mutate(t_start = start_ydays_growth)
      start <- all %>% 
        filter(t == t_start) %>% 
        mutate(TN_hm_start = (Ns + Nf) * hm * 5,
               TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw_start.mg = B_dw.mg) %>% 
        dplyr::select(cell_no, TN_hm_start, TN_start, B_dw_start.mg)
      end <- all %>% 
        filter(t == max(t)) %>% 
        rename(t_end = t) %>% 
        mutate(gr_mean = mean(growth_rate)/unname(spec_params['mu']),
               TN_hm = (Ns + Nf) * hm * 5,
               TN = (Ns + Nf) * 5) %>% 
        dplyr::select(t_start, t_end, cell_no, TN, TN_hm, gr_mean, B_dw.mg)
      merge(start, end, by = "cell_no") %>% 
        mutate(success = case_when(TN > TN_start ~ T, T ~ F),
               species = species_names)
    },
    pattern = map(cell_growth_extra, cross(map(species_data, species_names), start_ydays_growth, BARRA_C2_cell_nos_chunked))
  ),
  
  # Scenarios -----------------------------------------------------------------------------------------------------
  tar_target(
    name = scens,
    command = tibble::tribble(
      ~name,            ~num,   ~T_mod, ~bathy, ~kW_mod, ~S_mod, ~Ni_add,  ~Am_add,  ~UV_mod,
      "base",           1,      1,      50,     1,       1,      0,        0,        1,
      "fish_farm",      2,      1,      30,     1,       1,      Ni_add,   Am_add,   1,
      "deep_water",     3,      0.9,    100,    1.2,     1,      0,        0,        4,
      "shallow_bay",    4,      1.1,    15,     0.8,     1.1,    Ni_add,   0,        0.25,
      "estuary_mouth",  5,      1,      30,     0.8,     0.9,    Ni_add,   0,        1
    )
  ),
  tar_target(
    name = theo_names,
    command = unique(scens$name)
  ),
  tar_target(
    scen_growth,
    command = {
      ts <- scens[scens$name == theo_names, ]
      mat <- grow_macroalgae(
        start =        start_ydays_growth,
        grow_days =    42,
        temperature =  state_inputs$T_input_mean[start_ydays_growth:(start_ydays_growth + 42)] * ts$T_mod,
        salinity =     state_inputs$S_input_mean[start_ydays_growth:(start_ydays_growth + 42)] * ts$S_mod,
        light =        state_inputs$I_input_mean[start_ydays_growth:(start_ydays_growth + 42)],
        kW =           state_inputs$Kd_490_mean[start_ydays_growth:(start_ydays_growth + 42)] * ts$kW_mod,
        velocity =     state_inputs$UV_input_mean[start_ydays_growth:(start_ydays_growth + 42)] * ts$UV_mod,
        nitrate =      state_inputs$Ni_input_mean[start_ydays_growth:(start_ydays_growth + 42)] + ts$Ni_add,
        ammonium =     state_inputs$Am_input_mean[start_ydays_growth:(start_ydays_growth + 42)] + ts$Am_add,
        ni_uptake =    "linear",
        am_uptake =    "MM",
        site_params =  c(
          hc = 5,
          farmA = 50 * 50,
          turbulence = NA,
          d_top = 2.5,
          hz = ts$bathy
        ), 
        spec_params =  unlist(species_data),
        initials =     init_state
      )
      mat %>% 
        as.data.frame() %>% 
        tibble::remove_rownames() %>% 
        mutate(scen = theo_names,
               state = states,
               species = species_names)
    },
    pattern = cross(
      map(species_data, init_state, species_names),
      map(start_dates_growth, start_ydays_growth),
      theo_names,
      map(states, state_inputs)
    )
  ),
  tar_target(
    scen_growth_lims,
    command = {
      spec_params <- unlist(species_data)
      scen_growth %>% 
        mutate(
          lim = case_when(
            growth_rate == spec_params['mu'] ~ "No_limit",
            T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
            I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
            S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
            Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
            T ~ NA
          ),
          lim = as.factor(lim)
        )
    },
    pattern = map(scen_growth, cross(
      map(species_data, init_state, species_names),
      map(start_dates_growth, start_ydays_growth),
      theo_names,
      map(states, state_inputs)
    ))
  ),
  tar_target(
    scen_growth_end,
    command = {
      spec_params <- unlist(species_data)
      scen_growth %>% 
        as.data.frame() %>% 
        tibble::remove_rownames() %>%
        mutate(start = min(t)) %>% 
        filter(t == max(t)) %>% 
        mutate(gr_mean = mean(growth_rate)/unname(spec_params['mu']),
               TN_hm = (Ns + Nf) * hm * 5,
               TN = (Ns + Nf) * 5) %>% 
        dplyr::select(start, state, scen, t, TN, TN_hm, gr_mean, B_dw.mg, species)
    },
    pattern = scen_growth
  )
)

# nolint end

