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
  controller = crew_controller_local(workers = 10, seconds_idle = 120),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R"))
)

units::remove_unit("molN")
units::remove_unit("gN")
units::install_unit("gN")
units::install_unit(symbol = "molN", def = "14.007 gN")
Ni_add <- 3.5 %>% units::set_units("umolN L-1") %>% units::set_units("mgN m-3") %>% units::drop_units()
Am_add <- 6.5 %>% units::set_units("umolN L-1") %>% units::set_units("mgN m-3") %>% units::drop_units()

chunk_size <- 175

a_armata_gametophyte <- targets::tar_read(a_armata_gametophyte, store = "targets_outputs/_species")
a_taxiformis_gametophyte <- targets::tar_read(a_taxiformis_gametophyte, store = "targets_outputs/_species")

init_biomass <- 0.005 %>% units::set_units("g L-1") %>% units::set_units("mg m-3") %>% units::drop_units()
init_Q_rel <- 0.25
init_state_armata <- c(Q_rel = init_Q_rel)
init_state_armata <- c(init_state_armata, 
                       macrogrow::biomass_to_Nf(biomass = init_biomass, Q_rel = init_state_armata["Q_rel"], 
                                                spec_params = a_armata_gametophyte, dry = T))
init_state_taxiformis <- c(Q_rel = init_Q_rel)
init_state_taxiformis <- c(init_state_taxiformis, 
                           macrogrow::biomass_to_Nf(biomass = init_biomass, Q_rel = init_state_taxiformis["Q_rel"], 
                                                    spec_params = a_taxiformis_gametophyte, dry = T))


list(
  # Species data --------------------------------------------------------------------------------------------------
  tar_target(species_names, c("armata", "taxiformis")),

  tar_target(BARRA_C2_cell_coords_file, "data/processed_env_inputs/BARRA_C2_cell_coords.parquet", format = "file"),
  tar_target(BARRA_C2_cell_coords, arrow::read_parquet(BARRA_C2_cell_coords_file)),
  tar_target(states, unique(BARRA_C2_cell_coords$state)),
  
  # Cell data -----------------------------------------------------------------------------------------------------
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
      setdiff(cells_to_omit)
  ),
  tar_target(
    BARRA_C2_cell_nos_chunked, 
    split(BARRA_C2_cell_nos, ceiling(seq_along(BARRA_C2_cell_nos)/chunk_size))
  ),
  
  # Cell inputs ---------------------------------------------------------------------------------------------------
  tar_target(
    cell_inputs_files, 
    file.path("data", "processed_env_inputs", "cell_input_all") %>% 
      list.files(full.names = T) %>% 
      str_subset("cell_input_all"), 
    format = "file"
  ),
  tar_target(
    cell_inputs_chunked, 
    command = {
      df <- purrr::map_dfr(cell_inputs_files, function(file) {
        arrow::read_parquet(file) %>% 
          dplyr::filter(cell_no %in% BARRA_C2_cell_nos_chunked[[1]])
      })
      split(df, df$cell_no)
    },
    pattern = BARRA_C2_cell_nos_chunked
  ), 
  
  tar_target(state_inputs_file, "data/processed_env_inputs/state_input_timeseries.parquet", format = "file"),
  tar_target(
    state_inputs, 
    state_inputs_file %>% 
      purrr::map_dfr(arrow::read_parquet) %>% 
      dplyr::filter(state == states),
    pattern = states
  ),
  
  # Total growth --------------------------------------------------------------------------------------------------
  tar_target(months_growth, 1:12),
  tar_target(start_ydays_growth, lubridate::make_date(2023, months_growth, 1) %>% lubridate::yday()),
  tar_target(end_ydays_growth, c(start_ydays_growth[2:12], 366)),
  tar_target(
    cell_growth_armata,
    command = {
      cell_nos <- BARRA_C2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})

      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    end_ydays_growth-start_ydays_growth,
          temperature =  input$T_input,
          salinity =     input$S_input,
          light =        input$I_input,
          kW =           input$Kd_490,
          velocity =     input$UV_input,
          nitrate =      input$Ni_input,
          ammonium =     input$Am_input,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = unique(input$hz)), 
          spec_params =  a_armata_gametophyte,
          initials =     init_state_armata
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            cell_no = cell,
            start = start_ydays_growth,
            lim = case_when(
              T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
              I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
              S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
              Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
              Nf_loss >= Ns_to_Nf ~ "Nf_loss",
              Ns_loss >= (up_Am + up_Ni) ~ "Ns_loss",
              T ~ NA
            ),
            lim = as.factor(lim),
            species = as.factor("armata")
          )
      }) %>% dplyr::bind_rows()
    },
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked))
  ),
  
  tar_target(
    cell_growth_taxiformis,
    command = {
      cell_nos <- BARRA_C2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})
      
      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    end_ydays_growth-start_ydays_growth,
          temperature =  input$T_input,
          salinity =     input$S_input,
          light =        input$I_input,
          kW =           input$Kd_490,
          velocity =     input$UV_input,
          nitrate =      input$Ni_input,
          ammonium =     input$Am_input,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = unique(input$hz)), 
          spec_params =  a_taxiformis_gametophyte,
          initials =     init_state_taxiformis
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            cell_no = cell,
            start = start_ydays_growth,
            lim = case_when(
              T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
              I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
              S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
              Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
              Nf_loss >= Ns_to_Nf ~ "Nf_loss",
              Ns_loss >= (up_Am + up_Ni) ~ "Ns_loss",
              T ~ NA
            ),
            lim = as.factor(lim),
            species = as.factor("taxiformis")
          )
      }) %>% dplyr::bind_rows()
    },
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked))
  ),
  
  tar_target(
    cell_growth_end,
    command = {
      start <- cell_growth_armata %>% 
        filter(t == start) %>% 
        mutate(hm_start = hm, TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw.mg_start = B_dw.mg) %>% 
        dplyr::select(start, cell_no, hm_start, TN_start, B_dw.mg_start, species)
      end <- cell_growth_armata %>% 
        filter(t == max(t)) %>% 
        rename(end = t, B_dw.mg_end = B_dw.mg) %>% 
        mutate(gr_mean = mean(growth_rate), hm_end = hm, TN_end = (Ns + Nf) * 5) %>% 
        dplyr::select(start, end, cell_no, hm_end, TN_end, gr_mean, B_dw.mg_end)
      arma <- merge(start, end, by = c("cell_no", "start")) %>% 
        mutate(success = case_when(TN_end > TN_start ~ T, T ~ F))
      
      start <- cell_growth_taxiformis %>% 
        filter(t == start) %>% 
        mutate(hm_start = hm, TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw.mg_start = B_dw.mg) %>% 
        dplyr::select(start, cell_no, hm_start, TN_start, B_dw.mg_start, species)
      end <- cell_growth_taxiformis %>% 
        filter(t == max(t)) %>% 
        rename(end = t, B_dw.mg_end = B_dw.mg) %>% 
        mutate(gr_mean = mean(growth_rate), hm_end = hm, TN_end = (Ns + Nf) * 5) %>% 
        dplyr::select(start, end, cell_no, hm_end, TN_end, gr_mean, B_dw.mg_end)
      taxi <- merge(start, end, by = c("cell_no", "start")) %>% 
        mutate(success = case_when(TN_end > TN_start ~ T, T ~ F))
      
      rbind(arma, taxi)
      },
    pattern = map(cell_growth_armata, cell_growth_taxiformis)
  ),
  
  # Supplemented growth -------------------------------------------------------------------------------------------
  tar_target(
    cell_growth_supp_armata,
    command = {
      cell_nos <- BARRA_C2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})
      
      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    end_ydays_growth-start_ydays_growth,
          temperature =  input$T_input,
          salinity =     input$S_input,
          light =        input$I_input,
          kW =           input$Kd_490,
          velocity =     input$UV_input,
          nitrate =      input$Ni_input + Ni_add,
          ammonium =     input$Am_input + Am_add,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = unique(input$hz)), 
          spec_params =  a_armata_gametophyte,
          initials =     init_state_armata
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            cell_no = cell,
            start = start_ydays_growth,
            lim = case_when(
              T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
              I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
              S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
              Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
              Nf_loss >= Ns_to_Nf ~ "Nf_loss",
              Ns_loss >= (up_Am + up_Ni) ~ "Ns_loss",
              T ~ NA
            ),
            lim = as.factor(lim),
            species = as.factor("armata")
          )
      }) %>% dplyr::bind_rows()
    },
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked))
  ),
  tar_target(
    cell_growth_supp_taxiformis,
    command = {
      cell_nos <- BARRA_C2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})
      
      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    end_ydays_growth-start_ydays_growth,
          temperature =  input$T_input,
          salinity =     input$S_input,
          light =        input$I_input,
          kW =           input$Kd_490,
          velocity =     input$UV_input,
          nitrate =      input$Ni_input + Ni_add,
          ammonium =     input$Am_input + Am_add,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = unique(input$hz)), 
          spec_params =  a_taxiformis_gametophyte,
          initials =     init_state_taxiformis
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            cell_no = cell,
            start = start_ydays_growth,
            lim = case_when(
              T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
              I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
              S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
              Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
              Nf_loss >= Ns_to_Nf ~ "Nf_loss",
              Ns_loss >= (up_Am + up_Ni) ~ "Ns_loss",
              T ~ NA
            ),
            lim = as.factor(lim),
            species = as.factor("taxiformis")
          )
      }) %>% dplyr::bind_rows()
    },
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_C2_cell_nos_chunked, cell_inputs_chunked))
  ),
  
  tar_target(
    cell_growth_supp_end,
    command = {
      start <- cell_growth_supp_armata %>% 
        filter(t == start) %>% 
        mutate(hm_start = hm, TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw.mg_start = B_dw.mg) %>% 
        dplyr::select(start, cell_no, hm_start, TN_start, B_dw.mg_start, species)
      end <- cell_growth_supp_armata %>% 
        filter(t == max(t)) %>% 
        rename(end = t, B_dw.mg_end = B_dw.mg) %>% 
        mutate(gr_mean = mean(growth_rate), hm_end = hm, TN_end = (Ns + Nf) * 5) %>% 
        dplyr::select(start, end, cell_no, hm_end, TN_end, gr_mean, B_dw.mg_end)
      arma <- merge(start, end, by = c("cell_no", "start")) %>% 
        mutate(success = case_when(TN_end > TN_start ~ T, T ~ F))
      
      start <- cell_growth_supp_taxiformis %>% 
        filter(t == start) %>% 
        mutate(hm_start = hm, TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw.mg_start = B_dw.mg) %>% 
        dplyr::select(start, cell_no, hm_start, TN_start, B_dw.mg_start, species)
      end <- cell_growth_supp_taxiformis %>% 
        filter(t == max(t)) %>% 
        rename(end = t, B_dw.mg_end = B_dw.mg) %>% 
        mutate(gr_mean = mean(growth_rate), hm_end = hm, TN_end = (Ns + Nf) * 5) %>% 
        dplyr::select(start, end, cell_no, hm_end, TN_end, gr_mean, B_dw.mg_end)
      taxi <- merge(start, end, by = c("cell_no", "start")) %>% 
        mutate(success = case_when(TN_end > TN_start ~ T, T ~ F))
      
      rbind(arma, taxi)
    },
    pattern = map(cell_growth_supp_armata, cell_growth_supp_taxiformis)
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
    scen_growth_armata,
    command = {
      ts <- scens[scens$name == theo_names, ]
      st_ins <- split(state_inputs, state_inputs$state)
      
      purrr::map2(as.character(states), st_ins, function(st, ins) {
        grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    end_ydays_growth-start_ydays_growth,
          temperature =  ins$T_input_mean[start_ydays_growth:end_ydays_growth] * ts$T_mod,
          salinity =     ins$S_input_mean[start_ydays_growth:end_ydays_growth] * ts$S_mod,
          light =        ins$I_input_mean[start_ydays_growth:end_ydays_growth],
          kW =           ins$Kd_490_mean[start_ydays_growth:end_ydays_growth] * ts$kW_mod,
          velocity =     ins$UV_input_mean[start_ydays_growth:end_ydays_growth] * ts$UV_mod,
          nitrate =      ins$Ni_input_mean[start_ydays_growth:end_ydays_growth] + ts$Ni_add,
          ammonium =     ins$Am_input_mean[start_ydays_growth:end_ydays_growth] + ts$Am_add,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = ts$bathy), 
          spec_params =  a_armata_gametophyte,
          initials =     init_state_armata
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            scen = as.factor(theo_names),
            state = as.factor(st),
            start = start_ydays_growth,
            lim = case_when(
              T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
              I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
              S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
              Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
              Nf_loss >= Ns_to_Nf ~ "Nf_loss",
              Ns_loss >= (up_Am + up_Ni) ~ "Ns_loss",
              T ~ NA
            ),
            lim = as.factor(lim),
            species = as.factor("armata")
          )
      }) %>% dplyr::bind_rows()
    },
    pattern = cross(map(start_ydays_growth, end_ydays_growth), theo_names)
  ),
  tar_target(
    scen_growth_taxiformis,
    command = {
      ts <- scens[scens$name == theo_names, ]
      st_ins <- split(state_inputs, state_inputs$state)
      
      purrr::map2(as.character(states), st_ins, function(st, ins) {
        grow_macroalgae(
          start =        start_ydays_growth,
          grow_days =    end_ydays_growth-start_ydays_growth,
          temperature =  ins$T_input_mean[start_ydays_growth:end_ydays_growth] * ts$T_mod,
          salinity =     ins$S_input_mean[start_ydays_growth:end_ydays_growth] * ts$S_mod,
          light =        ins$I_input_mean[start_ydays_growth:end_ydays_growth],
          kW =           ins$Kd_490_mean[start_ydays_growth:end_ydays_growth] * ts$kW_mod,
          velocity =     ins$UV_input_mean[start_ydays_growth:end_ydays_growth] * ts$UV_mod,
          nitrate =      ins$Ni_input_mean[start_ydays_growth:end_ydays_growth] + ts$Ni_add,
          ammonium =     ins$Am_input_mean[start_ydays_growth:end_ydays_growth] + ts$Am_add,
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = ts$bathy), 
          spec_params =  a_taxiformis_gametophyte,
          initials =     init_state_taxiformis
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            scen = as.factor(theo_names),
            state = as.factor(st),
            start = start_ydays_growth,
            lim = case_when(
              T_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "T_lim",
              I_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "I_lim",
              S_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "S_lim",
              Q_lim == pmin(T_lim, I_lim, S_lim, Q_lim) ~ "Q_lim",
              Nf_loss >= Ns_to_Nf ~ "Nf_loss",
              Ns_loss >= (up_Am + up_Ni) ~ "Ns_loss",
              T ~ NA
            ),
            lim = as.factor(lim),
            species = as.factor("taxiformis")
          )
      }) %>% dplyr::bind_rows()
    },
    pattern = cross(map(start_ydays_growth, end_ydays_growth), theo_names)
  ),
  tar_target(
    scen_growth_end,
    command = {
      start <- scen_growth_armata %>% 
        filter(t == start) %>% 
        mutate(hm_start = hm, TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw.mg_start = B_dw.mg) %>% 
        dplyr::select(start, state, scen, hm_start, TN_start, B_dw.mg_start, species)
      end <- scen_growth_armata %>% 
        filter(t == max(t)) %>% 
        rename(end = t, B_dw.mg_end = B_dw.mg) %>% 
        mutate(gr_mean = mean(growth_rate), hm_end = hm, TN_end = (Ns + Nf) * 5) %>% 
        dplyr::select(start, end, state, scen, hm_end, TN_end, gr_mean, B_dw.mg_end)
      arma <- merge(start, end, by = c("state", "scen", "start")) %>% 
        mutate(success = case_when(TN_end > TN_start ~ T, T ~ F))
      
      start <- scen_growth_taxiformis %>% 
        filter(t == start) %>% 
        mutate(hm_start = hm, TN_start = (Ns + Nf) * 5) %>% 
        rename(B_dw.mg_start = B_dw.mg) %>% 
        dplyr::select(start, state, scen, hm_start, TN_start, B_dw.mg_start, species)
      end <- scen_growth_taxiformis %>% 
        filter(t == max(t)) %>% 
        rename(end = t, B_dw.mg_end = B_dw.mg) %>% 
        mutate(gr_mean = mean(growth_rate), hm_end = hm, TN_end = (Ns + Nf) * 5) %>% 
        dplyr::select(start, end, state, scen, hm_end, TN_end, gr_mean, B_dw.mg_end)
      taxi <- merge(start, end, by = c("state", "scen", "start")) %>% 
        mutate(success = case_when(TN_end > TN_start ~ T, T ~ F))
      
      rbind(arma, taxi)
    },
    pattern = map(scen_growth_armata, scen_growth_taxiformis)
  )
)

# nolint end

