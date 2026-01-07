# nolint start

# Setup -----------------------------------------------------------------------------------------------------------
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
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("dplyr", "macrogrow", "geosphere", "stringr"),
  format = "qs", 
  controller = crew_controller_local(workers = 12, seconds_idle = 120),
  workspace_on_error = T
)

tar_source(
  files = c(
    file.path("R_scripts", "00_functions.R")
  ))

# Since we are always working in N, mol and molN are interchangeable
units::remove_unit("molN")
units::remove_unit("gN")
units::install_unit("gN")
units::install_unit(symbol = "molN", def = "14.007 gN")

# Set standard supplements to nitrate and ammonium timeseries
Ni_add <- 3.5 %>% 
  units::set_units("umolN L-1") %>% 
  units::set_units("mgN m-3") %>% 
  units::drop_units()
Am_add <- 6.5 %>% 
  units::set_units("umolN L-1") %>% 
  units::set_units("mgN m-3") %>% 
  units::drop_units()

# Growth units (cell * time) are chunked for more efficient running
chunk_size <- 100

# Load previously-processed species data
armata <- c(
  V_am = 4.991988e+00,
  K_am = 1.402071e+02,
  M_ni = 0.000000e+00,
  C_ni = 1.138183e-01,
  Q_min = 1.943162e+01,
  Q_max = 4.303730e+01,
  K_c = 1.200000e+01,
  mu = 3.200000e-01,
  D_m = 1.148264e-02,
  D_ve = 1.562892e-02,
  a_cs = 1.530430e-04,
  I_o = 1.500000e+02,
  T_opt = 1.909091e+01,
  T_min = 8.895487e+00,
  T_max = 2.242924e+01,
  S_opt = 3.525000e+01,
  S_min = 2.115000e+01,
  S_max = 4.165000e+01,
  h_a = 3.850000e+03,
  h_b = 1.750000e+00,
  h_c = 1.000000e-03,
  h_max = 2.600000e-01,
  DWWW = 7.680000e+00
)
taxiformis <- armata
taxiformis["T_opt"] <- 24.6153846153846
taxiformis["T_min"] <- 16.0326603325416
taxiformis["T_max"] <- 29.2262470308789

biom <- macrogrow::biomass_to_Nf(
  biomass = 0.0555 %>% units::set_units("g L-1") %>% units::set_units("mg m-3") %>% units::drop_units(),
  Q_rel = 0.75,
  spec_params = armata,
  dry = T
)
init_state_armata <- init_state_taxiformis <- c(biom, Q_rel = 0.75)

list(
  # Species data --------------------------------------------------------------------------------------------------
  tar_target(species_names, c("armata", "taxiformis")),

  # Cell data -----------------------------------------------------------------------------------------------------
  tar_target(BARRA_R2_cell_coords_file, "data/processed_env_inputs/BARRA_R2_cell_coords.parquet", format = "file"),
  tar_target(BARRA_R2_cell_coords, arrow::read_parquet(BARRA_R2_cell_coords_file)),
  tar_target(states, unique(BARRA_R2_cell_coords$state)),
  tar_target(BARRA_R2_cell_nos_file, "data/processed_env_inputs/BARRA_R2_cell_nos.qs", format = "file"),
  tar_target(bathy_file, "data/processed_env_inputs/cell_bathy.parquet", format = "file"),
  tar_target(cell_bathy, arrow::read_parquet(bathy_file)),
  tar_target(
    cells_to_omit, 
    cell_bathy %>% 
      dplyr::filter(hz < 8) %>% 
      pull(cell_no)
  ),
  tar_target(
    BARRA_R2_cell_nos, 
    qs::qread(BARRA_R2_cell_nos_file) %>% 
      setdiff(cells_to_omit) #%>% sample(size = 1195)
  ),
  tar_target(
    BARRA_R2_cell_nos_chunked, 
    split(BARRA_R2_cell_nos, ceiling(seq_along(BARRA_R2_cell_nos)/chunk_size))
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
          dplyr::filter(cell_no %in% BARRA_R2_cell_nos_chunked[[1]])
      })
      split(df, df$cell_no)
    },
    pattern = BARRA_R2_cell_nos_chunked
  ), 
  
  tar_target(state_inputs_file, "data/processed_env_inputs/state_input_timeseries.parquet", format = "file"),
  tar_target(
    state_inputs, 
    state_inputs_file %>% 
      purrr::map_dfr(arrow::read_parquet) %>% 
      dplyr::filter(state == states),
    pattern = states
  ),
  
  # Test growth ---------------------------------------------------------------------------------------------------
  tar_target(
    test_growth_armata,
    command = {
      state_ins <- split(state_inputs, state_inputs$state)
      purrr::map(state_ins, function(input) {
        grow_macroalgae(
          t =            start_ydays_growth:end_ydays_growth,
          temperature =  input$T_input_mean[start_ydays_growth:end_ydays_growth],
          salinity =     input$S_input_mean[start_ydays_growth:end_ydays_growth],
          light =        input$I_input_mean[start_ydays_growth:end_ydays_growth],
          kW =           input$Kd_490_mean[start_ydays_growth:end_ydays_growth],
          velocity =     input$UV_input_mean[start_ydays_growth:end_ydays_growth],
          nitrate =      input$Ni_input_mean[start_ydays_growth:end_ydays_growth],
          ammonium =     input$Am_input_mean[start_ydays_growth:end_ydays_growth],
          ni_uptake =    "linear",
          am_uptake =    "MM",
          site_params =  c(hc = 5, farmA = 50 * 50, turbulence = NA, d_top = 2.5, hz = 30), 
          spec_params =  armata,
          initials =     init_state_armata
        ) %>% 
          as.data.frame() %>% 
          tibble::remove_rownames() %>% 
          mutate(
            state = unique(input$state),
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
    pattern = map(start_ydays_growth, end_ydays_growth)
  ),
  
  # Total growth --------------------------------------------------------------------------------------------------
  tar_target(months_growth, 1:12),
  tar_target(start_ydays_growth, lubridate::make_date(2023, months_growth, 1) %>% lubridate::yday()),
  tar_target(end_ydays_growth, c(start_ydays_growth[2:12], 366)),
  tar_target(
    cell_growth_armata,
    command = {
      cell_nos <- BARRA_R2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})

      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          t =            start_ydays_growth:end_ydays_growth,
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
          spec_params = armata,
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
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_R2_cell_nos_chunked, cell_inputs_chunked))
  ),
  
  tar_target(
    cell_growth_taxiformis,
    command = {
      cell_nos <- BARRA_R2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})
      
      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          t =            start_ydays_growth:end_ydays_growth,
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
          spec_params =  taxiformis,
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
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_R2_cell_nos_chunked, cell_inputs_chunked))
  ),
  
  tar_target(
    cell_growth_end,
    command = {
      st <- cell_growth_armata %>% 
        group_by(species, cell_no, start) %>% 
        slice_min(t) %>% 
        ungroup() %>% 
        mutate(TN_start = (Ns + Nf) * hm) %>% 
        dplyr::select(start, cell_no, TN_start, species)
      en <- cell_growth_armata %>% 
        group_by(species, cell_no, start) %>% 
        slice_max(t) %>% 
        ungroup() %>% 
        rename(end = t) %>% 
        mutate(TN_end = (Ns + Nf) * hm) %>% 
        dplyr::select(species, start, cell_no, TN_end)
      arma <- merge(st, en, by = c("species", "cell_no", "start")) %>% 
        mutate(TN_rem = TN_end - TN_start,
               success = case_when(TN_rem > 0 ~ T, T ~ F),
               start = as.integer(start))
      
      st <- cell_growth_taxiformis %>% 
        group_by(species, cell_no, start) %>% 
        slice_min(t) %>% 
        ungroup() %>% 
        mutate(TN_start = (Ns + Nf) * hm) %>% 
        dplyr::select(start, cell_no, TN_start, species)
      en <- cell_growth_taxiformis %>% 
        group_by(species, cell_no, start) %>% 
        slice_max(t) %>% 
        ungroup() %>% 
        rename(end = t) %>% 
        mutate(TN_end = (Ns + Nf) * hm) %>% 
        dplyr::select(species, start, cell_no, TN_end)
      taxi <- merge(st, en, by = c("species", "cell_no", "start")) %>% 
        mutate(TN_rem = TN_end - TN_start,
               success = case_when(TN_rem > 0 ~ T, T ~ F),
               start = as.integer(start))
      
      rbind(arma, taxi)
      },
    pattern = map(cell_growth_armata, cell_growth_taxiformis)
  ),

  # Supplemented growth -------------------------------------------------------------------------------------------
  tar_target(
    cell_growth_supp_armata,
    command = {
      cell_nos <- BARRA_R2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})
      
      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          t =            start_ydays_growth:end_ydays_growth,
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
          spec_params =  armata,
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
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_R2_cell_nos_chunked, cell_inputs_chunked))
  ),
  tar_target(
    cell_growth_supp_taxiformis,
    command = {
      cell_nos <- BARRA_R2_cell_nos_chunked[[1]]
      cell_ins <- lapply(cell_inputs_chunked, function(df) {df[df$yday %in% start_ydays_growth:end_ydays_growth, ]})
      
      purrr::map2(cell_nos, cell_ins, function(cell, input) {
        grow_macroalgae(
          t =            start_ydays_growth:end_ydays_growth,
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
          spec_params =  taxiformis,
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
    pattern = cross(map(start_ydays_growth, end_ydays_growth), map(BARRA_R2_cell_nos_chunked, cell_inputs_chunked))
  ),
  
  tar_target(
    cell_growth_supp_end,
    command = {
      st <- cell_growth_supp_armata %>% 
        group_by(species, cell_no, start) %>% 
        slice_min(t) %>% 
        ungroup() %>% 
        mutate(TN_start = (Ns + Nf) * hm) %>% 
        dplyr::select(start, cell_no, TN_start, species)
      en <- cell_growth_supp_armata %>% 
        group_by(species, cell_no, start) %>% 
        slice_max(t) %>% 
        ungroup() %>% 
        rename(end = t) %>% 
        mutate(TN_end = (Ns + Nf) * hm) %>% 
        dplyr::select(species, start, cell_no, TN_end)
      arma <- merge(st, en, by = c("species", "cell_no", "start")) %>% 
        mutate(TN_rem = TN_end - TN_start,
               success = case_when(TN_rem > 0 ~ T, T ~ F),
               start = as.integer(start))
      
      st <- cell_growth_supp_taxiformis %>% 
        group_by(species, cell_no, start) %>% 
        slice_min(t) %>% 
        ungroup() %>% 
        mutate(TN_start = (Ns + Nf) * hm) %>% 
        dplyr::select(start, cell_no, TN_start, species)
      en <- cell_growth_supp_taxiformis %>% 
        group_by(species, cell_no, start) %>% 
        slice_max(t) %>% 
        ungroup() %>% 
        rename(end = t) %>% 
        mutate(TN_end = (Ns + Nf) * hm) %>% 
        dplyr::select(species, start, cell_no, TN_end)
      taxi <- merge(st, en, by = c("species", "cell_no", "start")) %>% 
        mutate(TN_rem = TN_end - TN_start,
               success = case_when(TN_rem > 0 ~ T, T ~ F),
               start = as.integer(start))
      
      rbind(arma, taxi)
    },
    pattern = map(cell_growth_supp_armata, cell_growth_supp_taxiformis)
  )
)

# nolint end

