library(targets)
library(tarchetypes)

# Set up ----------------------------------------------------------------------------------------------------------
tar_option_set(
  packages = c("tidyverse", "ggplot2", "arrow", "lubridate", "macrogrow", "units", "macrogrow", "qs", "qs2"), 
  format = "qs",
  controller = crew::crew_controller_local(workers = 2, seconds_idle = 60),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_targets_functions.R")) 
)

list(
  ## Asparagopsis armata ------------------------------------------------------------------------------------------
  tar_target(
    name = a_armata_gametophyte_file,
    command = file.path("data", "asparagopsis_armata_gametophyte.csv"),
    format = "file"
  ),
  tar_target(
    name = a_armata_gametophyte,
    command = {
      vec <- read.csv(a_armata_gametophyte_file)$value
      names(vec) <- read.csv(a_armata_gametophyte_file)$parameter
      vec <- vec[!names(vec) %in% c("M_am", "C_am", "V_ni", "K_ni", "V_ot", "K_ot", "M_ot", "C_ot")]
    }
  ),
  
  ## Asparagopsis taxiformis --------------------------------------------------------------------------------------
  tar_target(
    name = a_taxiformis_gametophyte_file,
    command = file.path("data", "asparagopsis_taxiformis_gametophyte.csv"),
    format = "file"
  ),
  tar_target(
    name = a_taxiformis_gametophyte,
    command = {
      vec <- read.csv(a_taxiformis_gametophyte_file)$value
      names(vec) <- read.csv(a_taxiformis_gametophyte_file)$parameter
      vec <- vec[!names(vec) %in% c("M_am", "C_am", "V_ni", "K_ni", "V_ot", "K_ot", "M_ot", "C_ot")]
    }
  ),
  
  # Sensitivities -------------------------------------------------------------------------------------------------
  ## Standard values ----------------------------------------------------------------------------------------------
  tar_target(
    sens_conditions,
    data.frame(t_span = rep(1:732, 2)) %>% 
      mutate(
        Ni_input = 26.60 + 25.90 * sin((t_span * pi - 310) / 182.5),
        Am_input = 3.15 + 3 * sin((t_span * pi - 325) / 182.5),
        I_input  = 325 + 250 * sin((t_span * pi + 300) / 182.5),
        S_input  = 35.1 + 0.125 * sin((t_span * pi + 450) / 182.5),
        T_input  = c(26.75 + 3.5 * sin((1:732 * pi + 220) / 182.5), 
                     16.75 + 4.75 * sin((1:732 * pi + 120) / 182.5)),
        T_level  = rep(c(1,2), each = 732),
        U_input  = 0.3 + 0.25 * sin((t_span * pi + 180) / 182.5)
      )
  ),
  tar_target(T_levels, c(1,2)),
  
  ## Culture conditions -------------------------------------------------------------------------------------------
  tar_target(name = starts, command = seq(5, 365, 10)),
  tar_target(name = culture_depths, command = c(0.5, 2.5, 5)),
  tar_target(name = factors, command = c(0.95, 1, 1.05)),
  tar_target(
    name = site_params_sens,
    command = c(
      d_top = culture_depths,
      hc = 5,
      farmA = 50 * 50,
      hz = 50,
      kW = Kd_from_Secchi(15),
      turbulence = NA
    ),
    pattern = culture_depths,
    iteration = "list"
  ),
  
  tar_target(name = param_names, command = names(a_armata_gametophyte)),

  # Sensitivity armata --------------------------------------------------------------------------------------------
  tar_target(
    name = init_state,
    command = c(
      Nf = biomass_to_Nf(
        biomass = 0.0005 %>% units::set_units("g L-1") %>% units::set_units("mg m-3") %>% units::drop_units(),
        Q_rel = 0.5,
        spec_params = unlist(a_armata_gametophyte),
        dry = T
      ),
      Q_rel = 0.5
    )
  ),
  tar_target(
    name = spec_params_arma,
    command = adj_params(
      params = unlist(a_armata_gametophyte),
      focus_param = param_names,
      factor = factors
    ),
    pattern = cross(param_names, factors),
    iteration = "list"
  ),
  tar_target(
    name = sensitivity_run_arma,
    command = grow_macroalgae(
      start = as.Date(parse_date_time(x = paste0("2022-", starts), orders = "yj")), 
      grow_days = 42,
      temperature = sens_conditions$T_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      salinity = sens_conditions$S_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      light = sens_conditions$I_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      velocity = sens_conditions$U_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)],
      nitrate = sens_conditions$Ni_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      ammonium = sens_conditions$Am_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      ni_uptake = "linear", 
      am_uptake = "MM",
      site_params = unlist(site_params_sens),
      spec_params = unlist(spec_params_arma),
      initials = unlist(init_state)
    ),
    pattern = cross(starts, map(spec_params_arma, cross(param_names, factors)), map(culture_depths, site_params_sens), T_levels)
  ),
  
  tar_target(
    name = total_N_end_arma,
    command = cbind(sensitivity_run_arma[,'date'], sensitivity_run_arma[,'Nf'] + sensitivity_run_arma[,'Ns']) %>% 
      as.data.frame() %>% 
      remove_rownames() %>% 
      rename(date = V1, TN = V2) %>% 
      filter(date == max(date)) %>% 
      mutate(
        date = as.Date(date),
        species = "A. armata",
        cult_dep = culture_depths,
        param = param_names,
        factor = factors,
        T_level = T_levels
      ),
    pattern = map(sensitivity_run_arma, cross(starts, map(spec_params_arma, cross(param_names, factors)), map(culture_depths, site_params_sens), T_levels))
  ),

  # Sensitivity taxiformis -----------------------------------------------------------------------------------------
  tar_target(
    name = spec_params_taxi,
    command = adj_params(
      params = unlist(a_taxiformis_gametophyte),
      focus_param = param_names,
      factor = factors
    ),
    pattern = cross(param_names, factors),
    iteration = "list"
  ),
  tar_target(
    name = sensitivity_run_taxi,
    command = grow_macroalgae(
      start = as.Date(parse_date_time(x = paste0("2022-", starts), orders = "yj")), 
      grow_days = 42,
      temperature = sens_conditions$T_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      salinity = sens_conditions$S_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      light = sens_conditions$I_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      velocity = sens_conditions$U_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)],
      nitrate = sens_conditions$Ni_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      ammonium = sens_conditions$Am_input[sens_conditions$T_level == T_levels & sens_conditions$t_span %in% starts:(starts + 42)], 
      ni_uptake = "linear", 
      am_uptake = "MM",
      site_params = unlist(site_params_sens),
      spec_params = unlist(spec_params_taxi),
      initials = unlist(init_state)
    ),
    pattern = cross(starts, map(spec_params_taxi, cross(param_names, factors)), map(culture_depths, site_params_sens), T_levels)
  ),

  tar_target(
    name = total_N_end_taxi,
    command = cbind(sensitivity_run_taxi[,'date'], sensitivity_run_taxi[,'Nf'] + sensitivity_run_taxi[,'Ns']) %>% 
      as.data.frame() %>% 
      remove_rownames() %>% 
      rename(date = V1, TN = V2) %>% 
      filter(date == max(date)) %>% 
      mutate(
        date = as.Date(date),
        species = "A. taxiformis",
        cult_dep = culture_depths,
        param = param_names,
        factor = factors,
        T_level = T_levels
      ),
    pattern = map(sensitivity_run_taxi, cross(starts, map(spec_params_arma, cross(param_names, factors)), map(culture_depths, site_params_sens), T_levels))
  )
)

