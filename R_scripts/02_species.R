library(targets)
library(tarchetypes)

# Set up ----------------------------------------------------------------------------------------------------------
tar_option_set(
  packages = c("tidyverse", "ggplot2", "arrow", "lubridate", "macrogrow", "units", "macrogrow", "qs", "qs2"), 
  format = "qs",
  controller = crew::crew_controller_local(workers = 6, seconds_idle = 60),
  workspace_on_error = T
)

tar_source(
  files = c(file.path("R_scripts", "00_functions.R")) 
)

list(
  ## Asparagopsis armata ------------------------------------------------------------------------------------------
  tar_target(a_armata_file, file.path("data", "asparagopsis_armata.qs"), format = "file"),
  tar_target(
    name = a_armata,
    command = {
      vec <- qs::qread(a_armata_file)
      vec <- vec[!is.na(vec)]
      vec <- vec[!names(vec) %in% c("M_am", "C_am", "V_ot", "K_ot", "M_ot", "C_ot")]
    }
  ),
  
  ## Asparagopsis taxiformis --------------------------------------------------------------------------------------
  tar_target(a_taxiformis_file, file.path("data", "asparagopsis_taxiformis.qs"), format = "file"),
  tar_target(
    name = a_taxiformis,
    command = {
      vec <- qs::qread(a_taxiformis_file)
      vec <- vec[!is.na(vec)]
      vec <- vec[!names(vec) %in% c("M_am", "C_am", "V_ot", "K_ot", "M_ot", "C_ot")]
    }
  ),
  
  # Choose starting values ----------------------------------------------------------------------------------------
  
  tar_target(
    explore_starting_biomass,
    command = {
      biom_range <- seq(0.0005, 0.5, 0.0005) %>% 
        units::set_units("g L-1") %>% 
        units::set_units("mg m-3") %>% 
        units::drop_units()
      bf <- sapply(
        X = biom_range,
        FUN = biomass_to_Nf,
        Q_rel = 0.5,
        spec_params = a_armata,
        dry = T
      )
      init_states <- data.frame(Q_rel = 0.5, Nf = bf["Nf", ], Ns = bf["Ns", ])
      init_states$n = 1:nrow(init_states)
      
      res <- purrr::map(unique(init_states$n), function(nx) {
        init_state <- init_states[init_states$n == nx, ]
        df <- grow_macroalgae(
          t = 1:32,
          temperature = rep(a_armata["T_opt"], 32), 
          salinity = rep(a_armata["S_opt"], 32), 
          light = rep(a_armata["I_o"], 32), 
          kW = rep(0.0995, 32),
          velocity = rep(0.001, 32),
          nitrate = rep(40, 32), 
          ammonium = rep(6, 32), 
          ni_uptake = "linear",
          am_uptake = "MM",
          site_params = c(d_top = 1, hc = 5, farmA = 50 * 50, hz = 35, turbulence = NA),
          spec_params = a_armata,
          initials = c(Nf = init_state$Nf, Ns = init_state$Ns, Q_rel = init_state$Q_rel)
        ) %>% 
          as.data.frame() %>% 
          dplyr::select(Nf, Ns, t) %>% 
          mutate(biom_start = nx)
        merge(slice_min(df, t), slice_max(df, t), by = "biom_start") %>% 
          mutate(N_rem = 5*(Nf.y + Ns.y) - 5*(Nf.x + Ns.x)) %>% 
          dplyr::select(biom_start, N_rem)
      }) %>% bind_rows()
      
      merge(res, init_states, by.y = "n", by.x = "biom_start")
    }
  ),

  # Sensitivities -------------------------------------------------------------------------------------------------
  ## Standard values ----------------------------------------------------------------------------------------------
  tar_target(
    sens_conditions,
    data.frame(t_span = rep(1:732, 2)) %>% 
      mutate(
        Ni_input = 20.50 + 19.5 * sin((t_span * pi - 350) / 182.5),
        Am_input = 3.15 + 3 * sin((t_span * pi - 325) / 182.5),
        I_input  = 325 + 250 * sin((t_span * pi + 300) / 182.5),
        S_input  = 35.1 + 0.85 * sin((t_span * pi + 450) / 182.5),
        T_input  = c(26.75 + 3.5 * sin((1:732 * pi + 220) / 182.5), 
                     16.75 + 4.75 * sin((1:732 * pi + 120) / 182.5)),
        T_level  = rep(c(1,2), each = 732),
        U_input  = 0.3 + 0.25 * sin((t_span * pi + 180) / 182.5),
        kW       = Secchi_to_Kd(20)
      )
  ),
  tar_target(T_levels, c(1,2)),
  
  ## Culture conditions -------------------------------------------------------------------------------------------
  tar_target(name = starts, command = seq(5, 365, 15)),
  tar_target(name = culture_depths, command = c(0.5, 2.5, 5, 10)),
  tar_target(name = factors, command = c(0.95, 1, 1.05)),
  tar_target(name = param_names, command = names(a_armata)),

  # Sensitivity armata --------------------------------------------------------------------------------------------
  tar_target(
    name = sensitivity_run_arma,
    command = {
      init_state <- explore_starting_biomass[25,]
      init_state <- c(Nf = init_state$Nf, Ns = init_state$Ns, Q_rel = init_state$Q_rel)
      ydays <- starts:(starts + 31)
      
      purrr::map(factors, function(fct) {
        spec_params_arma <- adj_params(
          params = a_armata,
          focus_param = param_names,
          factor = fct
        )
        conds <- sens_conditions[sens_conditions$T_level == 1 & sens_conditions$t_span %in% ydays, ]
        grow_macroalgae(
          t = ydays,
          temperature = conds$T_input, 
          salinity = conds$S_input, 
          light = conds$I_input, 
          kW = rep(0.0995, length(ydays)),
          velocity = conds$U_input,
          nitrate = conds$Ni_input, 
          ammonium = conds$Am_input, 
          ni_uptake = "linear",
          am_uptake = "MM",
          site_params = c(d_top = culture_depths, hc = 5, farmA = 50 * 50, hz = 35, turbulence = NA),
          spec_params = spec_params_arma,
          initials = init_state
        ) %>% 
          as.data.frame() %>% 
          mutate(factor = as.factor(fct))
      }) %>% 
        bind_rows() %>% 
        pivot_longer(
          names_to = "output", names_transform = list(output = as.factor), 
          values_to = "value", values_transform = list(value = as.numeric), 
          cols = c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "T_lim", "S_lim", "I_top", "I_lim", "u_c")
        ) %>%
        mutate(
          start_t = as.integer(starts),
          adj_param = as.factor(param_names),
          cult_depth = as.factor(culture_depths)
        )
    },
    pattern = cross(starts, param_names, culture_depths)
  ),
  
  tar_target(
    name = sens_N_end_arma_0,
    command = {
      all <- sensitivity_run_arma %>% 
        filter(output %in% c("Nf", "Ns")) %>% 
        group_by(t, factor, start_t, adj_param, cult_depth) %>% 
        reframe(N = sum(value))
      start <- all %>% filter(t == start_t) %>% select(-t)
      end <- all %>% group_by(start_t) %>% slice_max(t) %>% ungroup() %>% select(-t)
      Nrem <- merge(start, end, by = c("start_t", "cult_depth", "adj_param", "factor")) %>% 
        mutate(N_rem = N.y - N.x) %>% 
        select(-c(N.y, N.x)) %>% 
        pivot_wider(names_from = factor, values_from = N_rem, names_prefix = "fact_") %>% 
        mutate(sens = (fact_1.05 - fact_0.95)/(0.1*fact_1),
               species = as.factor("armata")) %>% 
        select(-contains("fact"))
    },
    pattern = sensitivity_run_arma
  ),
  tar_target(
    sens_N_end_arma,
    sens_N_end_arma_0 %>% 
      group_by(species, cult_depth, adj_param) %>% 
      reframe(mean = meanna(sens),
              sd = sdna(sens))
  ),

  # Sensitivity taxiformis -----------------------------------------------------------------------------------------
  tar_target(
    name = sensitivity_run_taxi,
    command = {
      init_state <- explore_starting_biomass[25,]
      init_state <- c(Nf = init_state$Nf, Ns = init_state$Ns, Q_rel = init_state$Q_rel)
      ydays <- starts:(starts + 31)
      
      purrr::map(factors, function(fct) {
        spec_params_taxi <- adj_params(
          params = a_taxiformis,
          focus_param = param_names,
          factor = fct
        )
        conds <- sens_conditions[sens_conditions$T_level == 1 & sens_conditions$t_span %in% ydays, ]
        grow_macroalgae(
          t = ydays,
          temperature = conds$T_input, 
          salinity = conds$S_input, 
          light = conds$I_input, 
          kW = rep(0.0995, length(ydays)),
          velocity = conds$U_input,
          nitrate = conds$Ni_input, 
          ammonium = conds$Am_input, 
          ni_uptake = "linear",
          am_uptake = "MM",
          site_params = c(d_top = culture_depths, hc = 5, farmA = 50 * 50, hz = 35, turbulence = NA),
          spec_params = spec_params_taxi,
          initials = init_state
        ) %>% 
          as.data.frame() %>% 
          mutate(factor = as.factor(fct))
      }) %>% 
        bind_rows() %>% 
        pivot_longer(
          names_to = "output", names_transform = list(output = as.factor), 
          values_to = "value", values_transform = list(value = as.numeric), 
          cols = c("Nf", "Ns", "growth_rate", "Ns_to_Nf", "Ns_loss", "Nf_loss", "Q_int", "Q_rel", "Q_lim", "B_dw.mg", "B_ww.mg", "hm", "conc_nitrate", "up_Ni", "conc_ammonium", "up_Am", "T_lim", "S_lim", "I_top", "I_lim", "u_c")
        ) %>%
        mutate(
          start_t = as.integer(starts),
          adj_param = as.factor(param_names),
          cult_depth = as.factor(culture_depths)
        )
    },
    pattern = cross(starts, param_names, culture_depths)
  ),
  
  tar_target(
    name = sens_N_end_taxi_0,
    command = {
      all <- sensitivity_run_taxi %>% 
        filter(output %in% c("Nf", "Ns")) %>% 
        group_by(t, factor, start_t, adj_param, cult_depth) %>% 
        reframe(N = sum(value))
      start <- all %>% filter(t == start_t) %>% select(-t)
      end <- all %>% group_by(start_t) %>% slice_max(t) %>% ungroup() %>% select(-t)
      Nrem <- merge(start, end, by = c("start_t", "cult_depth", "adj_param", "factor")) %>% 
        mutate(N_rem = N.y - N.x) %>% 
        select(-c(N.y, N.x)) %>% 
        pivot_wider(names_from = factor, values_from = N_rem, names_prefix = "fact_") %>% 
        mutate(sens = (fact_1.05 - fact_0.95)/(0.1*fact_1),
               species = as.factor("taxiformis")) %>% 
        select(-contains("fact"))
    },
    pattern = sensitivity_run_taxi
  ),
  tar_target(
    sens_N_end_taxi,
    sens_N_end_taxi_0 %>% 
      group_by(species, cult_depth, adj_param) %>% 
      reframe(mean = meanna(sens),
              sd = sdna(sens))
  )
)

