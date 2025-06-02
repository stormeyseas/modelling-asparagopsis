# nolint start

# Setup -----------------------------------------------------------------------------------------------------------
library(targets)
library(tarchetypes)
library(crew)
library(tidyr)
library(conflicted)
library(qs)
library(qs2)
library(dplyr)
library(macrogrow)
library(geosphere)
library(stringr)
library(magrittr)
library(stringr)
library(geotargets)
library(raster)
library(sp)
library(units)
library(furrr)
library(tictoc)
library(arrow)
library(purrr)
library(lubridate)
library(terra)
conflicts_prefer(dplyr::filter(), dplyr::select(), magrittr::`%>%`, .quiet = T)
remove_unit("molN")
install_unit(symbol = "molN", def = "14.007 g")

# Load functions and global variables ---------------------------------------------------------------------------
source(file.path("R_scripts", "00_targets_functions.R"))
extract_yday <- function(filename) {
  yday_match <- regexpr("yday\\d+", filename) 
  if (yday_match != -1) {
    return(regmatches(filename, yday_match))
  } else {
    return(NA)
  }
}
extract_chunk <- function(filename) {
  chunk_match <- regexpr("chunk\\d+", filename) 
  if (chunk_match != -1) {
    return(regmatches(filename, chunk_match))
  } else {
    return(NA)
  }
}
Ni_add <- 3.5 %>% units::set_units("umolN L-1") %>% units::set_units("mg m-3") %>% units::drop_units()
Am_add <- 6.5 %>% units::set_units("umolN L-1") %>% units::set_units("mg m-3") %>% units::drop_units()

# Load data -----------------------------------------------------------------------------------------------------
# Check global environment first, then load from files
if(!exists("species_data")) {
  species_data <- list(
    tar_read(a_armata_gametophyte, store = "targets_outputs/_species"),
    tar_read(a_taxiformis_gametophyte, store = "targets_outputs/_species")
  )
}

if(!exists("BARRA_C2_cell_coords")) {
  BARRA_C2_cell_coords <- arrow::read_parquet("data/processed_env_inputs/BARRA_C2_cell_coords.parquet")
}

if(!exists("BARRA_C2_cell_nos")) {
  BARRA_C2_cell_nos <- qs::qread("data/processed_env_inputs/BARRA_C2_cell_nos.qs")
}

if(!exists("cell_inputs")) {
  cell_inputs_files <- list.files("data/processed_env_inputs", full.names = T) %>% str_subset("cell_input_all")
  cell_inputs <- cell_inputs_files %>% 
    purrr::map_dfr(arrow::read_parquet) %>% 
    mutate(state = as.factor(state)) %>% 
    dplyr::filter(yday <= 400)
  
  # Isolate cell bathymetry (static)
  bathy <- cell_inputs %>% 
    distinct(cell_no, hz) %>% 
    # After a certain point, depth is irrelevant
    mutate(hz = case_when(hz > 100 ~ 100, T ~ hz))
  
  # Find cells without bathymetry (don't know how that happened)
  idx <- which(is.na(bathy$hz))
  for (id in idx) {
    bathy$hz[id] <- meanna(c(bathy$hz[id-1], bathy$hz[id+1]))
  }
  idx <- which(is.na(bathy$hz)) # check again
  rm(idx)
  # hist(bathy$hz)
  
  # Exclude cells that won't grow because they're too shallow
  shallow_cells <- bathy$cell_no[bathy$hz < 8]
  cell_inputs <- cell_inputs[!cell_inputs$cell_no %in% shallow_cells, ] %>% dplyr::select(-hz)
  bathy <- bathy[!bathy$cell_no %in% shallow_cells, ]
  BARRA_C2_cell_nos <- BARRA_C2_cell_nos[!BARRA_C2_cell_nos %in% shallow_cells]
  gc()
}

if(!exists("state_inputs")) {
  state_inputs <- arrow::read_parquet("data/processed_env_inputs/state_input_timeseries.parquet") %>% 
    dplyr::filter(yday <= 400) %>% 
    mutate(state = as.factor(state))
  
  # Replace cell inputs with state inputs if needed
  prog <- as.integer(seq(0, length(unique(cell_inputs$cell_no)), length.out = 101))
  for (i in seq_along(unique(cell_inputs$cell_no))) {
    cn <- unique(cell_inputs$cell_no)[i]
    
    salinity <- cell_inputs$S_input[cell_inputs$cell_no == cn]
    velocity <- cell_inputs$UV_input[cell_inputs$cell_no == cn]
    kW <- cell_inputs$Kd_490[cell_inputs$cell_no == cn]
    temperature <- cell_inputs$T_input[cell_inputs$cell_no == cn]
    light <- cell_inputs$I_input[cell_inputs$cell_no == cn]
    nitrate <- cell_inputs$Ni_input[cell_inputs$cell_no == cn]
    ammonium <- cell_inputs$Am_input[cell_inputs$cell_no == cn]
    
    # Check for NA values in data from BRAN
    if (any(is.na(salinity))) {
      st <- unique(cell_inputs$state[cell_inputs$cell_no == cn])
      cell_inputs$S_input[cell_inputs$cell_no == cn] <- 
        state_inputs$S_input_mean[state_inputs$state == st]
    }
    if (any(is.na(velocity))) {
      cell_inputs$UV_input[cell_inputs$cell_no == cn] <- 
        state_inputs$UV_input_mean[state_inputs$state == st]
    }
    if (any(is.na(kW))) {
      cell_inputs$Kd_490[cell_inputs$cell_no == cn] <- 
        state_inputs$Kd_490_mean[state_inputs$state == st]
    }
    if (any(is.na(temperature))) {
      cell_inputs$T_input[cell_inputs$cell_no == cn] <- 
        state_inputs$T_input_mean[state_inputs$state == st]
    }
    if (any(is.na(light))) {
      cell_inputs$I_input[cell_inputs$cell_no == cn] <- 
        state_inputs$I_input_mean[state_inputs$state == st]
    }
    if (any(is.na(nitrate))) {
      cell_inputs$Ni_input[cell_inputs$cell_no == cn] <- 
        state_inputs$Ni_input_mean[state_inputs$state == st]
    }
    if (any(is.na(ammonium))) {
      cell_inputs$Am_input[cell_inputs$cell_no == cn] <- 
        state_inputs$Am_input_mean[state_inputs$state == st]
    }
    if (i %in% prog) {
      print(str_c("Finished ", i, " of ", length(unique(cell_inputs$cell_no)), " cells, ", 
                  round(i/length(unique(cell_inputs$cell_no)), 2), "% done"))
    }
  }
}

# Initialize parameters -------------------------------------------------------------------------------------------
init_biomass <- 0.005 %>% units::set_units("g L-1") %>% 
  units::set_units("mg m-3") %>% units::drop_units()

init_states <- purrr::map(species_data, ~{
  c(Nf = biomass_to_Nf(biomass = init_biomass, Q_rel = 0.5, spec_params = unlist(.x), dry = T),
    Q_rel = 0.5)
})

start_dates_growth <- lubridate::ceiling_date(lubridate::make_date(2023, 1:12, 15), 'month') - 
  lubridate::days(42+1)
start_ydays_growth <- lubridate::yday(start_dates_growth)

# Run full growth -------------------------------------------------------------------------------------------------
overwrite <- F

# This needs to be done in chunks or else it's too large
# Define chunks
chunks <- split(BARRA_C2_cell_nos, ceiling(seq_along(BARRA_C2_cell_nos)/1000))
all_site_params <- c(
  hc = 5,
  farmA = 50 * 50,
  turbulence = NA,
  d_top = 2.5,
  hz = NA
)

start_ydays_growth <- start_ydays_growth[6:7]

for(sp in 1:2) {
  for (ch in seq_along(chunks)) {
    # Restrict data to cells in chunk
    chunk_cells <- chunks[[ch]]
    # Pre-process input data by cell
    chunk_cell_inputs <- cell_inputs[cell_inputs$cell_no %in% chunk_cells, ]
    chunk_state_inputs <- state_inputs[state_inputs$state %in% unique(chunk_cell_inputs$state), ]
    chunk_cell_inputs <- split(chunk_cell_inputs, chunk_cell_inputs$cell_no)
    chunk_bathy <- bathy[bathy$cell_no %in% chunk_cells, ]
    
    # Cell growth calculations
    for (start_yday in start_ydays_growth) {
      fnms <- c(
        file.path("data", "processed_model_running", "chunks", str_c("cell_growth_lims_sp", sp, "_chunk", fixnum(ch,3), "_yday", fixnum(start_yday,3), ".qs")),
        file.path("data", "processed_model_running", "chunks", str_c("cell_growth_end_TN_sp", sp, "_chunk", fixnum(ch,3), "_yday", fixnum(start_yday,3), ".qs"))
      )
      
      tic()
      if (any(!file.exists(fnms)) | overwrite) {
        # Get results for all cells in the chunk for 1 start day
        results_1_month <- list()
        for (i in 1:length(chunk_cells)) {
          cell_no <- chunk_cells[i]
          cell_input <- chunk_cell_inputs[[i]]
          site_params <- all_site_params
          start <- start_yday
          state_input <- chunk_state_inputs[chunk_state_inputs$state == unique(cell_input$state), ]
          site_params['hz'] <- chunk_bathy$hz[chunk_bathy$cell_no == cell_no]

          mat <- grow_macroalgae(
            start =        start,
            grow_days =    42,
            temperature =  cell_input$T_input[start:(start + 42)],
            salinity =     cell_input$S_input[start:(start + 42)],
            light =        cell_input$I_input[start:(start + 42)],
            kW =           cell_input$Kd_490[start:(start + 42)],
            velocity =     cell_input$UV_input[start:(start + 42)],
            nitrate =      cell_input$Ni_input[start:(start + 42)],
            ammonium =     cell_input$Am_input[start:(start + 42)],
            ni_uptake =    "linear",
            am_uptake =    "MM",
            site_params =  site_params,
            spec_params =  species_data[[sp]],
            initials =     init_states[[sp]]
          )
          
          # Add cell number and total N to growth matrix
          lims <- get_limiting(mat, spec_params = unlist(species_data[[sp]]))
          
          # Output results as data frame
          results_1_month[[i]] <- mat %>%
            as.data.frame() %>%
            mutate(cell_no = cell_no,
                   lims = lims) %>% 
            select(cell_no, t, Nf, Ns, lims)
        }
        
        results_1_month <- bind_rows(results_1_month) %>% 
          mutate(start = start_yday,
                 TN = (Nf+Ns)*site_params['hc']) %>% 
          dplyr::select(-c(Nf, Ns))
        
        cell_growth_lims <- results_1_month %>%
          dplyr::select(-TN) %>%
          mutate(lims = as.factor(lims))
        
        cell_growth_end_TN <- results_1_month %>%
          dplyr::select(-lims) %>%
          group_by(cell_no) %>%
          slice_tail(n = 1) %>%
          ungroup()
        
        # Save intermediate results
        qs::qsave(cell_growth_lims, fnms[1])
        qs::qsave(cell_growth_end_TN, fnms[2])
      }
      t <- toc(quiet = T)
      print(str_c("Finished chunk ", ch, ", species ", sp, ", yday ", start_yday, ", ", t$callback_msg))
    }
  }
}

# Theoretical scenarios -------------------------------------------------------------------------------------------
overwrite <- T

# Define scenarios
theo_scens <- tibble::tribble(
  ~name,            ~num,   ~T_mod, ~bathy, ~kW_mod, ~S_mod, ~Ni_add,  ~Am_add,  ~UV_mod,
  "base",           1,      1,      50,     1,       1,      0,        0,        1,
  "fish_farm",      2,      1,      30,     1,       1,      Ni_add,   Am_add,   1,
  "deep_water",     3,      0.9,    100,    1.2,     1,      0,        0,        4,
  "shallow_bay",    4,      1.1,    15,     0.8,     1.1,    Ni_add,   0,        0.25,
  "estuary_mouth",  5,      1,      30,     0.8,     0.9,    Ni_add,   0,        1
)
states <- levels(state_inputs$state)

site_params <- c(
  hc = 5,
  farmA = 50 * 50,
  turbulence = NA,
  d_top = 2.5,
  hz = 40
)

for(sp in 1:2) {
  # Cell growth calculations
  for (start_yday in start_ydays_growth) {
    fnms <- c(
      file.path("data", "processed_model_running", str_c("scen_growth_lims_sp", sp, "_yday", fixnum(start_yday,3), ".qs")),
      file.path("data", "processed_model_running", str_c("scen_growth_end_TN_sp", sp, "_yday", fixnum(start_yday,3), ".qs"))
    )
    
    tic()
    if (any(!file.exists(fnms)) | overwrite) {
      # Get results for all states in the chunk for 1 start day
      results_1_month_1 <- results_1_month_2 <- list()
      for (i in 1:nrow(theo_scens)) {
        for (j in 1:length(states)) {
          input <- state_inputs %>% 
            filter(state == states[j]) %>% 
            mutate(
              T_input_mean = T_input_mean * theo_scens$T_mod[theo_scens$num == i],
              S_input_mean = S_input_mean * theo_scens$S_mod[theo_scens$num == i],
              Kd_490_mean = Kd_490_mean * theo_scens$kW_mod[theo_scens$num == i],
              UV_input_mean = UV_input_mean * theo_scens$UV_mod[theo_scens$num == i],
              Ni_input_mean = Ni_input_mean + theo_scens$Ni_add[theo_scens$num == i],
              Am_input_mean = Am_input_mean + theo_scens$Am_add[theo_scens$num == i]
              )
          
          mat <- grow_macroalgae(
            start =        start_yday,
            grow_days =    42,
            temperature =  input$T_input_mean[start_yday:(start_yday + 42)],
            salinity =     input$S_input_mean[start_yday:(start_yday + 42)],
            light =        input$I_input_mean[start_yday:(start_yday + 42)],
            kW =           input$Kd_490_mean[start_yday:(start_yday + 42)],
            velocity =     input$UV_input_mean[start_yday:(start_yday + 42)],
            nitrate =      input$Ni_input_mean[start_yday:(start_yday + 42)],
            ammonium =     input$Am_input_mean[start_yday:(start_yday + 42)],
            ni_uptake =    "linear",
            am_uptake =    "MM",
            site_params =  site_params,
            spec_params =  species_data[[sp]],
            initials =     init_states[[sp]]
          )
          mat <- cbind(
            mat,
            matrix(data = mat[, 'Ns_loss'] + mat[, 'Nf_loss'], nrow = nrow(mat), ncol = 1),
            matrix(data = (mat[, 'Nf'] + mat[, 'Ns']) * 5, nrow = nrow(mat), ncol = 1),
            matrix(data = j, nrow = nrow(mat), ncol = 1),
            matrix(data = i, nrow = nrow(mat), ncol = 1)
          )
          colnames(mat) <- c(colnames(mat)[1:24], "N_loss", "TN", "state", "scen")
          lims <- get_limiting(mat, spec_params = species_data[[sp]])
          
          results_1_month_1[[j]] <- mat %>%
            as.data.frame() %>%
            select(state, scen, t, TN) %>%
            mutate(lims = lims,
                   start = start_yday)
        }
        results_1_month_2[[i]] <- bind_rows(results_1_month_1)
      }
      results_1_month <- bind_rows(results_1_month_2) 
      results_1_month <- results_1_month %>% 
        mutate(state = factor(state, labels = states),
               scen = factor(scen, labels = theo_scens$name))

      scen_growth_lims <- results_1_month %>%
        dplyr::select(-TN) %>%
        mutate(lims = as.factor(lims))
      
      scen_growth_end_TN <- results_1_month %>%
        dplyr::select(-lims) %>%
        group_by(state, scen, start) %>%
        slice_tail(n = 1) %>%
        ungroup()
      
      # Save intermediate results
      qs::qsave(scen_growth_lims, fnms[1])
      qs::qsave(scen_growth_end_TN, fnms[2])
    }
    t <- toc(quiet = T)
    print(str_c("Finished species ", sp, ", yday ", start_yday, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  }
}

# Get ready for further analysis ---------------------------------------------------------------------------------
## A. armata -----------------------------------------------------------------------------------------------------
# Cell growth
file_names <- file.path("data", "processed_model_running", "chunks") %>% 
  list.files(full.names = T) %>% str_subset("sp1") %>% str_subset("cell_growth")

# Subset by files by chunk (all times)
result_list <- split(file_names, sapply(file_names, extract_chunk))
for (i in seq_along(result_list)) {
  cell_growth_lims <- result_list[[1]] %>% 
    str_subset("lims") %>% 
    purrr::map_dfr(qs::qread) %>% 
    merge(BARRA_C2_cell_coords, by = "cell_no") %>% 
    dplyr::select(-layer)
  
  cell_growth_endTN <- result_list[[1]] %>% 
    str_subset("end_TN") %>% 
    purrr::map_dfr(qs::qread) %>% 
    merge(BARRA_C2_cell_coords, by = "cell_no") %>% 
    dplyr::select(-layer)
}

# Scenario growth
file_names <- file.path("data", "processed_model_running") %>% 
  list.files(full.names = T) %>% str_subset("sp1") %>% str_subset("scen_growth")
scen_growth_lims <- file_names %>% 
  str_subset("lims") %>% 
  purrr::map_dfr(qs::qread)

scen_growth_lims %>% 
  group_by(state, scen, lims) %>% 
  reframe(n = n()) %>% 
  ggplot(aes(x = lims, y = n, fill = scen)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~state)

## A. taxiformis -------------------------------------------------------------------------------------------------
file_names <- file.path("data", "processed_model_running", "chunks") %>% 
  list.files(full.names = T) %>% str_subset("sp2") %>% str_subset("cell_growth")

# nolint end