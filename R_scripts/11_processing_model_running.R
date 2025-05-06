suppressMessages(suppressWarnings({
  library(targets)
  library(tarchetypes)
  library(crew)
  library(tidyr)
  library(dplyr)
  library(conflicted)
  library(units)
  library(qs)
  library(qs2)
  library(magrittr)
  library(stringr)
  library(geotargets)
  library(raster)
  library(sp)
  library(units)
  library(tictoc)
  library(streamMetabolizer)
}))
# devtools::install_github("https://github.com/stormeyseas/macrogrow.git", quiet = T)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)

source("R_scripts/00_targets_functions.R")

# This script takes files from the targets pipeline and saves them as files, which are easier to load
out_path <- file.path("data", "processed_cell_data")
out_path_2 <- file.path("data", "processed_model_running")

projects <- tibble::tribble(
  ~state, ~st,
  "QLD",  "targets_outputs/_model_running_QLD",
  "TAS",  "targets_outputs/_model_running_TAS",  
  "SAU",  "targets_outputs/_model_running_SAU",  
  "WAS",  "targets_outputs/_model_running_WAS",  
  "NSW",  "targets_outputs/_model_running_NSW",  
  "VIC",  "targets_outputs/_model_running_VIC",  
  "WAN",  "targets_outputs/_model_running_WAN",  
  "NTE",  "targets_outputs/_model_running_NTE"
)
sc <- "R_scripts/09_model_running.R"
# mani <- tar_manifest(script = sc)
# meta <- tar_meta(store = projects$st[6])

# Globals ---------------------------------------------------------------------------------------------------------
states <- tar_read(states, store = projects$st[1])
theo_scens <- tar_read(theo_scens, store = projects$st[1]) %>% 
  arrow::write_parquet(file.path(out_path_2, "theo_scens.parquet"))

tar_read(BARRA_C2_cell_coords, store = projects$st[1]) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "BARRA_C2_cell_coords.parquet"))

tar_read(refstation_coords, store = projects$st[1]) %>% 
  arrow::write_parquet(file.path(out_path, "refstation_coords.parquet"))

tar_read(outfall_locations, store = projects$st[1]) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_locations.parquet"))

# State targets ---------------------------------------------------------------------------------------------------
for (i in 1:nrow(projects)) {
  # Cell input data
  tar_read(cell_input_timeseries, store = projects$st[i]) %>% 
    dplyr::filter(yday <= 365) %>% 
    mutate(state = factor(state, levels = states)) %>% 
    arrow::write_parquet(file.path(out_path, paste0("cell_input_timeseries_", projects$state[i], ".parquet")))

  tar_read(state_input_timeseries, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("state_input_timeseries_", projects$state[i], ".parquet")))
  
  BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos, store = projects$st[i])
  refstation_cell_match <- tar_read(refstation_cell_match, store = projects$st[i])
  refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
  arrow::write_parquet(refstation_cell_match, file.path(out_path, paste0("refstation_cell_match_", projects$state[i], ".parquet")))
  rm(refstation_cell_match)
  
  tar_read(outfall_sites_cellpaired, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("outfall_sites_cellpaired_", projects$state[i], ".parquet")))

  tar_read(total_cell_growth, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("total_cell_growth_", projects$state[i], ".parquet")))

  tar_read(theo_scen_growth, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("total_cell_growth_", projects$state[i], ".parquet")))
  
  tar_read(total_cell_growth_end_TN, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("total_cell_growth_end_TN_", projects$state[i], ".parquet")))
  
  tar_read(total_theo_growth_end_TN, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("total_theo_growth_end_TN_", projects$state[i], ".parquet")))
  
  tar_read(cell_growth_lims, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("cell_growth_lims_", projects$state[i], ".parquet")))
  
  tar_read(theo_growth_lims, store = projects$st[i]) %>% 
    arrow::write_parquet(file.path(out_path, paste0("theo_growth_lims_", projects$state[i], ".parquet")))
}

