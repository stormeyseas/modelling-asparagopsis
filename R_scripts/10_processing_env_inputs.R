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
  library(arrow)
  library(purrr)
  library(future)
  library(future.apply)
  library(furrr)
}))
# devtools::install_github("https://github.com/stormeyseas/macrogrow.git", quiet = T)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)

source("R_scripts/00_targets_functions.R")

# This script takes files from the targets pipeline and saves them as files, which are easier to load
script <- "R_scripts/09_env_inputs.R"
store <- "targets_outputs/_env_inputs"
out_path <- "data/processed_env_inputs"

tar_read(states_bbox, store = store) %>% 
  qsave(file.path(out_path, "states_bbox.qs"))

tar_read(BARRA_C2_cell_coords, store = store) %>% 
  mutate(state = as.factor(state)) %>% 
  write_parquet(file.path(out_path, "BARRA_C2_cell_coords.parquet"))

tar_read(BARRA_C2_cell_nos, store = store) %>% 
  qsave(file.path(out_path, "BARRA_C2_cell_nos.qs"))

BARRA_C2_cell_nos_chunked <- tar_read(BARRA_C2_cell_nos_chunked, store = store)

tar_read(state_input_timeseries, store = store) %>% 
  write_parquet(file.path(out_path, "state_input_timeseries.parquet"))

# # Cell input timeseries processing ------------------------------------------------------------------------------
for (i in seq_along(BARRA_C2_cell_nos_chunked)) {
  tar_read(cell_input_chunked_gapfilled, store = store, branches = i) %>%
    write_parquet(file.path(out_path, "cell_input_all", paste0("cell_input_all_", fixnum(i), ".parquet")))
}
cell_inputs_all <- file.path(out_path, "cell_input_all") %>%
  list.files(full.names = T) %>%
  purrr::map(., function(fnm) {
    arrow::read_parquet() %>% distinct(cell_no, hz)
  }) %>% 
  dplyr::bind_rows() %>%
  write_parquet(file.path(out_path, "cell_bathy.parquet"))

# N data prioritisation outputs -----------------------------------------------------------------------------------
stations <- list()
for (i in seq_along(BARRA_C2_cell_nos_chunked)) {
  # Get combined data
  N_data_combined <- tar_read(N_data_combined, store = store, branches = i) 
  N_data_combined %>%
    write_parquet(file.path(out_path, "combined_N_data", paste0("N_data_combined_", fixnum(i), ".parquet")))
  
  tar_read(N_data_prioritised, store = store, branches = i) %>%
    write_parquet(file.path(out_path, "combined_N_data", paste0("N_data_parameters_", fixnum(i), ".parquet")))
  
  # Count how many outfalls/refstations per cell
  stations[[i]] <- rbind(
    tar_read(N_outfall_data, store = store, branches = i) %>% 
      distinct(cell_no, data_source, name),
    tar_read(N_refstation_data, store = store, branches = i) %>% 
      rename(name = StationName) %>% 
      distinct(cell_no, data_source, name)
  )
}

total_cells <- tar_read(BARRA_C2_cell_nos, store = store) %>% length()

# rbenchmark::benchmark(
#   purrr = {purrr::list_rbind(stations)},
#   dplyr = {dplyr::bind_rows(stations)}
# )

stations <- stations %>% 
  dplyr::bind_rows() %>% 
  write_parquet(file.path(out_path, paste0("N_data_stations_byname.parquet"))) 

# How many refstations/outfalls are influencing each cell?
stations %>% 
  group_by(cell_no, data_source) %>% 
  reframe(num = n()) %>% 
  write_parquet(file.path(out_path, paste0("N_data_stations_bycell.parquet"))) 

# How many cells are being influenced by each refstations/outfall?
stations %>% 
  group_by(data_source, name) %>% 
  reframe(num = n(),
          total_cells = total_cells,
          prop = num/total_cells) %>% 
  write_parquet(file.path(out_path, paste0("N_data_cells_bystation.parquet"))) 


