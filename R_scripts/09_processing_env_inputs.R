# This script takes files from the env_inputs targets pipeline and saves them as files, which are easier to load

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

source("R_scripts/00_functions.R")

Sys.setenv(TAR_PROJECT = "project_inputs")
out_path <- here("data/processed_env_inputs")

tar_load_globals()

tar_load(states_bbox) 
states_bbox %>% 
  qsave(file.path(out_path, "states_bbox.qs"))

tar_load(BARRA_R2_cell_coords) 
BARRA_R2_cell_coords %>% 
  mutate(state = as.factor(state)) %>% 
  write_parquet(file.path(out_path, "BARRA_R2_cell_coords.parquet"))

tar_load(BARRA_R2_cell_nos) 
BARRA_R2_cell_nos %>% 
  qsave(file.path(out_path, "BARRA_R2_cell_nos.qs"))

tar_load(BARRA_R2_cell_nos_chunked)

# # Cell input timeseries processing ------------------------------------------------------------------------------
hz_list <- list()
for (i in seq_along(BARRA_R2_cell_nos_chunked)) {
  tar_load(cell_input_chunked_gapfilled, branches = i)

  hz_list[[i]] <- cell_input_chunked_gapfilled %>% 
    distinct(cell_no, hz)

  cell_input_chunked_gapfilled %>%
    dplyr::select(-hz) %>% 
    write_parquet(file.path(out_path, "cell_input_all", paste0("cell_input_all_", fixnum(i), ".parquet")))

  cell_input_chunked_gapfilled %>%
    dplyr::select(-hz) %>% 
    mutate(
      Ni_input = Ni_input + Ni_add,
      Am_input = Am_input + Am_add
    ) %>% 
    write_parquet(file.path(out_path, "cell_input_supp", paste0("cell_input_supp_", fixnum(i), ".parquet")))
}
hz_list %>% 
  bind_rows() %>% 
  write_parquet(file.path(out_path, "cell_bathy.parquet"))

tar_load(missing_data)
write_parquet(missing_data, file.path(out_path, "missing_data.parquet"))

# N data prioritisation outputs -----------------------------------------------------------------------------------
stations <- list()
for (i in seq_along(BARRA_R2_cell_nos_chunked)) {
  tar_read(N_data_prioritised, branches = i) %>%
    write_parquet(file.path(out_path, "combined_N_data", paste0("N_data_parameters_", fixnum(i), ".parquet")))
  
  print(paste("Branch", i, "of", length(BARRA_R2_cell_nos_chunked), "done"))
}
tar_load(N_data_combined)
N_data_combined %>%
  write_parquet(file.path(out_path, "combined_N_data", paste0("N_data_combined.parquet")))

total_cells <- tar_read(BARRA_R2_cell_nos) %>% length()

# rbenchmark::benchmark(
#   purrr = {purrr::list_rbind(stations)},
#   dplyr = {dplyr::bind_rows(stations)}
# )

# Count how many outfalls/refstations per cell
stations <- tar_read(N_data_combined) %>% 
  distinct(cell_no, data_source, name)

stations %>% 
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

# State input timeseries ------------------------------------------------------------------------------------------
tar_load(state_input_timeseries) 
state_input_timeseries %>% 
  write_parquet(file.path(out_path, "state_input_timeseries.parquet"))
