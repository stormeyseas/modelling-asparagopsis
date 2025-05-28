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
plan(multisession, workers = 18)

source("R_scripts/00_targets_functions.R")

# This script takes files from the targets pipeline and saves them as files, which are easier to load
script <- "R_scripts/09_env_inputs.R"
store <- "targets_outputs/_env_inputs"
out_path <- "data/processed_env_inputs"

states_bbox <- tar_read(states_bbox, store = store)
qsave(states_bbox, file.path(out_path, "states_bbox.qs"))

tar_read(BARRA_C2_cell_coords, store = store) %>% 
  mutate(state = as.factor(state)) %>% 
  write_parquet(file.path(out_path, "BARRA_C2_cell_coords.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos, store = store)
qsave(BARRA_C2_cell_nos, file.path(out_path, "BARRA_C2_cell_nos.qs"))

tar_read(state_input_timeseries, store = store) %>% 
  write_parquet(file.path(out_path, "state_input_timeseries.parquet"))

# # Cell input timeseries processing
input_branches <- seq(0, length(BARRA_C2_cell_nos), length.out = 101) %>% as.integer()
for (i in 1:(length(input_branches)-1)) {
  tar_read(cell_input_gapfilled, store = store, branches = (input_branches[i]+1):input_branches[i+1]) %>%
    write_parquet(file.path(out_path, paste0("cell_input_all_", fixnum(i), ".parquet")))
}
cell_inputs_all <- out_path %>%
  list.files(full.names = T) %>%
  str_subset("cell_input_all") %>%
  purrr::map_dfr(arrow::read_parquet)
hz <- cell_inputs_all %>%
  distinct(cell_no, hz) %>%
  write_parquet(file.path(out_path, "cell_bathy.parquet"))

# N data prioritisation outputs
input_branches <- seq(0, length(BARRA_C2_cell_nos), length.out = 101) %>% as.integer()
for (i in 2:length(input_branches)) {
  ls <- tar_read(N_data_prioritised, store = store, branches = (input_branches[i-1]+1):input_branches[i])
  ls2 <- list()
  for (j in 1:length(ls)) {
    ls1 <- ls[[j]]
    cell_no <- unique(ls1$data$cell_no)
    ls2[[j]] <- list(
      data = ls1$data,
      mean_sd = data.frame(
        cell_no = cell_no,
        measure = names(ls1$mean),
        mean_value = as.numeric(ls1$mean),
        sd_value = as.numeric(ls1$sd),
        stringsAsFactors = FALSE
      ),
      parameters = data.frame(
        cell_no = cell_no,
        measure = names(ls1$a),
        a = as.numeric(ls1$a),
        b = as.numeric(ls1$b),
        c = as.integer(ls1$c),
        stringsAsFactors = FALSE
      )
    )
  }
  combined_data <- bind_rows(map(ls2, ~ .x$data))
  combined_mean_sd <- bind_rows(map(ls2, ~ .x$mean_sd))
  combined_parameters <- bind_rows(map(ls2, ~ .x$parameters))
  
  write_parquet(combined_data, file.path(out_path, paste0("combined_N_data_", fixnum(i-1, 4), ".parquet")))
  write_parquet(combined_mean_sd, file.path(out_path, paste0("combined_N_mean_sd_", fixnum(i-1, 4), ".parquet")))
  write_parquet(combined_parameters, file.path(out_path, paste0("combined_N_parameters_", fixnum(i-1, 4), ".parquet")))
}

# Station/outfall influence ---------------------------------------------------------------------------------------
# Count how many refstations and outfalls are influencing each cell
total_cells <- (1:length(BARRA_C2_cell_nos))
input_branches <- split(total_cells, ceiling(seq_along(total_cells)/250))

by_refstation <- future_map_dfr(input_branches, function(br) {
  tar_read(N_refstation_data, store = store, branches = br) %>% 
    mutate(StationName = as.factor(StationName)) %>% 
    distinct(cell_no, StationName) %>% 
    group_by(cell_no, StationName) %>% 
    reframe(n = n())
}, .progress = T, .options = furrr_options(seed = TRUE))

by_refstation %>% 
  group_by(StationName) %>% 
  reframe(nc = sum(n)) %>% 
  write_parquet(file.path(out_path, "cells_by_station_refstations.parquet"))

by_refstation %>% 
  group_by(cell_no) %>% 
  reframe(ns = sum(n)) %>% 
  write_parquet(file.path(out_path, "stations_by_cell_refstations.parquet"))


by_outfall <- future_map_dfr(input_branches, function(br) {
  tar_read(N_outfall_data, store = store, branches = br) %>% 
    mutate(name = as.factor(name)) %>% 
    distinct(cell_no, name) %>% 
    group_by(cell_no, name) %>% 
    reframe(n = n())
}, .progress = T, .options = furrr_options(seed = TRUE))
by_outfall %>% 
  group_by(name) %>% 
  reframe(nc = sum(n)) %>% 
  write_parquet(file.path(out_path, "cells_by_station_outfalls.parquet"))

by_outfall %>% 
  group_by(cell_no) %>% 
  reframe(ns = sum(n)) %>% 
  write_parquet(file.path(out_path, "stations_by_cell_outfalls.parquet"))

rm(by_refstation, by_outfall)

