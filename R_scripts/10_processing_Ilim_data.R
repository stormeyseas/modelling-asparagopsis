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
out_path <- file.path("data", "processed_model_running")

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

states <- tar_read(states, store = proj$st[1])

for (i in 1:nrow(projects)) {
  tar_read(Ilim_cell_PL, store = proj$st[i]) %>% 
    mutate(state = factor(state, levels = states),
           depth = as.numeric(depth)) %>% 
    arrow::write_parquet(file.path(out_path, paste0("Ilim_cell_PL_", proj$state[i], ".parquet")))
}

