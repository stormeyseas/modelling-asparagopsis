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

# Globals ---------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_NSW")
states <- tar_read(states)

# New South Wales -------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_NSW")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_NSW.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_NSW.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_NSW.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_NSW.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_NSW.parquet"))

# Tasmania --------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_TAS")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_TAS.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_TAS.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_TAS.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_TAS.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_TAS.parquet"))

# Victoria --------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_VIC")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_VIC.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_VIC.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_VIC.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_VIC.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_VIC.parquet"))

# Queensland ------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_QLD")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_QLD.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_QLD.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_QLD.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_QLD.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_QLD.parquet"))

# South Australia -------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_SAU")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_SAU.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_SAU.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_SAU.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_SAU.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_SAU.parquet"))

# Western Australia North -----------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_WAN")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_WAN.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_WAN.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_WAN.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_WAN.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_WAN.parquet"))

# Western Australia South -----------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_WAS")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_WAS.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_WAS.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_WAS.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_WAS.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_WAS.parquet"))

# Northern Territory ----------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_NTE")

tar_read(Ilim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_consec_NTE.parquet"))

tar_read(Tlim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -temperature) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_PL_NTE.parquet"))

tar_read(Tlim_cell_consec) %>% 
  mutate(state = factor(state, levels = states),
         species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Tlim_cell_consec_NTE.parquet"))

tar_read(Slim_cell_PL) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  dplyr::select(-t, -salinity) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_PL_NTE.parquet"))

tar_read(Slim_cell_consec) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path_2, "Slim_cell_consec_NTE.parquet"))

