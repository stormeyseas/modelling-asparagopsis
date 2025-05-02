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

BARRA_C2_cell_coords <- tar_read(BARRA_C2_cell_coords) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "BARRA_C2_cell_coords.parquet")) %>% 
  group_by(state) %>% 
  reframe(num_cells = n())

# New South Wales -------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_NSW")

tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_NSW.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_NSW.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_NSW.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_NSW.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_NSW.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_NSW.parquet"))

# Tasmania --------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_TAS")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_TAS.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_TAS.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_TAS.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_TAS.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_TAS.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_TAS.parquet"))

# Victoria --------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_VIC")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_VIC.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_VIC.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_VIC.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_VIC.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_VIC.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_VIC.parquet"))

# Queensland ------------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_QLD")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_QLD.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_QLD.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_QLD.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_QLD.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_QLD.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_QLD.parquet"))

# South Australia -------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_SAU")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_SAU.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_SAU.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_SAU.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_SAU.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_SAU.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_SAU.parquet"))

# Western Australia North -----------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_WAN")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_WAN.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_WAN.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_WAN.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_WAN.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_WAN.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_WAN.parquet"))

# Western Australia South -----------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_WAS")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_WAS.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_WAS.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_WAS.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_WAS.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_WAS.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_WAS.parquet"))

# Northern Territory ----------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_model_running_NTE")

cell_input_timeseries <- tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = factor(state, levels = states)) %>% 
  arrow::write_parquet(file.path(out_path, "cell_input_timeseries_NTE.parquet")) %>% 
  group_by(state, yday) %>% 
  reframe(
    I_input_mean = meanna(I_input), I_input_min = minna(I_input), I_input_max = maxna(I_input), I_input_sd = sdna(I_input),
    T_input_mean = meanna(T_input), T_input_min = minna(T_input), T_input_max = maxna(T_input), T_input_sd = sdna(T_input),
    S_input_mean = meanna(S_input), S_input_min = minna(S_input), S_input_max = maxna(S_input), S_input_sd = sdna(S_input),
    UV_input_mean = meanna(UV_input), UV_input_min = minna(UV_input), UV_input_max = maxna(UV_input), UV_input_sd = sdna(UV_input)
  ) %>% 
  arrow::write_parquet(file.path(out_path, "state_input_timeseries_NTE.parquet"))

BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos)
refstation_cell_match <- tar_read(refstation_cell_match)
refstation_cell_match$cell_no <- rep(BARRA_C2_cell_nos, each = 2)
arrow::write_parquet(refstation_cell_match, file.path(out_path, "refstation_cell_match_NTE.parquet"))
rm(refstation_cell_match)

tar_read(outfall_sites_cellpaired) %>% 
  arrow::write_parquet(file.path(out_path, "outfall_sites_cellpaired_NTE.parquet"))

tar_read(Ilim_cell_PL) %>% 
  mutate(state = factor(state, levels = states),
         depth = as.numeric(depth)) %>% 
  arrow::write_parquet(file.path(out_path, "Ilim_cell_PL_NTE.parquet")) %>% 
  arrow::write_parquet(file.path(out_path_2, "Ilim_cell_PL_NTE.parquet"))
