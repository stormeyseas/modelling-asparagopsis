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
  library(streamMetabolizer)
}))
# devtools::install_github("https://github.com/stormeyseas/macrogrow.git", quiet = T)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)

# This script takes files from the targets pipeline and saves them as files, which are easier to load
out_path <- file.path("data", "processed_model_running")

# Sys.setenv(TAR_PROJECT = "project_model_running")
# meta <- tar_meta(fields = c("names", "parent", "time")) %>%
#   mutate(name = case_when(is.na(parent) ~ name,
#                           T ~ sub("_[a-z0-9]+$", "", name)),
#          time = as.Date(time)) 
# meta <- meta %>% 
#   group_by(name) %>%
#   reframe(time = max(time, na.rm = T)) #%>% filter(time == as.Date(lubridate::now()))

this_state <- tar_read(this_state)

# Cells pathway ---------------------------------------------------------------------------------------------------
tar_read(BARRA_C2_cell_coords) %>% arrow::write_parquet(file.path(out_path, "BARRA_C2_cell_coords.parquet"))

# Cell inputs -----------------------------------------------------------------------------------------------------
# All inputs
tar_read(cell_input_timeseries) %>% 
  filter(yday <= 365) %>% 
  mutate(state = this_state) %>% 
  arrow::write_parquet(file.path(out_path, paste0("cell_input_timeseries_", this_state, ".parquet"))) %>% 
  group_by(state, yday) %>% 
  reframe(I_input = mean(I_input, na.rm = T),
          T_input = mean(T_input, na.rm = T),
          S_input = mean(S_input, na.rm = T),
          UV_input = mean(UV_input, na.rm = T)) %>% 
  arrow::write_parquet(file.path(out_path, paste0("state_input_timeseries_", this_state, ".parquet")))

# Partial lims ----------------------------------------------------------------------------------------------------
tar_read(N_loss_PL) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path, paste0("N_loss_PL_", this_state, ".parquet")))

tar_read(Tlim_cell_PL) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path, paste0("Tlim_cell_PL_", this_state, ".parquet")))

tar_read(Slim_cell_PL) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path, paste0("Slim_cell_PL_", this_state, ".parquet")))

tar_read(Ilim_cell_PL) %>% 
  mutate(species = as.factor(species)) %>% 
  arrow::write_parquet(file.path(out_path, paste0("Ilim_cell_PL_", this_state, ".parquet")))




