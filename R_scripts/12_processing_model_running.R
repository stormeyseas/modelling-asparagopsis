# Setup -----------------------------------------------------------------------------------------------------------
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
library(here)
library(future)
library(future.apply)
library(furrr)
# devtools::install_github("https://github.com/stormeyseas/macrogrow.git", quiet = T)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)
source("R_scripts/00_targets_functions.R")

Sys.setenv(TAR_PROJECT = "project_running")
out_path <- here() %>% file.path("data", "processed_model_running")
mani <- tar_manifest() %>% 
  select(-command)

# General cell info -----------------------------------------------------------------------------------------------
BARRA_C2_cell_nos <- tar_read(BARRA_C2_cell_nos) 
BARRA_C2_cell_nos %>% 
  qsave(file.path(out_path, "cells_included.qs"))

tar_read(cells_to_omit) %>% 
  qsave(file.path(out_path, "cells_omitted.qs"))

BARRA_C2_cell_coords <- tar_read(BARRA_C2_cell_coords) %>% 
  filter(cell_no %in% BARRA_C2_cell_nos)
BARRA_C2_cell_coords %>% 
  qsave(file.path(out_path, "BARRA_C2_cell_coords.parquet"))
  
# Cell input timeseries processing
BARRA_C2_cell_nos_chunked <- tar_read(BARRA_C2_cell_nos_chunked)
# for (i in seq_along(BARRA_C2_cell_nos_chunked)) {
#   tar_read(cell_inputs_chunked, branches = i)[[1]] %>% 
#     unname()
# } %>% dplyr::bind_rows()
# write_parquet(df, file.path(out_path, "cell_inputs", paste0("cell_inputs_chunked_", fixnum(i), ".parquet")))
gc()

# Cell growth -----------------------------------------------------------------------------------------------------
cell_branches <- length(tar_read(start_ydays_growth)) * length(tar_read(BARRA_C2_cell_nos_chunked))
growth_chunks <- seq(0, cell_branches, length.out = 101) %>% as.integer()

for (i in 1:(length(growth_chunks)-1)) {
  stt <- growth_chunks[i]+1
  end <- growth_chunks[i+1]
  
  # Normal growth
  rbind(
    tar_read(cell_growth_armata, branches = stt:end),
    tar_read(cell_growth_taxiformis, branches = stt:end)
  ) %>% 
    mutate(t = as.integer(t),
           start = as.integer(start),
           cell_no = as.integer(cell_no)) %>% 
    write_parquet(file.path(out_path, "cell_growth", paste0("cell_growth_", fixnum(i), ".parquet"))) %>% 
    dplyr::select(t, growth_rate, start, species, cell_no, lim) %>% 
    write_parquet(file.path(out_path, "cell_growth", paste0("cell_growth_abbrev_", fixnum(i), ".parquet")))
  
  # Supplemented growth
  rbind(
    tar_read(cell_growth_supp_armata, branches = stt:end),
    tar_read(cell_growth_supp_taxiformis, branches = stt:end)
  ) %>% 
    mutate(t = as.integer(t),
           start = as.integer(start),
           cell_no = as.integer(cell_no)) %>% 
    write_parquet(file.path(out_path, "cell_growth", paste0("cell_growth_supp_", fixnum(i), ".parquet"))) %>% 
    dplyr::select(t, growth_rate, start, species, cell_no, lim) %>% 
    write_parquet(file.path(out_path, "cell_growth", paste0("cell_growth_supp_abbrev_", fixnum(i), ".parquet")))
}

tar_read(cell_growth_end) %>% 
  mutate(end = as.integer(end),
         start = as.integer(start),
         cell_no = as.integer(cell_no)) %>% 
  write_parquet(file.path(out_path, "cell_growth_end.parquet"))

tar_read(cell_growth_supp_end) %>% 
  mutate(end = as.integer(end),
         start = as.integer(start),
         cell_no = as.integer(cell_no)) %>% 
  write_parquet(file.path(out_path, "cell_growth_supp_end.parquet"))

gc()

## Further process ------------------------------------------------------------------------------------------------
### Regular cell growth -------------------------------------------------------------------------------------------
cell_growth_lims_fnms <- file.path(out_path, "cell_growth") %>% 
  list.files(full.names = T) %>% 
  str_subset("cell_growth_abbrev")

dom_lims <- map_dfr(cell_growth_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet(col_select = c("start", "species", "cell_no", "lim")) %>% 
    mutate(start = (lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           start = case_when(start == 13 ~ 1, start == 14 ~ 2, T ~ start),
           start = as.integer(start)) %>% 
    group_by(cell_no, species, start) %>%
    count(lim) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup()
})
dom_lims %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_dom_lims.parquet")))

growth_rate <- map_dfr(cell_growth_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet(c("start", "species", "cell_no", "growth_rate")) %>% 
    mutate(start = (lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           start = case_when(start == 13 ~ 1, start == 14 ~ 2, T ~ start),
           start = as.integer(start)) %>% 
    group_by(cell_no, species, start) %>%
    reframe(growth_rate = meanna(growth_rate))
})
growth_rate %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_rate.parquet")))

### Supplemented growth -------------------------------------------------------------------------------------------
cell_growth_lims_fnms <- file.path(out_path, "cell_growth") %>% 
  list.files(full.names = T) %>% 
  str_subset("cell_growth_supp_abbrev")

dom_lims <- map_dfr(cell_growth_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet(col_select = c("start", "species", "cell_no", "lim")) %>% 
    mutate(start = (lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           start = case_when(start == 13 ~ 1, start == 14 ~ 2, T ~ start),
           start = as.integer(start)) %>% 
    group_by(cell_no, species, start) %>%
    count(lim) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup()
})
dom_lims %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_supp_dom_lims.parquet")))

growth_rate <- map_dfr(cell_growth_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet(c("start", "species", "cell_no", "growth_rate")) %>% 
    mutate(start = (lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           start = case_when(start == 13 ~ 1, start == 14 ~ 2, T ~ start),
           start = as.integer(start)) %>% 
    group_by(cell_no, species, start) %>%
    reframe(growth_rate = meanna(growth_rate))
})
growth_rate %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_supp_rate.parquet")))

# Scenarios -------------------------------------------------------------------------------------------------------
tar_read(scens) %>% 
  qsave(file.path(out_path, "theo_scens.qs"))

# Everything contained in here is also in lims
# tar_read(scen_growth) %>% 
#   qsave(file.path(out_path, "scen_growth_all.qs"))

rbind(
  tar_read(scen_growth_armata),
  tar_read(scen_growth_taxiformis)
  ) %>% 
  mutate(t = as.integer(t),
         start = as.integer(start)) %>% 
  write_parquet(file.path(out_path, paste0("scen_growth.parquet")))

tar_read(scen_growth_end) %>% 
  write_parquet(file.path(out_path, "scen_growth_end_all.parquet"))


