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
for (i in seq_along(BARRA_C2_cell_nos_chunked)) {
  tar_read(cell_inputs_chunked, branches = i)[[1]] %>% 
    unname() %>% 
    bind_rows() %>% 
    write_parquet(file.path(out_path, "cell_inputs", paste0("cell_inputs_chunked_", fixnum(i), ".parquet")))
}
gc()

# Cell growth -----------------------------------------------------------------------------------------------------
cell_branches <- 2 * length(tar_read(start_ydays_growth)) * length(tar_read(BARRA_C2_cell_nos_chunked))
growth_chunks <- c(seq(0, cell_branches, by = 100), last(cell_branches))

for (i in 1:(length(growth_chunks)-1)) {
  stt <- growth_chunks[i]+1
  end <- growth_chunks[i+1]
  
  tar_read(cell_growth_lims, branches = stt:end) %>% 
    mutate(species = as.factor(species)) %>% 
    write_parquet(file.path(out_path, "cell_growth", paste0("cell_growth_lims_", fixnum(i), ".parquet"))) %>% 
    dplyr::select(t, growth_rate, start, species, cell_no, lim) %>% 
    write_parquet(file.path(out_path, "cell_growth", paste0("cell_growth_lims_abbrev_", fixnum(i), ".parquet")))
}

tar_read(cell_growth_end) %>% 
  write_parquet(file.path(out_path, "cell_growth_end.parquet"))

gc()

## Further process ------------------------------------------------------------------------------------------------
# List all cell_growth_lims files
cell_growth_lims_fnms <- file.path(out_path, "cell_growth") %>% 
  list.files(full.names = T) %>% 
  str_subset("cell_growth_lims_abbrev")

dom_lims <- map_dfr(cell_growth_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet(col_select = c("start", "species", "cell_no", "lim")) %>% 
    mutate(start = as.integer(lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           start = case_when(start == 13 ~ 1, start == 14 ~ 2, T ~ start)) %>% 
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
    mutate(start = as.integer(lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           start = case_when(start == 13 ~ 1, start == 14 ~ 2, T ~ start)) %>% 
    group_by(cell_no, species, start) %>%
    reframe(growth_rate = meanna(growth_rate))
})
growth_rate %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_rate.parquet")))

# Cell growth extra -----------------------------------------------------------------------------------------------
for (i in 2:length(growth_chunks)) {
  #   Everything contained in here is also in lims
  #   tar_read(cell_growth_extra, branches = (growth_chunks[i-1]+1):growth_chunks[i]) %>%
  #     write_parquet(file.path(out_path, "cell_growth_extra", paste0("cell_growth_extra_", fixnum(i-1), ".parquet")))
  tar_read(cell_growth_extra_lims, branches = (growth_chunks[i-1]+1):growth_chunks[i]) %>% 
    mutate(species = as.factor(species)) %>% 
    write_parquet(file.path(out_path, "cell_growth_extra", paste0("cell_growth_extra_lims_", fixnum(i-1), ".parquet"))) %>% 
    dplyr::select(t, growth_rate, start, species, cell_no, lim) %>% 
    write_parquet(file.path(out_path, "cell_growth_extra", paste0("cell_growth_extra_lims_abbrev_", fixnum(i-1), ".parquet")))
}
tar_read(cell_growth_extra_end) %>%
  qsave(file.path(out_path, "cell_growth_extra_end.qs"))
gc()

## Further process ------------------------------------------------------------------------------------------------
cell_growth_extra_lims_fnms <- file.path(out_path, "cell_growth_extra") %>% 
  list.files(full.names = T) %>% 
  str_subset("abbrev")

purrr::map_dfr(cell_growth_extra_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet() %>% 
    mutate(lim = as.factor(lim),
           month = as.integer(lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           month = case_when(month == 13 ~ 1, month == 14 ~ 2, T ~ month)) %>% 
    group_by(cell_no, species, month) %>%
    count(lim) %>%
    slice_max(n, with_ties = FALSE) %>%
    select(cell_no, species, month, dom_lim = lim)
}) %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_extra_dom_lims.parquet")))

purrr::map_dfr(cell_growth_extra_lims_fnms, function(fnm) {
  fnm %>% 
    arrow::read_parquet() %>% 
    mutate(lim = as.factor(lim),
           month = as.integer(lubridate::month(parse_date_time(x = start, orders = "j"))+1),
           month = case_when(month == 13 ~ 1, month == 14 ~ 2, T ~ month)) %>% 
    group_by(cell_no, species, month) %>%
    reframe(growth_rate = meanna(growth_rate))
}) %>% 
  write_parquet(file.path(out_path, paste0("cell_growth_extra_rate.parquet")))

# Scenarios -------------------------------------------------------------------------------------------------------
tar_read(scens) %>% 
  qsave(file.path(out_path, "theo_scens.qs"))

# Everything contained in here is also in lims
# tar_read(scen_growth) %>% 
#   qsave(file.path(out_path, "scen_growth_all.qs"))

tar_read(scen_growth_end) %>% 
  qsave(file.path(out_path, "scen_growth_end_all.qs"))

tar_read(scen_growth_lims) %>% 
  qsave(file.path(out_path, "scen_growth_lims_all.qs"))

