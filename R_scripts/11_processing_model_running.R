# This script takes files from the model_running targets pipeline and saves them as files, which are easier to load

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
library(terra)
library(units)
library(tictoc)
library(streamMetabolizer)
library(arrow)
library(purrr)
library(here)
library(future)
library(future.apply)
library(furrr)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)
source("R_scripts/00_functions.R")

Sys.setenv(TAR_PROJECT = "project_running")
out_path <- here() %>% file.path("data", "processed_model_running")
mani <- tar_manifest() %>%
  select(-command)

month_ydays <- lapply(1:12, function(m) {
  first_day <- ymd(paste(2023, m, 1, sep = "-"))
  days_seq <- seq(first_day, first_day + days(days_in_month(first_day) - 1), by = "day")
  yday(days_seq)
}) %>% setNames(1:12)

# General cell info -----------------------------------------------------------------------------------------------
BARRA_R2_cell_nos <- tar_read(BARRA_R2_cell_nos) 
BARRA_R2_cell_nos %>% 
  qsave(file.path(out_path, "cells_included.qs"))

tar_read(cells_to_omit) %>% 
  qsave(file.path(out_path, "cells_omitted.qs"))

BARRA_R2_cell_coords <- tar_read(BARRA_R2_cell_coords) %>% 
  filter(cell_no %in% BARRA_R2_cell_nos)
BARRA_R2_cell_coords %>% 
  qsave(file.path(out_path, "BARRA_R2_cell_coords.parquet"))
  
# Cell input timeseries processing
BARRA_R2_cell_nos_chunked <- tar_read(BARRA_R2_cell_nos_chunked)
# for (i in seq_along(BARRA_R2_cell_nos_chunked)) {
#   tar_read(cell_inputs_chunked, branches = i)[[1]] %>% 
#     unname()
# } %>% dplyr::bind_rows()
# write_parquet(df, file.path(out_path, "cell_inputs", paste0("cell_inputs_chunked_", fixnum(i), ".parquet")))
gc()

# Cell growth -----------------------------------------------------------------------------------------------------
cell_branches <- length(tar_read(start_ydays_growth)) * length(tar_read(BARRA_R2_cell_nos_chunked))
growth_chunks <- seq(0, cell_branches, length.out = 101) %>% as.integer()

for (i in 1:(length(growth_chunks)-1)) {
  stt <- growth_chunks[i]+1
  end <- growth_chunks[i+1]
  
  # Normal growth
  normal_growth <- rbind(
    tar_read(cell_growth_armata, branches = stt:end),
    tar_read(cell_growth_taxiformis, branches = stt:end)
  ) %>% 
    mutate(t = as.integer(t),
           start = as.integer(start),
           cell_no = as.integer(cell_no))
  write_parquet(normal_growth, file.path(out_path, "cell_growth", paste0("cell_growth_", fixnum(i), ".parquet")))
  
  normal_growth <- normal_growth %>% 
    dplyr::select(t, growth_rate, start, species, cell_no, lim)
  write_parquet(normal_growth, file.path(out_path, "cell_growth", paste0("cell_growth_abbrev_", fixnum(i), ".parquet")))
  
  # Supplemented growth
  supplemented_growth <- rbind(
    tar_read(cell_growth_supp_armata, branches = stt:end),
    tar_read(cell_growth_supp_taxiformis, branches = stt:end)
  ) %>% 
    mutate(t = as.integer(t),
           start = as.integer(start),
           cell_no = as.integer(cell_no))
  write_parquet(supplemented_growth, file.path(out_path, "cell_growth", paste0("cell_growth_supp_", fixnum(i), ".parquet"))) 
  
  supplemented_growth <- supplemented_growth %>% 
    dplyr::select(t, growth_rate, start, species, cell_no, lim)
  write_parquet(supplemented_growth, file.path(out_path, "cell_growth", paste0("cell_growth_supp_abbrev_", fixnum(i), ".parquet")))
}

tar_load(cell_growth_end)
write_parquet(cell_growth_end, file.path(out_path, "cell_growth_end.parquet"))

tar_load(cell_growth_supp_end)
write_parquet(cell_growth_supp_end, file.path(out_path, "cell_growth_supp_end.parquet"))

gc()
rm(cell_growth_end, cell_growth_supp_end)

## Further process ------------------------------------------------------------------------------------------------
### Regular cell growth -------------------------------------------------------------------------------------------
cell_growth_lims_fnms <- file.path(out_path, "cell_growth") %>% 
  list.files(full.names = T) %>% 
  str_subset("cell_growth_abbrev")

# This gives a daily growth limiter for every cell
armata_growth <- cell_growth_lims_fnms %>%
  purrr::map(function(fnm) {
    arrow::read_parquet(fnm) %>%
      filter(species == "armata") %>%
      mutate(t = as.integer(case_when(t > 365 ~ t-365, T ~ t))) %>%
      distinct(cell_no, t, growth_rate) %>% 
      mutate(month = as.factor(month(parse_date_time(t, orders = "j"))))
  }) %>%
  dplyr::bind_rows() %>%
  merge(select(BARRA_R2_cell_coords, -layer), by = "cell_no") %>%
  write_parquet(file.path(out_path, "armata_growth.parquet"))

taxiformis_growth <- cell_growth_lims_fnms %>%
  purrr::map(function(fnm) {
    arrow::read_parquet(fnm) %>%
      filter(species == "taxiformis") %>%
      mutate(t = as.integer(case_when(t > 365 ~ t-365, T ~ t))) %>%
      distinct(cell_no, t, growth_rate) %>% 
      mutate(month = as.factor(month(parse_date_time(t, orders = "j"))))
  }) %>%
  dplyr::bind_rows() %>%
  merge(select(BARRA_R2_cell_coords, -layer), by = "cell_no") %>%
  write_parquet(file.path(out_path, "taxiformis_growth.parquet"))

# Save growth as rasters (for maps)
for (m in 1:12) {
  r <- armata_growth %>% 
    filter(month == m) %>% 
    group_by(latitude, longitude, month) %>%
    reframe(growth_rate = meanna(growth_rate))
  
  rr <- rast(r[ , c("longitude", "latitude", "growth_rate")])
  crs(rr) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(rr, file.path(out_path, str_c("armata_growth_rate_raster_", fixnum(m,2), ".tif")), overwrite = T)
  
  r <- taxiformis_growth %>% 
    filter(month == m) %>% 
    group_by(latitude, longitude, month) %>%
    reframe(growth_rate = meanna(growth_rate))

  r <- rast(r[ , c("longitude", "latitude", "growth_rate")])
  writeRaster(rr, file.path(out_path, str_c("taxiformis_growth_rate_raster_", fixnum(m,2), ".tif")), overwrite = T)
}

### Regular growth lims -------------------------------------------------------------------------------------------
cell_growth_lims_fnms <- file.path(out_path, "cell_growth") %>% 
  list.files(full.names = T) %>% 
  str_subset("cell_growth_abbrev")

# This gives a daily growth limiter for every cell
armata_lims <- cell_growth_lims_fnms %>%
  purrr::map(function(fnm) {
    arrow::read_parquet(fnm) %>%
      filter(species == "armata") %>%
      mutate(t = as.integer(case_when(t > 365 ~ t-365, T ~ t))) %>%
      distinct(cell_no, t, lim) %>% 
      mutate(month = as.factor(month(parse_date_time(t, orders = "j"))))
  }) %>%
  dplyr::bind_rows() %>%
  merge(select(BARRA_R2_cell_coords, -layer), by = "cell_no") %>%
  write_parquet(file.path(out_path, "armata_lims.parquet"))

armata_lims %>%
  dplyr::select(c(latitude, month, cell_no, lim)) %>%
  write_parquet(file.path(out_path, "forhist_armata.parquet"))

# This gives a daily growth limiter for every cell
taxiformis_lims <- cell_growth_lims_fnms %>%
  purrr::map(function(fnm) {
    arrow::read_parquet(fnm) %>%
      filter(species == "taxiformis") %>%
      mutate(t = as.integer(case_when(t > 365 ~ t-365, T ~ t))) %>%
      distinct(cell_no, t, lim) %>% 
      mutate(month = as.factor(month(parse_date_time(t, orders = "j"))))
  }) %>%
  dplyr::bind_rows() %>%
  merge(select(BARRA_R2_cell_coords, -layer), by = "cell_no") %>%
  write_parquet(file.path(out_path, "taxiformis_lims.parquet"))

taxiformis_lims %>% 
  dplyr::select(c(latitude, month, cell_no, lim)) %>%
  write_parquet(file.path(out_path, "forhist_taxiformis.parquet"))

# Calculate dominant limiters and save as rasters (for maps)
for (m in 1:12) {
  r <- armata_lims %>% 
    filter(month == m) %>% 
    group_by(latitude, longitude) %>%
    count(lim) %>%
    slice_max(n, with_ties = FALSE) %>% 
    ungroup() %>% 
    # Manually set levels to make double sure they're consistent
    mutate(lim = factor(lim, levels = c("T_lim", "I_lim", "Q_lim", "S_lim")), 
           lim = as.numeric(lim))
  
  rr <- rast(r[ , c("longitude", "latitude", "lim")])
  crs(rr) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(rr, file.path(out_path, str_c("armata_dom_lims_raster_", fixnum(m,2), ".tif")), overwrite = T)
  
  r <- taxiformis_lims %>% 
    filter(month == m) %>% 
    group_by(latitude, longitude) %>%
    count(lim) %>%
    slice_max(n, with_ties = FALSE) %>% 
    ungroup() %>% 
    # Manually set levels to make double sure they're consistent
    mutate(lim = factor(lim, levels = c("T_lim", "I_lim", "Q_lim", "S_lim")), 
           lim = as.numeric(lim))
  
  rr <- rast(r[ , c("longitude", "latitude", "lim")])
  crs(rr) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(rr, file.path(out_path, str_c("taxiformis_dom_lims_raster_", fixnum(m,2), ".tif")), overwrite = T)
}

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


