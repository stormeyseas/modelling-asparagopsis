library(ncdf4)
library(stringr)
library(lubridate)
library(reshape2)
library(magrittr)
library(dplyr)
library(arrow)
library(oceanmap)
library(ggplot2)
library(conflicted)
library(rasterVis)
conflicts_prefer(dplyr::filter(), dplyr::select())

fixnum <- function(n, digits = 4) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}
fixnum2 <- function(n, digits = 6) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}

# Land area fraction ----------------------------------------------------------------------------------------------
if (!file.exists(file.path(data_path, "BARRA2_cell_mask.Rdata")) | overwrite == T) {
  land_path <- file.path(data_path, "sftlf_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1.nc")

  nc <- nc_open(land_path)
  lat <- ncvar_get(nc, "lat")
  lat_ind <- which(lat > -45 & lat < -9)
  lon <- ncvar_get(nc, "lon")
  lon_ind <- which(lon > 108 & lon < 155)
  land <- sea <- ncvar_get(nc, "sftlf")
  land[land == 0] <- NA
  sea[sea == 100] <- NA
  
  # Create raster of land % values
  rast <- raster(list(x = lon[lon_ind], y = lat[lat_ind], z = land[lon_ind, lat_ind]))
  crs(rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  plot(rast)
  
  # Create raster of distance from land values, mask out cells too far from the coast
  dist <- distance(rast) # takes ages
  dist[dist$layer > 70*10^3] <- NA
  # plot(dist)
  cell_mask <- raster(list(x = lon[lon_ind], y = lat[lat_ind], z = sea[lon_ind, lat_ind]))
  crs(cell_mask) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  cell_mask <- mask(x = cell_mask, mask = dist)
  
  # Cut out some extraneous bits
  cutcells <- cellsFromExtent(cell_mask, extent(148, 155, -18, 0))
  cell_mask$layer[cutcells] <- NA
  cutcells <- cellsFromExtent(cell_mask, extent(115, 125.75, -13.2, 0))
  cell_mask$layer[cutcells] <- NA
  cutcells <- cellsFromExtent(cell_mask, extent(146, 150, -12, 0))
  cell_mask$layer[cutcells] <- NA
  cutcells <- cellsFromExtent(cell_mask, extent(125, 129, -11, 0))
  cell_mask$layer[cutcells] <- NA
  cutcells <- cellsFromExtent(cell_mask, extent(140, 141.4, -10.2, 0))
  cell_mask$layer[cutcells] <- NA
  cutcells <- cellsFromExtent(cell_mask, extent(143.25, 145, -10.6, 0))
  cell_mask$layer[cutcells] <- NA
  cutcells <- cellsFromExtent(cell_mask, extent(141, 144, -9.6, 0))
  cell_mask$layer[cutcells] <- NA
  
  plot(cell_mask)
  plot(crop(cell_mask, extent(139.5, 145, -12, 0)))
  
  as.data.frame(cell_mask, xy = TRUE) %>% 
    ggplot(aes(x = x, y = y, fill = layer)) +
    geom_raster()
  
  save(lat, file = file.path(data_path, "BARRA2_lats.Rdata"))
  save(lon, file = file.path(data_path, "BARRA2_lons.Rdata"))
  save(lat_ind, file = file.path(data_path, "BARRA2_lat_inds.Rdata"))
  save(lon_ind, file = file.path(data_path, "BARRA2_lon_inds.Rdata"))
  save(cell_mask, file = file.path(data_path, "BARRA2_cell_mask.Rdata"))
  
  nc_close(nc)
  rm(nc)
} else {
  load(file.path(data_path, "BARRA2_lats.Rdata"))
  load(file.path(data_path, "BARRA2_lons.Rdata"))
  load(file.path(data_path, "BARRA2_lat_inds.Rdata"))
  load(file.path(data_path, "BARRA2_lon_inds.Rdata"))
  load(file.path(data_path, "BARRA2_cell_mask.Rdata"))
}

## Non-NA cells ---------------------------------------------------------------------------------------------------
if (!file.exists(file.path(data_path, "BARRA2_cell_coords.Rdata")) | overwrite == T) {
  cell_coords <- cell_mask %>% as.data.frame()
  all_cells <- cell_mask
  all_cells[is.na(all_cells$layer)] <- 10
  cell_coords <- rasterToPoints(all_cells) %>% as.data.frame()
  cell_coords$cell_no <- 1:nrow(cell_coords)
  cell_coords <- cell_coords %>% filter(layer != 10)
  save(cell_coords, file = file.path(data_path, "BARRA2_cell_coords.Rdata"))
} else {
  load(file.path(data_path, "BARRA2_cell_coords.Rdata"))
}

# Extract data ----------------------------------------------------------------------------------------------------
# Shortwave radiation ---------------------------------------------------------------------------------------------
var <- "rsdsdir"
var_path <- file.path(data_path, var)
data_file <- list.files(var_path, full.names = T) %>% str_subset("raster", negate = T) %>% str_subset("cell_vals", negate = T)
dest_path <- file.path(var_path, "rasters")
dest_path_2 <- file.path(var_path, "raster_means")
dest_path_3 <- file.path("data_raw", "BARRA2", "cell_vals")

# Raster every timestep from the nc files, mask by cell_mask
for (i in 1:length(data_file)) {
  nc <- nc_open(data_file[i])
  time_unit <- nc$dim$time$units %>% 
    str_remove("days since ") %>% 
    as.Date()
  times <- ncvar_get(nc, 'time')
  times <- time_unit + lubridate::duration(times, "days")
  
  # Get data, constrain a bit
  var_data <- ncvar_get(nc, var, start = c(min(lon_ind), min(lat_ind), 1), count = c(length(lon_ind), length(lat_ind), length(times)))
  
  for (j in 1:length(times)) {
    var_fnm <- file.path(dest_path, str_c(var, "_", str_replace_all(str_replace_all(times[j], ":", "-"), " ", "-"), ".tif"))
    if (!file.exists(var_fnm) | overwrite == T) {
      var_rast <- raster(list(x = lon[lon_ind], y = lat[lat_ind], z = var_data[,,j]))
      crs(var_rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      var_rast <- crop(var_rast, extent(cell_mask))
      var_rast <- mask(var_rast, cell_mask)
      writeRaster(var_rast, , overwrite = T)
    }
  }
  nc_close(nc)
}
rm(nc)

# Mean over all doys
rast_files <- data.frame(
  fnm = list.files(file.path(dest_path), full.names = T),
  time = list.files(file.path(dest_path))
) %>% 
  mutate(time = str_remove(time, str_c(var, "_")),
         time = str_remove(time, ".tif"),
         time = lubridate::parse_date_time(time, orders = "%Y-%m-%d-%H-%M-%S"),
         date = format(time, "%Y-%m-%d"),
         leap_year = case_when(leap_year(time) ~ T, T ~ F),
         leap_day = case_when(month(time) == 2 & day(time) == 29 ~ yday(time)-1,
                              leap_year == T & yday(time) > yday(as.Date("2024-02-28")) ~ yday(time)-1,
                              T ~ yday(time)))
for (d in unique(rast_files$leap_day)) {
  d_fnm <- file.path(dest_path_2, str_c(var, "_doy_", fixnum(d), ".tif"))
  if (!file.exists(d_fnm) | overwrite == T) {
    rf <- rast_files[rast_files$leap_day == d, ]
    rast_stack <- terra::rast(rf$fnm)
    r_mean <- terra::app(rast_stack, mean)
    writeRaster(r_mean, d_fnm, overwrite = T)
  }
  print(str_c(var, " doy ", fixnum(d), " saved"))
}

# Extract each cell from mean year raster stack
fnms <- list.files(file.path(dest_path_2), full.names = T)
rast_stack <- terra::rast(fnms)
for (c in 1:length(cell_coords$cell_no)) {
  cname <- file.path(dest_path_3, str_c(var, "_cell_", fixnum(cell_coords$cell_no[c], 6), ".qs"))
  if (!file.exists(cname)) {
    cv <- raster::extract(rast_stack, cell_coords$cell_no[c]) %>% unlist() %>% unname()
    qsave(cv, cname)
  }
  if (c %in% as.integer(seq(0, length(cell_coords$cell_no), length.out = 101))) {
    print(str_c(c, " of ", length(cell_coords$cell_no), " saved"))
  }
}

# Surface temperature ---------------------------------------------------------------------------------------------
var <- "ts"
var_path <- file.path(data_path, var)
data_file <- list.files(var_path, full.names = T) %>% str_subset("raster", negate = T) %>% str_subset("cell_vals", negate = T)
dest_path <- file.path(var_path, "rasters")
dest_path_2 <- file.path(var_path, "raster_means")
dest_path_3 <- file.path("data_raw", "BARRA2", "cell_vals")

# Raster every timestep from the nc files, mask by cell_mask
for (i in 1:length(data_file)) {
  nc <- nc_open(data_file[i])
  time_unit <- nc$dim$time$units %>% 
    str_remove("days since ") %>% 
    as.Date()
  times <- ncvar_get(nc, 'time')
  times <- time_unit + lubridate::duration(times, "days")
  
  # Get data, constrain a bit
  var_data <- ncvar_get(nc, var, start = c(min(lon_ind), min(lat_ind), 1), count = c(length(lon_ind), length(lat_ind), length(times)))
  
  for (j in 1:length(times)) {
    var_fnm <- file.path(dest_path, str_c(var, "_", str_replace_all(str_replace_all(times[j], ":", "-"), " ", "-"), ".tif"))
    if (!file.exists(var_fnm) | overwrite == T) {
      var_rast <- raster(list(x = lon[lon_ind], y = lat[lat_ind], z = var_data[,,j]))
      crs(var_rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      var_rast <- crop(var_rast, extent(cell_mask))
      var_rast <- mask(var_rast, cell_mask)
      writeRaster(var_rast, , overwrite = T)
    }
  }
  nc_close(nc)
}
rm(nc)

# Mean over all doys
rast_files <- data.frame(
  fnm = list.files(file.path(dest_path), full.names = T),
  time = list.files(file.path(dest_path))
) %>% 
  mutate(time = str_remove(time, str_c(var, "_")),
         time = str_remove(time, ".tif"),
         time = lubridate::parse_date_time(time, orders = "%Y-%m-%d-%H-%M-%S"),
         date = format(time, "%Y-%m-%d"),
         leap_year = case_when(leap_year(time) ~ T, T ~ F),
         leap_day = case_when(month(time) == 2 & day(time) == 29 ~ yday(time)-1,
                              leap_year == T & yday(time) > yday(as.Date("2024-02-28")) ~ yday(time)-1,
                              T ~ yday(time)))
for (d in unique(rast_files$leap_day)) {
  d_fnm <- file.path(dest_path_2, str_c(var, "_doy_", fixnum(d), ".tif"))
  if (!file.exists(d_fnm) | overwrite == T) {
    rf <- rast_files[rast_files$leap_day == d, ]
    rast_stack <- terra::rast(rf$fnm)
    r_mean <- terra::app(rast_stack, mean)
    writeRaster(r_mean, d_fnm, overwrite = T)
  }
  print(str_c(var, " doy ", fixnum(d), " saved"))
}

# Extract each cell from mean year raster stack
fnms <- list.files(file.path(dest_path_2), full.names = T)
rast_stack <- terra::rast(fnms)
for (c in 1:length(cell_coords$cell_no)) {
  cname <- file.path(dest_path_3, str_c(var, "_cell_", fixnum(cell_coords$cell_no[c], 6), ".qs"))
  if (!file.exists(cname)) {
    cv <- raster::extract(rast_stack, cell_coords$cell_no[c]) %>% unlist() %>% unname()
    qsave(cv, cname)
  }
  if (c %in% as.integer(seq(0, length(cell_coords$cell_no), length.out = 101))) {
    print(str_c(c, " of ", length(cell_coords$cell_no), " saved, ", round(100*c/length(cell_coords$cell_no), 2), "% done"))
  }
}

