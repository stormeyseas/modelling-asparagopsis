suppressMessages(suppressWarnings({
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
  library(qs)
  library(qs2)
}))
conflicts_prefer(dplyr::filter(), dplyr::select())

clean_memory <- function() {
  gc()
  terra::terraOptions(memfrac=0.5)
}

# Land area fraction ----------------------------------------------------------------------------------------------
if (!file.exists(file.path(raw_data_path, "BARRA_C2_cell_mask.Rdata")) | overwrite == T) {
  land_path <- file.path(raw_data_path, "sftlf_AUST-04_ERA5_historical_hres_BOM_BARRA-C2_v1.nc")
  
  nc <- nc_open(land_path)
  lat <- ncvar_get(nc, "lat")
  lat_ind <- which(lat > -45 & lat < -9)
  lon <- ncvar_get(nc, "lon")
  lon_ind <- which(lon > 108 & lon < 155)
  land <- sea <- ncvar_get(nc, "sftlf")
  land[land == 0] <- NA
  sea[sea == 1] <- NA
  
  # Create raster of land % values
  rast <- raster(list(x = lon[lon_ind], y = lat[lat_ind], z = land[lon_ind, lat_ind]))
  crs(rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  # Cut out some extraneous bits before doing the distance calculation to reduce computational load
  cutcells <- cellsFromExtent(rast, extent(148, 155, -18, 0))
  rast[cutcells] <- NA
  cutcells <- cellsFromExtent(rast, extent(115, 125.75, -13.2, 0))
  rast[cutcells] <- NA
  cutcells <- cellsFromExtent(rast, extent(146, 150, -12, 0))
  rast[cutcells] <- NA
  cutcells <- cellsFromExtent(rast, extent(125, 129, -11, 0))
  rast[cutcells] <- NA
  cutcells <- cellsFromExtent(rast, extent(140, 141.4, -10.2, 0))
  rast[cutcells] <- NA
  cutcells <- cellsFromExtent(rast, extent(143.25, 145, -10.6, 0))
  rast[cutcells] <- NA
  cutcells <- cellsFromExtent(rast, extent(141, 144, -9.6, 0))
  rast[cutcells] <- NA
  plot(rast)
  
  # Create raster of distance from land values, mask out cells too far from the coast
  # dist <- distance(rast) # takes ages, load if you can
  # save(dist, file = file.path(raw_data_path, "BARRA_C2_fulldist.Rdata"))
  load(file.path(raw_data_path, "BARRA_C2_fulldist.Rdata"))
  dist[dist > 71*10^3] <- NA
  plot(dist)
  
  cell_mask <- raster(list(x = lon[lon_ind], y = lat[lat_ind], z = sea[lon_ind, lat_ind]))
  cell_mask <- mask(x = cell_mask, mask = dist)
  plot(cell_mask)
  # plot(crop(cell_mask, extent(139.5, 145, -12, 0))) # Zoom in for detail
  
  # cell_coords <- rasterToPoints(cell_mask) %>% as.data.frame()
  # as.data.frame(cell_mask, xy = TRUE) %>% 
  #   ggplot(aes(x = x, y = y, fill = layer)) +
  #   geom_raster()
  
  save(lat, file = file.path(raw_data_path, "BARRA_C2_lats.Rdata"))
  save(lon, file = file.path(raw_data_path, "BARRA_C2_lons.Rdata"))
  save(lat_ind, file = file.path(raw_data_path, "BARRA_C2_lat_inds.Rdata"))
  save(lon_ind, file = file.path(raw_data_path, "BARRA_C2_lon_inds.Rdata"))
  save(cell_mask, file = file.path(raw_data_path, "BARRA_C2_cell_mask.Rdata"))
  
  nc_close(nc)
  rm(nc)
} else {
  load(file.path(raw_data_path, "BARRA_C2_lats.Rdata")) # lat
  load(file.path(raw_data_path, "BARRA_C2_lons.Rdata")) # lon
  load(file.path(raw_data_path, "BARRA_C2_lat_inds.Rdata")) # lat_ind
  load(file.path(raw_data_path, "BARRA_C2_lon_inds.Rdata")) # lon_ind
  load(file.path(raw_data_path, "BARRA_C2_cell_mask.Rdata")) # cell_mask
}

# Extract raw data ------------------------------------------------------------------------------------------------
# This script takes the raw (nc) files and turns them into raster files.
## Shortwave radiation --------------------------------------------------------------------------------------------
var <- "rsdsdir"
var_path <- file.path(raw_data_path, var)
var_nc_files <- list.files(var_path, full.names = T, pattern = ".nc")
var_rasters <- file.path(var_path, "rasters")
doy_raster_means <- file.path(var_path, "raster_means")
dir.create(var_rasters)
dir.create(doy_raster_means)

for (i in 1:length(var_nc_files)) {
  # Open data file
  nc <- nc_open(var_nc_files[i])
  time_unit <- nc$dim$time$units %>% 
    str_remove("days since ") %>% 
    as.Date()
  times <- ncvar_get(nc, 'time')
  times <- time_unit + lubridate::duration(times, "days")
  # Create filenames for each timestep
  jfiles <- file.path(var_rasters, str_c(var, "_", str_replace_all(str_replace_all(times, ":", "-"), " ", "-"), ".tif"))
  
  if (any(!file.exists(jfiles)) | overwrite == T) {
    # Pre-allocate array for efficiency
    var_data <- array(NA, dim=c(length(lon_ind), length(lat_ind), length(times)))
    
    # Read data in chunks for memory efficiency
    chunk_size <- 5
    n_chunks <- ceiling(length(times)/chunk_size)
    
    for(chunk in 1:n_chunks) {
      start_idx <- (chunk-1)*chunk_size + 1
      end_idx <- min(chunk*chunk_size, length(times))
      var_data[,,start_idx:end_idx] <- ncvar_get(nc, var, 
                                                  start=c(min(lon_ind), min(lat_ind), start_idx),
                                                  count=c(length(lon_ind), length(lat_ind), end_idx-start_idx+1))
    }
    nc_close(nc)
    
    times_df <- data.frame(time = times, index = 1:length(times))
    for (j in seq_along(times_df$index)) {
      fname <- file.path(var_rasters, str_c(var, "_", str_replace_all(str_replace_all(times[j], ":", "-"), " ", "-"), ".tif"))
      if (!file.exists(fname) | overwrite == T) {
        var_rast <- rast(list(x = lon[lon_ind], y = lat[lat_ind], z = var_data[,,j]))
        crs(var_rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        writeRaster(var_rast, fname, overwrite = T)
      }
    }
    clean_memory()
  } else {
    nc_close(nc)
  }
print(str_c(i, " of ", length(var_nc_files), " raw files done at ", Sys.time()))
}

## Surface temperature ---------------------------------------------------------------------------------------------
var <- "ts"
var_path <- file.path(raw_data_path, var)
var_nc_files <- list.files(var_path, full.names = T, pattern = ".nc")
var_rasters <- file.path(var_path, "rasters")
dir.create(var_rasters)

for (i in 1:length(var_nc_files)) {
  nc <- nc_open(var_nc_files[i])
  time_unit <- nc$dim$time$units %>% 
    str_remove("days since ") %>% 
    as.Date()
  times <- ncvar_get(nc, 'time')
  times <- time_unit + lubridate::duration(times, "days")
  jfiles <- file.path(var_rasters, str_c(var, "_", str_replace_all(str_replace_all(times, ":", "-"), " ", "-"), ".tif"))
  
  if (any(!file.exists(jfiles)) | overwrite == T) {
    # Pre-allocate array for efficiency
    var_data <- array(NA, dim=c(length(lon_ind), length(lat_ind), length(times)))
    
    # Read data in chunks for memory efficiency
    chunk_size <- 5
    n_chunks <- ceiling(length(times)/chunk_size)
    
    for(chunk in 1:n_chunks) {
      start_idx <- (chunk-1)*chunk_size + 1
      end_idx <- min(chunk*chunk_size, length(times))
      var_data[,,start_idx:end_idx] <- ncvar_get(nc, var, 
                                                 start=c(min(lon_ind), min(lat_ind), start_idx),
                                                 count=c(length(lon_ind), length(lat_ind), end_idx-start_idx+1))
    }
    nc_close(nc)
    
    times_df <- data.frame(time = times, index = 1:length(times))
    for (j in seq_along(times_df$index)) {
      fname <- file.path(var_rasters, str_c(var, "_", str_replace_all(str_replace_all(times[j], ":", "-"), " ", "-"), ".tif"))
      if (!file.exists(fname) | overwrite == T) {
        var_rast <- rast(list(x = lon[lon_ind], y = lat[lat_ind], z = var_data[,,j]))
        crs(var_rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        writeRaster(var_rast, fname, overwrite = T)
      }
    }
    clean_memory()
  } else {
    nc_close(nc)
  }
  print(str_c(i, " of ", length(var_nc_files), " raw files done at ", Sys.time()))
}

# Mean over DOYs --------------------------------------------------------------------------------------------------
# This takes the raw rasters from all years and creates a mean year for each variable
vars <- c("rsdsdir", "ts")

for (var in vars) {
  var_path <- file.path(raw_data_path, var)
  var_nc_files <- list.files(var_path, full.names = T) %>% 
    str_subset("raster", negate = T) %>% 
    str_subset("cell_vals", negate = T)
  var_rasters <- file.path(var_path, "rasters")
  doy_raster_means <- file.path(var_path, "raster_means")
  dir.create(doy_raster_means)
  
  rast_files <- data.frame(fnm = list.files(file.path(var_rasters), full.names = T),
                           time = list.files(file.path(var_rasters))) %>% 
    mutate(time = str_remove(time, str_c(var, "_")),
           time = str_remove(time, ".tif"),
           time = lubridate::parse_date_time(time, orders = "%Y-%m-%d-%H-%M-%S"),
           date = format(time, "%Y-%m-%d"),
           leap_year = case_when(leap_year(time) ~ T, T ~ F),
           leap_day = case_when(month(time) == 2 & day(time) == 29 ~ yday(time)-1,
                                leap_year == T & yday(time) > yday(as.Date("2024-02-28")) ~ yday(time)-1,
                                T ~ yday(time)))
  
  for (d in unique(rast_files$leap_day)) {
    fname <- file.path(doy_raster_means, str_c(var, "_doy_", fixnum(d), ".tif"))
    if (!file.exists(fname) | overwrite == T) {
      rf <- rast_files[rast_files$leap_day == d, ]
      r_mean <- terra::app(terra::rast(rf$fnm), mean)
      writeRaster(x = r_mean, filename = fname, overwrite = T)
      clean_memory()
    }
    print(str_c(var, " file DOY ", fixnum(d), " saved"))
  }
}

# Final rasters ---------------------------------------------------------------------------------------------------
land_rast <- terra::rast(cell_mask)
writeRaster(x = land_rast, filename = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_land_rast.tif"), overwrite = T)

doy_raster_means <- file.path(raw_data_path, "rsdsdir", "raster_means") %>% 
  list.files(full.names = T)
BARRA_C2_rsdsdir_data <- terra::rast(doy_raster_means)
writeRaster(x = BARRA_C2_rsdsdir_data, filename = file.path(final_data_path, "BARRA_C2_rsdsdir_data.tif"), overwrite = T)
rm(BARRA_C2_rsdsdir_data)

doy_raster_means <- file.path(raw_data_path, "ts", "raster_means") %>% 
  list.files(full.names = T)
BARRA_C2_ts_data <- terra::rast(doy_raster_means)
writeRaster(x = BARRA_C2_ts_data, filename = file.path(final_data_path, "BARRA_C2_ts_data.tif"), overwrite = T)
rm(BARRA_C2_ts_data)



