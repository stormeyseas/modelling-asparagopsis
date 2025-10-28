# Load required packages
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
library(terra)
library(qs)
library(qs2)
library(future)
library(furrr)
conflicts_prefer(dplyr::filter(), dplyr::select())

# Setup parallel processing
# plan(multisession, workers = availableCores() - 1)

# Function to clear memory
clean_memory <- function() {
  gc()
  terra::terraOptions(memfrac=0.5)
}

# Land area fraction ----------------------------------------------------------------------------------------------
tmp <- file.path(raw_data_path, "ocean_salt") %>% list.files(full.names = T) %>% str_subset(".nc") %>% first() %>% nc_open()
# salt[xt_ocean,yt_ocean,st_ocean,Time]
dep <- ncvar_get(tmp, "st_ocean")
lon <- ncvar_get(tmp, "xt_ocean")
lat <- ncvar_get(tmp, "yt_ocean")

dep_ind <- which(dep <= 25)
lon_ind <- which(lon >=  108 & lon <= 155)
lat_ind <- which(lat >=  -44.99 & lat <= -8)

save(lat, file = file.path(raw_data_path, "BRAN2023_lats.Rdata"))
save(lon, file = file.path(raw_data_path, "BRAN2023_lons.Rdata"))
save(dep, file = file.path(raw_data_path, "BRAN2023_deps.Rdata"))
save(lat_ind, file = file.path(raw_data_path, "BRAN2023_lat_inds.Rdata"))
save(lon_ind, file = file.path(raw_data_path, "BRAN2023_lon_inds.Rdata"))
save(dep_ind, file = file.path(raw_data_path, "BRAN2023_dep_inds.Rdata"))

nc_close(tmp)
clean_memory()

# Get variables ---------------------------------------------------------------------------------------------------
# load(file.path(raw_data_path, "BRAN2023_cell_coords.Rdata")) # cell_coords
load(file = file.path(raw_data_path, "BRAN2023_lats.Rdata")) # lat
load(file = file.path(raw_data_path, "BRAN2023_lons.Rdata")) # lon
load(file = file.path(raw_data_path, "BRAN2023_deps.Rdata")) # dep
load(file = file.path(raw_data_path, "BRAN2023_lat_inds.Rdata")) # lat_ind
load(file = file.path(raw_data_path, "BRAN2023_lon_inds.Rdata")) # lon_ind
load(file = file.path(raw_data_path, "BRAN2023_dep_inds.Rdata")) # dep_ind

# Convert cell_mask from RasterLayer (package raster) to SpatRaster (package terra)
#load(file.path(raw_data_path, "BRAN2023_cell_mask.Rdata")) # cell_mask
#BARRA_cell_mask <- rast(cell_mask)
#rm(cell_mask)

# Function to extract raw data from BRAN files
extract_BRAN_data <- function(
  var,        # variable name for file naming purposes
  var1,       # variable name inside the data nc files
  var_path    # path to raw files
) {
  data_file <- list.files(var_path, full.names = T) %>% str_subset(".nc")
  dest_path <- file.path(var_path, "rasters")

  for (i in seq_along(data_file)) {
    nc <- nc_open(data_file[i])

    # Get times
    time_unit <- nc$dim$Time$units %>% 
      str_remove("days since ") %>% 
      as.Date()
    times <- ncvar_get(nc, 'Time')
    times <- time_unit + lubridate::duration(times, "days")

    # Pre-allocate array and define chunks for efficiency
    var_data <- array(NA, dim=c(length(lon_ind), length(lat_ind), length(dep_ind), length(times)))
    chunk_size <- 10
    n_chunks <- ceiling(length(times)/chunk_size)
    
    for(chunk in 1:n_chunks) {
      # Extract data (only specified indices of lat, lon, depth)
      start_idx <- (chunk-1)*chunk_size + 1
      end_idx <- min(chunk*chunk_size, length(times))
      var_data[,,,start_idx:end_idx] <- ncvar_get(
        nc, 
        var1, 
        start=c(min(lon_ind), min(lat_ind), min(dep_ind), start_idx),
        count=c(length(lon_ind), length(lat_ind), length(dep_ind), end_idx-start_idx+1)
      )
    }
    nc_close(nc)
    
    # Name new files by time
    jfiles <- file.path(dest_path, str_c(var, "_", str_replace_all(str_replace_all(times, ":", "-"), " ", "-"), ".tif"))

    for (j in seq_along(jfiles)) {
      # Extract one time, average over depths
      var_data_j <- apply(var_data[,,,j], c(1,2), mean)
      var_rast <- rast(list(x = lon[lon_ind], y = lat[lat_ind], z = var_data_j))
      crs(var_rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      writeRaster(var_rast, jfiles[j], overwrite = T)
    }

    clean_memory()
    print(str_c(i, " of ", length(data_file), " raw files done at ", Sys.time()))
  }
}

## Extract salinity data to rasters ------------------------------------------------------------------------------
extract_BRAN_data(
  var = "ocean_salt",
  var1 = "salt",
  var_path = file.path(raw_data_path, "ocean_salt")
)

## Raster & mask U --------------------------------------------------------------------------------------------------------
extract_BRAN_data(
  var = "ocean_u",
  var1 = "u",
  var_path = file.path(raw_data_path, "ocean_u")
)

## Raster & mask V ------------------------------------------------------------------------------------------------
extract_BRAN_data(
  var = "ocean_v",
  var1 = "v",
  var_path = file.path(raw_data_path, "ocean_v")
)

# Create raster means ---------------------------------------------------------------------------------------------
create_raster_mean <- function(
  var,
  var1,
  var_path
) {
  raster_path <- file.path(var_path, "rasters")
  dest_path <- file.path(var_path, "raster_means")

  # Create new filenames for each day (no leap years)
  rast_files <- data.frame(
    fnm = list.files(file.path(raster_path), full.names = T),
    time = list.files(file.path(raster_path))
  ) %>% 
    mutate(
      time = str_remove(time, str_c(var, "_")),
      time = str_remove(time, ".tif"),
      time = lubridate::parse_date_time(time, orders = "%Y-%m-%d-%H-%M-%S"),
      date = format(time, "%Y-%m-%d"),
      leap_year = case_when(leap_year(time) ~ T, T ~ F),
      leap_day = case_when(
        month(time) == 2 & day(time) == 29 ~ yday(time)-1,
        leap_year == T & yday(time) > yday(as.Date("2024-02-28")) ~ yday(time)-1,
        T ~ yday(time)
      )
    )

  # Create mean raster for all non-leap-year ydays
  for (d in unique(rast_files$leap_day)) {
    fname <- file.path(dest_path, str_c(var, "_doy_", fixnum(d), ".tif"))

    rf <- rast_files[rast_files$leap_day == d, ]
    r_mean <- terra::app(terra::rast(rf$fnm), mean)
    writeRaster(x = r_mean, filename = fname, overwrite = T)
    clean_memory()

    print(str_c(var, " file DOY ", fixnum(d), " saved"))
  }
}

## Mean salinity over DOYs ----------------------------------------------------------------------------------------
create_raster_mean(
  var = "ocean_salt",
  var1 = "salt",
  var_path = file.path(raw_data_path, "ocean_salt")
)

## Mean U over DOYs -----------------------------------------------------------------------------------------------
create_raster_mean(
  var = "ocean_u",
  var1 = "u",
  var_path = file.path(raw_data_path, "ocean_u")
)

## Mean V over DOYs ----------------------------------------------------------------------------------------------
create_raster_mean(
  var = "ocean_v",
  var1 = "v",
  var_path = file.path(raw_data_path, "ocean_v")
)

## Final raster for UV -------------------------------------------------------------------------------------------
fnames <- data.frame(
  u_files = file.path(raw_data_path, "ocean_u", "raster_means") %>% list.files(full.names = T),
  v_files = file.path(raw_data_path, "ocean_v", "raster_means") %>% list.files(full.names = T)
)
dest_path <- file.path(raw_data_path, "ocean_uv", "raster_means")
dir.create(dest_path, recursive = T)
dest_files <- list.files(file.path(raw_data_path, "ocean_u", "raster_means")) %>% 
  str_replace_all("ocean_u", "ocean_uv")

purrr::map2(1:nrow(fnames), dest_files, function(f, dest) {
  u <- terra::rast(fnames$u_files[f])
  v <- terra::rast(fnames$v_files[f])
  uv <- sqrt(u^2 + v^2)
  writeRaster(uv, file.path(dest_path, dest), overwrite = T)
})

# # Final rasters ---------------------------------------------------------------------------------------------------
final_data_path <- file.path("data_raw", "BRAN2023")

salt <- file.path(raw_data_path, "ocean_salt", "raster_means") %>% 
  list.files(full.names = T) %>% 
  terra::rast()
writeRaster(salt, file.path(final_data_path, "BRAN_ocean_salt.tif"), overwrite = T)

uv <- file.path(raw_data_path, "ocean_uv", "raster_means") %>% 
  list.files(full.names = T) %>% 
  terra::rast()
writeRaster(uv, file.path(final_data_path, "BRAN_ocean_uv.tif"), overwrite = T)

# load(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata")) # distance from coast cell_mask

# dest_path_2 <- file.path(raw_data_path, "ocean_salt", "raster_means")
# dest_path_3 <- file.path("data_raw", "BRAN2023")

# fnms <- list.files(dest_path_2, full.names = T)
# BRAN_salt_data <- terra::rast(fnms)
# BRAN_salt_data <- project(BRAN_salt_data, terra::rast(cell_mask), method = "average")
# writeRaster(x = BRAN_salt_data, filename = file.path(dest_path_3, "BRAN_salt_data.tif"), overwrite = T)
# rm(BRAN_salt_data)

# dest_path_2 <- file.path(raw_data_path, "ocean_u", "raster_means")
# dest_path_3 <- file.path("data_raw", "BRAN2023")

# fnms <- list.files(dest_path_2, full.names = T)
# BRAN_u_data <- terra::rast(fnms)
# BRAN_u_data <- project(BRAN_u_data, terra::rast(cell_mask), method = "average")
# writeRaster(x = BRAN_u_data, filename = file.path(dest_path_3, "BRAN_u_data.tif"), overwrite = T)
# rm(BRAN_u_data)

# dest_path_2 <- file.path(raw_data_path, "ocean_v", "raster_means")
# dest_path_3 <- file.path("data_raw", "BRAN2023")

# fnms <- list.files(dest_path_2, full.names = T)
# BRAN_v_data <- terra::rast(fnms)
# BRAN_v_data <- project(BRAN_v_data, terra::rast(cell_mask), method = "average")
# writeRaster(x = BRAN_v_data, filename = file.path(dest_path_3, "BRAN_v_data.tif"), overwrite = T)
# rm(BRAN_v_data)

# # Testing
# data_sm <- crop(terra::rast(fnms[1]), extent(147.241, 147.974, -43.509, -42.947))
# cell_sm <- crop(terra::rast(cell_mask), extent(147.241, 147.974, -43.509, -42.947))
# 
# data_sm_2 <- project(data_sm, cells_sm, method = "bilinear")
# data_sm_3 <- project(data_sm, cells_sm, method = "average")
# data_sm_4 <- project(data_sm, cells_sm, method = "cubic")
# 
# plot(data_sm)
# plot(data_sm_2)
# plot(data_sm_3)
# plot(data_sm_4)
