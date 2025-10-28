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

# Cell information ------------------------------------------------------------------------------------------------
if (!file.exists(file.path(data_path, "AquaMODIS_cell_coords.Rdata")) | overwrite == T) {
  assign("BARRA_cell_mask", get(load(file = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2", "BARRA_R2_cell_mask.Rdata")))) # cell_mask, RasterLayer
  rm(cell_mask)
    
  tmp <- file.path(data_path, "nc_downloads") %>% 
    list.files(full.names = T) %>% 
    first() %>% 
    nc_open()
  
  # Kd_490[lon,lat]
  lon <- ncvar_get(tmp, "lon")
  lat <- ncvar_get(tmp, "lat")
  
  lon_ind <- which(lon >=  108 & lon <= 155)
  lat_ind <- which(lat >=  -44.99 & lat <= -8.99)
  
  save(lat, file = file.path(data_path, "AquaMODIS_lats.Rdata"))
  save(lon, file = file.path(data_path, "AquaMODIS_lons.Rdata"))
  save(lat_ind, file = file.path(data_path, "AquaMODIS_lat_inds.Rdata"))
  save(lon_ind, file = file.path(data_path, "AquaMODIS_lon_inds.Rdata"))
  
  nc_close(tmp)
  rm(tmp)
}
clean_memory()

# Get data --------------------------------------------------------------------------------------------------------
load(file = file.path(data_path, "AquaMODIS_lats.Rdata")) # lat
load(file = file.path(data_path, "AquaMODIS_lons.Rdata")) # lon
load(file = file.path(data_path, "AquaMODIS_lat_inds.Rdata")) # lat_ind
load(file = file.path(data_path, "AquaMODIS_lon_inds.Rdata")) # lon_ind

nc_files <- file.path(data_path, "nc_downloads")
dest_path <- file.path(data_path, "raster_means")

# Create new filenames for each day (no leap years)
rast_files <- data.frame(
  fnm = list.files(nc_files, full.names = T),
  time = list.files(nc_files)
) %>% 
  mutate(
    time = str_split_i(time, "[.]", 2),
    time = lubridate::parse_date_time(time, orders = "%Y%m%d"),
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
  fname <- file.path(dest_path, str_c("Kd_490_doy_", fixnum(d), ".tif"))

  rf <- rast_files[rast_files$leap_day == d, ]
  
  ras_stack <- terra::rast(rf$fnm)
  ras_mean <- app(ras_stack, fun = mean, na.rm = TRUE) %>% 
    setNames(str_c("doy_", d))
  
  crs(ras_mean) <- "EPSG:4326"

  writeRaster(x = ras_mean, filename = fname, overwrite = T)
  clean_memory()
  
  print(str_c("Kd_490 file DOY ", fixnum(d), " saved"))
}

# Final raster ---------------------------------------------------------------------------------------------------
cell_mask_file <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2", "BARRA_R2_cell_mask.Rdata")
load(cell_mask_file) # loads 'cell_mask' object
cell_mask_raster <- terra::rast(cell_mask)

raster_means_path <- file.path(data_path, "raster_means")
raster_files <- list.files(
  raster_means_path,
  pattern = "Kd_490_doy_.*\\.tif$",
  full.names = TRUE
)

modis_raster_stack <- terra::rast(raster_files)
names(modis_raster_stack)

reprojected_raster <- terra::project(modis_raster_stack, cell_mask_raster, method = "average")

terra::writeRaster(
  x = reprojected_raster,
  filename = file.path(data_path, "AquaMODIS_Kd_490_data.tif"),
  overwrite = TRUE
)

rm(reprojected_raster, modis_raster_stack, cell_mask_raster)
clean_memory()
