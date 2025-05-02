# Load required packages
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
  library(terra)
  library(qs)
  library(qs2)
  library(future)
  library(furrr)
}))

# Setup parallel processing
plan(multisession, workers = 5)
# plan(multisession, workers = availableCores() - 1)

# Function to clear memory
clean_memory <- function() {
  gc()
  terra::terraOptions(memfrac=0.5)
}

conflicts_prefer(dplyr::filter(), dplyr::select())

fixnum <- function(n, digits = 4) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}

# Cell information ------------------------------------------------------------------------------------------------
if (!file.exists(file.path(data_path, "AquaMODIS_cell_coords.Rdata")) | overwrite == T) {
  assign("BARRA_cell_mask", get(load(file = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata")))) # cell_mask, RasterLayer
  rm(cell_mask)
    
  tmp <- file.path(data_path, "nc_downloads") %>% list.files(full.names = T) %>% str_subset(".nc") %>% first() %>% nc_open()
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

# Get variables ---------------------------------------------------------------------------------------------------
load(file = file.path(data_path, "AquaMODIS_lats.Rdata")) # lat
load(file = file.path(data_path, "AquaMODIS_lons.Rdata")) # lon
load(file = file.path(data_path, "AquaMODIS_lat_inds.Rdata")) # lat_ind
load(file = file.path(data_path, "AquaMODIS_lon_inds.Rdata")) # lon_ind

# Final rasters ---------------------------------------------------------------------------------------------------
load(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata")) # distance from coast cell_mask

dest_path <- file.path(data_path, "nc_downloads")

fnms <- list.files(dest_path, full.names = T)
var_rast <- rast(fnms)
crs(var_rast) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

var_rast <- project(var_rast, terra::rast(cell_mask), method = "average")
ext <- extent(cell_mask)
var_rast <- crop(var_rast, ext)
writeRaster(x = var_rast, filename = file.path(data_path, "AquaMODIS_Kd_490_data.tif"), overwrite = T)
rm(var_rast)

