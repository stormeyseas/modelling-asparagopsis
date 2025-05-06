suppressMessages(suppressWarnings({
  library(magrittr)
  library(arrow)
  library(dplyr)
  library(targets)
  library(stringr)
  library(geosphere)
  library(ncdf4)
  library(reshape2)
  library(terra)
  library(subsetnc)
  library(sf)
  library(ozmaps)
  library(raster)
  library(ggplot2)
  library(conflicted)
  library(qs)
  library(qs2)
}))
conflicts_prefer(dplyr::select(), dplyr::filter(), .quiet = T)

# Load BARRA-C2 cell_mask, RasterLayer
load(file = file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata"))
cell_mask <- terra::rast(cell_mask)

bathy <- file.path("D:", "FRDC-Seaweed-Raw-Data", "AusBathyTopo 2024", "AusBathyTopo__Australia__2024_250m_MSL_cog.tif") %>% terra::rast()
crs(bathy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
bathy <- crop(bathy, ext(cell_mask))
bathy[bathy >= 0] <- NA

bathy_projected <- project(bathy, cell_mask, method = "average")
writeRaster(x = bathy_projected, filename = file.path("data_raw", "AusBathyTopo 2024", "bathy_projected.tif"), overwrite = T)
rm(bathy_projected)

# # Testing
# bathy_sm <- crop(bathy, extent(147.241, 147.974, -43.509, -42.947))
# cells_sm <- crop(cell_mask, extent(147.241, 147.974, -43.509, -42.947))
# 
# bathy_sm_1 <- project(bathy_sm, cells_sm)
# bathy_sm_2 <- project(bathy_sm, cells_sm, method = "bilinear")
# bathy_sm_3 <- project(bathy_sm, cells_sm, method = "average")
# bathy_sm_4 <- project(bathy_sm, cells_sm, method = "cubic")
# 
# plot(bathy_sm)
# plot(bathy_sm_1)
# plot(bathy_sm_2)
# plot(bathy_sm_3)
# plot(bathy_sm_4)

