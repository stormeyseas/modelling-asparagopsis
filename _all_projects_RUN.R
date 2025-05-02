library(targets, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(stringr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(gitcreds, warn.conflicts = F)

# Species ---------------------------------------------------------------------------------------------------------
# The species project collates all the species data into iterable lists for easier model running, and also runs all the sensitivity analyses on species parameters. It doesn't use anything from the other targets projects.
  # - parameterisation/asparagopsis-parameterisation.qmd

## Targets run ----------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_species")

tar_validate()
tar_make(reporter = "balanced", seconds_meta_append = 300)
source("R_scripts/03_processing-species.R")

# Spatial cells ---------------------------------------------------------------------------------------------------
# The spatial-cells project processes the more ad-hoc data for nitrogen, and collates all the environmental data 
overwrite <- F

# into discrete cell-specific (and state-specific) objects that are easier to call for actual runs.
load(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_coords.Rdata")) # cell_coords for BARRA
cell_sample <- NA # cell_coords %>% 
  # filter(layer <= 97.5) %>% 
  # slice_sample(n = 2500)

# BARRA2-R data was downloaded from: [Thredds](https://thredds.nci.org.au/thredds/fileServer/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/day/catalog.html). It requires access to the ob53 project.
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA2")
source("R_scripts/04_extract-BARRA2-data.R")

# BARRA-C2 data was downloaded from: [Thredds](https://thredds.nci.org.au/thredds/catalog/ob53/output/reanalysis/AUST-04/BOM/ERA5/historical/hres/BARRA-C2/v1/day/catalog.html). It requires access to the ob53 project.
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2")
source("R_scripts/05_extract-BARRA-C2-data.R")

# BRAN2023 data was downloaded from: [Thredds](https://thredds.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2023/daily/catalog.html). It requires access to the gb6 project.
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BRAN2023")
source("R_scripts/06_extract-BRAN2023-data.R")

data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "Aqua_MODIS_KD")
source("R_scripts/06_extract-MODIS-data.R")

# AusBathyTopo 2023 data was downloaded from: [Metadata catalogue](https://doi.org/10.26186/148758).
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "AusBathyTopo 2024")
source("R_scripts/07_extract-AusBathyTopo-data.R")

## Targets cell data ----------------------------------------------------------------------------------------------
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.1_model_running_TAS.R", store = "targets_outputs/_model_running_TAS")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.2_model_running_VIC.R", store = "targets_outputs/_model_running_VIC")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.3_model_running_QLD.R", store = "targets_outputs/_model_running_QLD")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.4_model_running_SAU.R", store = "targets_outputs/_model_running_SAU")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.5_model_running_WAN.R", store = "targets_outputs/_model_running_WAN")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.6_model_running_WAS.R", store = "targets_outputs/_model_running_WAS")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.7_model_running_NTE.R", store = "targets_outputs/_model_running_NTE")
tar_make(cell_input_timeseries, seconds_meta_append = 300, script = "R_scripts/07.8_model_running_NSW.R", store = "targets_outputs/_model_running_NSW")

tar_prune()
tar_make(seconds_meta_append = 300)

source("R_scripts/08_processing_cell_data_1.R")


source("R_scripts/08_processing_cell_data_2.R")


# Renv files ------------------------------------------------------------------------------------------------------
library(gitcreds)
library(magrittr)

targets::tar_renv(
  extras = c(
    "bslib", "crew", "gt", "markdown", "rstudioapi", "shiny", "shinybusy", "shinyWidgets", "visNetwork", "qs",
    "qs2", "stormeyseas/macrogrow", "DOI-USGS/streamMetabolizer", "markwh/subsetnc"
  ),
  path = file.path("renv", "packages_spatial_cells.R"),
  script = file.path("R_scripts", "08_spatial_cells.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_model_running.R"),
  script = file.path("R_scripts", "12_model_running.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_species.R"),
  script = file.path("R_scripts", "02_species.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_model_running.R"),
  script = file.path("R_scripts", "12_model_running.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_theo_scens.R"),
  script = file.path("R_scripts", "12.5_theo_scens.R")
)

packs <- renv::dependencies()$Package %>% unique()
renv::install(
  packages = c(packs, "DOI-USGS/streamMetabolizer", "stormeyseas/macrogrow", "markwh/subsetnc", "appling/unitted"),
  exclude = c("streamMetabolizer", "macrogrow", "subsetnc", "unitted"),
  dependencies = T
)
renv::snapshot()



