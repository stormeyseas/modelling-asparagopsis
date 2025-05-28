library(targets, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(stringr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(gitcreds, warn.conflicts = F)
fixnum <- function(n, digits = 4) {
  vapply(n, function(x) {
    str_flatten(c(rep("0", digits-nchar(as.character(x))), as.character(x)))
  }, character(1))
}

# Species ---------------------------------------------------------------------------------------------------------
# The species project collates all the species data into iterable lists for easier model running, and also runs all the sensitivity analyses on species parameters. It doesn't use anything from the other targets projects.
  # - parameterisation/asparagopsis-parameterisation.qmd

## Targets run ----------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_species")

tar_outdated()
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
source("R_scripts/07_extract-MODIS-data.R")

# AusBathyTopo 2023 data was downloaded from: [Metadata catalogue](https://doi.org/10.26186/148758).
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "AusBathyTopo 2024")
source("R_scripts/08_extract-AusBathyTopo-data.R")

## Targets run ----------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_inputs")

tar_outdated()
tar_make(seconds_meta_append = 90)
source("R_scripts/10_processing_env_inputs.R")

# Model running ---------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_running")

tar_outdated()
tar_make(
  names = c("cell_growth_extra_lims", "cell_growth_extra", "cell_growth_extra_end"), 
  shortcut = T, 
  seconds_meta_append = 190
)
source("R_scripts/12_processing_model_running.R")

# Renv files ------------------------------------------------------------------------------------------------------
library(gitcreds)
library(magrittr)

targets::tar_renv(
  extras = c("bslib", "crew", "gt", "markdown", "rstudioapi", "shiny", "shinybusy", "shinyWidgets", "visNetwork", "qs", "qs2"),
  path = file.path("renv", "09_env_inputs.R"),
  script = file.path("R_scripts", "09_env_inputs.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_species.R"),
  script = file.path("R_scripts", "02_species.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_running.R"),
  script = file.path("R_scripts", "11_model_running.R")
)

packs <- renv::dependencies()$Package %>% unique()
renv::install(
  packages = c(packs, "stormeyseas/macrogrow", "DOI-USGS/streamMetabolizer", "markwh/subsetnc", "ropensci/geotargets"),
  exclude = c("macrogrow", "streamMetabolizer", "subsetnc", "geotargets"),
  dependencies = T
)
renv::snapshot()



