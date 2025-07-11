library(targets)
library(magrittr)
library(stringr)
library(dplyr)
library(gitcreds)
library(here)
library(macrogrow)

here() %>% file.path("R_scripts", "00_functions.R") %>% source()

# Species ---------------------------------------------------------------------------------------------------------
# The species markdown details the parameterisation process for _Asparagopsis armata_ and _Asparagopsis taxiformis_. 
# It collates all the species parameters into named number vectors for easier model running, and also outputs them 
# into human-readable .csv files in the "data" folder. 
  # terminal: quarto render parameterisation/asparagopsis-parameterisation.qmd

# The project_species targets pipeline collates the relevant species data into iterable lists for easier model 
# running later on, and also runs all the sensitivity analyses on species parameters. 
Sys.setenv(TAR_PROJECT = "project_species")
tar_make(seconds_meta_append = 300)
source("R_scripts/03_processing-species.R")

# Environmental data ----------------------------------------------------------------------------------------------
# The spatial-cells project processes the ad-hoc data for nitrogen, and collates all the environmental data into 
# discrete cell-specific (and state-specific) objects that are easier to call for actual runs.

# Running these scripts requires that the user has downloaded the relevant raw data into the following path:
big_path <- file.path("D:", "FRDC-Seaweed-Raw-Data")

load(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_coords.Rdata")) # cell_coords for BARRA
cell_sample <- NA

# BARRA2-R data was downloaded from Thredds: 
#   https://thredds.nci.org.au/thredds/fileServer/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/day/catalog.html
# It requires access to the ob53 project.
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA2")
source("R_scripts/04_extract-BARRA2-data.R")

# BARRA-C2 data was downloaded from Thredds: 
#   https://thredds.nci.org.au/thredds/catalog/ob53/output/reanalysis/AUST-04/BOM/ERA5/historical/hres/BARRA-C2/v1/day/catalog.html
# It requires access to the ob53 project.
raw_data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2")
final_data_path <- here() %>% file.path("data_raw", "BARRA-C2")
overwrite <- F
source("R_scripts/05_extract-BARRA-C2-data.R")

# BRAN2023 data was downloaded from: [Thredds](https://thredds.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2023/daily/catalog.html). It requires access to the gb6 project.
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BRAN2023")
source("R_scripts/06_extract-BRAN2023-data.R")

data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "Aqua_MODIS_KD")
source("R_scripts/07_extract-MODIS-data.R")

# AusBathyTopo 2023 data was downloaded from: [Metadata catalogue](https://doi.org/10.26186/148758).
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "AusBathyTopo 2024")
source("R_scripts/08_extract-AusBathyTopo-data.R")

## Spatial data targets -------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_inputs")

tar_outdated()
tar_make(seconds_meta_append = 90)
source("R_scripts/10_processing_env_inputs.R")

# Model running ---------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_running")

targets::tar_outdated()
targets::tar_make(seconds_meta_append = 190)
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



