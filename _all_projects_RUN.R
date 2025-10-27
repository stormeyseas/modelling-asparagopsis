# Setup -----------------------------------------------------------------------------------------------------------
library(targets)
library(magrittr)
library(stringr)
library(dplyr)
library(gitcreds)
library(here)
library(macrogrow)
library(ggplot2)

here() %>% file.path("R_scripts", "00_functions.R") %>% source()

# Species ---------------------------------------------------------------------------------------------------------
# The species markdown details the parameterisation process for _Asparagopsis armata_ and _Asparagopsis taxiformis_. 
# It collates all the species parameters into named number vectors for easier model running, and also outputs them 
# into human-readable .csv files in the "data" folder. 
  # terminal: quarto render parameterisation/asparagopsis-parameterisation.qmd

# The project_species targets pipeline collates the relevant species data into iterable lists for easier model 
# running later on, and also runs all the sensitivity analyses on species parameters. 
Sys.setenv(TAR_PROJECT = "project_species")
tar_make(
  explore_starting_biomass,
  seconds_meta_append = 300,
  callr_function = NULL
)
tar_load(explore_starting_biomass)
explore_starting_biomass %>% 
  ggplot(aes(x = Nf, y = N_rem)) +
  geom_line()
tar_prune()

tar_make(
  seconds_meta_append = 300,
  callr_function = NULL
)
source("R_scripts/03_processing-species.R")

# Environmental data ----------------------------------------------------------------------------------------------
# The spatial-cells project processes the ad-hoc data for nitrogen, and collates all the environmental data into 
# discrete cell-specific (and state-specific) objects that are easier to call for actual runs.

# Running these scripts requires that the user has downloaded the relevant raw data into the following path:
big_path <- file.path("D:", "FRDC-Seaweed-Raw-Data")

load(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA2", "BARRA2_cell_coords.Rdata")) # cell_coords for BARRA

# BARRA2 data was downloaded from Thredds: 
#   https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f9057_2475_0540_0329
# It requires access to the ob53 project.
raw_data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-R2")
final_data_path <- here() %>% file.path("data_raw", "BARRA-R2")
overwrite <- F
source("R_scripts/05_extract-BARRA-C2-data.R")

# BRAN2023 data was downloaded from: [Thredds](https://thredds.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2023/daily/catalog.html). It requires access to the gb6 project.
raw_data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "BRAN2023")
source("R_scripts/05_extract-BRAN2023-data.R")

# Aqua MODIS data was downloaded as png files from: [MODIS NASA](http://oceancolor.gsfc.nasa.gov/l3)
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "Aqua_MODIS_KD")
source("R_scripts/06_extract-MODIS-data.R")

# AusBathyTopo 2023 data was downloaded from: [Metadata catalogue](https://doi.org/10.26186/148758).
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "AusBathyTopo 2024")
source("R_scripts/07_extract-AusBathyTopo-data.R")

## Spatial data targets -------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_inputs")

tar_outdated()
tar_make(
  seconds_meta_append = 120,
  callr_function = NULL
)
tar_prune()
source("R_scripts/09_processing_env_inputs.R")

# Model running ---------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_running")

tar_outdated()

# Quicker test
tar_invalidate("test_growth_armata")
tar_make(test_growth_armata)
tar_load(test_growth_armata)
ggplot(test_growth_armata, aes(x = t, y = B_dw.mg, color = as.factor(start))) +
  geom_line(linewidth = 0.75) +
  facet_grid(rows = vars(state)) +
  theme_classic()
ggplot(test_growth_armata, aes(x = t, y = (Nf+Ns)*hm, color = as.factor(start))) +
  geom_line(linewidth = 0.75) +
  facet_grid(rows = vars(state)) +
  theme_classic()

tar_make(
  seconds_meta_append = 90,
  callr_function = NULL
)
tar_prune()

source("R_scripts/11_processing_model_running.R")

# Renv files ------------------------------------------------------------------------------------------------------
library(gitcreds)
library(magrittr)

targets::tar_renv(
  extras = c("bslib", "crew", "gt", "markdown", "rstudioapi", "shiny", "shinybusy", "shinyWidgets", "visNetwork", "qs", "qs2", "rbenchmark"),
  path = file.path("renv", "packages_inputs.R"),
  script = file.path("targets_pipelines", "env_inputs.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_species.R"),
  script = file.path("targets_pipelines", "species.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_running.R"),
  script = file.path("targets_pipelines", "model_running.R")
)

packs <- renv::dependencies()$Package %>% unique()
renv::install(
  packages = c(packs, "stormeyseas/macrogrow", "DOI-USGS/streamMetabolizer", "markwh/subsetnc", "ropensci/geotargets"),
  exclude = c("macrogrow", "streamMetabolizer", "subsetnc", "geotargets", "RPostgreSQL"),
  dependencies = T
)
renv::snapshot()



