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
source("R_scripts/07_extract-MODIS-data.R")

# AusBathyTopo 2023 data was downloaded from: [Metadata catalogue](https://doi.org/10.26186/148758).
data_path <- file.path("D:", "FRDC-Seaweed-Raw-Data", "AusBathyTopo 2024")
source("R_scripts/08_extract-AusBathyTopo-data.R")

## Targets cell data ----------------------------------------------------------------------------------------------
projects <- tibble::tribble(
  ~state, ~sc, ~st, ~upd,
  "QLD",  "R_scripts/07.3_model_running_QLD.R", "targets_outputs/_model_running_QLD", NA, 
  "TAS",  "R_scripts/07.1_model_running_TAS.R", "targets_outputs/_model_running_TAS", NA, 
  "SAU",  "R_scripts/07.4_model_running_SAU.R", "targets_outputs/_model_running_SAU", NA, 
  "WAS",  "R_scripts/07.6_model_running_WAS.R", "targets_outputs/_model_running_WAS", NA, 
  "NSW",  "R_scripts/07.8_model_running_NSW.R", "targets_outputs/_model_running_NSW", NA, 
  "VIC",  "R_scripts/07.2_model_running_VIC.R", "targets_outputs/_model_running_VIC", NA, 
  "WAN",  "R_scripts/07.5_model_running_WAN.R", "targets_outputs/_model_running_WAN", NA, 
  "NTE",  "R_scripts/07.7_model_running_NTE.R", "targets_outputs/_model_running_NTE", NA
)
sc <- "R_scripts/09_model_running.R"

### Prune ---------------------------------------------------------------------------------------------------------
for (i in 1:nrow(projects)) {
  projects$state[i] %>% qsave("targets_outputs/this_state.qs")
  try({
    tar_prune(script = sc, store = projects$st[i])
    projects$upd[i] <- 1
  })
}

### Ilim only -----------------------------------------------------------------------------------------------------
# Targets needed for 10_processing_Ilim_data.R
# Uncomment slice_sample(n = 2250) line in BARRA_C2_cell_nos
projects <- projects %>% dplyr::filter(is.na(upd))
for (i in 1:nrow(projects)) {
  projects$state[i] %>% qsave("targets_outputs/this_state.qs")
  try({
    tar_make(
      names = "Ilim_cell_PL", #shortcut = T, 
      seconds_meta_append = 90, script = sc, store = projects$st[i]
      )
    projects$upd[i] <- 1
    })
}
source("R_scripts/10_processing_Ilim_data.R") # I_lim depth testing

### All other targets ---------------------------------------------------------------------------------------------
# Comment out the slice_sample(n = 2250) line in BARRA_C2_cell_nos
# DO NOT run Ilim_cell_PL - maybe comment that out too

tar_workspace(total_cell_growth_234a709e9f8ad1f0, store = "targets_outputs/_model_running_QLD")

projects <- projects %>% dplyr::filter(is.na(upd))
for (i in 1:nrow(projects)) {
  projects$state[i] %>% qsave("targets_outputs/this_state.qs")
  try({
    tar_make(
      names = c("total_cell_growth", "growth_lims", "theo_scens", "theo_site_params"),
      seconds_meta_append = 90, shortcut = F, script = sc, store = projects$st[i]
    )
    projects$upd[i] <- 1
  })
}
source("R_scripts/11_processing_model_running.R")

# Renv files ------------------------------------------------------------------------------------------------------
library(gitcreds)
library(magrittr)

targets::tar_renv(
  extras = c("bslib", "crew", "gt", "markdown", "rstudioapi", "shiny", "shinybusy", "shinyWidgets", "visNetwork", "qs", "qs2"),
  path = file.path("renv", "packages_model_running.R"),
  script = file.path("R_scripts", "09_model_running.R")
)

targets::tar_renv(
  path = file.path("renv", "packages_species.R"),
  script = file.path("R_scripts", "02_species.R")
)

packs <- renv::dependencies()$Package %>% unique()
renv::install(
  packages = c(packs, "stormeyseas/macrogrow", "DOI-USGS/streamMetabolizer", "markwh/subsetnc", "ropensci/geotargets"),
  exclude = c("macrogrow", "streamMetabolizer", "subsetnc", "geotargets"),
  dependencies = T
)
renv::snapshot()



