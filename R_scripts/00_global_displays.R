library(here)
library(tidyverse)
library(qs2)

# This path is for my (Tormey's) personal longterm data storage
permdata_path <- file.path("D:", "modelling-asparagopsis")

spec_store <- file.path(permdata_path, "targets_outputs", "_species")
envi_store <- file.path(permdata_path, "targets_outputs", "_env_inputs")
runs_store <- file.path(permdata_path, "targets_outputs", "_model_running")
spec_data  <- file.path(permdata_path, "data", "processed_species")
envi_data  <- file.path(permdata_path, "data", "processed_env_inputs")
runs_data  <- file.path(permdata_path, "data", "processed_model_running")

bounds <- file.path(envi_data, "states_bbox.qs") %>% qs::qread()

prettyplot <- function() {
  theme_classic() +
    theme(text = element_text(family = "serif", size = 11, colour = "black"),
          legend.position = "none",
          axis.title.y = element_text(vjust = 1.5),
          axis.title.x = element_text(vjust = 1.5),
          legend.title = element_blank())
}
env_plot <- function() {
  prettyplot() + theme(aspect.ratio = 0.55)
}
rm.y <- function() {
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
}

find_read <- function(path, pattern) {
  fnm <- list.files(path, full.names = T) %>% str_subset(pattern)
  if (all(str_detect(fnm, ".qs"))) {dat <- purrr::map(fnm, qs::qread)}
  if (all(str_detect(fnm, ".parquet"))) {dat <- purrr::map_dfr(fnm, arrow::read_parquet)}
  return(dat)
}

nitrogen_lab <- expression("Concentration ( "*mu*"M N)")
nitrate_lab <- expression("Nitrate concentration ( "*mu*"M NO"[3]*")")
ammonium_lab <- expression("Ammonium concentration ( "*mu*"M NH"[4]^+")")

# Allows using "units" package to convert between moles and g of nitrogen ONLY
remove_unit("mol")
install_unit("mol", "14.0067 g")

prettynum <- function(num, sigfig=2) {format(round(as.numeric(num), sigfig), nsmall=1, big.mark=",")}

species_pal <- c("red3", "dodgerblue3")
depths_pal <- c(
  "0.5" = "green3", 
  "2.5" = "#FF8247", 
  "5" = "dodgerblue3", 
  "10" = "grey56"
)
scens_pal <- c(
  "base" = "#292928", 
  "deep_water" = "dodgerblue3", 
  "shallow_bay" = "#66a61e", 
  "fish_farm" = "#e7298a", 
  "estuary_mouth" = "#e6ab02"
)

states_pal <- c(
  "NTE" = "#d95f02", "Northern Territory" = "#d95f02",
  "QLD" = "#e6ab02", "Queensland" = "#e6ab02",
  "NSW" = "#7570b3", "New South Wales" = "#7570b3",
  "SAU" = "#1b9e77", "South Australia" = "#1b9e77",
  "TAS" = "#66a61e", "Tasmania" = "#66a61e",
  "VIC" = "#e7298a", "Victoria" = "#e7298a",
  "WAN" = "#a6761d", "Western Australia (N)" = "#a6761d",
  "WAS" = "#666666", "Western Australia (S)" = "#666666"
)
lims_pal <- c("No_limit" = "black", "T_lim" = "#e7298a", "I_lim" = "#e6ab02", "S_lim" = "#7570b3", 
              "Q_lim" = "#66a61e", "Nf_loss" = "#1b9e77", "Ns_loss" = "#e7298a")
states <- tar_read(states, store = runs_store)
states_ord <- c("NTE", "QLD", "NSW", "SAU", "TAS", "VIC", "WAN", "WAS")
states_lng <- c("Northern Territory", "Queensland", "New South Wales", "South Australia", "Tasmania", "Victoria", "Western Australia (N)", "Western Australia (S)")
