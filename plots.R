library(stringr)
library(units)
library(targets)
library(tarchetypes)
library(macrogrow)
library(ozmaps)
library(arrow)
library(gganimate)
library(magrittr)
library(ggplot2)
library(dplyr)
library(terra)
library(tidyterra)
library(here)
library(lubridate)
library(tidyr)
library(cowplot)
library(conflicted)
conflicts_prefer(dplyr::select, dplyr::filter)

prettyplot <- theme_classic() +
  theme(text = element_text(family = "sans", size = 12, colour = "black"),
        axis.title = element_text(vjust = 2),
        legend.position = "none")

nitrogen_lab <- expression("Concentration ("*mu*"M N)")
nitrate_lab <- expression("Nitrate concentration ("*mu*"M NO"[3]*")")
ammonium_lab <- expression("Ammonium concentration ("*mu*"M NH"[4]^+")")

# Allows using "units" package to convert between moles and g of nitrogen ONLY
remove_unit("mol")
install_unit("mol", "14.0067 g")

prettynum <- function(num, sigfig=2) {format(round(as.numeric(num), sigfig), nsmall=1, big.mark=",")}

species_pal <- c("red3", "dodgerblue3")
depths_pal <- c("green3", "#FF8247", "dodgerblue3", "#292928")
states_pal <- c("#d95f02", "#e6ab02", "#7570b3", "#1b9e77", "#66a61e", "#e7298a", "#a6761d", "#666666")
scens_pal_1 <- c("base" = "#292928", "deep_water" = "dodgerblue3", "shallow_bay" = "#66a61e", "fish_farm" = "#e7298a", "estuary_mouth" = "#e6ab02")
scens_pal_2 <- c("dodgerblue3", "#66a61e", "#e7298a", "#e6ab02")
states_pal_2 <- c(states_pal, "black")
states_ord <- c("NTE", "QLD", "NSW", "SAU", "TAS", "VIC", "WAN", "WAS")
states_lng <- c("Northern Territory", "Queensland", "New South Wales", "South Australia", "Tasmania", "Victoria", "Western Australia (N)", "Western Australia (S)")

cell_coords <- arrow::read_parquet("data/processed_env_inputs/BARRA_C2_cell_coords.parquet") %>% 
  mutate(state = factor(state, levels = states_ord)) %>% 
  dplyr::select(-layer)



# Growth scenarios ------------------------------------------------------------------------------------------------
scen_growth_sp1 <- "data/processed_model_running" %>% 
  list.files(full.names = T) %>% 
  str_subset("scen_growth") %>% 
  str_subset("sp1") %>% 
  purrr::map_dfr(qs::qread)
scen_growth_sp2 <- "data/processed_model_running" %>% 
  list.files(full.names = T) %>% 
  str_subset("scen_growth") %>% 
  str_subset("sp2") %>% 
  purrr::map_dfr(qs::qread)

scen_growth_sp1 %>% 
  group_by(state, scen) %>% 
  reframe(sd = sdna(TN),
          TN = meanna(TN)) %>% 
  ggplot(aes(fill = scen, y = TN, x = state)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = scens_pal_1)











