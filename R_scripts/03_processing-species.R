library(tidyr)
library(dplyr)
library(magrittr)
library(stringr)
library(units)
library(targets)
library(tarchetypes)
library(macrogrow)
library(ozmaps)
library(arrow)
library(gganimate)
library(conflicted)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)

remove_unit("mol")
install_unit("mol", def = "14.0067 g")

output_path <- file.path("data", "processed_species")

Sys.setenv(TAR_PROJECT = "project_species")
meta <- tar_meta(fields = c("names", "parent", "time")) %>% 
  group_by(parent) %>% 
  reframe(time = max(time),
          time = as.Date(time)) %>% 
  filter(time == as.Date(lubridate::now()))
mani <- tar_manifest()

arma <- tar_read(a_armata_gametophyte)
taxi <- tar_read(a_taxiformis_gametophyte)

spec_ni_uptake <- c("linear", "linear")
spec_am_uptake <- c("MM", "MM")
species_names <- c("A. armata", "A. taxiformis")

# Responses -------------------------------------------------------------------------------------------------------
## Nitrogen -----------------------------------------------------------------------------------------------------
conc <- data.frame(uM = seq(0, 25, 0.25)) %>%
  mutate(uM = set_units(uM, "umol L-1"),
         mgm3 = set_units(uM, "mg m-3")) %>%
  mutate(uM = drop_units(uM),
         mgm3 = drop_units(mgm3))

conc_1 <- conc %>% mutate(form = "nitrate")
conc_1$uptake <- sapply(X = conc_1$mgm3,
                        FUN = get_uptake,
                        Nform_abbr = "ni",
                        spec_params = arma,
                        uptake_shape = spec_ni_uptake[1])
conc_2 <- conc %>% mutate(form = "ammonium")
conc_2$uptake <- sapply(X = conc_1$mgm3,
                        FUN = get_uptake,
                        Nform_abbr = "am",
                        spec_params = arma,
                        uptake_shape = spec_am_uptake[1])

conc <- rbind(conc_1, conc_2) %>%
  mutate(form = factor(form, levels = c("nitrate", "ammonium")),
         uptake_mgm3 = set_units(uptake, "mg m-3 d-1"),
         uptake_uM = set_units(uptake_mgm3, "umol L-1 d-1")) %>%
  dplyr::select(-uptake) %>%
  mutate(uptake_mgm3 = drop_units(uptake_mgm3),
         uptake_uM = drop_units(uptake_uM)) %>%
  write_parquet(file.path(output_path, "N_uptake.parquet"))
rm(conc, conc_1, conc_2)

## Temperature --------------------------------------------------------------------------------------------------
temp <- data.frame(temp = seq(-1, 40, 0.5))
temp$arma <- sapply(X = temp$temp, FUN = T_lim, spec_params = arma)
temp$taxi <- sapply(X = temp$temp, FUN = T_lim, spec_params = taxi)
temp <- temp %>%
  pivot_longer(cols = c("arma", "taxi"), names_to = "species", values_to = "Tlim") %>%
  write_parquet(file.path(output_path, "T_response.parquet"))
rm(temp)

## Salinity --------------------------------------------------------------------------------------------------
sali <- data.frame(sali = seq(0, 65, 0.25))
sali$arma <- sapply(X = sali$sali, FUN = S_lim, spec_params = arma)
sali$taxi <- sapply(X = sali$sali, FUN = S_lim, spec_params = taxi)
sali <- sali %>%
  pivot_longer(cols = c("arma", "taxi"), names_to = "species", values_to = "Slim") %>%
  write_parquet(file.path(output_path, "S_response.parquet"))
rm(sali)

## Light --------------------------------------------------------------------------------------------------------
light <- data.frame(light = seq(0, 1600, 25))
culture_depths <- tar_read(culture_depths)
ls1 <- list()
for (d in 1:length(culture_depths)) {
  light$arma <- sapply(X = light$light, FUN = I_lim, Nf = 500, site_params = c(d_top = culture_depths[d], kW = 0.2), spec_params = arma)
  ls1[[d]] <- light %>%
    pivot_longer(cols = c("arma"), names_to = "species", values_to = "Ilim") %>%
    mutate(depth_m = culture_depths[d])
}

ls1 <- ls1 %>%
  bind_rows() %>%
  write_parquet(file.path(output_path, "I_response.parquet"))
rm(light, ls1)

# Sensitivities ---------------------------------------------------------------------------------------------------
factors <- tar_read(factors)
## Conditions -----------------------------------------------------------------------------------------------------
tar_read(sens_conditions) %>%
  rename(t = t_span) %>%
  write_parquet(file.path(output_path, "sensitivity_conditions.parquet"))

## End calculations -----------------------------------------------------------------------------------------------
### Asparagopsis armata -------------------------------------------------------------------------------------------
s_arma <- tar_read(total_N_end_arma) 
param_names <- names(arma)[!is.na(arma)]
param_names <- param_names[!param_names %in% c("D_lo", "D_mi", "D_hi", "M_am", "C_am", "V_ni", "K_ni", "M_ot", "C_ot")]

s_arma <- s_arma %>% 
  write_parquet(file.path(output_path, "sensitivity_arma_raw.parquet")) %>% 
  filter(param %in% param_names & !is.na(TN)) %>%
  mutate(param = factor(param, levels = param_names), 
         factor = as.character(factor), 
         T_level = as.factor(T_level),
         species = factor(species, levels = species_names)) 

s_arma <- s_arma %>% 
  pivot_wider(names_from = factor, names_prefix = "p_", values_from = TN) %>%
  mutate(sens = (p_1.05 - p_0.95) / (0.1 * p_1)) %>% 
  write_parquet(file.path(output_path, "sensitivity_arma_calc_full.parquet"))

# Summarise across dates within T_levels
s_arma <- s_arma %>% 
  group_by(species, cult_dep, param, T_level) %>% 
  reframe(sd = sd(sens, na.rm = T), sens = mean(sens, na.rm = T)) %>% 
  write_parquet(file.path(output_path, "sensitivity_arma_calc_sepTU.parquet"))

# Summarise across dates and T_levels
s_arma <- s_arma %>% 
  group_by(species, cult_dep, param) %>% 
  reframe(sd = sd(sens, na.rm = T), sens = mean(sens, na.rm = T)) %>% 
  write_parquet(file.path(output_path, "sensitivity_arma_calc_combTU.parquet"))

### Asparagopsis taxiformis ---------------------------------------------------------------------------------------
s_taxi <- tar_read(total_N_end_taxi)
s_taxi <- s_taxi %>% 
  write_parquet(file.path(output_path, "sensitivity_taxi_raw.parquet")) %>% 
  filter(param %in% param_names & !is.na(TN)) %>%
  mutate(param = factor(param, levels = param_names), 
         factor = as.character(factor), 
         T_level = as.factor(T_level),
         species = factor(species, levels = species_names)) 

s_taxi <- s_taxi %>% 
  pivot_wider(names_from = factor, names_prefix = "p_", values_from = TN) %>%
  mutate(sens = (p_1.05 - p_0.95) / (0.1 * p_1)) %>% 
  write_parquet(file.path(output_path, "sensitivity_taxi_calc_full.parquet"))

# Summarise across dates within T_levels
s_taxi <- s_taxi %>% 
  group_by(species, cult_dep, param, T_level) %>% 
  reframe(sd = sd(sens, na.rm = T), sens = mean(sens, na.rm = T)) %>% 
  write_parquet(file.path(output_path, "sensitivity_taxi_calc_sepTU.parquet"))

# Summarise across dates and T_levels
s_taxi <- s_taxi %>% 
  group_by(species, cult_dep, param) %>% 
  reframe(sd = sd(sens, na.rm = T), sens = mean(sens, na.rm = T)) %>% 
  write_parquet(file.path(output_path, "sensitivity_taxi_calc_combTU.parquet"))



