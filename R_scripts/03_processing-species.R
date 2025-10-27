# This script takes files from the species targets pipeline and saves them as files, which are easier to load

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
library(here)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)

remove_unit("mol")
install_unit("mol", def = "14.0067 g")

output_path <- here() %>% file.path("data", "processed_species")

Sys.setenv(TAR_PROJECT = "project_species")
tar_names <- tar_manifest()$name

spec_ni_uptake <- c("linear", "linear")
spec_am_uptake <- c("MM", "MM")
species_names <- c("A. armata", "A. taxiformis")

armata <- macrogrow::a_armata
armata["S_min"] <- 21.25
armata["S_max"] <- 41
taxiformis <- armata
taxiformis["T_opt"] <- 24.6153846153846
taxiformis["T_min"] <- 16.0326603325416
taxiformis["T_max"] <- 29.2262470308789

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
                        spec_params = armata,
                        uptake_shape = spec_ni_uptake[1])
conc_2 <- conc %>% mutate(form = "ammonium")
conc_2$uptake <- sapply(X = conc_1$mgm3,
                        FUN = get_uptake,
                        Nform_abbr = "am",
                        spec_params = armata,
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
temp$armata <- sapply(X = temp$temp, FUN = T_lim, spec_params = armata)
temp$taxiformis <- sapply(X = temp$temp, FUN = T_lim, spec_params = taxiformis)
temp <- temp %>%
  pivot_longer(cols = c("armata", "taxiformis"), names_to = "species", values_to = "Tlim") %>%
  write_parquet(file.path(output_path, "T_response.parquet"))
rm(temp)

## Salinity --------------------------------------------------------------------------------------------------
sali <- data.frame(sali = seq(0, 65, 0.25))
sali$armata <- sapply(X = sali$sali, FUN = S_lim, spec_params = armata)
sali$taxiformis <- sapply(X = sali$sali, FUN = S_lim, spec_params = taxiformis)
sali <- sali %>%
  pivot_longer(cols = c("armata", "taxiformis"), names_to = "species", values_to = "Slim") %>%
  write_parquet(file.path(output_path, "S_response.parquet"))
rm(sali)

## Light --------------------------------------------------------------------------------------------------------
light <- data.frame(light = seq(0, 1600, 25))
tar_load(culture_depths)
ls1 <- list()
for (d in 1:length(culture_depths)) {
  light$armata <- sapply(X = light$light, FUN = I_lim, kW = 0.1, Nf = 500, site_params = c(d_top = culture_depths[d], kW = 0.2), spec_params = armata)
  ls1[[d]] <- light %>%
    pivot_longer(cols = c("armata"), names_to = "species", values_to = "Ilim") %>%
    mutate(depth_m = culture_depths[d])
}

ls1 <- ls1 %>%
  bind_rows() %>%
  write_parquet(file.path(output_path, "I_response.parquet"))
rm(light, ls1)

# Sensitivities ---------------------------------------------------------------------------------------------------
tar_load(factors)
## Conditions -----------------------------------------------------------------------------------------------------
tar_load(sens_conditions) 
sens_conditions %>%
  rename(t = t_span) %>%
  write_parquet(file.path(output_path, "sensitivity_conditions.parquet"))

## End calculations -----------------------------------------------------------------------------------------------
### Asparagopsis armata -------------------------------------------------------------------------------------------
tar_load(sens_N_end_arma) 
param_names <- names(armata)[!is.na(armata)]
param_names <- param_names[!param_names %in% c("D_lo", "D_mi", "D_hi", "M_am", "C_am", "V_ni", "K_ni", "M_ot", "C_ot")]

sens_N_end_arma %>% 
  write_parquet(file.path(output_path, "sensitivity_arma_raw.parquet")) %>% 
  filter(adj_param %in% param_names) %>%
  mutate(adj_param = factor(adj_param, levels = param_names), 
         species = factor(species, levels = species_names)) %>% 
  write_parquet(file.path(output_path, "sensitivity_arma_calc_full.parquet"))

### Asparagopsis taxiformis ---------------------------------------------------------------------------------------
tar_load(sens_N_end_taxi)

sens_N_end_taxi %>% 
  write_parquet(file.path(output_path, "sensitivity_taxi_raw.parquet")) %>% 
  filter(adj_param %in% param_names) %>%
  mutate(adj_param = factor(adj_param, levels = param_names), 
         species = factor(species, levels = species_names)) %>% 
  write_parquet(file.path(output_path, "sensitivity_taxi_calc_full.parquet"))

