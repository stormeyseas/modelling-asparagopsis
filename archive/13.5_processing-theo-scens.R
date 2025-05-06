library(targets, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(stringr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(tidyr, warn.conflicts = F)
library(arrow, warn.conflicts = F)
library(conflicted, warn.conflicts = F)
library(ggplot2, warn.conflicts = F)
conflicts_prefer(dplyr::select(), dplyr::mutate(), dplyr::filter(), .quiet = T)

# This script takes files from the model_running project and saves them as files, which are a little easier to load than targets objects
# mani <- tar_manifest(script = file.path("_targets-scripts", "_model_running.R"))$name %>% 
#   str_subset("_file", negate = T)

base_path <- file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "FRDC-seaweed")
out_path <- file.path(base_path, "data", "processed-theo-scens")

Sys.setenv(TAR_PROJECT = "project_theo_scens")

# New -------------------------------------------------------------------------------------------------------------
theo_scen_grow_ulva <- tar_read(theo_scen_grow_ulva) %>% 
  select(-c(hm, note))

theo_scen_grow_ulva <- theo_scen_grow_ulva %>% 
  arrange(-value) %>% 
  group_by(start_date, species, state, scen) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(diff = value - start_value,
         date = as.Date(date),
         start_date = as.Date(start_date),
         days = as.integer(date - start_date),
         state = as.factor(state),
         scen = as.factor(scen)) %>% 
  relocate(start_date, .before = date) %>% 
  relocate(start_value, .after = start_date) %>% 
  write_parquet(file.path(out_path, "theo_scen_grow_ulva.parquet"))

theo_scen_grow_eckl <- tar_read(theo_scen_grow_eckl) %>% 
  select(-c(hm, note))

theo_scen_grow_eckl <- theo_scen_grow_eckl %>% 
  arrange(-value) %>% 
  group_by(start_date, species, state, scen) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(diff = value - start_value,
         date = as.Date(date),
         start_date = as.Date(start_date),
         days = as.integer(date - start_date),
         state = as.factor(state),
         scen = as.factor(scen)) %>% 
  relocate(start_date, .before = date) %>% 
  relocate(start_value, .after = start_date) %>% 
  write_parquet(file.path(out_path, "theo_scen_grow_eckl.parquet"))

theo_scen_grow_arma <- tar_read(theo_scen_grow_arma) %>% 
  select(-c(hm, note))

theo_scen_grow_arma <- theo_scen_grow_arma %>% 
  arrange(-value) %>% 
  group_by(start_date, species, state, scen) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(diff = value - start_value,
         date = as.Date(date),
         start_date = as.Date(start_date),
         days = as.integer(date - start_date),
         state = as.factor(state),
         scen = as.factor(scen)) %>% 
  relocate(start_date, .before = date) %>% 
  relocate(start_value, .after = start_date) %>% 
  write_parquet(file.path(out_path, "theo_scen_grow_arma.parquet"))

theo_scen_grow_taxi <- tar_read(theo_scen_grow_taxi) %>% 
  select(-c(hm, note))

theo_scen_grow_taxi <- theo_scen_grow_taxi %>% 
  arrange(-value) %>% 
  group_by(start_date, species, state, scen) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(diff = value - start_value,
         date = as.Date(date),
         start_date = as.Date(start_date),
         days = as.integer(date - start_date),
         state = as.factor(state),
         scen = as.factor(scen)) %>% 
  relocate(start_date, .before = date) %>% 
  relocate(start_value, .after = start_date) %>% 
  write_parquet(file.path(out_path, "theo_scen_grow_taxi.parquet"))

all <- rbind(theo_scen_grow_arma, theo_scen_grow_taxi, theo_scen_grow_eckl, theo_scen_grow_ulva) %>% 
  select(-c(value, start_value)) %>% 
  mutate(species = as.factor(species)) %>% 
  write_parquet(file.path(out_path, "theo_scen_grow_all.parquet"))


# Scenario conditions ---------------------------------------------------------------------------------------------
theo_T_input <- tar_read(theo_T_input) %>% 
  mutate(state = as.factor(state),
         theo_scen = as.factor(theo_scen)) %>% 
  write_parquet(file.path(out_path, "theo_T_input.parquet"))

# theo_T_input %>% 
#   group_by(state, theo_scen) %>% 
#   reframe(max = max(value),
#           min = min(value),
#           mean = mean(value))
# 
# ggplot(theo_T_input, aes(x = yday, y = value, color = theo_scen)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))


theo_S_input <- tar_read(theo_S_input) %>% 
  mutate(state = as.factor(state),
         theo_scen = as.factor(theo_scen)) %>% 
  write_parquet(file.path(out_path, "theo_S_input.parquet"))

# theo_S_input %>% 
#   group_by(state, theo_scen) %>% 
#   reframe(max = max(value),
#           min = min(value),
#           mean = mean(value))
# 
# ggplot(theo_S_input, aes(x = yday, y = value, color = theo_scen)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))

theo_I_input <- tar_read(theo_I_input) %>% 
  mutate(state = as.factor(state),
         theo_scen = as.factor(theo_scen)) %>% 
  write_parquet(file.path(out_path, "theo_I_input.parquet"))

# theo_I_input %>% 
#   group_by(state, theo_scen) %>% 
#   reframe(max = max(value),
#           min = min(value),
#           mean = mean(value))
# 
# ggplot(theo_I_input, aes(x = yday, y = value, color = theo_scen)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))

theo_UV_input <- tar_read(theo_UV_input) %>% 
  mutate(state = as.factor(state),
         theo_scen = as.factor(theo_scen)) %>% 
  write_parquet(file.path(out_path, "theo_UV_input.parquet"))

# theo_UV_input %>% 
#   group_by(state, theo_scen) %>% 
#   reframe(max = max(value),
#           min = min(value),
#           mean = mean(value))
# 
# ggplot(theo_UV_input, aes(x = yday, y = value, color = theo_scen)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))

theo_Ni_input <- tar_read(theo_Ni_input) %>% 
  mutate(state = as.factor(state),
         theo_scen = as.factor(theo_scen)) %>% 
  write_parquet(file.path(out_path, "theo_Ni_input.parquet"))
# 
# theo_Ni_input %>% 
#   group_by(state, theo_scen) %>% 
#   reframe(max = max(value),
#           min = min(value),
#           mean = mean(value))
# 
# ggplot(theo_Ni_input, aes(x = yday, y = value, color = theo_scen)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))

theo_Am_input <- tar_read(theo_Am_input) %>% 
  mutate(state = as.factor(state),
         theo_scen = as.factor(theo_scen)) %>% 
  write_parquet(file.path(out_path, "theo_Am_input.parquet"))

# theo_Am_input %>% 
#   group_by(state, theo_scen) %>% 
#   reframe(max = max(value),
#           min = min(value),
#           mean = mean(value))
# 
# ggplot(theo_Am_input, aes(x = yday, y = value, color = theo_scen)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))

# Scenario growth -------------------------------------------------------------------------------------------------
# theo_scen_grow <- tar_read(theo_scen_grow_arma)
# theo_scen_grow <- theo_scen_grow %>% 
#   distinct(start_date, date, species, state, scen, output, value)
# measures <- levels(theo_scen_grow$output)
# 
# for (lev in 1:length(measures)) {
#   df <- theo_scen_grow[theo_scen_grow$output == measures[lev], ]
#   write_parquet(df, file.path(out_path, str_c("theo_scen_grow_arma_", measures[lev], ".parquet")))
# }
# 
# theo_scen_grow <- tar_read(theo_scen_grow_taxi)
# theo_scen_grow <- theo_scen_grow %>% 
#   distinct(start_date, date, species, state, scen, output, value)
# 
# for (lev in 1:length(measures)) {
#   df <- theo_scen_grow[theo_scen_grow$output == measures[lev], ]
#   write_parquet(df, file.path(out_path, str_c("theo_scen_grow_taxi_", measures[lev], ".parquet")))
# }
# 
# theo_scen_grow <- tar_read(theo_scen_grow_eckl)
# theo_scen_grow <- theo_scen_grow %>% 
#   distinct(start_date, date, species, state, scen, output, value)
# 
# for (lev in 1:length(measures)) {
#   df <- theo_scen_grow[theo_scen_grow$output == measures[lev], ]
#   write_parquet(df, file.path(out_path, str_c("theo_scen_grow_eckl_", measures[lev], ".parquet")))
# }
# 
# theo_scen_grow <- tar_read(theo_scen_grow_ulva)
# theo_scen_grow <- theo_scen_grow %>% 
#   distinct(start_date, date, species, state, scen, output, value)
# 
# for (lev in 1:length(measures)) {
#   df <- theo_scen_grow[theo_scen_grow$output == measures[lev], ]
#   write_parquet(df, file.path(out_path, str_c("theo_scen_grow_ulva_", measures[lev], ".parquet")))
# }
# 
# rm(theo_scen_grow, df)
# gc()
# 
# # Growth rate -----------------------------------------------------------------------------------------------------
# files <- list.files(out_path, full.names = T) %>% 
#   str_subset("growth_rate") %>% 
#   str_subset("ALL", negate = T)
# fs <- list()
# for (i in 1:length(files)) {
#   df <- read_parquet(files[i]) %>% 
#     mutate(scen = as.factor(scen),
#            state = as.factor(state)) %>% 
#     select(-output)
#   fs[[i]] <- df
# }
# fs <- bind_rows(fs) %>% 
#   mutate(species = as.factor(species))
# levels(fs$species)
# 
# fs <- fs %>% 
#   mutate(start_date = as.Date(start_date), 
#          date = as.Date(date)) %>% 
#   write_parquet(file.path(out_path, str_c("total_growth_rate_ALL.parquet")))
# 
# # df <- fs %>% 
# #   filter(species == "A. armata" & state == "SAU")
# 
# ggplot(filter(df, state == "WAS"), aes(x = as.Date(date), y = value, colour = as.factor(scen))) +
#   geom_line(linewidth = 0.75) +
#   facet_wrap(facets = vars(start_date), ncol = 2) +
#   theme_classic()
# 
# # N removed -------------------------------------------------------------------------------------------------------
# files_ni <- list.files(out_path, full.names = T) %>% str_subset("conc_nit")
# files_am <- list.files(out_path, full.names = T) %>% str_subset("conc_amm")
# files_Nf <- list.files(out_path, full.names = T) %>% str_subset("Nf") %>% str_subset("Ns", negate = T) %>% str_subset("loss", negate = T)
# files_Ns <- list.files(out_path, full.names = T) %>% str_subset("Ns") %>% str_subset("Nf", negate = T) %>% str_subset("loss", negate = T)
# files_hm <- list.files(out_path, full.names = T) %>% str_subset("hm")
# 
# fs <- list()
# for (i in 1:length(files_ni)) {
#   N1 <- read_parquet(files_ni[i]) %>%
#     rename(nitrate = value) %>% 
#     select(-c(output))
#   
#   N2 <- read_parquet(files_am[i]) %>% 
#     rename(ammonium = value) %>% 
#     select(-c(output))
#   
#   N <- merge(N1, N2, by = c("date", "species", "state", "scen", "start_date")) %>% 
#     mutate(amb_total = (nitrate + ammonium) * 5, # canopy depth
#            amb_nit = nitrate * 5,
#            amb_amm = ammonium * 5) %>% 
#     pivot_longer(cols = c(amb_total, amb_nit, amb_amm), names_to = "pool", values_to = "N_mgm2") %>% 
#     select(-c(nitrate, ammonium)) #%>% as.data.frame()
#   
#   
#   N1 <- read_parquet(files_Nf[i]) %>%
#     rename(Nf = value) %>% 
#     select(-c(output))
# 
#   N2 <- read_parquet(files_Ns[i]) %>%
#     rename(Ns = value) %>% 
#     select(-c(output))
#   
#   N3 <- read_parquet(files_hm[i]) %>%
#     rename(hm = value) %>% 
#     select(-c(output))
# 
#   Nm <- merge(N1, N2, by = c("date", "species", "state", "scen", "start_date")) #%>% as_tibble()
#   Nm <- merge(Nm, N3, by = c("date", "species", "state", "scen", "start_date"))  %>% #as_tibble() %>% 
#     mutate(pool = "macroalgae",
#            N_mgm2 = (Nf + Ns) * 5) %>% 
#     select(-c(Nf, Ns, hm)) %>% 
#     as_tibble()
#   
#   fs[[i]] <- rbind(N, Nm) %>% 
#     mutate(pool = as.factor(pool))
# }
# rm(N, Nm)
# 
# fs <- bind_rows(fs) %>% 
#   mutate(start_date = as.Date(start_date),
#          date = as.Date(date),
#          species = as.factor(species),
#          scen = as.factor(scen)) 
# 
# fs <- fs %>% 
#   write_parquet(file.path(out_path, str_c("total_N_removed_ALL.parquet")))

# ggplot(fs, aes(x = date, y = N_mgm2, colour = as.factor(species), linetype = pool)) +
#   geom_line() +
#   facet_wrap(facets = vars(state))



