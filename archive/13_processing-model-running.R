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
out_path <- file.path(base_path, "data", "processed-model-running")

Sys.setenv(TAR_PROJECT = "project_model_running")

## Total growth all vars ------------------------------------------------------------------------------------------
tar_read(total_growth_arma) %>%
  write_parquet(file.path(out_path, "total_growth_arma.parquet"))
tar_read(total_growth_taxi) %>%
  write_parquet(file.path(out_path, "total_growth_taxi.parquet"))
tar_read(total_growth_eckl) %>%
  write_parquet(file.path(out_path, "total_growth_eckl.parquet"))
tar_read(total_growth_ulva) %>%
  write_parquet(file.path(out_path, "total_growth_ulva.parquet"))

# Total growth ----------------------------------------------------------------------------------------------------
file_list <- c(
  file.path(out_path, "total_growth_arma.parquet"),
  file.path(out_path, "total_growth_taxi.parquet"),
  file.path(out_path, "total_growth_eckl.parquet"),
  file.path(out_path, "total_growth_ulva.parquet")
)
file_list2 <- c(
  file.path(out_path, "total_N_removed_arma.parquet"),
  file.path(out_path, "total_N_removed_taxi.parquet"),
  file.path(out_path, "total_N_removed_eckl.parquet"),
  file.path(out_path, "total_N_removed_ulva.parquet")
)
species <- c("A. armata", "A. taxiformis", "Ecklonia", "Ulva")

for (i in 1:4) {
  tg <- read_parquet(file_list[i]) %>% 
    filter(output %in% c("Nf", "Ns")) %>% 
    mutate(output = droplevels(output))
  
  tg %>% 
    group_by(note) %>% 
    reframe(num_cells = n()) %>% 
    write_parquet(file.path(out_path, str_c("run_outcomes_", i, ".parquet"))) %>% 
    select(-note)
  
  tg <- tg %>% 
    group_by(start_date, date, state, cell_ID) %>% 
    reframe(value = sum(value, na.rm = T) * 5)
  
  start_values <- tg %>% 
    filter(date == start_date) %>% 
    distinct(start_date, cell_ID, value) %>% 
    rename(start_value = value)
  
  tg <- merge(tg, start_values, by = c("start_date", "cell_ID")) %>% 
    # mutate(removed = start_value - value) %>% 
    # select(-c(start_value, value)) %>% 
    filter(start_value != 0)

  tg <- tg %>% 
    arrange(-value) %>% 
    group_by(state, cell_ID, start_date) %>% 
    slice_head(n = 1) %>% 
    ungroup() %>% 
    rename(harv_date = date)
  
  tg <- tg %>% 
    mutate(diff = value - start_value,
           species = species[i]) %>% 
    write_parquet(file_list2[i])
  
  # summ <- tg %>% 
  #   group_by(cell_ID, start_date) %>% 
  #   reframe(n = n())
}

# Nf --------------------------------------------------------------------------------------------------------------
# files <- list.files(out_path, full.names = T) %>% 
#   str_subset("Nf.parquet") %>% 
#   str_subset("ulva") %>% 
#   str_subset("Ns", negate = T)
# 
# Nf <- files %>% 
#   read_parquet()
# 
# Nf %>% 
#   filter(state == "SAU") %>% 
# ggplot(aes(x = date, y = value, color = cell_ID)) +
#   geom_line(linewidth = 0.75) +
#   facet_grid(facets = vars(batch)) +
#   theme_classic() 
# rm(Nf)
# gc()

# Growth rate -----------------------------------------------------------------------------------------------------
# files <- list.files(out_path, full.names = T) %>% 
#   str_subset("growth_rate") %>% 
#   str_subset("ALL", negate = T)
# 
# fs <- list()
# for (i in 1:length(files)) {
#   df <- read_parquet(files[i]) %>% 
#     select(-output, -note) %>% 
#     group_by(state, cell_ID, species) %>% 
#     reframe(n = n()) %>% 
#     mutate(batch = cur_group_id()) %>% 
#     ungroup()
#   fs[[i]] 
# }
# fs <- bind_rows(fs) %>% 
#   mutate(note = as.factor(note),
#          species = as.factor(species))
# 
# levels(fs$note)
# fs_err <- fs[!is.na(fs$note), ]
# 
# fs <- fs %>% 
#   group_by(date, state, species) %>% 
#   reframe(num_cells = n(),
#           mean = mean(value, na.rm = T),
#           sd = sd(value, na.rm = T)) %>% 
#   mutate(date = as.Date(date)) %>% 
#   write_parquet(file.path(out_path, str_c("total_growth_rate_ALL.parquet")))
# 
# ggplot(filter(fs, state == "SAU"), 
#        aes(x = date, y = mean, color = species)) +
#   geom_line(linewidth = 0.75) +
#   # facet_wrap(facets = vars(state), ncol = 2) +
#   theme_classic()

# N removed -------------------------------------------------------------------------------------------------------
files <- list.files(out_path, full.names = T)
files_ni <- files %>% str_subset("conc_nit")
files_am <- files %>% str_subset("conc_amm")
files_Nf <- files %>% str_subset("Nf") %>% str_subset("Ns", negate = T) %>% str_subset("loss", negate = T)
files_Ns <- files %>% str_subset("Ns") %>% str_subset("Nf", negate = T) %>% str_subset("loss", negate = T)
files_hm <- files %>% str_subset("hm")

fs <- list()
for (i in 1:length(files_ni)) {
  Nf <- read_parquet(files_Nf[i]) %>% select(-c(note)) %>% mutate(output = droplevels(output))
  Ns <- read_parquet(files_Ns[i]) %>% select(-c(note)) %>% mutate(output = droplevels(output))
  
Nm <- merge(Nf, Ns, by = c("harv_date", "start_date", "state", "cell_ID", "species"))

    bind_rows() %>%
    select(-c(note)) %>% 
    pivot_wider(names_from = output, values_from = harv_value, id_cols = c(start_date, harv_date, state, cell_ID, species)) %>%
    mutate(pool = "macroalgae",
           N_mgm2 = (Nf + Ns) * 5) %>%
    select(-c(Nf, Ns, hm))


  fs[[i]] <- Nm %>% 
    # rbind(N, Nm) %>%
    mutate(pool = as.factor(pool))
}
rm(Nm)

fs <- bind_rows(fs) %>%
  mutate(note = as.factor(note),
         species = as.factor(species)) %>%
  write_parquet(file.path(out_path, str_c("total_N_removed_cells_taxi.parquet")))

fs <- fs %>%
  group_by(start_date, note, state, species, pool) %>%
  reframe(num_cells = n(),
          mean = mean(N_mgm2, na.rm = T),
          sd = sd(N_mgm2, na.rm = T)) %>%
  write_parquet(file.path(out_path, str_c("total_N_removed_states_taxi.parquet")))





