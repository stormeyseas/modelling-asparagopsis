library(magrittr)
library(dplyr)
library(units)
library(gganimate)
library(ggplot2)
library(gifski)
library(ozmaps)
library(arrow)
library(targets)
library(stringr)
library(here)

runs_data <- here() %>% file.path("data", "processed-model-running")
out_path <- here() %>% file.path("data", "gifs")
if (!dir.exists(out_path)) {dir.create(out_path)}

sc <- "R_scripts/09_model_running.R"
projects <- tibble::tribble(
  ~state, ~store, ~state_long,
  "QLD",  "targets_outputs/_model_running_QLD", "Queensland", 
  "TAS",  "targets_outputs/_model_running_TAS",  "Tasmania", 
  "SAU",  "targets_outputs/_model_running_SAU",  "South Australia", 
  "WAS",  "targets_outputs/_model_running_WAS",  "Western Australia (S)",
  "NSW",  "targets_outputs/_model_running_NSW",  "New South Wales", 
  "VIC",  "targets_outputs/_model_running_VIC",  "Victoria", 
  "WAN",  "targets_outputs/_model_running_WAN",  "Western Australia (N)", 
  "NTE",  "targets_outputs/_model_running_NTE", "Northern Territory"
)

states_bbox <- tar_read(states_bbox, store = projects$st[1])
states <- tar_read(states, store = projects$st[1])

remove_unit("molN")
install_unit("molN", "14.0067 g")

# I want to create gifs of the factors limiting growth in the model for each species and state.
# Tlim ------------------------------------------------------------------------------------------------------
T_clims <- c(0, 1)
T_cbreaks <- seq(T_clims[1], T_clims[2], 0.2)

Tlim_arma <- Tlim_cell[Tlim_cell$species == "A. armata", ]
Tlim_arma <- merge(Tlim_arma, cell_coords, by = c("cell_ID")) %>% 
  rename(state = state.y) %>% 
  select(-state.x)

for (st in c("TAS")){
df <- df %>% 
  filter(state == st) %>% 
  mutate(date = parse_date_time(x = yday, orders = "j"))

pgif <- ggplot(data = df, aes(x = lon, y = lat, fill = T_lim, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", 
                       limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature limitation")) +
  coord_sf(xlim = c(states_bbox[[st]]["lonmin"], states_bbox[[st]]["lonmax"]),
           ylim = c(states_bbox[[st]]["latmin"], states_bbox[[st]]["latmax"])) +
  scale_x_continuous(breaks = seq(100, 160, 2.5)) +
  scale_y_continuous(breaks = seq(-5, -45, -2.5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1] / 6, "inches"),
        legend.key.height = unit(0.2, "inches")) +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(title = "Temperature limitation at DOY: {closest_state}")

# Animate and save
anim <- gganimate::animate(pgif, fps = 4, width = 14.5, height = 14, units = "cm", res = 300, nframes = length(unique(df$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_lim-arma-", st, ".gif")))
}

# Temperature -----------------------------------------------------------------------------------------------------
T_input_cells <- file.path(cell_data, "T_input_cell_av.parquet") %>% 
  read_parquet() %>% 
  # filter(yday %% 2 == 0) %>% # To only show every few time steps (for speed)
  mutate(state_short = factor(state, levels = states_abbr_2),
         state_long = factor(state, levels = states_abbr_2, labels = states_long))

# T_clims <- c(min(T_input_cells$value), max(T_input_cells$value))
T_clims <- c(0, 40)
T_cbreaks <- seq(T_clims[1], T_clims[2], 5)

## Australia plot -------------------------------------------------------------------------------------------------
pgif <- ggplot(data = T_input_cells, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", 
                       limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  coord_sf(xlim = c(bounds[["AUS"]]["lonmin"], bounds[["AUS"]]["lonmax"]),
           ylim = c(bounds[["AUS"]]["latmin"], bounds[["AUS"]]["latmax"])) +
  scale_x_continuous(breaks = seq(100, 160, 5)) +
  scale_y_continuous(breaks = seq(-5, -45, -5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1] / 6, "inches"),
        legend.key.height = unit(0.2, "inches")) +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(title = "Surface temperature at DOY: {closest_state}")

# Animate and save
anim <- gganimate::animate(pgif, fps = 4, width = 14.5, height = 14, units = "cm", res = 300,
                           nframes = length(unique(T_input_cells$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-AUS.gif")))

# Salinity --------------------------------------------------------------------------------------------------------
S_input_cells <- file.path(cell_data, "S_input_cell_av.parquet") %>% 
  read_parquet() %>% 
  # filter(yday %% 2 == 0) %>% # To only show every few time steps (for speed)
  mutate(state_short = factor(state, levels = states_abbr_2),
         state_long = factor(state, levels = states_abbr_2, labels = states_long))

# S_clims <- c(min(S_input_cells$value), max(S_input_cells$value))
S_clims <- c(30, 40)
S_cbreaks <- seq(S_clims[1], S_clims[2], 2)

## Australia plot -------------------------------------------------------------------------------------------------
pgif <- ggplot(data = S_input_cells, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", 
                       limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  coord_sf(xlim = c(bounds[["AUS"]]["lonmin"], bounds[["AUS"]]["lonmax"]),
           ylim = c(bounds[["AUS"]]["latmin"], bounds[["AUS"]]["latmax"])) +
  scale_x_continuous(breaks = seq(100, 160, 5)) +
  scale_y_continuous(breaks = seq(-5, -45, -5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1] / 6, "inches"),
        legend.key.height = unit(0.2, "inches")) +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(title = "Salinity at DOY: {closest_state}")

# Animate and save
anim <- gganimate::animate(pgif, fps = 4, width = 14.5, height = 14, units = "cm", res = 300,
                           nframes = length(unique(S_input_cells$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-AUS.gif")))

# Water velocity --------------------------------------------------------------------------------------------------
U_input_cells <- file.path(cell_data, "U_input_cell_av.parquet") %>% 
  read_parquet() %>% 
  # filter(yday %% 2 == 0) %>% # To only show every few time steps (for speed)
  mutate(state_short = factor(state, levels = states_abbr_2),
         state_long = factor(state, levels = states_abbr_2, labels = states_long))

# U_clims <- c(min(U_input_cells$value), max(U_input_cells$value))
U_clims <- c(0, 1.4)
U_cbreaks <- seq(U_clims[1], U_clims[2], 0.2)

## Australia plot -------------------------------------------------------------------------------------------------
pgif <- ggplot(data = U_input_cells, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", 
                       limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  coord_sf(xlim = c(bounds[["AUS"]]["lonmin"], bounds[["AUS"]]["lonmax"]),
           ylim = c(bounds[["AUS"]]["latmin"], bounds[["AUS"]]["latmax"])) +
  scale_x_continuous(breaks = seq(100, 160, 5)) +
  scale_y_continuous(breaks = seq(-5, -45, -5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1] / 6, "inches"),
        legend.key.height = unit(0.2, "inches")) +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(title = "Water velocity at DOY: {closest_state}")

# Animate and save
anim <- gganimate::animate(pgif, fps = 4, width = 14.5, height = 14, units = "cm", res = 300,
                           nframes = length(unique(U_input_cells$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-AUS.gif")))

# Nitrate concentration -------------------------------------------------------------------------------------------
Ni_input_cells <- file.path(cell_data, "Ni_input_cell_av.parquet") %>% 
  read_parquet() %>% 
  mutate(value = set_units(value, "mg m-3"),
         value = set_units(value, "umol L-1"),
         value = drop_units(value)) %>% 
  # filter(yday %% 2 == 0) %>% # To only show every few time steps (for speed)
  mutate(state_short = factor(state, levels = states_abbr_2),
         state_long = factor(state, levels = states_abbr_2, labels = states_long))

# Ni_clims <- c(min(Ni_input_cells$value), max(Ni_input_cells$value))
Ni_clims <- c(0, 3.5)
Ni_cbreaks <- seq(Ni_clims[1], Ni_clims[2], 0.5)

## Australia plot -------------------------------------------------------------------------------------------------
pgif <- ggplot(data = Ni_input_cells, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", 
                       limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate concentration ("*mu*"M)")) +
  coord_sf(xlim = c(bounds[["AUS"]]["lonmin"], bounds[["AUS"]]["lonmax"]),
           ylim = c(bounds[["AUS"]]["latmin"], bounds[["AUS"]]["latmax"])) +
  scale_x_continuous(breaks = seq(100, 160, 5)) +
  scale_y_continuous(breaks = seq(-5, -45, -5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1] / 6, "inches"),
        legend.key.height = unit(0.2, "inches")) +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(title = "Nitrate concentration at DOY: {closest_state}")

# Animate and save
anim <- gganimate::animate(pgif, fps = 4, width = 14.5, height = 14, units = "cm", res = 300,
                           nframes = length(unique(Ni_input_cells$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-AUS.gif")))

# Ammonium concentration ------------------------------------------------------------------------------------------
Am_input_cells <- file.path(cell_data, "Am_input_cell_av.parquet") %>% 
  read_parquet() %>% 
  mutate(value = set_units(value, "mg m-3"),
         value = set_units(value, "umol L-1"),
         value = drop_units(value)) %>% 
  # filter(yday %% 2 == 0) %>% # To only show every few time steps (for speed)
  mutate(state_short = factor(state, levels = states_abbr_2),
         state_long = factor(state, levels = states_abbr_2, labels = states_long))

# Am_clims <- c(min(Am_input_cells$value), max(Am_input_cells$value))
Am_clims <- c(0, 0.7)
Am_cbreaks <- seq(Am_clims[1], Ni_clims[2], 0.1)

## Australia plot -------------------------------------------------------------------------------------------------
pgif <- ggplot(data = Am_input_cells, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", 
                       limits = Am_clims, 
                       breaks = seq(Am_clims[1], Am_clims[2], 0.1), 
                       name = expression("Nitrate concentration ("*mu*"M)")) +
  coord_sf(xlim = c(bounds[["AUS"]]["lonmin"], bounds[["AUS"]]["lonmax"]),
           ylim = c(bounds[["AUS"]]["latmin"], bounds[["AUS"]]["latmax"])) +
  scale_x_continuous(breaks = seq(100, 160, 5)) +
  scale_y_continuous(breaks = seq(-5, -45, -5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1] / 6, "inches"),
        legend.key.height = unit(0.2, "inches")) +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(title = "Nitrate concentration at DOY: {closest_state}")

# Animate and save
anim <- gganimate::animate(pgif, fps = 4, width = 14.5, height = 14, units = "cm", res = 300,
                           nframes = length(unique(Am_input_cells$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-AUS.gif")))

# State plots -----------------------------------------------------------------------------------------------------
state_theme <- theme_classic() +
  theme(text = element_text(family = "sans", size = 7, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title.position = "bottom", 
        legend.key.width = unit(dev.size()[1]/12, "inches"),
        legend.key.height = unit(0.2, "inches"))

## Tasmania -------------------------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "TAS")
T_input_cells_st <- T_input_cells %>% filter(state == "TAS")
S_input_cells_st <- S_input_cells %>% filter(state == "TAS")
U_input_cells_st <- U_input_cells %>% filter(state == "TAS")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "TAS")
Am_input_cells_st <- Am_input_cells %>% filter(state == "TAS")
xlims <- c(143, 149.25)
ylims <- c(-44.25, -39.05)
xbreaks <- seq(143, 149, 1)
ybreaks <- seq(-45, -38, 1)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-TAS.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-TAS.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-TAS.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-TAS.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-TAS.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-TAS.gif")))

## Victoria -------------------------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "VIC")
T_input_cells_st <- T_input_cells %>% filter(state == "VIC")
S_input_cells_st <- S_input_cells %>% filter(state == "VIC")
U_input_cells_st <- U_input_cells %>% filter(state == "VIC")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "VIC")
Am_input_cells_st <- Am_input_cells %>% filter(state == "VIC")
xlims <- c(141, 150.5)
ylims <- c(-39.75, -37.25)
xbreaks <- seq(xlims[1], xlims[2], 1)
ybreaks <- seq(ylims[1], ylims[2], 0.5)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-VIC.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-VIC.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-VIC.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-VIC.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-VIC.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-VIC.gif")))

## South Australia ------------------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "SAU")
T_input_cells_st <- T_input_cells %>% filter(state == "SAU")
S_input_cells_st <- S_input_cells %>% filter(state == "SAU")
U_input_cells_st <- U_input_cells %>% filter(state == "SAU")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "SAU")
Am_input_cells_st <- Am_input_cells %>% filter(state == "SAU")
xlims <- c(128.5, 141.25)
ylims <- c(-39, -31)
xbreaks <- seq(128.5, 141.5, 1)
ybreaks <- seq(-39, -31, 1)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  transition_states(yday, transition_length = 1, state_length = 30) +
  ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-SAU.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-SAU.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-SAU.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-SAU.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-SAU.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-SAU.gif")))

## Western Australia (North) --------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "WAN")
T_input_cells_st <- T_input_cells %>% filter(state == "WAN")
S_input_cells_st <- S_input_cells %>% filter(state == "WAN")
U_input_cells_st <- U_input_cells %>% filter(state == "WAN")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "WAN")
Am_input_cells_st <- Am_input_cells %>% filter(state == "WAN")
xlims <- c(112, 128)
ylims <- c(-26, -13)
xbreaks <- seq(112, 128, 2)
ybreaks <- seq(-26, -13, 2)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  # transition_states(yday, transition_length = 1, state_length = 30) +
  # ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")
pgif

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-WAN.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-WAN.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-WAN.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-WAN.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-WAN.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-WAN.gif")))

## Western Australia (South) --------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "WAS")
T_input_cells_st <- T_input_cells %>% filter(state == "WAS")
S_input_cells_st <- S_input_cells %>% filter(state == "WAS")
U_input_cells_st <- U_input_cells %>% filter(state == "WAS")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "WAS")
Am_input_cells_st <- Am_input_cells %>% filter(state == "WAS")
xlims <- c(112, 130)
ylims <- c(-36, -24)
xbreaks <- seq(112, 130, 2)
ybreaks <- seq(-36, -24, 2)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  # transition_states(yday, transition_length = 1, state_length = 30) +
  # ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")
pgif

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-WAS.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-WAS.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-WAS.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-WAS.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-WAS.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-WAS.gif")))

## Queensland -----------------------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "QLD")
T_input_cells_st <- T_input_cells %>% filter(state == "QLD")
S_input_cells_st <- S_input_cells %>% filter(state == "QLD")
U_input_cells_st <- U_input_cells %>% filter(state == "QLD")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "QLD")
Am_input_cells_st <- Am_input_cells %>% filter(state == "QLD")
xlims <- c(112, 130)
ylims <- c(-36, -24)
xbreaks <- seq(112, 130, 2)
ybreaks <- seq(-36, -24, 2)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  # transition_states(yday, transition_length = 1, state_length = 30) +
  # ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")
pgif

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-QLD.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-QLD.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-QLD.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-QLD.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-QLD.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-QLD.gif")))

## Northern Territory ---------------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "NTE")
T_input_cells_st <- T_input_cells %>% filter(state == "NTE")
S_input_cells_st <- S_input_cells %>% filter(state == "NTE")
U_input_cells_st <- U_input_cells %>% filter(state == "NTE")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "NTE")
Am_input_cells_st <- Am_input_cells %>% filter(state == "NTE")
xlims <- c(112, 130)
ylims <- c(-36, -24)
xbreaks <- seq(112, 130, 2)
ybreaks <- seq(-36, -24, 2)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  # transition_states(yday, transition_length = 1, state_length = 30) +
  # ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")
pgif

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-NTE.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-NTE.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-NTE.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-NTE.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-NTE.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-NTE.gif")))

## New South Wales ------------------------------------------------------------------------------------------------
I_input_cells_st <- I_input_cells %>% filter(state == "NSW")
T_input_cells_st <- T_input_cells %>% filter(state == "NSW")
S_input_cells_st <- S_input_cells %>% filter(state == "NSW")
U_input_cells_st <- U_input_cells %>% filter(state == "NSW")
Ni_input_cells_st <- Ni_input_cells %>% filter(state == "NSW")
Am_input_cells_st <- Am_input_cells %>% filter(state == "NSW")
xlims <- c(112, 130)
ylims <- c(-36, -24)
xbreaks <- seq(112, 130, 2)
ybreaks <- seq(-36, -24, 2)

pgif <- ggplot(I_input_cells_st, aes(x = lon, y = lat, fill = value, frame = yday)) +
  geom_raster() +
  geom_sf(data = ozmap_data(data = "states"), inherit.aes = F) +
  scale_fill_viridis_c(guide = "colorbar", limits = I_clims, breaks = I_cbreaks, 
                       name = expression("Irradiance (photons m"^-2*" s"^-1*")")) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) +
  state_theme +
  # transition_states(yday, transition_length = 1, state_length = 30) +
  # ease_aes('linear') +   # Smooth transitions
  labs(x = "Longitude", y = "Latitude", title = "Surface irradiance at DOY: {closest_state}")
pgif

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(I_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("I_input-NSW.gif")))

pgif <- pgif %+% T_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = T_clims, breaks = T_cbreaks, 
                       name = expression("Temperature ("*degree*"C)")) +
  labs(x = "Longitude", y = "Latitude", title = "Surface temperature at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(T_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("T_input-NSW.gif")))

pgif <- pgif %+% S_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = S_clims, breaks = S_cbreaks, 
                       name = expression("Salinity (g L"^-1*")")) +
  labs(title = "Salinity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(S_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("S_input-NSW.gif")))

pgif <- pgif %+% U_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = U_clims, breaks = U_cbreaks, 
                       name = expression("Water velocity (m s"^-1*")")) +
  labs(x = "Longitude", y = "Latitude", title = "Water velocity at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(U_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("U_input-NSW.gif")))

pgif <- pgif %+% Ni_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Ni_clims, breaks = Ni_cbreaks, 
                       name = expression("Nitrate ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Nitrate concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Ni_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Ni_input-NSW.gif")))

pgif <- pgif %+% Am_input_cells_st +
  scale_fill_viridis_c(guide = "colorbar", limits = Am_clims, breaks = Am_cbreaks, 
                       name = expression("Ammonium ("*mu*"M)")) +
  labs(x = "Longitude", y = "Latitude", title = "Ammonium concentration at DOY: {closest_state}")

anim <- gganimate::animate(pgif, fps = 4, width = 10, height = 9.5, units = "cm", res = 300, nframes = length(unique(Am_input_cells_st$yday)))
anim_save(animation = anim, filename = file.path(out_path, str_c("Am_input-NSW.gif")))
