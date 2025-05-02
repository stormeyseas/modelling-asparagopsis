library(GGally)
library(tidyverse)
library(ggplot2)

scens <- data.frame(
  name = c("Shallow bay", "Estuary mouth", "Offshore", "Near salmon"), 
  ammonium = c(0,0,-1,1),
  nitrate = c(0,1,-1,0),
  exposure_depth = c(-1,0,1,0),
  salinity = c(1,-1,0,0),
  light = c(0,-1,1,0)
  # temperature = c(1,0,0,0)
) %>% 
  pivot_longer(names_to = "factor", values_to = "level", cols = c(ammonium, nitrate, exposure_depth, salinity, light))

ggplot(scens, aes(x = factor, y = level, color = name)) +
  geom_point(size = 3, position = position_dodge(0.25)) +
  scale_y_continuous(breaks = c(-1,0,1), labels = c("Low", "Medium", "High")) +
  theme_classic() +
  coord_flip()

(
  p <- ggparcoord(scens,
                  columns = 2:6,
                  scale = "globalminmax",
                  groupColumn = "name") +
    geom_point(size = 1) +
    theme_classic() +
    scale_y_continuous(breaks = c(-1,0,1), labels = c("Low", "Moderate", "High")) +
    scale_x_discrete(breaks = c("ammonium", "nitrate", "exposure_depth", "salinity", "light"#, "temperature"
                                ), 
                     labels = c("Ammonium", "Nitrate", "Exposure/Depth", "Salinity", "Light"#, "Temperature"
                                )) +
    theme(element_text(family = "sans", size = 14, colour = "black"),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none"
          ) +
    theme(aspect.ratio = 0.65)
)

ggsave(plot = p, width = 6, height = 4,
       filename = file.path("C:", "Users", "treimer", "Documents", "R-temp-files", "FRDC-seaweed", "model-outputs", "parcoordplot.png"))


