graphics.off()
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(deSolve)
library(mvtnorm)
library(openxlsx)
library(ggpubr)
library(cowplot)
library(gridExtra)

####################
### DATA LOADING ###
####################

################
### Figure 1 ###
################
library(ggmap)
library(osmdata)
library(raster)
library(viridis)
library(ggnewscale)

### Get poligons of all municipalities
italia_comuni=getData("GADM",country="ITA",level=3)
emilia_con_comuni=italia_comuni[which(italia_comuni$NAME_1=="Emilia-Romagna"),]
emilia_con_comuni$NAME_3 <- toupper(emilia_con_comuni$NAME_3)
marche_con_comuni <- italia_comuni[which(italia_comuni$NAME_1=="Marche"),]
nuovi_comuni <- c("Casteldelci", "Maiolo", "Pennabilli", "Novafeltria", "San Leo", "Sant' Agata Feltria", "Talamello")
emilia.data <- rbind(fortify(emilia_con_comuni), fortify(marche_con_comuni[which(marche_con_comuni$NAME_3 %in% nuovi_comuni),]))

# Get bird data clustered by municipality
birds.by.town <- read.xlsx("birds_by_town.xlsx")
emilia.borders <- merge(birds.by.town[,c("id", "n_specimens", "prev", "subregion")], cbind(n = 1:nrow(emilia_data), emilia_data), all=T)
emilia.borders <- emilia.borders[order(emilia.borders$n),]
emilia.borders$n_specimens[is.na(emilia.borders$n_specimens)] <- 0

### Get Map
#an API key to download Stadia maps can be obtained at https://client.stadiamaps.com/signup/
register_stadiamaps(key = "your-API-key")
map.emilia <- get_map(location = getbb("emilia-romagna"), maptype = "stamen_terrain_lines", source = "stadia")

### Bubble Plot
p1 <- ggmap(map_emilia) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = cluster), data = emilia_borders,
               alpha = .4, size = 0.5, color = "grey60") +
  scale_fill_manual(values = colconf[1:3], name = "Subregion") +
  geom_point(data = comuni, aes(x=long, y=lat, size=esemplari_reali, color=prev)) +
  scale_size_continuous(range=c(1, 5), trans = "log10", name = "Sample size", breaks = c(1, 10, 100, 1000, 10000)) +
  scale_color_viridis(name = "WNV prevalence", option = "plasma", trans = "sqrt") +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
ggsave("bubble_plot.png", p1, dpi = 1200, width = 9, height = 6)
