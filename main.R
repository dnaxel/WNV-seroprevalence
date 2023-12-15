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

### Get poligons of all municipalities
italia_comuni=getData("GADM",country="ITA",level=3)
emilia_con_comuni=italia_comuni[which(italia_comuni$NAME_1=="Emilia-Romagna"),]
bird.map.data <- uccelli[uccelli$ANNO!=2012,]
emilia_con_comuni$NAME_3 <- toupper(emilia_con_comuni$NAME_3)
#to get ID: emilia_con_comuni@polygons[[1]]@ID
marche_con_comuni <- italia_comuni[which(italia_comuni$NAME_1=="Marche"),]
nuovi_comuni <- c("Casteldelci", "Maiolo", "Pennabilli", "Novafeltria", "San Leo", "Sant' Agata Feltria", "Talamello")

### Get Map
#an API key to download Stadia maps can be obtained at https://client.stadiamaps.com/signup/
register_stadiamaps(key = "your-API-key")
map_emilia <- get_map(location = getbb("emilia-romagna"), maptype = "stamen_terrain_lines", source = "stadia")
