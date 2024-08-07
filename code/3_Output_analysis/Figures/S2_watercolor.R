# Plots to show effects of NAO, ENSO across the Americas

#Get at predicted arrival for particular species, using both models

library(ggplot2)
library(MCMCvis)
library(rgl)
library(gridExtra)
library(grid)
library(gridBase)
library(dplyr)
library(dggridR)
library(ebirdst)
library(sf)
library(geosphere)
library(terra)
library(tidyterra)

color1 <- "palevioletred"
color2 <- "mediumturquoise"

colfunc <- colorRampPalette(c(color1, color2))
colfunc(100)
plot(rep(1,100),col=alpha(colfunc(100),1),pch="|",cex=10)

d2 <- readRDS("data/mdl_data/ms_data2.rds")
mdl_out <- readRDS("data/output/sp_mdl_params.rds")

co_corr <- rast("co_corr.nc")


sp_enso_nao_eff <- data.frame(species = unique(d2$species))
sp_enso_nao_eff$mu_phi1 <- mdl_out[[4]]$mean
sp_enso_nao_eff$mu_phi2 <- mdl_out[[5]]$mean

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,long < -40)
lakes <- map_data("lakes") %>% filter(long < -40)

p1 <- ggplot() + coord_map("mercator",xlim=c(-145,-65),ylim=c(-25,65)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5))) 
p1

p2 <- ggplot() + coord_map("mercator",xlim=c(-145,-65),ylim=c(-25,65)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5))) 
p2

count <- 0
for (each_species in sp_enso_nao_eff$species){
  count <- count + 1
  print(each_species)
  path <- ebirdst_download(species = each_species,force=FALSE,tifs_only = TRUE,pattern="range_")
  ranges <- load_ranges(path, resolution = "lr")
  nb_range <- ranges$geom[2]
  if (sp_enso_nao_eff$mu_phi1[count] > 0){
    p1 <- p1 + geom_sf(data = nb_range,fill=color1,
                     alpha=abs(sp_enso_nao_eff$mu_phi1[count])/20,
                     color="transparent")
  } else {
    p1 <- p1 + geom_sf(data = nb_range,fill=color2,
                       alpha=abs(sp_enso_nao_eff$mu_phi1[count])/20,
                       color="transparent")
  }
  
  if (sp_enso_nao_eff$mu_phi2[count] > 0){
    p2 <- p2 + geom_sf(data = nb_range,fill=color1,
                       alpha=abs(sp_enso_nao_eff$mu_phi2[count])/20,
                       color="transparent")
  } else {
    p2 <- p2 + geom_sf(data = nb_range,fill=color2,
                       alpha=abs(sp_enso_nao_eff$mu_phi2[count])/20,
                       color="transparent")
  }
}

p1_cropped <- p1 + coord_sf(crs = st_crs("ESRI:54009"),
              xlim=c(-12040095.7,-3040095.7),ylim=c(-6020047.85,7520047.85))

p2_cropped <- p2 + coord_sf(crs = st_crs("ESRI:54009"),
              xlim=c(-12040095.7,-3040095.7),ylim=c(-6020047.85,7520047.85))


ggsave("p1_cropped.png",plot=p1_cropped)
ggsave("p2_cropped.png",plot=p2_cropped)

enso_temp <- co_corr$ENSO_Temperature
nao_temp <- co_corr$NAO_Temperature

crs(enso_temp) <- "+proj=longlat +datum=WGS84 +no_defs"
enso_temp <- project(enso_temp,"ESRI:54009")

crs(nao_temp) <- "+proj=longlat +datum=WGS84 +no_defs"
nao_temp <- project(nao_temp,"ESRI:54009")

enso_temp_plot <- ggplot() + 
  geom_spatraster(data=enso_temp,interpolate = TRUE,maxcell = 100000000) + 
  coord_sf(xlim=c(-12040095.7,-3040095.7),ylim=c(-6020047.85,7520047.85)) +
  scale_fill_gradient2(low="dodgerblue3", mid="grey",high="firebrick3",
                       midpoint=0,lim=c(-1,1),na.value=NA,name="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = alpha("white",.5))) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(2, 'cm'),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24),
        legend.position = c(.25,.28))
enso_temp_plot
nao_temp_plot <- ggplot() + 
  geom_spatraster(data=nao_temp,interpolate = TRUE,maxcell = 100000000) + 
  coord_sf(xlim=c(-12040095.7,-3040095.7),ylim=c(-6020047.85,7520047.85)) +
  scale_fill_gradient2(low="dodgerblue3", mid="grey",high="firebrick3",
                       midpoint=0,lim=c(-1,1),na.value=NA,name="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = alpha("white",.5))) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(2, 'cm'),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24),
        legend.position = c(.25,.28))
nao_temp_plot

ggsave("data/output/figures/nao_temp_plot.png",plot=nao_temp_plot)
ggsave("data/output/figures/enso_temp_plot.png",plot=enso_temp_plot)
