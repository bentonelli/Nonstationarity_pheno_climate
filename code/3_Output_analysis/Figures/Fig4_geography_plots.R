library(dggridR)
library(dplyr)

sp_cl_df <- readRDS("data/output/sp_cl_mdl_params.rds")
snow_p_data <- readRDS("data/mdl_data/ms_data2.rds") %>% select(cl_sp,snow_p)

water_col = "aliceblue"
country_col = "grey80"
midpoint_col = "white"
grid_outline_col = "grey60"

hexgrid6 <- dggridR::dgconstruct(res = 6) 
countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,long < -40)
lakes <- map_data("lakes") %>% filter(long < -40)

median_eff <- sp_cl_df %>% group_by(cell) %>% summarise(eff_med = median(b1s_median),count = n()) %>% filter(count >=10)

grid  <- dgcellstogrid(hexgrid6,median_eff$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,median_eff,by.x="cell",by.y="cell")

p1 <- ggplot() + coord_map("mollweide",xlim=c(-145,-65),ylim=c(25,65)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group),alpha=1,fill=country_col) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill=water_col) +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=eff_med), alpha=1)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color=grid_outline_col) +
  geom_path(data=lakes, aes(x=long, y=lat, group=group),alpha=1,color="black") +
  geom_path(data=countries, aes(x=long, y=lat, group=group),alpha=1,color="black") +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), alpha=.4,fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(3.25, 'cm'),
        legend.title=element_text(size=32), 
        legend.text=element_text(size=32),
        legend.position = c(.15,.4),
        legend.background = element_rect(fill=alpha(water_col,))) +
  theme(panel.background = element_rect(fill = alpha(water_col))) +
  scale_fill_gradient2(low="firebrick4", mid=midpoint_col,high="dodgerblue4",
                       name = "",midpoint=0,lim=c(-2.5,0.5),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
p1
ggsave("temp_geo.png", width = 20.58,height=11.29, device='png', dpi=500)

sp_eff_wo_no_snow <- merge(sp_cl_df,snow_p_data,by="cl_sp")

median_eff <- sp_eff_wo_no_snow %>% filter(snow_p == 1) %>% group_by(cell) %>% summarise(eff_med = median(b2s_median),count = n()) %>% filter(count >=10)

grid  <- dgcellstogrid(hexgrid6,median_eff$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,median_eff,by.x="cell",by.y="cell")

p2 <- ggplot() + coord_map("mollweide",xlim=c(-145,-65),ylim=c(25,65)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group),alpha=1,fill=country_col) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill=water_col) +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=eff_med), alpha=1)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color=grid_outline_col) +
  geom_path(data=lakes, aes(x=long, y=lat, group=group),alpha=1,color="black") +
  geom_path(data=countries, aes(x=long, y=lat, group=group),alpha=1,color="black") +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), alpha=.4,fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(3.25, 'cm'),
        legend.title=element_text(size=32), 
        legend.text=element_text(size=32),
        legend.position = c(.15,.4),
        legend.background = element_rect(fill=water_col)) +
  theme(panel.background = element_rect(fill = alpha(water_col))) +
  scale_fill_gradient2(low="wheat3", mid=midpoint_col,high="skyblue4",name = "",
                       midpoint=0,lim=c(-0.25,1.5),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
p2
ggsave("snow_geo.png", width = 20.58,height=11.29, device='png', dpi=500)


median_eff <- sp_cl_df %>% group_by(cell) %>% summarise(eff_med = median(b3s_median),count = n()) %>% filter(count >=10)

grid  <- dgcellstogrid(hexgrid6,median_eff$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,median_eff,by.x="cell",by.y="cell")
p3 <- ggplot() + coord_map("mollweide",xlim=c(-145,-65),ylim=c(25,65)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group),alpha=1,fill=country_col) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill=water_col) +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=eff_med), alpha=1)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color=grid_outline_col) +
  geom_path(data=lakes, aes(x=long, y=lat, group=group),alpha=1,color="black") +
  geom_path(data=countries, aes(x=long, y=lat, group=group),alpha=1,color="black") +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), alpha=.4,fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(3.25, 'cm'),
        legend.title=element_text(size=32), 
        legend.text=element_text(size=32),
        legend.position = c(.15,.4),
        legend.background = element_rect(fill=water_col)) +
  theme(panel.background = element_rect(fill = alpha(water_col))) +
  scale_fill_gradient2(low="wheat3", mid=midpoint_col,high="darkgreen",
                       name = "",midpoint=0,lim=c(-0.35,0.35),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
#scale_color_gradient2(low="tan4", mid="grey",high="forestgreen",midpoint=0,lim=c(min(cellcenters$eff),max(cellcenters$eff)))
p3
ggsave("precip_geo.png", width = 20.58,height=11.29, device='png', dpi=500)
