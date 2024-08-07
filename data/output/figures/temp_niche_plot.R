library(dggridR)
library(dplyr)

sp_cl_df <- readRDS("sp_cl_eff_df_10_26_23.rds")
d2 <- readRDS("ms_data2_7_26_23.rds") %>% select(c("species","cell","mean_tmin_29","mean_prcp_29","mean_swe_29","sd_tmin_c29","sd_prcp_c29","sd_swe_c29"))

hexgrid6 <- dggridR::dgconstruct(res = 6) 

sp_cl_w_cc <- sp_cl_df
colnames(sp_cl_w_cc)[1] <- "species"
sp_cl_w_cc <- merge(sp_cl_w_cc,d2,by=c("species","cell")) %>% unique()

cellcenters <- as.data.frame(dgSEQNUM_to_GEO(hexgrid6,sp_cl_w_cc$cell)) 
sp_cl_w_cc$latitude <- cellcenters$lat_deg
sp_cl_w_cc$longitude <- cellcenters$lon_deg

sp_temp <- sp_cl_w_cc %>% select(mean_tmin_29,sd_tmin_c29,b1s_median,latitude) %>% unique()
sp_temp$days_per_degree <- sp_temp$b1s_median/sp_temp$sd_tmin_c29


colfunc <- viridis::plasma(40,direction=-1)
lat_cols <- round(sp_temp$latitude)-min(round(sp_temp$latitude)) + 1
color_lats_ind <- colfunc[lat_cols]
        
par(las=1)
par(pty="s")
plot(NULL,xlim=c(-20,20),ylim=c(-3,1.5),xlab="Mean Temperature, Arrival", ylab="Arrival Change, Days/ +1 Degree",
     cex.lab=2,cex.axis=2)
abline(h=0,lty=2,col="grey20")
abline(v=0,lty=2,col="grey20")
points(sp_temp$mean_tmin_29,sp_temp$days_per_degree,
       pch=21,cex=.5,bg=alpha(color_lats_ind,.85),col="transparent")
legend_image <- as.raster(matrix(colfunc, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
#text(x=1.5, y = seq(0,1,l=5), labels = seq(min(round(sp_temp$latitude)),max(round(sp_temp$latitude)),by=8)) #Not working right
rasterImage(legend_image, 0, 0, 1,1)

temp_plot <- ggplot(sp_temp,aes(x=mean_tmin_29,y=days_per_degree,fill=latitude)) +
  geom_hline(yintercept = 0,lty=4,col="darkgrey") + 
  geom_point(pch=21,color="transparent",size=1.5) +
  scale_fill_viridis_c("Latitude (Degree N\u00B0)",direction=-1,option="inferno",alpha=.9) +
  guides(fill = guide_colorbar(
    direction="horizontal",
    title.position = "top",
    title.hjust = .5)) +
  theme_classic() + 
  xlab("Mean Minimum Temperature (\u00B0C)") +
  ylab("Sensitivity (days per +1\u00B0C)") +
  theme(aspect.ratio=1,
        text = element_text(size = 20,family = "sans"),
        axis.text = element_text(size = 20,colour = "black"),
        legend.position = c(0.76, 0.125),
        legend.key.size = unit(1.7, 'cm'), #change legend key size
        legend.key.height = unit(1.7, 'cm'), #change legend key height
        legend.key.width = unit(1.7, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14))
temp_plot
ggsave("temp_plot.png", width = 9,height=8, device='png', dpi=500)




