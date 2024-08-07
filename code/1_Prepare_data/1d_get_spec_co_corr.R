#Script to get the impact of temp, prcp on wintering areas of species of interest

library(raster)
library(ncdf4)
library(RColorBrewer)
library(ebirdst)
library(sf)

#Read in species names
spec_names <- readRDS("data/mdl_data/ms_data2.rds")
spec_names <- unique(spec_names$species)

list_vars <- c("ENSO_Temperature","ENSO_Precipitation",
               "NAO_Temperature","NAO_Precipitation")

all_co_corr <- nc_open("co_corr.nc")
my.palette <- brewer.pal(n = 11, name = "BrBG")

sp_co_rec <- c()
pdf("wb_wintering_range_co.pdf")
for (each_spec in spec_names){
  print(each_spec)
  # Species of interest
  #each_spec <- "Purple Martin"
  path <- ebirdst_download(species = each_spec,pattern = "range",force=FALSE,tifs_only = TRUE)
  ranges <- load_ranges(path, resolution = "lr")
  nb_ind <- which(ranges$season == "nonbreeding")

  if (length(nb_ind) > 0){
    #Sample points
    points_in_range <- st_sample(ranges$geom[nb_ind],1000,by_polygon=TRUE,exact = FALSE)
    
    #Get corr values
    points_in_range <- st_coordinates(points_in_range)
    
    #Only keep points in NA and SA
    points_in_range <- points_in_range[which(points_in_range[,1] >= -170 & points_in_range[,1] <= -30),]
    
    for (targ_metric in list_vars){
      # Metric of interest
      #targ_metric <- "NAO_Precipitation"
      tmp_raster <- brick("data/climate_oscillations/co_corr.nc", varname=targ_metric)
      tmp_raster <- rotate(tmp_raster)
      
      # Plot
      plot(tmp_raster,zlim=c(-1,1),ylim=c(-25,80),xlim=c(-150,-20),main=targ_metric,col=my.palette)
      
      #Plot NB range
      plot(ranges$geom[nb_ind],add=TRUE,lwd=1,lty=3)
      
      vals <- raster::extract(tmp_raster$layer,points_in_range)
      
      sp_co_var_abs_median <- median(abs(vals),na.rm=TRUE)
      sp_co_var_median <- median(vals,na.rm=TRUE)
      to_add <- c(each_spec,targ_metric,sp_co_var_median,sp_co_var_abs_median)
      sp_co_rec <- rbind(to_add,sp_co_rec)
      }
    } else {
      print(paste(each_spec,"has no nb data,sorry!"))
  }
}
dev.off()
sp_co_rec <- as.data.frame(sp_co_rec)
colnames(sp_co_rec) <- c("spec","co_var","median_impact","abs_median_impact")
sp_co_rec$median_impact <- as.numeric(sp_co_rec$median_impact)
sp_co_rec$abs_median_impact <- as.numeric(sp_co_rec$abs_median_impact)

saveRDS(sp_co_rec,"data/mdl_data/sp_co_rec.rds")
