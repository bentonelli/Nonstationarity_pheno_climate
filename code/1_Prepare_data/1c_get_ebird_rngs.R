library(ebirdst)
library(sf)
library(dggridR)
library(readr)

spec_list <- readRDS("data/mdl_data/ms_data2.rds")
spec_list <- spec_list$species
for (trg_sp in spec_list){
  #trg_sp <- "Clay-colored Sparrow"
  print(trg_sp)
  #path <- ebirdst_download(species = trg_sp,force=TRUE,tifs_only = TRUE
  if(trg_sp %in% ebirdst_runs$common_name){
    path <- ebirdst_download(species = trg_sp,force=TRUE,tifs_only = TRUE,pattern="range_")
    ranges <- load_ranges(path, resolution = "lr")
    
    #Plot breeding range
    pdf(paste("data/spec_breed_plots/",trg_sp,".pdf",sep=""))
    plot(ranges$geom[1])
    
    #Sample ~10000 draws in the range. This is optimized for speed 
    points_in_range <- st_sample(ranges$geom[1],10000,by_polygon=TRUE,exact = FALSE)
    plot(points_in_range,cex=.1)
    
    as.data.frame(points_in_range)
    hexgrid6 <- dggridR::dgconstruct(res = 6) 
    
    points_as_df <- as.data.frame(st_coordinates(points_in_range))
    
    colnames(points_as_df) <- c("lon","lat")
    
    points_as_df$cell <- dgGEO_to_SEQNUM(hexgrid6, points_as_df$lon,points_as_df$lat)$seqnum
    
    grid  <- dgcellstogrid(hexgrid6,points_as_df$cell,frame=TRUE,wrapcells=TRUE)
    grid  <- merge(grid,points_as_df,by.x="cell",by.y="cell")
    
    countries <- map_data("world")
    states <- map_data("state")
    countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                                "Belize","Honduras","El Salvador","Nicaragua",
                                                "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                                "Dominican Republic","Puerto Rico","Colombia",
                                                "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
    
    p <- ggplot() + coord_map("mollweide",xlim=c(-150,-55),ylim=c(20,70)) + ggtitle(trg_sp) +
      geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),axis.text.y = element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "null"),
            panel.spacing = unit(c(0, 0, 0, 0), "null")) +
      theme(panel.background = element_rect(fill = alpha("white",.5))) +
      geom_polygon(data=grid,      aes(x=long, y=lat.x, group=group), alpha=0.8)  +
      geom_path   (data=grid,      aes(x=long, y=lat.x, group=group), alpha=0.4, color="white")# +
    #geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=sp_yr_data$count)) +
    #scale_fill_gradient2(low="dodgerblue4", mid="grey",high="firebrick4",midpoint=0)
    print(p)
    dev.off()
    saveRDS(unique(points_as_df$cell),file=paste("data/spec_breed_cells/",trg_sp,".rds",sep=""))
  } else {
    print(paste(trg_sp,"is not modeled by ebirdst"))
  }
}
