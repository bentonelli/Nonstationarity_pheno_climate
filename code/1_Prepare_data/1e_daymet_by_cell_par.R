#Download daymet by cell
library(daymetr)
library(dggridR)
library(dplyr)
library(sf)
library(terra)
library(zoo)
library(foreach)
library(parallel)

#Get cell data
all_data <- readRDS("data/mdl_data/ms_data2.rds")

#Environmental variables
ev_list <- c("tmin","swe","prcp")

#Get list of cells
cl_list <- unique(all_data$cell)

for (param_targ in ev_list){
  #Set up to run in parallel - do this 3 times, one for each enviro variable
  n.cores <- parallel::detectCores() - 2
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  x <- foreach(
    cl_target = cl_list, 
    .combine = 'c'
  ) %dopar% {
    library(daymetr)
    library(dggridR)
    library(dplyr)
    library(sf)
    library(terra)
    library(zoo)
    library(foreach)
    library(parallel)
    # Construct grid
    hexgrid6 <- dggridR::dgconstruct(res = 6) 
    
    grid  <- dgcellstogrid(hexgrid6,cl_list,frame=TRUE,wrapcells=TRUE)
    
    cl_details <- filter(grid,cell==cl_target)
    #Point to folder
    new_folder_name <- paste(param_targ,"/",cl_target,sep="")
    
    folder_path <- paste("data/enviro_var_daily/temp/",new_folder_name,sep="")
    
    #Create bounding box with max/min long + lat
    bb_cl <- c(max(cl_details$lat),min(cl_details$long),min(cl_details$lat),max(cl_details$long))
    
    #For some reason, the daymetr package doesn't like cells 32 and 58. To fix
    #this, adjust the bounding box to be slightly larger.
    if (cl_target == 32 | cl_target == 58){
      bb_cl[4] <- -157
    }
    
    
    #For each year, download daymet data, process, and add to file to be saved.
    cell_data_all <- c()
    for (targ_year in 2001:2022){
      #targ_year <- 2001
      
      #Download daymet data for bounding box
      download_daymet_ncss(
        location = bb_cl,
        start = targ_year,
        end = targ_year,
        param = param_targ,
        frequency = "daily",
        mosaic = "na",
        path = folder_path,
        silent = FALSE,
        force = FALSE,
        ssl = TRUE
      )
      
      #Read in file
      t1 = terra::rast(paste(folder_path,"/",param_targ,"_daily_",targ_year,"_ncss.nc",sep=""))
      
      #Sample points, then limit to extent of cell
      ss_points <- spatSample(t1,5000,method = "random",replace=TRUE,na.rm=TRUE,as.points=TRUE,xy=TRUE)
      #points(ss_points$x,ss_points$y)
      ev_vals <- as.data.frame(ss_points[,3:367])
      
      lonlats_in <- crds(project(ss_points,"epsg:4326"),df=TRUE)
      lonlats_in$cell <- dgGEO_to_SEQNUM(hexgrid6, lonlats_in$x,lonlats_in$y)$seqnum
      lonlats_in <- cbind(lonlats_in,ev_vals)
      lonlats_in <- filter(lonlats_in,cell==cl_target)
      
      cell_data_add <- as.data.frame(1:365)
      cell_data_add$year <- targ_year
      cell_data_add$enviro_avg <- colMeans(lonlats_in,na.rm=TRUE)[4:368]
      
      colnames(cell_data_add) <- c("day","year",param_targ)
      
      #Delete temporary file
      unlink(paste(folder_path,"/",param_targ,"_daily_",targ_year,"_ncss.nc",sep=""))
      cell_data_all <- rbind(cell_data_all,cell_data_add)
      #print(nrow(cell_data_all))
    }
    #Save
    saveRDS(cell_data_all,paste("data/enviro_var_daily/",param_targ,"_cell/","cl_",cl_target,"_",param_targ,".rds",sep=""))
  }
  parallel::stopCluster(cl = my.cluster)
}
