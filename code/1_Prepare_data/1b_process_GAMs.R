# Code to get arrival estimates for species-cell-years
# Adapted from Youngflesh, 2021

library(nloptr)
library(raster)
library(lubridate)
library(dggridR)
library(rstan)
library(rstanarm)
library(MCMCvis)
library(dplyr)

cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

print(args)
year_list <- as.numeric(args[[2]][2])
year_list <- 2012:2021

spec_list <- c("Yellow-bellied Sapsucker")
for (year in year_list){
  #Create df for record keeping
  print(year)
  #Read in year file
  f_name <- paste("data/year_files/",year,"_all.txt",sep="")
  data_in <- read.delim(f_name,sep = "|",na.strings = "")
  
  #only keep needed data
  data_in <- dplyr::select(data_in,c(COMMON.NAME,LATITUDE,LONGITUDE,OBSERVATION.DATE,SAMPLING.EVENT.IDENTIFIER,DURATION.MINUTES,GROUP.IDENTIFIER))
  
  #Set up grids
  hexgrid6 <- dggridR::dgconstruct(res = 6) 
  data_in$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                           in_lon_deg = data_in$LONGITUDE, 
                                           in_lat_deg = data_in$LATITUDE)[[1]]
  
  #Read in unique checklists
  grid_cell_checklist <- unique(dplyr::select(data_in,c("SAMPLING.EVENT.IDENTIFIER","cell","DURATION.MINUTES")))
  
  #Change effort from minutes to hours
  grid_cell_checklist$DURATION.MINUTES <- grid_cell_checklist$DURATION.MINUTES/60
  colnames(grid_cell_checklist)[c(1,3)] <- c("checklist_num","EH")
  
  #Get unique sampking events by removing duplicates
  sampling_events <- select(data_in,c("SAMPLING.EVENT.IDENTIFIER","GROUP.IDENTIFIER"))
  sampling_events <- sampling_events[!duplicated(sampling_events[,"GROUP.IDENTIFIER"], incomparables = NA),]
  sampling_events <- unique(select(sampling_events,"SAMPLING.EVENT.IDENTIFIER"))
  
  #Add julian days, lat long
  julian_days <- unique(dplyr::select(data_in,c("SAMPLING.EVENT.IDENTIFIER","OBSERVATION.DATE","LATITUDE","LONGITUDE")))
  julian_days$julian_day <- yday(julian_days$OBSERVATION.DATE)
  sampling_events <- merge(sampling_events,julian_days,by="SAMPLING.EVENT.IDENTIFIER")
  
  #Filter to checklists where julian day < 200, rename columns
  sampling_events <- sampling_events %>%
    filter(julian_day < 200) %>%
    dplyr::select(c("SAMPLING.EVENT.IDENTIFIER","julian_day","LATITUDE","LONGITUDE"))
  colnames(sampling_events) <- c("checklist_num","julian_day","LATITUDE","LONGITUDE") 
  
  #Merge together with cell identities
  spec_matrix <- merge(sampling_events,grid_cell_checklist,by="checklist_num")
  
  #Get elevation
  lat_longs <- dplyr::select(spec_matrix,c("LONGITUDE","LATITUDE"))
  coordinates(lat_longs)= ~ LONGITUDE+LATITUDE
  # Read in elevation file (can be downloaded here: https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation)
  fff <- raster(x = "data/GMTED_2010/medians/NA_30AS_Elev_Med.tif")
  spec_matrix$elev <- raster::extract(fff,lat_longs)
  
  # There are some NAs in the great lakes region from the raster - 
  # it should only happen for a minuscule number of checklists (if any).
  # In these cases, we want to set to sea level.
  if (sum(is.na(spec_matrix$elev)>0)){
    #Print out if there are NAs
    print(paste("NA cells:",sum(is.na(spec_matrix$elev)),sep = " "))
    spec_matrix$elev[which(is.na(spec_matrix$elev))] <- 0
  }
  
  #Remove Lat-long, no longer needed
  spec_matrix <- dplyr::select(spec_matrix,c("checklist_num","julian_day","cell","EH","elev"))
  
  #For each species, add columns, and P/NP information for each checklist
  for (spec_name in spec_list){
    #Add species columns
    spec_data <- data_in %>% 
      filter(COMMON.NAME == spec_name) %>%
      summarise(checklist_num = unique(SAMPLING.EVENT.IDENTIFIER),p = 1)
    
    spec_matrix <- merge(spec_matrix,spec_data,by="checklist_num", all.x = TRUE)
    spec_matrix$p[which(is.na(spec_matrix$p))] <- 0
    colnames(spec_matrix)[ncol(spec_matrix)] <- spec_name
  }
  
  #Save total number of species
  num_species <- length(spec_list)
  
  #Set up stan model conditions
  DELTA <- .95
  TREE_DEPTH <- 10
  ITER <- 1500
  CHAINS <- 4
  
  #Clear memory
  lat_longs <- fff <- julian_days <- sampling_events <- grid_cell_checklist <- data_in <- grid_cell_checklist <- NULL
  
  for (each_spec in 1:num_species){
    arrival_df_all <- c()
    #Counter for cell-species-year combos
    counter <- 0
    #Index for the column where species information is given
    spec_ind <- each_spec + 5
    
    print(paste(spec_list[each_spec],each_spec,sep=" "))
    print(spec_ind)
    
    # For each cell, check:
    # Is species reported on >= 1% of checklists?
    # Are <= 2% of detections from before julian_day 60?
    # Are 20+ days of non-detection present before first detection?
    # Detections on >=20 days
    
    unique_cells <- unique(spec_matrix$cell)
    GAM_cells <- c()
    low_prob_cells <- c()
    for (each_cell in unique_cells){
      cell_data <- spec_matrix %>% filter(cell == each_cell) %>% dplyr::select(c(1,2,3,4,5,spec_ind))
      
      total_detections <- sum(cell_data[,6])
      
      if (total_detections > 0){
        
        perc_checklists_reported <- sum(cell_data[,6])/nrow(cell_data)
        
        first_detection <- min(cell_data$julian_day[which(cell_data[,6]==1)])
        
        non_detections <- cell_data$julian_day[which(cell_data[,6] == 0)]
        
        num_detection_dates <- length(unique(cell_data$julian_day[which(cell_data[,6]==1)]))
        
        detections_before_60 <- sum(cell_data$julian_day[which(cell_data[,6]==1)] < 60)/total_detections
        
        nd_before_first <- sum(non_detections < first_detection)
        
        if(detections_before_60 <= .02 & nd_before_first >= 20 & num_detection_dates >= 20){
          #Add extra if layer to save cells that are eliminated because they have a low overall probability
          if(perc_checklists_reported >= .01){
            GAM_cells <- c(GAM_cells,each_cell)
          } else {
            low_prob_cells <- c(low_prob_cells,each_cell)
          }
        }
      }
    }
    #Check how many cells for given species
    print(paste("Number of cells: ",length(GAM_cells),sep=""))
    print(paste("Low prob cells removed: ",length(low_prob_cells),sep=""))
    #As long as there are cells to run, create arrival df to save model output
    if (length(GAM_cells) > 0 ){
      hm_mat <- matrix(data = NA, nrow = length(GAM_cells), ncol = ((ITER/2)*CHAINS))
      max_mat <- matrix(data = NA, nrow = length(GAM_cells), ncol = ((ITER/2)*CHAINS))
      colnames(hm_mat) <- paste0('hm_iter_', 1:((ITER/2)*CHAINS))
      colnames(max_mat) <- paste0('max_iter_', 1:((ITER/2)*CHAINS))
      arrival_df <- data.frame(species = spec_list[each_spec], 
                               year = rep(year, each = length(GAM_cells)), 
                               cell = GAM_cells,
                               num_og_records = NA,
                               max_Rhat = NA,
                               min_neff = NA,
                               mlmax = NA,
                               plmax = NA,
                               num_diverge = NA,
                               num_tree = NA,
                               num_BFMI = NA,
                               delta = NA,
                               tree_depth = NA,
                               t_iter = NA,
                               elev = NA,
                               arr_GAM_mean = NA,
                               arr_GAM_sd = NA,
                               hm_mat,
                               max_mat)
      
      #For each cell to run a model in
      for(each_cell in GAM_cells){
        
        
        counter <- counter + 1
        print(paste("Running #",counter," of ",length(GAM_cells), " cells"))
        
        DELTA <- 0.95
        TREE_DEPTH <- 15
        
        #Filter data to sp-cell-yr
        cell_data <- filter(spec_matrix,cell==each_cell) %>% dplyr::select(c(1,2,3,4,5,spec_ind))
        colnames(cell_data)[6] <- "detect"
        og_number_rec <- nrow(cell_data)
        #If there are >500 detection records, thin to ~500 positive records.
        #This drastically reduces runtimes for high volume cell-year combos
        
        if (nrow(cell_data) > 5000){
          
          print(paste("Thinning from:",og_number_rec,"records"))
          
          targ_records <- 5000
          
          cell_data <- cell_data[sample(1:nrow(cell_data),targ_records),]
          reduction <- (1-(nrow(cell_data)/og_number_rec))*100
          print(paste("Thinned to",nrow(cell_data),"records,","a",round(reduction,2),"% reduction"))
          print(paste("Total 1s:",sum(cell_data[,6])))
        }
        
        #Scale effort
        cell_data$shr <- scale(cell_data$EH, scale = FALSE)[,1]
        
        # Get elevation as 500m intervals: 
        # 0) <500m,1) 500m-1km,2) 1km-1.5km,
        # 3) 1.5km-2km,4) 2km-2.5km,5) 2.5-3km,
        # 6) >3km
        cell_data$elev_factor <- floor(cell_data$elev/500)
        cell_data$elev_factor[which(cell_data$elev_factor < 0)] <- 0
        cell_data$elev_factor[which(cell_data$elev_factor > 6)] <- 6
        
        #Save as factor
        cell_data$elev_factor <- as.factor(cell_data$elev_factor)
        
        print(paste("Elev. levels: ",length(table(cell_data$elev_factor)),sep=""))
        
        #As long as there are checklists from >1 factor level
        if (length(table(cell_data$elev_factor)) > 1){
          fit2 <- rstanarm::stan_gamm4(detect ~ s(julian_day, k = 30, by = elev_factor) + shr + elev_factor,
                                       data = cell_data,
                                       family = binomial(link = "logit"),
                                       algorithm = 'sampling',
                                       iter = ITER,
                                       chains = CHAINS,
                                       cores = CHAINS,
                                       adapt_delta = .95,
                                       control = list(max_treedepth = 10))
          #If not, we can just run a simpler version of the model
        } else {
          fit2 <- rstanarm::stan_gamm4(detect ~ s(julian_day, k = 30) + shr,
                                       data = cell_data,
                                       family = binomial(link = "logit"),
                                       algorithm = 'sampling',
                                       iter = ITER,
                                       chains = CHAINS,
                                       cores = CHAINS,
                                       adapt_delta = .95,
                                       control = list(max_treedepth = 10))
        }
        
        #Save model output
        num_diverge <- rstan::get_num_divergent(fit2$stanfit)
        num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
        num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
        max_Rhat <- round(max(summary(fit2)[, 'Rhat']), 2)
        min_neff <- min(summary(fit2)[, 'n_eff'])
        
        #If there are divergences, try, try again
        while (num_diverge > 0 & DELTA <= 0.98){
          DELTA <- DELTA + 0.01
          
          if (length(table(cell_data$elev_factor)) > 1){
            fit2 <- rstanarm::stan_gamm4(detect ~ s(julian_day, k = 30, by = elev_factor) + shr + elev_factor,
                                         data = cell_data,
                                         family = binomial(link = "logit"),
                                         algorithm = 'sampling',
                                         iter = ITER,
                                         chains = CHAINS,
                                         cores = CHAINS,
                                         adapt_delta = .95,
                                         control = list(max_treedepth = 10))
          } else {
            fit2 <- rstanarm::stan_gamm4(detect ~ s(julian_day, k = 30) + shr,
                                         data = cell_data,
                                         family = binomial(link = "logit"),
                                         algorithm = 'sampling',
                                         iter = ITER,
                                         chains = CHAINS,
                                         cores = CHAINS,
                                         adapt_delta = .95,
                                         control = list(max_treedepth = 10))
          }
          
          num_diverge <- rstan::get_num_divergent(fit2$stanfit)
          num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
          num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
        }
        arrival_df$num_og_records[counter] <- og_number_rec
        arrival_df$num_diverge[counter] <- num_diverge
        arrival_df$num_tree[counter] <- num_tree
        arrival_df$num_BFMI[counter] <- num_BFMI
        arrival_df$delta[counter] <- DELTA
        arrival_df$tree_depth[counter] <- TREE_DEPTH
        arrival_df$t_iter[counter] <- ITER
        arrival_df$max_Rhat[counter] <- max_Rhat
        arrival_df$min_neff[counter] <- min_neff
        
        #generate predict data
        predictDays <- range(cell_data$julian_day)[1]:range(cell_data$julian_day)[2]
        
        # Get the most common elevation factor, use that to make predictions
        elev_tbl <- table(cell_data$elev_factor)
        elev_mode <- as.factor(names(which(elev_tbl == max(elev_tbl))))
        
        #If there is an even split, pick the lower elevation
        if (length(elev_mode) > 1){
          elev_mode <- elev_mode[1]
        }
        
        arrival_df$elev[counter] <- as.character(elev_mode)
        
        if (length(table(cell_data$elev_factor)) > 1){
          #Want to get arrival at most common elevation
          newdata <- data.frame(julian_day = predictDays, shr = 0,elev_factor=elev_mode)
        } else {
          #If elevation is only in one bucket:
          newdata <- data.frame(julian_day = predictDays, shr = 0)
        }
        
        
        #predict response
        dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
        halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
        max_fit <- rep(NA, ((ITER/2)*CHAINS))
        tlmax <- rep(NA, ((ITER/2)*CHAINS))
        #day at which probability of occurence is half local maximum value
        for (L in 1:((ITER/2)*CHAINS)){
          rowL <- as.vector(dfit[L,])
          #first detection
          fd <- min(cell_data$julian_day[which(cell_data$detect == 1)])
          #local maximum(s)
          #from: stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
          lmax_idx <- which(diff(sign(diff(rowL))) == -2) + 1
          lmax <- predictDays[lmax_idx]
          #first local max to come after first detection
          flm <- which(lmax > fd)
          if (length(flm) > 0)
          {
            #first local max to come after first detection
            lmax2_idx <- lmax_idx[min(flm)]
            lmax2 <- lmax[min(flm)]
            tlmax[L] <- TRUE
          } else {
            #no local max
            lmax2_idx <- which.max(rowL)
            lmax2 <- predictDays[which.max(rowL)]
            tlmax[L] <- FALSE
          }
          #store local max
          max_fit[L] <- lmax2 
          #local mins before max (global and local mins)
          lmin_idx <- c(which.min(rowL[1:lmax2_idx]), 
                        which(diff(sign(diff(rowL[1:lmax2_idx]))) == 2) + 1)
          lmin <- predictDays[lmin_idx]
          #local min nearest to local max
          lmin2_idx <- lmin_idx[which.min(lmax2 - lmin)]
          lmin2 <- predictDays[lmin2_idx]
          
          #value at local max - value at min (typically 0)
          dmm <- rowL[lmax2_idx] - rowL[lmin2_idx]
          #all positions less than or equal to half diff between max and min + value min
          tlm <- which(rowL <= ((dmm/2) + rowL[lmin2_idx]))
          #which of these come before max and after or at min
          
          vgm <- tlm[which(tlm < lmax2_idx & tlm >= lmin2_idx)]
          #insert halfmax (first day for situations where max is a julian_day = 1)
          if (length(vgm) > 0)
          {
            halfmax_fit[L] <- predictDays[max(vgm)]
          } else {
            halfmax_fit[L] <- predictDays[1]
          }
        }
        
        #number of iterations that had local max
        arrival_df$plmax[counter] <- round(sum(tlmax)/((ITER/2)*CHAINS), 3)
        
        #model fit
        mn_dfit <- apply(dfit, 2, mean)
        LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
        UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
        
        #check whether local max exists for mean model fit
        mlmax <- sum(diff(sign(diff(mn_dfit))) == -2)
        if (mlmax > 0){
          arrival_df$mlmax[counter] <- TRUE
        } else {
          arrival_df$mlmax[counter] <- FALSE
        }
        
        #estimated halfmax
        mn_hm <- mean(halfmax_fit)
        LCI_hm <- quantile(halfmax_fit, probs = 0.025)
        UCI_hm <- quantile(halfmax_fit, probs = 0.975)
        
        #estimated max
        mn_max <- mean(max_fit)
        LCI_max <- quantile(max_fit, probs = 0.025)
        UCI_max <- quantile(max_fit, probs = 0.975)
        
        #fill df with halfmax iter
        cndf <- colnames(arrival_df)
        hm_iter_ind <- grep('hm_iter', cndf)
        arrival_df[counter, hm_iter_ind] <- halfmax_fit
        
        #fill df with max iter
        max_iter_ind <- grep('max_iter', cndf)
        arrival_df[counter, max_iter_ind] <- max_fit
        
        arrival_df[counter,"arr_GAM_mean"] <- mean(halfmax_fit)
        arrival_df[counter,"arr_GAM_sd"] <- sd(halfmax_fit)
        
        ########################
        #PLOT MODEL FIT AND DATA
        
        pdf(paste0("data/spec_hm_pdfs/",spec_list[each_spec], '_', year, '_',each_cell, '_arrival.pdf'))
        plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
             ylim = c(0, max(UCI_dfit)),
             main = paste0(spec_list[each_spec], ' - ', year, ' - ', each_cell,' - ',as.character(elev_mode)),
             xlab = 'Julian Day', ylab = 'Probability of occurrence')
        lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
        lines(predictDays, mn_dfit, lwd = 2)
        dd <- cell_data$detect
        dd[which(dd == 1)] <- max(UCI_dfit)
        points(cell_data$julian_day, dd, col = rgb(0,0,0,0.25))
        abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
        abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
        abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
        abline(v = mn_max, col = rgb(0,1,0,0.5), lwd = 2)
        abline(v = LCI_max, col = rgb(0,1,0,0.5), lwd = 2, lty = 2)
        abline(v = UCI_max, col = rgb(0,1,0,0.5), lwd = 2, lty = 2)
        legend('topleft',
               legend = c('Model fit', 'CI fit', 'Half max', 'CI HM'),
               col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
               lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
        dev.off()
      }
      #Add to full record
      arrival_df_all <- rbind(arrival_df_all,arrival_df)
    }
    arrival_df_all[,1:16]
    spec_name_gsub <- gsub(" ","_",spec_list[each_spec])
    #If any GAMs were succesfull, save
    if (!is.null(nrow(arrival_df_all))){
      f_name <- paste("data/spec_yr_files/",spec_name_gsub,"_",year,".rds",sep="")
      saveRDS(object = arrival_df_all,f_name)
    }
  }
}

