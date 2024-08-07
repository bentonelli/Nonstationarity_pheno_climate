library(readr)
library(zoo)
library(dplyr)

#Script to read what cells have full data, then compile environmental and phenology data

all_data <- readRDS("data/sp_yr_files/sp_yr_cl_dt.rds")
cl_list <- unique(all_data$cell)

tmin_lst <- list.files("data/enviro_var_daily/tmin_cell/")
tmin_lst_cl <- gsub("cl_","",tmin_lst)
tmin_lst_cl <- as.numeric(gsub("_tmin.rds","",tmin_lst_cl))

cl_list[!(cl_list %in% tmin_lst_cl)]

prcp_lst <- list.files("data/enviro_var_daily/prcp_cell/")
prcp_lst_cl <- gsub("cl_","",prcp_lst)
prcp_lst_cl <- as.numeric(gsub("_prcp.rds","",prcp_lst_cl))

cl_list[!(cl_list %in% prcp_lst_cl)]

swe_lst <- list.files("data/enviro_var_daily/swe_cell/")
swe_lst_cl <- gsub("cl_","",swe_lst)
swe_lst_cl <- as.numeric(gsub("_swe.rds","",swe_lst_cl))

cl_list[!(cl_list %in% swe_lst_cl)]

### Combine enviro data ####
#Get cells for each file
tmin_cl_names <- Sys.glob("data/enviro_var_daily/tmin_cell/*.rds")
tmin_cl_names <- gsub("data/enviro_var_daily/tmin_cell/cl_","",tmin_cl_names)
tmin_cl_names <- as.numeric(gsub("_tmin.rds","",tmin_cl_names))

#Get list of files
tmin_lst_set <- lapply(Sys.glob("data/enviro_var_daily/tmin_cell/*.rds"), readRDS)
for (each_lst in 1:length(tmin_lst_set)){
  #Add cells to each file
  tmin_lst_set[[each_lst]]$cell <- tmin_cl_names[each_lst]
  #Add rolling average - ignore NAs and center rolling average
  tmin_lst_set[[each_lst]]$tmin_c29 <- rollapply(tmin_lst_set[[each_lst]]$tmin, 
                                                 width=29, FUN=function(x) mean(x, na.rm=TRUE),
                                                 by=1, partial=TRUE, fill=NA, align="center")
  #Filter to study period
  tmin_lst_set[[each_lst]] <- filter(tmin_lst_set[[each_lst]],year > 2001)
  
  tmin_anom_cl <- tmin_lst_set[[each_lst]] %>% 
    select(day,tmin_c29) %>%
    group_by(day) %>% 
    summarise(mean_tmin_29 = mean(tmin_c29),sd_tmin_c29 = sd(tmin_c29))
  tmin_lst_set[[each_lst]] <- merge(tmin_lst_set[[each_lst]],tmin_anom_cl,by="day")
  tmin_lst_set[[each_lst]]$anom_tmin_c29 <- (tmin_lst_set[[each_lst]]$tmin_c29 - 
                                               tmin_lst_set[[each_lst]]$mean_tmin_29)/tmin_lst_set[[each_lst]]$sd_tmin_c29
}
#Combine
tmin_lst_df <- do.call(bind_rows,tmin_lst_set)

# Do for each enviro variable
#Get cells for each file
prcp_cl_names <- Sys.glob("data/enviro_var_daily/prcp_cell/*.rds")
prcp_cl_names <- gsub("data/enviro_var_daily/prcp_cell/cl_","",prcp_cl_names)
prcp_cl_names <- as.numeric(gsub("_prcp.rds","",prcp_cl_names))

#Get list of files
prcp_lst_set <- lapply(Sys.glob("data/enviro_var_daily/prcp_cell/*.rds"), readRDS)
for (each_lst in 1:length(prcp_lst_set)){
  #Add cells to each file
  prcp_lst_set[[each_lst]]$cell <- prcp_cl_names[each_lst]
  #Add rolling average - ignore NAs and center rolling average
  prcp_lst_set[[each_lst]]$prcp_c29 <- rollapply(prcp_lst_set[[each_lst]]$prcp, 
                                                 width=29, FUN=function(x) mean(x, na.rm=TRUE),
                                                 by=1, partial=TRUE, fill=NA, align="center")
  
  #Filter to study period
  prcp_lst_set[[each_lst]] <- filter(prcp_lst_set[[each_lst]],year > 2001)
  
  prcp_anom_cl <- prcp_lst_set[[each_lst]] %>% 
    select(day,prcp_c29) %>%
    group_by(day) %>% 
    summarise(mean_prcp_29 = mean(prcp_c29),sd_prcp_c29 = sd(prcp_c29))
  prcp_lst_set[[each_lst]] <- merge(prcp_lst_set[[each_lst]],prcp_anom_cl,by="day")
  prcp_lst_set[[each_lst]]$anom_prcp_c29 <- (prcp_lst_set[[each_lst]]$prcp_c29 - 
                                               prcp_lst_set[[each_lst]]$mean_prcp_29)/prcp_lst_set[[each_lst]]$sd_prcp_c29
}
#Combine
prcp_lst_df <- do.call(bind_rows,prcp_lst_set)

#Get cells for each file
swe_cl_names <- Sys.glob("data/enviro_var_daily/swe_cell/*.rds")
swe_cl_names <- gsub("data/enviro_var_daily/swe_cell/cl_","",swe_cl_names)
swe_cl_names <- as.numeric(gsub("_swe.rds","",swe_cl_names))

#Get list of files
swe_lst_set <- lapply(Sys.glob("data/enviro_var_daily/swe_cell/*.rds"), readRDS)
for (each_lst in 1:length(swe_lst_set)){
  #Add cells to each file
  swe_lst_set[[each_lst]]$cell <- swe_cl_names[each_lst]
  #Add rolling average - ignore NAs and center rolling average
  swe_lst_set[[each_lst]]$swe_c29 <- rollapply(swe_lst_set[[each_lst]]$swe, 
                                               width=29, FUN=function(x) mean(x, na.rm=TRUE),
                                               by=1, partial=TRUE, fill=NA, align="center")
  
  #Filter to study period
  swe_lst_set[[each_lst]] <- filter(swe_lst_set[[each_lst]],year > 2001)
  
  swe_anom_cl <- swe_lst_set[[each_lst]] %>% 
    select(day,swe_c29) %>%
    group_by(day) %>% 
    summarise(mean_swe_29 = mean(swe_c29),sd_swe_c29 = sd(swe_c29))
  swe_lst_set[[each_lst]] <- merge(swe_lst_set[[each_lst]],swe_anom_cl,by="day")
  swe_lst_set[[each_lst]]$anom_swe_c29 <- (swe_lst_set[[each_lst]]$swe_c29 - 
                                             swe_lst_set[[each_lst]]$mean_swe_29)/swe_lst_set[[each_lst]]$sd_swe_c29
}

#Combine
swe_lst_df <- do.call(bind_rows,swe_lst_set)


comb_enviro_df <- merge(tmin_lst_df,prcp_lst_df,by=c("cell","year","day"))
comb_enviro_df <- merge(comb_enviro_df,swe_lst_df,by=c("cell","year","day"))

#We want to trim the 2001 data, as the rolling mean algo filled those in
#without much data

comb_enviro_df_trim <- filter(comb_enviro_df,year != 2001)
saveRDS(comb_enviro_df_trim,"data/mdl_data/enviro_by_cell.rds")

### Calculate average arrival for each sp_cl combo, add enviro data based on time window ####

avg_sp_cl <- all_data %>% 
  group_by(species,cell) %>% 
  summarise(mean_arr_date = mean(arr_GAM_mean))

all_data_w_mean <- merge(all_data,avg_sp_cl,by=c("species","cell"))

#About 99% of all arrival dates fall within +- 12 days (14 for waterbirds).
quantile((all_data_w_mean$arr_GAM_mean-all_data_w_mean$mean_arr_date),probs=c(.01,.99))

#Round to calendar day
all_data_w_mean$day <- round(all_data_w_mean$mean_arr_date)

#Now add enviro data
envir_phen_comb <- merge(all_data_w_mean,comb_enviro_df_trim,by=c("cell","year","day"),)
envir_phen_comb <- envir_phen_comb %>% 
  arrange(day) %>% 
  arrange(year) %>% 
  arrange(species)

#Remove the few cells with no data
envir_phen_comb <- envir_phen_comb[which(!is.na(envir_phen_comb$tmin_c29)),]

#For places with no snow, change NAs to 0
envir_phen_comb$anom_swe_c29[which(envir_phen_comb$mean_swe_29==0)] <- 0
#Save
saveRDS(envir_phen_comb,"data/mdl_data/envir_phen_comb.rds")
