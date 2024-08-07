# Code to prepare data for a combined stan file that simultaneously
# estimates breeding weather effects and ENSO/NAO. 

library(dplyr)
library(readr)
library(rstan)
library(MCMCvis)
library(shinystan)
library(dggridR)
library(rsoi)

#Read in WB and LB data
all_data <- readRDS("data/mdl_data/envir_phen_comb.rds")

nrow(all_data)

all_data <- all_data[order(all_data$species),]

#Group by
all_data_trim <- c()
for (each_spec in unique(all_data$species)){
  #targ_spec <- "Warbling Vireo"
  targ_spec <- each_spec
  
  #Filter to specific species, filter out estimates from Jan.
  spec_data <- filter(all_data,species==targ_spec & arr_GAM_mean > 31) %>% 
    select(c("species","cell","year","arr_GAM_mean","arr_GAM_sd",
             "mean_tmin_29","mean_prcp_29","mean_swe_29",
             "anom_tmin_c29","anom_prcp_c29","anom_swe_c29",
             "tmin_c29","prcp_c29","swe_c29",
             "sd_tmin_c29","sd_prcp_c29","sd_swe_c29"))
  
  # Get the range species niche data by comparing the average condition in each cell in comparison to
  # all other cells
  
  #First get the unique cells and 
  cl_env_dt <- select(spec_data,c("cell","mean_tmin_29","mean_prcp_29","mean_swe_29")) %>% unique()
  
  cl_env_dt$tmin_avg_std <- (cl_env_dt$mean_tmin_29 - mean(cl_env_dt$mean_tmin_29))/sd(cl_env_dt$mean_tmin_29)
  cl_env_dt$prcp_avg_std <- (cl_env_dt$mean_prcp_29 - mean(cl_env_dt$mean_prcp_29))/sd(cl_env_dt$mean_prcp_29)
  cl_env_dt$swe_avg_std <- (cl_env_dt$mean_swe_29 - mean(cl_env_dt$mean_swe_29))/sd(cl_env_dt$mean_swe_29)
  cl_env_dt <- select(cl_env_dt,-c("mean_tmin_29","mean_prcp_29","mean_swe_29"))
  
  spec_data <- merge(spec_data,cl_env_dt,by="cell")
  
  #Filter to species >50 yr/cl estimates
  if (nrow(spec_data) >= 50){
    all_data_trim <- rbind(all_data_trim,spec_data)
  }
}

### Now download NAO ####
#nao_data <- download_nao() # broken for some reason
download.file("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table","nao_dt.txt")
nao_data <- read.table("nao_dt.txt",fill = NA)

nao_data$Year <- as.numeric(rownames(nao_data))

nao_data <- reshape2::melt(nao_data,id.vars="Year",variable.names = "Month",value.name = "NAO")
colnames(nao_data)[2] <- "Month"

#Add one to all decembers to account for cross year mark
nao_data$Year[which(nao_data$Month=="Dec")] <- nao_data$Year[which(nao_data$Month=="Dec")] + 1

nao_months <- nao_data %>% filter(Month %in% c("Dec","Jan","Feb")) %>%
  group_by(Year) %>%
  summarise(NAO = mean(NAO))

all_data_trim <- merge(all_data_trim,nao_months,by.x="year",by.y="Year")

#Get cell correlations to enviro variables
#nao_cl_cor <- readRDS("nao_cl_cor.rds")
#colnames(nao_cl_cor) <- c("cell","cor_prcp_nao","cor_tmax_nao","cor_swe_nao")

#all_data_trim <- merge(all_data_trim,nao_cl_cor,by="cell")

### Now download ENSO ####
#enso_data <- download_enso()

download.file("https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/nino34.long.anom.data","enso_dt.txt")
enso_data <- read.table("enso_dt.txt",fill = NA,skip = 1)
colnames(enso_data) <- c("Year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
enso_data <- filter(enso_data,Year <= 2022 & Year >= 1979)
enso_data <- enso_data %>% mutate_if(is.character, as.numeric)
enso_data <- reshape2::melt(enso_data,id.vars="Year",variable.names = "Month",value.name = "ENSO")
colnames(enso_data)[2] <- "Month"

#Add one to all decembers to account for cross year mark
enso_data$Year[which(enso_data$Month=="Dec")] <- enso_data$Year[which(enso_data$Month=="Dec")] + 1

enso_months <- enso_data %>% filter(Month %in% c("Dec","Jan","Feb")) %>%
  group_by(Year) %>%
  summarise(ENSO = mean(ENSO))

all_data_trim <- merge(all_data_trim,enso_months,by.x="year",by.y="Year")

#Get cell correlations to enviro variables
#enso_cl_cor <- readRDS("enso_cl_cor.rds")
#colnames(enso_cl_cor) <- c("cell","cor_prcp_enso","cor_tmax_enso","cor_swe_enso")

#all_data_trim <- merge(all_data_trim,enso_cl_cor,by="cell")

### Get winter associations ####
#Get species wintering ranges associations with CO variables
sp_co_rec <- readRDS("data/mdl_data/sp_co_rec.rds")

#Filter to ENSO
enso_temp <- sp_co_rec %>% filter(co_var=="ENSO_Temperature") %>% select(c(spec,median_impact))
colnames(enso_temp)[2] <- "enso_temp_wint"
enso_prcp <- sp_co_rec %>% filter(co_var=="ENSO_Precipitation") %>% select(c(spec,median_impact))
colnames(enso_prcp)[2] <- "enso_prcp_wint"

all_data_trim <- merge(all_data_trim,enso_temp,by.x="species",by.y="spec")
all_data_trim <- merge(all_data_trim,enso_prcp,by.x="species",by.y="spec")

#Filter to NAO
nao_temp <- sp_co_rec %>% filter(co_var=="NAO_Temperature") %>% select(c(spec,median_impact))
colnames(nao_temp)[2] <- "nao_temp_wint"
nao_prcp <- sp_co_rec %>% filter(co_var=="NAO_Precipitation") %>% select(c(spec,median_impact))
colnames(nao_prcp)[2] <- "nao_prcp_wint"

all_data_trim <- merge(all_data_trim,nao_temp,by.x="species",by.y="spec")
all_data_trim <- merge(all_data_trim,nao_prcp,by.x="species",by.y="spec")

#all_data_trim$cor_swe_nao[which(is.na(all_data_trim$cor_swe_nao))] <- 0
#all_data_trim$cor_swe_enso[which(is.na(all_data_trim$cor_swe_enso))] <- 0

# Get unique cell, species combos
all_data_trim <- all_data_trim %>% arrange(species,cell)
cl_sp_combos <- unique(data.frame(cell=all_data_trim$cell,species=all_data_trim$species))
cl_sp_combos$cl_sp <- 1:nrow(cl_sp_combos)

sp_cl_in <- as.numeric(as.factor(cl_sp_combos$species))
cl_in <- as.numeric(as.factor(cl_sp_combos$cell))

all_data_trim <- merge(all_data_trim,cl_sp_combos,by=c("cell","species"))
all_data_trim <- all_data_trim %>% arrange(species)

temp_niche <- all_data_trim %>% select(c("cell","species","tmin_avg_std")) %>% unique()
prcp_niche <- all_data_trim %>% select(c("cell","species","prcp_avg_std")) %>% unique()
swe_niche <- all_data_trim %>% select(c("cell","species","swe_avg_std")) %>% unique()

nrow(all_data_trim)
table(all_data_trim$species)

#ENSO Variables
sp_wint_enso_temp <- as.numeric(unique(cbind(all_data_trim$species,all_data_trim$enso_temp_wint))[,2])
sp_wint_enso_prcp <- as.numeric(unique(cbind(all_data_trim$species,all_data_trim$enso_prcp_wint))[,2])

enso_pc <- prcomp(cbind(sp_wint_enso_temp,sp_wint_enso_prcp),retx = TRUE)
enso_pca1 <- -enso_pc$x[,1]
cor(enso_pca1,sp_wint_enso_temp)^2
cor(enso_pca1,sp_wint_enso_prcp)^2


#NAO variables
sp_wint_nao_temp <- as.numeric(unique(cbind(all_data_trim$species,all_data_trim$nao_temp_wint))[,2])
sp_wint_nao_prcp <- as.numeric(unique(cbind(all_data_trim$species,all_data_trim$nao_prcp_wint))[,2])

nao_pc <- prcomp(cbind(sp_wint_nao_temp,sp_wint_nao_prcp),retx = TRUE)
nao_pca1 <- nao_pc$x[,1]
cor(nao_pca1,sp_wint_nao_temp)^2
cor(nao_pca1,sp_wint_nao_prcp)^2

# Add "snow present" variable

snow_p <- 1 * (all_data_trim$mean_swe_29 > 0.1)
sum(snow_p)/length(snow_p)

all_data_trim$snow_p <- snow_p

#Add latitudes of cell centers - centered within species
grid_size<- dgconstruct(res=6)
all_data_trim$cell_lat  <- dgSEQNUM_to_GEO(grid_size,all_data_trim$cell)$lat_deg
sp_cl_lat <- all_data_trim %>% 
  group_by(species) %>% 
  select(species,cl_sp,cell_lat) %>% 
  summarise(mean_sp_cl_lat = mean(cell_lat))

all_data_trim <- merge(all_data_trim,sp_cl_lat,by="species")
all_data_trim$cell_lat_centered <- all_data_trim$cell_lat - all_data_trim$mean_sp_cl_lat

#### Run Model ####
data_in1 <- list(
  N = nrow(all_data_trim),
  
  yy = (all_data_trim$arr_GAM_mean-mean(all_data_trim$arr_GAM_mean)),
  yy_unc = all_data_trim$arr_GAM_sd,
  
  temp = all_data_trim$anom_tmin_c29,
  precip = all_data_trim$anom_prcp_c29,
  swe = all_data_trim$anom_swe_c29,
  
  snow_p = snow_p,
  
  enso = all_data_trim$ENSO,
  nao = all_data_trim$NAO,
  
  cell_lat = all_data_trim$cell_lat_centered,

  Ncell = length(unique(all_data_trim$cell)),
  cc = as.numeric(as.factor(all_data_trim$cell)),
  
  Nsp = length(unique(all_data_trim$species)),
  ii = as.numeric(as.factor(all_data_trim$species)),
  
  N_sp_cl = length(unique(all_data_trim$cl_sp)),
  sp_cl = as.numeric(as.factor(all_data_trim$cl_sp)),
  sp_cl_in = sp_cl_in,
  cl_in = cl_in,
  
  temp_niche = temp_niche$tmin_avg_std,
  prcp_niche = prcp_niche$prcp_avg_std,
  swe_niche = swe_niche$swe_avg_std,
  
  sp_enso_pca1 = enso_pca1,
  sp_nao_pca1 = nao_pca1
  
)

saveRDS(data_in1,"data/mdl_data/ms_data1.rds")
saveRDS(all_data_trim,"data/mdl_data/ms_data2.rds")

