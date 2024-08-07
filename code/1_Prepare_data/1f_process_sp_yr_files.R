library(dplyr)
library(dggridR)
library(mapdata)
library(ggplot2)
library(rstan)
library(MCMCvis)
library(shinystan)
library(readr)

#Bring in all species, year files
all_yr_fits <- lapply(Sys.glob("data/sp_yr_files/*.rds"), readRDS)
compiled_yrs <- do.call(bind_rows,all_yr_fits)
all_yr_fits <- c()
nrow(compiled_yrs) #154k

table(compiled_yrs$species,compiled_yrs$year)
table(compiled_yrs$year)

#Filter out bad GAM fits
compiled_yrs <- filter(compiled_yrs,max_Rhat < 1.02 & min_neff >= 400 & 
                         mlmax == TRUE & plmax >=.99 & num_diverge == 0 & arr_GAM_sd <= 15)
nrow(compiled_yrs) #107k

#Remove >95% water cells (assessed visually), cells outside climate data range
water_cells <- c(3667,3699,3726,3754,812,785,734,736,768,767,307)
water_list <- compiled_yrs$cell %in% water_cells
compiled_yrs <- compiled_yrs[!water_list,]

nrow(compiled_yrs) #106k

table(compiled_yrs$species)
table(compiled_yrs$species,compiled_yrs$year)

#Filter to breeding ranges - here we have individual files for species in the analysis
spec_list_br <- gsub(".rds","",list.files("data/spec_breed_cells/"))

compiled_yrs_trim<- c()
for (each_spec in unique(compiled_yrs$species)){
  #targ_spec <- "Rufous Hummingbird"
  targ_spec <- each_spec
  
  # Only read in species in the list 
  if (targ_spec %in% spec_list_br){
    #Read in breeding range data, limit to breeding range cells
    breed_rng_cls <- readRDS(paste("data/spec_breed_cells/",targ_spec,".rds",sep=""))
    
    spec_data <- filter(compiled_yrs,species==targ_spec) %>% select(c("species","cell","year","arr_GAM_mean","arr_GAM_sd"))
    
    spec_data <- filter(spec_data, cell %in% breed_rng_cls)
    
    #Remove cells with < 3 years of data
    cell_count <- spec_data %>% group_by(cell) %>% summarise(count=n()) %>% filter(count > 2)
    spec_data <- filter(spec_data,cell %in% cell_count$cell)
    
    if (nrow(spec_data) >= 1){
      compiled_yrs_trim <- rbind(compiled_yrs_trim,spec_data)
    }
  } else {
    print(targ_spec)
    print("dropped")
  }
  
}

compiled_yrs <- compiled_yrs_trim
compiled_yrs_trim <- NULL

nrow(compiled_yrs) #92k

#Limit to species with at least 4 cells represented
cl_count <- compiled_yrs %>% group_by(species) %>%
  summarise(num_cls = length(unique(cell)))

compiled_yrs <- merge(cl_count,compiled_yrs,by="species")
compiled_yrs <- filter(compiled_yrs,num_cls >= 4)

compiled_yrs <- select(compiled_yrs,-c("num_cls"))

nrow(compiled_yrs) #92k

saveRDS(compiled_yrs,"data/mdl_data/sp_yr_cl_dt.rds")

table(compiled_yrs$species,compiled_yrs$year)

#unique_cells <- unique(compiled_yrs$cell)
#write_csv(as.data.frame(unique_cells),"unique_cells_wb.csv")

