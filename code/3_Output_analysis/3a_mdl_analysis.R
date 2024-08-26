#Get effect sizes, cl, sp
library(MCMCvis)
library(dplyr)

model_in <- readRDS("data/output/PhenoClimate_mdl.rds")
data_in <- readRDS("data/mdl_data/ms_data1.rds")
d2 <- readRDS("data/mdl_data/ms_data2.rds")

MCMCdiag(model_in,
         "PhenoClimate_mdl_diag.txt",
         summary=FALSE
         )
as <- MCMCsummary(model_in,"alpha")
b1s <- MCMCsummary(model_in,"beta1")
b2s <- MCMCsummary(model_in,"beta2")
b3s <- MCMCsummary(model_in,"beta3")
p1s <- MCMCsummary(model_in,"phi1")
p2s <- MCMCsummary(model_in,"phi2")

#Save effect sizes, mean arrival dates by species/cell
sp_cl_df <- d2 %>% select(species,cell,cl_sp,cell_lat,snow_p) %>% unique()

sp_cl_df <- data.frame(sp=sp_cl_df$species,
                       cell=sp_cl_df$cell,
                       cl_sp=sp_cl_df$cl_sp,
                       cell_lat=sp_cl_df$cell_lat,
                       snow_p=sp_cl_df$snow_p,
                       as_median = as$`50%`,
                       b1s_median = b1s$`50%`,
                       b2s_median = b2s$`50%`,
                       b3s_median = b3s$`50%`,
                       p1s_median = p1s$`50%`,
                       p2s_median = p2s$`50%`)
sp_cl_df$as_og <- sp_cl_df$as_median + mean(d2$arr_GAM_mean)

saveRDS(sp_cl_df,"sp_cl_mdl_params.rds")

