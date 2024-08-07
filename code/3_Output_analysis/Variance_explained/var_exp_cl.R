#Get at predicted arrival for particular species, using both models
library(MCMCvis)
library(dplyr)

model_in <- readRDS("data/output/PhenoClimate_mdl.rds")
data_in <- readRDS("data/mdl_data/ms_data1.rds")
d2 <- readRDS("data/mdl_data/ms_data2.rds")

arr_mean_cross_sp <- mean(d2$arr_GAM_mean)
spec_names <- unique(d2$species)
cell_names <- unique(d2$cell)

snow_p <- d2$snow_p

yy_trues <- MCMCchains(model_in,params="yy_true")

ma1 <- MCMCchains(model_in,params="alpha")
mb1 <- MCMCchains(model_in,params="beta1")
mb2 <- MCMCchains(model_in,params="beta2")
mb3 <- MCMCchains(model_in,params="beta3")
ph1 <- MCMCchains(model_in,params="phi1")
ph2 <- MCMCchains(model_in,params="phi2")


### Variance explained, by species ####
#Pick a species, get chains
#cl_ind <- which(d2$species == "Eastern Phoebe")
var_exp_cm_c <- c()
for(each_cell in unique(d2$cell)){
  print(each_cell)
  cl_ind <- which(d2$cell == each_cell)
  
  unq_sp <- length(unique(d2$species[cl_ind]))
  samp_ind <- sample(nrow(yy_trues),200,replace=FALSE)
  c_var_exp_temp <- c()
  c_var_exp_prcp <- c()
  c_var_exp_swe <- c()
  c_var_exp_enso <- c()
  c_var_exp_nao <- c()
  c_var_exp_all <- c()
  c_var_all <- c()
  c_arr_mean_all <- c()
  for (nn in samp_ind){
    #print(nn)
    
    # TMIN
    grep_temp <-  as.numeric(mb1[nn,d2$cl_sp[cl_ind]] * d2$anom_tmin_c29[cl_ind])
    e_temp <- (yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]]) - grep_temp
    var_exp_temp <- var(grep_temp)/(var(grep_temp) + var(e_temp))
    c_var_exp_temp <- c(c_var_exp_temp,var_exp_temp)
    
    # SWE
    grep_swe <-  as.numeric(mb2[nn,d2$cl_sp[cl_ind]] * snow_p[cl_ind] * d2$anom_swe_c29[cl_ind])
    e_swe <- (yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]]) - grep_swe
    var_exp_swe <- var(grep_swe)/(var(grep_swe) + var(e_swe))
    c_var_exp_swe <- c(c_var_exp_swe,var_exp_swe)
    
    # Precip
    grep_prcp <-  as.numeric(mb3[nn,d2$cl_sp[cl_ind]] * d2$anom_prcp_c29[cl_ind])
    e_prcp <- (yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]]) - grep_prcp
    var_exp_prcp <- var(grep_prcp)/(var(grep_prcp) + var(e_prcp))
    c_var_exp_prcp <- c(c_var_exp_prcp,var_exp_prcp)
    
    # ENSO
    grep_enso <-  as.numeric(ph1[nn,d2$cl_sp[cl_ind]] * d2$ENSO[cl_ind])
    e_enso <- (yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]]) - grep_enso
    var_exp_enso <- var(grep_enso)/(var(grep_enso) + var(e_enso))
    c_var_exp_enso <- c(c_var_exp_enso,var_exp_enso)
    
    #NAO
    grep_nao <-  as.numeric(ph2[nn,d2$cl_sp[cl_ind]] * d2$NAO[cl_ind])
    e_nao <- (yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]]) - grep_nao
    var_exp_nao <- var(grep_nao)/(var(grep_nao) + var(e_nao))
    c_var_exp_nao <- c(c_var_exp_nao,var_exp_nao)
    
    # EVERYTHING
    grep_all <- as.numeric(mb1[nn,d2$cl_sp[cl_ind]] * d2$anom_tmin_c29[cl_ind] + 
                             mb2[nn,d2$cl_sp[cl_ind]] * d2$anom_swe_c29[cl_ind] + 
                             mb3[nn,d2$cl_sp[cl_ind]] * d2$anom_prcp_c29[cl_ind] +
                             ph1[nn,d2$cl_sp[cl_ind]] * d2$ENSO[cl_ind] +
                             ph2[nn,d2$cl_sp[cl_ind]] * d2$NAO[cl_ind])
    e_all <- (yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]]) - grep_all
    var_exp_all <- var(grep_all)/(var(grep_all) + var(e_all))
    c_var_exp_all <- c(c_var_exp_all,var_exp_all)
    
    # Variance, explained or otherwise
    var_all <- var(yy_trues[nn,cl_ind] - ma1[nn,d2$cl_sp[cl_ind]])
    c_var_all <- c(c_var_all,var_all)
    
    # Mean arrival across all cells, years
    arr_mean_all <- mean(yy_trues[nn,cl_ind]+arr_mean_cross_sp)
    c_arr_mean_all <- c(c_arr_mean_all,arr_mean_all)
    
    var_exp <- data.frame(tmin=c_var_exp_temp,
                          prcp=c_var_exp_prcp,
                          swe=c_var_exp_swe,
                          enso = c_var_exp_enso,
                          nao = c_var_exp_nao,
                          total_exp=c_var_exp_all,
                          rel_var = c_var_all,
                          mean_arr_true = c_arr_mean_all,
                          nn_count = unq_sp)
  }
  
  var_exp_cm <- colMeans(var_exp)
  var_exp_cm_c <- rbind(var_exp_cm_c,var_exp_cm)
}

var_exp_cm_c <- as.data.frame(var_exp_cm_c)
var_exp_cm_c$cell <- unique(d2$cell)

write.csv(var_exp_cm_c,"data/output/var_exp/var_exp_cell.csv")

