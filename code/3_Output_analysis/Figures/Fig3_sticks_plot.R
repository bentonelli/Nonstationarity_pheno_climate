# Effect over observed temp niche
library(MCMCvis)
library(dplyr)
library(readr)
library(ggplot2)

model_in <- readRDS("data/output/PhenoClimate_mdl.rds")
data_in <- readRDS("data/mdl_data/ms_data1.rds")
d2 <- readRDS("data/mdl_data/ms_data2.rds")
spec_names <- read_csv("data/mdl_data/spec_groups.csv")

mb1 <- MCMCsummary(model_in,"mu_beta1")
mb2 <- MCMCsummary(model_in,"mu_beta2")
mb3 <- MCMCsummary(model_in,"mu_beta3")

t1 <- MCMCsummary(model_in,"theta1")
t2 <- MCMCsummary(model_in,"theta2")
t3 <- MCMCsummary(model_in,"theta3")

sp_cl_df <- readRDS("data/output/sp_cl_mdl_params.rds")
d2_eff <- merge(d2,sp_cl_df) %>% select(-c(sp))

head(d2_eff)

par(mfrow=c(2,1))
par(mar=c(5,1,1,.5))
par(las=1)
par(pty="s")
plot(NULL,xlim=c(-20,25),ylim=c(-4,2),xlab="Mean Arrival Temperature (°C)",ylab="Effect +1SD Temp.",
     cex.lab=1.5,cex.axis=1.25)
abline(h=0,lty=2)
for (each_spec in 1:222){
  spec_name <- unique(d2$species)[each_spec]
  sp_d2_eff <- d2_eff %>% filter(species==spec_name) %>% select(cl_sp,mean_tmin_29,tmin_avg_std) %>% unique()
  min_x <- min(sp_d2_eff$tmin_avg_std)
  max_x <- max(sp_d2_eff$tmin_avg_std)
  xx_range <- seq(min_x,max_x,by=.1)*sd(sp_d2_eff$mean_tmin_29) + mean(sp_d2_eff$mean_tmin_29)
  calc_y <- mb1$`50%`[each_spec] + t1$`50%`[each_spec] * seq(min_x,max_x,by=.1)
  if (spec_name == "Western Kingbird"){
    xx_range_sp_interest <- xx_range
    calc_y_sp_interest <- calc_y
  } else {
    points(xx_range,calc_y,pch=19,col=alpha("firebrick4",.2),type="l",lwd=2)
  }
}
points(xx_range_sp_interest,calc_y_sp_interest,pch=19,col=alpha("firebrick4",1),type="l",lwd=5)

#plot(NULL,xlim=c(.1,800),ylim=c(-3,3),xlab="Mean Arrival Snow",ylab="Effect +1SD Snow",
     #cex.lab=1.75,cex.axis=1.75, log='x')
#abline(h=0,lty=2)
for (each_spec in 1:222){
  spec_name <- unique(d2$species)[each_spec]
  sp_d2_eff <- d2_eff %>% filter(species==spec_name) %>% select(cl_sp,mean_swe_29,swe_avg_std) %>% unique()
  min_x <- min(sp_d2_eff$swe_avg_std)
  max_x <- max(sp_d2_eff$swe_avg_std)
  xx_range <- seq(min_x,max_x,by=.1)*sd(sp_d2_eff$mean_swe_29) + mean(sp_d2_eff$mean_swe_29)
  calc_y <- mb2$`50%`[each_spec] + t2$`50%`[each_spec] * seq(min_x,max_x,by=.1)
  if (spec_name == "Black-headed Grosbeak"){
    xx_range_sp_interest <- xx_range
    calc_y_sp_interest <- calc_y
  } else {
    #points(xx_range,calc_y,pch=19,col=alpha("skyblue3",.25),type="l",lwd=2)
  }
}
#points(xx_range_sp_interest,calc_y_sp_interest,pch=19,col=alpha("grey30",1),type="l",lwd=3)

plot(NULL,xlim=c(0,9),ylim=c(-1,1.5),xlab="Mean Arrival Precip. (mm/day)",ylab="Effect +1SD Precipitation",
     cex.lab=1.5,cex.axis=1.25)
abline(h=0,lty=2)
for (each_spec in 1:222){
  spec_name <- unique(d2$species)[each_spec]
  sp_d2_eff <- d2_eff %>% filter(species==spec_name) %>% select(cl_sp,mean_prcp_29,prcp_avg_std) %>% unique()
  min_x <- min(sp_d2_eff$prcp_avg_std)
  max_x <- max(sp_d2_eff$prcp_avg_std)
  xx_range <- seq(min_x,max_x,by=.1)*sd(sp_d2_eff$mean_prcp_29) + mean(sp_d2_eff$mean_prcp_29)
  calc_y <- mb3$`50%`[each_spec] + t3$`50%`[each_spec] * seq(min_x,max_x,by=.1)
  if (spec_name == "Greater Yellowlegs"){
    xx_range_sp_interest <- xx_range
    calc_y_sp_interest <- calc_y
  } else {
    points(xx_range,calc_y,pch=19,col=alpha("darkgreen",.2),type="l",lwd=2)
  }
}
points(xx_range_sp_interest,calc_y_sp_interest,pch=19,col=alpha("darkgreen",1),type="l",lwd=5)


