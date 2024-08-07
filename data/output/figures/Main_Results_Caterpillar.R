library(MCMCvis)
library(ggplot2)
library(ggplot2)
model_in <- readRDS("PhenoClimate_mdl.rds")

sum_data <- MCMCsummary(model_in,
         params=c("eta1","eta2","eta3",
                  "nu1","nu2","nu3",
                  "omega1","omega2",
                  "kappa1","kappa2"),
         probs = c(0.055,0.5,0.945),
         pg0 = TRUE
         )
saveRDS(sum_data,"sum_data.rds")

MCMCplot(model_in,
         params=c("eta1","eta2","eta3",
                             "nu1","nu2","nu3",
                             "omega1","omega2",
                             "kappa1","kappa2"),
        labels = c("Global Temp","Global SWE","Global Precip",
                   "Niche Temp","Niche SWE","Niche Precip",
                   "Global ENSO","Global NAO",
                   "Winter ENSO","Winter NAO"),
        ci = c(50, 89),
        xlab="", sz_med = 1.5,sz_thick=4,sz_thin = 2)

mm_sum <- MCMCsummary(model_in,
         params=c("eta1","eta2","eta3",
                  "nu1","nu2","nu3",
                  "omega1","omega2",
                  "kappa1","kappa2"),
         probs = c(0.055, 0.5, 0.945),
         hpd_prob=.89,
         pg0 = TRUE)

mm_beta1 <- MCMCsummary(model_in,
                        params=c("mu_beta1"),
                        probs = c(0.055, 0.5, 0.945),
                        hpd_prob=.89)

mm_beta2 <- MCMCsummary(model_in,
                        params=c("mu_beta2"),
                        probs = c(0.055, 0.5, 0.945),
                        hpd_prob=.89)

mm_beta3 <- MCMCsummary(model_in,
                        params=c("mu_beta3"),
                        probs = c(0.055, 0.5, 0.945),
                        hpd_prob=.89)

mm_theta1 <- MCMCsummary(model_in,
                        params=c("theta1"),
                        probs = c(0.055, 0.5, 0.945),
                        hpd_prob=.89)

mm_theta2 <- MCMCsummary(model_in,
                         params=c("theta2"),
                         probs = c(0.055, 0.5, 0.945),
                         hpd_prob=.89)

mm_theta3 <- MCMCsummary(model_in,
                         params=c("theta3"),
                         probs = c(0.055, 0.5, 0.945),
                         hpd_prob=.89)

mm_phi1 <- MCMCsummary(model_in,
                         params=c("mu_phi1"),
                         probs = c(0.055, 0.5, 0.945),
                         hpd_prob=.89)

mm_phi2 <- MCMCsummary(model_in,
                         params=c("mu_phi2"),
                         probs = c(0.055, 0.5, 0.945),
                         hpd_prob=.89)



disp_factor <- 1
plot(NULL,xlim=c(-4,4),ylim=c(0,11))
abline(v=0,lwd=.5,col="grey",lty=2)

points(mm_beta1$`50%`,rep(10,nrow(mm_beta1))+rnorm(nrow(mm_beta1),0,disp_factor),cex=.75,col=alpha("firebrick3",.1),pch=19)
points(mm_beta2$`50%`,rep(9,nrow(mm_beta2))+rnorm(nrow(mm_beta2),0,disp_factor),cex=.75,col=alpha("tan3",.1),pch=19)
points(mm_beta3$`50%`,rep(8,nrow(mm_beta3))+rnorm(nrow(mm_beta3),0,disp_factor),cex=.75,col=alpha("skyblue3",.1),pch=19)

points(mm_theta1$`50%`,rep(7,nrow(mm_theta1))+rnorm(nrow(mm_theta1),0,disp_factor),cex=.75,col=alpha("firebrick3",.1),pch=19)
points(mm_theta2$`50%`,rep(6,nrow(mm_theta2))+rnorm(nrow(mm_theta2),0,disp_factor),cex=.75,col=alpha("tan3",.1),pch=19)
points(mm_theta3$`50%`,rep(5,nrow(mm_theta3))+rnorm(nrow(mm_theta3),0,disp_factor),cex=.75,col=alpha("skyblue3",.1),pch=19)

points(mm_phi1$`50%`,rep(4,nrow(mm_phi1))+rnorm(nrow(mm_phi1),0,disp_factor),cex=.75,col=alpha("chocolate4",.1),pch=19)
points(mm_phi2$`50%`,rep(3,nrow(mm_phi1))+rnorm(nrow(mm_phi1),0,disp_factor),cex=.75,col=alpha("darkolivegreen",.1),pch=19)

points(mm_sum$`50%`,nrow(mm_sum):1,cex=2.5,ylim=c(0,(nrow(mm_sum)+1)),pch="|",xlim=c(-4,4),col=alpha("black",.9))

for (nn in 1:nrow(mm_sum)){
  points(c(mm_sum$`5.5%`[nn],mm_sum$`94.5%`[nn]),1+c(nrow(mm_sum)-nn,nrow(mm_sum)-nn),type="l",col="grey50",lwd=4)
}


### simpler, split plot ####
mm_sum_simple <- MCMCsummary(model_in,
                      params=c("eta1","eta2","eta3",
                               "omega1","omega2"),
                      probs = c(0.055, 0.5, 0.945),
                      pg0 = TRUE,
                      hpd_prob=.89)

# Version 1
par(pty="s")
par(mar=c(8,12,2,2))
disp_factor1 <- -.4
disp_factor2 <- .4
trans_fact <- .6
sp_cex <- 1.25
plot(NULL,xlim=c(-3.5,3.5),ylim=c(0.5,5.5),xlab="Effect on Arrival Phenology (days)",ylab="",yaxt="n",
     cex.lab=2,cex.axis = 2)
abline(v=0,lwd=3,col="grey70",lty=2)

points(mm_beta1$`50%`,rep(5,nrow(mm_beta1))+runif(nrow(mm_beta1),disp_factor1,disp_factor2),cex=sp_cex,bg=alpha("firebrick4",trans_fact),col="transparent",pch=21)
points(mm_beta2$`50%`,rep(4,nrow(mm_beta2))+runif(nrow(mm_beta2),disp_factor1,disp_factor2),cex=sp_cex,bg=alpha("skyblue4",trans_fact),col="transparent",pch=21)
points(mm_beta3$`50%`,rep(3,nrow(mm_beta3))+runif(nrow(mm_beta3),disp_factor1,disp_factor2),cex=sp_cex,bg=alpha("darkgreen",trans_fact),col="transparent",pch=21)

points(mm_phi1$`50%`,rep(2,nrow(mm_phi1))+runif(nrow(mm_phi1),disp_factor1,disp_factor2),cex=sp_cex,bg=alpha("orchid2",trans_fact),col="transparent",pch=21)
points(mm_phi2$`50%`,rep(1,nrow(mm_phi2))+runif(nrow(mm_phi2),disp_factor1,disp_factor2),cex=sp_cex,bg=alpha("chocolate3",trans_fact),col="transparent",pch=21)

for (nn in 1:nrow(mm_sum_simple)){
  if (mm_sum_simple$`p>0`[nn] >= .89 | mm_sum_simple$`p>0`[nn] <= .11){
    points(mm_sum_simple$`50%`[nn],(nrow(mm_sum_simple):1)[nn],cex=4,pch="|",col=alpha("black",1))
    #points(c(mm_sum_simple$`5.5%`[nn],mm_sum_simple$`94.5%`[nn]),1+c(nrow(mm_sum_simple)-nn,nrow(mm_sum_simple)-nn),type="l",col="black",lwd=4)
  } else {
    points(mm_sum_simple$`50%`[nn],(nrow(mm_sum_simple):1)[nn],cex=4,pch="|",col=alpha("black",1))
    #points(c(mm_sum_simple$`5.5%`[nn],mm_sum_simple$`94.5%`[nn]),1+c(nrow(mm_sum_simple)-nn,nrow(mm_sum_simple)-nn),type="l",col="grey40",lwd=4)
  }
  
}
ytick<-seq(1, 5, by=1)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1], ytick,  
     labels = rev(c("Temperature","Snowpack","Precipitation","ENSO","NAO")), srt = 0, pos = 2, xpd = TRUE,cex=2)


#text(labels=c("Temp.","Snow","Precip.","ENSO","NAO"),x=c(3,3,3,3,3),y=rev(c(1,2,3,4,5)),cex=2)
