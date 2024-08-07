library(readr)
library(dplyr)
library(ggplot2)
library(viridis)

var_exp <- read_csv("data/output/var_exp/var_exp_sp.csv")

library(pals)
#col_list <- rev(parula(n=nrow(spec_names)))
var_exp$arr_round <- round(var_exp$mean_arr_true,0)
arr_dates <- unique(var_exp$arr_round )
min_save <- min(arr_dates)
arr_dates <- arr_dates - min(arr_dates) + 1
col_list <- rev(parula(n=max(arr_dates)))

var_exp$arr_round_col <- var_exp$arr_round - min_save + 1            
point_size <- 1.25
point_alpha <- .6
par(pty="s")
par(mar=c(5,5,5,5))
sc_y <- rnorm(nrow(var_exp),0,.15) 
plot(var_exp$total_exp*100,rep(6,nrow(var_exp))+sc_y,xlim=c(.15,100),ylim=c(.25,6.75),xaxt="n",
     col=alpha(col_list[var_exp$arr_round_col],point_alpha),pch=19,log="x",xlab="Variance Explained (%)",yaxt="n",ylab="",cex=point_size,
     cex.lab=1.5,cex.axis=1.5)
axis(1, at = c(.2,.5,1,2,5,10,20,50,100),las=1,cex.axis=1.5,
     cex.lab=1.5)
axis(2, at = c(6,5,4,3,2,1),
     cex.axis=1.5,
     cex.lab=1.5,
     labels = c("All","Temp","Snow","Precip.","ENSO","NAO"),las=2)
abline(v=.5,lty=2,lwd=2.5,col="grey20")
abline(v=5,lty=2,lwd=2.5,col="grey20")
abline(v=50,lty=2,lwd=2.5,col="grey20")
points(var_exp$tmin*100,rep(5,nrow(var_exp))+sc_y,col=alpha(col_list[var_exp$arr_round_col],point_alpha),pch=19,cex=point_size)
points(var_exp$swe*100,rep(4,nrow(var_exp))+sc_y,col=alpha(col_list[var_exp$arr_round_col],point_alpha),pch=19,cex=point_size)
points(var_exp$prcp*100,rep(3,nrow(var_exp))+sc_y,col=alpha(col_list[var_exp$arr_round_col],point_alpha),pch=19,cex=point_size)
points(var_exp$enso*100,rep(2,nrow(var_exp))+sc_y,col=alpha(col_list[var_exp$arr_round_col],point_alpha),pch=19,cex=point_size)
points(var_exp$nao*100,rep(1,nrow(var_exp))+sc_y,col=alpha(col_list[var_exp$arr_round_col],point_alpha),pch=19,cex=point_size)

legend_image <- rev(as.raster(matrix(col_list, ncol=1)))
plot(c(0,2),c(0,max(var_exp$arr_round_col)),type = 'n', axes = F,xlab = '', ylab = '')
rasterImage(legend_image, 0, 0, 1,max(var_exp$arr_round_col))
text(x=1.25, y = seq(min(var_exp$arr_round_col),max(var_exp$arr_round_col),by=10), 
     labels = seq(min(var_exp$arr_round),max(var_exp$arr_round),by=10),cex=2)

### ggplot version ####

#var_exp$swe[which(var_exp$swe==0)] <- .00005
ggplot() + 
  #ylim(-.5,5.5) + 
  scale_y_continuous(breaks=c(-0.5,0,1,2,3,4,5,5.5),labels=c("","NAO","ENSO","Precipitation","Snowpack","Temperature","All","")) +
  scale_x_continuous(trans='pseudo_log',breaks=c(.1,1,5,10,20,50,100)) +
  geom_vline(xintercept=c(.1,1,5,10,20,50,100),lty=2,col="grey") +
  xlab("Variance Explained (%)") +
  ylab("") +
  geom_point(data=var_exp, aes(x=total_exp*100, y=5+runif(nrow(var_exp),-.25,.25),color=arr_round),size=2) +
  geom_point(data=var_exp, aes(x=tmin*100, y=4+runif(nrow(var_exp),-.25,.25),color=arr_round),size=2) +
  geom_point(data=var_exp, aes(x=swe*100, y=3+runif(nrow(var_exp),-.25,.25),color=arr_round),size=2) +
  geom_point(data=var_exp, aes(x=prcp*100, y=2+runif(nrow(var_exp),-.25,.25),color=arr_round),size=2) +
  geom_point(data=var_exp, aes(x=enso*100, y=1+runif(nrow(var_exp),-.25,.25),color=arr_round),size=2) +
  geom_point(data=var_exp, aes(x=nao*100, y=0+runif(nrow(var_exp),-.25,.25),color=arr_round),size=2) +
  scale_color_viridis_c(option="viridis",direction = 1,begin=.1,end=1,name="Arrival Date",alpha=.9) + 
  theme(legend.position = c(.85,.2)) +
  guides(color = guide_colourbar(title.position="top", title.hjust = 0.5)) + 
  theme_classic(base_size = 18) +
  theme(legend.position = c(.17,.9), 
        legend.direction = "horizontal",
        legend.key.width = unit(.85,"cm"),
        legend.text = element_text(size=11),
        legend.background = element_rect(fill=alpha("white",.5),
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        aspect.ratio = 1)




