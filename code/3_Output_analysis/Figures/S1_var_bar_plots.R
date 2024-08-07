library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(RColorBrewer)

var_exp_data <- read_csv("data/output/var_exp/var_exp_sp.csv")

var_exp_by_group <- var_exp_data %>% 
  group_by(Group) %>% 
  summarise(mean_rl_var = mean(rel_var),
            mean_tmin = mean(tmin),
            mean_swe = mean(swe),
            mean_prcp = mean(prcp),
            mean_enso = mean(enso),
            mean_nao = mean(nao),
            mean_exp = mean(total_exp),
            unexp_var = (1-mean(total_exp)),
            count = n()) %>% 
  filter(count >= 5) %>% arrange(mean_exp)

var_exp_by_group$Group[which(var_exp_by_group$Group == "Tyrant Flycatchers: Pewees, Kingbirds, and Allies")] <- "Tyrant Flycatchers" 
rs_var_exp <- reshape2::melt(var_exp_by_group[,c(1,3:7)],value.name = "var_exp")
rs_var_exp$var_exp <- rs_var_exp$var_exp * 100
rs_var_exp$variable_grouped <- factor(rs_var_exp$variable, levels = rev(c("mean_tmin","mean_swe","mean_prcp","mean_enso","mean_nao")))
rs_var_exp$group_grouped <- factor(rs_var_exp$Group, levels = (unique(rs_var_exp$Group)))

# Stacked
ggplot(rs_var_exp, aes(fill=variable_grouped, y=group_grouped, x=var_exp)) + 
  geom_vline(xintercept = 5,lty=2,col="darkgrey") +
  geom_vline(xintercept = 10,lty=2,col="darkgrey") +
  geom_vline(xintercept = 15,lty=2,col="darkgrey") +
  geom_vline(xintercept = 20,lty=2,col="darkgrey") +
  geom_bar(position="dodge", stat="identity",color="black") +
  scale_fill_manual(labels = rev(c("Temp.", "Snow","Precip.","ENSO","NAO")),
                    values=rev(c("firebrick3", "skyblue2", "forestgreen","orchid3","tan3"))) +
  guides(fill = guide_legend(title=element_blank(),reverse = TRUE)) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        legend.position = c(0.85, 0.2),
        legend.key.size = unit(.85, 'cm'),
        legend.background = element_rect(fill="white",
        size=0.5, linetype="solid", colour ="black")) +
  labs(x=("Variance Explained (%)"),y=(""))

