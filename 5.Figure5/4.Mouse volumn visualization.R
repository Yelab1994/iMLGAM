library(ggquickeda)
library(ggprism)
library(ggsci)
library(ggpubr)
library(ggplot2)
##Read mouse tumor volumn data
narrow_df=read.csv('mouse tumor volumn data.csv')
##Draw the line charts
colors=ggsci::pal_d3("category20")(20)
ggplot(data=narrow_df,
       aes(x=days,y=value,
           group=tumor))+
  scale_y_continuous(limits = c(0, 1000),breaks = seq(0,1000,200),expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 21),breaks = seq(0,21,3),expand = c(0, 0))+
  labs(x="days", y="Tumor volume (mm³)")+
  geom_smooth(aes(color=Group),method = "loess", se = FALSE,size=0.6) +
  theme_prism(base_size = 10)+
  facet_wrap("Group", scales = "free_y", ncol =4 )+
  scale_color_manual(values = colors)+theme(legend.position = "none")
