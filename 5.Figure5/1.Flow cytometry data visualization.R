library(ggprism)
library(ggsci)
library(ggpubr)
library(ggplot2)
###Read flow cytometry data
dat <- read.csv("flow cytometry data.csv",header = T,row.names = 1)
###Draw the bar chart.
mycolors <-ggsci::pal_d3("category20")(4)
group=levels(factor(dat$group))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
p1 =ggplot(data = dat,   
       aes(x = factor(group), y = exp, fill = group)) +  
  stat_summary(fun = "mean",  
               geom = "bar",  
               width = 0.6,
               fill='white',
               colour = mycolors,  
               size = 0.9) +  
  geom_jitter(size = 3, 
              aes(colour = group, shape = group),  
              fill = "white",  
              stroke = 1,  
              width = 0.2) +  
  scale_shape_manual(values = rep(21, length(unique(dat$group)))) +  
  scale_y_continuous(limits = c(0, 42), expand = c(0, 0)) +  
  stat_summary(fun = "mean",  
               fun.max = function(x) mean(x) + sd(x),  
               fun.min = function(x) mean(x) - sd(x),  
               geom = "errorbar",
               color=mycolors,
               width = 0.4,  
               size = 0.9) +  
  theme_prism(axis_text_angle = 45) +  
  theme(legend.position = '',
        axis.title.y = element_text(size = 13),  
        axis.text.y  = element_text(size = 11.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +  
  scale_fill_manual(values = mycolors) +  
  scale_color_manual(values = mycolors) +  
  labs(x = element_blank(), y = "CD8+ cells(% of CD3+ T cells)", title = '') 

color_value <- c(mycolors[2:4], mycolors[3:4], mycolors[4])

x_value <- rep(1:3, 3:1)
y_value_group1 <- rep(max(dat$exp[dat$group==group1 ]),3)+0.1
y_value_group2 <- rep(max(dat$exp[dat$group==group2 ]),2)+0.1
y_value_group3 <- rep(max(dat$exp[dat$group==group3 ]),1)+0.1

y_value <- c(y_value_group1,y_value_group2,y_value_group3)



y_value <- y_value + c(0.1*1:3, 0.1*1:2, 0.1*1:1)
color_value <- c(mycolors [2:4], mycolors [3:4], mycolors [4])



for (i in 1:nrow(stat.test)) {
  if (stat.test$p.adj.signif[i] != "ns") {
    y_tmp <- y_value[i]
    p1 <- p1+annotate(geom = "text",
                      label = stat.test$p.adj.signif[i],
                      x = x_value[i],
                      y = y_tmp,
                      color = color_value[i])
  }
}
p1
