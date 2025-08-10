library(ggquickeda)
library(ggprism)
library(ggsci)
library(ggpubr)
library(ggplot2)
##Read mouse survival data
data=read.csv('mouse survival data.csv')
##Draw survival curves
surv_object <- Surv(data$futime, data$fustat)  
fit <- survfit(surv_object ~ group, data = data)  
surv_diff <- survdiff(surv_object ~ group, data = data)  
p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)  
ggplot(data, aes(time = futime, status = fustat, color = factor(group))) + 
  geom_km() + 
  theme_classic(base_size = 18)+
  scale_color_manual(values = colors)+
  theme_prism(base_size = 10)+
  scale_y_continuous(limits = c(0, 1.1),breaks = seq(0,1.1,0.25),expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 100),breaks = seq(0,100,25),expand = c(0, 0))+
  labs(x = 'time(days)', y = 'Percent survival(%)') + 
  annotate("text", x = 25, y = 0.2,   
           label = paste("P =", format(p_value, digits = 2)), size = 4.5)+
  theme(legend.position = "")
