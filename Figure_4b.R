library(tidyverse)
library(latex2exp)
library(ggpubr)
library(scales)

### plot type 1 error
input_dir  <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'
# load type 1 data
setwd(input_dir)
load(file='./Type_1_2D_GFL_stop_criteria_K_grid_8_level_2_0_sim_times_1000_random_seed_2021.RData')

df_p0_plot <- data.frame(naive = sort(unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                         hyun = sort(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                         union = sort(unlist(lapply(p_val_result,function(x)x[["Union"]]))))

df_p0_plot <- df_p0_plot %>% 
  mutate(theoretical = c(1:nrow(df_p0_plot))/nrow(df_p0_plot)) %>%
  pivot_longer(-theoretical, names_to = "p_type", values_to = "p_val") %>%
  mutate(p_type = factor(as.factor(p_type), levels=c("naive","hyun","union")))

type_1_2d <- ggplot(data = df_p0_plot) +  
  geom_point(aes(x = theoretical, y=p_val, colour=p_type),size=0.6) +
  ylab('P-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_discrete('Tests',
                       labels = unname(TeX(c('$p_{Naive}$','$p_{Hyun}$','$p_{C_{1},C_{2}}'))))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_4_b.png'),
    width = 6,height = 6, res=200,units='in')

print(type_1_2d)

dev.off()






