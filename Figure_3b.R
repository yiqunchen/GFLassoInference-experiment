# load libraries for data processing and plotting
library(tidyverse)
library(scales)
library(latex2exp)
library(ggpubr)

### plot type 1 error
input_dir <-  "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/"

setwd(input_dir)
load(file='./Type_I_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0_sim_times_1000_random_seed_2021.RData')

df_p0_plot <- data.frame(naive = sort(unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                         hyun = sort(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                         union = sort(unlist(lapply(p_val_result,function(x)x[["Union"]]))))


df_p0_plot <- df_p0_plot %>% 
  mutate(theoretical = c(1:nrow(df_p0_plot))/nrow(df_p0_plot)) %>%
  pivot_longer(-theoretical, names_to = "p_type", values_to = "p_val") %>%
  mutate(p_type = factor(as.factor(p_type), levels=c("naive","hyun","union")))


type_1 <- ggplot(data = df_p0_plot) +  
  geom_point(aes(x = theoretical, y=p_val, colour=p_type),size=0.6) +
  ylab('P-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_discrete('Tests',
                       labels = unname(TeX(c('$p_{Naive}$','$p_{Hyun}$',
                                             '$p_{\\hat{C}_{1},\\hat{C}_{2}}'))))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_3_b.png'),
    width = 7,height = 7, res=200,units='in')
type_1
dev.off()




