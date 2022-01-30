library(tidyverse)
library(scales)
library(latex2exp)
library(ggpubr)

alpha_thresh <- 0.05
### plot type 1 error
input_dir <-  "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'

setwd(input_dir)
load(file='./estimated_sigma_MAD_Type_I_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0_sim_times_1000_random_seed_2021.RData')

df_p0_plot <- data.frame(naive = (unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                         hyun = (unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                         union = (unlist(lapply(p_val_result,function(x)x[["Union"]]))))

manual_color_pal <- scales::hue_pal()(3)[2:3]

p_sample_type1 <- df_p0_plot %>% 
  mutate(theoretical = c(1:nrow(df_p0_plot))/nrow(df_p0_plot)) %>%
  pivot_longer(-theoretical, names_to = "p_type", values_to = "p_val") %>%
  mutate(p_type = factor(as.factor(p_type), levels=c("naive","hyun","union"))) %>%
  filter(p_type!="naive") %>%
  group_by(p_type) %>%
  arrange(p_val) %>%
  mutate(theoretical = ecdf(p_val)(p_val)) %>%
  ggplot(aes(y = p_val, x = theoretical,color = p_type)) +
  geom_point(size=0.6) +
  ylab('Selective p-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  #ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_manual('Tests', 
                     values = manual_color_pal, 
                     labels = unname(TeX(c('$p_{Hyun}(\\hat{\\sigma}_{Sample}$)','$p_{\\hat{C}_{1},\\hat{C}_{2}}(\\hat{\\sigma}_{Sample}$)'))))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_13_c.png'),
    width = 7,height = 7, res=200,units='in')
p_sample_type1
dev.off()



