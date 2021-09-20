
library(tidyverse)
library(scales)
### plot type 1 error
input_dir <-  "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'

setwd(input_dir)
load(file='./estimated_sigma_Type_I_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0_sim_times_1000_random_seed_2021.RData')
library(latex2exp)
library(ggpubr)

df_p0_plot <- data.frame(naive = sort(unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                         hyun = sort(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                         union = sort(unlist(lapply(p_val_result,function(x)x[["Union"]]))))


df_p0_plot <- df_p0_plot %>% 
  mutate(theoretical = c(1:nrow(df_p0_plot))/nrow(df_p0_plot)) %>%
  pivot_longer(-theoretical, names_to = "p_type", values_to = "p_val") %>%
  mutate(p_type = factor(as.factor(p_type), levels=c("naive","hyun","union"))) %>%
  filter(p_type!="naive")

manual_color_pal <- scales::hue_pal()(3)

type_1_estimated_sigma <- ggplot(data = df_p0_plot) +  
  geom_point(aes(x = theoretical, y=p_val, colour=p_type),size=0.6) +
  ylab('P-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_manual('Tests', values = manual_color_pal[2:3], 
                     labels = unname(TeX(c('$p_{Hyun}$','$p_{\\hat{C}_{1},\\hat{C}_{2}}'))))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_7_a.png'),
    width = 6,height = 6, res=200,units='in')
type_1_estimated_sigma
dev.off()


alpha_thresh <- 0.05

sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("estimated_sigma_New_Power_1D_GFL*500*.*"),full.names = T)

concat_p_value_list <- list()

for (current_file in sim_files){
  load(file=current_file)
  concat_p_value_list <- c(concat_p_value_list,p_val_result)
  rm(p_val_result)
}

hyun_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Hyun"]]))
union_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Union"]]))
effect_size <- abs(unlist(lapply(concat_p_value_list, function(x)x[["mu_diff"]])))
nu_norm <- unlist(lapply(concat_p_value_list, function(x)x[["sd"]]))


power_new_df <- data.frame(effect_size=effect_size,
                           hyun_power = as.numeric(hyun_p_value<=alpha_thresh), 
                           union_power = as.numeric(union_p_value<=alpha_thresh)) %>%
  pivot_longer(-effect_size, names_to = "p_type", values_to = "cond_power")


manual_color_pal <- hue_pal()(3)[2:3]

power_by_cut <- power_new_df %>%
  mutate(break_up_mu = cut_interval(effect_size, n=8,dig.lab=0)) %>%
  group_by(break_up_mu,p_type) %>%
  summarise(mean_power = mean(cond_power),
            sd_power = sqrt(mean(cond_power)*(1-mean(cond_power))/n()))


plot_power <- power_by_cut %>%
  ggplot(aes(x=break_up_mu,y=mean_power,color=p_type,group=p_type))+ 
  geom_point()+
  geom_pointrange(aes(ymin = mean_power+1.96*sd_power, ymax = mean_power-1.96*sd_power))+
  geom_line()+
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('$| \\nu^T\\beta | / \\sigma$'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  scale_color_manual( "Tests", values = c(manual_color_pal[1],manual_color_pal[2]),
                      labels = unname(TeX(c('$p_{Hyun}$',"$p_{\\hat{C}_{1},\\hat{C}_{2}}$"))))


png(paste0(plot_output_dir,'Figure_7_b.png'),
    width = 6,height = 6, res=200,units='in')
plot_power
dev.off()




