library(tidyverse)
library(scales)
library(latex2exp)


### plot type 1 error
input_dir  <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
setwd(input_dir)
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'
num_rej_hyun <- c()
num_rej_union <- c()
delta_list <- c()
alpha_thresh <- 0.05
num_sims <- 500


load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0.5_sim_times_500_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 0.5)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_1_sim_times_500_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 1)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_1.5_sim_times_501_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 1.5)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_2_sim_times_501_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 2)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_2.5_sim_times_500_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 2.5)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_3_sim_times_500_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 3)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_3.5_sim_times_501_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 3.5)

load(file='./cond_power_1D_GFL_middle_stop_criteria_K_grid_200_level_2_4_sim_times_501_random_seed_2021.RData')

num_rej_hyun <- c(num_rej_hyun,mean(unlist(lapply(p_val_result,function(x)x[["Hyun"]]))<=0.05))
num_rej_union <- c(num_rej_union,mean(unlist(lapply(p_val_result,function(x)x[["Union"]]))<=0.05))
delta_list <- c(delta_list, 4)

manual_color_pal <- hue_pal()(3)[2:3]

cond_power_curve <- data.frame(delta=delta_list, 
                               hyun = num_rej_hyun, 
                               union=num_rej_union) %>%
  pivot_longer(-delta, names_to = "p_type", values_to = "cond_power") %>%
  mutate(se_cond_power = sqrt(cond_power*(1-cond_power)/num_sims))

plot_1d_middle_power <- ggplot(cond_power_curve, aes(x=delta, y=cond_power, colour=p_type))+
  geom_point()+
  geom_line()+
  geom_pointrange(aes(ymin=cond_power-1.96*se_cond_power, ymax=cond_power+1.96*se_cond_power), width=.2,
                  position=position_dodge(0.0))+
  theme_bw()+
  ylab(TeX('Conditional power at $\\alpha$=$0.05$'))+
  xlab(TeX('Difference in means: $\\delta$'))+
  scale_y_continuous(limits = c(0,1))+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  scale_color_manual( "Tests", values = manual_color_pal, breaks = c("hyun","union"),
                      labels = unname(TeX(c('$p_{Hyun}$','$p_{\\hat{C}_{1},\\hat{C}_{2}}'))))

png(paste0(plot_output_dir,'Figure_7_b.png'),
    width = 6,height = 6, res=200,units='in')
plot_1d_middle_power
dev.off()













