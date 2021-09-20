library(tidyverse)
library(latex2exp)
library(scales)
library(intervals)
manual_color_pal <- scales::hue_pal()(3)[2:3]
### plot power 2D
input_dir  <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'
alpha_thresh <- 0.05
# load power data
setwd(input_dir)
sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("*2D*500*.*"),full.names = T)

concat_p_value_list <- list()
concat_rand_list <- c()

for (current_file in sim_files){
  load(file=current_file)
  concat_p_value_list <- c(concat_p_value_list,p_val_result)
  concat_rand_list <- c(concat_rand_list,rand_list)
  rm(p_val_result)
}

concat_rand_list <- unlist(concat_rand_list)
hyun_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Hyun"]]))
union_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Union"]]))
effect_size <- abs(unlist(lapply(concat_p_value_list, function(x)x[["mu_diff"]])))
nu_norm <- unlist(lapply(concat_p_value_list, function(x)x[["sd"]]))

power_new_df <- data.frame(effect_size=effect_size,
                           rand_index = concat_rand_list,
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value_rerun)) %>%
  pivot_longer(-c(rand_index,effect_size), names_to = "p_type", 
               values_to = "cond_power") %>%
  mutate(cond_power = cond_power<=alpha_thresh)

power_by_cut <- power_new_df %>%
  mutate(break_up_mu = cut_interval(effect_size, n=8)) %>%
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
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=12,hjust = 0),
        axis.title=element_text(size=15))+
  scale_color_manual( "Tests", values = c(manual_color_pal[1],manual_color_pal[2]),
                      labels = unname(TeX(c('$p_{Hyun}$',"$p_{C_{1},C_{2}}$"))))

png(paste0(plot_output_dir,'Figure_4_c.png'),
    width = 7,height = 7, res=200,units='in')
plot_power
dev.off()






