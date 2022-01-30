library(tidyverse)
library(scales)
library(latex2exp)
library(ggpubr)

alpha_thresh <- 0.05
### plot type 1 error
input_dir <-  "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/"

setwd(input_dir)

sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("K_init_4*New_Power*_1500_*_2000.RData"),full.names = T)


concat_p_value_list <- list()
concat_p_value_list <- list()
concat_rand_list <- c()
concat_delta_list <- c()
concat_sigma_list <- c()

for (current_file in sim_files){
  load(file=current_file)
  # extract info 
  sigma_strsplit <- as.numeric(strsplit(strsplit(current_file,"sigma_")[[1]][2],
                                        "_sim_times")[[1]][1])
  delta_strsplit <- as.numeric(strsplit(strsplit(current_file,"level_2_")[[1]][2],
                                        "_sigma")[[1]][1])
  concat_p_value_list <- c(concat_p_value_list,p_val_result)
  concat_rand_list <- c(concat_rand_list, rand_list)
  concat_delta_list <- c(concat_delta_list,
                         rep((delta_strsplit), times = length(p_val_result)))
  concat_sigma_list <- c(concat_sigma_list,
                         (rep((sigma_strsplit), times = length(p_val_result))))
  
  rm(p_val_result)
  rm(rand_list)
}

naive_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Naive"]]))
hyun_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Hyun"]]))
union_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Union"]]))
effect_size <- abs(unlist(lapply(concat_p_value_list, function(x)x[["mu_diff"]])))
nu_norm <- unlist(lapply(concat_p_value_list, function(x)x[["sd"]]))
vdelta_list <- unlist(concat_delta_list)
sigma_list <- unlist(concat_sigma_list)

power_new_df <- data.frame(effect_size = effect_size,
                           naive_power = as.numeric(naive_p_value),
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value),
                           rand_index = unlist(concat_rand_list),
                           sigma = sigma_list,
                           delta = delta_list)



power_by_cut <- power_new_df %>%
  filter(effect_size>1e-5) %>%
  mutate(break_up_mu=cut_interval(effect_size, n=7,dig.lab=2)) %>%
  pivot_longer(-c(effect_size,delta,break_up_mu,sigma,rand_index), names_to = "p_type", 
               values_to = "p_value") %>%
  filter(p_type!="naive_power") %>%
  group_by(break_up_mu,p_type,sigma) %>%
  summarise(mean_power = mean(p_value<=alpha_thresh,na.rm=TRUE),
            sd_power = sqrt(mean_power*(1-mean_power)/n()))

manual_color_pal <- scales::hue_pal()(3)[2:3]

plot_power_one_d_K_4 <- power_by_cut %>%
  ggplot(aes(x=break_up_mu,y=mean_power,
             color=p_type, group = p_type))+ 
  geom_point()+
  geom_pointrange(aes(ymin = mean_power+1.96*sd_power, ymax = mean_power-1.96*sd_power))+
  geom_line()+
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('|$\\nu^{T}\\beta$|'))+
  theme_bw()+
  scale_linetype_discrete(name=TeX('$\\sigma$'))+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=12),
        legend.key.size = unit(2.5, "lines"),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  scale_color_manual( "Tests", 
                      values = c(manual_color_pal[1],manual_color_pal[2]),
                      labels = unname(TeX(c('$p_{Hyun}$',"$p_{C_{1},C_{2}}$"))))


png(paste0(plot_output_dir,'Figure_12_b.png'),
    width = 7,height = 7, res=300,units='in')
plot_power_one_d_K_4
dev.off()


#### plot Type I error control
load("./K_init_4_Type_I_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0_sim_times_1000_random_seed_2021.RData")

df_p0_plot <- data.frame(naive = (unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                         hyun = (unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                         union = (unlist(lapply(p_val_result,function(x)x[["Union"]]))))


df_p0_plot_type_1 <- df_p0_plot %>% 
  pivot_longer(everything(), names_to = "p_type", values_to = "p_val") %>%
  mutate(p_type = factor(as.factor(p_type), levels=c("naive","hyun","union"))) %>%
  group_by(p_type) %>%
  arrange(p_val) %>%
  mutate(theoretical = ecdf(p_val)(p_val)) %>%
  ungroup()

p_type_1_K_4 <- df_p0_plot_type_1 %>%
  ggplot(aes(y = p_val, x = theoretical,color = p_type)) +
  geom_point(size=0.6) +
  ylab('P-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  #ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_manual('Tests', 
                     values = scales::hue_pal()(3), 
                     labels = unname(TeX(c("$p_{Naive}$",'$p_{Hyun}','$p_{\\hat{C}_{1},\\hat{C}_{2}}'))))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))


png(paste0(plot_output_dir,'Figure_12_a.png'),
    width = 7,height = 7, res=300,units='in')
p_type_1_K_4
dev.off()

