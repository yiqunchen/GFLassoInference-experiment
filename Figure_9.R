library(tidyverse)
library(latex2exp)
library(ggpubr)
library(scales)
library(intervals)

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

png(paste0(plot_output_dir,'Figure_9_a.png'),
    width = 6,height = 6, res=200,units='in')

print(type_1_2d)

dev.off()

manual_color_pal <- scales::hue_pal()(3)[2:3]
### plot power 2D
input_dir  <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'
alpha_thresh <- 0.05
# load power data
setwd(input_dir)
sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("*2D*stop_criteria_CC*500*.*"),full.names = T)

concat_p_value_list <- list()
concat_rand_list <- c()
concat_delta_list <- c()
concat_sigma_list <- c()

for (current_file in sim_files){
  load(file=current_file)
  delta_strsplit <- strsplit(strsplit(current_file,"level_2_")[[1]][2],
                             "_sim_times")[[1]][1]
  sigma_strsplit <- strsplit(strsplit(current_file,"sigma_")[[1]][2],
                             "_random_seed")[[1]][1]
  concat_p_value_list <- c(concat_p_value_list,p_val_result)
  concat_rand_list <- c(concat_rand_list,rand_list)
  concat_delta_list <- c(concat_delta_list,
                         rep(as.numeric(delta_strsplit), times = length(rand_list)))
  concat_sigma_list <- c(concat_sigma_list,
                         (rep(as.numeric(sigma_strsplit), times = length(rand_list))))
  rm(p_val_result)
}

concat_rand_list <- unlist(concat_rand_list)
hyun_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Hyun"]]))
union_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Union"]]))
effect_size <- (unlist(lapply(concat_p_value_list, function(x)x[["mu_diff"]])))
nu_norm <- unlist(lapply(concat_p_value_list, function(x)x[["sd"]]))
delta_list <- unlist(concat_delta_list)
sigma_list <- unlist(concat_sigma_list)

power_new_df <- data.frame(effect_size=effect_size,
                           rand_index = concat_rand_list,
                           naive_power = as.numeric(naive_p_value),
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value),
                           sigma = sigma_list,
                           delta = delta_list) %>%
  pivot_longer(-c(rand_index,effect_size,delta,sigma), names_to = "p_type", 
               values_to = "p_value") 

power_by_cut <- power_new_df %>%
  filter(effect_size>1e-5) %>%
  mutate(break_up_mu = cut_interval(effect_size, n=7,dig.lab=2)) %>%
  group_by(break_up_mu,p_type,sigma) %>%
  summarise(mean_power = mean(p_value<=alpha_thresh),
            sd_power = sqrt(mean_power*(1-mean_power))/n())


plot_power <- power_by_cut %>%
  filter(p_type!="naive_power") %>%
  ggplot(aes(x=break_up_mu,y=mean_power, 
             color=p_type,
             group=interaction(p_type,sigma),
             linetype=as.factor(sigma)))+ 
  geom_point(size=2)+
  geom_pointrange(aes(ymin = mean_power+1.96*sd_power, 
    ymax = mean_power-1.96*sd_power))+
  geom_line()+
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('$| \\nu^T\\beta |$'))+
  theme_bw()+
  scale_linetype_discrete(name=TeX('$\\sigma$'))+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=15,hjust = 0),
        legend.key.size = unit(2.5, "lines"),
        axis.title=element_text(size=15))+
  scale_color_manual( "Tests", 
                      values = c(manual_color_pal[1],manual_color_pal[2]),
                      labels = unname(TeX(c('$p_{Hyun}$',"$p_{C_{1},C_{2}}$"))))

png(paste0(plot_output_dir,'Figure_9_b.png'),
    width = 7,height = 7, res=200,units='in')
plot_power
dev.off()

# alternative definition of power: consider detection prob and cond. power separately
df_rand_index <- power_new_df %>%
  group_by(delta,sigma) %>%
  summarise(avg_ARI = mean(rand_index), 
            sd_ARI = sd(rand_index)/sqrt(n()),
            mean_perfect_recovery = mean(rand_index==1),
            sd_perfect_recovery = 
              sqrt(mean_perfect_recovery*(1-mean_perfect_recovery)/n()))

p_detect_prob <- df_rand_index %>%
  ggplot(aes(x=delta,y=mean_perfect_recovery,
             linetype=as.factor(sigma),
             group=as.factor(sigma)))+ 
  geom_point()+
  geom_pointrange(aes(ymin = mean_perfect_recovery+1.96*sd_perfect_recovery, 
                      ymax = mean_perfect_recovery-1.96*sd_perfect_recovery))+
  geom_line()+
  ylab(TeX('Detection Probability'))+
  xlab(TeX('Difference in means: $\\delta$'))+
  theme_bw(base_size=18) + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        legend.position="bottom",
        legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.key.size = unit(2.5, "lines"),
        legend.text = element_text(size=18,hjust = 0),
        axis.title=element_text(size=18))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_linetype_discrete(name=TeX('$\\sigma$'))+
  scale_y_continuous(limits = c(0,1))
  
png(paste0(plot_output_dir,'Figure_9_c.png'),
    width = 7,height = 7, res=300,units='in')
p_detect_prob
dev.off()

cond_power_by_cut <- power_new_df %>%
  filter(effect_size>1e-5) %>%
  filter(rand_index==1) %>%
  group_by(delta,p_type,sigma) %>%
  summarise(mean_power = mean(p_value<=alpha_thresh),
            sd_power = sqrt(mean_power*(1-mean_power))/n()) %>%
  ungroup()%>%
  add_row(delta=0.5,p_type="hyun_power",sigma=1,mean_power=0,sd_power=0) %>%
  add_row(delta=0.5,p_type="union_power",sigma=1,mean_power=0,sd_power=0) %>%
  add_row(delta=1,p_type="union_power",sigma=2,mean_power=0,sd_power=0)%>%
  add_row(delta=1,p_type="hyun_power",sigma=2,mean_power=0,sd_power=0)

plot_cond_power <- cond_power_by_cut %>%
  filter(p_type!="naive_power") %>%
  ggplot(aes(x=delta,y=mean_power,
             color=p_type,
             linetype=as.factor(sigma)))+ 
  geom_point()+
  geom_pointrange(aes(ymin = mean_power+1.96*sd_power, ymax = mean_power-1.96*sd_power))+
  geom_line()+
  ylab(TeX('Conditional Power at $\\alpha = 0.05$'))+
  xlab(TeX('Difference in means: $\\delta$'))+
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

png(paste0(plot_output_dir,'Figure_9_d.png'),
    width = 7,height = 7, res=300,units='in')
plot_cond_power
dev.off()





