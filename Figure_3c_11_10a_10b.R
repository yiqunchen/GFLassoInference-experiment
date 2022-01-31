library(tidyverse)
library(latex2exp)
library(scales)

manual_color_pal <- scales::hue_pal()(3)[2:3]

### plot power
input_dir <-  "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/"
alpha_thresh <- 0.05
num_sims <- 1500

sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("New_Power_1D_GFL*1500*.*"),full.names = T)

concat_p_value_list <- list()
concat_delta_list <- c()
concat_sigma_list <- c()

for (current_file in sim_files){
  # extract delta and sigma sequences
  sigma_strsplit <- as.numeric(strsplit(strsplit(current_file,"sigma_")[[1]][2],
                      "_sim_times")[[1]][1])
  delta_strsplit <- as.numeric(strsplit(strsplit(current_file,"level_2_")[[1]][2],
                                          "_sigma")[[1]][1])
  load(file=current_file)
  concat_delta_list <- c(concat_delta_list,
                         rep((delta_strsplit), times = length(p_val_result)))
  concat_sigma_list <- c(concat_sigma_list,
                         (rep((sigma_strsplit), times = length(p_val_result))))
  concat_p_value_list <- c(concat_p_value_list,p_val_result)
  rm(p_val_result)
}

hyun_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Hyun"]]))
union_p_value <- unlist(lapply(concat_p_value_list, function(x)x[["Union"]]))
effect_size <- (unlist(lapply(concat_p_value_list, function(x)x[["mu_diff"]])))
nu_norm <- unlist(lapply(concat_p_value_list, function(x)x[["sd"]]))
rand_index <- unlist(lapply(concat_p_value_list, 
  function(x)x[["current_adj_rand"]]))
delta_list <- unlist(concat_delta_list)
sigma_list <- unlist(concat_sigma_list)

power_new_df <- data.frame(effect_size=effect_size,
                           naive_power = as.numeric(naive_p_value),
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value),
                           rand_index = unlist(rand_index),
                           sigma = sigma_list,
                           delta = delta_list)


# compute power in each bucket of effect size 
power_by_cut <- power_new_df %>% 
  filter(effect_size>1e-5) %>%
  pivot_longer(-c(effect_size,delta,sigma), names_to = "p_type", 
               values_to = "p_value")  %>%
  mutate(break_up_mu=cut_interval(effect_size, n=7,dig.lab=2)) %>%
  group_by(break_up_mu,p_type,sigma) %>%
  summarise(mean_power = mean(p_value<=alpha_thresh,na.rm=TRUE),
            sd_power = sqrt(mean_power*(1-mean_power)/n()))

plot_power <- power_by_cut %>%
  filter(p_type!="naive_power") %>%
  ggplot(aes(x=break_up_mu,y=mean_power,
             color=p_type,
             group=interaction(p_type,sigma),
             linetype=as.factor(sigma)))+
  geom_point()+
  geom_pointrange(aes(ymin = mean_power+1.96*sd_power,
                      ymax = mean_power-1.96*sd_power))+
  geom_line()+
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('$| \\nu^T\\beta |$'))+
  theme_bw()+
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
    breaks = c("hyun_power","union_power"),
    labels = unname(TeX(c('$p_{Hyun}$',"$p_{\\hat{C}_{1},\\hat{C}_{2}}$"))))+
  scale_linetype_discrete(name=TeX('$\\sigma$'))

png(paste0(plot_output_dir,'Figure_3_c.png'),
    width = 7,height = 7, res=200,units='in')
plot_power
dev.off()


plot_1d_smooth_function <- power_new_df %>% 
  pivot_longer(-c(effect_size,delta,sigma), names_to = "p_type", 
               values_to = "p_value")  %>%
  filter(sigma>=1,p_type!="naive_power")%>%
  drop_na() %>%
  ggplot(aes(x=effect_size,  
             y=as.numeric(p_value<=alpha_thresh),
             linetype=as.factor(sigma),
             colour = p_type)) +
  geom_smooth()+
  #geom_line() +
  theme_bw() +
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('Effect size: $| \\nu^T\\beta |$'))+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))+
  scale_color_manual( "Tests", values = manual_color_pal,
    breaks = c("hyun_power","union_power"),
    labels = unname(TeX(c('$p_{Hyun}$','$p_{\\hat{C}_{1},\\hat{C}_{2}}'))) )+
  scale_linetype_discrete(name=TeX('$\\sigma$'))

png(paste0(plot_output_dir,'Figure_11.png'),
    width = 7,height = 7, res=200,units='in')
plot_1d_smooth_function
dev.off()

p_detect_ggplot <- power_new_df %>% group_by(sigma,delta) %>%
  summarise(avg_detect_p = mean(rand_index==1),
            sd_detect_p = sqrt(avg_detect_p*(1-avg_detect_p)/n())) %>%
ggplot(aes(x=delta,y=avg_detect_p,
           linetype=as.factor(sigma),
           group=as.factor(sigma)))+ 
  geom_point()+
  geom_pointrange(aes(ymin = avg_detect_p+1.96*sd_detect_p, 
                      ymax = avg_detect_p-1.96*sd_detect_p))+
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


png(paste0(plot_output_dir,'Figure_10a.png'),
    width = 7,height = 7, res=300,units='in')
p_detect_ggplot
dev.off()

cond_power_by_cut <- power_new_df %>%
  filter(effect_size>1e-5) %>%
  filter(rand_index==1) %>%
  pivot_longer(-c(effect_size,delta,sigma,rand_index), names_to = "p_type", 
               values_to = "p_value") %>%
  filter(p_type!="naive_power") %>%
  group_by(delta,p_type,sigma) %>%
  summarise(mean_power = mean(p_value<=alpha_thresh,na.rm=TRUE),
            sd_power = sqrt(mean_power*(1-mean_power))/n()) %>%
  ungroup()


plot_cond_power_one_d <- cond_power_by_cut %>%
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

png(paste0(plot_output_dir,'Figure_10b.png'),
    width = 7,height = 7, res=300,units='in')
plot_cond_power_one_d
dev.off()




