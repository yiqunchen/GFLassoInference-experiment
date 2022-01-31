library(tidyverse)
library(scales)
library(latex2exp)
library(ggpubr)

alpha_thresh <- 0.05
input_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/more_powerful_GFL_experiments/plot_output/'

# collect data for sample variance
sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("estimated_sigma_sample_New_Power*.*"),full.names = T)

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
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value)) %>%
  pivot_longer(-effect_size, names_to = "p_type", values_to = "cond_power") %>%
  mutate(sigma_type = 'sample_variance')

write_csv(power_new_df, file = paste0(plot_output_dir, "one_d_power_sample_variance.csv"))

# collect data for MAD variance
sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("estimated_sigma_MAD_New_Power*.*"),full.names = T)

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
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value)) %>%
  pivot_longer(-effect_size, names_to = "p_type", values_to = "cond_power") %>%
  mutate(sigma_type = 'MAD')

write_csv(power_new_df, file = paste0(plot_output_dir, "one_d_power_MAD.csv"))

alpha_thresh <- 0.05

sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("estimated_sigma_residual_New_Power*.*"),full.names = T)

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
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value)) %>%
  pivot_longer(-effect_size, names_to = "p_type", values_to = "cond_power") %>%
  mutate(sigma_type = 'residual')

write_csv(power_new_df, file = paste0(plot_output_dir, "one_d_power_residual.csv"))

# true sigma
sim_files <- list.files(path = input_dir,
                        pattern = glob2rx("New_Power*_1500_*_2000.RData"),full.names = T)

sim_files <- sim_files[!grepl("sigma",sim_files)]

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
                           hyun_power = as.numeric(hyun_p_value), 
                           union_power = as.numeric(union_p_value)) %>%
  pivot_longer(-effect_size, names_to = "p_type", values_to = "cond_power") %>%
  mutate(sigma_type = 'True')

write_csv(power_new_df, file = paste0(plot_output_dir, "one_d_power_true_sigma.csv"))




alpha_thresh <- 0.05
files <- dir(plot_output_dir, pattern = "one_d_power_*.csv") # get file names

power_data <- files %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(file.path(data_path, .))) %>% 
  reduce(rbind)


power_by_cut <- power_data %>%
  filter(effect_size>0) %>%
  mutate(break_up_mu = cut_interval(effect_size, n=7,dig.lab=0)) %>%
  group_by(break_up_mu,p_type, sigma_type) %>%
  summarise(mean_power = mean(cond_power<=alpha_thresh),
            sd_power = sqrt(mean(cond_power)*(1-mean(cond_power))/n()))

manual_color_pal <- scales::hue_pal()(3)[2:3]

plot_power_all_estimator <- power_by_cut %>% 
  ggplot(aes(x=break_up_mu,y=mean_power,color=p_type,linetype = sigma_type,
             group = interaction(p_type,sigma_type)))+ 
  facet_wrap(.~sigma_type, nrow=2, scales="fixed")+
  geom_point()+
  geom_pointrange(aes(ymin = mean_power+1.96*sd_power, ymax = mean_power-1.96*sd_power))+
  geom_line()+
  ylab(TeX('Power at $\\alpha = 0.05$'))+
  xlab(TeX('$| \\nu^T\\beta |$'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=17),
        legend.position="bottom",
        legend.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.key.size = unit(2.5, "lines"),
        axis.text.x = element_text(size=13),
        legend.text = element_text(size=20,hjust = 0),
        axis.title=element_text(size=20),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  scale_color_manual("Tests", values = c(manual_color_pal[1],manual_color_pal[2]),
                      breaks = c("hyun_power","union_power"),
                      labels = unname(TeX(c('$p_{Hyun}$',
                                            "$p_{\\hat{C}_{1},\\hat{C}_{2}}$"))))+
  scale_linetype_manual("Variance", 
                        breaks = c("TRUE","MAD","residual","sample_variance"),
                        labels = unname(TeX(c("$\\sigma$", 
                                              "$\\hat{\\sigma}_{MAD}$",
                                              "$\\hat{\\sigma}_{Residual}$",
                                              "$\\hat{\\sigma}_{Sample}$"))), 
                        values = c("solid","twodash","longdash","dotted"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE, order=1),
    linetype=guide_legend(nrow=2,byrow=TRUE,order=0))


#### think about how to plot them nicely
png(paste0(data_path,'Figure_14.png'),
    width = 15,height = 13, res=300,units='in')
plot_power_all_estimator
dev.off()







