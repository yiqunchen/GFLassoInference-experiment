library(tidyverse)
library(latex2exp)
library(scales)
library(intervals)

### plot type 1 error
input_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
setwd(input_dir)
plot_output_dir <- '~/Dropbox/dissertation_papers/GFLasso_manuscript/figure/revision_plot/'

alpha_thresh <- 0.05

load("./early_stopping_K_init_2_Type_I_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0_sim_times_1000_random_seed_2021.RData")

df_p0_plot_early_stop <- data.frame(naive = (unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                         hyun = (unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                         union = (unlist(lapply(p_val_result,function(x)x[["Union"]]))))


load("./Type_I_1D_GFL_middle_stop_criteria_K_grid_200_level_2_0_sim_times_1000_random_seed_2021.RData")

df_p0_plot_full_stop <- data.frame(naive = (unlist(lapply(p_val_result,function(x)x[["Naive"]]))),
                                    hyun = (unlist(lapply(p_val_result,function(x)x[["Hyun"]]))),
                                    union = (unlist(lapply(p_val_result,function(x)x[["Union"]]))))


early_stop_pval <- df_p0_plot_early_stop$union
full_stop_pval <- df_p0_plot_full_stop$union
# the same under type 1??
all.equal(early_stop_pval,full_stop_pval)
# "Mean relative difference: 2.373242e-06"

p_type_1_early_stop <- data.frame(early_stop = early_stop_pval,
           full_stop = full_stop_pval) %>%
  ggplot(aes(x=-log10(early_stop), y=-log10(full_stop))) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  geom_point(size=0.6) +
  ylab(TeX('$-log_{10}(p_{\\hat{C}_{1},\\hat{C}_{2}})$'))+
  xlab(TeX(paste0(('$-log_{10}(p_{\\hat{C}_{1},\\hat{C}_{2}})(\\delta)$'), ' with ',
           '$\\delta = max (0,10 \\sigma || \\nu ||-| \\nu^{T}y |$)') ))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))

png(paste0(plot_output_dir,'Figure8a.png'),
    width = 7,height = 7, res=300,units='in')
p_type_1_early_stop
dev.off()

sim_files_sigma_early <- list.files(path = input_dir,
                                  pattern = glob2rx("early_stopping_K_init_2_New_Power*"),
                                  full.names = TRUE)

concat_p_value_list_early_stop <- list()

for (current_file in sim_files_sigma_early){
  load(current_file)
  concat_p_value_list_early_stop <- c(concat_p_value_list_early_stop,p_val_result)
  rm(p_val_result)
}


# sigma sim files 
sim_files_sigma_full <- list.files(path = input_dir,
                                    pattern = glob2rx("*K_init_2*New_Power_1D_*_sigma_1*"),
                                    full.names = TRUE)

concat_p_value_list_full <- list()

for (current_file in sim_files_sigma_full){
  load(file=current_file)
  concat_p_value_list_full <- c(concat_p_value_list_full,p_val_result)
  rm(p_val_result)
}

power_p_value_early <- unlist(lapply(concat_p_value_list_early_stop, function(x)x[["Union"]]))
power_p_value_full <- unlist(lapply(concat_p_value_list_full, function(x)x[["Union"]]))


p_power_early_stop <- data.frame(early_stop = power_p_value_early,
                                  full_stop = power_p_value_full) %>%
  ggplot(aes(x=-log10(early_stop), y=-log10(full_stop))) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  geom_point(size=0.6) +
  ylab(TeX('$-log_{10}(p_{\\hat{C}_{1},\\hat{C}_{2}})$'))+
  xlab(TeX(paste0(('$-log_{10}(p_{\\hat{C}_{1},\\hat{C}_{2}}(\\delta))$'), ' with ',
                  '$\\delta = max (0,10 \\sigma || \\nu ||-| \\nu^{T}y |$)') ))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))


png(paste0(plot_output_dir,'Figure8b.png'),
    width = 7,height = 7, res=300,units='in')
p_power_early_stop
dev.off()

