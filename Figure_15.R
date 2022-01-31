library(GFLassoInference)
library(igraph)
library(mclust)
library(genlasso)
library(ggplot2)

plot_output_dir <- '~/Desktop/dissertation/gflasso_project/more_powerful_GFL_experiments/plot_output/'

# simulation times 
sim_target_times <- 1000
p_val_result <- vector('list',length = sim_target_times)
rand_list <- vector('list',length = sim_target_times)
elasped_time <- rep(NA, times=sim_target_times)

# initialize the setup for the simulation
random_seed <- 2021
set.seed(random_seed)
delta <- 0
sigma <- 1
A <- matrix(c(rep(0,100),rep(delta, 40), rep(0, 60)),nrow=1)
n <- length(A)
# create signal vector beta
beta <- c(t(A))
# create D matrix for 1D
Dmat <- genlasso::getD1d(n)
counter_valid <- 0
counter <- 0
# Number of steps for the dual path algorithms
K <- 2
# L <- 3
stop_criteria <- "K"


while(counter_valid < sim_target_times){
  cat("counter_valid", counter_valid,"\n")
  # observation from the data generating mechanism
  y <- beta + rnorm(n,mean=0,sd=sigma) 
  complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
  resulting_cc <- complete_sol$pathobjs$i
  counter <- counter + 1
  # convert back to cps
  estimated_cp <- match(unique(resulting_cc),resulting_cc)-1
  estimated_cp <- estimated_cp[estimated_cp>0]
  cc_order <- unique(resulting_cc)
  
  all_pairs <- combn(length(cc_order),2)
  for (i in c(1:ncol(all_pairs)) ){
    # ith estimated changepoint
    # test only if we get a correct detection
    CC_to_test <- list(which(resulting_cc==all_pairs[1,i]),
                       which(resulting_cc==all_pairs[2,i]))
    counter_valid <- counter_valid+1
    contrast <- rep(0, times = length(y))
    contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
    contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
    mu_diff <- (sum(contrast*beta))
    curr_time <- system.time(p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                         D=Dmat,
                                                         c1=all_pairs[1,i], 
                                                         c2=all_pairs[2,i],
                                                         method="K",
                                                         sigma=sigma,
                                                         K=K))
    
    p_union <- p_val_segment_cc$Union
    p_Hyun <- p_val_segment_cc$Hyun
    # add vTbeta, for our defn of power
    p_val_segment_cc$mu_diff <- mu_diff 
    # store the result
    p_val_result[[counter_valid]] <- p_val_segment_cc
    elasped_time[counter_valid] <- curr_time[3]
    
  }
}

p_timing <- ggplot(data.frame(time = elasped_time[1:1000]), aes(x=time)) + 
  geom_histogram(color="black", fill="white",aes(y=stat(count)/sum(count)))+
  theme(plot.title = element_text(hjust = 0.5,size=17),
        legend.title = element_text(size=17),
        axis.text = element_text(size= 17),
        legend.text = element_text(size=17,hjust = 0),
        axis.title=element_text(size=17))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw(base_size=17) +
  xlab("Time (seconds)")+
  ylab("Frequency")+
  scale_y_continuous(labels = scales::percent_format()) 
  

png(paste0(plot_output_dir,'Figure_15b.png'),
    width = 6,height = 4, res=300,units='in')
p_timing
dev.off()

I_card <- unlist(lapply(p_val_result, function(x)x[["sum_intervals"]]))

all_pvals <- unlist(lapply(p_val_result, function(x)x[["pval"]]))

p_aggregate_type_1 <- data.frame(p_values = all_pvals) %>%
  filter(p_values!=1) %>%
  mutate(theoretical = ecdf(p_values)(p_values)) %>%
  ggplot() +
  geom_point(aes(y = p_values, x = (theoretical)),size=0.6) +
  ylab('P-value Quantiles')+
  xlab('Uniform(0,1) Quantiles')+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  geom_abline(intercept = 0, slope = 1)
  

I_1d_plot_type_1 <- ggplot(data.frame(I_card = I_card), aes(x=I_card)) + 
  geom_histogram(color="black", fill="white",aes(y=stat(count)/sum(count)))+
  theme(plot.title = element_text(hjust = 0.5,size=17),
        #legend.position="bottom",
        legend.title = element_text(size=17),
        axis.text = element_text(size= 17),
        legend.text = element_text(size=17,hjust = 0),
        axis.title=element_text(size=17))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw(base_size=17) +
  xlab("Number of Graph Fused Lasso instances")+
  scale_fill_discrete(name=TeX('$\\sigma$'), labels=c(0.5,1,2))+
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position='none')+
  ylab("Frequency")


png(paste0(plot_output_dir,'Figure_15a.png'),
    width = 6,height = 4, res=300,units='in')
I_1d_plot_type_1
dev.off()



