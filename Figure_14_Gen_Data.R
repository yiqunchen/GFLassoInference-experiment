library(GFLassoInference)
library(igraph)
library(mclust)
library(genlasso)

output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"

sim_target_times <- 1500
# set random seed
random_seed <- 2000
set.seed(random_seed)
# initial parameters
delta_seq <- seq(0.5,5,by=0.5)
sigma_seq <- c(1)

estimate_sigma <- function(y, cc_y){
  s <- split(y, cc_y)
  s_mean <- sapply(s,mean)
  rss <- sum((unlist(sapply(seq_along(s), function(i) s[[i]]-s_mean[[i]])))^2)
  s_sigma <- sqrt(rss/(length(y)-length(s)))
  return(s_sigma)
}

# using the residual based variance estimator 
for (sigma in sigma_seq){
  for (delta in delta_seq){
    # simulation times
    p_val_result <- vector('list',length = sim_target_times)
    rand_list <- vector('list',length = sim_target_times)
    
    # initialize signal
    A <- matrix(c(rep(0,100),rep(delta, 40), rep(0, 60)),nrow=1)
    beta <- c(t(A))
    # create D matrix for 1D
    n <- length(beta)
    Dmat <- genlasso::getD1d(n)
    # true piecewise constant regions
    true_cc <- as.factor(c(rep(0,100),rep(1,40),rep(2,60)))
    
    counter_valid <- 0
    counter <- 0
    K <- 2
    # L <- 3
    stop_criteria <- "K"
    
    while(counter_valid < sim_target_times){
      # sample noisy observations
      y <- beta + rnorm(n,mean=0,sd=sigma) 
      complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
      resulting_cc <- complete_sol$pathobjs$i
      N_CC <- length(unique(resulting_cc))
      # compute the rand index 
      current_adj_rand <- mclust::adjustedRandIndex(resulting_cc,true_cc)
      cc_order <- unique(resulting_cc)
      # sample the pair of connected components
      index_i <- sample(c(1:(N_CC-1)), replace = F,size = 1)
      # connected components to test
      CC_to_test <- list(which(resulting_cc==cc_order[index_i]), 
                         which(resulting_cc==cc_order[index_i+1]))
      counter_valid <- counter_valid+1
      contrast <- rep(0, times = length(y))
      contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
      contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
      mu_diff <- (sum(contrast*beta))
      # compute p-values
      sigma_hat <- estimate_sigma(y,resulting_cc)
      p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                           D=Dmat,
                                                           c1=cc_order[index_i], 
                                                           c2=cc_order[index_i+1],
                                                           method="K",
                                                           sigma=sigma_hat,
                                                           K=K)
      cat("counter", counter_valid, p_val_segment_cc$Union,"\n")
      p_val_segment_cc$mu_diff <- mu_diff 
      p_val_segment_cc$ <- current_adj_rand
      rand_list[[counter_valid]] <- current_adj_rand
      p_val_result[[counter_valid]] <- p_val_segment_cc
      
    }
  }
    
  save(p_val_result,rand_list,
       file =paste0(output_dir,'estimated_sigma_residual_New_Power_1D_GFL_middle_stop_criteria_',stop_criteria,'_grid_',n,
                    '_level_2_',delta,'_sigma_',sigma,'_sim_times_',
                    counter_valid,'_random_seed_',random_seed,'.RData'))
  
}

set.seed(random_seed)

# using the sample variance estimator 
for (sigma in sigma_seq){
  for (delta in delta_seq){
    # simulation times
    p_val_result <- vector('list',length = sim_target_times)
    rand_list <- vector('list',length = sim_target_times)
    
    # initialize signal
    A <- matrix(c(rep(0,100),rep(delta, 40), rep(0, 60)),nrow=1)
    beta <- c(t(A))
    # create D matrix for 1D
    n <- length(beta)
    Dmat <- genlasso::getD1d(n)
    # true piecewise constant regions
    true_cc <- as.factor(c(rep(0,100),rep(1,40),rep(2,60)))
    
    counter_valid <- 0
    counter <- 0
    K <- 2
    # L <- 3
    stop_criteria <- "K"
    
    while(counter_valid < sim_target_times){
      # sample noisy observations
      y <- beta + rnorm(n,mean=0,sd=sigma) 
      complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
      resulting_cc <- complete_sol$pathobjs$i
      N_CC <- length(unique(resulting_cc))
      # compute the rand index 
      current_adj_rand <- mclust::adjustedRandIndex(resulting_cc,true_cc)
      cc_order <- unique(resulting_cc)
      # sample the pair of connected components
      index_i <- sample(c(1:(N_CC-1)), replace = F,size = 1)
      # connected components to test
      CC_to_test <- list(which(resulting_cc==cc_order[index_i]), 
                         which(resulting_cc==cc_order[index_i+1]))
      counter_valid <- counter_valid+1
      contrast <- rep(0, times = length(y))
      contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
      contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
      mu_diff <- (sum(contrast*beta))
      # compute p-values
      sigma_hat <- sd(y)
      p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                           D=Dmat,
                                                           c1=cc_order[index_i], 
                                                           c2=cc_order[index_i+1],
                                                           method="K",
                                                           sigma=sigma_hat,
                                                           K=K)
      cat("counter", counter_valid, p_val_segment_cc$Union,"\n")
      p_val_segment_cc$mu_diff <- mu_diff 
      p_val_segment_cc$ <- current_adj_rand
      rand_list[[counter_valid]] <- current_adj_rand
      p_val_result[[counter_valid]] <- p_val_segment_cc
      
    }
  }
    
  save(p_val_result,rand_list,
       file =paste0(output_dir,'estimated_sigma_sample_New_Power_1D_GFL_middle_stop_criteria_',stop_criteria,'_grid_',n,
                    '_level_2_',delta,'_sigma_',sigma,'_sim_times_',
                    counter_valid,'_random_seed_',random_seed,'.RData'))
  
}

set.seed(random_seed)

# using the MAD estimator 
for (sigma in sigma_seq){
  for (delta in delta_seq){
    # simulation times
    p_val_result <- vector('list',length = sim_target_times)
    rand_list <- vector('list',length = sim_target_times)
    
    # initialize signal
    A <- matrix(c(rep(0,100),rep(delta, 40), rep(0, 60)),nrow=1)
    beta <- c(t(A))
    # create D matrix for 1D
    n <- length(beta)
    Dmat <- genlasso::getD1d(n)
    # true piecewise constant regions
    true_cc <- as.factor(c(rep(0,100),rep(1,40),rep(2,60)))
    
    counter_valid <- 0
    counter <- 0
    K <- 2
    # L <- 3
    stop_criteria <- "K"
    
    while(counter_valid < sim_target_times){
      # sample noisy observations
      y <- beta + rnorm(n,mean=0,sd=sigma) 
      complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
      resulting_cc <- complete_sol$pathobjs$i
      N_CC <- length(unique(resulting_cc))
      # compute the rand index 
      current_adj_rand <- mclust::adjustedRandIndex(resulting_cc,true_cc)
      cc_order <- unique(resulting_cc)
      # sample the pair of connected components
      index_i <- sample(c(1:(N_CC-1)), replace = F,size = 1)
      # connected components to test
      CC_to_test <- list(which(resulting_cc==cc_order[index_i]), 
                         which(resulting_cc==cc_order[index_i+1]))
      counter_valid <- counter_valid+1
      contrast <- rep(0, times = length(y))
      contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
      contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
      mu_diff <- (sum(contrast*beta))
      # compute p-values
      sigma_hat <- mad(diff(y,1))/sqrt(2)
      p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                           D=Dmat,
                                                           c1=cc_order[index_i], 
                                                           c2=cc_order[index_i+1],
                                                           method="K",
                                                           sigma=sigma_hat,
                                                           K=K)
      cat("counter", counter_valid, p_val_segment_cc$Union,"\n")
      p_val_segment_cc$mu_diff <- mu_diff 
      p_val_segment_cc$ <- current_adj_rand
      rand_list[[counter_valid]] <- current_adj_rand
      p_val_result[[counter_valid]] <- p_val_segment_cc
      
    }
  }
    
  save(p_val_result,rand_list,
       file =paste0(output_dir,'estimated_sigma_MAD_New_Power_1D_GFL_middle_stop_criteria_',stop_criteria,'_grid_',n,
                    '_level_2_',delta,'_sigma_',sigma,'_sim_times_',
                    counter_valid,'_random_seed_',random_seed,'.RData'))
  
}

