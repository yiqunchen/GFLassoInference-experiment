library(GFLassoInference)
library(igraph)
library(mclust)
library(genlasso)

output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"

# simulation times 
sim_target_times <- 1000
p_val_result <- vector('list',length = sim_target_times)
rand_list <- vector('list',length = sim_target_times)

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

estimate_sigma <- function(y, cc_y){
  s <- split(y, cc_y)
  s_mean <- sapply(s,mean)
  rss <- sum((unlist(sapply(seq_along(s), function(i) s[[i]]-s_mean[[i]])))^2)
  s_sigma <- sqrt(rss/(length(y)-length(s)))
  return(s_sigma)
}

while(counter_valid < sim_target_times){
  
  # observation from the data generating mechanism
  y <- beta + rnorm(n,mean=0,sd=sigma) 
  complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
  resulting_cc <- complete_sol$pathobjs$i
  counter <- counter + 1
  # convert back to cps
  estimated_cp <- match(unique(resulting_cc),resulting_cc)-1
  estimated_cp <- estimated_cp[estimated_cp>0]
  cc_order <- unique(resulting_cc)
  
  for (i in seq_along(estimated_cp)){
    # ith estimated changepoint
    # test only if we get a correct detection
    CC_to_test <- list(which(resulting_cc==cc_order[i]), which(resulting_cc==cc_order[i+1]))
    counter_valid <- counter_valid+1
    contrast <- rep(0, times = length(y))
    contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
    contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
    mu_diff <- (sum(contrast*beta))
    sigma_hat <- estimate_sigma(y,resulting_cc)
    p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                         D=Dmat,
                                                         c1=cc_order[i], 
                                                         c2=cc_order[i+1],
                                                         method="K",
                                                         sigma=sigma_hat,
                                                         K=K)
    
    p_union <- p_val_segment_cc$Union
    p_Hyun <- p_val_segment_cc$Hyun
    # add vTbeta, for our defn of power
    p_val_segment_cc$mu_diff <- mu_diff 
    # store the result
    p_val_result[[counter_valid]] <- p_val_segment_cc
    
  }
}

save(p_val_result,
     file =paste0(output_dir,'estimated_sigma_residual_Type_I_1D_GFL_middle_stop_criteria_',
                  stop_criteria,'_grid_',n,'_level_2_',delta,'_sim_times_',
                  counter_valid,'_random_seed_',random_seed,'.RData'))




