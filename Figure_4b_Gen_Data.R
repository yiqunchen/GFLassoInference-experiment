library(GFLassoInference)
library(igraph)
library(mclust)

random_seed <- 2021
set.seed(random_seed)

output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
delta <- 0
stop_criteria <- 'K' 
lev1 <- 0 # mean for group 1
lev2 <- 0 # delta
sigma <- 1 # level of noise
nn <- 8 # grid size; vary 5 or 10
Dmat <- genlasso::getD2d(nn, nn) # create D matrix 
# create three strips
A <- matrix(lev1, ncol=nn, nrow = nn)
A[1:(nn/3+1),1:(nn/3+1)] <- 1*delta 
A[(nn-2):(nn),(nn-2):(nn)] <- -1*delta
# create signal beta
beta <- c(t(A)) # create true mean vector
true_cc <- as.factor(beta)

K <- 15
# L <- 3
# number of simulations
sim_target_times <- 1000
p_val_result <- vector('list',length = sim_target_times)
rand_list <- vector('list',length = sim_target_times)
counter_valid <- 0
counter <- 0

while(counter_valid<sim_target_times){
  # sample a noisy observation
  y <- beta + rnorm(nn^2,mean=0,sd=sigma) 
  # solve fused lasso
  complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
  resulting_cc <- complete_sol$pathobjs$i
  # k_cc: number of estimated connected components 
  K_cc <- complete_sol$pathobjs$q
  # only perform inference if we estimated more than one CC
  if (K_cc>1){
    cat("counter_valid", counter_valid,"\n")
    counter_valid <- counter_valid+1
    current_adj_rand <- mclust::adjustedRandIndex(resulting_cc,true_cc)
    # randomly choose a pair to test, and construct contrast
    contrast <- rep(0, times = length(y))
    pair_to_test <- sample(c(1:K_cc), replace = F,size = 2)
    CC_to_test <- list(which(resulting_cc==pair_to_test[1]), 
                       which(resulting_cc==pair_to_test[2]))
    contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
    contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
    mu_diff <- (sum(contrast*beta))
    
    p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                         D=Dmat,
                                                         c1=pair_to_test[1], 
                                                         c2=pair_to_test[2],
                                                         method="K",
                                                         sigma=sigma,
                                                         K=K)

    p_val_segment_cc$mu_diff <- mu_diff 
    cat(p_val_segment_cc$Hyun,p_val_segment_cc$Union,"\n")
    p_val_result[[counter_valid]] <- p_val_segment_cc
    rand_list[[counter_valid]] <- current_adj_rand
    
  }
}


save(p_val_result,rand_list,
     file =paste0( output_dir,'Type_1_2D_GFL_stop_criteria_',
                   stop_criteria,'_grid_',nn,'_level_2_',lev2,
                   '_sim_times_',counter_valid,'_random_seed_',random_seed,'.RData'))




