library(PGInference)
library(igraph)
library(mclust)
library(Matrix)
library(MASS)

# set output dir
output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
# set random seed
random_seed <- 2021
set.seed(random_seed)
# set parameters for the power experiment
stop_criteria <- "K"
lev1 <- 0 # mean for group 1
lev2 <- delta # mean for groyp 2
sigma <- 1 # level of noise
nn <- 8 # grid size; vary 5 or 10
Dmat <- genlasso::getD2d(nn, nn) # create D matrix 
# create three strips
A <- matrix(lev1,ncol=nn, nrow = nn)
# create three strips
A[1:(nn/3+1),1:(nn/3+1)] <- 1*delta 
A[(nn-2):(nn),(nn-2):(nn)] <- -1*delta
# create true mean vector
beta <- c(t(A))
true_cc <- as.factor(beta)
# more parameters for running FLasso
K <- 15
L <- 3
sim_target_times <- 500
delta_seq <- seq(0.5,5,by=0.5)

for (delta in delta_seq){
  p_val_result <- vector('list',length = sim_target_times)
  rand_list <- vector('list',length = sim_target_times)
  counter_valid <- 0
  counter <- 0
  while(counter_valid<sim_target_times){
    # create the noise vector
    A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma) 
    # sample noisy observation
    y <- c(t(A.noisy))
    # fused lasso estimator
    complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
    resulting_cc <- complete_sol$pathobjs$i
    # number of estimated CC
    K_cc <- complete_sol$pathobjs$q
    if (K_cc>1){
      #cat("counter_valid", counter_valid,"\n")
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
      p_val_result[[counter_valid]] <- p_val_segment_cc
      rand_list[[counter_valid]] <- current_adj_rand
    }
  }

  save(p_val_result,rand_list,
       file =paste0(output_dir,'New_Power_2D_GFL_stop_criteria_',
                    stop_criteria,'_grid_',nn,'_level_2_',delta,
                    '_sim_times_',counter_valid,'_random_seed_',random_seed,'.RData'))
  

}


