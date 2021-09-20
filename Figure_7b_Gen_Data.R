library(GFLassoInference)
library(igraph)
library(mclust)

output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
random_seed <- 2021
set.seed(random_seed)
# simulation set-up
delta_seq <- seq(0.5,5,by=0.5)
sigma <- 1
true_cp <- c(100,140)
true_cc <- as.factor(c(rep(0,100),rep(1,40),rep(2,60)))
cp_detection_thresh <- 2
alpha_threshold <- 0.05
sim_target_times <- 500
# sim loop
for (delta in delta_seq){
  # initialize signal
  A <- matrix(c(rep(0,100),rep(delta, 40), rep(0, 60)),nrow=1)
  beta <- c(t(A))
  # create D matrix for 1D
  n <- length(beta)
  Dmat <- genlasso::getD1d(n)
  
  p_val_result <- vector('list',length = sim_target_times)
  rand_list <- vector('list',length = sim_target_times)
  counter_valid <- 0
  counter <- 0
  #sim_target_times <- 1000
  K <- 2
  stop_criteria <- "K"
  
  num_estimation <- 0
  num_correct_estimation <- 0
  num_correct_rejection_hyun <- 0
  num_correct_rejection_union <- 0
  
  while(counter_valid < sim_target_times){
    
    cat("delta",delta,"counter_valid",counter_valid,"\n")
    A.noisy <- A + rnorm(n,mean=0,sd=sigma) # create the noise vector
    y_vec <- do.call(c,lapply(1:nrow(A.noisy),function(irow)A.noisy[irow,]))
    y <- beta + rnorm(n,mean=0,sd=sigma) 
    complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
    resulting_cc <- complete_sol$pathobjs$i
    counter <- counter + 1
    # convert back to cps
    estimated_cp <- match(unique(resulting_cc),resulting_cc)-1
    estimated_cp <- estimated_cp[estimated_cp>0]
    #  only care about the first cp?
    #  maybe 2nd one is fine too
    correct_estimated_cp <- (abs(estimated_cp-true_cp)<=cp_detection_thresh)
    
    cc_order <- unique(resulting_cc)
    num_estimation <- num_estimation + length(estimated_cp)
    for (i in seq_along(correct_estimated_cp)){
      # ith estimated changepoint
      correct_detect <- correct_estimated_cp[i]
      # test only if we get a correct detection
      if(correct_detect){
        CC_to_test <- list(which(resulting_cc==cc_order[i]), which(resulting_cc==cc_order[i+1]))
        num_correct_estimation <- num_correct_estimation+1
        counter_valid <- counter_valid+1
        contrast <- rep(0, times = length(y_vec))
        contrast[CC_to_test[[1]]] = 1/(length(CC_to_test[[1]]))
        contrast[CC_to_test[[2]]] = -1/(length(CC_to_test[[2]]))
        mu_diff <- (sum(contrast*beta))
        # compute p-vals
        p_val_segment_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                            D=Dmat,
                                                            c1=cc_order[i], 
                                                            c2=cc_order[i+1],
                                                            method="K",
                                                            sigma=sigma,
                                                            K=K)
        p_union <- p_val_segment_cc$Union
        p_Hyun <- p_val_segment_cc$Hyun
        # add vTbeta, for our defn of power
        p_val_segment_cc$mu_diff <- mu_diff 
        
        if(p_union<=alpha_threshold){
          num_correct_rejection_union <- num_correct_rejection_union+1
        }
        if(p_Hyun<=alpha_threshold){
          num_correct_rejection_hyun <- num_correct_rejection_hyun+1
        }
        p_val_result[[counter_valid]] <- p_val_segment_cc
      }
    }
  }
  
  save(p_val_result,num_estimation,num_correct_estimation,
       num_correct_rejection_hyun,num_correct_rejection_union,
       file =paste0(output_dir,'cond_power_1D_GFL_middle_stop_criteria_',stop_criteria,'_grid_',n,'_level_2_',delta,'_sim_times_',
                    counter_valid,'_random_seed_',random_seed,'.RData'))
  
}




