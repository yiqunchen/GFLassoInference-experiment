
### load  libraries
library(truncnorm)
library(latex2exp)
library(ggpubr)
library(GFLassoInference)
library(genlasso)
library(ggplot2)
library(dplyr)

### Set color scales
library(lattice)
library(RColorBrewer)


## utility functions for plotting the density in Figure 1d
log_sum_exp <- function(logx, logy) {
  if (logx > logy) {
    a = logx;
  } else {
    a = logy;
  }
  if (abs(a) == Inf) {
    a = 0;
  }
  out = exp(logx - a);
  out = out + exp(logy - a);
  return(log(out) + a)
}


log1mexp <- function(a) {
  if (a >= 0 && a <= log(2)) {
    return(log(-expm1(-a)));
  } else if (a > log(2)) {
    return(log1p(-exp(-a)));
  } else {
    stop(paste0("trying to log1mexp with a = ", a))
  }
}
log_subtract <- function(x, y) {
  if (x < y) {
    stop("log_subtract:: cannot take log of (-)ve number");
  }
  return(x + log1mexp(abs(y - x)));
}

g_double_truncation <- function(x, mean = 0, sd = 1, a, b, c, d) {
  normalization_1 <-  exp(log_subtract(pnorm(b, mean = mean, sd = sd, log.p = TRUE) ,
                                       pnorm(a, mean = mean, sd = sd, log.p = TRUE)))
  
  normalization_2 <-  exp(log_subtract(pnorm(d, mean = mean, sd = sd, log.p = TRUE) ,
                                       pnorm(c, mean = mean, sd = sd, log.p = TRUE)))
  
  new_density <- ifelse((x > a & x <= b) | (x > c & x <= d),
                        dnorm(x, mean = mean, sd = sd) /
                          (normalization_1+normalization_2),
                        0)
  
  return(new_density)
}

g_single_truncation <- function(x, mean = 0, sd = 1, a, b) {
  ifelse((x > a & x <= b),
         dnorm(x, mean = mean, sd = sd) /
           exp(log_subtract(pnorm(b, mean = mean, sd = sd, log.p = TRUE) ,
                            pnorm(a, mean = mean, sd = sd, log.p = TRUE))),
         0)
}


input_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/"
plot_output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/"

# test the dual path fused and svd functions
lev1 <- 0 # mean for group 1
delta <- 3 # mean for group 2; null simulation
sigma <- 1 # level of noise
nn <- 8 # grid size
Dmat <- genlasso::getD2d(nn, nn) # create D matrix
A <- matrix(lev1,ncol=nn, nrow = nn)
# create three strips
# fix the size of boundary - unfair comparison otw
A[1:round(nn/3),1:round(nn/3)] <- 1*delta
A[(nn-2):(nn),(nn-2):(nn)] <- -1*delta
# create true mean vector
beta <- c(t(A))


set.seed(2005)
A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma)
A_noisy_value_plot <- data.frame(expand.grid(c(1:8),c(1:8)), z = c(t(A.noisy)))
# create the noise vector
y <- c(t(A.noisy))

A_true_value_plot <- data.frame(expand.grid(c(1:8),c(1:8)), z = c(t(A)))


#### 
p_A_true_ggplot <- ggplot(A_true_value_plot, aes(x=Var1, y=Var2, fill = z)) +
  geom_tile() +
  xlab("") + ylab("") +
  ggtitle("True signal")+
  scale_fill_distiller(palette = "PiYG")+
  theme_bw()+
  theme(plot.title=element_text(size=15,hjust = 0.5),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.title=element_text(size=15),
        strip.text.y = element_text(size=15,hjust=0,vjust = 1,angle=180,face="bold"))+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  annotate("text",x = 2 ,y = 2, label = TeX("C_1"),size=6)+
  annotate("text",x = 7 ,y = 7, label = TeX("C_3"),size=6)+
  annotate("text",x = 4.5 ,y = 4.5, label = TeX("C_2"),size=6)+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))

png(paste0(plot_output_dir,'Figure_1_a.png'),
           width = 3 ,height=4,res=200,units='in')
ggarrange(
  p_A_true_ggplot, ncol=1, nrow=1,
  common.legend = FALSE, legend="bottom")
dev.off()

p_A_noisy_ggplot <- ggplot(A_noisy_value_plot, aes(Var1, Var2, fill = z)) +
  geom_tile() +
  xlab("") + ylab("") +
  ggtitle("Noisy realization")+
  scale_fill_distiller(palette = "PiYG")+
  theme_bw()+
  theme(plot.title=element_text(size=15,hjust = 0.5),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.title=element_text(size=15),
        strip.text.y = element_text(size=15,hjust=0,vjust = 1,angle=180,face="bold"))+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))+  annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 5.5, ymax = 8.5,
                                               alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
                                               size=3)+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.5, ymax = 3.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
  annotate("rect", xmin = 0.5, xmax = 8.5, ymin = 0.5, ymax = 8.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
  scale_x_discrete(expand = c(0,0.04)) +
  scale_y_discrete(expand = c(0,0.04))+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))


png(paste0(plot_output_dir,'Figure_1_b.png'),
    width = 3 ,height=4,res=200,units='in')
ggarrange(
  p_A_noisy_ggplot, ncol=1, nrow=1,
  common.legend = FALSE, legend="bottom")
dev.off()


K <- 13
complete_sol <- genlasso::fusedlasso(y=y0,D=Dmat,maxsteps=(K))
beta_hat <- complete_sol$beta[,(K+1)]

neg_current_inf_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                       D=Dmat,
                                                       c1=2, c2=3,
                                                       method="K", 
                                                       sigma=sigma,
                                                       K=K)

pos_current_inf_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                       D=Dmat,
                                                       c1=1, c2=2,
                                                       method="K",
                                                       sigma=sigma,
                                                       K=K)

pos_neg_current_inf_cc <- GFLassoInference::fusedlasso_inf(y=y,
                                                           D=Dmat,
                                                           c1=1, c2=3,
                                                           method="K",
                                                           sigma=sigma,
                                                           K=K)


# get the p-values in Figure 1 c
neg_current_inf_cc
pos_current_inf_cc
pos_neg_current_inf_cc


# wrangle data for Figure 1d
set.seed(123)
plot_x <- seq(3.3, 3.45, by = 0.0001)
hyun_set <- as.vector(neg_current_inf_cc$hyun_set)
plot_y_hyun <- g_single_truncation(plot_x, mean = 0, sd = pos_current_inf_cc$sd,
                                   a = hyun_set[1], b = hyun_set[2])
hyun_plot_df <- data.frame(plot_x,plot_y_hyun)

# plot union density
plot_x <- seq(1.5, 4, by = 0.001)
union_set_1 <- as.vector(neg_current_inf_cc$truncation_set[2])
union_set_2 <- as.vector(neg_current_inf_cc$truncation_set[3])
plot_y_union <- g_double_truncation(plot_x, mean = 0, sd = neg_current_inf_cc$sd,
                                    a = union_set_1[1], b = union_set_1[2],
                                    c = union_set_2[1], d = union_set_2[2])
union_plot_df <- data.frame(plot_x, plot_y_union)


manual_color_pal <- scales::hue_pal()(3)[2:3]

p_cond_density <- ggplot() +
  ### add hyun set
  geom_line(data=hyun_plot_df%>%
              filter(plot_x>=hyun_set[1]-0.001,plot_x<=hyun_set[2]+0.001),
            aes(x= plot_x, y = plot_y_hyun))+
  geom_ribbon(data=hyun_plot_df,
              aes(x= plot_x, y = plot_y_hyun,ymax=plot_y_hyun,ymin=0,
                  fill="red"),alpha=0.4)+
  geom_line(data=union_plot_df,
            aes(x= plot_x, y = plot_y_union))+
  geom_ribbon(data=union_plot_df,
              aes(x= plot_x, y = plot_y_union,ymax=plot_y_union,ymin=0,
                  fill="blue"),alpha=0.4)+
  ylab('Density')+
  xlab(TeX('$\\nu^TY$'))+
  ggtitle(TeX('Conditional null distributions'))+
  theme_bw()+
  theme(plot.title=element_text(size=15,hjust = 0.5),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        legend.position="bottom",
        # axis.ticks.y=element_blank(),
        legend.title=element_text(size=15),
        legend.text = element_text(size=15),
        axis.title=element_text(size=15),
        strip.text.y = element_text(size=15,hjust=0,vjust = 1,angle=180,face="bold"))  +
  geom_vline(xintercept =  neg_current_inf_cc$test_stats,
             color = 'black', size= 1, linetype = 2)+
  xlim(1.5,4)+
  scale_y_sqrt()+
  scale_fill_manual(name = " ", values=rev(manual_color_pal),
                    labels = unname(rev(TeX(c('$p_{Hyun}$','$p_{\\hat{C}_{1},\\hat{C}_{2}}')))) )+
  theme(legend.position="right")


png(paste0(plot_output_dir,'Figure_1_d.png'),
    width = 9, height=2,res=300,units='in')
p_cond_density
dev.off()

