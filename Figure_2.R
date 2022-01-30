### load relevant libraries
library(genlasso)
library(latex2exp)
library(GFLassoInference)
library(lattice)
library(ggpubr)
library(RColorBrewer)

# Remember to change your input and output dir!
input_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/"
plot_output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/"

lev1 <- 0 # mean for group 1
delta <- 3 # mean for group 2; null simulation
sigma <- 1 # level of noise
nn <- 8 # grid size; vary 5 or 10
Dmat <-  genlasso::getD2d(nn, nn) # create D matrix 
A <- matrix(lev1, ncol=nn, nrow = nn)
# create three strips
A[1:(nn/3+1),1:(nn/3+1)] <- 1*delta 
A[(nn-2):(nn),(nn-2):(nn)] <- -1*delta
beta <- c(t(A))
K <- 13
set.seed(2005)
A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma) # create the noise vector
y <- c(t(A.noisy))

complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
beta_hat <- complete_sol$beta[,K]
resulting_cc <- complete_sol$pathobjs$i

A_plot <- data.frame(expand.grid(c(1:8),c(1:8)), z = y)

pos_segment_list <- list(which(resulting_cc==2),
                         which(resulting_cc==3)) 
pos_contrast_vec <- rep(0, times = length(resulting_cc))
pos_contrast_vec[pos_segment_list[[1]]] = 1/(length(pos_segment_list[[1]]))
pos_contrast_vec[pos_segment_list[[2]]] = -1/(length(pos_segment_list[[2]]))
v <- pos_contrast_vec
phi <- sum(v*y) #2.362945
new_y_no_change <- as.numeric(y-(sum(v*y)-phi)*v/sum(v*v))

# less change
phi <- 0
new_y_less_change <- as.numeric(y-(sum(v*y)-phi)*v/sum(v*v))
A_less_change <- data.frame(expand.grid(c(1:8),c(1:8)), z = new_y_less_change)
complete_sol_less_change <- genlasso::fusedlasso(y=new_y_less_change,D=Dmat,maxsteps=(K))
resulting_cc_less_change  <- complete_sol_less_change$pathobjs$i

# more pronounced change
phi <- -5
new_y_more_change <- as.numeric(y-(sum(v*y)-phi)*v/sum(v*v))
A_more_change <- data.frame(expand.grid(c(1:8),c(1:8)), z = new_y_more_change)
complete_sol_more_change <- genlasso::fusedlasso(y=new_y_more_change,D=Dmat,maxsteps=(K))
resulting_cc_more_change <-  complete_sol_more_change$pathobjs$i


A_plot$estimated_CC <- as.factor(resulting_cc)
A_less_change$estimated_CC <- as.factor(resulting_cc_less_change)
A_more_change$estimated_CC <- as.factor(resulting_cc_more_change)
table(A_more_change$estimated_CC, A_plot$estimated_CC)


### original A
p0_ggplot <- ggplot(A_plot, aes(Var1, Var2, fill = z)) +
  geom_tile() +
  xlab("") + ylab("") +
  ggtitle(TeX('Original data ($\\phi$=$3.36$)'))+
  scale_fill_distiller(palette = "PiYG")+
  theme_bw()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(expand = c(0,0.04)) +
  scale_y_discrete(expand = c(0,0.04))+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))+
  annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 5.5, ymax = 8.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.5, ymax = 3.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
  annotate("rect", xmin = 0.5, xmax = 8.5, ymin = 0.5, ymax = 8.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)

# perturbed with less change
p1_ggplot <- ggplot(A_less_change, aes(Var1, Var2, fill = z)) + geom_tile() +
  xlab("") + ylab("") +
  ggtitle(TeX('Peturbed data ($\\phi$=0)'))+
  scale_fill_distiller(palette = "PiYG")+
  theme_bw()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))+
annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 5.5, ymax = 8.5,
         alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
         size=3)+
  annotate("rect", xmin = 0.5, xmax = 8.5, ymin = 0.5, ymax = 8.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
  scale_x_discrete(expand = c(0,0.04)) +
  scale_y_discrete(expand = c(0,0.04))

p2_ggplot <-  ggplot(A_more_change, aes(Var1, Var2, fill = z)) + geom_tile() +
  xlab("") + ylab("") +
  ggtitle(TeX('Peturbed data ($\\phi$=-5)'))+
  scale_fill_distiller(palette = "PiYG")+
  theme_bw()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))+
  annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 5.5, ymax = 8.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 0.5, ymax = 3.5,
           alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
           size=3)+
annotate("rect", xmin = 0.5, xmax = 8.5, ymin = 0.5, ymax = 8.5,
         alpha = 0.1, fill = alpha("grey",0), colour = "grey50",
         size=3)+
  scale_x_discrete(expand = c(0,0.04)) +
  scale_y_discrete(expand = c(0,0.04))


pos_current_inf_cc <- fusedlasso_inf(y=y,
                                     D=Dmat,
                                     c1=2, 
                                     c2=3, 
                                     method="K", 
                                     sigma=sigma, 
                                     K=K) 


### Plot the panel (d)
union_set_1 <- as.vector(pos_current_inf_cc$truncation_set[1])
union_set_2 <- as.vector(pos_current_inf_cc$truncation_set[2])
union_set_3 <- as.vector(pos_current_inf_cc$truncation_set[3])

new_x_mid <- c(c(-10,union_set_1[2],union_set_2[1],union_set_2[2],union_set_3[1])+
                 c(union_set_1[2],union_set_2[1],union_set_2[2],union_set_3[1],70))/2

plot_rec <- data.frame(x=new_x_mid,
                       y=rep(c(1),each=5),
                       z=factor(c(1,2,1,2,1)),
                       w = c(diff(c(-10,union_set_1[2],union_set_2[1],union_set_2[2],union_set_3[1])),
                             70-union_set_3[1]))

cbPalette <- c( "#56B4E9","#E69F00", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")

phi_values <- ggplot(plot_rec,
                      aes(xmin = x - w/2 , xmax = x + w/2 , ymin = y, ymax = y + 1)) +
  geom_rect(aes(fill = z), colour = NA)+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  ylab("")+
  xlab(TeX('$\\phi$'))+
  ggtitle(TeX('Set of $\\phi$ for which $\\hat{C}_1, \\hat{C}_2 \\in CC_{13}(y\'(\\phi))$'))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.title.x=element_text(size=12),
        title=element_text(size=15))


png(paste0(plot_output_dir,'Figure_2_a.png'),
    width = 3 ,height=3.7,res=200,units='in')
ggarrange( 
  p0_ggplot, ncol=1, nrow=1,
  common.legend = FALSE, legend="bottom")
dev.off()

png(paste0(plot_output_dir,'Figure_2_b.png'),
    width = 3,height=3.7,res=200,units='in')
ggarrange( 
  p1_ggplot, ncol=1, nrow=1,
  common.legend = FALSE, legend="bottom")
dev.off()


png(paste0(plot_output_dir,'Figure_2_c.png'),
    width = 3 ,height=3.7,res=200,units='in')
ggarrange( 
  p2_ggplot, ncol=1, nrow=1,
  common.legend = FALSE, legend="bottom")
dev.off()


png(paste0(plot_output_dir,'Figure_2_d.png'),
    width = 9 ,height=1.5,res=300,units='in')

print(phi_values)

dev.off()





