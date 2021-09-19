input_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/"
plot_output_dir <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/"

lev1 <- 0 # mean for group 1
delta <- 3 # mean for group 2
sigma <- 1 # level of noise
A <- matrix(c(rep(0,100),rep(delta, 40), rep(0, 60)),nrow=1)
n <- length(A)
# create true mean vector
beta <- c(t(A))
# sample observation
set.seed(1234)
y <- beta+rnorm(n)
# create plot
plot_toy <- data.frame(x=seq(1,n,by=1), beta=beta, y=y)
toy_plot_1 <- ggplot(plot_toy)+
  scale_colour_grey()+
  geom_point(aes(x=x,y=y),colour='grey')+
  geom_line(aes(x=x,y=beta),colour='black')+
  theme_classic()+
  ylab('')+
  xlab('')+
  #ggtitle('Simulated example for the 1D fused lasso problem')+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        legend.position="bottom",
        legend.title = element_text(size=15),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=12),
        legend.text = element_text(size=15,hjust = 0),
        axis.title=element_text(size=15))

png(paste0(plot_output_dir,'Figure_3_a.png'),
    width = 7,height=7,res=200,units='in')

toy_plot_1

dev.off()
