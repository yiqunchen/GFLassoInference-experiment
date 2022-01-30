library(spdep)
library(maps)
library(maptools)
library(RColorBrewer)
library(latex2exp)
library(tidyverse)
library(igraph)
library(GFLassoInference)
library(genlasso)
library(scales)
manual_color_pal <- scales::hue_pal()(3)

## Input dir
input_dir  <- "~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/input/"
plot_output_dir <- '~/Desktop/dissertation/gflasso_project/GFLassoInference-experiment/plot_output/'

estimate_sigma <- function(y, cc_y){
  s <- split(y, cc_y)
  s_mean <- sapply(s,mean)
  rss <- sum((unlist(sapply(seq_along(s), function(i) s[[i]]-s_mean[[i]])))^2)
  s_sigma <- sqrt(rss/(length(y)-length(s)))
  return(s_sigma)
}

getDmat.from.adjmat <- function(adjmat, sparseMatrix = FALSE) {
  n = ncol(adjmat)
  if (sparseMatrix) {
    Dmat = Matrix(0, nrow = n^2, ncol = n, sparse = TRUE)
  }
  else {
    Dmat = matrix(0, nrow = n^2, ncol = n)
  }
  count = 1
  for (jj in 1:n) {
    for (kk in jj:n) {
      if (adjmat[jj, kk] == 1) {
        Dmat[count, jj] = 1
        Dmat[count, kk] = -1
        count = count + 1
      }
    }
  }
  Dmat = Dmat[1:(count - 1), ]
  return(Dmat)
}

## Create an adjacency matrix for the states in the US
usa.state <- maps::map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"), function(x) x[1])
usa.poly <- map2SpatialPolygons(usa.state, IDs=state.ID)
usa.nb <- poly2nb(usa.poly)
usa.adj.mat <- as.matrix(nb2mat(usa.nb, style="B"))
colnames(usa.adj.mat)<- rownames(usa.adj.mat)
# we don't do dc
dc_index <- which(colnames(usa.adj.mat)=="district of columbia")
usa.adj.mat <- usa.adj.mat[-dc_index, -dc_index]
# apparently 4 corners don't count
usa.adj.mat[which(rownames(usa.adj.mat)=='utah'),
            which(rownames(usa.adj.mat)=='new mexico')] <- 0
usa.adj.mat[which(rownames(usa.adj.mat)=='new mexico'),
            which(rownames(usa.adj.mat)=='utah')] <- 0
usa.adj.mat[which(rownames(usa.adj.mat)=='arizona'),
            which(rownames(usa.adj.mat)=='colorado')] <- 0
usa.adj.mat[which(rownames(usa.adj.mat)=='colorado'),
            which(rownames(usa.adj.mat)=='arizona')] <- 0

teen_birth <- read.csv(paste0(input_dir,'Teen Birth Rate by State.csv'))
teen_birth <- teen_birth %>%
  dplyr::filter(YEAR==2018) %>%
  dplyr::select(c(STATE,RATE))%>%
  mutate(log_rate=log(RATE))

colnames(teen_birth) <- c('states','rate','log_rate')

state_name <- read.csv(paste0(input_dir,'state_name_abbreviation.csv'))
state_name$State <- tolower(state_name$State)
teen_birth <- teen_birth %>%
  inner_join(state_name, by = c('states'='Abbreviation'))

all_states <- map_data("state")
teen_birth_plot <- teen_birth %>%
  inner_join(all_states, by = c("State"="region") )

# color ramp
yellow <- colorRampPalette(c("yellow","orange"))(30)
red <- colorRampPalette(c("orangered","darkred"))(30)

#### observed rate
obs_rate_teen <- ggplot(data = teen_birth_plot) +
  geom_polygon(aes(x = long, y = lat,
                   fill = exp(log_rate), group = group), color = "white") +
  coord_fixed(1.3) +
  guides(fill=FALSE)+
  theme_bw()+
  ylab('')+
  xlab('')+
  scale_fill_gradientn(colours=c(yellow, red))+
  guides(fill = guide_colourbar())+
  theme_void()+
  labs(fill='Birth Rate')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour = "white"),
        panel.background = element_blank(),
        legend.position="bottom",plot.margin=unit(c(-0.0,0,0,0), "null"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0))+
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1))


png(paste0(plot_output_dir,'Figure_6_a.png'),
    width = 9,height=5,res=300,units='in')
print(obs_rate_teen)
dev.off()


####
teen_birth <- teen_birth %>%
  filter(!State%in%c('hawaii','alaska'))
K <- 30
log_rate <- teen_birth$log_rate
# aux function
US_Dmat <- getDmat.from.adjmat(usa.adj.mat)
# solve for solution
GFL_sol <- genlasso::fusedlasso(y=log_rate,D=US_Dmat,maxsteps=K)
gfl_sol_staes <- data.frame(
  log_rate = teen_birth$log_rate,
  State = teen_birth$State,
  sol_step_K = GFL_sol$beta[,K]
)


teen_birth_plot_sol <- teen_birth_plot %>%
  inner_join(gfl_sol_staes, by = 'State')

fused_k_rate  <- ggplot(data = teen_birth_plot_sol) +
  geom_polygon(aes(x = long, y = lat,
                   fill = exp(sol_step_K), group = group), color = "white") +
  coord_fixed(1.3) +
  guides(fill=FALSE)+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour = "white"),
        panel.background = element_blank())+
  scale_fill_gradientn(colours=c(yellow, red))+
  labs(fill='Birth Rate')+
  guides(fill = guide_colourbar())+
  annotate(geom="label", x=-110, y=42, label=TeX("$\\hat{C}_5$"),
           fill="black", colour="white",size=7)+
  annotate(geom="label", x=-85, y=32, label=TeX("$\\hat{C}_3$"),
           fill="black", colour="white",size=6)+
  annotate(geom="label", x=-72, y=44,label=TeX("$\\hat{C}_1$"),
           fill="black", colour="white",size=5)+
  annotate(geom="label", x=-78, y=40.5, label=TeX("$\\hat{C}_4$"),
           fill="black", colour="white",size=3)+
  annotate(geom="label", x=-75, y=38, label=TeX("$\\hat{C}_2$"),
           fill="black", colour="white",size=3)+
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1))+
  theme(legend.position="bottom",plot.margin=unit(c(-0.0,0,0,0), "null"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0))

png(paste0(plot_output_dir,'Figure_6_b.png'),
    width = 9,height=5,res=300,units='in')
print(fused_k_rate)
dev.off()

#### Figure 5c and 5d -- test #### 
sigma_hat_K <- estimate_sigma(gfl_sol_staes$log_rate, 
                               as.factor(gfl_sol_staes$sol_step_K))


p_val_c1_c2 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=1,
                                                c2=2,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c1_c3 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=1, 
                                                c2=3,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c1_c4 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=1, 
                                                c2=4,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c1_c5 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=1, 
                                                c2=5,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c2_c3 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=2, 
                                                c2=3,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c2_c4 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=2, 
                                                c2=4,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c2_c5 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=2, 
                                                c2=5,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c3_c4 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=3, 
                                                c2=4,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c3_c5 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=3, 
                                                c2=5,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

p_val_c4_c5 <- GFLassoInference::fusedlasso_inf(y=gfl_sol_staes$log_rate,
                                                D=US_Dmat,
                                                c1=4, 
                                                c2=5,
                                                method="K",
                                                sigma=sigma_hat_K,
                                                K=K,
                                                compute_ci=TRUE)

aggr_ci_info <- function(p_val_segment_cc){
  # ah hoc; compute CI correspondong to p Hyun
  ci_p_hyun <- GFLassoInference::compute_CI(
    vTy= p_val_segment_cc$test_stats,
    vTv=(p_val_segment_cc$sd)^2, sigma=1,
    truncation = data.frame(p_val_segment_cc$hyun_set, ncol=2),
    alpha=0.05)
  # aggregate results
  result_aggr <- data.frame("p_union" = p_val_segment_cc$Union, "p_naive" = p_val_segment_cc$Naive, 
                            "p_hyun" = p_val_segment_cc$Hyun, "test_stat" = p_val_segment_cc$test_stats,
                            "hyun_lcb" = ci_p_hyun[1], "hyun_ucb" = ci_p_hyun[2],
                            "union_lcb" = p_val_segment_cc$CI_result[1], "union_ucb" = p_val_segment_cc$CI_result[2],
                            "naive_lcb" = p_val_segment_cc$test_stats-1.96*p_val_segment_cc$sd,
                            "naive_ucb" = p_val_segment_cc$test_stats+1.96*p_val_segment_cc$sd)
  return(result_aggr)
}
# Plot CI collection
list_results <- list(p_val_c1_c2, p_val_c1_c3, p_val_c1_c4, p_val_c1_c5, p_val_c2_c3,
                     p_val_c2_c4, p_val_c2_c5, p_val_c3_c4, p_val_c3_c5, p_val_c4_c5)

result_list_df <-  dplyr::bind_rows(lapply(list_results,function(x)aggr_ci_info(x)))
result_list_df$x_plot <- c(1:10)

#### aggregate results
result_list_df_ci_plot <- result_list_df %>% 
  select(-c(p_union,p_naive,p_hyun)) %>%
  pivot_longer(cols=-c(x_plot,test_stat)) %>% 
  separate(name, c("type", "cb_type"), "_") %>%
  pivot_wider(names_from=cb_type,values_from = value) %>%
  mutate(sig_0 = ((lcb<0)&(ucb<0))|((lcb>0)&(ucb>0))) %>%
  mutate(type = factor(type,levels=c("naive","hyun","union")))

### create x axis ticks 
all_comb_names <- t(utils::combn(5, 2))
tick_labels <- (paste0("$\\bar{\\beta}_{",all_comb_names[,1],"}","-",
                       "\\bar{\\beta}_{",all_comb_names[,2],"}$"))

p_ci_teen <- ggplot(data=result_list_df_ci_plot %>%
                          mutate(type=factor(type, levels = c("naive","union","hyun"))),
 aes(x=as.factor(x_plot), y=test_stat, colour=type))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_pointrange(aes(ymin = lcb, ymax = ucb),
    fatten = 0.5, size = 0.7,
    position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0))+
  coord_cartesian(ylim=c(-4,4))+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="bottom",
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15,hjust = 0),
        legend.margin=margin(-10, 0, 0, 0))+
  scale_color_manual( " ", values = c("#F8766D","#619CFF", "#00BA38"),
                      breaks=c("naive","union","hyun"),
                      labels = unname(TeX(c('$p_{Naive}$-based CI',
                                            '$p_{\\hat{C}_{1},\\hat{C}_{2}}$-based CI',
                                            '$p_{Hyun}$-based CI' )))) +
  scale_x_discrete(breaks=c(1:10), labels=unname(TeX(tick_labels)))

png(paste0(plot_output_dir,'Figure_6_d.png'),
    width = 12,height=3,res=200,units='in')
print(p_ci_teen)
dev.off()

